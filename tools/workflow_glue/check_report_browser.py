"""Fail when a generated HTML report produces browser console errors."""

import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import subprocess
import threading


def argparser():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("report", type=Path)
    parser.add_argument("--browser", required=True)
    return parser


def main(args):
    report = args.report.resolve()
    if not report.is_file():
        raise RuntimeError(f"report does not exist: {report}")

    handler = partial(SimpleHTTPRequestHandler, directory=str(report.parent))
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        url = f"http://127.0.0.1:{server.server_port}/{report.name}"
        result = subprocess.run(
            [
                args.browser,
                "--headless",
                "--no-sandbox",
                "--disable-dev-shm-usage",
                "--enable-logging=stderr",
                "--v=1",
                "--virtual-time-budget=5000",
                "--dump-dom",
                url,
            ],
            capture_output=True,
            check=False,
            text=True,
        )
    finally:
        server.shutdown()
        thread.join()
        server.server_close()

    if result.returncode != 0:
        raise RuntimeError(
            f"headless browser exited with status {result.returncode}:\n{result.stderr}")

    console_errors = [
        line for line in result.stderr.splitlines()
        if "ERROR:CONSOLE" in line
    ]
    if console_errors:
        raise RuntimeError(
            "browser console errors:\n" + "\n".join(console_errors))
    return 0
