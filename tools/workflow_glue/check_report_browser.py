"""Fail on report console errors or broken interactive plot bindings."""

import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import subprocess
import tempfile
import threading


BROWSER_CHECK_SCRIPT = r"""
<script>
window.addEventListener("load", () => {
    window.setTimeout(() => {
        const root = document.documentElement;
        const finish = (status) => {
            root.setAttribute("data-fastcat-browser-check", status);
        };
        try {
            const roots = Array.from(
                document.querySelectorAll("[data-root-id]"));
            if (roots.length === 0) {
                throw new Error("report contains no Bokeh roots");
            }
            if (roots.some((plotRoot) => plotRoot.childElementCount === 0)) {
                throw new Error("one or more Bokeh roots did not render");
            }

            const doc = Bokeh.documents[0];
            if (doc == null) {
                throw new Error("Bokeh document is missing");
            }
            const coverageTabs = Array.from(
                document.querySelectorAll("[data-fastcat-coverage-target]"));
            if (coverageTabs.length < 2) {
                finish("not-applicable");
                return;
            }

            const source = doc.get_model_by_name("fastcat_coverage_source");
            const hoverSource = doc.get_model_by_name("fastcat_coverage_hover_source");
            const targets = doc.get_model_by_name("fastcat_coverage_targets");
            const subtitle = doc.get_model_by_name("fastcat_coverage_subtitle");
            const xRange = doc.get_model_by_name("fastcat_coverage_x_range");
            const yRange = doc.get_model_by_name("fastcat_coverage_y_range");
            const hover = doc.get_model_by_name("fastcat_coverage_hover");
            if ([source, hoverSource, targets, subtitle, xRange, yRange, hover].some((model) => model == null)) {
                throw new Error("coverage selector models are missing");
            }
            const accuracyTab = document.querySelector(
                '[data-fastcat-accuracy-tab="1"]');
            const alignmentSection = document.getElementById("fastcat-alignment-section");
            const hasAccuracyPlot = accuracyTab != null
                && alignmentSection != null
                && alignmentSection.style.display !== "none";
            if (hasAccuracyPlot) {
                accuracyTab.click();
            }
            const hoverHasDataValues =
                hover.renderers.length === 1
                && hover.renderers[0].data_source === hoverSource
                && JSON.stringify(hover.tooltips) === JSON.stringify([
                    ["Minimum Depth", "@x{0}"],
                    ["Reference Positions Covered / %", "@y{0.0}"],
                ]);
            if (!hoverHasDataValues) {
                throw new Error("coverage hover is not bound to data values");
            }

            const index = 1;
            coverageTabs[index].click();
            window.setTimeout(() => {
                try {
                    const sameValues = (actual, expected) =>
                        JSON.stringify(Array.from(actual))
                        === JSON.stringify(Array.from(expected));
                    const valid =
                        sameValues(source.data.x, targets.data.x[index])
                        && sameValues(source.data.y, targets.data.y[index])
                        && sameValues(hoverSource.data.x, targets.data.hover_x[index])
                        && sameValues(hoverSource.data.y, targets.data.hover_y[index])
                        && sameValues(hoverSource.data.width, targets.data.hover_width[index])
                        && sameValues(hoverSource.data.height, targets.data.hover_height[index])
                        && subtitle.text === targets.data.subtitle[index]
                        && xRange.end === targets.data.x_end[index]
                        && yRange.end === targets.data.y_end[index];
                    if (!valid) {
                        throw new Error(
                            "coverage selector did not update its plot models");
                    }
                    if (hasAccuracyPlot) {
                        const accuracyFigure = doc.get_model_by_name(
                            "fastcat_accuracy_qscore");
                        const accuracyView = accuracyFigure == null
                            ? null : Bokeh.index[accuracyFigure.id];
                        if (accuracyView == null) {
                            throw new Error(
                                "alignment-concordance plot view is missing");
                        }
                    }
                    finish("ok");
                } catch (error) {
                    console.error("fastcat browser check failed", error);
                    finish("failed");
                }
            }, 100);
        } catch (error) {
            console.error("fastcat browser check failed", error);
            finish("failed");
        }
    }, 1000);
});
</script>
"""


def argparser():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("report", type=Path)
    parser.add_argument("--browser", required=True)
    return parser


def main(args):
    report = args.report.resolve()
    if not report.is_file():
        raise RuntimeError(f"report does not exist: {report}")

    report_html = report.read_text(encoding="utf-8")
    body_end = report_html.rfind("</body>")
    if body_end < 0:
        raise RuntimeError(f"report has no closing body tag: {report}")
    report_html = (
        report_html[:body_end]
        + BROWSER_CHECK_SCRIPT
        + report_html[body_end:])

    with tempfile.TemporaryDirectory() as temp_dir:
        instrumented = Path(temp_dir, report.name)
        instrumented.write_text(report_html, encoding="utf-8")
        handler = partial(SimpleHTTPRequestHandler, directory=temp_dir)
        server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            url = f"http://127.0.0.1:{server.server_port}/{instrumented.name}"
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
        if "CONSOLE" in line
        and any(token in line for token in ("Uncaught", "Error", "ERROR"))
    ]
    if console_errors:
        raise RuntimeError(
            "browser console errors:\n" + "\n".join(console_errors))
    markers = (
        'data-fastcat-browser-check="ok"',
        'data-fastcat-browser-check="not-applicable"',
    )
    if not any(marker in result.stdout for marker in markers):
        raise RuntimeError(
            "browser report check did not complete successfully")
    return 0
