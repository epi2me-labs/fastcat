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
            const models = Array.from(doc._all_models.values());
            const numericScales = models.filter(
                (model) => model.type === "LinearAxis" || model.type === "ColorBar");
            if (numericScales.length === 0
                    || numericScales.some(
                        (scale) => scale.formatter?.type !== "CustomJSTickFormatter")) {
                throw new Error("one or more numeric report scales lack adaptive formatting");
            }
            const formatter = numericScales[0].formatter;
            const expectLabels = (ticks, expected) => {
                const actual = formatter.doFormat(ticks, {loc: 0});
                if (JSON.stringify(actual) !== JSON.stringify(expected)) {
                    throw new Error(
                        "adaptive tick formatting mismatch for "
                        + JSON.stringify(ticks) + ": "
                        + JSON.stringify(actual) + " != "
                        + JSON.stringify(expected));
                }
            };
            expectLabels(
                [0, 500, 1000],
                ["0", "500", "1000"]);
            expectLabels(
                [0, 0.5, 1, 1.5, 2],
                ["0.0", "0.5", "1.0", "1.5", "2.0"]);
            expectLabels(
                [0, 0.05, 0.1, 0.15, 0.2],
                ["0.00", "0.05", "0.10", "0.15", "0.20"]);
            expectLabels(
                [0, 10000, 15000, 20000],
                ["0.0 × 10⁴", "1.0 × 10⁴", "1.5 × 10⁴", "2.0 × 10⁴"]);
            expectLabels(
                [0, NaN, 1],
                ["0", "", "1"]);

            const heatmapDetails = [
                ["fastcat_length_qscore", "Read length / kb"],
                ["fastcat_speed_qscore", "Translocation speed / bp/s"],
            ];
            const accuracyFigure = doc.get_model_by_name("fastcat_accuracy_qscore");
            if (accuracyFigure != null) {
                heatmapDetails.push([
                    "fastcat_accuracy_qscore", "Alignment concordance / Q"]);
            }
            for (const [name, xLabel] of heatmapDetails) {
                const figure = doc.get_model_by_name(name);
                const renderer = figure?.renderers?.[0];
                const source = renderer?.data_source;
                const hover = figure?.toolbar?.tools?.find(
                    (tool) => tool.renderers?.length === 1 && tool.renderers[0] === renderer);
                const fields = source == null ? [] : Object.keys(source.data).sort();
                const expectedTooltips = [
                    [xLabel, "@x{custom}"],
                    ["Quality", "@y{custom}"],
                    ["Count", "@image{0}"],
                ];
                if (
                        figure == null || figure.renderers.length !== 1
                        || renderer.glyph?.type !== "Image"
                        || JSON.stringify(fields) !== JSON.stringify(["dh", "dw", "image", "x", "y"])
                        || hover == null
                        || JSON.stringify(hover.tooltips) !== JSON.stringify(expectedTooltips)) {
                    throw new Error(`${name} does not use direct image hover data`);
                }
                const cell = {index: 0, i: 0, j: 0, flat_index: 0};
                const image = source.data.image[0];
                if (image == null) {
                    continue;
                }
                const x = Number(hover.formatters["@x"].format(
                    source.data.x[0], "custom", {image_index: cell}));
                const y = Number(hover.formatters["@y"].format(
                    source.data.y[0], "custom", {image_index: cell}));
                const expectedX = source.data.x[0] + source.data.dw[0] / image.shape[1] / 2;
                const expectedY = source.data.y[0] + source.data.dh[0] / image.shape[0] / 2;
                if (Math.abs(x - expectedX) > 0.001 || Math.abs(y - expectedY) > 0.1) {
                    throw new Error(`${name} hover coordinates do not use image cells`);
                }
                if (Math.abs(source.data.dh[0] / image.shape[0] - 0.5) > 1e-9) {
                    throw new Error(`${name} embeds quality bins outside its visible range`);
                }
                if (name === "fastcat_accuracy_qscore"
                        && Math.abs(source.data.dw[0] / image.shape[1] - 0.5) > 1e-9) {
                    throw new Error(`${name} embeds concordance bins outside its visible range`);
                }
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
