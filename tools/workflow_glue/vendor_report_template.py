"""Convert a wf-template report into embedded fastcat C report assets."""

import argparse
import gzip
from pathlib import Path
import re


# ezcharts emits a complete static Bokeh document with example plot data baked
# into the JSON blob. fastcat keeps that plot/layout JSON and swaps only the
# data sources at report-open time, avoiding surgery on Bokeh internals during
# the C runtime. This patch function is injected into the Bokeh embed_document().
PATCH_FUNCTION = r"""
        function fc_fastcat_map(entries) {
          return {"type": "map", "entries": entries};
        }

        function fc_fastcat_set_source(fig, entries) {
          if (!fig || !fig.attributes || !fig.attributes.renderers || fig.attributes.renderers.length === 0) {
            return;
          }
          const renderer = fig.attributes.renderers[0];
          const source = renderer.attributes && renderer.attributes.data_source;
          if (source && source.attributes) {
            source.attributes.data = fc_fastcat_map(entries);
          }
        }

        function fc_fastcat_set_subtitle(fig, text) {
          const above = fig && fig.attributes && fig.attributes.above;
          if (above && above.length > 0 && above[0].attributes) {
            above[0].attributes.text = text;
          }
        }

        function fc_patch_fastcat_report(docs, data) {
          const doc = docs[Object.keys(docs)[0]];
          if (!doc || !Array.isArray(doc.roots)) {
            return;
          }
          const figures = doc.roots.filter((root) => root && root.type === "object" && root.name === "Figure");
          if (figures.length < 3) {
            return;
          }
          fc_fastcat_set_source(figures[0], [
            ["top", data.quality.top],
            ["left", data.quality.left],
            ["right", data.quality.right],
          ]);
          fc_fastcat_set_subtitle(figures[0], data.quality.subtitle);

          fc_fastcat_set_source(figures[1], [
            ["top", data.length.top],
            ["left", data.length.left],
            ["right", data.length.right],
          ]);
          fc_fastcat_set_subtitle(figures[1], data.length.subtitle);

          fc_fastcat_set_source(figures[2], [
            ["x", data.yield.x],
            ["y", data.yield.y],
          ]);
          fc_fastcat_set_subtitle(figures[2], data.yield.subtitle);
        }
"""


def argparser():
    parser = argparse.ArgumentParser(
        description="patch wf-template HTML and emit embedded C assets",
        add_help=False)
    parser.add_argument("--input", required=True, help="example wf-template report HTML")
    parser.add_argument("--template-output", help="optional patched runtime template HTML")
    parser.add_argument("--header", required=True, help="output header path")
    parser.add_argument("--source", required=True, help="output C source path")
    parser.add_argument("--symbol", default="fc_report_template_html", help="symbol base name")
    return parser


def replace_embed_call(html):
    needle = "        function embed_document(root) {\n"
    if needle not in html:
        raise RuntimeError("could not find Bokeh embed_document function")
    html = html.replace(needle, PATCH_FUNCTION + needle, 1)

    # patch the Bokeh document immediately before Bokeh embeds it, while the
    # document is still plain JSON rather than live browser objects
    old = "        root.Bokeh.embed.embed_items(docs_json, render_items);"
    new = (
        "        const fastcat_report_data = __FASTCAT_REPORT_DATA__;\n"
        "        const fastcat_docs = JSON.parse(docs_json);\n"
        "        fc_patch_fastcat_report(fastcat_docs, fastcat_report_data);\n"
        "        root.Bokeh.embed.embed_items(JSON.stringify(fastcat_docs), render_items);"
    )
    if old not in html:
        raise RuntimeError("could not find Bokeh embed_items call")
    return html.replace(old, new, 1)


def remove_tag_block_containing(html, tag, marker):
    marker_pos = html.find(marker)
    if marker_pos < 0:
        return html

    tag_open = f"<{tag}"
    tag_close = f"</{tag}>"
    start = html.rfind(tag_open, 0, marker_pos)
    if start < 0:
        raise RuntimeError(f"could not find <{tag}> containing {marker!r}")

    pos = start
    depth = 0
    while True:
        next_open = html.find(tag_open, pos)
        next_close = html.find(tag_close, pos)
        if next_close < 0:
            raise RuntimeError(f"could not find closing </{tag}> for {marker!r}")
        if next_open >= 0 and next_open < next_close:
            depth += 1
            pos = next_open + len(tag_open)
            continue
        depth -= 1
        pos = next_close + len(tag_close)
        if depth == 0:
            return html[:start] + html[pos:]


def remove_versions_panel(html):
    html = html.replace(
        '<div class=""><a class="dropdown-item" href="#versions">Versions</a></div>',
        "",
    )
    return remove_tag_block_containing(html, "section", 'id="versions"')


def replace_params_table(html):
    pattern = re.compile(
        r'(<section[^>]*id="parameters"[^>]*>.*?<h2[^>]*>)Workflow parameters'
        r'(</h2>.*?<tbody>\s*).*?(\s*</tbody>)',
        flags=re.DOTALL,
    )
    html, count = pattern.subn(
        r"\1Run configuration\2__FASTCAT_RUN_CONFIG_ROWS__\3",
        html,
        count=1,
    )
    if count != 1:
        raise RuntimeError("could not find workflow parameters table")
    return html


def apply_fastcat_template_wording(html):
    html = html.replace("fastcat sequencing report", "Fastcat Sequencing QC Report")
    html = html.replace(
        "Results generated through the fastcat workflow provided by Oxford Nanopore Technologies.",
        "Sequencing QC generated by fastcat.",
    )
    html = remove_versions_panel(html)
    html = replace_params_table(html)
    html = html.replace(" nextflow workflow (", " (")
    return html


def bytes_to_c_array(data, width=16):
    parts = [f"0x{b:02x}" for b in data]
    lines = []
    for i in range(0, len(parts), width):
        lines.append("    " + ", ".join(parts[i : i + width]))
    return ",\n".join(lines)


def emit_c_assets(html, header_path, source_path, symbol):
    raw = html.encode("utf-8")
    compressed = gzip.compress(raw, compresslevel=9, mtime=0)
    array_payload = bytes_to_c_array(compressed)

    guard = f"{symbol.upper()}_DATA_H"
    header_text = f"""#ifndef {guard}
#define {guard}

#include <stddef.h>

extern const unsigned char {symbol}_gz[];
extern const size_t {symbol}_gz_len;
extern const size_t {symbol}_len;

#endif
"""

    source_text = f"""#include "{header_path.name}"

/*
 * this template is emitted as gzip-compressed bytes rather than a C string so
 * arbitrary HTML, JavaScript, CSS, quotes, backslashes, and embedded data do
 * not need fragile C string escaping, and the executable does not carry the
 * full uncompressed report template in .rodata
 */

const size_t {symbol}_gz_len = {len(compressed)};
const size_t {symbol}_len = {len(raw)};

const unsigned char {symbol}_gz[] = {{
{array_payload}
}};
"""

    header_path.parent.mkdir(parents=True, exist_ok=True)
    source_path.parent.mkdir(parents=True, exist_ok=True)
    header_path.write_text(header_text, encoding="utf-8")
    source_path.write_text(source_text, encoding="utf-8")
    print(f"wrote {header_path}")
    print(f"wrote {source_path}")


def main(args):
    html = Path(args.input).read_text(encoding="utf-8")
    html = apply_fastcat_template_wording(html)
    html = replace_embed_call(html)

    if args.template_output:
        template_path = Path(args.template_output)
        template_path.parent.mkdir(parents=True, exist_ok=True)
        template_path.write_text(html, encoding="utf-8")

    emit_c_assets(html, Path(args.header), Path(args.source), args.symbol)
    return 0