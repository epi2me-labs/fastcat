"""Local Bokeh-backed replacement for the fastcat seqsummary component."""

import os

from bokeh.models import CustomJSTickFormatter, NumeralTickFormatter, Title
import ezcharts as ezc
from ezcharts.components.ezchart import EZChart
from ezcharts.layout.base import Snippet
from ezcharts.layout.snippets import Grid, Tabs
from ezcharts.plots.util import empty_plot
import numpy as np
import pandas as pd
from .bokeh_style import ReportFontSizes, ReportFonts


QUALITY_BINWIDTH = 0.2
QUALITY_XMIN = 0.0
SUMMARY_PLOT_HEIGHT = "500px"
VISIBLE_QUANTILE = 0.98
VISIBLE_PADDING = 0.10
READ_LENGTH_BINS = 100
YIELD_THINNING = 100
QUALITY_XROUND = 10.0
LENGTH_XROUND_KB = 50.0


SCIENTIFIC_SUPERSCRIPT_JS = """
const superscripts = {
  "-": "\u207b",
  "0": "\u2070",
  "1": "\u00b9",
  "2": "\u00b2",
  "3": "\u00b3",
  "4": "\u2074",
  "5": "\u2075",
  "6": "\u2076",
  "7": "\u2077",
  "8": "\u2078",
  "9": "\u2079",
};

if (!Number.isFinite(tick)) {
  return "";
}

const value = Number(tick);
const abs_value = Math.abs(value);
if (abs_value === 0) {
  return "0";
}

const max_abs_tick = ticks
  .map((v) => Math.abs(Number(v)))
  .filter((v) => Number.isFinite(v))
  .reduce((acc, v) => Math.max(acc, v), 0);
const use_scientific = max_abs_tick >= 1e4;

if (!use_scientific) {
  if (Math.abs(value - Math.round(value)) < 1e-9) {
    return Math.round(value).toString();
  }
  return value.toFixed(2).replace(/\\.0+$/, "").replace(/(\\.\\d*[1-9])0+$/, "$1");
}

const exponent = Math.floor(Math.log10(max_abs_tick));
const mantissa = value / Math.pow(10, exponent);
let mantissa_text = mantissa
  .toFixed(1);
const exponent_text = exponent
  .toString()
  .split("")
  .map((ch) => superscripts[ch] || ch)
  .join("");
return `${mantissa_text} \u00d7 10 ${exponent_text}`;
"""


def load_histogram(hist_dir, dtype):
    """Load one fastcat histogram file."""
    dt = float
    if "length" in dtype:
        dt = int

    return pd.read_csv(
        os.path.join(hist_dir, f"{dtype}.hist"),
        sep="\t",
        names=["start", "end", "count"],
        dtype={"start": dt, "end": dt, "count": int},
    )


def sum_hists(hists):
    """Sum histogram dataframes based on their intervals."""
    return (
        pd.concat(hists)
        .reset_index(drop=True)
        .sort_values(by="start")
        .groupby(["start", "end"], as_index=False)
        .sum()
    )


def load_summary_histograms(path):
    """Load quality and length histograms, including unmapped if available."""
    quality = load_histogram(path, "quality")
    length = load_histogram(path, "length")

    try:
        quality = sum_hists((quality, load_histogram(path, "quality.unmap")))
    except Exception:
        pass

    try:
        length = sum_hists((length, load_histogram(path, "length.unmap")))
    except Exception:
        pass

    return length, quality


def histogram_quantile(hist, q):
    """Calculate a weighted quantile from a histogram using bin starts."""
    if hist.empty:
        return 0.0
    cumsum = np.cumsum(hist["count"].to_numpy(dtype=float))
    total = cumsum[-1]
    if not (total > 0.0):
        return 0.0
    threshold = q * total
    idx = int(np.searchsorted(cumsum, threshold, side="left"))
    idx = min(max(idx, 0), len(hist) - 1)
    return float(hist.iloc[idx]["start"])


def histogram_median(hist):
    return histogram_quantile(hist, 0.5)


def histogram_mean_midpoint(hist):
    if hist.empty:
        return 0.0
    mid = 0.5 * (hist["start"].to_numpy(dtype=float) + hist["end"].to_numpy(dtype=float))
    weights = hist["count"].to_numpy(dtype=float)
    total = weights.sum()
    if not (total > 0.0):
        return 0.0
    return float(np.average(mid, weights=weights))


def histogram_visible_max(hist):
    """Use the padded visible quantile cutoff for the x-range."""
    if hist.empty:
        return 1.0
    q95 = histogram_quantile(hist, VISIBLE_QUANTILE)
    xmax = q95 * (1.0 + VISIBLE_PADDING)
    fallback = float(hist["end"].max())
    if not np.isfinite(xmax) or xmax <= 0.0:
        xmax = fallback
    if not np.isfinite(xmax) or xmax <= 0.0:
        xmax = 1.0
    return xmax


def _round_up_to(value, step):
    if not np.isfinite(value) or value <= 0.0:
        return step
    return float(np.ceil(value / step) * step)


def _format_sigfig(value):
    return f"{value:.3g}"


def _format_length_stat(value):
    value = float(value)
    if value < 1000.0:
        return f"{int(round(value))}"
    return f"{value / 1000.0:.3g} kb"


def apply_scientific_count_formatter(axis):
    """Format large count axes as mantissa x 10^n with unicode superscripts."""
    axis.formatter = CustomJSTickFormatter(code=SCIENTIFIC_SUPERSCRIPT_JS)


def apply_fixed_precision_formatter(axis, fmt):
    """Apply a fixed numeral formatter."""
    axis.formatter = NumeralTickFormatter(format=fmt)


def read_quality_plot(hist, color=None):
    """Create the read-quality histogram with a quantile-based x-range."""
    xmax = _round_up_to(histogram_visible_max(hist), QUALITY_XROUND)
    mean_val = histogram_mean_midpoint(hist)
    median_val = histogram_median(hist)
    plot = ezc.histplot(
        data=hist["start"],
        weights=hist["count"],
        binwidth=QUALITY_BINWIDTH,
        binrange=(QUALITY_XMIN, xmax),
        color=color,
    )
    plot._fig.xaxis.axis_label = "Quality"
    plot._fig.yaxis.axis_label = "Number of reads"
    plot._fig.x_range.start = QUALITY_XMIN
    plot._fig.x_range.end = xmax
    apply_scientific_count_formatter(plot._fig.yaxis[0])
    plot._fig.add_layout(
        Title(text=f"Mean: {mean_val:.1f}. Median: {median_val:.1f}", text_font_size=ReportFontSizes.SUBTITLE, text_font=ReportFonts.PLOT),
        "above",
    )
    return plot


def read_length_plot(hist, color=None):
    """Create the read-length histogram with a quantile-based x-range."""
    xmax_kb = _round_up_to(histogram_visible_max(hist) / 1000.0, LENGTH_XROUND_KB)
    mean_length = histogram_mean_midpoint(hist)
    median_length = histogram_median(hist)
    min_length = int(hist["start"].min()) if not hist.empty else 0
    max_length = int(hist["start"].max()) if not hist.empty else 0

    plot = ezc.histplot(
        data=hist["start"] / 1000.0,
        weights=hist["count"],
        bins=READ_LENGTH_BINS,
        binrange=(0.0, xmax_kb),
        color=color,
    )
    plot._fig.xaxis.axis_label = "Read length / kb"
    plot._fig.yaxis.axis_label = "Number of reads"
    plot._fig.x_range.start = 0.0
    plot._fig.x_range.end = xmax_kb
    apply_scientific_count_formatter(plot._fig.yaxis[0])
    plot._fig.add_layout(
        Title(
            text=(
                f"Mean: {_format_length_stat(mean_length)}. "
                f"Median: {_format_length_stat(median_length)}. "
                f"Min: {_format_length_stat(min_length)}. "
                f"Max: {_format_length_stat(max_length)}"
            ),
            text_font_size=ReportFontSizes.SUBTITLE,
            text_font=ReportFonts.PLOT,
        ),
        "above",
    )
    return plot


def compute_n50_kb(hist):
    """Calculate N50 directly from the histogram."""
    if hist.empty:
        return 0.0
    lengths = hist["start"].to_numpy(dtype=float)
    weights = hist["count"].to_numpy(dtype=float)
    total_yield = float(np.sum(lengths * weights))
    if not (total_yield > 0.0):
        return 0.0
    cumsum = np.cumsum((lengths * weights)[::-1])[::-1]
    idx = len(lengths) - int(np.searchsorted(cumsum[::-1], total_yield * 0.5, side="left")) - 1
    idx = min(max(idx, 0), len(lengths) - 1)
    return lengths[idx] / 1000.0


def base_yield_plot(hist, color=None):
    """Create the anti-cumulative yield plot with the same x-range as length."""
    lengths = hist["start"].to_numpy(dtype=float)
    weights = hist["count"].to_numpy(dtype=float)
    cumsum = np.cumsum((lengths * weights)[::-1])[::-1]
    total_yield = cumsum[0] if cumsum.size else 0.0
    xmax_kb = _round_up_to(histogram_visible_max(hist) / 1000.0, LENGTH_XROUND_KB)

    df = pd.DataFrame({
        "Read length / kb": lengths / 1000.0,
        "Yield above length / Gbases": cumsum / 1e9,
    })
    df = df[df["Read length / kb"] <= xmax_kb].copy()
    if df.empty:
        df = pd.DataFrame({
            "Read length / kb": [0.0],
            "Yield above length / Gbases": [0.0],
        })
    elif len(df) > YIELD_THINNING:
        step = max(len(df) // YIELD_THINNING, 1)
        df = pd.concat((df.iloc[::step, :], df.iloc[[-1]]), axis=0).drop_duplicates()

    plot = ezc.lineplot(
        data=df,
        x="Read length / kb",
        y="Yield above length / Gbases",
        hue=None,
        color=color,
        marker=False,
    )
    plot._fig.xaxis.axis_label = "Read length / kb"
    plot._fig.yaxis.axis_label = "Yield above length / Gbases"
    plot._fig.x_range.start = 0.0
    plot._fig.x_range.end = xmax_kb
    plot._fig.add_layout(
        Title(
            text=f"Total yield: {_format_sigfig(total_yield / 1e9)} Gb.  N50: {_format_sigfig(compute_n50_kb(hist))}kb",
            text_font_size=ReportFontSizes.SUBTITLE,
            text_font=ReportFonts.PLOT,
        ),
        "above",
    )
    return plot


class SeqSummary(Snippet):
    """Generate run-statistics plots from fastcat histogram directories."""

    def __init__(self, seq_summary=None, sample_names=None, color=None, height=SUMMARY_PLOT_HEIGHT):
        super().__init__(styles=None, classes=None)
        self.color = color
        self.height = height

        if seq_summary is None:
            raise ValueError("`seq_summary` must be provided.")

        with self:
            if isinstance(seq_summary, tuple):
                tabs = Tabs()
                for sample_name, data in zip(sample_names, seq_summary):
                    with tabs.add_tab(sample_name):
                        self._draw_summary_plots(data)
            else:
                self._draw_summary_plots(seq_summary)

    def _draw_summary_plots(self, data):
        try:
            ldata, qdata = load_summary_histograms(data)
        except Exception:
            with Grid(columns=3):
                EZChart(empty_plot(), height=self.height)
                EZChart(empty_plot(), height=self.height)
                EZChart(empty_plot(), height=self.height)
            return

        with Grid(columns=3):
            EZChart(read_quality_plot(qdata, color=self.color), height=self.height)
            EZChart(read_length_plot(ldata, color=self.color), height=self.height)
            EZChart(base_yield_plot(ldata, color=self.color), height=self.height)
