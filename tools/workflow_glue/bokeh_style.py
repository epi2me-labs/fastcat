"""Shared Bokeh plot styling used by fastcat workflow report helpers."""

from bokeh.models import CustomJSTickFormatter, LinearAxis


class ReportFontSizes:
    TITLE = "13pt"
    SUBTITLE = "9pt"
    AXIS_LABEL = "12pt"
    TICK_LABEL = "10pt"
    COLORBAR_TITLE = "12pt"
    COLORBAR_TICK = "10pt"


class ReportFonts:
    PLOT = "Noto Sans"


ADAPTIVE_TICK_FORMATTER_JS = r"""
const superscripts = {
  "-": "⁻",
  "0": "⁰",
  "1": "¹",
  "2": "²",
  "3": "³",
  "4": "⁴",
  "5": "⁵",
  "6": "⁶",
  "7": "⁷",
  "8": "⁸",
  "9": "⁹",
};

if (!Number.isFinite(tick)) {
  return "";
}

const values = ticks
  .map(Number)
  .filter((value) => Number.isFinite(value));
if (values.length === 0) {
  return "";
}

function shared_precision(numbers) {
  for (let precision = 0; precision <= 6; ++precision) {
    const scale = Math.pow(10, precision);
    const representable = numbers.every((number) => {
      const scaled = number * scale;
      return Math.abs(scaled - Math.round(scaled))
        <= 1e-8 * Math.max(1, Math.abs(scaled));
    });
    if (representable) {
      return precision;
    }
  }
  return 6;
}

const value = Number(tick);
const max_abs_tick = values.reduce(
  (maximum, number) => Math.max(maximum, Math.abs(number)), 0);
if (max_abs_tick < 1e4) {
  return value.toFixed(shared_precision(values));
}

const exponent = Math.floor(Math.log10(max_abs_tick));
const scale = Math.pow(10, exponent);
const precision = shared_precision(values.map((number) => number / scale));
const exponent_text = exponent
  .toString()
  .split("")
  .map((character) => superscripts[character] || character)
  .join("");
return `${(value / scale).toFixed(precision)} × 10 ${exponent_text}`;
"""


def adaptive_tick_formatter():
    """Create the shared formatter for numeric report scales."""
    return CustomJSTickFormatter(code=ADAPTIVE_TICK_FORMATTER_JS)


def patch_bokeh_plot_class(bokeh_plot_cls):
    """Monkey patch ezcharts BokehPlot once to increase axis typography."""
    if getattr(bokeh_plot_cls, "_fastcat_axis_fonts_patched", False):
        return

    original_init = bokeh_plot_cls.__init__

    def patched_init(self, *args, **kwargs):
        original_init(self, *args, **kwargs)
        for axis in list(self._fig.xaxis) + list(self._fig.yaxis):
            axis.axis_label_text_font_size = ReportFontSizes.AXIS_LABEL
            axis.axis_label_text_font = ReportFonts.PLOT
            axis.major_label_text_font_size = ReportFontSizes.TICK_LABEL
            axis.major_label_text_font = ReportFonts.PLOT
            if isinstance(axis, LinearAxis):
                axis.formatter = adaptive_tick_formatter()

    bokeh_plot_cls.__init__ = patched_init
    bokeh_plot_cls._fastcat_axis_fonts_patched = True


def apply_colorbar_fonts(color_bar):
    """Use readable colorbar title and tick fonts."""
    color_bar.title_text_font_size = ReportFontSizes.COLORBAR_TITLE
    color_bar.title_text_font = ReportFonts.PLOT
    color_bar.major_label_text_font_size = ReportFontSizes.COLORBAR_TICK
    color_bar.major_label_text_font = ReportFonts.PLOT
    color_bar.formatter = adaptive_tick_formatter()
