"""Shared Bokeh plot styling used by fastcat workflow report helpers."""


class ReportFontSizes:
    TITLE = "13pt"
    SUBTITLE = "9pt"
    AXIS_LABEL = "12pt"
    TICK_LABEL = "10pt"
    COLORBAR_TITLE = "12pt"
    COLORBAR_TICK = "10pt"


class ReportFonts:
    PLOT = "Noto Sans"


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

    bokeh_plot_cls.__init__ = patched_init
    bokeh_plot_cls._fastcat_axis_fonts_patched = True


def apply_colorbar_fonts(color_bar):
    """Use readable colorbar title and tick fonts."""
    color_bar.title_text_font_size = ReportFontSizes.COLORBAR_TITLE
    color_bar.title_text_font = ReportFonts.PLOT
    color_bar.major_label_text_font_size = ReportFontSizes.COLORBAR_TICK
    color_bar.major_label_text_font = ReportFonts.PLOT
