"""Helpers for an in-figure null-data state on Bokeh plots."""

from bokeh.models import Label
from ezcharts.plots import BokehPlot
from ezcharts.plots.util import ReportFontSizes, ReportFonts


NULL_LABEL_NAME = "fastcat_null_label"
NULL_XAXIS_NAME = "fastcat_null_xaxis"
NULL_YAXIS_NAME = "fastcat_null_yaxis"
NULL_XGRID_NAME = "fastcat_null_xgrid"
NULL_YGRID_NAME = "fastcat_null_ygrid"
NULL_COLORBAR_NAME = "fastcat_null_colorbar"


def attach_null_state(plot):
    """Decorate a plot with a hidden centered null-data label."""
    fig = plot._fig
    if getattr(fig, "_fastcat_null_state_attached", False):
        return plot

    fig._fastcat_real_toolbar_location = fig.toolbar_location

    if fig.xaxis:
        fig.xaxis[0].name = NULL_XAXIS_NAME
    if fig.yaxis:
        fig.yaxis[0].name = NULL_YAXIS_NAME
    if fig.xgrid:
        fig.xgrid[0].name = NULL_XGRID_NAME
    if fig.ygrid:
        fig.ygrid[0].name = NULL_YGRID_NAME
    for item in fig.right:
        if item.__class__.__name__ == "ColorBar":
            item.name = NULL_COLORBAR_NAME

    label = Label(
        x=0.5,
        y=0.5,
        x_units="data",
        y_units="data",
        text="",
        text_align="center",
        text_baseline="middle",
        text_font=ReportFonts.PLOT,
        text_font_size=ReportFontSizes.AXIS_LABEL,
        visible=False,
        name=NULL_LABEL_NAME)
    fig.add_layout(label)
    fig._fastcat_null_label = label
    fig._fastcat_null_state_attached = True
    return plot


def set_null_state(plot, message):
    """Show a centered null-data message and hide plot chrome."""
    attach_null_state(plot)
    fig = plot._fig
    fig._fastcat_null_label.text = message
    fig._fastcat_null_label.visible = True
    for axis in fig.xaxis:
        axis.visible = False
    for axis in fig.yaxis:
        axis.visible = False
    for grid in fig.xgrid:
        grid.visible = False
    for grid in fig.ygrid:
        grid.visible = False
    for item in fig.right:
        if getattr(item, "name", None) == NULL_COLORBAR_NAME:
            item.visible = False
    fig.toolbar_location = None
    fig.x_range.start = 0.0
    fig.x_range.end = 1.0
    fig.y_range.start = 0.0
    fig.y_range.end = 1.0
    return plot


def clear_null_state(plot):
    """Restore normal plot chrome after a null-data state."""
    attach_null_state(plot)
    fig = plot._fig
    fig._fastcat_null_label.visible = False
    for axis in fig.xaxis:
        axis.visible = True
    for axis in fig.yaxis:
        axis.visible = True
    for grid in fig.xgrid:
        grid.visible = True
    for grid in fig.ygrid:
        grid.visible = True
    for item in fig.right:
        if getattr(item, "name", None) == NULL_COLORBAR_NAME:
            item.visible = True
    fig.toolbar_location = fig._fastcat_real_toolbar_location
    return plot


def null_plot(message, x_label="", y_label=""):
    """Create a placeholder figure showing an in-plot null-data message."""
    plot = BokehPlot(title="", output_backend="canvas")
    plot._fig.xaxis.axis_label = x_label
    plot._fig.yaxis.axis_label = y_label
    attach_null_state(plot)
    set_null_state(plot, message)
    return plot
