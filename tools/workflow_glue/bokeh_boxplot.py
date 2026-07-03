"""Bokeh-backed box plots for workflow reports."""

from bokeh.models import ColumnDataSource, HoverTool
from ezcharts.plots import BokehPlot
import numpy as np
import pandas as pd


def _unsupported(**kwargs):
    for name, value in kwargs.items():
        if value is not None:
            raise NotImplementedError(
                f"The parameter '{name}' with value {value} is not yet implemented.")


def boxplot(
        data=None, *, x=None, y=None, hue=None, order=None, hue_order=None,
        orient=None, color=None, palette=None, saturation=0.65, fill=True,
        dodge=None, width=None, gap=None, linecolor='auto', linewidth=1,
        hue_norm=None, native_scale=None, log_scale=None,
        formatter=None, legend=None, ax=None, lower=None, q1=None, q2=None,
        q3=None, upper=None, cap_width=0.3, nonnegative_y=True, hover_items=None,
        **kwargs):
    """Draw a Bokeh box plot from precomputed quantiles.

    The API follows ezcharts/seaborn `boxplot` where possible, but expects
    `lower`, `q1`, `q2`, `q3`, and `upper` to name existing columns.
    """
    _unsupported(
        hue=hue, order=order, hue_order=hue_order, orient=orient, dodge=dodge, gap=gap,
        hue_norm=hue_norm, native_scale=native_scale, log_scale=log_scale,
        formatter=formatter, legend=legend, ax=ax)

    if kwargs:
        raise NotImplementedError(
            f"Unsupported parameters: {', '.join(sorted(kwargs))}")

    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame with precomputed quantile columns")
    if x is None:
        x = data.columns[0]
    if y is None:
        y = q2
    if not all(col is not None for col in (lower, q1, q2, q3, upper)):
        raise ValueError("lower, q1, q2, q3, and upper must all be provided")

    required = [x, lower, q1, q2, q3, upper]
    missing = [col for col in required if col not in data]
    if missing:
        raise ValueError(f"missing required columns: {', '.join(missing)}")
    stats = data.copy()
    if y is not None and y not in stats:
        stats[y] = stats[q2]
    if "mean" not in stats:
        stats["mean"] = stats[q2]

    stats[x] = stats[x].astype(str)
    stats = stats.replace([np.inf, -np.inf], np.nan).dropna(
        subset=[lower, q1, q2, q3, upper])
    if stats.empty:
        return None

    groups = list(stats[x].unique())
    if width is None:
        width = 0.7

    if color is None:
        palette = list(BokehPlot.colors if palette is None else palette)
        if not palette:
            raise ValueError("palette must contain at least one color")
        color = palette[0]
    if not fill:
        color = None
        linecolor = "black"
    elif linecolor == 'auto':
        linecolor = "black"

    cap_half_width = cap_width * 0.5
    stats["cap_left"] = [(value, -cap_half_width) for value in stats[x]]
    stats["cap_right"] = [(value, cap_half_width) for value in stats[x]]

    y_min = float(stats[lower].min())
    y_max = float(stats[upper].max())
    y_padding = (y_max - y_min) * 0.1
    if y_padding == 0:
        y_padding = abs(y_max) * 0.1 or 1.0
    y_start = y_min - y_padding
    if nonnegative_y and y_min >= 0:
        y_start = max(0.0, y_start)
    y_range = (y_start, y_max + y_padding)

    source = ColumnDataSource(stats)
    plot = BokehPlot(x_range=groups, y_range=y_range)
    fig = plot._fig

    fig.segment(
        x0=x, y0=lower, x1=x, y1=q1, source=source,
        line_width=linewidth, line_color=linecolor)
    fig.segment(
        x0=x, y0=q3, x1=x, y1=upper, source=source,
        line_width=linewidth, line_color=linecolor)
    fig.segment(
        x0="cap_left", y0=lower, x1="cap_right", y1=lower, source=source,
        line_width=linewidth, line_color=linecolor)
    fig.segment(
        x0="cap_left", y0=upper, x1="cap_right", y1=upper, source=source,
        line_width=linewidth, line_color=linecolor)

    fig.vbar(
        x=x, width=width, top=q3, bottom=q2, source=source,
        fill_color=color, line_color=linecolor,
        line_width=linewidth, fill_alpha=saturation)
    fig.vbar(
        x=x, width=width, top=q2, bottom=q1, source=source,
        fill_color=color, line_color=linecolor,
        line_width=linewidth, fill_alpha=saturation)

    fig.xgrid.grid_line_color = None
    fig.xaxis.axis_label = str(x).capitalize()
    fig.yaxis.axis_label = str(y).capitalize() if y is not None else ""
    fig.xaxis.major_label_orientation = 1.0

    hover = fig.select(dict(type=HoverTool))
    hover.tooltips = hover_items or [
        ("Q1", f"@{{{q1}}}"),
        ("Median", f"@{{{q2}}}"),
        ("Q3", f"@{{{q3}}}"),
        ("Mean", "@mean"),
        ("Min", f"@{{{lower}}}"),
        ("Max", f"@{{{upper}}}"),
    ]

    return plot
