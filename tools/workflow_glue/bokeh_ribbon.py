"""Bokeh-backed quantile ribbon plots for workflow reports."""

from bokeh.models import ColumnDataSource, HoverTool, Title
from ezcharts.plots import BokehPlot
import numpy as np
import pandas as pd


def _unsupported(**kwargs):
    for name, value in kwargs.items():
        if value is not None:
            raise NotImplementedError(
                f"The parameter '{name}' with value {value} is not yet implemented.")


def ribbonplot(
        data=None, *, x=None, y=None, hue=None, order=None, hue_order=None,
        orient=None, color=None, palette=None, saturation=0.32, fill=True,
        dodge=None, width=None, gap=None, linecolor="#222222", linewidth=1.5,
        hue_norm=None, native_scale=None, log_scale=None,
        formatter=None, legend=None, ax=None, lower=None, q1=None, q2=None,
        q3=None, upper=None, outer_alpha=0.18, subtitle=None,
        nonnegative_y=True, hover_items=None, **kwargs):
    """Draw a Bokeh ribbon plot from precomputed quantiles."""
    _unsupported(
        hue=hue, order=order, hue_order=hue_order, orient=orient, dodge=dodge, gap=gap,
        hue_norm=hue_norm, native_scale=native_scale, log_scale=log_scale,
        formatter=formatter, legend=legend, ax=ax, width=width)

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
    stats = stats.replace([np.inf, -np.inf], np.nan).dropna(
        subset=[x, lower, q1, q2, q3, upper])
    if stats.empty:
        return None

    stats[x] = pd.to_numeric(stats[x], errors="coerce")
    stats = stats.dropna(subset=[x]).sort_values(x).reset_index(drop=True)
    if stats.empty:
        return None

    if y is not None and y not in stats:
        stats[y] = stats[q2]

    if color is None:
        if palette:
            color = palette[0]
        else:
            color = "#6aa6f8"

    x_min = float(stats[x].min())
    x_max = float(stats[x].max())
    x_padding = (x_max - x_min) * 0.02
    if x_padding == 0:
        x_padding = abs(x_max) * 0.02 or 1.0

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
    plot = BokehPlot(
        title="",
        x_range=(x_min - x_padding, x_max + x_padding),
        y_range=y_range,
        output_backend="canvas")
    fig = plot._fig

    if fill:
        fig.varea(
            x=x,
            y1=lower,
            y2=upper,
            source=source,
            fill_color=color,
            fill_alpha=outer_alpha)
        fig.varea(
            x=x,
            y1=q1,
            y2=q3,
            source=source,
            fill_color=color,
            fill_alpha=saturation)

    fig.line(
        x=x, y=lower, source=source,
        line_color=linecolor, line_alpha=0.6, line_width=1.0)
    fig.line(
        x=x, y=q1, source=source,
        line_color=color, line_alpha=0.9, line_width=1.2)
    fig.line(
        x=x, y=q2, source=source,
        line_color=linecolor, line_width=linewidth)
    fig.line(
        x=x, y=q3, source=source,
        line_color=color, line_alpha=0.9, line_width=1.2)
    fig.line(
        x=x, y=upper, source=source,
        line_color=linecolor, line_alpha=0.6, line_width=1.0)
    hover_renderer = fig.scatter(
        x=x, y=q2, source=source,
        size=10, line_alpha=0.0, fill_alpha=0.0)

    fig.xgrid.grid_line_alpha = 0.25
    fig.ygrid.grid_line_alpha = 0.25
    fig.xaxis.axis_label = str(x).capitalize()
    fig.yaxis.axis_label = str(y).capitalize() if y is not None else ""

    hover = fig.select(dict(type=HoverTool))
    hover.renderers = [hover_renderer]
    hover.tooltips = hover_items or [
        (str(x), f"@{{{x}}}{{0.00}}"),
        ("Q05", f"@{{{lower}}}{{0.00}}"),
        ("Q25", f"@{{{q1}}}{{0.00}}"),
        ("Q50", f"@{{{q2}}}{{0.00}}"),
        ("Q75", f"@{{{q3}}}{{0.00}}"),
        ("Q95", f"@{{{upper}}}{{0.00}}"),
    ]

    if subtitle:
        fig.add_layout(Title(text=subtitle, text_font_size="0.8em"), "above")

    return plot
