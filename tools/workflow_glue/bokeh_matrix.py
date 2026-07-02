"""Bokeh-backed matrix plots for workflow reports."""

from bokeh.models import BasicTicker, ColorBar, ColumnDataSource, HoverTool
from bokeh.models import LinearColorMapper, PrintfTickFormatter, Title
from bokeh.palettes import Viridis256
import numpy as np
import pandas as pd

from ezcharts.plots import BokehPlot


def heatmap(
        data, *, vmin=None, vmax=None, cmap=None, center=None, robust=False,
        annot=True, fmt='.2g', annot_kws=None, linewidths=0, linecolor='white',
        cbar=True, cbar_kws=None, cbar_ax=None, square=False, xticklabels='auto',
        yticklabels='auto', mask=None, ax=None, vmax_quantile=None,
        title=None, subtitle=None, **kwargs):
    """Plot rectangular data as a Bokeh heatmap.

    ezcharts heatmap hasn't been ported to Bokeh yet, so this is an implementation
    to get us going. it should find its way back into ezcharts eventually.
    """
    def _step(values):
        # compute the minimum step size between unique values in a sequence
        vals = np.unique(np.asarray(values, dtype=float))
        diffs = np.ediff1d(vals)
        diffs = diffs[diffs > 0.0]
        if diffs.size:
            return float(diffs.min())
        return 1.0

    values = pd.DataFrame(data).sort_index().sort_index(axis=1)

    image = values.to_numpy(dtype=float).T
    image = np.array(image, copy=True)
    if mask is not None:
        mask_values = np.asarray(mask, dtype=bool)
        if mask_values.shape != values.shape:
            raise ValueError("mask must have the same shape as data")
        image[mask_values.T] = np.nan
    image[~np.isfinite(image)] = np.nan

    x_name = values.index.name or "x"
    y_name = values.columns.name or "y"
    x_step = _step(values.index)
    y_step = _step(values.columns)

    # compute the image bounds, with a half-step padding on each side
    x0 = float(values.index.min()) - x_step * 0.5 if len(values.index) else 0.0
    x1 = float(values.index.max()) + x_step * 0.5 if len(values.index) else 1.0
    y0 = float(values.columns.min()) - y_step * 0.5 if len(values.columns) else 0.0
    y1 = float(values.columns.max()) + y_step * 0.5 if len(values.columns) else 1.0

    image_source = ColumnDataSource(data={
        "image": [image],
        "x": [x0],
        "y": [y0],
        "dw": [x1 - x0],
        "dh": [y1 - y0],
    })

    # figure out the color mapping range, using quantiles if requested
    finite = image[np.isfinite(image)]
    if finite.size:
        if vmax is None and vmax_quantile is not None:
            vmax = float(np.quantile(finite, vmax_quantile, method="nearest"))
        low = float(np.min(finite)) if vmin is None else vmin
        high = float(np.max(finite)) if vmax is None else vmax
    else:
        low, high = 0.0, 1.0

    if not np.isfinite(low):
        low = 0.0
    if not np.isfinite(high) or high <= low:
        high = low + 1.0

    palette = cmap if cmap is not None else Viridis256
    mapper = LinearColorMapper(palette=palette, low=low, high=high)

    # get plotted
    plt = BokehPlot(title="", output_backend="canvas")
    image_renderer = plt._fig.image(
        image="image", x="x", y="y", dw="dw", dh="dh",
        source=image_source,
        color_mapper=mapper)

    plt._fig.xaxis.axis_label = x_name
    plt._fig.yaxis.axis_label = y_name
    plt._fig.grid.grid_line_color = None
    plt._fig.x_range.start = x0
    plt._fig.x_range.end = x1
    plt._fig.y_range.start = y0
    plt._fig.y_range.end = y1

    if xticklabels is False:
        plt._fig.xaxis.major_label_text_font_size = "0pt"
    if yticklabels is False:
        plt._fig.yaxis.major_label_text_font_size = "0pt"

    # add tools and colorbar
    hover = plt._fig.select(dict(type=HoverTool))
    hover.renderers = [image_renderer]
    hover.tooltips = [
        (x_name, "$x{0.0}"),
        (y_name, "$y{0.000}"),
    ]

    if cbar:
        colorbar_options = dict(cbar_kws or {})
        colorbar_title = colorbar_options.pop(
            "title", colorbar_options.pop("label", None))
        colorbar_ticker = colorbar_options.pop("ticker", BasicTicker())
        colorbar_format = colorbar_options.pop("format", "%.2g")
        colorbar_formatter = colorbar_options.pop(
            "formatter", PrintfTickFormatter(format=colorbar_format))
        color_bar = ColorBar(
            color_mapper=mapper,
            ticker=colorbar_ticker,
            formatter=colorbar_formatter,
            title=colorbar_title)
        for attr, value in colorbar_options.items():
            setattr(color_bar, attr, value)
        plt._fig.add_layout(color_bar, "right")

    # add decorations
    if subtitle:
        plt._fig.add_layout(Title(text=subtitle, text_font_size="0.8em"), 'above')
    if title:
        plt._fig.add_layout(Title(text=title, text_font_size="1.5em"), 'above')

    return plt
