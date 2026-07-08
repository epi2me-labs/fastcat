"""Relational plots with optional stepped lines."""

import numbers

from bokeh.models import HoverTool, Title
from seaborn.relational import _LinePlotter

from ezcharts.plots import BokehPlot, util


DEFAULT_PALETTE = util.choose_palette()
DEFAULT_COLOR = DEFAULT_PALETTE[0]


class Mixin:
    """Mixin class to define how to plot lines."""

    def plot(self, plt, kws):
        """Plot the plot."""
        data = self.plot_data.dropna()
        if data.empty:
            return plt, []

        if "hue" not in data:
            data["hue"] = "values"

        x_name = self.variables.get("x", None)
        y_name = self.variables.get("y", None)

        symbol_size = kws.get("s", 8)
        if not isinstance(symbol_size, numbers.Number):
            raise NotImplementedError(
                "symbol size (s) currently only accepts a constant value")

        line_width = kws.get("linewidth", 3)
        if not isinstance(line_width, numbers.Number):
            raise NotImplementedError(
                "linewidth currently only accepts a constant value")

        marker = kws.get("marker", False)

        y_min = 0
        renderers = []
        for hue_name, df in data.groupby("hue"):
            relational_kwargs = {}
            if data["hue"].nunique() > 1:
                relational_kwargs["legend_label"] = hue_name

            y_min = min(y_min, df.y.min())
            color = kws.get("color")
            if self._hue_map.levels:
                color = self._hue_map(hue_name)

            if self.series_type == "line":
                renderer = plt._fig.line(
                    df.x, df.y, line_color=color,
                    line_width=line_width, **relational_kwargs)
            elif self.series_type == "steps-pre":
                renderer = plt._fig.step(
                    df.x, df.y, mode="before", line_color=color,
                    line_width=line_width, **relational_kwargs)
            elif self.series_type == "steps-mid":
                renderer = plt._fig.step(
                    df.x, df.y, mode="center", line_color=color,
                    line_width=line_width, **relational_kwargs)
            elif self.series_type == "steps-post":
                renderer = plt._fig.step(
                    df.x, df.y, mode="after", line_color=color,
                    line_width=line_width, **relational_kwargs)
            else:
                raise NotImplementedError(f"Unsupported series type {self.series_type}")
            renderers.append(renderer)

            if marker is not False:
                if marker is None or marker is True:
                    marker = "circle"
                plt._fig.scatter(
                    df.x, df.y, size=symbol_size, marker=marker,
                    color=color, **relational_kwargs)

            plt._fig.xaxis.axis_label = x_name
            plt._fig.yaxis.axis_label = y_name

        plt._fig.y_range.start = y_min
        return plt, renderers


class LinePlotter(Mixin, _LinePlotter):
    """Making line plots."""

    series_type = "line"


def lineplot(
        data=None, *,
        x=None, y=None, hue=None, size=None, style=None, units=None,
        palette=None, hue_order=None, hue_norm=None,
        sizes=None, size_order=None, size_norm=None,
        dashes=True, markers=False, style_order=None,
        estimator="mean", errorbar=("ci", 95), n_boot=1000, seed=None,
        orient="x", sort=True, err_style="band", err_kws=None,
        legend="auto", ci="deprecated", ax=None, bokeh_kwargs=None,
        hover_items=None, subtitle=None, drawstyle=None, **kwargs):
    """Draw a line plot with possibility of several semantic groupings."""
    if palette is None:
        palette = DEFAULT_PALETTE

    if not kwargs.get("color"):
        kwargs["color"] = DEFAULT_COLOR

    variables = _LinePlotter.get_semantics(locals())
    p = LinePlotter(
        data=data, variables=variables, estimator=estimator, n_boot=n_boot,
        seed=seed, errorbar=errorbar, sort=sort, orient=orient,
        err_style=err_style, err_kws=err_kws, legend=legend)
    p.map_hue(palette=palette, order=hue_order, norm=hue_norm)
    p.map_size(sizes=sizes, order=size_order, norm=size_norm)
    p.map_style(markers=markers, dashes=dashes, order=style_order)

    if drawstyle is None:
        p.series_type = "line"
    elif drawstyle in ("steps-pre", "steps-mid", "steps-post"):
        p.series_type = drawstyle
    else:
        raise NotImplementedError(f"Unsupported drawstyle: {drawstyle}")

    title = kwargs.pop("title", "")
    plt = BokehPlot(title=title, **(bokeh_kwargs or {}))
    hover = plt._fig.select(dict(type=HoverTool))
    hover.tooltips = hover_items or [(x, "@x"), (y, "@y")]

    plt, renderers = p.plot(plt, kwargs)
    if renderers:
        hover.renderers = renderers

    if subtitle:
        plt._fig.add_layout(Title(text=subtitle, text_font_size="0.8em"), "above")

    return plt
