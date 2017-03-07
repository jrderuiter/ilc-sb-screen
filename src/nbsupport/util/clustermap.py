from collections import OrderedDict

import numpy as np
import pandas as pd

from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize, rgb2hex


def color_annotation(df, colors, vlims=None, bg_color='#ffffff'):
    """Converts a data frame of annotations to colors."""

    if vlims is None:
        vlims = [None] * len(colors)

    # TODO: Fix boolean/nan case.

    colored_cols = OrderedDict()
    color_maps = OrderedDict()

    for (col_name, values), color, vlim in zip(df.items(), colors, vlims):
        colored, color_map = _color_column(values, color, bg_color, vlim=vlim)
        colored_cols[col_name] = colored
        color_maps[col_name] = color_map

    return pd.DataFrame(colored_cols), color_maps


def _color_column(series, color, bg_color, vlim=None):
    """Converts an annotation column to colors."""

    if series.dtype == bool:
        # Boolean case.
        return _color_bool(series, color, bg_color)
    elif series.dtype.kind in 'biufc':
        # Numeric case.
        vlim = vlim or (None, None)
        return _color_numeric(
            series, color, bg_color, vmin=vlim[0], vmax=vlim[1])
    elif series.dtype.kind == 'O':
        # Strings and categorical case.
        return _color_categorical(series, color, bg_color)
    else:
        raise ValueError('Unsupported dtype: {}'.format(series.dtype))


def _color_bool(series, color, bg_color):
    """Converts a boolean annotation column to colors."""

    mapped = series.map({True: color, False: bg_color}, na_action='ignore')
    return mapped, color


def _color_numeric(series, color, bg_color, n=200, vmax=None, vmin=None):
    """Converts a numeric annotation column to colors."""

    vmax = vmax or series.max()
    vmin = vmin or series.min()

    cmap = LinearSegmentedColormap.from_list(
        name='custom_cmap', colors=[bg_color, color], N=n)
    norm = Normalize(vmin=vmin, vmax=vmax)
    mappable = ScalarMappable(norm=norm, cmap=cmap)

    rgba_colors = mappable.to_rgba(series.values)
    hex_colors = (rgb2hex(rgba) for rgba in rgba_colors)

    mapped = pd.Series(hex_colors, index=series.index, name=series.name)

    return mapped, color


def _color_categorical(series, colors, bg_color):
    """Converts a categorical annotation column to colors."""

    if series.dtype.name == 'category':
        str_values = series.cat.categories
    else:
        str_values = set(series.dropna())

    if isinstance(colors, list):
        colors = [_rgb_to_hex(c, normalized=True) for c in colors]
        color_map = OrderedDict(zip(str_values, colors))
        mapped = series.map(color_map).fillna(bg_color)
    else:
        numeric_values = np.linspace(0, 1, len(str_values) + 1)[1:]
        numeric_map = OrderedDict(zip(str_values, numeric_values))

        mapped = _color_numeric(
            series.map(numeric_map), colors, bg_color=bg_color)[0]
        color_map = dict(zip(series.values, mapped.values))

    return mapped, color_map


def _rgb_to_hex(rgb, normalized=True):
    if normalized:
        rgb = tuple(map(lambda x: int(x * 255), rgb))
    return '#%02x%02x%02x' % rgb


def sort_matrix(df, sort_columns=True):
    """Sorts a DataFrame, first by its columns and then by its rows."""

    if sort_columns:
        freqs = (df.astype(float) > 0).sum(axis=0)
        order = list(freqs.sort_values(ascending=False).index)
    else:
        order = list(df.columns)

    df_sorted = df[order]
    df_sorted = df_sorted.sort_values(by=order, ascending=False)

    return df_sorted
