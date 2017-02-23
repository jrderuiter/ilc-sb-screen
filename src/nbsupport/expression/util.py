from ..clustermap import clustermap as _clustermap

from itertools import chain

import matplotlib.pyplot as plt
from matplotlib import patches as mpatches

from scipy.spatial import distance
from scipy.cluster import hierarchy


def calc_zscore(data2d, axis=1):
    """Standarize the mean and variance of the data axis."""
    other_axis = 0 if axis == 1 else 1
    return (data2d.subtract(
        data2d.mean(axis=other_axis), axis=axis).divide(
            data2d.std(axis=other_axis), axis=axis))


def calc_linkage(data2d,
                 axis=0,
                 z_score=None,
                 metric='euclidean',
                 method='complete'):
    if z_score is not None:
        data2d = calc_zscore(data2d, axis=z_score)

    if axis == 0:
        data2d = data2d.T

    dist = distance.pdist(data2d, metric=metric)
    linkage = hierarchy.linkage(dist, method=method)

    return linkage


def clustermap(*args, legend_cmaps=None, legend_kws=None, **kwargs):
    """Plots PAM50 genemap using sns.clustermap."""

    cm = _clustermap(*args, **kwargs)

    cm.ax_heatmap.set_xticks([])
    plt.setp(cm.ax_heatmap.get_yticklabels(), rotation=0, style='italic')

    cm.ax_heatmap.set_xlabel('Samples')
    cm.ax_heatmap.set_ylabel('')

    if legend_cmaps is not None:
        draw_legends(legend_cmaps, ax=cm.ax_heatmap, **(legend_kws or {}))

    return cm

# def clustermap(data,
#                *args,
#                col_annotation=None,
#                col_palette=None,
#                row_annotation=None,
#                row_palette=None,
#                legend_position=(1.2, 1),
#                legend_offset=0.2,
#                **kwargs):

#     # Color annotation if given.
#     color_maps = []

#     if col_annotation is not None:
#         kwargs['col_colors'], col_color_map = \
#             color_annotation(col_annotation, col_palette)
#         color_maps.append(col_color_map)

#     if row_annotation is not None:
#         kwargs['row_colors'], row_color_map = \
#             color_annotation(row_annotation, row_palette)
#         color_maps.append(row_color_map)

#     # Draw clustermap.
#     g = sns.clustermap(data, *args, **kwargs)

#     # Draw legends.
#     draw_legends(
#         color_maps,
#         ax=g.ax_heatmap,
#         position=legend_position,
#         offset=legend_offset)

#     return g


def draw_legends(color_maps,
                 ax,
                 position=(1.2, 1),
                 offset=0.2,
                 loc='upper right',
                 **kwargs):
    # Flatten list of dicts to items.
    color_maps = chain.from_iterable((cmap.items() for cmap in color_maps))

    # Draw legends.
    for i, ((name, color_map), not_last) in enumerate(lookahead(color_maps)):
        bbox_position = (position[0], position[1] - (offset * i))
        legend = draw_legend(
            color_map,
            ax=ax,
            name=name,
            loc=loc,
            bbox_to_anchor=bbox_position,
            **kwargs)

        if not_last:
            ax.add_artist(legend)


def draw_legend(color_map, ax, name=None, **kwargs):
    patches = [mpatches.Patch(
        color=color, label=label) for label, color in color_map.items()]
    legend = ax.legend(handles=patches, title=name, **kwargs)
    return legend


def lookahead(iterable):
    """Pass through all values from the given iterable, augmented by the
    information if there are more values to come after the current one
    (True), or if it is the last value (False).
    """

    # Get an iterator and pull the first value.
    it = iter(iterable)
    last = next(it)

    # Run the iterator to exhaustion (starting from the second value).
    for val in it:
        # Report the *previous* value (more to come).
        yield last, True
        last = val

    # Report the last value.
    yield last, False
