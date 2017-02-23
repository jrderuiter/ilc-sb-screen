import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns
import toolz

from ..seaborn import clustermap as _clustermap
from ..util import color_annotation


def plot_nmf_coefficients(nmf_coef,
                          sort=True,
                          names=None,
                          order=None,
                          **kwargs):
    default_kws = dict(cmap='YlOrRd')

    # Assign names if given.
    if names is not None:
        nmf_coef = nmf_coef.copy()
        nmf_coef.columns = names

    # Reorder if needed.
    if order is not None:
        nmf_coef = nmf_coef[order]

    if sort:
        nmf_coef = nmf_sort_by_cluster(nmf_coef)

    # Create cluster map.
    g = _clustermap(
        nmf_coef.T,
        row_cluster=False,
        col_cluster=not sort,
        **toolz.merge(default_kws, kwargs))

    # Remove unneeded labels.
    g.ax_heatmap.set_xticklabels([])

    # Correct rotation.
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    # g.ax_heatmap.set_yticklabels([])

    # For some reason we need to correct rotation?
    if g.ax_col_colors is not None:
        ticks = g.ax_col_colors.get_yticklabels()
        g.ax_col_colors.set_yticklabels(ticks, rotation=0)

    # Move colorbar.
    # g.cax.set_position([.85, .2, .01, .45])

    return g


def _reindex_colors(colors, index):
    if colors is None or not isinstance(colors, pd.DataFrame):
        return colors
    else:
        missing = set(index) - set(colors.index)
        if len(missing) > 0:
            raise ValueError('Missing rows {}'.format(missing))
        return colors.ix[index]


def nmf_sort_by_cluster(nmf_coef):
    # Assign to clusters.
    clust = nmf_assign_clusters(nmf_coef)

    # Sort using assignments.
    clust['sort_value'] = clust['cluster'].map(
        dict(zip(nmf_coef.columns, range(nmf_coef.shape[0]))))
    clust.sort_values(
        ['sort_value', 'max_coef'], ascending=[True, False], inplace=True)

    # Apply sorting to coefficients.
    nmf_coef = nmf_coef.ix[clust.index]

    return nmf_coef


def nmf_assign_clusters(nmf_coef):
    return pd.DataFrame({'cluster': nmf_coef.idxmax(axis=1),
                         'max_coef': nmf_coef.max(axis=1)})


def plot_gene_boxplot(data, gene, cluster, ax=None, **kwargs):
    """Plots expression of given gene as boxplot, stratified by clusters."""

    # Draw boxplot.
    ax = sns.boxplot(data=data, x=cluster, y=gene, ax=ax, **kwargs)

    # Style lines.
    plt.setp(ax.lines, color='0.2')
    plt.setp(ax.artists, edgecolor='0.2')
    sns.despine(ax=ax)

    # Style x-labels.
    plt.setp(ax.get_xticklabels(), rotation=25)
    ax.set_xlabel('')

    return ax


def plot_gene_boxplots(data, genes, cluster, axes=None, **kwargs):
    """Plots expression of multiple genes as boxplot, stratified by clusters."""

    if axes is None:
        _, axes = plt.subplots(ncols=len(genes))
        axes = [axes] if len(genes) == 1 else axes.flatten()

    for gene, ax in zip(genes, axes):
        plot_gene_boxplot(data, gene, cluster, ax=ax, **kwargs)
        ax.set_title(gene, fontstyle='italic')

    # Style y-labels.
    axes[0].set_ylabel('Expression')
    for ax in axes[1:]:
        ax.set_ylabel('')

    return axes


def plot_gene_heatmap(data,
                      genes,
                      cluster,
                      order=None,
                      palette=None,
                      col_cluster=True,
                      **kwargs):
    palette = palette or sns.color_palette()
    order = order or list(data[cluster].unique())

    clusters = data[cluster].astype('category', categories=order, ordered=True)
    clusters.name = cluster.capitalize()

    cluster_colors = color_annotation(clusters.to_frame(), [palette])[0]

    if not col_cluster:
        sample_order = clusters.sort_values().index
        data = data.ix[sample_order]

    cm = _clustermap(
        data[genes].T,
        z_score=0,
        row_cluster=False,
        col_colors=cluster_colors,
        col_cluster=col_cluster,
        **kwargs)

    if not col_cluster:
        breaks = list(clusters.value_counts().ix[order].cumsum())
        breaks = breaks[:-1]

        for break_ in breaks:
            cm.ax_heatmap.axvline(break_, color='black', lw=0.8)
            cm.ax_col_colors.axvline(break_, color='black', lw=0.8)

    cm.ax_heatmap.set_xticklabels([])
    cm.ax_heatmap.set_ylabel('')
    plt.setp(cm.ax_heatmap.get_yticklabels(), rotation=0, fontstyle='italic')

    return cm

    # col_cluster=False, colorbar=False,
#, col_ratios={'side_colors': 0.12})
