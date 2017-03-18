import math

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import toolz

from nbsupport.morphology import MORPHOLOGY_ORDER, MORPHOLOGY_COLORS
from nbsupport.util.clustermap import color_annotation
from nbsupport.util.seaborn import clustermap

from .util import tidy_expression

SUBTYPE_ORDER = ['ILC-1', 'ILC-2', 'Spindle cell-like', 'Squamous-like']
SUBTYPE_COLORS = [sns.color_palette()[i] for i in list(range(3)) + [4]]


def plot_nmf_subtypes(nmf_coef, morphology, **kwargs):
    morphology_colors = color_annotation(
        morphology[list(MORPHOLOGY_ORDER)], colors=MORPHOLOGY_COLORS)[0]

    g = plot_nmf_coefficientmap(
        nmf_coef,
        sort=True,
        order=SUBTYPE_ORDER,
        col_colors=morphology_colors,
        col_ratios=dict(side_colors=0.5),
        row_colors=SUBTYPE_COLORS,
        row_ratios=dict(side_colors=0.011),
        cmap='Purples',
        colorbar=False,
        **kwargs)

    g.ax_heatmap.set_yticklabels([])

    g.ax_row_colors.set_yticks(np.arange(nmf_coef.shape[1], 0, -1) - 0.5)
    g.ax_row_colors.set_yticklabels(nmf_coef.columns)

    g.ax_row_colors.set_ylabel('NMF clusters')
    g.ax_col_colors.set_title('NMF subtypes')

    return g


def plot_nmf_coefficientmap(nmf_coef,
                            sort=True,
                            names=None,
                            order=None,
                            **kwargs):
    default_kws = dict(cmap='Purples')

    # Assign names if given.
    if names is not None:
        nmf_coef = nmf_coef.copy()
        nmf_coef.columns = names

    # Reorder if needed.
    if order is not None:
        nmf_coef = nmf_coef[order]

    if sort:
        nmf_coef = _sort_by_subtype(nmf_coef)

    # Create cluster map.
    g = clustermap(
        nmf_coef.T,
        row_cluster=False,
        col_cluster=not sort,
        **toolz.merge(default_kws, kwargs))

    # Remove unneeded labels.
    g.ax_heatmap.set_xticklabels([])

    # Correct rotation.
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    # For some reason we need to correct rotation?
    if g.ax_col_colors is not None:
        ticks = g.ax_col_colors.get_yticklabels()
        g.ax_col_colors.set_yticklabels(ticks, rotation=0)

    # Move colorbar.
    # g.cax.set_position([.85, .2, .01, .45])

    return g


def _sort_by_subtype(nmf_coef):
    # Assign to clusters.
    clust = assign_nmf_subtypes(nmf_coef)

    # Sort using assignments.
    clust['sort_value'] = clust['subtype'].map(
        dict(zip(nmf_coef.columns, range(nmf_coef.shape[0]))))
    clust.sort_values(
        ['sort_value', 'max_coef'], ascending=[True, False], inplace=True)

    # Apply sorting to coefficients.
    nmf_coef = nmf_coef.ix[clust.index]

    return nmf_coef


def assign_nmf_subtypes(nmf_coef):
    return pd.DataFrame({
        'subtype': nmf_coef.idxmax(axis=1),
        'max_coef': nmf_coef.max(axis=1)
    })


def plot_boxplot(expr, design, gene, ax=None, **kwargs):
    if ax is None:
        _, ax = plt.subplots()

    data = tidy_expression(expr.ix[[gene]], design)

    sns.boxplot(
        data=data,
        x='subtype',
        y='value',
        ax=ax,
        order=SUBTYPE_ORDER,
        palette=SUBTYPE_COLORS,
        **kwargs)

    ax.set_xlabel('Subtype')
    ax.set_ylabel('Expression')
    ax.set_title(gene, fontstyle='italic')

    # Style xtick labels.
    plt.setp(ax.get_xticklabels(), rotation=25)
    ax.set_xlabel('')

    sns.despine(ax=ax)

    return ax


def plot_boxplots(expr, design, genes, ncols=None, figsize=None, **kwargs):
    if ncols is None:
        ncols, nrows = len(genes), 1
    else:
        nrows = math.ceil(len(genes) / ncols)

    fig, axes = plt.subplots(
        ncols=ncols, nrows=nrows, sharex=True, figsize=figsize)

    if nrows == 1:
        axes = np.array([axes])

    for gene, ax in zip(genes, axes.flatten()):
        plot_boxplot(expr, design, gene, ax=ax, **kwargs)

    for row in axes:
        for ax in row[1:]:
            ax.set_ylabel('')

    plt.tight_layout()

    return fig, axes


def plot_heatmap(expr,
                 design,
                 genes,
                 order=None,
                 palette=None,
                 col_cluster=True,
                 **kwargs):

    data = pd.concat([expr.ix[genes].T, design[['subtype']]], axis=1)

    palette = palette or sns.color_palette()
    order = order or list(data['subtype'].unique())

    clusters = data['subtype'].astype(
        'category', categories=order, ordered=True)
    clusters.name = 'Subtype'

    cluster_colors = color_annotation(clusters.to_frame(), [palette])[0]

    if not col_cluster:
        sample_order = clusters.sort_values().index
        data = data.ix[sample_order]

    cm = clustermap(
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
