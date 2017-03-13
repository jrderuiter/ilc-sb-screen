from collections import OrderedDict

import bisect

import matplotlib.pyplot as plt
import matplotlib.patches as mpl_patches
import matplotlib.collections as mpl_coll
import numpy as np
import pandas as pd
import seaborn as sns

import pybiomart
import toolz

from geneviz.tracks import BiomartTrack, FeatureTrack, RugTrack, plot_tracks

from .util.clustermap import color_annotation


def annotate_with_clonality(insertions):
    """Annotates an insertion dataframe with the clonality of insertions.

    The clonality of each insertion is calculated by dividing the support
    of the insertion by the support score of the strongest insertion in the
    same tumor.

    Parameters
    ----------
    insertions : pd.DataFrame
        Insertion dataframe.

    Returns
    -------
    pd.DataFrame
        Annotated insertion dataframe, containing an extra 'clonality' column.

    """

    grps = []
    for _, grp in insertions.groupby('sample'):
        grp = grp.copy()
        grp['clonality'] = grp['support'] / grp['support'].max()

        grps.append(grp)
    return pd.concat(grps, axis=0)


def plot_insertion_stats(insertions, fig_kws=None, suffix='', color=None):
    """Plots an overview of various statistics for the given insertions.

    Creates the following for plots for the given insertions:
        - Histogram of insertion depth
        - Histogram of insertion clonality
        - Histogram of the maximum insertion depth in each sample
        - Histogram of the number of insertions in each sample

    Parameters
    ----------
    insertions : pd.DataFrame
        Dataframe of insertions.
    fig_kws : dict[str, Any]
        Keyword arguments
    suffix : str
        Optional suffix to add to the subplot titles.
    color
        A valid matplotlib color, used to color the histograms.

    Returns
    -------
    Tuple[matplotlib.Figure, matplotlib.Axes]
        Returns a tuple containing the drawn figure and axes.

    """

    fig, axes = plt.subplots(ncols=2, nrows=2, **(fig_kws or {}))
    axes = axes.flatten()

    # Insertion depth.
    ins_depth = insertions.groupby('id')['depth_unique'].max()
    ins_depth.hist(
        bins=np.arange(0, np.max(ins_depth), 10),
        ax=axes[0],
        grid=False,
        color=color)
    axes[0].set_title('Insertion depth' + suffix)
    axes[0].set_xlabel('Insertion depth')
    axes[0].set_ylabel('Number of insertions')

    # Maximum depth per sample.
    ins_depth_sample = insertions.groupby('sample')['depth_unique'].max()
    ins_depth_sample.hist(
        bins=np.arange(0, 180, 10), ax=axes[2], grid=False, color=color)
    axes[2].set_title('Maximum insertion depth per sample' + suffix)
    axes[2].set_xlabel('Maximum insertion depth')
    axes[2].set_ylabel('Number of samples')

    # Insertion clonality
    ins_clon = insertions.groupby('id')['clonality'].max()
    ins_clon.hist(bins=20, ax=axes[1], grid=False, color=color)
    axes[1].set_title('Insertion clonality' + suffix)
    axes[1].set_xlabel('Insertion clonality')
    axes[1].set_ylabel('Number of insertions')

    # Number of insertion per sample.
    ins_sample = insertions.groupby('sample')['id'].nunique()
    ins_sample.hist(bins=15, ax=axes[3], grid=False, color=color)
    axes[3].axvline(ins_sample.median(), linestyle='dashed', color='red')
    axes[3].set_title('Number of insertions per sample' + suffix)
    axes[3].set_xlabel('Number of insertions')
    axes[3].set_ylabel('Number of samples')

    sns.despine(fig)
    fig.tight_layout()

    return fig, axes


def gene_statistics(insertions):
    """Calculates some basic statistics for genes in the given insertion set.

    Includes the following statistics:
        - n_samples: The number of samples in which a gene has an insertion
        - sense_fraction: The fraction of insertions that are integrated in
            a sense orientation with respect to their target gene.
        - sense_fraction_weighted: Same as sense_fraction, but insertions
            are weighted by their clonality when calculating the fraction.
        - mean_clonality: Mean clonality of insertions in the gene.

    Parameters
    ----------
    insertions : pd.DataFrame
        Insertion dataframe.

    Returns
    -------
    pd.DataFrame
        Dataframe containing four columns with the various statistics.

    """

    return pd.DataFrame(
        {
            'n_samples': gene_sample_count(insertions),
            'sense_fraction': gene_sense_fraction(insertions),
            'sense_fraction_weighted':
            gene_sense_fraction_weighted(insertions),
            'mean_clonality': gene_clonality(insertions)
        },
        columns=[
            'n_samples', 'mean_clonality', 'sense_fraction',
            'sense_fraction_weighted'
        ])


def gene_sample_count(insertions, name='n_samples'):
    """Calculates the number of samples with insertions for each gene.

    Parameters
    ----------
    insertions : pd.DataFrame
        Insertion dataframe.
    name : str
        Name to use for the returned series.

    Returns
    -------
    pd.Series
        Gene sample counts for the given insertion set.
    """

    count = insertions.groupby('gene_name')['sample'].nunique()
    count.name = name
    return count


def gene_sense_fraction(insertions):
    """Calculates the fraction of sense insertions per gene.

    Parameters
    ----------
    insertions : pd.DataFrame
        Insertion dataframe.

    Returns
    -------
    pd.Series
        Gene sense fractions for the given insertion set.
    """

    def _sense_frac(x):
        is_sense = x['gene_orientation'] == 'sense'
        return np.sum(is_sense) / len(x)

    return insertions.groupby(['gene_name']).apply(_sense_frac)


def gene_sense_fraction_weighted(insertions):
    """Calculates the weighted fraction of sense insertions per gene.

    Parameters
    ----------
    insertions : pd.DataFrame
        Insertion dataframe.

    Returns
    -------
    pd.Series
        Gene weighted sense fractions for the given insertion set.
    """

    def _sense_frac_weighted(x):
        is_sense = x['gene_orientation'] == 'sense'
        return np.sum(is_sense * x['clonality']) / np.sum(x['clonality'])

    return insertions.groupby(['gene_name']).apply(_sense_frac_weighted)


def gene_clonality(insertions):
    """Calculates the average clonality of insertions per gene.

    Parameters
    ----------
    insertions : pd.DataFrame
        Insertion dataframe.

    Returns
    -------
    pd.Series
        Gene clonalities for the given insertion set.
    """

    return clonality_matrix(insertions).mean(axis=1)


def clonality_matrix(insertions, value='clonality', fill_value=None):
    """Summarizes insertion clonality values in genes-by-samples matrix.

    Parameters
    ----------
    insertions : pd.DataFrame
        Insertion dataframe.
    value : str
        Name of the column to use for clonality values.
    fill_value : float
        Value to use for missing entries.

    Returns
    -------
    pd.DataFrame
        Matrix of genes-by-samples containing clonality values.
    """

    # Summarize gene clonality per sample.
    gene_clon = pd.pivot_table(
        insertions,
        index='gene_name',
        columns='sample',
        values=value,
        aggfunc='max',
        fill_value=fill_value)

    return gene_clon


def sort_matrix(mat):
    """Sorts a 2D matrix, first by its rows and then by its columns.

    Parameters
    ----------
    mat : pd.DataFrame
        Matrix to sort.

    Returns
    -------
    pd.DataFrame
        Sorted matrix.

    """

    freqs = (mat > 0).sum(axis=1)
    order = list(freqs.sort_values(ascending=False).index)

    mat_sorted = mat.ix[order]
    mat_sorted = mat_sorted.T.sort_values(by=order, ascending=False).T

    return mat_sorted


def plot_insertion_matrix(insertions,
                          genes=None,
                          value='clonality',
                          ax=None,
                          **kwargs):
    """Plots a 2D sample-by-gene heatmap of an insertion matrix.

    Parameters
    ----------
    insertions : pd.DataFrame
        Insertion dataframe.
    genes : List[str]
        List of genes to plot.
    value : str
        Name of column containing values to use in the heatmap.
    ax : matplotlib.Axis
        Axis to draw on.
    kwargs:
        Any keyword arguments are passed to seaborns heatmap function.

    Returns
    -------
    matplotlib.Axis
        Axis that was drawn on.

    """

    if ax is None:
        _, ax = plt.subplots()

    # Convert to matrix.
    ins_matrix = clonality_matrix(insertions, value=value)

    # Subset if needed.
    if genes is not None:
        ins_matrix = ins_matrix.loc[genes]
        ins_matrix = ins_matrix.dropna(how='all')

    # Sort matrix by columns.
    ins_matrix = sort_matrix(ins_matrix)

    # Draw heatmap.
    cmap = sns.light_palette('#376CA3', as_cmap=True)
    #cbar_ax = fig.add_axes([0.92, 0.35, 0.01, 0.3])
    sns.heatmap(
        ins_matrix,
        ax=ax,
        cmap=cmap,
        cbar_kws={'label': 'Clonality'},
        **kwargs)

    ax.set_xticks([])
    plt.setp(ax.get_yticklabels(), rotation=0, style='italic')

    ax.set_xlabel('Samples')
    ax.set_ylabel('Genes')

    for loc in ['top', 'right', 'left', 'bottom']:
        ax.spines[loc].set_visible(True)

    return ax


def plot_subtype_counts(insertions,
                        subtypes,
                        color='lightgrey',
                        highlight=None,
                        highlight_color='#c44e52',
                        highlight_labels=None,
                        legend_kws=None,
                        annotate_kws=None,
                        **kwargs):
    """Draws bar factorplot counting gene occurrences per subtype.

    Parameters
    ----------
    insertions : pd.DataFrame
        Insertion dataframe.
    subtypes : pd.Series
        Series containing the subtype assignment. Index should reflect
        samples in the given insertion set.
    color : str
        Color to use for the bars.
    highlight : Dict[str, List[str]]
        Dict of subtypes -> genes, indicating which genes should be highlighted
        per subtype (if any). Highlighted genes are given a different color.
    highlight_color : str
        Color to use for the bars of highlighted genes.
    legend_kws : Dict[str, Any]
        Dictionary of keyword arguments to pass to ax.legend.
    annotate_kws : Dict[str, Any]
        Dictionary of keyword arguments to pass to ax.annotate when
        labeling gene counts.
    **kwargs
        Dictionary of other keyword arguments to pass to seaborn.factorplot.

    Returns
    -------
    g : FacetGrid
        Returns the FacetGrid object with the plot on it for further tweaking.
    """

    legend_kws = toolz.merge(dict(loc='lower right', frameon=True),
                             legend_kws or {}) # yapf: disable

    # Summarize gene clonality per sample.
    subtype_map = dict(zip(subtypes.index, subtypes.values))

    merged = (clonality_matrix(insertions).unstack().dropna()
              .reset_index(name='clonality'))
    merged['subtype'] = merged['sample'].map(subtype_map)

    # Count occurrences / fraction.
    merged_count = (pd.crosstab(merged['gene_name'], merged['subtype'])
                    .unstack().reset_index(name='count'))

    subtype_size = subtypes.value_counts()
    merged_count['fraction'] = merged_count.apply(
        lambda r: r['count'] / subtype_size.ix[r['subtype']], axis=1)

    # Draw barplots.
    g = sns.factorplot(
        data=merged_count,
        x='fraction',
        y='gene_name',
        col='subtype',
        kind='bar',
        color=color,
        **kwargs)

    # Add counts.
    _annotate_with_counts(g, subtype_size, **(annotate_kws or {}))

    # Highlight entries if needed.
    if highlight is not None:
        _highlight_bars(g, highlight, highlight_color)

        if highlight_labels is None:
            highlight_labels = ('Highlighted genes', 'Other genes')

        cmap = OrderedDict(zip(highlight_labels, (highlight_color, color)))
        _draw_legend(ax=g.axes[0, 0], color_map=cmap, **(legend_kws or {}))

    # Adjust titles to remove prefix.
    for col_name, ax in zip(g.col_names, g.axes.flatten()):
        ax.set_title(col_name)

    # Clean up axes.
    for ax in g.axes.flatten():
        ax.set_xlim(0, 1)
        ax.set_xlabel('Fraction of samples')

    g.axes[0, 0].set_ylabel('Genes')

    # Make gene labels italic.
    plt.setp(g.axes[0, 0].get_yticklabels(), style='italic')

    return g


def _annotate_with_counts(g, col_total, **kwargs):
    """Annotates fractional horizontal bar factor plot with counts."""

    for clust, ax in zip(g.col_names, g.axes.flatten()):
        total = col_total.ix[clust]
        for patch in ax.patches:
            fraction = patch.get_width()

            if np.isnan(fraction):
                fraction = 0.0

            count = int(round(fraction * total))

            ax.text(
                x=fraction + 0.02,
                y=(patch.get_y() + (patch.get_height() / 2)),
                s='{}'.format(count),
                va='center',
                **kwargs)


def _highlight_bars(g, values, color):
    """Highlights specific entries of a bar factorplot."""

    # Determine row entry order.
    row_names = [t.get_text() for t in g.axes[0, 0].get_yticklabels()]

    # If values is not a dict, replicate for col_names.
    if not isinstance(values, dict):
        values = {col_name: values for col_name in g.col_names}

    # Color patches.
    for col_name, ax in zip(g.col_names, g.axes.flatten()):
        # If we have an entry for col_name, iterate over its values.
        if col_name in values:
            for item in values[col_name]:
                try:
                    # Highlight corresponding patch.
                    idx = row_names.index(item)
                    ax.patches[idx].set_facecolor(color)
                except ValueError:
                    # Skip any entries not in plot.
                    pass


def _draw_legend(color_map, ax, **kwargs):
    patches = [
        mpl_patches.Patch(
            color=color, label=label) for label, color in color_map.items()
    ]
    legend = ax.legend(handles=patches, **kwargs)
    return legend


def plot_gene_clonality(insertions,
                        ax=None,
                        label_min_freq=None,
                        label_offsets=None,
                        label_kws=None,
                        **kwargs):
    if ax is None:
        _, ax = plt.subplots()

    # Calculate frequency and mean clonality.
    data = gene_statistics(insertions).reset_index()

    # Draw plot.
    ax.plot(data['n_samples'], data['mean_clonality'], '.', **kwargs)

    ax.set_xlabel('Number of samples')
    ax.set_ylabel('Mean clonality')
    ax.set_title('Candidate frequency vs. clonality')

    ax.set_xlim(0, None)
    # ax.set_ylim(0, 1)
    sns.despine(ax=ax)

    # Draw labels.
    if label_min_freq is None:
        label_data = data
    else:
        label_data = data.ix[data['n_samples'] >= label_min_freq]
        ax.axvline(
            label_min_freq - 0.5,
            color='darkgrey',
            linestyle='dashed',
            zorder=0)

    _draw_labels(label_data, ax=ax, offsets=label_offsets, **(label_kws or {}))

    return ax


def _draw_labels(data,
                 ax,
                 offsets=None,
                 x='n_samples',
                 y='mean_clonality',
                 label='gene_name',
                 **kwargs):

    offsets = offsets or {}

    for _, row in data.iterrows():
        offset = offsets.get(row[label], None)

        if offset is None:
            offset, arrowprops = (5, 0), dict()
        else:
            arrowprops = dict(arrowstyle="-", lw=0.5)

        ax.annotate(
            row[label],
            va='center',
            xy=(row[x], row[y]),
            xycoords='data',
            xytext=offset,
            textcoords='offset points',
            arrowprops=arrowprops,
            **kwargs)


def plot_gene_clonality_boxplot(insertions, order=None, ax=None, **kwargs):
    if ax is None:
        _, ax = plt.subplots()

    if order is None:
        # Order genes by their frequency.
        order = list(gene_sample_count(insertions)
                     .sort_values(ascending=False).index) # yapf: disable

    # Plot clonality of insertions.
    sns.boxplot(
        data=insertions,
        x='gene_name',
        y='clonality',
        ax=ax,
        color='lightgrey',
        order=order,
        **kwargs)

    ax.set_xlabel('Genes')
    ax.set_ylabel('Clonality')
    ax.set_ylim(0, 1.02)

    plt.setp(ax.get_xticklabels(), rotation=90, fontstyle='italic')
    sns.despine(ax=ax)

    return ax


def plot_orientation_bias(insertions,
                          label_offsets=None,
                          min_samples=5,
                          figsize=(12, 5.5)):
    label_offsets = label_offsets or {}

    stats = gene_statistics(insertions).reset_index()

    fig, ax = plt.subplots(figsize=figsize)
    sns.regplot(
        data=stats,
        x='n_samples',
        y='sense_fraction_weighted',
        ax=ax,
        fit_reg=False,
        scatter_kws={'alpha': 1})

    for _, row in stats.query('n_samples > 5').iterrows():
        if row.n_samples > min_samples:
            offset = label_offsets.get(row['gene_name'], (5, 0))

            if offset is None:
                arrow_props = {}
            else:
                arrow_props = {'arrowstyle': "-", 'lw': 1}

            ax.annotate(
                row['gene_name'],
                va='center',
                xy=(row.n_samples, row.sense_fraction_weighted),
                xycoords='data',
                textcoords='offset points' if offset is not None else None,
                xytext=offset,
                arrowprops=arrow_props,
                fontstyle='italic')

    ax.set_xlim(0, None)
    ax.set_ylim(-0.05, 1.05)

    ax.set_xlabel('Number of samples')
    ax.set_ylabel('Sense fraction (weighted)')

    ax.axvline(min_samples, color='darkgrey', linestyle='dashed', zorder=-1)

    ax.axhline(
        0.5, color=sns.color_palette()[2], alpha=0.8,
        linestyle='dashed', zorder=-1) # yapf: disable

    return fig


def plot_insertion_track(insertions,
                         region,
                         gene=None,
                         transcript_id=None,
                         ins_ratio=1 / 25,
                         linewidth=0.5,
                         **kwargs):
    ori_order = ['sense', 'antisense']
    ori_palette = [sns.color_palette()[0], sns.color_palette()[2]]

    width = (region[2] - region[1]) * ins_ratio

    ins_track = FeatureTrack.from_position(
        data=insertions,
        width=width,
        height=0.25,
        hue='gene_orientation',
        hue_order=ori_order,
        palette=ori_palette,
        patch_kws={'edgecolor': 'black',
                   'linewidth': linewidth})

    rug_track = RugTrack(
        data=insertions,
        height=0.25,
        hue='gene_orientation',
        hue_order=ori_order,
        palette=ori_palette)

    if gene is not None:
        filter = 'gene_name == {!r}'.format(gene)
    elif transcript_id is not None:
        filter = 'transcript_id == {!r}'.format(transcript_id)
    else:
        filter = None

    gene_track = BiomartTrack(
        dataset='mmusculus_gene_ensembl',
        height=0.4,
        collapse='transcript',
        filter=filter,
        gene_id='gene_name',
        patch_kws={'linewidth': linewidth},
        line_kws={'lw': 1},
        label_kws={'fontstyle': 'italic'})

    fig = plot_tracks(
        [ins_track, rug_track, gene_track],
        region=region,
        despine=True,
        **kwargs)

    fig.axes[0].set_title('Insertions')
    fig.axes[-1].set_xlabel('Chromosome {}'.format(region[0]))

    return fig


def plot_exon_expression(coverage,
                         gene,
                         z_score=None,
                         insertions=None,
                         **kwargs):

    # Subset counts to gene.
    coverage = coverage.loc[[gene]]

    if insertions is not None:
        # Summarize insertions.
        ins_summ = _assign_exon_positions(insertions, coverage, gene=gene)
        ins_summ = ins_summ.ix[ins_summ['index'] > 0]
        ins_summ = ins_summ.sort_values(['index', 'clonality'])

        # Subset expression to samples with insertions.
        coverage = coverage[ins_summ['sample'].unique()]

        # Determine row colors based on clonality.
        ins_clon = ins_summ.groupby('sample')['clonality'].max()
        row_colors = color_annotation(
            ins_clon.to_frame(), colors=['black'], vlims=[(0, 1)])[0]
        row_colors = row_colors.rename(columns=lambda s: s.capitalize())
    else:
        ins_summ = None
        row_colors = None

    # Plot expression.
    g = _plot_exon_coverage(
        coverage, gene, z_score, row_colors=row_colors, **kwargs)

    # Annotate insertion sites.
    if ins_summ is not None:
        samples = list(coverage.columns)[::-1]

        segments = []
        for tup in ins_summ.itertuples():
            i = samples.index(tup.sample)
            segments.append([(tup.index, i), (tup.index, i + 1)])

    line_segments = mpl_coll.LineCollection(
        segments, linewidths=(2, ), linestyles='solid', color='black')
    g.ax_heatmap.add_collection(line_segments)

    return g


def _assign_exon_positions(insertions, coverage, gene):
    # Get sample, position and clonality fields.
    ins_gene = insertions.query('gene_name == {!r}'.format(gene))
    summary = ins_gene[['sample', 'position', 'clonality']].copy()

    # Determine orientation of gene.
    first = ins_gene.iloc[0]
    relative_ori = 1 if first.gene_orientation == 'sense' else -1
    gene_strand = first.strand * relative_ori

    # Determine index position of insertion.
    coverage = coverage.ix[[gene]]
    exon_ends = sorted(coverage.index.get_level_values('end'))

    summary['index'] = summary['position'].map(
        lambda pos: bisect.bisect(exon_ends, pos))

    if gene_strand == -1:
        summary['index'] = len(exon_ends) - summary['index']

    return summary


def _plot_exon_coverage(coverage, transcript_id, z_score=None, **kwargs):
    # Transform expression, taking z-score if requested.
    expr = np.log2(coverage.loc[transcript_id] + 1)

    if z_score is not None:
        expr = _calc_zscore(expr, axis=z_score)

    # Sort by ascending exons.
    strand = expr.index.get_level_values(3)[0]
    expr.sort_index(ascending=strand == '+', inplace=True)

    # Draw heatmap.
    g = sns.clustermap(expr.T, row_cluster=False, col_cluster=False, **kwargs)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    # Draw xticks.
    xticks = [1] + list(range(5, len(expr) + 1, 5))
    g.ax_heatmap.set_xticks(np.array(xticks) - 0.5)
    g.ax_heatmap.set_xticklabels(map(str, xticks), rotation=0)
    g.ax_heatmap.set_yticks([])

    # Add axis labels and title.
    g.ax_heatmap.set_xlabel('Exons')
    g.ax_heatmap.set_ylabel('Samples')
    g.ax_heatmap.set_title('{} exon expression'.format(transcript_id))

    # Draw border around heatmap.
    _add_axis_border(g.ax_heatmap, linewidth=0.5)

    # Draw border around annotation.
    if kwargs.get('row_colors', None) is not None:
        _add_axis_border(g.ax_row_colors, linewidth=0.5)

    return g


def _calc_zscore(data2d, axis=1):
    """Standarize the mean and variance of the data axis."""
    other_axis = 0 if axis == 1 else 1
    return (data2d.subtract(
        data2d.mean(axis=other_axis), axis=axis).divide(
            data2d.std(axis=other_axis), axis=axis))


def _add_axis_border(ax, linewidth=1):
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(linewidth)
