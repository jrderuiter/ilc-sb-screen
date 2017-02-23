from collections import OrderedDict

import matplotlib.pyplot as plt
import matplotlib.patches as mpl_patches
import numpy as np
import pandas as pd
import seaborn as sns

import pybiomart
import toolz

from geneviz.tracks import (GtfTrack, BiomartTrack, RugTrack, FeatureTrack,
                            plot_tracks)


def annotate_with_clonality(insertions):
    grps = []
    for _, grp in insertions.groupby('sample'):
        grp = grp.copy()
        grp['clonality'] = grp['support'] / grp['support'].max()

        grps.append(grp)
    return pd.concat(grps, axis=0)


def plot_insertion_stats(insertions, fig_kws=None, suffix='', color=None):
    # Plot some more advanced statistics of dataset.
    fig, axes = plt.subplots(ncols=2, nrows=2, **(fig_kws or {}))
    axes = axes.flatten()

    # Insertion depth.
    ins_depth = insertions.groupby('id')['depth_unique'].max()
    ins_depth.hist(bins=np.arange(0, np.max(ins_depth), 10), ax=axes[0],
                   grid=False, color=color)
    axes[0].set_title('Insertion depth' + suffix)
    axes[0].set_xlabel('Insertion depth')
    axes[0].set_ylabel('Number of insertions')

    # Maximum depth per sample.
    ins_depth_sample = insertions.groupby('sample')['depth_unique'].max()
    ins_depth_sample.hist(bins=np.arange(0, 180, 10), ax=axes[1],
                          grid=False, color=color)
    axes[1].set_title('Maximum insertion depth per sample' + suffix)
    axes[1].set_xlabel('Maximum insertion depth')
    axes[1].set_ylabel('Number of samples')

    # Insertion clonality
    ins_clon = insertions.groupby('id')['clonality'].max()
    ins_clon.hist(bins=20, ax=axes[2], grid=False, color=color)
    axes[2].set_title('Insertion clonality' + suffix)
    axes[2].set_xlabel('Insertion clonality')
    axes[2].set_ylabel('Number of insertions')

    # Number of insertion per sample.
    ins_sample = insertions.groupby('sample')['id'].nunique()
    ins_sample.hist(bins=15, ax=axes[3], grid=False, color=color)
    axes[3].axvline(ins_sample.median(), linestyle='dashed',
                    color='red')
    axes[3].set_title('Number of insertions per sample' + suffix)
    axes[3].set_xlabel('Number of insertions')
    axes[3].set_ylabel('Number of samples')

    sns.despine(fig)
    fig.tight_layout()

    return fig, axes


def gene_statistics(insertions):
    return pd.DataFrame(
        {
            'n_samples': gene_sample_count(insertions),
            'sense_fraction': gene_sense_fraction(insertions),
            'sense_fraction_weighted': gene_sense_fraction_weighted(insertions),
            'mean_clonality': gene_clonality(insertions)
        },
        columns=['n_samples', 'mean_clonality', 'sense_fraction',
                 'sense_fraction_weighted'])


def gene_sample_count(insertions, name='n_samples'):
    count = insertions.groupby('gene_name')['sample'].nunique()
    count.name = name
    return count


def gene_sense_fraction(insertions):
    def _sense_frac(x):
        is_sense = x['gene_orientation'] == 'sense'
        return np.sum(is_sense) / len(x)

    return insertions.groupby(['gene_name']).apply(_sense_frac)


def gene_sense_fraction_weighted(insertions):
    def _sense_frac_weighted(x):
        is_sense = x['gene_orientation'] == 'sense'
        return np.sum(is_sense * x['clonality']) / np.sum(x['clonality'])

    return insertions.groupby(['gene_name']).apply(_sense_frac_weighted)


def gene_clonality(insertions):
    return clonality_matrix(insertions).mean(axis=1)


def clonality_matrix(insertions, value='clonality', fill_value=None):
    """Summarizes insertions in gene/sample matrix."""

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
    """Sorts a 2d matrix, first by its rows and then by its columns."""

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
    """Plots a 2d sample/by gene overview of an insertion matrix."""

    if ax is None:
        _, ax = plt.subplots()

    # Convert to matrix.
    ins_matrix = clonality_matrix(insertions, value=value)

    # Subset if needed.
    if genes is not None:
        ins_matrix = ins_matrix.ix[genes]
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
    """Draws bar factorplot counting gene occurrences per subtype."""

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
        for p in ax.patches:
            fraction = p.get_width()
            count = int(round(fraction * total))

            ax.text(
                x=fraction + 0.02,
                y=(p.get_y() + (p.get_height() / 2)),
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
    patches = [mpl_patches.Patch(
        color=color, label=label) for label, color in color_map.items()]
    legend = ax.legend(handles=patches, **kwargs)
    return legend


def plot_gene_clonality(insertions,
                        ax=None,
                        label_min_freq=0,
                        label_extra=None,
                        label_offsets=None,
                        label_kws=None):
    if ax is None:
        _, ax = plt.subplots()

    # Calculate frequency and mean clonality.
    data = gene_statistics(insertions).reset_index()

    # Draw plot.
    sns.regplot(
        data=data, x='n_samples', y='mean_clonality', ax=ax, fit_reg=False,
        scatter_kws={'alpha': 1})

    ax.set_xlabel('Number of samples')
    ax.set_ylabel('Mean clonality')
    ax.set_title('Candidate frequency vs. clonality')

    ax.set_xlim(0, None)
    sns.despine(ax=ax)

    # Draw labels.
    label_data = data.ix[data['n_samples'] >= label_min_freq]
    _draw_labels(label_data, ax=ax, offsets=label_offsets, **(label_kws or {}))

    # Draw labels for extra genes.
    if label_extra is not None:
        label_data_extra = data.ix[data['gene_name'].isin(label_extra)]
        _draw_labels(
            label_data_extra,
            ax=ax,
            offsets=label_offsets,
            color='grey',
            **(label_kws or {}))

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


DEFAULT_GENE_KWS = {
    'collapse': 'transcript',
    'height': 0.3,
    'plot_kws': {'lw': 0.5},
    'arrow_size': 70,
    'arrow_spacing': 40,
    'stack_kws': {'spacing': 0.1}
}

DEFAULT_INS_KWS = dict(
    hue='gene_orientation',
    palette={
        'sense': sns.color_palette()[0],
        'antisense': sns.color_palette()[2]
    },
    height=0.2,
    plot_kws={'lw': 0.5},
    stack_kws=dict(spacing=0.01))

DEFAULT_RUG_KWS = dict(
    hue='gene_orientation',
    palette={'sense': sns.color_palette()[0],
             'antisense': sns.color_palette()[2]},
    height=0.25)


def plot_insertion_track(insertions, gene_name, track_kws=None, **kwargs):
    """Draws insertions above a gene track from Biomart."""

    track_kws = track_kws or {}

    # Get gene position.
    dataset = pybiomart.Dataset(
        name='mmusculus_gene_ensembl', host='http://www.ensembl.org')

    genes = dataset.query(
        ['external_gene_name', 'chromosome_name', 'start_position',
         'end_position', 'strand'],
        use_attr_names=True)

    genes = genes.rename(columns={
        'start_position': 'start',
        'end_position': 'end',
        'chromosome_name': 'contig'
    })

    gene = genes.set_index('external_gene_name').ix[gene_name]

    # Setup gene track.
    gene_kws = toolz.merge(DEFAULT_GENE_KWS, track_kws.get('gene', {}))
    gene_track = BiomartTrack(
        dataset='mmusculus_gene_ensembl',
        filter='gene_name == {!r}'.format(gene.name),
        **gene_kws)

    # Draw insertions together with this gene track.
    return _plot_insertions(
        insertions, gene_track, gene, track_kws=track_kws, **kwargs)


def _plot_insertions(insertions,
                     gene_track,
                     gene,
                     track_kws=None,
                     ins_size=1 / 120,
                     **kwargs):
    """Shared function for drawing non-gene tracks."""

    track_kws = track_kws or {}

    rug_kws = toolz.merge(DEFAULT_RUG_KWS, track_kws.get('rug', {}))
    rug_track = RugTrack(
        data=insertions.rename(columns={'chromosome': 'seqname'}), **rug_kws)

    # Determine insertions width.
    padding = kwargs.get('padding', (0, 0))
    insertion_width = (abs(gene.end - gene.start) + sum(padding)) * ins_size

    # Build insertion track.
    ins_data = _preprocess_insertions(insertions, insertion_width)

    ins_kws = toolz.merge(DEFAULT_INS_KWS, track_kws.get('insertions', {}))
    ins_track = FeatureTrack(data=ins_data, **ins_kws)

    reverse = gene.strand != 1

    return plot_tracks(
        [ins_track, rug_track, gene_track],
        seqname=gene.contig,
        start=gene.start,
        end=gene.end,
        despine=True,
        tick_top=False,
        reverse=reverse,
        **kwargs)


# def plot_insertion_track_gtf(insertions,
#                              gtf_path,
#                              gene_name,
#                              track_kws=None,
#                              **kwargs):
#     """Draws insertions using a GTF file for gene definitions."""

#     track_kws = track_kws or {}

#     # Get gene position.
#     gtf_file = GtfFile(gtf_path)
#     gene = gtf_file.get_gene(gene_name, field_name='gene_name')
#     gene.strand = 1 if gene.strand == '+' else -1

#     # Setup tracks.
#     gene_kws = toolz.merge(DEFAULT_GENE_KWS, track_kws.get('gene', {}))
#     gene_track = GtfTrack.from_path(
#         gtf_path, filters='gene_name == {!r}'.format(gene_name), **gene_kws)

#     return _plot_insertions(
#         insertions, gene_track, gene, track_kws=track_kws, **kwargs)


def _preprocess_insertions(insertions, insertion_width):
    return (insertions.assign(
        start=lambda df: df['position'] - insertion_width,
        end=lambda df: df['position'] + insertion_width)
            .rename(columns={'chromosome': 'seqname'}))