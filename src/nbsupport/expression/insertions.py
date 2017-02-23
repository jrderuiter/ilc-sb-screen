import bisect

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mpl_coll
import seaborn as sns

# from imfusion.expression.de_test import de_exon
# from imfusion.util.tabix import GtfFile as IMFGtfFile

from ..util import color_annotation


def plot_de(insertions,
            dexseq_gtf,
            exon_counts_path,
            gene_id,
            gene_name=None,
            ax=None):
    """
    Plots differential expression over a genes insertion sites using IM-Fusion.

    Parameters
    ----------
    insertions : pandas.DataFrame
        Insertions in DataFrame format.
    dexseq_gtf : GtfFile or str/Path
        (Path to) DEXSeq gtf file containing gene exons for which expression
        counts were generated. Note that genes should not have been aggregated
        in the gtf file.
    exon_counts_path : str or Path
        Path to a TSV containing exon expression counts. Exon definitions are
        expected to correspond with the given gtf file.
    gene_id : str
        (Ensembl) ID of the gene to plot for.
    gene_name : str, optional
        Name to use for the gene in the title of the plot. If not given,
        the value of gene_id is used.
    ax : matplotlib.Axes
        Axis on which the plot is drawn.

    Returns
    -------
    matplotlib.Axes
        Axis on which was drawn.

    """

    if ax is None:
        _, ax = plt.subplots()

    if not isinstance(dexseq_gtf, IMFGtfFile):
        dexseq_gtf = IMFGtfFile(dexseq_gtf)

    de_result = de_exon(
        insertions.rename(columns={'sample': 'sample_id'}),
        gene_id=gene_id,
        dexseq_gtf=dexseq_gtf,
        exon_counts=exon_counts_path)

    de_result.plot_boxplot(log=True, ax=ax, color='lightgrey')
    ax.set_title('{} (p = {:.2e})'.format(gene_name or gene_id,
                                          de_result.p_value))

    ax.set_xticklabels(['With\ninsertion', 'Without\ninsertion'])
    ax.set_xlabel('')

    sns.despine(ax=ax)

    return ax


def plot_exon_expression(exon_counts,
                         gene_id,
                         z_score=None,
                         insertions=None,
                         **kwargs):

    # Subset counts to gene.
    exon_counts = exon_counts.ix[[gene_id]]

    if insertions is not None:
        # Summarize insertions.
        ins_summ = _summarize_insertions(insertions, exon_counts, gene_id)
        ins_summ = ins_summ.ix[ins_summ['index'] > 0]
        ins_summ = ins_summ.sort_values(['index', 'clonality'])

        # Subset expression to samples with insertions.
        exon_counts = exon_counts[ins_summ['sample'].unique()]

        # Determine row colors based on clonality.
        ins_clon = ins_summ.groupby('sample')['clonality'].max()
        row_colors = color_annotation(
            ins_clon.to_frame(), colors=['black'], vlims=[(0, 1)])[0]
        row_colors = row_colors.rename(columns=lambda s: s.capitalize())
    else:
        ins_summ = None
        row_colors = None

    # Plot expression.
    g = _plot_exon_expression(
        exon_counts, gene_id, z_score, row_colors=row_colors, **kwargs)

    # Annotate insertion sites.
    if ins_summ is not None:
        samples = list(exon_counts.columns)[::-1]

        segments = []
        for tup in ins_summ.itertuples():
            i = samples.index(tup.sample)
            segments.append([(tup.index, i), (tup.index, i + 1)])

    line_segments = mpl_coll.LineCollection(
        segments, linewidths=(2, ), linestyles='solid', color='black')
    g.ax_heatmap.add_collection(line_segments)

    return g


def _summarize_insertions(insertions, exon_counts, gene_id):
    # Get sample, position and clonality fields.
    ins_gene = insertions.query('gene_id == {!r}'.format(gene_id))
    summary = ins_gene[['sample', 'position', 'clonality']].copy()

    # Determine orientation of gene.
    first = ins_gene.iloc[0]
    relative_ori = 1 if first.gene_orientation == 'sense' else -1
    gene_strand = first.strand * relative_ori

    # Determine index position of insertion.
    exon_counts = exon_counts.ix[[gene_id]]
    exon_ends = sorted(exon_counts.index.get_level_values('end'))

    summary['index'] = summary['position'].map(
        lambda pos: bisect.bisect(exon_ends, pos))

    if gene_strand == -1:
        summary['index'] = len(exon_ends) - summary['index']

    return summary


def _plot_exon_expression(exon_counts, gene_id, z_score=None, **kwargs):
    # Transform expression, taking z-score if requested.
    expr = np.log2(exon_counts.ix[gene_id] + 1)

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
    g.ax_heatmap.set_title('{} exon expression'.format(gene_id))

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
