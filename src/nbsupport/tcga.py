import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from scipy import stats


def parse_barcode(barcode):
    """Parse a tcga barcode into its constituents."""

    # Define fields and their locations in the barcode.
    fields = {
        'project': (0, ),
        'tss': (1, ),
        'participant': (2, ),
        'sample': (3, (0, 2)),
        'vial': (3, (2, )),
        'portion': (4, (0, 2)),
        'analyte': (4, (2, )),
        'plate': (5, ),
        'center': (6, )
    }

    # Split barcode and extract fields.
    split = barcode.split('-')

    values = {
        field: _extract_field(split, indices)
        for field, indices in fields.items()
    }

    return values


def _extract_field(split, indices):
    try:
        if len(indices) == 1:
            value = split[indices[0]]
        elif len(indices) == 2:
            ind1, ind2 = indices
            if len(ind2) == 1:
                value = split[ind1][ind2[0]]
            else:
                value = split[ind1][ind2[0]:ind2[1]]
        return value
    except IndexError:
        return None


def extract_sample(barcode):
    parsed = parse_barcode(barcode)
    return '-'.join([
        parsed['project'], parsed['tss'], parsed['participant'],
        parsed['sample']
    ])


def extract_particpant(barcode):
    parsed = parse_barcode(barcode)
    return '-'.join([parsed['project'], parsed['tss'], parsed['participant']])


def plot_expr_vs_cn(expr,
                    cnv,
                    gene_name,
                    show_points=True,
                    ax=None,
                    boxplot_kws=None,
                    strip_kws=None,
                    label_kws=None):
    merged = pd.merge(
        cnv.ix[gene_name].to_frame('copy_number'),
        expr.ix[gene_name].to_frame('expr'),
        left_index=True,
        right_index=True)

    ax = ax or plt.subplots()[1]

    palette = {
        -1: (0.65, 0.81, 0.89),
        -2: (0.13, 0.47, 0.71),
        1: (0.98, 0.60, 0.59),
        2: (0.89, 0.10, 0.11),
        0: 'lightgrey'
    }

    order = [-2, -1, 0, 1, 2]

    sns.boxplot(
        data=merged,
        x='copy_number',
        y='expr',
        ax=ax,
        palette=palette,
        showfliers=False,
        order=order,
        **(boxplot_kws or {}))

    if show_points:
        sns.swarmplot(
            data=merged,
            x='copy_number',
            y='expr',
            ax=ax,
            color='black',
            order=order,
            **(strip_kws or {}))

    ax.set_xlabel('Copy number status (Gistic)')
    ax.set_ylabel('Expression (log2)')

    ax.set_title('{} expression vs. copy number'.format(gene_name))
    sns.despine(ax=ax)

    test = stats.spearmanr(merged['copy_number'], merged['expr'])
    label = ('\u03C1 = {:3.2f}\np-value = {:3.2e}'
             .format(test.correlation, test.pvalue))

    ax.annotate(
        s=label,
        xy=(0.05, 0.88),
        xycoords='axes fraction',
        ha='left',
        **(label_kws or {}))

    labels = ['Deletion', 'Het loss', 'Neutral', 'Gain', 'Amplification']
    ax.set_xticklabels(labels)

    return ax
