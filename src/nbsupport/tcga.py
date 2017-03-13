import collections

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from scipy import stats


def parse_barcode(barcode):
    """Parse a tcga barcode into its constituents.

    Parameters
    ----------
    barcode : str
        TCGA barcode to parse.

    Returns
    -------
    dict[str, str]
        Dictionary containing the parsed barcode elements.

    """

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


def extract_sample_barcode(barcode):
    """Extracts the sample barcode from a TCGA barcode.

    Extracts the sample barcode from TCGA barcodes, which is comprised
    of the first four elements of the barcode.

    Parameters
    ----------
    barcode : str
        TCGA barcode to parse.

    Returns
    -------
    str
        Parsed sample barcode.

    """

    parsed = parse_barcode(barcode)
    return '-'.join([
        parsed['project'], parsed['tss'], parsed['participant'],
        parsed['sample']
    ])


def extract_participant_barcode(barcode):
    """Extracts the participant barcode from a TCGA barcode.

    Extracts the sample participant from TCGA barcodes, which is comprised
    of the first three elements of the barcode.

    Parameters
    ----------
    barcode : str
        TCGA barcode to parse.

    Returns
    -------
    str
        Parsed participant barcode.

    """
    parsed = parse_barcode(barcode)
    return '-'.join([parsed['project'], parsed['tss'], parsed['participant']])


def select_tumor_samples(data):
    """Selects columns corresponding to TCGA tumor samples.

    Parameters
    ----------
    data : pd.DataFrame
        Dataframe containing TCGA data, with samples along the columns.
        Samples should be named using TCGA barcodes.

    Returns
    -------
    pd.DataFrame
        Subsetted dataframe, containing only tumor samples.

    """

    def _is_tumor(barcode):
        return parse_barcode(barcode)['sample'].startswith('0')

    return data.get([c for c in data.columns if _is_tumor(c)])


def drop_duplicate_columns(data):
    """Drops duplicate columns from a dataframe.

    Parameters
    ----------
    data : pd.DataFrame
        Dataframe to filter.

    Returns
    -------
    pd.DataFrame
        Subsetted dataframe, without any duplicate column names.

    """

    counts = collections.Counter(data.columns)
    duplicates = [k for k, v in counts.items() if v > 1]
    return data.drop(duplicates, axis=1)


def plot_cnv_expr_corr(expr,
                       cnv,
                       gene_name,
                       show_points=True,
                       ax=None,
                       boxplot_kws=None,
                       strip_kws=None,
                       label_kws=None):
    """Plots correlation between CNV calls and gene expression values.

    """

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
