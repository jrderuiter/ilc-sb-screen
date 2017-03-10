from matplotlib.colors import ListedColormap
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import seaborn as sns
from statsmodels.sandbox.stats.multicomp import multipletests

from .util.clustermap import sort_matrix

MORPHOLOGY_ORDER = ['ILC', 'Spindle cell', 'Squamous']
MORPHOLOGY_COLORS = ['#66c2a5', '#fa8d62', '#e68ac3']

METASTASIS_COLOR = '#8d9fca'


def parse_morphology(samples, pathology='pathology_type'):
    """Parses pathology labels into predefined morphologies.

    Parses the pathology column of the sample overview into a 2D matrix,
    indicating the membership of each sample to three predefined morphologies:
    'ILC', 'Spindle cell' and 'Squamous'. Note that samples can have multiple
    morphologies.

    Parameters
    ----------
    samples : pd.DataFrame
        Sample overview. Is expected to contain a 'sample' column and the
        pathology column (which defaults to 'pathology_type').
    pathology : str
        Name of the column containing sample pathology.

    Returns
    -------
    pd.DataFrame
        Boolean matrix indicating the membership of each sample to each
        of the three morphologies.

    """

    pathology = samples[pathology].str.lower()

    morphology = pd.DataFrame(
        {
            'sample': samples['sample'],
            'ILC': pathology.str.contains('ilc'),
            'Spindle cell': pathology.str.contains('spindle cell'),
            'Squamous': pathology.str.contains('squamous')
        },
        columns=['sample'] + MORPHOLOGY_ORDER)

    morphology = morphology.set_index('sample').fillna(False)

    return sort_matrix(morphology, sort_columns=False)


def plot_morphology(morphology, sample_metastases=None, ax=None, **kwargs):
    """Plots morphology matrix as 2D heatmap.

    Plots a morphology matrix (typically obtained from the parse_morphology
    function) as a 2D heatmap. Matrix is expected to correspond with the
    three categories returned by parse_morphology (ILC, Spindle cell
    and Squamous).

    Parameters
    ----------
    morphology : pd.DataFrame
        Boolean matrix of samples-by-morphologies.
    sample_metastases : pd.DataFrame
        Optional dataframe (single column) indicating which samples have
        a metastasis. Used to draw metastases as an extra row in the heatmap.
    ax : matplotlib.Axis
        Axis to use for plotting.
    **kwargs
        Any kwargs are passed to seaborns heatmap function.

    Returns
    -------
    matplotlib.Axis
        Axis that was used for plotting.

    """

    if ax is None:
        _, ax = plt.subplots()

    # Sort by rows/columns.
    morphology = morphology[MORPHOLOGY_ORDER]

    if sample_metastases is not None:
        morphology = pd.concat([morphology, sample_metastases], axis=1)
        palette = MORPHOLOGY_COLORS + [METASTASIS_COLOR]
    else:
        palette = MORPHOLOGY_COLORS

    morphology = sort_matrix(morphology, sort_columns=False)

    # Convert to numeric matrix (for heatmap).
    num_matrix = pd.DataFrame(
        {
            col: morphology[col].astype(int) * (i + 1)
            for i, col in enumerate(morphology)
        },
        columns=morphology.columns)

    # Draw heatmap.
    cmap = ListedColormap(['#f6f6f6'] + palette)
    sns.heatmap(
        num_matrix.T,
        ax=ax,
        cbar=False,
        cmap=cmap,
        vmin=0,
        vmax=len(palette),
        **kwargs)

    ax.set_xticks([])
    ax.set_xlim(0, num_matrix.shape[0])

    ax.set_title('Tumor morphology')
    ax.set_xlabel('Samples')

    # Add counts to labels.
    ax.set_yticklabels(
        ['{} ({})'.format(k, v) for k, v in morphology.sum().items()][::-1],
        rotation=0)

    return ax


def plot_metastases(metastases, ax=None, **kwargs):
    """Plots matrix of metastasis sites as 2D heatmap.

    Plots a matrix of samples-by-metastasis-sites as a 2D heatmap,
    indicating which samples have metastases to which sites.

    Parameters
    ----------
    metastases : pd.DataFrame
        Boolean 2D matrix of samples-by-met-sites.
    ax : matplotlib.Axis
        Axis to use for plotting.
    **kwargs
        Any kwargs are passed to seaborns heatmap function.

    Returns
    -------
    ax : matplotlib.Axis
        Axis that was used for plotting.

    """

    if ax is None:
        _, ax = plt.subplots()

    # Sort by frequency.
    metastases = sort_matrix(metastases).T

    # Draw heatmap.
    cmap = ListedColormap(['#F6F6F6', sns.color_palette()[0]])
    sns.heatmap(metastases, ax=ax, cbar=False, cmap=cmap, **kwargs)

    # Set ticks/labels/title.
    ax.set_xticks([])
    ax.set_xlabel('Mice with metastases ({})'.format(metastases.shape[0]))
    ax.set_title('Metastasis sites')

    # Add counts to labels.
    ax.set_yticklabels([
        '{} ({})'.format(k, v) for k, v in metastases.sum(axis=1).items()
    ][::-1])

    return ax


def test_strain_bias(data,
                     value,
                     sample='sample',
                     strain='t2onc_type',
                     samples=None,
                     incl_neg=True):
    """Tests if a given property is biased towards a specific mouse strain.

    Test is performed using a Fishers Exact test, corrected for multiple
    testing using Benjamini-Hochberg multiple testing correction.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing (among others) the value to be tested,
        together with sample names and the sample strain.
    value : str
        Name of the value column, the (categorical) values of which
        will be tested for association with the strain.
    sample : str
        Name of the column containing sample names.
    strain : str
        Name of the column containing the strain names.
    samples : pd.DataFrame
        Optional dataframe containing the sample overview. If given,
        this frame is used to calculate the total number of samples
        per strain instead of the 'data' dataframe. This is useful in
        cases where data does not necessarily contain all samples.
    incl_neg : bool
        Whether to return the negative counts (number of samples without
        the corresponding value) in the result DataFrame.

    Returns
    -------
    pd.DataFrame
        Dataframe containing the test results, with one row for
        each distinct value in the value column.

    """
    assert data[strain].nunique() == 2

    # Build contingency table.
    if samples is None:
        strain_counts = data.groupby(strain)[sample].nunique()
    else:
        strain_counts = samples.groupby(strain)[sample].nunique()

    with_counts = (data.groupby([strain, value])[sample].nunique()
                   .reset_index(name='count').pipe(
                       pd.pivot_table,
                       index=value,
                       columns=strain,
                       values='count').fillna(0.0))

    without_counts = strain_counts - with_counts

    contingency_table = pd.concat(
        [
            with_counts.rename(columns=lambda c: 'pos_' + str(c)),
            without_counts.rename(columns=lambda c: 'neg_' + str(c))
        ],
        axis=1).astype(int)

    # Calculate p-values.
    result = contingency_table.copy()
    result.columns.name = None

    result['p_value'] = [
        fisher_exact(np.reshape(tup[1:], [2, 2]))[1]
        for tup in result.itertuples()
    ]

    # Calculate q-values.
    result['q_value'] = multipletests(result['p_value'], method='fdr_bh')[1]

    if not incl_neg:
        result = result.drop(
            [c for c in result.columns if c.startswith('neg')], axis=1)

    return result
