from matplotlib.colors import ListedColormap
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns

from .util.clustermap import sort_matrix

MORPHOLOGY_ORDER = ['ILC', 'Spindle cell', 'Squamous']
MORPHOLOGY_COLORS = ['#66c2a5', '#fa8d62', '#e68ac3']

METASTASIS_COLOR = '#8d9fca'


def parse_morphology(samples):
    """Parses morphology labels into predefined categories."""

    pathology = samples['pathology_type'].str.lower()

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


def plot_morphology(morphology, sample_metastases=None, ax=None):
    """Plots morphology labels as 2D matrix."""

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
    sns.heatmap(num_matrix.T, ax=ax, cbar=False, cmap=cmap)

    ax.set_xticks([])
    ax.set_xlim(0, num_matrix.shape[0])

    ax.set_title('Tumor morphology')
    ax.set_xlabel('Samples')

    # Add counts to labels.
    ax.set_yticklabels(
        ['{} ({})'.format(k, v) for k, v in morphology.sum().items()][::-1])

    return ax


def plot_metastases(metastases, ax=None, **kwargs):
    """Plots metastasis sites as 2d boolean matrix."""

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
