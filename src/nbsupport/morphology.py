
from matplotlib.colors import ListedColormap
import pandas as pd
import seaborn as sns

from .util import sort_matrix


def plot_morphology(morphology, palette, ax=None, sort_columns=True):
    if ax is None:
        _, ax = plt.subplots()

    # Sort by rows/columns.
    morphology = sort_matrix(morphology, sort_columns=sort_columns)

    # Convert to numeric matrix (for heatmap).
    num_matrix = pd.DataFrame(
        {col: morphology[col].astype(int) * (i + 1)
         for i, col in enumerate(morphology)},
        columns=morphology.columns)

    # Draw heatmap.
    cmap = ListedColormap(['#f6f6f6'] + palette[:num_matrix.shape[1]])
    sns.heatmap(num_matrix.T, ax=ax, cbar=False, cmap=cmap)

    ax.set_xticks([])
    ax.set_xlim(0, num_matrix.shape[0])

    ax.set_title('Tumor morphology')
    ax.set_xlabel('Samples')

    # Add counts to labels.
    ax.set_yticklabels(['{} ({})'.format(k, v)
                       for k, v in morphology.sum().items()][::-1])

    return ax


def plot_metastases(metastases, ax=None, **kwargs):
    """Plots metastasis sites as 2d boolean matrix."""

    if ax is None:
        _, ax = plt.subplots()

    # Sort by frequency.
    metastases = sort_matrix(metastases).T

    # Draw heatmap.
    cmap = ListedColormap(['#F6F6F6', sns.color_palette()[0]])
    sns.heatmap(metastases, ax=ax, cbar=False,
                cmap=cmap, **kwargs)

    # Set ticks/labels/title.
    ax.set_xticks([])
    ax.set_xlabel('Mice with metastases ({})'
                  .format(metastases.shape[0]));
    ax.set_title('Metastasis sites')

    # Add counts to labels.
    ax.set_yticklabels(['{} ({})'.format(k, v) for k, v in
                        metastases.sum(axis=1).items()][::-1])

    return ax
