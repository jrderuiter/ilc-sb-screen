import pandas as pd
from matplotlib import pyplot as plt
from scipy.spatial import distance
from scipy.cluster import hierarchy
import seaborn as sns

from nbsupport.util.clustermap import color_annotation
from nbsupport.util.seaborn import clustermap

from .combat import combat


def plot_pam50(combined_expr_combat,
               combined_design,
               genes,
               subset=None,
               metric='euclidean',
               method='average',
               extra_annotation=None,
               **kwargs):

    # Perform batch correction with ComBat.
    # combined_expr_combat = combat(
    #     combined_expr, batch=combined_design['organism'])

    # Prepare annotation.
    pam50_order = ['LumA', 'LumB', 'Her2', 'Basal', 'Normal']
    model_order = ['Basal-like', 'Luminal', 'SB']

    annotation = pd.DataFrame(
        {
            'Mouse model': pd.Categorical(
                combined_design['mouse_model'], categories=model_order),
            'PAM 50 (Human)': pd.Categorical(
                combined_design['subtype'], categories=pam50_order)
        },
        index=combined_design.index,
        columns=['PAM 50 (Human)', 'Mouse model'])

    # Color annotation.
    model_palette = sns.color_palette(
        sns.xkcd_palette(['greyish', 'midnight', 'tangerine']))
    pam50_palette = sns.color_palette(
        ['#a8cee2', '#2777b1', '#f89998', '#e02025', '#87ae73'])

    colored_annotation, _ = color_annotation(
        annotation, colors=[pam50_palette, model_palette])

    if extra_annotation is not None:
        colored_annotation = pd.concat(
            [colored_annotation, extra_annotation], axis=1)

    # Calculate column clustering.
    col_linkage = _calc_linkage(
        combined_expr_combat.ix[genes],
        z_score=0,
        metric=metric,
        method=method)

    # Draw heatmap.
    fig = clustermap(
        combined_expr_combat.ix[subset or genes],
        z_score=0,
        col_linkage=col_linkage,
        method=method,
        col_colors=colored_annotation,
        **kwargs)

    plt.setp(fig.ax_heatmap.get_yticklabels(), rotation=0, fontstyle='italic')
    fig.ax_heatmap.set_xticklabels([])

    return fig


def _calc_zscore(data2d, axis=1):
    """Standarize the mean and variance of the data axis."""
    other_axis = 0 if axis == 1 else 1
    return (data2d.subtract(
        data2d.mean(axis=other_axis), axis=axis).divide(
            data2d.std(axis=other_axis), axis=axis))


def _calc_linkage(data2d,
                  axis=0,
                  z_score=None,
                  metric='euclidean',
                  method='complete'):
    if z_score is not None:
        data2d = _calc_zscore(data2d, axis=z_score)

    if axis == 0:
        data2d = data2d.T

    dist = distance.pdist(data2d, metric=metric)
    linkage = hierarchy.linkage(dist, method=method)

    return linkage
