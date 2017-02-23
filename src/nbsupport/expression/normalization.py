import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()

r_base = None
r_sva = None


def combat(expr, batch):
    """Normalizes expression data for batch effects using ComBat.

    Parameters
    ----------
    expr : pandas.DataFrame
        Dataframe of samples-by-genes containing expression values.

    """

    assert expr.shape[1] == batch.shape[0]

    global r_base
    global r_sva

    # Perform imports if needed.
    if r_base is None:
        r_base = importr('base')

    if r_sva is None:
        r_sva = importr('sva')

    # Filter for minimum variance.
    expr_filt = expr.ix[expr.var(axis=1) > 1e-12]

    # Run ComBat.
    expr_combat = r_sva.ComBat(
        r_base.as_matrix(expr_filt), batch=r_base.as_character(batch))

    # Convert back to frame and return.
    return pd.DataFrame(
        np.array(expr_combat),
        index=expr_filt.index,
        columns=expr_filt.columns)


def normalize_counts(counts):
    """Normalizes expression counts using DESeq's median-of-ratios approach.

    Parameters
    ----------
    counts : pandas.DataFrame
        Dataframe of genes-by-samples containing raw expression counts.

    Returns
    -------
    pandas.DataFrame
        Dataframe containing expression values that are normalized for
        differences in sequencing depth between samples.

    """

    size_factors = estimate_size_factors(counts)
    return counts.divide(size_factors, axis='columns')


def estimate_size_factors(counts):
    """Estimates size factors for median-of-ratios count normalization.

    Parameters
    ----------
    counts : pandas.DataFrame
        Dataframe of genes-by-samples containing raw expression counts.

    Returns
    -------


    """
    with np.errstate(divide='ignore'):
        log_geo_means = np.mean(np.log(counts), axis=1)
    sf = np.apply_along_axis(
        _estimate_size_factors_col,
        axis=0,
        arr=counts,
        log_geo_means=log_geo_means)
    return sf


def _estimate_size_factors_col(counts, log_geo_means):
    with np.errstate(divide='ignore'):
        log_counts = np.log(counts)
    mask = np.isfinite(log_geo_means) & (counts > 0)
    return np.exp(np.median((log_counts - log_geo_means)[mask]))


def plot_distribution(expr, is_log=True, ax=None):
    """Plots a simple expression value distribution."""

    if ax is None:
        _, ax = plt.subplots()

    if not is_log:
        expr = np.log2(expr + 1)

    for _, values in expr.items():
        sns.distplot(values, hist=False, ax=ax)

    ax.set_xlabel('Expression (log2)')
    ax.set_ylabel('Density')
    sns.despine(ax=ax)

    return ax
