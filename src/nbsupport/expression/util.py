import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.spatial import distance
from scipy.cluster import hierarchy
import seaborn as sns
import toolz


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
    size_factors = np.apply_along_axis(
        _estimate_size_factors_col,
        axis=0,
        arr=counts,
        log_geo_means=log_geo_means)
    return size_factors


def _estimate_size_factors_col(counts, log_geo_means):
    with np.errstate(divide='ignore'):
        log_counts = np.log(counts)
    mask = np.isfinite(log_geo_means) & (counts > 0)
    return np.exp(np.median((log_counts - log_geo_means)[mask]))


def plot_count_distribution(expr, is_log=True, ax=None):
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


def tidy_expression(expr, design=None):
    """Converts expression matrix into a tidy 'long' format."""

    df_long = pd.melt(
        _reset_index(
            expr, name='gene'), id_vars=['gene'], var_name='sample')

    if design is not None:
        df_long = pd.merge(
            df_long,
            _reset_index(
                design, name='sample'),
            on='sample',
            how='left')

    return df_long


def _reset_index(df, name=None):
    """Resets the index of the dataframe, optionally renaming the index."""

    if name is None:
        return df.reset_index()
    else:
        index_name = df.index.name or 'index'
        return df.reset_index().rename(columns={index_name: name})