import itertools

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def contingency_table(a, b):
    """Calculates a 2x2 contigency table for two boolean Series.

    Args:
        a (pandas.Series): Boolean series a.
        b (pandas.Series): Boolean series b.

    Returns:
        pd.DataFrame: A 2x2 contingency table for a and b.

    """

    if a.name == b.name:
        raise ValueError('a and b should not have the same name ({})'
                         .format(a.name))

    # Subset to common indices and convert to categoricals.
    common_indices = list(set(a.index) & set(b.index))
    a_cat = a.ix[common_indices].astype('category', categories=[True, False])
    b_cat = b.ix[common_indices].astype('category', categories=[True, False])

    # Create 2D contingency table.
    return pd.crosstab(a_cat, b_cat)


def test_association(a, b, alternative='two-sided'):
    """Tests for an association between series a and b using fishers exact.

    Args:
        a (pandas.Series): Boolean series a.
        b (pandas.Series): Boolean series b.
        alternative (str): Alternative hypothesis to test.

    Returns:
        float: p-value from the fishers exact test, reflecting the degree
            of association between series a and b.

    """

    table = contingency_table(a, b)
    return fisher_exact(table, alternative=alternative)


def test_associations(a,
                      b,
                      alternative='two-sided',
                      corr_method='fdr_bh',
                      labels=('a', 'b')):
    """Tests for significant associations between columns in dataframes a and b.

    Args:
        a (pandas.DataFrame): DataFrame a. Should either be boolean or contain
            categorical values that can be expanded to a boolean matrix.
        b (pandas.DataFrame): DataFrame b, same format as a.
        alternative (str): Alternative hypothesis to test.
        corr_method (str): Correction method to use for to correct p-values
            for multiple testing.
        labels (str): Labels to use to in result.

    Returns:
        pd.DataFrame: Result dataframe containing associations and
            corresponding (corrected) p-values.

    """

    # Ensure a and b are boolean.
    a = convert_to_boolean(a)
    b = convert_to_boolean(b)

    # Combinations to try.
    combinations = itertools.product(a.iteritems(), b.iteritems())

    # Generate result frame from combinations.
    rows = ((a, b, test_association(
        a_values, b_values, alternative=alternative)[1])
            for (a, a_values), (b, b_values) in combinations)

    col_names = labels + ('p_value', )
    result = pd.DataFrame.from_records(rows, columns=col_names)

    # Apply multiple testing correction.
    result['p_value_corr'] = multipletests(
        result['p_value'], method=corr_method)[1]

    return result


def convert_to_boolean(data):
    """Converts categorical series or dataframe to boolean."""

    if isinstance(data, pd.DataFrame):
        return _convert_frame_to_boolean(data)
    elif isinstance(data, pd.Series):
        return _convert_series_to_boolean(data)
    else:
        raise ValueError('Unknown data type ({})'.format(type(data)))


def _convert_series_to_boolean(series):
    """Converts categorical series to boolean matrix."""

    if series.dtype.name == 'bool':
        return series
    else:
        expanded = pd.pivot_table(
            series.reset_index(),
            index=series.index.name or 'index',
            columns=series.name,
            aggfunc=len)
        return expanded > 0


def _convert_frame_to_boolean(frame):
    """Converts categorical DataFrame to boolean matrix."""

    converted = []
    for col_name, col_values in frame.items():
        if col_values.dtype.name != 'bool':
            conv = convert_to_boolean(col_values)
            conv.columns = ['{}_{}'.format(col_name, c) for c in conv.columns]
            converted.append(conv)
        else:
            converted.append(col_values)

    return pd.concat(converted, axis=1)


def test_set_associations(a, b, sets_a=None, sets_b=None, **kwargs):
    """Tests for associations between a and b, after grouping by sets.

    Tests for associations between a and b, after combining columns
    in a and/or b by the given genesets. For example, if a is a matrix
    of samples-by-genes and the given sets represent sets of genes, then
    columns of a are combined (using a column-wise any) into a single
    column reflecting that gene set. This reduced matrix is then used
    to test for associations. Note that columns not in any set are dropped.

    """

    # Convert a and b to boolean matrices.
    a = convert_to_boolean(a)
    b = convert_to_boolean(b)

    # Summarize a and b into sets (based on column values).
    if sets_a is not None:
        a = _summarize_sets(a, sets_a)

    if sets_b is not None:
        b = _summarize_sets(b, sets_b)

    # Test for associations.
    return test_associations(a, b, **kwargs)


def _summarize_sets(bool_mat, sets):
    """Summarizes boolean matrix for values in sets."""

    return pd.DataFrame({set_name: bool_mat[list(set_values)].any(axis=1)
                         for set_name, set_values in sets.items()})
