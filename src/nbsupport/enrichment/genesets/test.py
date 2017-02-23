# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

from itertools import chain

import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


def test_sets(selected_genes,
              gene_sets,
              all_genes=None,
              corr_method='fdr_bh',
              include_overlap=False):
    """Tests a list of selected genes for enriched gene sets.

    Args:
        selected_genes (list[str]): List of gene identifiers to test for
            enrichment. IDs should correspond with the ids used in gene_sets.
        gene_sets (dict[str, set[str]]): Dict of gene sets, with keys
            representing the names of the set and the set values containing
            the corresponding genes.
        all_genes (set[str]): Set representing the universe of all genes.
            This is typically the union of genes over all gene sets (which
            is used as the default if not given).
        corr_method (str): The method used to calculate p-values that are
            corrected for multiple testing.
        include_overlap (bool): Whether to include an extra column detailing
            the genes that overlapped with the given set.

    Returns:
        pandas.DataFrame: DataFrame with the names and (corrected) p-values
        for each of the tested sets.

    """

    # Use union of sets as default for all genes.
    if all_genes is None:
        all_genes = set(chain.from_iterable(gene_sets.values()))

    # Reduce selected genes to those within all_genes.
    selected_genes = all_genes.intersection(selected_genes)

    # Perform tests.
    result = pd.DataFrame.from_records(
        ((set_name, test_set(selected_genes, set_genes, all_genes))
         for set_name, set_genes in gene_sets.items()),
        columns=['gene_set', 'p_value'])

    # Add corrected p-value.
    result['p_value_corr'] = multipletests(
        result['p_value'], method=corr_method)[1]

    # Add genes in an extra column.
    if include_overlap:
        result['overlap'] = [
            ', '.join(sorted(selected_genes & gene_sets[set_name]))
            for set_name in result['gene_set']
        ]

    # Sort by corrected p-value.
    result.sort_values(by='p_value_corr', inplace=True, ascending=True)

    return result


def test_set(selected_genes, set_genes, all_genes):
    """Tests a list of genes for enrichment w.r.t. the given geneset.

    Args:
        selected_genes (list[str]): List of gene identifiers to test for
            enrichment. IDs should correspond with the ids used in gene_sets.
        gene_sets (dict[str, set[str]]): Dict of gene sets, with keys
            representing the names of the set and the set values containing
            the corresponding genes.
        all_genes (set[str]): Set representing the universe of all genes.
            This is typically the union of genes over all gene sets (which
            is used as the default if not given).

    Returns:
        float: A p-value between 0.0 and 1.0, reflecting the significance
            of the enrichment of the given gene set. The p-value is calculcated
            using the fishers exact test.

    """

    # Reduce the gene_set to the universe of all_genes,
    # as we can only sample from this set.
    all_genes = all_genes
    set_genes = set_genes & all_genes
    selected_genes = selected_genes & all_genes

    # Calculate overlap of the gene_set with selected genes.
    selected_set_genes = set_genes & selected_genes

    # Calculate p-value using the hyper-geometric test.
    p_val = hypergeom.sf(M=len(all_genes),
                         n=len(set_genes),
                         N=len(selected_genes),
                         k=len(selected_set_genes),
                         loc=1)

    return p_val
