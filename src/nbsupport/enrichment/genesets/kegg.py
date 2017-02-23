# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import io
import re

import pandas as pd
import requests
import requests_cache

KEGG_HOST = 'http://rest.kegg.jp'

KEGG_PATHWAY_URL = '/list/pathway/{species}'
KEGG_GENE_URL = '/list/{species}'
KEGG_MAPPING_URL = '/link/pathway/{species}'
KEGG_CONV_URL = '/conv/ncbi-geneid/{species}'

requests_cache.install_cache('_kegg_cache')


def get_gene_sets(species, use_kegg_ids=False, use_name=True):
    """Gets mapping from gene ids to pathways as genesets.

    Args:
        species (str): Name of the species to query. Use or example 'hsa'
            for human or 'mmu' for mouse.
        use_kegg_ids (bool): Whether to return gene ids as entrez
            ids (False), or as KEGG gene ids (True).
        use_name (bool): Whether to use the names of the pathways in the
            geneset dict (True), or the pathway ids (False).

    Returns:
        dict of sets: Dict mapping pathways (the keys) to gene ids (the sets).

    """
    # Fetch mapping frame.
    pathway_map = get_pathway_mapping(
        species, use_kegg_ids=use_kegg_ids, use_name=use_name)

    # Group into genesets.
    groupby_key = 'description' if use_name else 'pathway_id'
    genesets = {key: set(grp['gene_id'])
                for key, grp in pathway_map.groupby(groupby_key)}

    return genesets


def get_genes(species):
    """Gets KEGG genes and their descriptions."""
    response = _get(KEGG_HOST + KEGG_GENE_URL.format(species=species))
    return _read_response(response, sep='\t', names=['gene_id', 'description'])


def get_pathways(species):
    """Gets KEGG pathways and their descriptions."""
    response = _get(KEGG_HOST + KEGG_PATHWAY_URL.format(species=species))
    frame = _read_response(response, names=['pathway_id', 'description'])
    frame['description'] = _remove_species(frame['description'])
    return frame


def get_pathway_mapping(species, use_kegg_ids=False, use_name=False):
    """Gets the mapping from gene ids to KEGG pathways.

    Args:
        species (str): Name of the species to query. Use or example 'hsa'
            for human or 'mmu' for mouse.
        use_kegg_ids (bool): Whether to return gene ids as entrez
            ids (False), or as KEGG gene ids (True).
        use_name (bool): Whether to include the name of the KEGG
            pathway as an extra column in the DataFrame.

    Returns:
        pd.DataFrame: Pandas DataFrame containing 'gene_id' and
            'pathway_id' columns for mapping genes to KEGG pathways.
            A 'description' column is included if with_name is True.

    """

    response = _get(KEGG_HOST + KEGG_MAPPING_URL.format(species=species))
    frame = _read_response(response, names=['gene_id', 'pathway_id'])

    if not use_kegg_ids:
        # Convert gene ids to (numeric) entrez identifiers.
        entrez_map = _get_entrez_mapping(species)
        entrez_map = {kegg: entrez
                      for kegg, entrez in zip(entrez_map['kegg_id'],
                                              entrez_map['entrez_id'])}
        frame['gene_id'] = frame['gene_id'].map(entrez_map)

    if use_name:
        # Add pathway description to mapping.
        pathways = get_pathways(species)
        frame = pd.merge(frame, pathways, on='pathway_id')

    return frame


def _get_entrez_mapping(species):
    """Gets the mapping of KEGG gene ids to entrez gene ids."""
    response = _get(KEGG_HOST + KEGG_CONV_URL.format(species=species))
    frame = _read_response(response, names=['kegg_id', 'entrez_id'])

    frame['entrez_id'] = frame['entrez_id'].str.replace('ncbi-geneid:', '')
    frame['entrez_id'] = frame['entrez_id'].astype(str)

    return frame


def _get(url):
    """Internal function for getting requests."""
    r = requests.get(url)
    r.raise_for_status()
    return r


def _read_response(response, sep='\t', header=None, **kwargs):
    """Reads response content into a DataFrame using pd.read_csv."""
    bytes_obj = io.BytesIO(response.content)
    return pd.read_csv(bytes_obj, sep=sep, header=header, **kwargs)


def _remove_species(series):
    """Removes species label from pathway description."""

    regex = re.compile(r'(\s-\s.+)')
    match = regex.search(series[0])

    if match is not None:
        species_str = match.groups()[0]
        series = series.apply(lambda x: x.replace(species_str, ''))

    return series
