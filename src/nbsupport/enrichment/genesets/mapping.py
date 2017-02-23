# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *

# pylint: enable=W0622,W0614,W0401


def map_sets(gene_sets, from_type, to_type, mapper, **kwargs):
    """Maps genesets to different ids/organisms.

    Args:
        gene_sets (dict[str, set[str]]): Dict of gene sets to be mapped.
            Keys should represent names of gene sets, whilst values
            should be sets of gene ids.
        from_type (str): The type of id used in the given gene sets.
        to_type (str): The type of id the sets should be mapped to.
        mapper (str): The name of the mapper used to map gene sets.
        **kwargs: Any kwargs are passed to the used mapper.

    Returns:
        dict[str, set[str]]: Dict of gene sets with mapped ids.

    """

    import genemap

    # Fetch mapping and convert to dict.
    mapping = genemap.get_map(
        from_type=from_type, to_type=to_type, mapper=mapper, **kwargs)
    mapping_dict = dict(zip(mapping.ix[:, 0], mapping.ix[:, 1]))

    # Map genesets using mapping.
    mapped_sets = {k: {mapping_dict[g]
                       for g in v if g in mapping_dict}
                   for k, v in gene_sets.items()}

    return mapped_sets
