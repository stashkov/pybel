# -*- coding: utf-8 -*-

from collections import defaultdict

from .collapse import collapse_nodes
from ...pipeline import in_place_mutator
from ....constants import RELATION, TRANSCRIBED_TO, TRANSLATED_TO

__all__ = [
    'build_central_dogma_collapse_dict',
    'collapse_by_central_dogma',
]


def build_central_dogma_collapse_dict(graph):
    """Builds a dictionary to direct the collapsing on the central dogma

    :param pybel.BELGraph graph: A BEL graph
    :return: A dictionary of {node: set of nodes}
    :rtype: dict[tuple,set[tuple]]
    """
    collapse_dict = defaultdict(set)
    r2p = {}

    for rna_node, protein_node, d in graph.edges_iter(data=True):
        if d[RELATION] != TRANSLATED_TO:
            continue

        collapse_dict[protein_node].add(rna_node)
        r2p[rna_node] = protein_node

    for gene_node, rna_node, d in graph.edges_iter(data=True):
        if d[RELATION] != TRANSCRIBED_TO:
            continue

        if rna_node in r2p:
            collapse_dict[r2p[rna_node]].add(gene_node)
        else:
            collapse_dict[rna_node].add(gene_node)

    return collapse_dict


@in_place_mutator
def collapse_by_central_dogma(graph):
    """Collapses all nodes from the central dogma (GENE, RNA, PROTEIN) to PROTEIN, or most downstream possible entity,
    in place. This function wraps :func:`collapse_nodes` and :func:`build_central_dogma_collapse_dict`.

    :param pybel.BELGraph graph: A BEL graph

    Equivalent to:

    >>> collapse_nodes(graph, build_central_dogma_collapse_dict(graph))
    """
    collapse_dict = build_central_dogma_collapse_dict(graph)
    collapse_nodes(graph, collapse_dict)
