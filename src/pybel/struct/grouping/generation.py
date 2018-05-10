# -*- coding: utf-8 -*-

"""Unbiased candidate mechanism generation.

The Unbiased Candidate Mechanism Generation workflow addresses the inconsistency in the definitions of the
boundaries of pathways, mechanisms, subgraphs, etc. in networks and systems biology that are introduced during curation
due to a variety of reasons.

A simple approach for generating unbiased candidate mechanisms is to take the upstream controlles

This module provides functions for generating subgraphs based around a single node, most likely a biological process.

Subgraphs induced around biological processes should prove to be subgraphs of the NeuroMMSig/canonical mechanisms
and provide an even more rich mechanism inventory.

Examples
~~~~~~~~
This method has been applied in the following Jupyter Notebooks:

- `Generating Unbiased Candidate Mechanisms <http://nbviewer.jupyter.org/github/pybel/pybel-notebooks/blob/master/algorithms/Generating%20Candidate%20Mechanisms.ipynb>`_
"""

from ..filters import data_missing_key_builder, filter_nodes, get_nodes_by_function, is_upstream_leaf
from ..mutation import expand_upstream_causal_subgraph, get_upstream_causal_subgraph
from ..mutation.collapse import collapse_consistent_edges, remove_inconsistent_edges
from ..pipeline import in_place_mutator, mutator, splitter
from ...constants import BIOPROCESS

__all__ = [
    'generate_mechanism',
    'generate_bioprocess_mechanisms',
]


def get_unweighted_upstream_leaves(graph, key):
    """Gets all leaves of the graph with no incoming edges, one outgoing edge, and without the given key in
    its data dictionary

    .. seealso :: :func:`data_does_not_contain_key_builder`

    :param pybel.BELGraph graph: A BEL graph
    :param str key: The key in the node data dictionary representing the experimental data
    :return: An iterable over leaves (nodes with an in-degree of 0) that don't have the given annotation
    :rtype: iter[tuple]
    """
    return filter_nodes(graph, [is_upstream_leaf, data_missing_key_builder(key)])


@in_place_mutator
def remove_unweighted_leaves(graph, key):
    """

    :param pybel.BELGraph graph: A BEL graph
    :param str key: The key in the node data dictionary representing the experimental data
    """
    unweighted_leaves = list(get_unweighted_upstream_leaves(graph, key))
    graph.remove_nodes_from(unweighted_leaves)


def is_unweighted_source(graph, node, key):
    """
    
    :param pybel.BELGraph graph: A BEL graph
    :param tuple node: A BEL node
    :param str key: The key in the node data dictionary representing the experimental data
    """
    return graph.in_degree(node) == 0 and key not in graph.node[node]


def get_unweighted_sources(graph, key):
    """Gets unannotated nodes on the periphery of the subgraph

    :param pybel.BELGraph graph: A BEL graph
    :param str key: The key in the node data dictionary representing the experimental data
    :return: An iterator over BEL nodes that are unannotated and on the periphery of this subgraph
    :rtype: iter[BaseEntity]
    """
    for node in graph:
        if is_unweighted_source(graph, node, key):
            yield node


@in_place_mutator
def remove_unweighted_sources(graph, key):
    """Prunes unannotated nodes on the periphery of the subgraph

    :param pybel.BELGraph graph: A BEL graph
    :param str key: The key in the node data dictionary representing the experimental data
    """
    nodes = list(get_unweighted_sources(graph, key))
    graph.remove_nodes_from(nodes)


@in_place_mutator
def prune_mechanism_by_data(graph, key):
    """Removes all leaves and source nodes that don't have weights. Is a thin wrapper around 
    :func:`remove_unweighted_leaves` and :func:`remove_unweighted_sources`

    :param pybel.BELGraph graph: A BEL Graph
    :param str key: The key in the node data dictionary representing the experimental data. If none, does not prune
                    unannotated nodes after generation

    Equivalent to:
    
    >>> remove_unweighted_leaves(graph, key)
    >>> remove_unweighted_sources(graph, key)
    """
    remove_unweighted_leaves(graph, key)
    remove_unweighted_sources(graph, key)


@mutator
def generate_mechanism(graph, node, key=None, prune_inconsistent=True):
    """Generates a mechanistic subgraph upstream of the given node

    :param pybel.BELGraph graph: A BEL Graph
    :param BaseEntity node: The target BEL node for generation
    :param str key: The key in the node data dictionary representing the experimental data. If none, does not prune
                    unannotated nodes after generation
    :return: A subgraph grown around the target BEL node
    :rtype: pybel.BELGraph
    """
    subgraph = get_upstream_causal_subgraph(graph, node)
    expand_upstream_causal_subgraph(graph, subgraph)

    if prune_inconsistent:
        remove_inconsistent_edges(subgraph)
        collapse_consistent_edges(subgraph)

    if key is not None:
        prune_mechanism_by_data(subgraph, key)

    return subgraph


@splitter
def generate_bioprocess_mechanisms(graph, key=None):
    """Generates a mechanistic subgraph for each biological process in the graph using :func:`generate_mechanism`

    :param pybel.BELGraph graph: A BEL Graph
    :param str key: The key in the node data dictionary representing the experimental data. If none, does not prune
                unannotated nodes after generation
    :return: A dictionary from {tuple bioprocess node: BELGraph candidate mechanism}
    :rtype: dict[bioprocess,pybel.BELGraph]
    """
    return {
        bp: generate_mechanism(graph, bp, key=key)
        for bp in get_nodes_by_function(graph, BIOPROCESS)
    }
