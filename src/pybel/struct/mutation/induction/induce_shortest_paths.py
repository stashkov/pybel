# -*- coding: utf-8 -*-

import itertools as itt
import logging
from itertools import product

import networkx as nx
from networkx import all_shortest_paths

from .induce import get_subgraph_by_induction
from ...pipeline import mutator
from ....constants import FUNCTION, PATHOLOGY

log = logging.getLogger(__name__)

__all__ = [
    'get_nodes_in_all_shortest_paths',
    'get_subgraph_by_all_shortest_paths',
]


def _get_nodes_in_all_shortest_paths_helper(graph, nodes, weight=None, remove_pathologies=True):
    if remove_pathologies:
        graph = graph.copy()
        for node, data in graph.nodes(data=True):
            if data[FUNCTION] == PATHOLOGY:
                graph.remove_node(node)

    for u, v in product(nodes, repeat=2):
        try:
            yield from all_shortest_paths(graph, u, v, weight=weight)
        except nx.exception.NetworkXNoPath:
            continue


def get_nodes_in_all_shortest_paths(graph, nodes, weight=None, remove_pathologies=True):
    """Gets all shortest paths from all nodes to all other nodes in the given list and returns the set of all nodes
    contained in those paths using :func:`networkx.all_shortest_paths`.

    :param pybel.BELGraph graph: A BEL graph
    :param iter[tuple] nodes: The list of nodes to use to use to find all shortest paths
    :param str weight: Edge data key corresponding to the edge weight. If none, uses unweighted search.
    :param bool remove_pathologies: Should pathology nodes be removed first?
    :return: A set of nodes appearing in the shortest paths between nodes in the BEL graph
    :rtype: set[tuple]

    .. note:: This can be trivially parallelized using :func:`networkx.single_source_shortest_path`
    """
    shortest_paths_nodes_iterator = _get_nodes_in_all_shortest_paths_helper(graph, nodes, weight=weight,
                                                                            remove_pathologies=remove_pathologies)

    return set(itt.chain.from_iterable(shortest_paths_nodes_iterator))


@mutator
def get_subgraph_by_all_shortest_paths(graph, nodes, weight=None, remove_pathologies=True):
    """Induces a subgraph over the nodes in the pairwise shortest paths between all of the nodes in the given list

    :param pybel.BELGraph graph: A BEL graph
    :param set[tuple] nodes: A set of nodes over which to calculate shortest paths
    :param str weight: Edge data key corresponding to the edge weight. If None, performs unweighted search
    :param bool remove_pathologies: Should the pathology nodes be deleted before getting shortest paths?
    :return: A BEL graph induced over the nodes appearing in the shortest paths between the given nodes
    :rtype: Optional[pybel.BELGraph]
    """
    query_nodes = []

    for node in nodes:
        if node not in graph:
            log.debug('%s not in %s', node, graph)
            continue
        query_nodes.append(node)

    if not query_nodes:
        return

    induced_nodes = get_nodes_in_all_shortest_paths(graph, query_nodes, weight=weight,
                                                    remove_pathologies=remove_pathologies)

    if not induced_nodes:
        return

    return get_subgraph_by_induction(graph, induced_nodes)
