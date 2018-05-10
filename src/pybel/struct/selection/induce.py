# -*- coding: utf-8 -*-

from ..filters import filter_edges
from ..filters.edge_predicates import is_causal_relation
from ..graph import BELGraph
from ..pipeline import mutator

__all__ = [
    'get_subgraph_by_induction',
    'get_subgraph_by_edge_filter',
    'get_causal_subgraph',
]


@mutator
def get_subgraph_by_induction(graph, nodes):
    """Induces a graph over the given nodes. Returns None if none of the nodes are in the given graph.

    :param pybel.BELGraph graph: A BEL graph
    :param iter[BaseEntity] nodes: A list of BEL nodes in the graph
    :rtype: Optional[pybel.BELGraph]
    """
    if all(node not in graph for node in nodes):
        return

    return graph.subgraph(nodes)


@mutator
def get_subgraph_by_edge_filter(graph, edge_predicates=None):
    """Induces a subgraph on all edges that pass the given filters

    :param pybel.BELGraph graph: A BEL graph
    :param edge_predicates: A predicate or list of predicates (graph, node, node, key, data) -> bool
    :type edge_predicates: (pybel.BELGraph, tuple, tuple, int) -> bool or list[(pybel.BELGraph, tuple, tuple, int) -> bool]
    :return: A BEL subgraph induced over the edges passing the given filters
    :rtype: pybel.BELGraph
    """
    result = BELGraph()

    for u, v, key in filter_edges(graph, edge_predicates=edge_predicates):
        result.add_edge(u, v, key=key, **graph.edge[u][v][key])

    return result


@mutator
def get_causal_subgraph(graph):
    """Builds a new subgraph induced over all edges that are causal

    :param pybel.BELGraph graph: A BEL graph
    :return: A subgraph of the original BEL graph
    :rtype: pybel.BELGraph
    """
    return get_subgraph_by_edge_filter(graph, is_causal_relation)
