# -*- coding: utf-8 -*-

from .subgraph import get_subgraph_by_edge_filter
from ...filters.edge_predicates import is_causal_relation
from ...pipeline import mutator

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
def get_causal_subgraph(graph):
    """Builds a new subgraph induced over all edges that are causal

    :param pybel.BELGraph graph: A BEL graph
    :return: A subgraph of the original BEL graph
    :rtype: pybel.BELGraph
    """
    return get_subgraph_by_edge_filter(graph, is_causal_relation)
