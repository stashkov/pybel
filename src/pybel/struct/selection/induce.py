# -*- coding: utf-8 -*-

from .. import pipeline

__all__ = [
    'get_subgraph_by_induction',
]


@pipeline.mutator
def get_subgraph_by_induction(graph, nodes):
    """Induces a graph over the given nodes. Returns None if none of the nodes are in the given graph.

    :param pybel.BELGraph graph: A BEL graph
    :param iter[BaseEntity] nodes: A list of BEL nodes in the graph
    :rtype: Optional[pybel.BELGraph]
    """
    if all(node not in graph for node in nodes):
        return

    return graph.subgraph(nodes)
