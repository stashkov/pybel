# -*- coding: utf-8 -*-

from ...pipeline import in_place_mutator
from ....constants import RELATION, unqualified_edges

__all__ = [
    'collapse_pair',
    'collapse_nodes',
]


@in_place_mutator
def collapse_pair(graph, *, survivor, from_node):
    """Rewires all edges from the synonymous node to the survivor node, then deletes the synonymous node.

    Does not keep edges between the two nodes.

    :param pybel.BELGraph graph: A BEL graph
    :param BaseEntity survivor: The BEL node to collapse all edges on the synonym to
    :param BaseEntity from_node: The BEL node to collapse into the surviving node
    """

    for _, successor, key, data in graph.out_edges(from_node, keys=True, data=True):
        if successor == survivor:
            continue

        if data[RELATION] in unqualified_edges:
            graph.add_unqualified_edge(survivor, successor, data[RELATION])
        else:
            graph.add_edge(survivor, successor, key=key, **data)

    for predecessor, _, key, data in graph.in_edges(from_node, keys=True, data=True):
        if predecessor == survivor:
            continue

        if data[RELATION] in unqualified_edges:
            graph.add_unqualified_edge(predecessor, survivor, data[RELATION])
        else:
            graph.add_edge(predecessor, survivor, key=key, **data)

    graph.remove_node(from_node)


@in_place_mutator
def collapse_nodes(graph, dict_of_sets_of_nodes):
    """Collapses all nodes in values to the key nodes, in place

    :param pybel.BELGraph graph: A BEL graph
    :param dict[BaseEntity,set[BaseEntity]] dict_of_sets_of_nodes: A dictionary of {node: set of nodes}
    """
    for key_node, value_nodes in dict_of_sets_of_nodes.items():
        for value_node in value_nodes:
            collapse_pair(graph, from_node=value_node, survivor=key_node)

    # Remove self edges
    for u, v, k in graph.edges(keys=True):
        if u == v:
            graph.remove_edge(u, v, k)
