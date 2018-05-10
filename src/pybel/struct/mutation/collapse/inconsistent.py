# -*- coding: utf-8 -*-

from ...pipeline import in_place_mutator
from ....constants import RELATION

__all__ = [
    'remove_inconsistent_edges',
    'collapse_consistent_edges',
]


def all_edges_iter(graph, u, v):
    """Lists all edges between the given nodes

    :param pybel.BELGraph graph: A BEL Graph
    :param BaseEntity u: A BEL node
    :param BaseEntity v: A BEL node
    :return: Yields tuples of (node, node, key)
    :rtype: iter[tuple[BaseEntity,BaseEntity,str]]
    """
    if u not in graph or v not in graph[u]:
        raise ValueError(f'Graph has no edges between {u} and {v}')

    for k in graph[u][v].keys():
        yield u, v, k


def get_all_relations(graph, u, v):
    """Returns the set of all relations between a given pair of nodes

    :param pybel.BELGraph graph: A BEL graph
    :param BaseEntity u: The source BEL node
    :param BaseEntity v: The target BEL node
    :rtype: set[str]
    """
    if v not in graph[u]:
        return set()

    return {
        data[RELATION]
        for data in graph[u][v].values()
    }


def pair_is_consistent(graph, u, v):
    """Returns if the edges between the given nodes are consistent, meaning they all have the same relation

    :param pybel.BELGraph graph: A BEL graph
    :param tuple u: The source BEL node
    :param tuple v: The target BEL node
    :return: If the edges aren't consistent, return false, otherwise return the relation type
    :rtype: Optional[str]
    """
    relations = get_all_relations(graph, u, v)

    if 1 != len(relations):
        return

    return list(relations)[0]


def _iter_pairs(graph):
    """Iterates over unique node-node pairs in the graph

    :param pybel.BELGraph graph: A BEL graph
    :rtype: iter
    """
    for u, v in set(graph.edges()):
        yield u, v


def get_inconsistent_edges(graph):
    """Returns an iterator over inconsistent edges

    :param pybel.BELGraph graph: A BEL graph
    :return: An iterator over (source, target) node pairs corresponding to edges with many inconsistent relations
    :rtype: iter
    """
    for u, v in _iter_pairs(graph):
        if not pair_is_consistent(graph, u, v):
            yield u, v


def get_consistent_edges(graph):
    for u, v in _iter_pairs(graph):
        relation = pair_is_consistent(graph, u, v)
        if relation:
            yield u, v, relation


@in_place_mutator
def remove_inconsistent_edges(graph):
    """Remove all edges between node paris with consistent edges.

    This is the all-or-nothing approach. It would be better to do more careful investigation of the evidences during
    curation.

    :param pybel.BELGraph graph: A BEL graph
    """
    for u, v in get_inconsistent_edges(graph):
        edges = list(all_edges_iter(graph, u, v))
        graph.remove_edges_from(edges)


@in_place_mutator
def collapse_consistent_edges(graph):
    """Collapses consistent edges together

    .. warning:: This operation doesn't preserve evidences or other annotations

    :param pybel.BELGraph graph: A BEL Graph
    """
    for u, v in graph.edges():
        if not pair_is_consistent(graph, u, v):
            continue

        edges = list(all_edges_iter(graph, u, v))

        if len(edges) == 1:
            continue

        graph.remove_edges_from(edges[1:])
