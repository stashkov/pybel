# -*- coding: utf-8 -*-

import logging

from ...filters.edge_filters import filter_edges
from ...filters.edge_predicate_builders import (
    build_annotation_dict_all_filter, build_annotation_dict_any_filter,
)
from ...pipeline import mutator

__all__ = [
    'get_subgraph_by_edge_filter',
    'get_subgraph_by_annotations',
    'get_subgraph_by_annotation_value',
]

log = logging.getLogger(__name__)


@mutator
def get_subgraph_by_edge_filter(graph, edge_filters):
    """Induces a subgraph on all edges that pass the given filters

    :param pybel.BELGraph graph: A BEL graph
    :param edge_filters: A predicate or list of predicates (graph, node, node, key, data) -> bool
    :type edge_filters: (pybel.BELGraph, tuple, tuple, int) -> bool or list[(pybel.BELGraph, tuple, tuple, int) -> bool]
    :return: A BEL subgraph induced over the edges passing the given filters
    :rtype: pybel.BELGraph
    """
    result = graph.fresh_copy()

    for u, v, k in filter_edges(graph, edge_filters):
        result.add_entity(u)
        result.add_entity(v)
        result.add_edge(u, v, key=k, **graph[u][v][k])

    return result


@mutator
def get_subgraph_by_annotations(graph, annotations, or_=None):
    """Returns the subgraph given an annotations filter.

    :param graph: pybel.BELGraph graph: A BEL graph
    :param dict[str,set[str]] annotations: Annotation filters (match all with :func:`pybel.utils.subdict_matches`)
    :param boolean or_: if True any annotation should be present, if False all annotations should be present in the
                        edge. Defaults to True.
    :return: A subgraph of the original BEL graph
    :rtype: pybel.BELGraph
    """
    edge_filter_builder = (
        build_annotation_dict_any_filter
        if (or_ is None or or_) else build_annotation_dict_all_filter
    )

    return get_subgraph_by_edge_filter(graph, edge_filter_builder(annotations))


@mutator
def get_subgraph_by_annotation_value(graph, annotation, value):
    """Builds a new subgraph induced over all edges whose annotations match the given key and value

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation: The annotation to group by
    :param str value: The value for the annotation
    :return: A subgraph of the original BEL graph
    :rtype: pybel.BELGraph
    """
    return get_subgraph_by_annotations(graph, {annotation: {value}})


@mutator
def get_subgraph_by_neighborhood(graph, nodes):
    """Gets a BEL graph around the neighborhoods of the given nodes. Returns none if no nodes are in the graph

    :param pybel.BELGraph graph: A BEL graph
    :param iter[BaseEntity] nodes: An iterable of BEL nodes
    :return: A BEL graph induced around the neighborhoods of the given nodes
    :rtype: Optional[pybel.BELGraph]
    """
    result = graph.fresh_copy()

    node_set = set(nodes)

    if all(node not in graph for node in node_set):
        return

    result.add_edges_from(graph.in_edges(nodes, keys=True, data=True))
    result.add_edges_from(graph.out_edges(nodes, keys=True, data=True))

    return result


