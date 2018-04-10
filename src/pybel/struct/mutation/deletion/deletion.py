# -*- coding: utf-8 -*-

from ... import pipeline
from ...filters.edge_filters import filter_edges
from ...filters.edge_predicates import is_associative_relation
from ...filters.node_filters import get_nodes
from ...filters.node_predicate_builders import function_inclusion_filter_builder
from ....constants import PATHOLOGY

__all__ = [
    'remove_filtered_edges',
    'remove_associations',
    'remove_pathologies',
]


@pipeline.in_place_mutator
def remove_filtered_edges(graph, edge_filters):
    """Removes edges passing the given edge filters

    :param pybel.BELGraph graph: A BEL graph
    :param edge_filters: An edge filter or list of edge filters (graph, node, node, key, data)-> bool
    :type edge_filters: types.FunctionType or iter[types.FunctionType]
    :return:
    """
    edges = list(filter_edges(graph, edge_filters))
    graph.remove_edges_from(edges)


@pipeline.in_place_mutator
def remove_pathologies(graph):
    """Remove pathology nodes

    :param pybel.BELGraph graph: A BEL graph
    """
    nodes = get_nodes(graph, function_inclusion_filter_builder(PATHOLOGY))
    graph.remove_nodes_from(nodes)


@pipeline.in_place_mutator
def remove_associations(graph):
    """Removes all associative relationships from the graph

    :param pybel.BELGraph graph: A BEL graph
    """
    remove_filtered_edges(graph, is_associative_relation)
