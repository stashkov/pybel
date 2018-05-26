# -*- coding: utf-8 -*-

import logging

from .expansion import expand_all_node_neighborhoods, expand_nodes_neighborhoods
from .induction.induce import get_subgraph_by_induction
from .induction.induce_shortest_paths import get_subgraph_by_all_shortest_paths
from .induction.subgraph import (
    get_subgraph_by_annotations, get_subgraph_by_neighborhood,
)
from ..pipeline import mutator
from ...constants.induce_subgraph_keys import *

__all__ = [
    'get_subgraph_by_second_neighbors',
    'get_subgraph',
]

log = logging.getLogger(__name__)


@mutator
def get_subgraph_by_second_neighbors(graph, nodes, filter_pathologies=False):
    """Gets a BEL graph around the neighborhoods of the given nodes, and expands to the neighborhood of those nodes

    :param pybel.BELGraph graph: A BEL graph
    :param iter[tuple] nodes: An iterable of BEL nodes
    :param bool filter_pathologies: Should expansion take place around pathologies?
    :return: A BEL graph induced around the neighborhoods of the given nodes
    :rtype: Optional[pybel.BELGraph]
    """
    result = get_subgraph_by_neighborhood(graph, nodes)

    if result is None:
        return

    expand_all_node_neighborhoods(graph, result, filter_pathologies=filter_pathologies)
    return result


@mutator
def get_subgraph(graph, seed_method=None, seed_data=None, expand_nodes=None, remove_nodes=None):
    """Runs pipeline query on graph with multiple subgraph filters and expanders.

    Order of Operations:

    1. Seeding by given function name and data
    2. Add nodes
    3. Remove nodes

    :param pybel.BELGraph graph: A BEL graph
    :param str seed_method: The name of the get_subgraph_by_* function to use
    :param seed_data: The argument to pass to the get_subgraph function
    :param list[tuple] expand_nodes: Add the neighborhoods around all of these nodes
    :param list[tuple] remove_nodes: Remove these nodes and all of their in/out edges
    :rtype: Optional[pybel.BELGraph]
    """

    # Seed by the given function
    if seed_method == SEED_TYPE_INDUCTION:
        result = get_subgraph_by_induction(graph, seed_data)

    elif seed_method == SEED_TYPE_PATHS:
        result = get_subgraph_by_all_shortest_paths(graph, seed_data)

    elif seed_method == SEED_TYPE_NEIGHBORS:
        result = get_subgraph_by_neighborhood(graph, seed_data)

    elif seed_method == SEED_TYPE_DOUBLE_NEIGHBORS:
        result = get_subgraph_by_second_neighbors(graph, seed_data)

    elif seed_method == SEED_TYPE_UPSTREAM:
        raise NotImplementedError
        result = get_multi_causal_upstream(graph, seed_data)

    elif seed_method == SEED_TYPE_DOWNSTREAM:
        raise NotImplementedError
        result = get_multi_causal_downstream(graph, seed_data)

    elif seed_method == SEED_TYPE_PUBMED:
        raise NotImplementedError
        result = get_subgraph_by_pubmed(graph, seed_data)

    elif seed_method == SEED_TYPE_AUTHOR:
        raise NotImplementedError
        result = get_subgraph_by_authors(graph, seed_data)

    elif seed_method == SEED_TYPE_ANNOTATION:
        result = get_subgraph_by_annotations(graph, seed_data['annotations'], or_=seed_data.get('or'))

    elif seed_method == SEED_TYPE_SAMPLE:
        raise NotImplementedError
        result = get_random_subgraph(
            graph,
            number_edges=seed_data.get('number_edges'),
            seed=seed_data.get('seed')
        )

    elif not seed_method:  # Otherwise, don't seed a subgraph
        result = graph.copy()
        log.debug('no seed function - using full network: %s', result.name)

    else:
        raise ValueError('Invalid seed method: {}'.format(seed_method))

    if result is None:
        log.debug('query returned no results')
        return

    log.debug('original graph has (%s nodes / %s edges)', result.number_of_nodes(), result.number_of_edges())

    # Expand around the given nodes
    if expand_nodes:
        expand_nodes_neighborhoods(graph, result, expand_nodes)
        log.debug('graph expanded to (%s nodes / %s edges)', result.number_of_nodes(), result.number_of_edges())

    # Delete the given nodes
    if remove_nodes:
        for node in remove_nodes:
            if node not in result:
                log.debug('%s is not in graph %s', node, graph.name)
                continue
            result.remove_node(node)
        log.debug('graph contracted to (%s nodes / %s edges)', result.number_of_nodes(), result.number_of_edges())

    log.debug(
        'Subgraph coming from %s (seed type) %s (data) contains %d nodes and %d edges',
        seed_method,
        seed_data,
        result.number_of_nodes(),
        result.number_of_edges()
    )

    return result