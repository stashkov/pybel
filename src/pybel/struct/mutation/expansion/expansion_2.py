# -*- coding: utf-8 -*-

import itertools as itt
import logging
from collections import defaultdict

from ... import pipeline
from ...filters.edge_filters import and_edge_predicates
from ....constants import RELATION

__all__ = [
    'expand_internal',
]

log = logging.getLogger(__name__)


# TODO should this bother checking multiple relationship types?
@pipeline.uni_in_place_mutator
def expand_internal(universe, graph, edge_filters=None):
    """Edges between entities in the subgraph that pass the given filters

    :param pybel.BELGraph universe: The full graph
    :param pybel.BELGraph graph: A subgraph to find the upstream information
    :param edge_filters: Optional list of edge filter functions (graph, node, node, key, data) -> bool
    :type edge_filters: list or lambda
    """
    edge_filter = and_edge_predicates(edge_filters)

    for u, v in itt.product(graph, repeat=2):
        if graph.has_edge(u, v) or not universe.has_edge(u, v):
            continue

        rs = defaultdict(list)
        for k in universe.edge[u][v]:
            if not edge_filter(universe, u, v, k):
                continue
            d = universe[u][v][k]
            rs[d[RELATION]].append(d)

        if 1 == len(rs):
            relation = list(rs)[0]
            for d in rs[relation]:
                graph.add_edge(u, v, attr_dict=d)
        else:
            log.debug('Multiple relationship types found between %s and %s', u, v)
