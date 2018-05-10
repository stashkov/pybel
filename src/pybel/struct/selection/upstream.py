# -*- coding: utf-8 -*-

import logging

from ..graph import BELGraph
from ..pipeline import mutator
from ...constants import CAUSAL_RELATIONS, RELATION

__all__ = [
    'get_upstream_causal_subgraph',
]

log = logging.getLogger(__name__)


@mutator
def get_upstream_causal_subgraph(graph, nbunch):
    """Induces a subgraph from all of the upstream causal entities of the nodes in the nbunch

    :param pybel.BELGraph graph: A BEL graph
    :param BaseEntity or iter[BaseEntity] nbunch: A BEL node or iterable of BEL nodes
    :rtype: pybel.BELGraph
    """
    rv = BELGraph()
    rv.add_edges_from(
        (u, v, key, data)
        for u, v, key, data in graph.in_edges(nbunch, keys=True, data=True)
        if data[RELATION] in CAUSAL_RELATIONS

    )
    return rv
