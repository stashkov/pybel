# -*- coding: utf-8 -*-

import logging

from ..induction import get_upstream_causal_subgraph
from ...graph import left_full_join
from ...pipeline import uni_in_place_mutator

__all__ = [
    'expand_upstream_causal_subgraph',
]

log = logging.getLogger(__name__)


@uni_in_place_mutator
def expand_upstream_causal_subgraph(universe, graph):
    """Adds the upstream causal relations to the given subgraph

    :param pybel.BELGraph universe: A BEL graph representing the universe of all knowledge
    :param pybel.BELGraph graph: The target BEL graph to enrich with upstream causal controllers of contained nodes
    """
    upstream = get_upstream_causal_subgraph(universe, graph)
    left_full_join(graph, upstream)
