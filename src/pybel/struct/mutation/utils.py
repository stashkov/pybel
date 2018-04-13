# -*- coding: utf-8 -*-

from ..pipeline import uni_in_place_mutator


@uni_in_place_mutator
def ensure_node_from_universe(universe, graph, node, raise_for_missing=False):
    """Makes sure that the subgraph has the given node

    :param pybel.BELGraph universe: The universe of all knowledge
    :param pybel.BELGraph graph: The target BEL graph
    :param BaseEntity node: A BEL node
    :param bool raise_for_missing: Should an error be thrown if the given node is not in the universe?
    """
    if raise_for_missing and node not in universe:
        raise IndexError('{} not in {}'.format(node, universe.name))

    if node not in graph:
        graph.add_entity(node)
