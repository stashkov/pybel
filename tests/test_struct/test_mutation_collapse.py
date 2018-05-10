# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.constants import DIRECTLY_INCREASES, INCREASES
from pybel.dsl import protein
from pybel.struct.mutation.collapse.inconsitent import pair_is_consistent


class TestPredicateBuilders(unittest.TestCase):

    def test_pair_is_consistent(self):
        """Test the data_missing_key_builder function."""
        graph = BELGraph()
        a, b, c = [protein(namespace='test', name=str(i)) for i in range(3)]

        graph.add_qualified_edge(a, b, INCREASES, 'asfa', 'asa')
        graph.add_qualified_edge(a, b, INCREASES, 'asfaasf', 'asfaasa')
        graph.add_qualified_edge(a, b, INCREASES, 'asfasfaasf', 'asfasfaasa')

        self.assertTrue(pair_is_consistent(graph, a, b))

        graph.add_qualified_edge(a, b, DIRECTLY_INCREASES, 'asfasfasfaasf', 'asfasfaasfasa')

        self.assertFalse(pair_is_consistent(graph, a, b))

    def test_pair_is_consistent_missing_edge(self):
        graph = BELGraph()
        a, b = [protein(namespace='test', name=str(i)) for i in range(2)]
        graph.add_entity(a)
        graph.add_entity(b)
        self.assertFalse(pair_is_consistent(graph, a, b))
