# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.dsl import protein
from pybel.struct.filters.node_predicate_builders import data_missing_key_builder


class TestPredicateBuilders(unittest.TestCase):

    def test_data_missing_key_builder(self):
        """Test the data_missing_key_builder function."""
        data_missing_key = data_missing_key_builder('test_key')

        graph = BELGraph()
        node = protein(namespace='test', name='test')
        graph.add_entity(node)

        self.assertTrue(data_missing_key(graph, node), msg='Should not have any stuff in')

        graph.node[node]['test_key'] = 714
        self.assertFalse(data_missing_key(graph, node))
