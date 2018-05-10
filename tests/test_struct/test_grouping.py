# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.constants import INCREASES
from pybel.dsl import protein
from pybel.struct.grouping import get_subgraphs_by_annotation


class TestGrouping(unittest.TestCase):

    def test_get_subgraphs_by_annotation(self):
        graph = BELGraph()
        a, b, c, d = [protein(namespace='test', name=str(i)) for i in range(4)]
        citation, evidence = '', ''
        graph.add_qualified_edge(a, b, INCREASES, citation, evidence, {'subgraph': {'1', '2'}})
        graph.add_qualified_edge(a, c, INCREASES, citation, evidence, {'subgraph': {'1'}})
        graph.add_qualified_edge(b, d, INCREASES, citation, evidence, {'subgraph': {'1', '2'}})
        graph.add_qualified_edge(a, d, INCREASES, citation, evidence, {'subgraph': {'2'}})
        graph.add_qualified_edge(c, d, INCREASES, citation, evidence)

        subgraphs = get_subgraphs_by_annotation(graph, annotation='subgraph', keep_undefined=False)

        self.assertIn('1', subgraphs)
        self.assertIn('2', subgraphs)

        subgraph_1 = subgraphs['1']

        self.assertIn(a, subgraph_1)
        self.assertIn(b, subgraph_1)
        self.assertIn(c, subgraph_1)
        self.assertIn(d, subgraph_1)

        self.assertIn(b, subgraph_1[a])
        self.assertIn(c, subgraph_1[a])
        self.assertIn(d, subgraph_1[b])
        self.assertNotIn(d, subgraph_1[a])
        self.assertNotIn(d, subgraph_1[c])

        subgraph_2 = subgraphs['2']

        self.assertIn(a, subgraph_2)
        self.assertIn(b, subgraph_2)
        self.assertNotIn(c, subgraph_2)
        self.assertIn(d, subgraph_2)

        self.assertIn(b, subgraph_2[a])
        self.assertNotIn(c, subgraph_2[a])
        self.assertIn(d, subgraph_2[b])
        self.assertIn(d, subgraph_2[a])
