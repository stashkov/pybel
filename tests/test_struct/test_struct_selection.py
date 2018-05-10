import unittest

from pybel import BELGraph
from pybel.constants import ASSOCIATION, INCREASES
from pybel.dsl import protein
from pybel.struct import get_causal_subgraph, get_upstream_causal_subgraph


class TestSelection(unittest.TestCase):
    def test_get_causal_subgraph(self):
        graph = BELGraph()
        a, b, c, d = [protein(namespace='test', name=str(i)) for i in range(4)]
        citation, evidence = '', ''

        graph.add_qualified_edge(a, b, INCREASES, citation, evidence)
        graph.add_qualified_edge(b, c, INCREASES, citation, evidence)
        graph.add_qualified_edge(a, c, ASSOCIATION, citation, evidence)
        graph.add_qualified_edge(d, c, ASSOCIATION, citation, evidence)
        self.assertEqual(4, graph.number_of_nodes())
        self.assertEqual(4, graph.number_of_edges())

        subgraph = get_causal_subgraph(graph)

        self.assertIn(a, subgraph)
        self.assertIn(b, subgraph)
        self.assertIn(c, subgraph)
        self.assertNotIn(d, subgraph)
        self.assertEqual(3, subgraph.number_of_nodes())

        self.assertIn(b, subgraph[a])
        self.assertIn(c, subgraph[b])
        self.assertEqual(2, subgraph.number_of_edges())

    def test_get_upstream_causal_subgraph(self):
        graph = BELGraph()
        a, b, c, d = [protein(namespace='test', name=str(i)) for i in range(4)]
        citation, evidence = '', ''

        graph.add_qualified_edge(a, b, INCREASES, citation, evidence)
        graph.add_qualified_edge(b, c, INCREASES, citation, evidence)
        graph.add_qualified_edge(a, c, ASSOCIATION, citation, evidence)
        graph.add_qualified_edge(d, c, ASSOCIATION, citation, evidence)
        self.assertEqual(4, graph.number_of_nodes())
        self.assertEqual(4, graph.number_of_edges())

        subgraph = get_upstream_causal_subgraph(graph, a)
        self.assertEqual(0, subgraph.number_of_nodes())

        subgraph = get_upstream_causal_subgraph(graph, b)
        self.assertIn(a, subgraph)
        self.assertIn(b, subgraph)
        self.assertNotIn(c, subgraph)
        self.assertNotIn(d, subgraph)
        self.assertEqual(2, subgraph.number_of_nodes())

        self.assertIn(b, subgraph[a])
        self.assertEqual(1, subgraph.number_of_edges())

        subgraph = get_upstream_causal_subgraph(graph, c)
        self.assertNotIn(a, subgraph)
        self.assertIn(b, subgraph)
        self.assertIn(c, subgraph)
        self.assertNotIn(d, subgraph)
        self.assertEqual(2, subgraph.number_of_nodes())

        self.assertIn(c, subgraph[b])
        self.assertEqual(1, subgraph.number_of_edges())
