# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.constants import ANNOTATIONS, INCREASES
from pybel.dsl import protein
from pybel.examples.statin_example import (
    avorastatin, ec_11134, ec_11188, fluvastatin, hmgcr, hmgcr_inhibitor, mevinolinic_acid, statin, statin_graph,
    synthetic_statin,
)
from pybel.struct.mutation import infer_child_relations, strip_annotations
from pybel.struct.mutation.transfer import iter_children


class TestMutations(unittest.TestCase):
    def test_strip_annotations(self):
        u = protein(namespace='HGNC', name='U')
        v = protein(namespace='HGNC', name='V')

        annotations = {
            'A': {'B': True}
        }

        graph = BELGraph()
        key = graph.add_qualified_edge(
            u,
            v,
            relation=INCREASES,
            citation='123456',
            evidence='Fake',
            annotations=annotations,
        )

        self.assertIn(ANNOTATIONS, graph[u][v][key])

        self.assertEqual(annotations, graph.get_edge_annotations(u, v, key))
        strip_annotations(graph)
        self.assertNotIn(ANNOTATIONS, graph[u][v][key])


class TestTransfer(unittest.TestCase):
    def test_get_children(self):
        children = list(iter_children(statin_graph, hmgcr_inhibitor))

        self.assertNotEqual(0, len(children), msg='no children found')
        self.assertIn(mevinolinic_acid, children, msg='direct child not found')

    def test_infer(self):
        graph = statin_graph.copy()
        self.assertEqual(9, graph.number_of_nodes())
        self.assertEqual(8, graph.number_of_edges())

        self.assertNotIn(ec_11134, graph.edge[fluvastatin])
        self.assertNotIn(ec_11188, graph.edge[fluvastatin])
        self.assertNotIn(ec_11134, graph.edge[avorastatin])
        self.assertNotIn(ec_11188, graph.edge[avorastatin])
        self.assertNotIn(ec_11134, graph.edge[synthetic_statin])
        self.assertNotIn(ec_11188, graph.edge[synthetic_statin])
        self.assertNotIn(ec_11134, graph.edge[statin])
        self.assertNotIn(ec_11188, graph.edge[statin])
        self.assertNotIn(ec_11134, graph.edge[mevinolinic_acid])
        self.assertNotIn(ec_11188, graph.edge[mevinolinic_acid])
        self.assertIn(ec_11134, graph.edge[hmgcr_inhibitor])
        self.assertIn(ec_11188, graph.edge[hmgcr_inhibitor])

        infer_child_relations(graph, hmgcr_inhibitor)

        self.assertIn(ec_11134, graph.edge[fluvastatin])
        self.assertIn(ec_11188, graph.edge[fluvastatin])
        self.assertIn(ec_11134, graph.edge[avorastatin])
        self.assertIn(ec_11188, graph.edge[avorastatin])
        self.assertIn(ec_11134, graph.edge[synthetic_statin])
        self.assertIn(ec_11188, graph.edge[synthetic_statin])
        self.assertIn(ec_11134, graph.edge[statin])
        self.assertIn(ec_11188, graph.edge[statin])
        self.assertIn(ec_11134, graph.edge[mevinolinic_acid])
        self.assertIn(ec_11188, graph.edge[mevinolinic_acid])
        self.assertIn(ec_11134, graph.edge[hmgcr_inhibitor])
        self.assertIn(ec_11188, graph.edge[hmgcr_inhibitor])

        self.assertEqual(9, graph.number_of_nodes())
        self.assertEqual(18, graph.number_of_edges())

        infer_child_relations(graph, ec_11134)

        self.assertIn(hmgcr, graph.edge[fluvastatin])
        self.assertIn(hmgcr, graph.edge[avorastatin])
        self.assertIn(hmgcr, graph.edge[synthetic_statin])
        self.assertIn(hmgcr, graph.edge[statin])
        self.assertIn(hmgcr, graph.edge[mevinolinic_acid])
        self.assertIn(hmgcr, graph.edge[hmgcr_inhibitor])

        self.assertEqual(9, graph.number_of_nodes())
        self.assertEqual(24, graph.number_of_edges())

        self.assertEqual(9, statin_graph.number_of_nodes(), msg='original graph nodes should not be modified')
        self.assertEqual(8, statin_graph.number_of_edges(), msg='original graph edges should not be modified')

    @unittest.skip('not yet sure if this is necessary')
    def test_does_not_redo(self):
        """Tests that :func:`propagate_node_relations` does not add the same edges twice"""
        graph = statin_graph.copy()
        self.assertEqual(9, graph.number_of_nodes())
        self.assertEqual(8, graph.number_of_edges())

        infer_child_relations(graph, hmgcr_inhibitor)
        self.assertEqual(9, graph.number_of_nodes())
        self.assertEqual(18, graph.number_of_edges())

        infer_child_relations(graph, hmgcr_inhibitor)
        self.assertEqual(9, graph.number_of_nodes())
        self.assertEqual(18, graph.number_of_edges(), msg='edges should not be added again')


if __name__ == '__main__':
    unittest.main()
