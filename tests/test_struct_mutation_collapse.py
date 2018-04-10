# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel.constants import unqualified_edge_code
from pybel.dsl import *
from pybel.struct.mutation.collapse import collapse_by_central_dogma, collapse_nodes
from tests.utils import n

HGNC = 'HGNC'
GOBP = 'GOBP'
CHEBI = 'CHEBI'

g1 = gene(HGNC, '1')
r1 = rna(HGNC, '1')
p1 = protein(HGNC, '1')

g2 = gene(HGNC, '2')
r2 = rna(HGNC, '2')
p2 = protein(HGNC, '2')

g3 = gene(HGNC, '3')
r3 = rna(HGNC, '3')
p3 = protein(HGNC, '3')

g4 = gene(HGNC, '4')
m4 = mirna(HGNC, '4')

a5 = abundance(CHEBI, '5')
p5 = pathology(GOBP, '5')


class TestCollapseDownstream(unittest.TestCase):
    def test_collapse_1(self):
        graph = BELGraph()

        graph.add_entity(p1)
        graph.add_entity(p2)
        graph.add_entity(p3)

        graph.add_qualified_edge(p1, p3, relation=INCREASES, citation=n(), evidence=n())
        graph.add_qualified_edge(p2, p3, relation=DIRECTLY_INCREASES, citation=n(), evidence=n())

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        d = {
            p1: {p2}
        }

        collapse_nodes(graph, d)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges(), msg=graph.edges(data=True, keys=True))

    def test_collapse_dogma_1(self):
        graph = BELGraph()

        graph.add_entity(p1)
        graph.add_entity(r1)

        graph.add_translation(r1, p1)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(1, graph.number_of_edges())

        collapse_by_central_dogma(graph)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

    def test_collapse_dogma_2(self):
        graph = BELGraph()

        graph.add_entity(p1)
        graph.add_entity(r1)
        graph.add_entity(g1)

        graph.add_edge(r1, p1, key=unqualified_edge_code[TRANSLATED_TO], **{RELATION: TRANSLATED_TO})
        graph.add_edge(g1, r1, key=unqualified_edge_code[TRANSCRIBED_TO], **{RELATION: TRANSCRIBED_TO})

        self.assertEqual(3, graph.number_of_nodes())
        self.assertEqual(2, graph.number_of_edges())

        collapse_by_central_dogma(graph)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())

    def test_collapse_dogma_3(self):
        graph = BELGraph()

        graph.add_entity(r1)
        graph.add_entity(g1)

        graph.add_edge(g1, r1, key=unqualified_edge_code[TRANSCRIBED_TO], **{RELATION: TRANSCRIBED_TO})

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(1, graph.number_of_edges())

        collapse_by_central_dogma(graph)

        self.assertEqual(1, graph.number_of_nodes())
        self.assertEqual(0, graph.number_of_edges())
