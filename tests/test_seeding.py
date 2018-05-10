# -*- coding: utf-8 -*-

import unittest
from copy import deepcopy

from pybel import BELGraph
from pybel.constants import IDENTIFIER
from pybel.examples import sialic_acid_graph
from pybel.examples.sialic_acid_example import cd33_phosphorylated, shp2, syk, trem2
from pybel.manager.models import Edge, Namespace, Network
from pybel.manager.query_manager import graph_from_edges
from tests.constants import TemporaryCacheClsMixin, TestGraphMixin
from tests.mocks import mock_bel_resources

chebi_url = 'https://arty.scai.fraunhofer.de/artifactory/bel/namespace/chebi/chebi-20170725.belns'
hgnc_url = 'https://arty.scai.fraunhofer.de/artifactory/bel/namespace/hgnc-human-genes/hgnc-human-genes-20170725.belns'
gobp_url = 'https://arty.scai.fraunhofer.de/artifactory/bel/namespace/go-biological-process/go-biological-process-20170725.belns'

class TestNodes(unittest.TestCase):
    def test_identifier_missing(self):
        trem2_copy = deepcopy(trem2)
        del trem2_copy[IDENTIFIER]

        self.assertNotIn(IDENTIFIER, trem2_copy)
        self.assertEqual('p(HGNC:TREM2)', trem2_copy.as_bel())


class TestSeeding(TemporaryCacheClsMixin, TestGraphMixin):
    """This module tests the seeding functions in the query manager"""

    @classmethod
    def setUpClass(cls):
        """Adds the sialic acid subgraph for all query tests"""
        super(TestSeeding, cls).setUpClass()

        @mock_bel_resources
        def insert(mock):
            """Inserts the Sialic Acid Subgraph using the mock resources"""
            cls.manager.insert_graph(sialic_acid_graph, store_parts=True)

        insert()

    def test_hgnc_namespace_existence(self):
        ns = self.manager.get_namespace_by_url(hgnc_url)
        self.assertIsNotNone(ns)
        self.assertEqual(hgnc_url, ns.url)

    def test_chebi_namespace_existence_b(self):
        ns = self.manager.get_namespace_by_url(chebi_url)
        self.assertIsNotNone(ns)
        self.assertEqual(chebi_url, ns.url)

    def test_namespace_existence_c(self):
        ns = self.manager.get_namespace_by_url(gobp_url)
        self.assertIsNotNone(ns)
        self.assertEqual(gobp_url, ns.url)


    def test_sialic_acid_in_node_store(self):
        name = 'sialic acid'
        n = self.manager.get_namespace_entry(chebi_url, name)
        self.assertIsNotNone(n)
        self.assertEqual(name, n.name)

    def test_network_existence(self):
        networks = self.manager.session.query(Network).all()

        self.assertEqual(1, len(networks))

    def test_edge_existence(self):
        self.assertEqual(11, self.manager.count_edges())

    def test_seed_by_pmid(self):
        pmids = ['26438529']

        edges = self.manager.query_edges_by_pubmed_identifiers(pmids)

        self.assertLess(0, len(edges))

    def test_seed_by_pmid_no_result(self):
        edges = self.manager.query_edges_by_pubmed_identifiers('11111')

        self.assertEqual(0, len(edges))

    def test_seed_by_induction_raise(self):
        with self.assertRaises(ValueError):
            self.manager.query_induction([])

    def test_seed_by_induction(self):
        trem2_copy = deepcopy(trem2)
        del trem2_copy[IDENTIFIER]

        syk_copy = deepcopy(syk)
        del syk_copy[IDENTIFIER]

        shp2_copy = deepcopy(shp2)
        del shp2_copy[IDENTIFIER]

        shp2_model = self.manager.get_node_by_dict(shp2)
        self.assertIsNotNone(shp2_model)
        self.assertEqual(shp2_copy, shp2_model.as_bel())

        syk_model = self.manager.get_node_by_dict(syk)
        self.assertIsNotNone(syk_model)
        self.assertEqual(syk_copy, syk_model.as_bel())

        trem2_model = self.manager.get_node_by_dict(trem2)
        self.assertIsNotNone(trem2_model)
        self.assertEqual(trem2_copy, trem2_model.as_bel())

        edges = self.manager.query_induction([shp2_model, syk_model, trem2_model])
        self.assertEqual(2, len(edges))
        for edge in edges:
            self.assertIsInstance(edge, Edge)

        graph = graph_from_edges(edges)
        self.assertIsInstance(graph, BELGraph)
        self.assertEqual(2, graph.number_of_edges())
        self.assertEqual(3, graph.number_of_nodes(), msg='Nodes: {}'.format(graph.nodes()))

        self.assertHasNode(graph, trem2_copy)
        self.assertHasNode(graph, syk_copy)
        self.assertHasNode(graph, shp2_copy)

        self.assertIn(trem2_copy, graph[syk_copy])
        self.assertIn(syk_copy, graph[shp2_copy])

    def test_seed_by_neighbors(self):
        node = self.manager.get_node_by_dict(shp2)
        edges = self.manager.query_neighbors([node])
        self.assertEqual(2, len(edges))

        graph = graph_from_edges(edges)

        syk_copy = deepcopy(syk)
        del syk_copy[IDENTIFIER]
        shp2_copy = deepcopy(shp2)
        del shp2_copy[IDENTIFIER]
        cd33_phosphorylated_copy = deepcopy(cd33_phosphorylated)
        del cd33_phosphorylated_copy[IDENTIFIER]

        self.assertIn(cd33_phosphorylated_copy, graph)
        self.assertIn(syk_copy, graph)
        self.assertIn(shp2_copy, graph)

        self.assertEqual(3, graph.number_of_nodes(), msg='Nodes:\n{}'.format('\n'.join(map(str, graph))))
        self.assertEqual(2, graph.number_of_edges())
