# -*- coding: utf-8 -*-

"""
Reference for testing Flask

- Flask Documentation http://flask.pocoo.org/docs/0.12/testing/
- Flask Cookbook: http://flask.pocoo.org/docs/0.12/tutorial/testing/
"""

import logging
import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel.dsl import *
from pybel.examples import egf_graph
from pybel.examples.homology_example import (
    homology_graph, mouse_csf1_protein, mouse_csf1_rna,
    mouse_mapk1_protein, mouse_mapk1_rna,
)
from pybel.examples.sialic_acid_example import dap12, shp1, shp2, sialic_acid_graph, syk, trem2
from pybel.struct.mutation import infer_central_dogma
from pybel.struct.mutation.collapse import collapse_by_central_dogma_to_genes
from pybel.struct.mutation.expansion import expand_internal
from pybel.struct.pipeline import Pipeline
from pybel.struct.query import Query
from pybel.struct.selection import get_subgraph_by_annotation_value
from tests.struct_mocks import MockQueryManager

log = logging.getLogger(__name__)
log.setLevel(10)

HGNC = 'HGNC'

protein_a = protein(namespace=HGNC, name='a')
protein_b = protein(namespace=HGNC, name='b')
gene_c = gene(namespace=HGNC, name='c')
rna_d = rna(namespace=HGNC, name='d')
protein_e = protein(namespace=HGNC, name='e')
gene_f = gene(namespace=HGNC, name='f')
protein_g = protein(namespace=HGNC, name='g')
protein_h = protein(namespace=HGNC, name='h')
protein_i = protein(namespace=HGNC, name='i')
protein_j = protein(namespace=HGNC, name='j')


def make_graph_1():
    """
    :rtype: pybel.BELGraph
    """
    graph = BELGraph(
        name='PyBEL Tools Example Network 1',
        version='1.1.0',
        description='Example Network for PyBEL Tools Tests',
        authors='Daniel Domingo-Fernández and Charles Tapley Hoyt',
        contact='charles.hoyt@scai.fraunhofer.de',
    )

    graph.add_node_from_data(protein_a)
    graph.add_node_from_data(protein_b)
    graph.add_node_from_data(gene_c)
    graph.add_node_from_data(rna_d)

    graph.add_qualified_edge(protein_a, protein_b, relation=INCREASES, citation='1', evidence='Evidence 1',
                             annotations={'Annotation': 'foo'})

    graph.add_qualified_edge(rna_d, protein_a, relation=INCREASES, citation='2', evidence='Evidence 2',
                             annotations={'Annotation': 'foo'})

    graph.add_qualified_edge(gene_c, protein_b, relation=DECREASES, citation='3', evidence='Evidence 3',
                             annotations={'Annotation': 'foo'})

    return graph


class TestMockManager(unittest.TestCase):
    def test_make(self):
        manager = MockQueryManager()
        self.assertEqual(0, manager.count_networks())

    def test_make_with_graph(self):
        graph_1 = make_graph_1()
        manager = MockQueryManager(graphs=[graph_1])
        self.assertEqual(1, manager.count_networks())


class QueryTest(unittest.TestCase):
    def setUp(self):
        super(QueryTest, self).setUp()
        self.manager = MockQueryManager()
        self.graph_1 = make_graph_1()

    def test_pipeline(self):
        test_graph_1 = self.graph_1  # Defined in test.constants.TestNetworks

        infer_central_dogma(test_graph_1)

        self.assertEqual(9, test_graph_1.number_of_nodes())  # 4 nodes already there +  2*2 proteins + 1 (rna)
        self.assertEqual(8, test_graph_1.number_of_edges())  # 3 already there + 2*2 proteins + 1 (rna)

        network = self.manager.insert_graph(test_graph_1)

        pipeline = Pipeline()
        pipeline.append(collapse_by_central_dogma_to_genes)

        query = Query(
            network_ids=[network.id],
            pipeline=pipeline
        )
        result_graph = query.run(self.manager)

        self.assertEqual(4, result_graph.number_of_nodes())  # same number of nodes than there were
        self.assertEqual(3, result_graph.number_of_edges())  # same number of edges than there were

    def test_pipeline_2(self):
        test_network = self.graph_1  # Defined in test.constants.TestNetworks

        infer_central_dogma(test_network)

        network = self.manager.insert_graph(test_network)
        network_id = network.id

        pipeline = Pipeline()
        pipeline.append(get_subgraph_by_annotation_value, 'Annotation', 'foo')

        query = Query(network_ids=[network_id])
        query.append_seeding_neighbors([protein_a])
        query.pipeline = pipeline

        result = query.run(self.manager, in_place=False)
        self.assertIsNotNone(result, msg='Query returned none')

        self.assertEqual(3, result.number_of_nodes())  # only expanded to node protein_a and gene_c
        self.assertEqual(2, result.number_of_edges())  # three nodes with two relationships

    def test_query_multiple_networks(self):
        egf_network = self.manager.insert_graph(sialic_acid_graph.copy())
        sialic_acid_network = self.manager.insert_graph(egf_graph.copy())

        query = Query()
        query.append_network(egf_network.id)
        query.append_network(sialic_acid_network.id)
        query.append_seeding_neighbors([syk])
        query.append_pipeline(infer_central_dogma)

        result = query.run(self.manager, in_place=False)
        self.assertIsNotNone(result, msg='Query returned none')

        self.assertIn(shp1, result)
        self.assertIn(shp2, result)
        self.assertIn(trem2, result)
        self.assertIn(dap12, result)

        self.assertEqual(15, result.number_of_nodes())
        self.assertEqual(14, result.number_of_edges())

    def test_get_subgraph_by_annotation_value(self):
        graph = homology_graph.copy()

        result = get_subgraph_by_annotation_value(graph, 'Species', '10090')
        self.assertIsNotNone(result, msg='Query returned none')
        self.assertIsInstance(result, BELGraph)

        self.assertIn(mouse_mapk1_protein, result)
        self.assertIn(mouse_csf1_protein, result)

        self.assertEqual(2, result.number_of_nodes())
        self.assertEqual(1, result.number_of_edges())

    def test_seeding_1(self):
        test_network_1 = self.manager.insert_graph(homology_graph.copy())

        query = Query(network_ids=[test_network_1.id])
        query.append_seeding_neighbors([mouse_csf1_rna, mouse_mapk1_rna])

        result = query.run(self.manager, in_place=False)
        self.assertIsNotNone(result, msg='Query returned none')
        self.assertIsInstance(result, BELGraph)

        self.assertIn(mouse_mapk1_rna, result)
        self.assertIn(mouse_csf1_rna, result)
        self.assertIn(mouse_mapk1_protein, result)
        self.assertIn(mouse_csf1_protein, result)

        self.assertEqual(6, result.number_of_nodes())
        self.assertEqual(4, result.number_of_edges())

    def test_seeding_with_pipeline(self):
        test_network_1 = self.manager.insert_graph(homology_graph.copy())

        query = Query(network_ids=[test_network_1.id])
        query.append_seeding_neighbors([mouse_csf1_rna, mouse_mapk1_rna])
        query.append_pipeline(expand_internal)
        result = query.run(self.manager, in_place=False)
        self.assertIsNotNone(result, msg='Query returned none')
        self.assertIsInstance(result, BELGraph)

        self.assertIn(mouse_mapk1_rna, result)
        self.assertIn(mouse_csf1_rna, result)
        self.assertIn(mouse_mapk1_protein, result)
        self.assertIn(mouse_csf1_protein, result)

        self.assertEqual(6, result.number_of_nodes())
        self.assertEqual(5, result.number_of_edges())

    def test_query_multiple_networks_with_api(self):
        test_network_1 = self.manager.insert_graph(homology_graph.copy())

        pipeline = Pipeline()
        pipeline.append(expand_internal)
        pipeline.append(get_subgraph_by_annotation_value, 'Species', '10090')

        query = Query(
            network_ids=[test_network_1.id],
            pipeline=pipeline
        )
        query.append_seeding_neighbors([mouse_csf1_rna, mouse_mapk1_rna])

        result = query.run(self.manager, in_place=False)
        self.assertIsNotNone(result, msg='Query returned none')

        self.assertEqual(2, result.number_of_nodes())
        self.assertIn(mouse_mapk1_protein, result)
        self.assertIn(mouse_csf1_protein, result)

        self.assertEqual(1, result.number_of_edges())
