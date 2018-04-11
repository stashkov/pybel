# -*- coding: utf-8 -*-

import logging
import unittest

from pybel.examples.egf_example import nfkb_complex, rela, nfkb1, nfkb2, egf_graph
from pybel.struct.mutation import build_delete_node_by_hash, build_expand_node_neighborhood_by_hash, infer_central_dogma
from pybel.struct.pipeline import MissingPipelineFunctionError, Pipeline, assert_is_mapped_to_pipeline, mapped
from tests.struct_mocks import MockQueryManager


log = logging.getLogger(__name__)
log.setLevel(10)


class TestEgfExample(unittest.TestCase):
    """Random test for mutation functions"""

    def setUp(self):
        self.graph = egf_graph.copy()
        self.original_number_nodes = self.graph.number_of_nodes()
        self.original_number_edges = self.graph.number_of_edges()

    def check_original_unchanged(self):
        self.assertEqual(self.original_number_nodes, self.graph.number_of_nodes(),
                         msg='original graph nodes should remain unchanged')
        self.assertEqual(self.original_number_edges, self.graph.number_of_edges(),
                         msg='original graph edges should remain unchanged')


class TestPipelineFailures(unittest.TestCase):

    def test_assert_failure(self):
        with self.assertRaises(MissingPipelineFunctionError):
            assert_is_mapped_to_pipeline('missing function')

    def test_assert_success(self):
        m = list(mapped)
        self.assertLess(0, len(m))
        m = m[0]
        assert_is_mapped_to_pipeline(m)

    def test_get_function_failure(self):
        pass

    def test_get_function_success(self):
        pass

    def test_fail_add(self):
        pipeline = Pipeline()

        with self.assertRaises(MissingPipelineFunctionError):
            pipeline.append('missing function')

    def test_fail_build(self):
        protocol = [{'function': 'missing function'}]
        with self.assertRaises(MissingPipelineFunctionError):
            Pipeline(protocol=protocol)

    def test_fail_from_json(self):
        protocol = [{'function': 'missing function'}]
        with self.assertRaises(MissingPipelineFunctionError):
            Pipeline.from_json(protocol)


class TestPipeline(TestEgfExample):
    def test_central_dogma_is_registered(self):
        self.assertIn('infer_central_dogma', mapped)

    def test_pipeline_by_string(self):
        pipeline = Pipeline()
        pipeline.append('infer_central_dogma')
        result = pipeline.run(self.graph, in_place=False)

        self.assertEqual(32, result.number_of_nodes())

        for node in self.graph:
            self.assertIn(node, result)

        self.check_original_unchanged()

    def test_pipeline_by_function(self):
        pipeline = Pipeline()
        pipeline.append(infer_central_dogma)
        result = pipeline.run(self.graph, in_place=False)

        self.assertEqual(32, result.number_of_nodes())

        for node in self.graph:
            self.assertIn(node, result)

        self.check_original_unchanged()


class TestBoundMutation(TestEgfExample):
    """Random test for mutation functions"""

    def setUp(self):
        super(TestBoundMutation, self).setUp()

        self.manager = MockQueryManager([self.graph])

        self.delete_node_by_hash = build_delete_node_by_hash(self.manager)
        self.expand_node_neighborhood_by_hash = build_expand_node_neighborhood_by_hash(self.manager)

    def test_mock_contents(self):
        self.assertIn(nfkb_complex, self.manager.graphs[0], msg='Graph missing NFKB complex')
        self.assertIn(rela, self.manager.graphs[0], msg='Graph missing RELA')

        self.assertIn(nfkb_complex, self.manager.sha512_to_node.values(), msg='NFKB is unindexed')
        self.assertIn(rela, self.manager.sha512_to_node.values(), msg='RELA is unindexed')

        self.assertIn(nfkb1.as_sha512(), self.manager.sha512_to_node, msg='NFKB is unindexed')
        self.assertIn(rela.as_sha512(), self.manager.sha512_to_node, msg='RELA is unindexed')

    def test_functions_registered(self):
        self.assertIn('delete_node_by_hash', mapped)
        self.assertIn('expand_node_neighborhood_by_hash', mapped)

    def test_bound_mutation(self):
        """Tests when a node is deleted then re-expanded"""
        pipeline = Pipeline(universe=self.graph)
        pipeline.append(self.delete_node_by_hash, nfkb_complex.as_sha512())
        pipeline.append(self.expand_node_neighborhood_by_hash, rela.as_sha512())

        result = pipeline.run(self.graph, in_place=False)

        self.check_original_unchanged()

        self.assertEqual(self.original_number_nodes, result.number_of_nodes())
        self.assertGreater(self.original_number_edges, result.number_of_edges())

    def test_bound_mutation_named(self):
        """Tests when a node is deleted then re-expanded"""
        pipeline = Pipeline(universe=self.graph)
        pipeline.append('delete_node_by_hash', nfkb_complex.as_sha512())
        pipeline.append('expand_node_neighborhood_by_hash', rela.as_sha512())

        result = pipeline.run(self.graph, in_place=False)

        self.check_original_unchanged()

        self.assertEqual(self.original_number_nodes, result.number_of_nodes())
        self.assertGreater(self.original_number_edges, result.number_of_edges())
