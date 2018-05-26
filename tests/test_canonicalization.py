# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.canonicalize import _to_bel_lines_body, canonicalize_edge
from pybel.constants import BEL_DEFAULT_NAMESPACE, INCREASES, MODIFIER
from pybel.dsl import *
from pybel.dsl.edges import extracellular, intracellular
from tests.utils import n


class TestCanonicalizeEdge(unittest.TestCase):
    """This class houses all testing for the canonicalization of edges such that the relation/modifications can be used
    as a second level hash"""

    def setUp(self):
        """Builds a graph and adds some nodes to it"""
        self.graph = BELGraph()
        self.u = protein(name='u', namespace='TEST')
        self.v = protein(name='v', namespace='TEST')
        self.graph.add_entity(self.u)
        self.graph.add_entity(self.v)

    def add_edge(self, subject_modifier=None, object_modifier=None, annotations=None):
        """Wraps adding a sample edge

        :param subject_modifier:
        :param object_modifier:
        :param annotations:
        :rtype: tuple
        """
        key = self.graph.add_qualified_edge(
            self.u,
            self.v,
            relation=INCREASES,
            evidence=n(),
            citation=n(),
            subject_modifier=subject_modifier,
            object_modifier=object_modifier,
            annotations=annotations,
        )

        data = self.graph[self.u][self.v][key]

        return canonicalize_edge(data)

    def test_failure(self):
        with self.assertRaises(ValueError):
            self.add_edge(subject_modifier={MODIFIER: 'nope'})

    def test_canonicalize_edge_info(self):
        c1 = self.add_edge(
            annotations={
                'Species': '9606'
            }
        )

        c2 = self.add_edge(
            annotations={
                'Species': '9606'
            }
        )

        c3 = self.add_edge(
            subject_modifier=activity('tport'),
        )

        c4 = self.add_edge(
            subject_modifier=activity('tport', namespace=BEL_DEFAULT_NAMESPACE),
        )

        self.assertEqual(c1, c2)
        self.assertNotEqual(c1, c3)
        self.assertEqual(c3, c4)

    def test_subject_degradation_location(self):
        self.assertEqual(
            self.add_edge(
                subject_modifier=degradation()
            ),
            self.add_edge(
                subject_modifier=degradation()
            )
        )

        self.assertEqual(
            self.add_edge(
                subject_modifier=degradation(location=entity(name='somewhere', namespace='GOCC'))
            ),
            self.add_edge(
                subject_modifier=degradation(location=entity(name='somewhere', namespace='GOCC'))
            )
        )

        self.assertNotEqual(
            self.add_edge(
                subject_modifier=degradation()
            ),
            self.add_edge(
                subject_modifier=degradation(location=entity(name='somewhere', namespace='GOCC'))
            )
        )

    def test_translocation(self):
        self.assertEqual(
            self.add_edge(subject_modifier=secretion()),
            self.add_edge(subject_modifier=secretion()),
        )

        self.assertEqual(
            self.add_edge(subject_modifier=secretion()),
            self.add_edge(subject_modifier=translocation(from_loc=intracellular, to_loc=extracellular)),
        )


class TestSerializeBEL(unittest.TestCase):
    def setUp(self):
        self.citation = n()
        self.evidence = n()
        self.url = n()
        self.graph = BELGraph()
        self.graph.namespace_url['HGNC'] = self.url

    def help_check_lines(self, lines):
        """Checks the given lines match the graph built during the tests

        :type lines: list[str]
        """
        self.assertEqual(lines, list(_to_bel_lines_body(self.graph)))

    def test_simple(self):
        """Tests a scenario with a qualified edge, but no annotaitons"""
        self.graph.add_qualified_edge(
            protein(namespace='HGNC', name='YFG1'),
            protein(namespace='HGNC', name='YFG'),
            relation=INCREASES,
            citation=self.citation,
            evidence=self.evidence
        )

        self.assertEqual(2, self.graph.number_of_nodes())
        self.assertEqual(1, self.graph.number_of_edges())

        expected_lines = [
            '#' * 80,
            'SET Citation = {{"PubMed", "{}"}}'.format(self.citation),
            'SET SupportingText = "{}"'.format(self.evidence),
            'p(HGNC:YFG1) increases p(HGNC:YFG)',
            'UNSET SupportingText',
            'UNSET Citation'
        ]

        self.help_check_lines(expected_lines)

    def test_single_annotation(self):
        """Tests a scenario with a qualified edge, but no annotaitons"""
        a1, v1 = map(lambda _: n(), range(2))

        self.graph.add_qualified_edge(
            protein(namespace='HGNC', name='YFG1'),
            protein(namespace='HGNC', name='YFG'),
            relation=INCREASES,
            citation=self.citation,
            evidence=self.evidence,
            annotations={
                a1: {v1}
            }
        )

        self.assertEqual(2, self.graph.number_of_nodes())
        self.assertEqual(1, self.graph.number_of_edges())

        expected_lines = [
            '#' * 80,
            'SET Citation = {{"PubMed", "{}"}}'.format(self.citation),
            'SET SupportingText = "{}"'.format(self.evidence),
            'SET {} = "{}"'.format(a1, v1),
            'p(HGNC:YFG1) increases p(HGNC:YFG)',
            'UNSET {}'.format(a1),
            'UNSET SupportingText',
            'UNSET Citation'
        ]

        self.help_check_lines(expected_lines)

    def test_multiple_annotations(self):
        a1, v1, v2 = map(lambda _: n(), range(3))
        v1, v2 = sorted([v1, v2])

        self.graph.add_qualified_edge(
            protein(namespace='HGNC', name='YFG1'),
            protein(namespace='HGNC', name='YFG'),
            relation=INCREASES,
            citation=self.citation,
            evidence=self.evidence,
            annotations={
                a1: {v1, v2}
            }
        )

        self.assertEqual(2, self.graph.number_of_nodes())
        self.assertEqual(1, self.graph.number_of_edges())

        expected_lines = [
            '#' * 80,
            'SET Citation = {{"PubMed", "{}"}}'.format(self.citation),
            'SET SupportingText = "{}"'.format(self.evidence),
            'SET {} = {{"{}", "{}"}}'.format(a1, v1, v2),
            'p(HGNC:YFG1) increases p(HGNC:YFG)',
            'UNSET {}'.format(a1),
            'UNSET SupportingText',
            'UNSET Citation',
        ]

        self.help_check_lines(expected_lines)
