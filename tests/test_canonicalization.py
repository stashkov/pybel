# -*- coding: utf-8 -*-

import unittest

from pybel import BELGraph
from pybel.canonicalize import canonicalize_edge
from pybel.constants import BEL_DEFAULT_NAMESPACE, INCREASES, MODIFIER
from pybel.dsl import *
from pybel.dsl.edges import extracellular, intracellular
from tests.utils import n


class TestCanonicalizeEdge(unittest.TestCase):
    """This class houses all testing for the canonicalization of edges such that the relation/modifications can be used
    as a second level hash"""

    def setUp(self):
        """Builds a graph and adds some nodes to it"""
        self.g = BELGraph()
        self.u = protein(name='u', namespace='TEST')
        self.v = protein(name='v', namespace='TEST')
        self.g.add_node_from_data(self.u)
        self.g.add_node_from_data(self.v)
        self.key = 0

    def get_data(self, k):
        """Gets the data associated with the sample key"""
        return self.g.edge[self.u][self.v][k]

    def add_edge(self, subject_modifier=None, object_modifier=None, annotations=None):
        """Wraps adding a sample edge

        :param subject_modifier:
        :param object_modifier:
        :param annotations:
        :rtype: str
        """
        self.key += 1

        self.g.add_qualified_edge(
            self.u,
            self.v,
            relation=INCREASES,
            evidence=n(),
            citation=n(),
            subject_modifier=subject_modifier,
            object_modifier=object_modifier,
            annotations=annotations,
            key=self.key
        )

        return canonicalize_edge(self.get_data(self.key))

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
