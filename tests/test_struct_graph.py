# -*- coding: utf-8 -*-

import unittest

from six import string_types

from pybel import BELGraph
from pybel.constants import (
    CITATION_REFERENCE, CITATION_TYPE, CITATION_TYPE_PUBMED, FUNCTION, HAS_VARIANT, IDENTIFIER,
    INCREASES, NAME, NAMESPACE, PROTEIN, unqualified_edge_code,
)
from pybel.dsl import *
from pybel.examples import sialic_acid_graph
from tests.utils import n


class TestStruct(unittest.TestCase):
    def test_example_metadata(self):
        self.assertIsNotNone(sialic_acid_graph.name)
        self.assertIsNotNone(sialic_acid_graph.version)
        self.assertIsNotNone(sialic_acid_graph.description)
        self.assertIsNone(sialic_acid_graph.disclaimer)
        self.assertIsNone(sialic_acid_graph.copyright)
        self.assertIsNone(sialic_acid_graph.license)

    def test_add_simple(self):
        g = BELGraph()

        namespace, name = n(), n()

        g.add_node_from_data(protein(namespace=namespace, name=name))
        self.assertEqual(1, g.number_of_nodes())

        g.add_node_from_data(protein(namespace=namespace, name=name))
        self.assertEqual(1, g.number_of_nodes())

    def test_str_kwargs(self):
        (
            name,
            version,
            description,
            authors,
            contact,
            licenses,
            copyrights,
            disclaimer
        ) = [n() for _ in range(8)]

        g = BELGraph(
            name=name,
            version=version,
            description=description,
            authors=authors,
            contact=contact,
            license=licenses,
            copyright=copyrights,
            disclaimer=disclaimer
        )

        self.assertEqual(name, g.name)
        self.assertEqual(version, g.version)
        self.assertEqual(description, g.description)
        self.assertEqual(authors, g.authors)
        self.assertEqual(contact, g.contact)
        self.assertEqual(licenses, g.license)
        self.assertEqual(copyrights, g.copyright)
        self.assertEqual(disclaimer, g.disclaimer)

        self.assertEqual('{name} v{version}'.format(name=name, version=version), str(g))

    def test_name(self):
        (
            name,
            version,
            description,
            authors,
            contact,
            licenses,
            copyrights,
            disclaimer
        ) = [n() for _ in range(8)]

        g = BELGraph()

        g.name = name
        g.version = version
        g.description = description
        g.authors = authors
        g.contact = contact
        g.license = licenses
        g.copyright = copyrights
        g.disclaimer = disclaimer

        self.assertEqual(name, g.name)
        self.assertEqual(version, g.version)
        self.assertEqual(description, g.description)
        self.assertEqual(authors, g.authors)
        self.assertEqual(contact, g.contact)
        self.assertEqual(licenses, g.license)
        self.assertEqual(copyrights, g.copyright)
        self.assertEqual(disclaimer, g.disclaimer)

        self.assertEqual('{name} v{version}'.format(name=name, version=version), str(g))

    def test_citation_type_error(self):
        g = BELGraph()

        with self.assertRaises(TypeError):
            g.add_qualified_edge(
                protein(namespace='TEST', name='YFG1'),
                protein(namespace='TEST', name='YFG2'),
                relation=INCREASES,
                evidence=n(),
                citation=5
            )


class TestDSL(unittest.TestCase):
    def test_add_robust_node(self):
        g = BELGraph()
        p = protein(name='yfg', namespace='test', identifier='1')
        g.add_node_from_data(p)
        self.assertIn(p, g)

    def test_add_identified_node(self):
        """What happens when a node with only an identifier is added to a graph"""
        g = BELGraph()
        p = protein(namespace='test', identifier='1')
        g.add_node_from_data(p)
        self.assertIn(p, g)


    def test_missing_information(self):
        """Checks that entity and abundance functions raise on missing name/identifier"""
        with self.assertRaises(ValueError):
            entity(namespace='test')

        with self.assertRaises(ValueError):
            protein(namespace='test')


class TestGetGraphProperties(unittest.TestCase):
    """The tests in this class check the getting and setting of node properties"""

    def setUp(self):
        self.graph = BELGraph()

    def test_get_qualified_edge(self):
        test_source = protein(namespace='TEST', name='YFG')
        test_target = protein(namespace='TEST', name='YFG2')
        test_evidence = n()
        test_pmid = n()

        test_key = self.graph.add_qualified_edge(
            test_source,
            test_target,
            relation=INCREASES,
            citation=test_pmid,
            evidence=test_evidence,
            annotations={
                'Species': '9606',
                'Confidence': 'Very High'
            },
        )

        citation = self.graph.get_edge_citation(test_source, test_target, test_key)

        self.assertIsNotNone(citation)
        self.assertIsInstance(citation, dict)
        self.assertIn(CITATION_TYPE, citation)
        self.assertEqual(CITATION_TYPE_PUBMED, citation[CITATION_TYPE])
        self.assertIn(CITATION_REFERENCE, citation)
        self.assertEqual(test_pmid, citation[CITATION_REFERENCE])

        evidence = self.graph.get_edge_evidence(test_source, test_target, test_key)

        self.assertIsNotNone(evidence)
        self.assertIsInstance(evidence, string_types)
        self.assertEqual(test_evidence, evidence)

        annotations = self.graph.get_edge_annotations(test_source, test_target, test_key)
        self.assertIsNotNone(annotations)
        self.assertIsInstance(annotations, dict)
        self.assertIn('Species', annotations)
        self.assertIn('9606', annotations['Species'])
        self.assertTrue(annotations['Species']['9606'])
        self.assertIn('Confidence', annotations)
        self.assertIn('Very High', annotations['Confidence'])
        self.assertTrue(annotations['Confidence']['Very High'])

    def test_get_unqualified_edge(self):
        test_source = protein(namespace='TEST', name='YFG')
        test_target = protein(namespace='TEST', name='YFG2')

        key = self.graph.add_unqualified_edge(
            test_source,
            test_target,
            relation=HAS_VARIANT,
        )

        citation = self.graph.get_edge_citation(test_source, test_target, key)
        self.assertIsNone(citation)

        evidence = self.graph.get_edge_evidence(test_source, test_target, key)
        self.assertIsNone(evidence)

        annotations = self.graph.get_edge_annotations(test_source, test_target, key)
        self.assertIsNone(annotations)


    def test_get_node_properties(self):
        test_name = n()
        test_identifier = n()

        node = protein(namespace='TEST', name=test_name, identifier=test_identifier)
        self.graph.add_node_from_data(node)

        self.assertIsNone(self.graph.get_node_description(node))

        test_description = n()
        self.graph.set_node_description(node, test_description)
        self.assertEqual(test_description, self.graph.get_node_description(node))


if __name__ == '__main__':
    unittest.main()
