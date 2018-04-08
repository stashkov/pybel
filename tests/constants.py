# -*- coding: utf-8 -*-

import logging
import os
import tempfile
import unittest
from json import dumps
from pathlib import Path

from requests.compat import urlparse

from pybel import BELGraph
from pybel.constants.bel_graph_keywords import *
from pybel.constants.citation_data_keys import *
from pybel.constants.edge_data_keys import *
from pybel.constants.edge_modifiers import *
from pybel.constants.edge_types import *
from pybel.constants.misc import (
    BEL_DEFAULT_NAMESPACE, FRAUNHOFER_RESOURCES, OPENBEL_ANNOTATION_RESOURCES,
    OPENBEL_NAMESPACE_RESOURCES,
)
from pybel.constants.node_data_keys import *
from pybel.dsl.edges import translocation
from pybel.dsl.nodes import (
    BaseEntity, abundance, bioprocess, complex_abundance, composite_abundance, fragment,
    fusion_range, gene, gene_fusion, gmod, hgvs, hgvs_reference, hgvs_unspecified, mirna, named_complex_abundance,
    pathology, pmod, protein, protein_fusion, reaction, rna, rna_fusion,
)
from pybel.manager import Manager
from pybel.parser.exc import *
from pybel.parser.parse_bel import BelParser
from pybel.utils import hash_edge, subdict_matches

log = logging.getLogger(__name__)
logging.getLogger('pybel.io').setLevel(logging.CRITICAL)
logging.getLogger('pybel.parser').setLevel(logging.CRITICAL)
logging.getLogger('pybel.manager').setLevel(logging.CRITICAL)
logging.getLogger('pybel.resources').setLevel(logging.CRITICAL)

dir_path = os.path.dirname(os.path.realpath(__file__))
owl_dir_path = os.path.join(dir_path, 'owl')
bel_dir_path = os.path.join(dir_path, 'bel')
belns_dir_path = os.path.join(dir_path, 'belns')
belanno_dir_path = os.path.join(dir_path, 'belanno')
beleq_dir_path = os.path.join(dir_path, 'beleq')

test_bel_simple = os.path.join(bel_dir_path, 'test_bel.bel')
test_bel_extensions = os.path.join(bel_dir_path, 'test_bel_owl_extension.bel')
test_bel_slushy = os.path.join(bel_dir_path, 'slushy.bel')
test_bel_thorough = os.path.join(bel_dir_path, 'thorough.bel')
test_bel_isolated = os.path.join(bel_dir_path, 'isolated.bel')
test_bel_misordered = os.path.join(bel_dir_path, 'misordered.bel')
test_bel_no_identifier_valiation = os.path.join(bel_dir_path, 'no_identifier_validation_test.bel')

test_owl_pizza = os.path.join(owl_dir_path, 'pizza_onto.owl')
test_owl_wine = os.path.join(owl_dir_path, 'wine.owl')
test_owl_ado = os.path.join(owl_dir_path, 'ado.owl')

test_an_1 = os.path.join(belanno_dir_path, 'test_an_1.belanno')

test_ns_1 = os.path.join(belns_dir_path, 'test_ns_1.belns')
test_ns_2 = os.path.join(belns_dir_path, 'test_ns_1_updated.belns')
test_ns_nocache = os.path.join(belns_dir_path, 'test_nocache.belns')
test_ns_empty = os.path.join(belns_dir_path, 'test_ns_empty.belns')

test_ns_nocache_path = Path(test_ns_nocache).as_uri()

test_eq_1 = os.path.join(beleq_dir_path, 'disease-ontology.beleq')
test_eq_2 = os.path.join(beleq_dir_path, 'mesh-diseases.beleq')

test_citation_dict = {
    CITATION_TYPE: 'PubMed',
    CITATION_NAME: 'TestName',
    CITATION_REFERENCE: '1235813'
}
SET_CITATION_TEST = 'SET Citation = {{"{type}","{name}","{reference}"}}'.format(**test_citation_dict)
test_evidence_text = 'I read it on Twitter'
test_set_evidence = 'SET Evidence = "{}"'.format(test_evidence_text)

CHEBI_KEYWORD = 'CHEBI'
CHEBI_URL = OPENBEL_NAMESPACE_RESOURCES + 'chebi.belns'
CELL_LINE_KEYWORD = 'CellLine'
CELL_LINE_URL = OPENBEL_ANNOTATION_RESOURCES + 'cell-line.belanno'
HGNC_KEYWORD = 'HGNC'
HGNC_URL = OPENBEL_NAMESPACE_RESOURCES + 'hgnc-human-genes.belns'
MESH_DISEASES_KEYWORD = 'MeSHDisease'
MESH_DISEASES_URL = OPENBEL_ANNOTATION_RESOURCES + "mesh-diseases.belanno"

pizza_iri = 'http://www.lesfleursdunormal.fr/static/_downloads/pizza_onto.owl'
wine_iri = 'http://www.w3.org/TR/2003/PR-owl-guide-20031209/wine'

test_connection = os.environ.get('PYBEL_TEST_CONNECTION')


def update_provenance(control_parser):
    """Sticks provenance in a BEL parser
    
    :param pybel.parser.parse_control.ControlParser bel_parser:
    :return: 
    """
    control_parser.citation.update(test_citation_dict)
    control_parser.evidence = test_evidence_text


def assertHasNode(self, node, graph, **kwargs):
    """A helper function for checking if a node with the given properties is contained within a graph

    :param self: A Test Case
    :type self: unittest.TestCase
    :param BaseEntity node:
    :param pybel.BELGraph graph:
    :param kwargs:
    """
    self.assertIsInstance(node, BaseEntity, msg='not a BaseEntity: {} {}'.format(node.__class__.__name__, node))

    self.assertTrue(graph.has_node(node),
                    msg='{} not found in graph. Actual contents:\n{}'.format(node, '\n'.join(
                        node.as_bel() for node in graph)))
    if kwargs:
        missing = set(kwargs) - set(graph.node[node])
        self.assertFalse(missing, msg="Missing {} in node data".format(', '.join(sorted(missing))))
        self.assertTrue(all(kwarg in graph.node[node] for kwarg in kwargs),
                        msg="Missing kwarg in node data")
        self.assertEqual(kwargs, {k: graph.node[node][k] for k in kwargs},
                         msg="Wrong values in node data")


def any_dict_matches(dict_of_dicts, query_dict):
    """

    :param dict_of_dicts:
    :param query_dict:
    :return:
    """
    return any(
        query_dict == sd
        for sd in dict_of_dicts.values()
    )


def any_subdict_matches(dict_of_dicts, query_dict):
    """Checks if dictionary target_dict matches one of the subdictionaries of a

    :param dict[any,dict] dict_of_dicts: dictionary of dictionaries
    :param dict query_dict: dictionary
    :return: if dictionary target_dict matches one of the subdictionaries of a
    :rtype: bool
    """
    return any(
        subdict_matches(sub_dict, query_dict)
        for sub_dict in dict_of_dicts.values()
    )


def assertHasEdge(self, u, v, graph, **kwargs):
    """A helper function for checking if an edge with the given properties is contained within a graph

    :param unittest.TestCase self: A TestCase
    :param BaseEntity u: source node
    :param BaseEntity v: target node
    :param BELGraph graph: underlying graph
    """
    self.assertIsInstance(u, BaseEntity)
    self.assertIsInstance(v, BaseEntity)

    self.assertIn(u, graph, msg='source not found in graph: {}'.format(u))
    self.assertIn(v, graph, msg='target not found in graph: {}'.format(v))
    self.assertIn(v, graph[u], msg='edge ({}, {}) not in graph'.format(u, v))

    self.assertIn(RELATION, kwargs)
    relation = kwargs[RELATION]

    if relation in unqualified_edges:
        self.assertTrue(graph.has_unqualified_edge(u, v, relation))

    else:
        key = hash_edge(u, v, kwargs)
        self.assertIn(key, graph[u][v],
                      msg='could not find edge from {} to {}:\n{}'.format(u, v,
                                                                          dumps(kwargs, indent=2, sort_keys=True)))

        # equivalent to prior test
        self.assertTrue(graph.has_edge(u, v, key))


class TestGraphMixin(unittest.TestCase):
    def assertHasNode(self, g, n, **kwargs):
        """Helper for asserting node membership
        
        :param pybel.BELGraph g: Graph
        :param BaseEntity n: Node
        :param kwargs: 
        """
        assertHasNode(self, n, g, **kwargs)

    def assertHasEdge(self, g, u, v, **kwargs):
        """Helper for asserting edge membership
        
        :param pybel.BELGraph g: Graph
        :param BaseEntity u: Source node
        :param BaseEntity v: Target node
        :param kwargs: 
        """
        assertHasEdge(self, u, v, g, **kwargs)


class TemporaryCacheMixin(unittest.TestCase):
    def setUp(self):
        self.fd, self.path = tempfile.mkstemp()
        self.connection = 'sqlite:///' + self.path
        log.info('Test generated connection string %s', self.connection)

        self.manager = Manager(connection=self.connection)
        self.manager.create_all()

    def tearDown(self):
        self.manager.session.close()
        os.close(self.fd)
        os.remove(self.path)


class TemporaryCacheClsMixin(unittest.TestCase):
    """Facilitates generating a database in a temporary file on a class-by-class basis"""

    manager = None
    fd = None
    path = None

    @classmethod
    def setUpClass(cls):
        if test_connection:
            cls.connection = test_connection
        else:
            cls.fd, cls.path = tempfile.mkstemp()
            cls.connection = 'sqlite:///' + cls.path
            log.info('Test generated connection string %s', cls.connection)

        cls.manager = Manager(connection=cls.connection)
        cls.manager.create_all()

    @classmethod
    def tearDownClass(cls):
        cls.manager.session.close()

        if not test_connection:
            os.close(cls.fd)
            os.remove(cls.path)
        else:
            cls.manager.drop_all()


class FleetingTemporaryCacheMixin(TemporaryCacheClsMixin):
    """This class makes a manager available for the entire existence of the class but deletes everything that gets
    stuck in it after each test"""

    def setUp(self):
        super(FleetingTemporaryCacheMixin, self).setUp()

        self.manager.drop_networks()
        self.manager.drop_edges()
        self.manager.drop_nodes()
        self.manager.drop_namespaces()
        self.manager.drop_annotations()


class TestTokenParserBase(unittest.TestCase):
    graph = None
    parser = None

    @classmethod
    def setUpClass(cls):
        """Makes a parser once"""
        cls.graph = BELGraph()
        cls.parser = BelParser(cls.graph, autostreamline=False)

    def setUp(self):
        """Clears the parser before doing anything"""
        self.parser.clear()

    def assertHasNode(self, member):
        """Asserts a node exists in the raph"""
        assertHasNode(self, member, self.graph, **member)

    def assertHasEdge(self, u, v, **kwargs):
        """Asserts an edge exists between the two nodes"""
        self.assertIn(RELATION, kwargs)
        assertHasEdge(self, u, v, self.graph, **kwargs)

    def assertNumberNodes(self, n):
        """
        :param int n: The number of expected nodes
        """
        self.assertEqual(n, self.graph.number_of_nodes(),
                         msg='Wrong number of nodes:\n{}'.format('\n'.join(map(str, self.graph))))

    def assertNumberEdges(self, n):
        """
        :param int n: The number of expected edges
        """
        self.assertEqual(n, self.graph.number_of_edges(),
                         msg='Wrong number of edges:\n{}'.format('\n'.join(
                             str((str(u), str(v), str(key)[:10], str(data))) for u, v, key, data in
                             self.graph.edges(keys=True, data=True))))

    def add_default_provenance(self):
        """Adds the default citation and evidence to the parser"""
        update_provenance(self.parser.control_parser)

    def help_test_parent_in_graph(self, node):
        """Checks that the parent is also in the graph

        :param pybel.dsl.nodes.CentralDogma node:
        """
        parent = node.get_parent()
        self.assertIsNotNone(parent)
        self.assertHasNode(parent)
        self.assertHasEdge(parent, node, relation=HAS_VARIANT)


expected_test_simple_metadata = {
    METADATA_NAME: "PyBEL Test Simple",
    METADATA_DESCRIPTION: "Made for testing PyBEL parsing",
    METADATA_VERSION: "1.6.0",
    METADATA_COPYRIGHT: "Copyright (c) Charles Tapley Hoyt. All Rights Reserved.",
    METADATA_AUTHORS: "Charles Tapley Hoyt",
    METADATA_LICENSES: "WTF License",
    METADATA_CONTACT: "charles.hoyt@scai.fraunhofer.de",
    METADATA_PROJECT: 'PyBEL Testing',
}

expected_test_thorough_metadata = {
    METADATA_NAME: "PyBEL Test Thorough",
    METADATA_DESCRIPTION: "Statements made up to contain many conceivable variants of nodes from BEL",
    METADATA_VERSION: "1.0.0",
    METADATA_COPYRIGHT: "Copyright (c) Charles Tapley Hoyt. All Rights Reserved.",
    METADATA_AUTHORS: "Charles Tapley Hoyt",
    METADATA_LICENSES: "WTF License",
    METADATA_CONTACT: "charles.hoyt@scai.fraunhofer.de"
}

expected_test_bel_4_metadata = {
    METADATA_NAME: "PyBEL Test OWL Extension",
    METADATA_DESCRIPTION: "Tests the use of OWL ontologies as namespaces",
    METADATA_VERSION: "1.6.0",
    METADATA_COPYRIGHT: "Copyright (c) Charles Tapley Hoyt. All Rights Reserved.",
    METADATA_AUTHORS: "Charles Tapley Hoyt",
    METADATA_LICENSES: "WTF License",
    METADATA_CONTACT: "charles.hoyt@scai.fraunhofer.de"
}

expected_test_slushy_metadata = {
    METADATA_NAME: "Worst. BEL Document. Ever.",
    METADATA_DESCRIPTION: "This document outlines all of the evil and awful work that is possible during BEL curation",
    METADATA_VERSION: "0.0",
    METADATA_AUTHORS: "Charles Tapley Hoyt",
    METADATA_LICENSES: "WTF License",
}


def get_uri_name(url):
    """Gets the file name from the end of the URL. Only useful for PyBEL's testing though since it looks specifically
    if the file is from the weird owncloud resources distributed by Fraunhofer"""
    url_parsed = urlparse(url)

    if url.startswith(FRAUNHOFER_RESOURCES):
        return url_parsed.query.split('=')[-1]
    else:
        url_parts = url_parsed.path.split('/')
        return url_parts[-1]


def help_check_hgnc(self, namespace_dict):
    """
    :param unittest.TestCase self: A unittest test case
    :param namespace_dict: A dictionary of values for testing a namespace
    """
    self.assertIn(HGNC_KEYWORD, namespace_dict)

    self.assertIn('MHS2', namespace_dict[HGNC_KEYWORD])
    self.assertEqual(set('G'), set(namespace_dict[HGNC_KEYWORD]['MHS2']))

    self.assertIn('MIATNB', namespace_dict[HGNC_KEYWORD])
    self.assertEqual(set('GR'), set(namespace_dict[HGNC_KEYWORD]['MIATNB']))

    self.assertIn('MIA', namespace_dict[HGNC_KEYWORD])
    self.assertEqual(set('GRP'), set(namespace_dict[HGNC_KEYWORD]['MIA']))


AKT1 = protein('HGNC', 'AKT1')
EGFR = protein('HGNC', 'EGFR')
FADD = protein('HGNC', 'FADD')
CASP8 = protein('HGNC', 'CASP8')
cftr = protein('HGNC', 'CFTR')
mia = protein('HGNC', 'MIA')

il6 = protein('HGNC', 'IL6')

adgrb1 = protein(namespace='HGNC', name='ADGRB1')
adgrb2 = protein(namespace='HGNC', name='ADGRB2')
adgrb_complex = complex_abundance([adgrb1, adgrb2])
achlorhydria = pathology(namespace='MESHD', name='Achlorhydria')

akt1_gene = gene('HGNC', 'AKT1')
akt1_rna = rna('HGNC', 'AKT1')
oxygen_atom = abundance('CHEBI', 'oxygen atom')
akt_methylated = gene('HGNC', 'AKT1', variants=gmod('Me'))
akt1_phe_508_del = gene('HGNC', 'AKT1', variants=hgvs('p.Phe508del'))

cftr_protein_unspecified_variant = protein('HGNC', 'CFTR', variants=hgvs_unspecified())
cftr_protein_phe_508_del = protein('HGNC', 'CFTR', variants=hgvs('p.Phe508del'))
adenocarcinoma = pathology('MESHD', 'Adenocarcinoma')
interleukin_23_complex = named_complex_abundance('GOCC', 'interleukin-23 complex')

hydrogen_peroxide = abundance('CHEBI', 'hydrogen peroxide')

tmprss2_erg_gene_fusion = gene_fusion(
    partner5p=gene('HGNC', 'TMPRSS2'),
    range5p=fusion_range('c', 1, 79),
    partner3p=gene('HGNC', 'ERG'),
    range3p=fusion_range('c', 312, 5034)
)

bcr_jak2_gene_fusion = gene_fusion(
    partner5p=gene('HGNC', 'BCR'),
    range5p=fusion_range('c', '?', 1875),
    partner3p=gene('HGNC', 'JAK2'),
    range3p=fusion_range('c', 2626, '?')
)

chchd4_aifm1_gene_fusion = gene_fusion(
    partner5p=gene('HGNC', 'CHCHD4'),
    partner3p=gene('HGNC', 'AIFM1')
)

tmprss2_erg_protein_fusion = protein_fusion(
    partner5p=protein('HGNC', 'TMPRSS2'),
    range5p=fusion_range('p', 1, 79),
    partner3p=protein('HGNC', 'ERG'),
    range3p=fusion_range('p', 312, 5034)
)

bcr_jak2_protein_fusion = protein_fusion(
    partner5p=protein('HGNC', 'BCR'),
    range5p=fusion_range('p', '?', 1875),
    partner3p=protein('HGNC', 'JAK2'),
    range3p=fusion_range('p', 2626, '?')
)

chchd4_aifm1_protein_fusion = protein_fusion(
    protein('HGNC', 'CHCHD4'),
    protein('HGNC', 'AIFM1')
)

bcr_jak2_rna_fusion = rna_fusion(
    partner5p=rna('HGNC', 'BCR'),
    range5p=fusion_range('r', '?', 1875),
    partner3p=rna('HGNC', 'JAK2'),
    range3p=fusion_range('r', 2626, '?')
)

chchd4_aifm1_rna_fusion = rna_fusion(
    partner5p=rna('HGNC', 'CHCHD4'),
    partner3p=rna('HGNC', 'AIFM1')
)

tmprss2_erg_rna_fusion = rna_fusion(
    partner5p=rna('HGNC', 'TMPRSS2'),
    range5p=fusion_range('r', 1, 79),
    partner3p=rna('HGNC', 'ERG'),
    range3p=fusion_range('r', 312, 5034)
)
tmprss2_erg_rna_fusion_unspecified = rna_fusion(
    partner5p=rna('HGNC', 'TMPRSS2'),
    partner3p=rna('HGNC', 'ERG')
)

BEL_THOROUGH_NODES = {
    oxygen_atom,
    tmprss2_erg_rna_fusion,
    tmprss2_erg_rna_fusion_unspecified,
    akt_methylated,
    bcr_jak2_rna_fusion,
    chchd4_aifm1_rna_fusion,
    akt1_gene,
    akt1_phe_508_del,
    AKT1,
    gene('HGNC', 'AKT1', variants=hgvs('c.308G>A')),
    tmprss2_erg_gene_fusion,
    gene('HGNC', 'AKT1', variants=[hgvs('c.1521_1523delCTT'), hgvs('c.308G>A'), hgvs('p.Phe508del')]),
    mirna('HGNC', 'MIR21'),
    bcr_jak2_gene_fusion,
    gene('HGNC', 'CFTR', variants=hgvs('c.1521_1523delCTT')),
    gene('HGNC', 'CFTR'),
    gene('HGNC', 'CFTR', variants=hgvs('g.117199646_117199648delCTT')),
    gene('HGNC', 'CFTR', variants=hgvs('c.1521_1523delCTT')),
    protein('HGNC', 'AKT1', variants=pmod('Ph', 'Ser', 473)),
    mirna('HGNC', 'MIR21', variants=hgvs('p.Phe508del')),
    protein('HGNC', 'AKT1', variants=hgvs('p.C40*')),
    protein('HGNC', 'AKT1', variants=[hgvs('p.Ala127Tyr'), pmod('Ph', 'Ser')]),
    chchd4_aifm1_gene_fusion,
    tmprss2_erg_protein_fusion,
    protein('HGNC', 'AKT1', variants=hgvs('p.Arg1851*')),
    bcr_jak2_protein_fusion,
    protein('HGNC', 'AKT1', variants=hgvs('p.40*')),
    chchd4_aifm1_protein_fusion,
    protein('HGNC', 'CFTR', variants=hgvs_reference()),
    cftr,
    EGFR,
    cftr_protein_unspecified_variant,
    adenocarcinoma,
    cftr_protein_phe_508_del,
    protein('HGNC', 'MIA', variants=fragment(5, 20)),
    mia,
    interleukin_23_complex,
    protein('HGNC', 'MIA', variants=fragment(1, '?')),
    protein('HGNC', 'MIA', variants=fragment()),
    protein('HGNC', 'MIA', variants=fragment(description='55kD')),
    protein('HGNC', 'CFTR', variants=hgvs('p.Gly576Ala')),
    akt1_rna,
    rna('HGNC', 'AKT1', variants=[hgvs('c.1521_1523delCTT'), hgvs('p.Phe508del')]),
    gene('HGNC', 'NCF1'),
    complex_abundance([
        gene('HGNC', 'NCF1'),
        protein('HGNC', 'HBP1')
    ]),
    protein('HGNC', 'HBP1'),
    complex_abundance([protein('HGNC', 'FOS'), protein('HGNC', 'JUN')]),
    protein('HGNC', 'FOS'),
    protein('HGNC', 'JUN'),
    rna('HGNC', 'CFTR', variants=hgvs('r.1521_1523delcuu')),
    rna('HGNC', 'CFTR'),
    rna('HGNC', 'CFTR', variants=hgvs('r.1653_1655delcuu')),
    composite_abundance([
        interleukin_23_complex,
        il6
    ]),
    il6,
    bioprocess('GOBP', 'cell cycle arrest'),
    hydrogen_peroxide,
    protein('HGNC', 'CAT'),
    gene('HGNC', 'CAT'),
    protein('HGNC', 'HMGCR'),
    bioprocess('GOBP', 'cholesterol biosynthetic process'),
    gene('HGNC', 'APP', variants=hgvs('c.275341G>C')),
    gene('HGNC', 'APP'),
    pathology('MESHD', 'Alzheimer Disease'),
    complex_abundance([protein('HGNC', 'F3'), protein('HGNC', 'F7')]),
    protein('HGNC', 'F3'),
    protein('HGNC', 'F7'),
    protein('HGNC', 'F9'),
    protein('HGNC', 'GSK3B', variants=pmod('Ph', 'Ser', 9)),
    protein('HGNC', 'GSK3B'),
    pathology('MESHD', 'Psoriasis'),
    pathology('MESHD', 'Skin Diseases'),
    reaction(
        reactants=[
            abundance('CHEBI', '(3S)-3-hydroxy-3-methylglutaryl-CoA'),
            abundance('CHEBI', 'NADPH'),
            abundance('CHEBI', 'hydron')
        ], products=[
            abundance('CHEBI', 'NADP(+)'),
            abundance('CHEBI', 'mevalonate')
        ]),
    abundance('CHEBI', '(3S)-3-hydroxy-3-methylglutaryl-CoA'),
    abundance('CHEBI', 'NADPH'),
    abundance('CHEBI', 'hydron'),
    abundance('CHEBI', 'mevalonate'),
    abundance('CHEBI', 'NADP(+)'),
    abundance('CHEBI', 'nitric oxide'),
    complex_abundance([protein('HGNC', 'ITGAV'), protein('HGNC', 'ITGB3')]),
    protein('HGNC', 'ITGAV'),
    protein('HGNC', 'ITGB3'),
    protein('HGNC', 'FADD'),
    abundance('TESTNS2', 'Abeta_42'),
    protein('TESTNS2', 'GSK3 Family'),
    protein('HGNC', 'PRKCA'),
    protein('HGNC', 'CDK5'),
    protein('HGNC', 'CASP8'),
    protein('HGNC', 'AKT1', variants=pmod(namespace='TESTNS2', name='PhosRes', code='Ser', position=473)),
    protein('HGNC', 'HRAS', variants=pmod('Palm')),
    bioprocess('GOBP', 'apoptotic process'),
    composite_abundance([
        abundance('TESTNS2', 'Abeta_42'),
        protein('HGNC', 'CASP8'),
        protein('HGNC', 'FADD')
    ]),
    reaction(
        reactants=protein('HGNC', 'CDK5R1'),
        products=protein('HGNC', 'CDK5')
    ),
    protein('HGNC', 'PRKCB'),
    named_complex_abundance('TESTNS2', 'AP-1 Complex'),
    protein('HGNC', 'PRKCE'),
    protein('HGNC', 'PRKCD'),
    protein('TESTNS2', 'CAPN Family'),
    gene('TESTNS2', 'AKT1 ortholog'),
    protein('HGNC', 'HRAS'),
    protein('HGNC', 'CDK5R1'),
    protein('TESTNS2', 'PRKC'),
    bioprocess('GOBP', 'neuron apoptotic process'),
    protein('HGNC', 'MAPT', variants=pmod('Ph')),
    protein('HGNC', 'MAPT'),
    gene('HGNC', 'ARRDC2'),
    gene('HGNC', 'ARRDC3'),
    gene('dbSNP', 'rs123456')
}

citation_1 = {
    CITATION_TYPE: 'PubMed',
    CITATION_NAME: 'That one article from last week',
    CITATION_REFERENCE: '123455'
}

citation_2 = {
    CITATION_TYPE: 'PubMed',
    CITATION_NAME: 'That one article from last week #2',
    CITATION_REFERENCE: '123456'
}

evidence_1 = "Evidence 1"
dummy_evidence = 'These are mostly made up'

bel_thorough_graph = BELGraph()
bel_thorough_graph.add_qualified_edge(
    oxygen_atom,
    akt_methylated,
    relation=INCREASES,
    evidence=dummy_evidence,
    citation=citation_1,
    annotations={
        'TESTAN1': {'1': True, '2': True},
        'TestRegex': {'9000': True}
    }
)
bel_thorough_graph.add_has_variant(akt1_gene, akt_methylated)
bel_thorough_graph.add_qualified_edge(
    akt1_gene,
    oxygen_atom,
    evidence=dummy_evidence,
    citation=citation_1,
    relation=DECREASES,
    subject_modifier={LOCATION: {NAMESPACE: 'GOCC', NAME: 'intracellular'}},
    object_modifier={LOCATION: {NAMESPACE: 'GOCC', NAME: 'intracellular'}},
)
bel_thorough_graph.add_has_variant(akt1_gene, gene('HGNC', 'AKT1', variants=hgvs('p.Phe508del')))
bel_thorough_graph.add_has_variant(akt1_gene, gene('HGNC', 'AKT1', variants=hgvs('c.308G>A')))
bel_thorough_graph.add_has_variant(akt1_gene, gene('HGNC', 'AKT1',
                                                   variants=[hgvs('c.1521_1523delCTT'), hgvs('c.308G>A'),
                                                             hgvs('p.Phe508del')]))
bel_thorough_graph.add_transcription(akt1_gene, akt1_rna)
BEL_THOROUGH_EDGES = [
    (gene('HGNC', 'AKT1', variants=hgvs('p.Phe508del')), AKT1, {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: DIRECTLY_DECREASES,
    }),
    (AKT1, protein('HGNC', 'AKT1', variants=pmod('Ph', 'Ser', 473)), {
        RELATION: HAS_VARIANT,
    }),
    (AKT1, protein('HGNC', 'AKT1', variants=hgvs('p.C40*')), {
        RELATION: HAS_VARIANT,
    }),
    (AKT1,
     protein('HGNC', 'AKT1', variants=[hgvs('p.Ala127Tyr'), pmod('Ph', 'Ser')]), {
         RELATION: HAS_VARIANT,
     }),
    (AKT1,
     protein('HGNC', 'AKT1', variants=[hgvs('p.Ala127Tyr'), pmod('Ph', 'Ser')]), {
         EVIDENCE: dummy_evidence,
         CITATION: citation_1,
         RELATION: DIRECTLY_DECREASES,
         SUBJECT: {LOCATION: {NAMESPACE: 'GOCC', NAME: 'intracellular'}},
         OBJECT: {LOCATION: {NAMESPACE: 'GOCC', NAME: 'intracellular'}},
     }),
    (AKT1, protein('HGNC', 'AKT1', variants=hgvs('p.Arg1851*')), {
        RELATION: HAS_VARIANT,
    }),
    (AKT1, protein('HGNC', 'AKT1', variants=hgvs('p.40*')), {
        RELATION: HAS_VARIANT,
    }),
    (AKT1, protein('HGNC', 'MIA', variants=fragment()), {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: INCREASES,
        SUBJECT: {MODIFIER: DEGRADATION},
    }),
    (AKT1, protein('HGNC', 'CFTR', variants=hgvs('p.Gly576Ala')), {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: INCREASES,
    }),
    (AKT1, rna('HGNC', 'CFTR', variants=hgvs('r.1521_1523delcuu')), {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: INCREASES,
        SUBJECT: {MODIFIER: ACTIVITY, EFFECT: {NAMESPACE: BEL_DEFAULT_NAMESPACE, NAME: 'kin'}},
    }),
    (AKT1, rna('HGNC', 'CFTR', variants=hgvs('r.1653_1655delcuu')), {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: INCREASES,
        SUBJECT: {MODIFIER: ACTIVITY},
    }),
    (AKT1, EGFR, {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: INCREASES,
        SUBJECT: {
            MODIFIER: ACTIVITY,
            EFFECT: {
                NAMESPACE: BEL_DEFAULT_NAMESPACE,
                NAME: 'cat'
            }
        },
        OBJECT: {MODIFIER: DEGRADATION},
    }),
    (AKT1, EGFR, {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: INCREASES,
        SUBJECT: {
            MODIFIER: ACTIVITY,
            EFFECT: {NAME: 'kin', NAMESPACE: BEL_DEFAULT_NAMESPACE}
        },
        OBJECT: translocation(
            {NAMESPACE: 'GOCC', NAME: 'intracellular'},
            {NAMESPACE: 'GOCC', NAME: 'extracellular space'}
        ),
    }),
    (gene('HGNC', 'AKT1', variants=hgvs('c.308G>A')), tmprss2_erg_gene_fusion, {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: CAUSES_NO_CHANGE,
    }),
    (gene('HGNC', 'AKT1', variants=hgvs('c.308G>A')),
     gene('HGNC', 'AKT1', variants=[hgvs('c.1521_1523delCTT'), hgvs('c.308G>A'), hgvs('p.Phe508del')]), {
         EVIDENCE: dummy_evidence,
         CITATION: citation_1,
         RELATION: INCREASES,
         SUBJECT: {LOCATION: {NAMESPACE: 'GOCC', NAME: 'intracellular'}},
     }),
    (mirna('HGNC', 'MIR21'), bcr_jak2_gene_fusion,
     {
         EVIDENCE: dummy_evidence,
         CITATION: citation_1,
         RELATION: DIRECTLY_INCREASES,
     }),
    (mirna('HGNC', 'MIR21'), protein('HGNC', 'AKT1', variants=pmod('Ph', 'Ser', 473)),
     {
         EVIDENCE: dummy_evidence,
         CITATION: citation_1,
         RELATION: DECREASES,
         SUBJECT: {LOCATION: {NAMESPACE: 'GOCC', NAME: 'intracellular'}},
     }),
    (mirna('HGNC', 'MIR21'), mirna('HGNC', 'MIR21', variants=hgvs('p.Phe508del')), {
        RELATION: HAS_VARIANT,
    }),
    (gene('HGNC', 'CFTR', variants=hgvs('c.1521_1523delCTT')), AKT1,
     {
         EVIDENCE: dummy_evidence,
         CITATION: citation_1,
         RELATION: INCREASES,
         OBJECT: {MODIFIER: DEGRADATION},
     }),
    (gene('HGNC', 'CFTR'), gene('HGNC', 'CFTR', variants=hgvs('c.1521_1523delCTT')), {
        RELATION: HAS_VARIANT,
    }),
    (gene('HGNC', 'CFTR'), gene('HGNC', 'CFTR', variants=hgvs('g.117199646_117199648delCTT')), {
        RELATION: HAS_VARIANT,
    }),
    (gene('HGNC', 'CFTR'), gene('HGNC', 'CFTR', variants=hgvs('c.1521_1523delCTT')), {
        RELATION: HAS_VARIANT,
    }),
    (gene('HGNC', 'CFTR', variants=hgvs('g.117199646_117199648delCTT')),
     gene('HGNC', 'CFTR', variants=hgvs('c.1521_1523delCTT')), {
         EVIDENCE: dummy_evidence,
         CITATION: citation_1,
         RELATION: INCREASES,
     }),
    (mirna('HGNC', 'MIR21', variants=hgvs('p.Phe508del')), protein('HGNC', 'AKT1', variants=hgvs('p.C40*')),
     {
         EVIDENCE: dummy_evidence,
         CITATION: citation_1,
         RELATION: INCREASES,
         SUBJECT: {LOCATION: {NAMESPACE: 'GOCC', NAME: 'intracellular'}},
     }),
    (chchd4_aifm1_gene_fusion, tmprss2_erg_protein_fusion, {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: INCREASES,
    }),
    (protein('HGNC', 'AKT1', variants=hgvs('p.Arg1851*')), bcr_jak2_protein_fusion, {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: INCREASES,
    }),
    (protein('HGNC', 'AKT1', variants=hgvs('p.40*')), chchd4_aifm1_protein_fusion,
     {
         EVIDENCE: dummy_evidence,
         CITATION: citation_1,
         RELATION: INCREASES,
     }),
    (protein('HGNC', 'CFTR', variants=hgvs_reference()), EGFR, {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: INCREASES,
        OBJECT: {
            MODIFIER: TRANSLOCATION,
            EFFECT: {
                FROM_LOC: {NAMESPACE: 'GOCC', NAME: 'intracellular'},
                TO_LOC: {NAMESPACE: 'GOCC', NAME: 'cell surface'}
            }
        },
    }),
    (cftr, protein('HGNC', 'CFTR', variants=hgvs('=')), {
        RELATION: HAS_VARIANT,
    }),
    (cftr, protein('HGNC', 'CFTR', variants=hgvs('?')), {
        RELATION: HAS_VARIANT,
    }),
    (cftr, protein('HGNC', 'CFTR', variants=hgvs('p.Phe508del')), {
        RELATION: HAS_VARIANT,
    }),
    (cftr, protein('HGNC', 'CFTR', variants=hgvs('p.Gly576Ala')), {
        RELATION: HAS_VARIANT,
    }),
    (mia, protein('HGNC', 'MIA', variants=fragment(5, 20)), {
        RELATION: HAS_VARIANT,
    }),
    (mia, protein('HGNC', 'MIA', variants=fragment(1, '?')), {
        RELATION: HAS_VARIANT,
    }),
    (mia, protein('HGNC', 'MIA', variants=fragment()), {
        RELATION: HAS_VARIANT,
    }),
    (mia, protein('HGNC', 'MIA', variants=fragment(description='55kD')), {
        RELATION: HAS_VARIANT,
    }),
    (akt1_rna, rna('HGNC', 'AKT1', variants=[hgvs('c.1521_1523delCTT'), hgvs('p.Phe508del')]), {
        RELATION: HAS_VARIANT,
    }),
    (akt1_rna, AKT1, {
        RELATION: TRANSLATED_TO,
    }),
    (gene('HGNC', 'APP'), gene('HGNC', 'APP', variants=hgvs('c.275341G>C')), {
        RELATION: HAS_VARIANT,
    }),
    (complex_abundance([protein('HGNC', 'F3'), protein('HGNC', 'F7')]), protein('HGNC', 'F3'), {
        RELATION: HAS_COMPONENT,
    }),
    (complex_abundance([protein('HGNC', 'F3'), protein('HGNC', 'F7')]), protein('HGNC', 'F7'), {
        RELATION: HAS_COMPONENT,
    }),
    (protein('HGNC', 'GSK3B'), protein('HGNC', 'GSK3B', variants=pmod('Ph', 'Ser', 9)), {
        RELATION: HAS_VARIANT,
    }),
    (pathology('MESHD', 'Psoriasis'), pathology('MESHD', 'Skin Diseases'), {
        RELATION: IS_A,
    }),

    (complex_abundance([gene('HGNC', 'NCF1'), protein('HGNC', 'HBP1')]), protein('HGNC', 'HBP1'), {
        RELATION: HAS_COMPONENT,
    }),
    (complex_abundance([gene('HGNC', 'NCF1'), protein('HGNC', 'HBP1')]), gene('HGNC', 'NCF1'), {
        RELATION: HAS_COMPONENT,
    }),

    (complex_abundance([protein('HGNC', 'FOS'), protein('HGNC', 'JUN')]), protein('HGNC', 'FOS'), {
        RELATION: HAS_COMPONENT,
    }),
    (complex_abundance([protein('HGNC', 'FOS'), protein('HGNC', 'JUN')]), protein('HGNC', 'JUN'), {
        RELATION: HAS_COMPONENT,
    }),
    (rna('HGNC', 'CFTR'), rna('HGNC', 'CFTR', variants=hgvs('r.1521_1523delcuu')), {
        RELATION: HAS_VARIANT,
    }),
    (rna('HGNC', 'CFTR'), rna('HGNC', 'CFTR', variants=hgvs('r.1653_1655delcuu')), {
        RELATION: HAS_VARIANT,
    }),
    (composite_abundance([interleukin_23_complex, il6]), il6, {
        RELATION: HAS_COMPONENT,
    }),
    (composite_abundance([interleukin_23_complex, il6]), interleukin_23_complex, {
        RELATION: HAS_COMPONENT,
    }),
    (protein('HGNC', 'CFTR', variants=hgvs('?')), pathology('MESHD', 'Adenocarcinoma'), {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: INCREASES,
    }),
    (rna('HGNC', 'AKT1', variants=[hgvs('c.1521_1523delCTT'), hgvs('p.Phe508del')]),
     tmprss2_erg_rna_fusion, {
         EVIDENCE: dummy_evidence,
         CITATION: citation_1,
         RELATION: DIRECTLY_INCREASES,
     }),
    (rna_fusion(rna('HGNC', 'TMPRSS2'), rna('HGNC', 'ERG')),
     complex_abundance([gene('HGNC', 'NCF1'), protein('HGNC', 'HBP1')]), {
         EVIDENCE: dummy_evidence,
         CITATION: citation_1,
         RELATION: INCREASES,
     }),
    (protein('HGNC', 'MIA', variants=fragment(5, 20)), named_complex_abundance('GOCC', 'interleukin-23 complex'), {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: INCREASES,
        OBJECT: {
            MODIFIER: TRANSLOCATION,
            EFFECT: {
                FROM_LOC: {NAMESPACE: 'GOCC', NAME: 'intracellular'},
                TO_LOC: {NAMESPACE: 'GOCC', NAME: 'extracellular space'}
            }
        },
    }),
    (protein('HGNC', 'MIA', variants=fragment(1, '?')), EGFR, {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: INCREASES,
        OBJECT: {
            MODIFIER: TRANSLOCATION,
            EFFECT: {
                FROM_LOC: {NAMESPACE: 'GOCC', NAME: 'cell surface'},
                TO_LOC: {NAMESPACE: 'GOCC', NAME: 'endosome'}
            }
        },
    }),
    (akt1_rna, EGFR, {
        EVIDENCE: dummy_evidence,
        CITATION: citation_1,
        RELATION: INCREASES,
        OBJECT: {
            MODIFIER: TRANSLOCATION, EFFECT: {
                FROM_LOC: {NAMESPACE: 'GOCC', NAME: 'cell surface'},
                TO_LOC: {NAMESPACE: 'GOCC', NAME: 'endosome'}
            }
        },
    }),
    (rna_fusion(rna('HGNC', 'CHCHD4'), rna('HGNC', 'AIFM1'), ),
     complex_abundance([protein('HGNC', 'FOS'), protein('HGNC', 'JUN')]),
     {
         EVIDENCE: dummy_evidence,
         CITATION: citation_1,
         RELATION: INCREASES,
     }),
    (composite_abundance([interleukin_23_complex, il6]),
     bioprocess('GOBP', 'cell cycle arrest'), {
         EVIDENCE: dummy_evidence,
         CITATION: citation_1,
         RELATION: DECREASES,
     }),
    (protein('HGNC', 'CAT'), hydrogen_peroxide, {
        EVIDENCE: 'These were all explicitly stated in the BEL 2.0 Specification',
        CITATION: citation_2,
        RELATION: DIRECTLY_DECREASES,
        SUBJECT: {LOCATION: {NAMESPACE: 'GOCC', NAME: 'intracellular'}},
    }),
    (gene('HGNC', 'CAT'), hydrogen_peroxide, {
        EVIDENCE: 'These were all explicitly stated in the BEL 2.0 Specification',
        CITATION: citation_2,
        RELATION: DIRECTLY_DECREASES,
        SUBJECT: {LOCATION: {NAMESPACE: 'GOCC', NAME: 'intracellular'}},
    }),
    (protein('HGNC', 'HMGCR'), bioprocess('GOBP', 'cholesterol biosynthetic process'), {
        EVIDENCE: 'These were all explicitly stated in the BEL 2.0 Specification',
        CITATION: citation_2,
        RELATION: RATE_LIMITING_STEP_OF,
        SUBJECT: {MODIFIER: ACTIVITY, EFFECT: {NAMESPACE: BEL_DEFAULT_NAMESPACE, NAME: 'cat'}},
    }),
    (gene('HGNC', 'APP', variants=hgvs('c.275341G>C')), pathology('MESHD', 'Alzheimer Disease'),
     {
         EVIDENCE: 'These were all explicitly stated in the BEL 2.0 Specification',
         CITATION: citation_2,
         RELATION: CAUSES_NO_CHANGE,
     }),

    (complex_abundance([protein('HGNC', 'F3'), protein('HGNC', 'F7')]), protein('HGNC', 'F9'), {
        EVIDENCE: 'These were all explicitly stated in the BEL 2.0 Specification',
        CITATION: citation_2,
        RELATION: REGULATES,
        SUBJECT: {MODIFIER: ACTIVITY, EFFECT: {NAME: 'pep', NAMESPACE: BEL_DEFAULT_NAMESPACE}},
        OBJECT: {MODIFIER: ACTIVITY, EFFECT: {NAME: 'pep', NAMESPACE: BEL_DEFAULT_NAMESPACE}},
    }),
    (protein('HGNC', 'GSK3B', variants=pmod('Ph', 'Ser', 9)), protein('HGNC', 'GSK3B'), {
        EVIDENCE: 'These were all explicitly stated in the BEL 2.0 Specification',
        CITATION: citation_2,
        RELATION: POSITIVE_CORRELATION,
        OBJECT: {MODIFIER: ACTIVITY, EFFECT: {NAMESPACE: BEL_DEFAULT_NAMESPACE, NAME: 'kin'}},
    }),

    (protein('HGNC', 'GSK3B'), protein('HGNC', 'GSK3B', variants=pmod('Ph', 'Ser', 9)), {
        EVIDENCE: 'These were all explicitly stated in the BEL 2.0 Specification',
        CITATION: citation_2,
        RELATION: POSITIVE_CORRELATION,
        SUBJECT: {MODIFIER: ACTIVITY, EFFECT: {NAMESPACE: BEL_DEFAULT_NAMESPACE, NAME: 'kin'}},
    }),

    (reaction(
        reactants=(
            abundance('CHEBI', '(3S)-3-hydroxy-3-methylglutaryl-CoA'),
            abundance('CHEBI', 'NADPH'),
            abundance('CHEBI', 'hydron')
        ),
        products=(
            abundance('CHEBI', 'NADP(+)'),
            abundance('CHEBI', 'mevalonate')
        )
    ),
     abundance('CHEBI', '(3S)-3-hydroxy-3-methylglutaryl-CoA'), {
         RELATION: HAS_REACTANT,
     }),
    (reaction(
        reactants=(
            abundance('CHEBI', '(3S)-3-hydroxy-3-methylglutaryl-CoA'),
            abundance('CHEBI', 'NADPH'),
            abundance('CHEBI', 'hydron')
        ),
        products=(
            abundance('CHEBI', 'NADP(+)'),
            abundance('CHEBI', 'mevalonate')
        )
    ),
     abundance('CHEBI', 'NADPH'), {
         RELATION: HAS_REACTANT,
     }),
    (
        reaction(
            reactants=(
                abundance('CHEBI', '(3S)-3-hydroxy-3-methylglutaryl-CoA'),
                abundance('CHEBI', 'NADPH'),
                abundance('CHEBI', 'hydron')
            ),
            products=(
                abundance('CHEBI', 'NADP(+)'),
                abundance('CHEBI', 'mevalonate')
            )
        ),
        abundance('CHEBI', 'hydron'), {
            RELATION: HAS_REACTANT,
        }),
    (reaction(
        reactants=(
            abundance('CHEBI', '(3S)-3-hydroxy-3-methylglutaryl-CoA'),
            abundance('CHEBI', 'NADPH'),
            abundance('CHEBI', 'hydron')
        ),
        products=(
            abundance('CHEBI', 'NADP(+)'),
            abundance('CHEBI', 'mevalonate'))
    ),
     abundance('CHEBI', 'mevalonate'), {
         RELATION: HAS_PRODUCT,
     }),
    (reaction(
        reactants=(
            abundance('CHEBI', '(3S)-3-hydroxy-3-methylglutaryl-CoA'),
            abundance('CHEBI', 'NADPH'),
            abundance('CHEBI', 'hydron')
        ),
        products=(
            abundance('CHEBI', 'NADP(+)'),
            abundance('CHEBI', 'mevalonate')
        )
    ),
     abundance('CHEBI', 'NADP(+)'), {
         RELATION: HAS_PRODUCT,
     }),
    (reaction(
        reactants=(
            abundance('CHEBI', '(3S)-3-hydroxy-3-methylglutaryl-CoA'),
            abundance('CHEBI', 'NADPH'),
            abundance('CHEBI', 'hydron')
        ),
        products=(
            abundance('CHEBI', 'NADP(+)'),
            abundance('CHEBI', 'mevalonate')
        )
    ),
     bioprocess('GOBP', 'cholesterol biosynthetic process'),
     {
         EVIDENCE: 'These were all explicitly stated in the BEL 2.0 Specification',
         CITATION: citation_2,
         RELATION: SUBPROCESS_OF,
     }),
    (abundance('CHEBI', 'nitric oxide'),
     complex_abundance([protein('HGNC', 'ITGAV'), protein('HGNC', 'ITGB3')]), {
         EVIDENCE: 'These were all explicitly stated in the BEL 2.0 Specification',
         CITATION: citation_2,
         RELATION: INCREASES,
         OBJECT: {
             MODIFIER: TRANSLOCATION,
             EFFECT: {
                 FROM_LOC: {NAMESPACE: 'GOCC', NAME: 'intracellular'},
                 TO_LOC: {NAMESPACE: 'GOCC', NAME: 'cell surface'}
             }
         },
     }),
    (complex_abundance([protein('HGNC', 'ITGAV'), protein('HGNC', 'ITGB3')]), protein('HGNC', 'ITGAV'), {
        RELATION: HAS_COMPONENT,
    }),
    (complex_abundance([protein('HGNC', 'ITGAV'), protein('HGNC', 'ITGB3')]), protein('HGNC', 'ITGB3'), {
        RELATION: HAS_COMPONENT,
    }),
    (gene('HGNC', 'ARRDC2'), gene('HGNC', 'ARRDC3'), {
        RELATION: EQUIVALENT_TO,
    }),
    (gene('HGNC', 'ARRDC3'), gene('HGNC', 'ARRDC2'), {
        RELATION: EQUIVALENT_TO,
        CITATION: citation_2,
        EVIDENCE: 'These were all explicitly stated in the BEL 2.0 Specification'
    }),
    (gene('dbSNP', 'rs123456'), gene('HGNC', 'CFTR', variants=hgvs('c.1521_1523delCTT')), {
        RELATION: ASSOCIATION,
        CITATION: citation_2,
        EVIDENCE: 'These were all explicitly stated in the BEL 2.0 Specification'
    }),
    (gene('HGNC', 'CFTR', variants=hgvs('c.1521_1523delCTT')), gene('dbSNP', 'rs123456'), {
        RELATION: ASSOCIATION,
        CITATION: citation_2,
        EVIDENCE: 'These were all explicitly stated in the BEL 2.0 Specification'
    }),
]


class BelReconstitutionMixin(TestGraphMixin):
    def bel_simple_reconstituted(self, graph, check_metadata=True):
        """Checks that test_bel.bel was loaded properly

        :param BELGraph graph: A BEL grpah
        :param bool check_metadata: Check the graph's document section is correct
        """
        self.assertIsNotNone(graph)
        self.assertIsInstance(graph, BELGraph)

        if check_metadata:
            self.assertIsNotNone(graph.document)
            self.assertEqual(expected_test_simple_metadata[METADATA_NAME], graph.name)
            self.assertEqual(expected_test_simple_metadata[METADATA_VERSION], graph.version)

        self.assertEqual(4, graph.number_of_nodes())

        # FIXME this should work, but is getting 8 for the upgrade function
        # self.assertEqual(6, graph.number_of_edges(), msg='Edges:\n{}'.format('\n'.join(map(str, graph.edges(keys=True, data=True)))))

        self.assertTrue(graph.has_node(protein(namespace='HGNC', name='AKT1')))
        self.assertTrue(graph.has_node(protein(namespace='HGNC', name='EGFR')))
        self.assertTrue(graph.has_node(protein(namespace='HGNC', name='FADD')))
        self.assertTrue(graph.has_node(protein(namespace='HGNC', name='CASP8')))

        bel_simple_citation_1 = {
            CITATION_NAME: "That one article from last week",
            CITATION_REFERENCE: "123455",
            CITATION_TYPE: "PubMed"
        }

        bel_simple_citation_2 = {
            CITATION_NAME: "That other article from last week",
            CITATION_REFERENCE: "123456",
            CITATION_TYPE: "PubMed"
        }

        evidence_1_extra = "Evidence 1 w extra notes"
        evidence_2 = 'Evidence 2'
        evidence_3 = 'Evidence 3'

        assertHasEdge(self, AKT1, EGFR, graph, **{
            RELATION: INCREASES,
            CITATION: bel_simple_citation_1,
            EVIDENCE: evidence_1_extra,
            ANNOTATIONS: {
                'Species': {'9606': True}
            }
        })
        assertHasEdge(self, EGFR, FADD, graph, **{
            RELATION: DECREASES,
            ANNOTATIONS: {
                'Species': {'9606': True},
                'CellLine': {'10B9 cell': True}
            },
            CITATION: bel_simple_citation_1,
            EVIDENCE: evidence_2
        })
        assertHasEdge(self, EGFR, CASP8, graph, **{
            RELATION: DIRECTLY_DECREASES,
            ANNOTATIONS: {
                'Species': {'9606': True},
                'CellLine': {'10B9 cell': True}
            },
            CITATION: bel_simple_citation_1,
            EVIDENCE: evidence_2,
        })
        assertHasEdge(self, FADD, CASP8, graph, **{
            RELATION: INCREASES,
            ANNOTATIONS: {
                'Species': {'10116': True}
            },
            CITATION: bel_simple_citation_2,
            EVIDENCE: evidence_3,
        })
        assertHasEdge(self, AKT1, CASP8, graph, **{
            RELATION: ASSOCIATION,
            ANNOTATIONS: {
                'Species': {'10116': True}
            },
            CITATION: bel_simple_citation_2,
            EVIDENCE: evidence_3,
        })
        assertHasEdge(self, CASP8, AKT1, graph, **{
            RELATION: ASSOCIATION,
            ANNOTATIONS: {
                'Species': {'10116': True}
            },
            CITATION: bel_simple_citation_2,
            EVIDENCE: evidence_3,
        })

    def bel_thorough_reconstituted(self, graph, check_metadata=True, check_warnings=True, check_provenance=True,
                                   check_citation_name=True):
        """Checks that thorough.bel was loaded properly

        :param BELGraph graph: A BEL graph
        :param bool check_metadata: Check the graph's document section is correct
        :param bool check_warnings: Check the graph produced the expected warnings
        :param bool check_provenance: Check the graph's definition section is correct
        :param bool check_citation_name: Check that the names in the citations get reconstituted. This isn't strictly
                                         necessary since this data can be looked up
        """
        self.assertIsNotNone(graph)
        self.assertIsInstance(graph, BELGraph)

        if check_warnings:
            self.assertEqual(0, len(graph.warnings),
                             msg='Document warnings:\n{}'.format('\n'.join(map(str, graph.warnings))))

        if check_metadata:
            self.assertLessEqual(set(expected_test_thorough_metadata), set(graph.document))
            self.assertEqual(expected_test_thorough_metadata[METADATA_NAME], graph.name)
            self.assertEqual(expected_test_thorough_metadata[METADATA_VERSION], graph.version)
            self.assertEqual(expected_test_thorough_metadata[METADATA_DESCRIPTION], graph.description)

        if check_provenance:
            self.assertEqual({'CHEBI', 'HGNC', 'GOBP', 'GOCC', 'MESHD', 'TESTNS2'}, set(graph.namespace_url))
            self.assertEqual(set(), set(graph.namespace_owl))
            self.assertEqual({'dbSNP'}, set(graph.namespace_pattern))
            self.assertEqual(set(), set(graph.annotation_owl))
            self.assertEqual({'TESTAN1', 'TESTAN2'}, set(graph.annotation_list))
            self.assertEqual({'TestRegex'}, set(graph.annotation_pattern))

        self.assertEqual(set(BEL_THOROUGH_NODES), set(graph))

        # FIXME
        # self.assertEqual(set((u, v) for u, v, _ in e), set(g.edges()))

        self.assertLess(0, graph.number_of_edges())

        for u, v, d in BEL_THOROUGH_EDGES:

            if not check_citation_name and CITATION in d and CITATION_NAME in d[CITATION]:
                d[CITATION] = d[CITATION].copy()
                del d[CITATION][CITATION_NAME]

            if d[RELATION] in unqualified_edges and CITATION in d:
                log.warning('shouldnt have citation with unqualified edge!')

                del d[CITATION]
                del d[EVIDENCE]

            assertHasEdge(self, u, v, graph, **d)

        for u, v, d in bel_thorough_graph.edges(data=True):
            assertHasEdge(self, u, v, graph, **d)

    def bel_slushy_reconstituted(self, graph, check_metadata=True, check_warnings=True):
        """Check that slushy.bel was loaded properly
        
        :param BELGraph graph: A BEL graph
        :param bool check_metadata: Check the graph's document section is correct
        :param bool check_warnings: Check the graph produced the expected warnings
        """
        self.assertIsNotNone(graph)
        self.assertIsInstance(graph, BELGraph)

        if check_metadata:
            self.assertIsNotNone(graph.document)
            self.assertIsInstance(graph.document, dict)
            self.assertEqual(expected_test_slushy_metadata[METADATA_NAME], graph.name)
            self.assertEqual(expected_test_slushy_metadata[METADATA_VERSION], graph.version)
            self.assertEqual(expected_test_slushy_metadata[METADATA_DESCRIPTION], graph.description)

        if check_warnings:
            expected_warnings = [
                (0, MissingMetadataException),
                (3, VersionFormatWarning),
                (26, MissingAnnotationKeyWarning),
                (29, MissingAnnotationKeyWarning),
                (34, InvalidCitationLengthException),
                (37, InvalidCitationType),
                (40, InvalidPubMedIdentifierWarning),
                (43, MissingCitationException),
                (48, MissingAnnotationKeyWarning),
                (51, MissingAnnotationKeyWarning),
                (54, MissingSupportWarning),
                (59, NakedNameWarning),
                (62, UndefinedNamespaceWarning),
                (65, MissingNamespaceNameWarning),
                (68, UndefinedAnnotationWarning),
                (71, MissingAnnotationKeyWarning),
                (74, IllegalAnnotationValueWarning),
                (77, MissingAnnotationRegexWarning),
                (80, MissingNamespaceRegexWarning),
                (83, MalformedTranslocationWarning),
                (86, PlaceholderAminoAcidWarning),
                (89, NestedRelationWarning),
                (92, InvalidFunctionSemantic),
                # (95, Exception),
                (98, BelSyntaxError),
            ]

            for (el, ew), (l, _, w, _) in zip(expected_warnings, graph.warnings):
                self.assertEqual(el, l, msg="Expected different error on line {}. Check line {}".format(el, l))
                self.assertIsInstance(w, ew, msg='Line: {}'.format(el))

        self.assertTrue(graph.has_node(protein(namespace='HGNC', name='AKT1')))
        self.assertTrue(graph.has_node(protein(namespace='HGNC', name='EGFR')))

        self.assertLess(0, graph.number_of_edges())

        assertHasEdge(self, AKT1, EGFR, graph, **{
            RELATION: INCREASES,
            CITATION: citation_1,
            EVIDENCE: evidence_1,
        })

    def bel_isolated_reconstituted(self, graph):
        """Runs the isolated node test

        :type graph: BELGraph
        """
        self.assertIsNotNone(graph)
        self.assertIsInstance(graph, BELGraph)

        self.assertHasNode(graph, adgrb1)
        self.assertHasNode(graph, adgrb2)
        self.assertHasNode(graph, adgrb_complex)
        self.assertHasNode(graph, achlorhydria)

        assertHasEdge(self, adgrb_complex, adgrb1, graph, relation=HAS_COMPONENT)
        assertHasEdge(self, adgrb_complex, adgrb2, graph, relation=HAS_COMPONENT)
