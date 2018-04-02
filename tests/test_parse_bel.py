# -*- coding: utf-8 -*-

import logging
import unittest

from pybel import BELGraph
from pybel.constants import *
from pybel.dsl import cell_surface_expression, entity, secretion, translocation
from pybel.dsl.nodes import (
    BaseAbundance, abundance, bioprocess, complex_abundance, composite_abundance, fragment, fusion_range, gene,
    gene_fusion, gmod, hgvs, hgvs_reference, hgvs_unspecified, mirna, named_complex_abundance, pathology, pmod, protein,
    protein_deletion, protein_fusion, protein_substitution, reaction, rna, rna_fusion,
)
from pybel.parser import BelParser
from pybel.parser.exc import MalformedTranslocationWarning
from pybel.parser.parse_bel import modifier_po_to_dict

from tests.constants import TestTokenParserBase, assertHasEdge, assertHasNode, update_provenance

log = logging.getLogger(__name__)

TEST_GENE_VARIANT = 'c.308G>A'
TEST_PROTEIN_VARIANT = 'p.Phe508del'


class TestAbundance(TestTokenParserBase):
    """2.1.1"""

    def setUp(self):
        self.parser.clear()
        self.parser.general_abundance.setParseAction(self.parser.handle_term)

        self.expected_node = abundance(namespace='CHEBI', name='oxygen atom')
        self.expected_node_tuple = self.expected_node.as_tuple()

        self.expected_canonical_bel = 'a(CHEBI:"oxygen atom")'

    def test_short_abundance(self):
        """small molecule"""
        statement = 'a(CHEBI:"oxygen atom")'

        self.parser.general_abundance.parseString(statement)

        nodes = list(self.graph)
        self.assertEqual(1, len(nodes))
        node = nodes[0]

        self.assertIsInstance(node, BaseAbundance)
        self.assertEqual(self.expected_node, node)
        self.assertEqual(self.expected_node_tuple, node.as_tuple())
        self.assertEqual(self.expected_canonical_bel, node.as_bel())

    def test_long_abundance(self):
        """small molecule"""
        statement = 'abundance(CHEBI:"oxygen atom", loc(GOCC:intracellular))'

        result = self.parser.general_abundance.parseString(statement)

        expected_result = {
            FUNCTION: ABUNDANCE,
            NAMESPACE: 'CHEBI',
            NAME: 'oxygen atom',
            LOCATION: {
                NAMESPACE: 'GOCC',
                NAME: 'intracellular'
            }
        }

        self.assertEqual(expected_result, result.asDict())

        modifier = modifier_po_to_dict(result)
        expected_modifier = {
            LOCATION: {NAMESPACE: 'GOCC', NAME: 'intracellular'}
        }
        self.assertEqual(expected_modifier, modifier)

        self.assertIn(self.expected_node, self.graph)


class TestGene(TestTokenParserBase):
    """2.1.4 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#XgeneA"""

    def setUp(self):
        self.parser.clear()
        self.parser.gene.setParseAction(self.parser.handle_term)

    def test_214a(self):
        statement = 'g(HGNC:AKT1)'

        result = self.parser.gene.parseString(statement)
        expected_list = [GENE, 'HGNC', 'AKT1']
        self.assertEqual(expected_list, result.asList())

        expected_dict = {
            FUNCTION: GENE,
            NAMESPACE: 'HGNC',
            NAME: 'AKT1'
        }
        self.assertEqual(expected_dict, result.asDict())

        node = gene('HGNC', 'AKT1')

        self.assertEqual(1, self.parser.graph.number_of_nodes())
        self.assertHasNode(node, **node)

    def test_214b(self):
        statement = 'g(HGNC:AKT1, loc(GOCC:intracellular))'

        result = self.parser.gene.parseString(statement)

        expected_dict = {
            FUNCTION: GENE,
            NAMESPACE: 'HGNC',
            NAME: 'AKT1',
            LOCATION: {
                NAMESPACE: 'GOCC',
                NAME: 'intracellular'
            }
        }

        self.assertEqual(expected_dict, result.asDict())

        expected_node = gene('HGNC', 'AKT1')
        self.assertHasNode(expected_node, **expected_node)

    def test_214c(self):
        """Test variant"""
        statement = 'g(HGNC:AKT1, var(p.Phe508del))'
        result = self.parser.gene.parseString(statement)

        expected_result = {
            FUNCTION: GENE,
            NAMESPACE: 'HGNC',
            NAME: 'AKT1',
            VARIANTS: [hgvs(TEST_PROTEIN_VARIANT)]
        }
        self.assertEqual(expected_result, result.asDict())

        expected_node = gene('HGNC', 'AKT1', variants=hgvs(TEST_PROTEIN_VARIANT))
        self.assertEqual('g(HGNC:AKT1, var(p.Phe508del))', expected_node.as_bel())

        parent = expected_node.get_parent()
        self.assertIsNotNone(parent)

        self.assertHasNode(expected_node, **expected_node)
        self.assertHasNode(parent, **parent)
        self.assertHasEdge(parent, expected_node, relation=HAS_VARIANT)

    def test_gmod(self):
        """Test Gene Modification"""
        statement = 'geneAbundance(HGNC:AKT1,gmod(M))'
        result = self.parser.gene.parseString(statement)

        expected_result = {
            FUNCTION: GENE,
            NAMESPACE: 'HGNC',
            NAME: 'AKT1',
            VARIANTS: [gmod('Me')]
        }
        self.assertEqual(expected_result, result.asDict())

        expected_node = gene('HGNC', 'AKT1', variants=gmod('Me'))
        self.assertEqual('g(HGNC:AKT1, gmod(Me))', expected_node.as_bel())

        self.assertHasNode(expected_node, **{FUNCTION: GENE})

        parent = expected_node.get_parent()
        self.assertIsNotNone(parent)
        self.assertHasNode(parent, **parent)
        self.assertHasEdge(parent, expected_node, relation=HAS_VARIANT)

    def test_214d(self):
        """Test BEL 1.0 gene substitution"""
        statement = 'g(HGNC:AKT1,sub(G,308,A))'
        result = self.parser.gene.parseString(statement)

        expected_result = gene(
            name='AKT1',
            namespace='HGNC',
            variants=[hgvs(TEST_GENE_VARIANT)]
        )

        self.assertEqual(expected_result, result.asDict())
        self.assertEqual('g(HGNC:AKT1, var(c.308G>A))', expected_result.as_bel())

        parent = expected_result.get_parent()
        self.assertIsNotNone(parent)
        self.assertHasNode(parent, **parent)

        self.assertHasEdge(parent, expected_result, relation=HAS_VARIANT)

    def test_variant_location(self):
        """Test BEL 1.0 gene substitution with location tag"""
        statement = 'g(HGNC:AKT1,sub(G,308,A),loc(GOCC:intracellular))'
        result = self.parser.gene.parseString(statement)

        expected_result = {
            FUNCTION: GENE,
            NAMESPACE: 'HGNC',
            NAME: 'AKT1',
            VARIANTS: [
                {
                    KIND: HGVS,
                    IDENTIFIER: TEST_GENE_VARIANT
                }
            ],
            LOCATION: {
                NAMESPACE: 'GOCC',
                NAME: 'intracellular'
            }
        }
        self.assertEqual(expected_result, result.asDict())

        expected_node = gene('HGNC', 'AKT1', variants=hgvs(TEST_GENE_VARIANT))
        self.assertEqual('g(HGNC:AKT1, var(c.308G>A))', expected_node.as_bel())
        self.assertHasNode(expected_node, **expected_node)

        parent = expected_node.get_parent()
        self.assertIsNotNone(parent)
        self.assertHasNode(parent, **parent)

        self.assertHasEdge(parent, expected_node, relation=HAS_VARIANT)

    def test_multiple_variants(self):
        """Test multiple variants"""
        statement = 'g(HGNC:AKT1, var(p.Phe508del), sub(G,308,A), var(c.1521_1523delCTT))'
        result = self.parser.gene.parseString(statement)

        expected_result = {
            FUNCTION: GENE,
            NAMESPACE: 'HGNC',
            NAME: 'AKT1',
            VARIANTS: [
                hgvs(TEST_PROTEIN_VARIANT),
                hgvs(TEST_GENE_VARIANT),
                hgvs('c.1521_1523delCTT')
            ]
        }
        self.assertEqual(expected_result, result.asDict())

        node = gene('HGNC', 'AKT1', variants=[
            hgvs('c.1521_1523delCTT'),
            hgvs(TEST_GENE_VARIANT),
            hgvs(TEST_PROTEIN_VARIANT)
        ])
        self.assertEqual('g(HGNC:AKT1, var(c.1521_1523delCTT), var(c.308G>A), var(p.Phe508del))', node.as_bel())
        self.assertHasNode(node, **node)
        self.help_test_parent_in_graph(node)

    def test_gene_fusion_1(self):
        self.maxDiff = None
        statement = 'g(fus(HGNC:TMPRSS2, c.1_79, HGNC:ERG, c.312_5034))'
        result = self.parser.gene.parseString(statement)
        expected_dict = {
            FUNCTION: GENE,
            FUSION: {
                PARTNER_5P: {NAMESPACE: 'HGNC', NAME: 'TMPRSS2'},
                PARTNER_3P: {NAMESPACE: 'HGNC', NAME: 'ERG'},
                RANGE_5P: {
                    FUSION_REFERENCE: 'c',
                    FUSION_START: 1,
                    FUSION_STOP: 79

                },
                RANGE_3P: {
                    FUSION_REFERENCE: 'c',
                    FUSION_START: 312,
                    FUSION_STOP: 5034
                }
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = gene_fusion(
            partner_5p=gene('HGNC', 'TMPRSS2'),
            range_5p=fusion_range('c', 1, 79),
            partner_3p=gene('HGNC', 'ERG'),
            range_3p=fusion_range('c', 312, 5034)
        )
        self.assertEqual(statement, expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_gene_fusion_2(self):
        self.maxDiff = None
        statement = 'g(fus(HGNC:TMPRSS2, c.1_?, HGNC:ERG, c.312_5034))'
        result = self.parser.gene.parseString(statement)
        expected_dict = {
            FUNCTION: GENE,
            FUSION: {
                PARTNER_5P: {NAMESPACE: 'HGNC', NAME: 'TMPRSS2'},
                PARTNER_3P: {NAMESPACE: 'HGNC', NAME: 'ERG'},
                RANGE_5P: {
                    FUSION_REFERENCE: 'c',
                    FUSION_START: 1,
                    FUSION_STOP: '?'

                },
                RANGE_3P: {
                    FUSION_REFERENCE: 'c',
                    FUSION_START: 312,
                    FUSION_STOP: 5034
                }
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = gene_fusion(
            partner_5p=gene('HGNC', 'TMPRSS2'),
            range_5p=fusion_range('c', 1, '?'),
            partner_3p=gene('HGNC', 'ERG'),
            range_3p=fusion_range('c', 312, 5034)
        )
        self.assertEqual(statement, expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_gene_fusion_3(self):
        self.maxDiff = None
        statement = 'g(fus(HGNC:TMPRSS2, ?, HGNC:ERG, c.312_5034))'
        result = self.parser.gene.parseString(statement)
        expected_dict = {
            FUNCTION: GENE,
            FUSION: {
                PARTNER_5P: {NAMESPACE: 'HGNC', NAME: 'TMPRSS2'},
                PARTNER_3P: {NAMESPACE: 'HGNC', NAME: 'ERG'},
                RANGE_5P: {
                    FUSION_MISSING: '?'
                },
                RANGE_3P: {
                    FUSION_REFERENCE: 'c',
                    FUSION_START: 312,
                    FUSION_STOP: 5034
                }
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = gene_fusion(
            partner_5p=gene('HGNC', 'TMPRSS2'),
            partner_3p=gene('HGNC', 'ERG'),
            range_3p=fusion_range('c', 312, 5034)
        )
        self.assertEqual(statement, expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_gene_fusion_legacy_1(self):
        self.maxDiff = None
        statement = 'g(HGNC:BCR, fus(HGNC:JAK2, 1875, 2626))'
        result = self.parser.gene.parseString(statement)

        expected_dict = {
            FUNCTION: GENE,
            FUSION: {
                PARTNER_5P: {NAMESPACE: 'HGNC', NAME: 'BCR'},
                PARTNER_3P: {NAMESPACE: 'HGNC', NAME: 'JAK2'},
                RANGE_5P: {
                    FUSION_REFERENCE: 'c',
                    FUSION_START: '?',
                    FUSION_STOP: 1875

                },
                RANGE_3P: {
                    FUSION_REFERENCE: 'c',
                    FUSION_START: 2626,
                    FUSION_STOP: '?'
                }
            }
        }

        self.assertEqual(expected_dict, result.asDict())

        expected_node = gene_fusion(
            partner_5p=gene('HGNC', 'BCR'),
            range_5p=fusion_range('c', '?', 1875),
            partner_3p=gene('HGNC', 'JAK2'),
            range_3p=fusion_range('c', 2626, '?')
        )
        self.assertEqual('g(fus(HGNC:BCR, c.?_1875, HGNC:JAK2, c.2626_?))', expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_gene_fusion_legacy_2(self):
        statement = 'g(HGNC:CHCHD4, fusion(HGNC:AIFM1))'
        result = self.parser.gene.parseString(statement)

        expected_dict = {
            FUNCTION: GENE,
            FUSION: {
                PARTNER_5P: {NAMESPACE: 'HGNC', NAME: 'CHCHD4'},
                PARTNER_3P: {NAMESPACE: 'HGNC', NAME: 'AIFM1'},
                RANGE_5P: {FUSION_MISSING: '?'},
                RANGE_3P: {FUSION_MISSING: '?'}
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = gene_fusion(
            partner_5p=gene('HGNC', 'CHCHD4'),
            partner_3p=gene('HGNC', 'AIFM1'),
        )
        self.assertEqual('g(fus(HGNC:CHCHD4, ?, HGNC:AIFM1, ?))', expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_gene_variant_snp(self):
        """2.2.2 SNP"""
        statement = 'g(SNP:rs113993960, var(c.1521_1523delCTT))'
        result = self.parser.gene.parseString(statement)

        expected_result = [GENE, 'SNP', 'rs113993960', [HGVS, 'c.1521_1523delCTT']]
        self.assertEqual(expected_result, result.asList())

        expected_node = gene('SNP', 'rs113993960', variants=hgvs('c.1521_1523delCTT'))
        self.assertEqual('g(SNP:rs113993960, var(c.1521_1523delCTT))', expected_node.as_bel())
        self.assertHasNode(expected_node)

        parent = expected_node.get_parent()
        self.assertIsNotNone(parent)
        self.assertHasNode(parent)

        self.assertHasEdge(parent, expected_node, relation=HAS_VARIANT)

    def test_gene_variant_chromosome(self):
        """2.2.2 chromosome"""
        statement = 'g(REF:"NC_000007.13", var(g.117199646_117199648delCTT))'
        result = self.parser.gene.parseString(statement)

        expected_result = [GENE, 'REF', 'NC_000007.13', [HGVS, 'g.117199646_117199648delCTT']]
        self.assertEqual(expected_result, result.asList())

        gene_node = gene('REF', 'NC_000007.13', variants=hgvs('g.117199646_117199648delCTT'))
        self.assertEqual(statement, gene_node.as_bel())
        self.assertHasNode(gene_node, **gene_node)

        parent = gene_node.get_parent()
        self.assertIsNotNone(parent)
        self.assertHasNode(parent, **parent)

        self.assertHasEdge(parent, gene_node, relation=HAS_VARIANT)

    def test_gene_variant_deletion(self):
        """2.2.2 gene-coding DNA reference sequence"""
        statement = 'g(HGNC:CFTR, var(c.1521_1523delCTT))'
        result = self.parser.gene.parseString(statement)

        expected_result = {
            FUNCTION: GENE,
            NAMESPACE: 'HGNC',
            NAME: 'CFTR',
            VARIANTS: [
                {KIND: HGVS, IDENTIFIER: 'c.1521_1523delCTT'}
            ]
        }
        self.assertEqual(expected_result, result.asDict())

        expected_node = gene('HGNC', 'CFTR', variants=hgvs('c.1521_1523delCTT'))
        self.assertEqual(statement, expected_node.as_bel())
        self.assertHasNode(expected_node, **expected_node)

        parent = expected_node.get_parent()
        self.assertIsNotNone(parent)
        self.assertHasNode(parent, **parent)

        self.assertHasEdge(parent, expected_node, relation=HAS_VARIANT)


class TestMiRNA(TestTokenParserBase):
    """2.1.5 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#XmicroRNAA"""

    def setUp(self):
        self.parser.clear()
        self.parser.mirna.setParseAction(self.parser.handle_term)

    def test_short(self):
        statement = 'm(HGNC:MIR21)'
        result = self.parser.mirna.parseString(statement)
        expected_result = [MIRNA, 'HGNC', 'MIR21']
        self.assertEqual(expected_result, result.asList())

        expected_dict = {
            FUNCTION: MIRNA,
            NAMESPACE: 'HGNC',
            NAME: 'MIR21'
        }

        self.assertEqual(expected_dict, result.asDict())

        node = mirna('HGNC', 'MIR21')
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node, **node)

    def test_long(self):
        statement = 'microRNAAbundance(HGNC:MIR21)'
        result = self.parser.mirna.parseString(statement)
        expected_result = [MIRNA, 'HGNC', 'MIR21']
        self.assertEqual(expected_result, result.asList())

        expected_dict = {
            FUNCTION: MIRNA,
            NAMESPACE: 'HGNC',
            NAME: 'MIR21',
        }

        self.assertEqual(expected_dict, result.asDict())

        expected_node = mirna('HGNC', 'MIR21')
        self.assertEqual('m(HGNC:MIR21)', expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_mirna_location(self):
        statement = 'm(HGNC:MIR21,loc(GOCC:intracellular))'
        result = self.parser.mirna.parseString(statement)

        expected_dict = {
            FUNCTION: MIRNA,
            NAMESPACE: 'HGNC',
            NAME: 'MIR21',
            LOCATION: {
                NAMESPACE: 'GOCC',
                NAME: 'intracellular'
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = mirna('HGNC', 'MIR21')
        self.assertEqual('m(HGNC:MIR21)', expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_mirna_variant(self):
        statement = 'm(HGNC:MIR21,var(p.Phe508del))'
        result = self.parser.mirna.parseString(statement)

        expected_dict = {
            FUNCTION: MIRNA,
            NAMESPACE: 'HGNC',
            NAME: 'MIR21',
            VARIANTS: [
                {
                    KIND: HGVS,
                    IDENTIFIER: TEST_PROTEIN_VARIANT
                },
            ]
        }
        self.assertEqual(expected_dict, result.asDict())

        self.assertEqual(2, self.graph.number_of_nodes())

        expected_node = mirna('HGNC', 'MIR21', variants=hgvs(TEST_PROTEIN_VARIANT))
        self.assertEqual('m(HGNC:MIR21, var(p.Phe508del))', expected_node.as_bel())
        self.assertHasNode(expected_node)

        parent = expected_node.get_parent()
        self.assertIsNotNone(parent)
        self.assertHasNode(parent)

        self.assertHasEdge(parent, expected_node, relation=HAS_VARIANT)

    def test_mirna_variant_location(self):
        statement = 'm(HGNC:MIR21,var(p.Phe508del),loc(GOCC:intracellular))'
        result = self.parser.mirna.parseString(statement)

        expected_dict = {
            FUNCTION: MIRNA,
            NAMESPACE: 'HGNC',
            NAME: 'MIR21',
            VARIANTS: [
                {
                    KIND: HGVS,
                    IDENTIFIER: TEST_PROTEIN_VARIANT
                },
            ],
            LOCATION: {
                NAMESPACE: 'GOCC',
                NAME: 'intracellular'
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        self.assertEqual(2, self.graph.number_of_nodes())

        expected_node = mirna('HGNC', 'MIR21', variants=hgvs(TEST_PROTEIN_VARIANT))
        self.assertEqual('m(HGNC:MIR21, var(p.Phe508del))', expected_node.as_bel())
        self.assertHasNode(expected_node)

        parent = expected_node.get_parent()
        self.assertIsNotNone(parent)
        self.assertHasNode(parent)

        self.assertHasEdge(parent, expected_node, relation=HAS_VARIANT)


class TestProtein(TestTokenParserBase):
    """2.1.6 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#XproteinA"""

    def setUp(self):
        self.parser.clear()

        # Just activate and streamline the protein parser
        self.parser.protein.setParseAction(self.parser.handle_term)

    def test_216a(self):
        statement = 'p(HGNC:AKT1)'

        result = self.parser.protein.parseString(statement)
        expected_result = [PROTEIN, 'HGNC', 'AKT1']
        self.assertEqual(expected_result, result.asList())

        expected_dict = {
            FUNCTION: PROTEIN,
            NAMESPACE: 'HGNC',
            NAME: 'AKT1',
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = protein('HGNC', 'AKT1')
        self.assertEqual(statement, expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_protein_wLocation(self):
        statement = 'p(HGNC:AKT1, loc(GOCC:intracellular))'

        result = self.parser.protein.parseString(statement)

        expected_dict = {
            FUNCTION: PROTEIN,
            NAMESPACE: 'HGNC',
            NAME: 'AKT1',
            LOCATION: {
                NAMESPACE: 'GOCC',
                NAME: 'intracellular'
            }
        }

        self.assertEqual(expected_dict, result.asDict())

        expected_node = protein('HGNC', 'AKT1')
        self.assertEqual('p(HGNC:AKT1)', expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_multiVariant(self):
        statement = 'p(HGNC:AKT1,sub(A,127,Y),pmod(Ph, Ser),loc(GOCC:intracellular))'

        result = self.parser.protein.parseString(statement)

        expected_dict = {
            FUNCTION: PROTEIN,
            NAMESPACE: 'HGNC',
            NAME: 'AKT1',
            LOCATION: {
                NAMESPACE: 'GOCC',
                NAME: 'intracellular'
            },
            VARIANTS: [
                hgvs('p.Ala127Tyr'),
                pmod(name='Ph', code='Ser')
            ]
        }

        self.assertEqual(expected_dict, result.asDict())

        node = protein('HGNC', 'AKT1', variants=[protein_substitution('Ala', 127, 'Tyr'), pmod(name='Ph', code='Ser')])
        self.assertEqual('p(HGNC:AKT1, pmod(Ph, Ser), var(p.Ala127Tyr))', node.as_bel())
        self.assertHasNode(node, **node)

        parent = node.get_parent()
        self.assertIsNotNone(parent)
        self.assertHasNode(parent, **parent)
        self.assertHasEdge(parent, node, relation=HAS_VARIANT)

    def test_protein_fusion_1(self):
        statement = 'p(fus(HGNC:TMPRSS2, p.1_79, HGNC:ERG, p.312_5034))'
        result = self.parser.protein.parseString(statement)
        expected_dict = {
            FUNCTION: PROTEIN,
            FUSION: {
                PARTNER_5P: {NAMESPACE: 'HGNC', NAME: 'TMPRSS2'},
                PARTNER_3P: {NAMESPACE: 'HGNC', NAME: 'ERG'},
                RANGE_5P: {
                    FUSION_REFERENCE: 'p',
                    FUSION_START: 1,
                    FUSION_STOP: 79

                },
                RANGE_3P: {
                    FUSION_REFERENCE: 'p',
                    FUSION_START: 312,
                    FUSION_STOP: 5034

                }
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = protein_fusion(
            partner_5p=protein('HGNC', 'TMPRSS2'),
            partner_3p=protein('HGNC', 'ERG'),
            range_5p=fusion_range('p', 1, 79),
            range_3p=fusion_range('p', 312, 5034)
        )
        self.assertEqual(statement, expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_protein_fusion_legacy_1(self):
        """Tests parsing of a BEL v1.0 protein fusion that contains a 5P stop adn 3P start"""
        self.maxDiff = None
        statement = 'p(HGNC:BCR, fus(HGNC:JAK2, 1875, 2626))'
        result = self.parser.protein.parseString(statement)

        expected_dict = {
            FUNCTION: PROTEIN,
            FUSION: {
                PARTNER_5P: {NAMESPACE: 'HGNC', NAME: 'BCR'},
                PARTNER_3P: {NAMESPACE: 'HGNC', NAME: 'JAK2'},
                RANGE_5P: {
                    FUSION_REFERENCE: 'p',
                    FUSION_START: '?',
                    FUSION_STOP: 1875

                },
                RANGE_3P: {
                    FUSION_REFERENCE: 'p',
                    FUSION_START: 2626,
                    FUSION_STOP: '?'

                }
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = protein_fusion(
            partner_5p=protein('HGNC', 'BCR'),
            partner_3p=protein('HGNC', 'JAK2'),
            range_5p=fusion_range('p', '?', 1875),
            range_3p=fusion_range('p', 2626, '?')
        )
        self.assertEqual('p(fus(HGNC:BCR, p.?_1875, HGNC:JAK2, p.2626_?))', expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_protein_fusion_legacy_2(self):
        statement = 'p(HGNC:CHCHD4, fusion(HGNC:AIFM1))'
        result = self.parser.protein.parseString(statement)

        expected_dict = {
            FUNCTION: PROTEIN,
            FUSION: {
                PARTNER_5P: {NAMESPACE: 'HGNC', NAME: 'CHCHD4'},
                PARTNER_3P: {NAMESPACE: 'HGNC', NAME: 'AIFM1'},
                RANGE_5P: {FUSION_MISSING: '?'},
                RANGE_3P: {FUSION_MISSING: '?'}
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = protein_fusion(partner_5p=protein('HGNC', 'CHCHD4'), partner_3p=protein('HGNC', 'AIFM1'))
        self.assertEqual('p(fus(HGNC:CHCHD4, ?, HGNC:AIFM1, ?))', expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_protein_trunc_1(self):
        statement = 'p(HGNC:AKT1, trunc(40))'
        result = self.parser.protein.parseString(statement)

        node = protein('HGNC', 'AKT1', variants=hgvs('p.40*'))  # FIXME replace with truncation?
        self.assertEqual('p(HGNC:AKT1, var(p.40*))', node.as_bel())
        self.assertHasNode(node, **{FUNCTION: PROTEIN})
        self.help_test_parent_in_graph(node)

    def test_protein_trunc_2(self):
        statement = 'p(HGNC:AKT1, var(p.Cys40*))'
        result = self.parser.protein.parseString(statement)

        expected_result = [PROTEIN, 'HGNC', 'AKT1', [HGVS, 'p.Cys40*']]
        self.assertEqual(expected_result, result.asList())

        node = protein('HGNC', 'AKT1', variants=hgvs('p.Cys40*'))
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node, **{FUNCTION: PROTEIN})
        self.help_test_parent_in_graph(node)

    def test_protein_trunc_3(self):
        statement = 'p(HGNC:AKT1, var(p.Arg1851*))'
        result = self.parser.protein.parseString(statement)

        expected_result = [PROTEIN, 'HGNC', 'AKT1', [HGVS, 'p.Arg1851*']]
        self.assertEqual(expected_result, result.asList())

        node = protein('HGNC', 'AKT1', variants=hgvs('p.Arg1851*'))
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node, **{FUNCTION: PROTEIN})
        self.help_test_parent_in_graph(node)

    def test_protein_pmod_1(self):
        """2.2.1 Test default BEL namespace and 1-letter amino acid code:"""
        statement = 'p(HGNC:AKT1, pmod(Ph, S, 473))'
        self.parser.protein.parseString(statement)

        node = protein('HGNC', 'AKT1', variants=pmod('Ph', code='Ser', position=473))
        self.assertEqual('p(HGNC:AKT1, pmod(Ph, Ser, 473))', node.as_bel())
        self.assertHasNode(node, **node)
        self.help_test_parent_in_graph(node)

    def test_protein_pmod_2(self):
        """2.2.1 Test default BEL namespace and 3-letter amino acid code:"""
        statement = 'p(HGNC:AKT1, pmod(Ph, Ser, 473))'
        result = self.parser.protein.parseString(statement)

        node = protein('HGNC', 'AKT1', variants=pmod('Ph', code='Ser', position=473))
        self.assertEqual('p(HGNC:AKT1, pmod(Ph, Ser, 473))', node.as_bel())
        self.assertHasNode(node, **node)
        self.help_test_parent_in_graph(node)

    def test_protein_pmod_3(self):
        """2.2.1 Test PSI-MOD namespace and 3-letter amino acid code:"""
        statement = 'p(HGNC:AKT1, pmod(MOD:PhosRes, Ser, 473))'
        result = self.parser.protein.parseString(statement)

        node = protein('HGNC', 'AKT1', variants=pmod(namespace='MOD', name='PhosRes', code='Ser', position=473))
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node)
        self.help_test_parent_in_graph(node)

    def test_protein_pmod_4(self):
        """2.2.1 Test HRAS palmitoylated at an unspecified residue. Default BEL namespace"""
        statement = 'p(HGNC:HRAS, pmod(Palm))'
        result = self.parser.protein.parseString(statement)

        node = protein('HGNC', 'HRAS', variants=pmod('Palm'))
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node, **node)
        self.help_test_parent_in_graph(node)

    def test_protein_variant_reference(self):
        """2.2.2 Test reference allele"""
        statement = 'p(HGNC:CFTR, var(=))'
        result = self.parser.protein.parseString(statement)
        expected_result = [PROTEIN, 'HGNC', 'CFTR', [HGVS, '=']]
        self.assertEqual(expected_result, result.asList())

        node = protein('HGNC', 'CFTR', variants=hgvs_reference())
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node, **node)
        self.help_test_parent_in_graph(node)

    def test_protein_variant_unspecified(self):
        """2.2.2 Test unspecified variant"""
        statement = 'p(HGNC:CFTR, var(?))'
        result = self.parser.protein.parseString(statement)

        expected_result = [PROTEIN, 'HGNC', 'CFTR', [HGVS, '?']]
        self.assertEqual(expected_result, result.asList())

        node = protein('HGNC', 'CFTR', variants=hgvs_unspecified())
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node, **node)
        self.help_test_parent_in_graph(node)

    def test_protein_variant_substitution(self):
        """2.2.2 Test substitution"""
        statement = 'p(HGNC:CFTR, var(p.Gly576Ala))'
        result = self.parser.protein.parseString(statement)
        expected_result = [PROTEIN, 'HGNC', 'CFTR', [HGVS, 'p.Gly576Ala']]
        self.assertEqual(expected_result, result.asList())

        node = protein('HGNC', 'CFTR', variants=protein_substitution('Gly', 576, 'Ala'))
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node, **node)
        self.help_test_parent_in_graph(node)

    def test_protein_variant_deletion(self):
        """2.2.2 deletion"""
        statement = 'p(HGNC:CFTR, var(p.Phe508del))'
        result = self.parser.protein.parseString(statement)

        expected_result = [PROTEIN, 'HGNC', 'CFTR', [HGVS, TEST_PROTEIN_VARIANT]]
        self.assertEqual(expected_result, result.asList())

        node = protein('HGNC', 'CFTR', variants=protein_deletion('Phe', 508))
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node, **node)
        self.help_test_parent_in_graph(node)

    def test_protein_fragment_known(self):
        """2.2.3 fragment with known start/stop"""
        statement = 'p(HGNC:YFG, frag(5_20))'
        result = self.parser.protein.parseString(statement)

        node = protein('HGNC', 'YFG', variants=fragment(5, 20))
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node, **node)
        self.help_test_parent_in_graph(node)

    def test_protein_fragment_unbounded(self):
        """2.2.3 amino-terminal fragment of unknown length"""
        statement = 'p(HGNC:YFG, frag(1_?))'
        result = self.parser.protein.parseString(statement)

        node = protein('HGNC', 'YFG', variants=fragment(start=1))
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node, **node)
        self.help_test_parent_in_graph(node)

    def test_protein_fragment_unboundTerminal(self):
        """2.2.3 carboxyl-terminal fragment of unknown length"""
        statement = 'p(HGNC:YFG, frag(?_*))'
        result = self.parser.protein.parseString(statement)

        node = protein('HGNC', 'YFG', variants=fragment(stop='*'))
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node, **node)
        self.help_test_parent_in_graph(node)

    def test_protein_fragment_unknown(self):
        """2.2.3 fragment with unknown start/stop"""
        statement = 'p(HGNC:YFG, frag(?))'
        result = self.parser.protein.parseString(statement)

        node = protein('HGNC', 'YFG', variants=fragment())
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node, **node)
        self.help_test_parent_in_graph(node)

    def test_protein_fragment_descriptor(self):
        """2.2.3 fragment with unknown start/stop and a descriptor"""
        statement = 'p(HGNC:YFG, frag(?, "55kD"))'
        result = self.parser.protein.parseString(statement)

        node = protein('HGNC', 'YFG', variants=fragment(description='55kD'))
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node, **node)
        self.help_test_parent_in_graph(node)

    def test_ensure_no_dup_edges(self):
        """Ensure node and edges aren't added twice, even if from different statements and has origin completion"""
        s1 = 'p(HGNC:AKT1)'
        s2 = 'deg(p(HGNC:AKT1))'

        self.parser.bel_term.parseString(s1)
        self.parser.bel_term.parseString(s2)

        node = protein('HGNC', 'AKT1')

        self.assertEqual(1, self.parser.graph.number_of_nodes())
        self.assertEqual(0, self.parser.graph.number_of_edges())
        self.assertHasNode(node)


class TestRna(TestTokenParserBase):
    """2.1.7 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#XrnaA"""

    def setUp(self):
        self.parser.clear()
        self.parser.rna.setParseAction(self.parser.handle_term)

    def test_217a(self):
        statement = 'r(HGNC:AKT1)'

        result = self.parser.rna.parseString(statement)
        expected_result = [RNA, 'HGNC', 'AKT1']
        self.assertEqual(expected_result, result.asList())

        expected_dict = {
            FUNCTION: RNA,
            NAMESPACE: 'HGNC',
            NAME: 'AKT1'
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = rna('HGNC', 'AKT1')
        self.assertEqual(statement, expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_214e(self):
        """Test multiple variants"""
        statement = 'r(HGNC:AKT1, var(p.Phe508del), var(c.1521_1523delCTT))'
        result = self.parser.rna.parseString(statement)

        expected_result = {
            FUNCTION: RNA,
            NAMESPACE: 'HGNC',
            NAME: 'AKT1',
            VARIANTS: [
                hgvs(TEST_PROTEIN_VARIANT),
                hgvs('c.1521_1523delCTT')
            ]
        }
        self.assertEqual(expected_result, result.asDict())

        node = rna('HGNC', 'AKT1', variants=[hgvs('c.1521_1523delCTT'), hgvs(TEST_PROTEIN_VARIANT)])
        self.assertEqual('r(HGNC:AKT1, var(c.1521_1523delCTT), var(p.Phe508del))', node.as_bel())
        self.assertHasNode(node, **node)
        self.help_test_parent_in_graph(node)

    def test_rna_fusion_1(self):
        """2.6.1 RNA abundance of fusion with known breakpoints"""
        statement = 'r(fus(HGNC:TMPRSS2, r.1_79, HGNC:ERG, r.312_5034))'
        result = self.parser.rna.parseString(statement)

        expected_dict = {
            FUNCTION: RNA,
            FUSION: {
                PARTNER_5P: {NAMESPACE: 'HGNC', NAME: 'TMPRSS2'},
                PARTNER_3P: {NAMESPACE: 'HGNC', NAME: 'ERG'},
                RANGE_5P: {
                    FUSION_REFERENCE: 'r',
                    FUSION_START: 1,
                    FUSION_STOP: 79

                },
                RANGE_3P: {
                    FUSION_REFERENCE: 'r',
                    FUSION_START: 312,
                    FUSION_STOP: 5034
                }
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = rna_fusion(
            partner_5p=rna('HGNC', 'TMPRSS2'),
            partner_3p=rna('HGNC', 'ERG'),
            range_5p=fusion_range('r', 1, 79),
            range_3p=fusion_range('r', 312, 5034)
        )
        self.assertEqual(statement, expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_rna_fusion_2(self):
        """2.6.1 RNA abundance of fusion with unspecified breakpoints"""
        statement = 'r(fus(HGNC:TMPRSS2, ?, HGNC:ERG, ?))'
        result = self.parser.rna.parseString(statement)

        expected_dict = {
            FUNCTION: RNA,
            FUSION: {
                PARTNER_5P: {NAMESPACE: 'HGNC', NAME: 'TMPRSS2'},
                PARTNER_3P: {NAMESPACE: 'HGNC', NAME: 'ERG'},
                RANGE_5P: {
                    FUSION_MISSING: '?',
                },
                RANGE_3P: {
                    FUSION_MISSING: '?',
                }
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = rna_fusion(
            partner_5p=rna('HGNC', 'TMPRSS2'),
            partner_3p=rna('HGNC', 'ERG'),
        )
        self.assertEqual(statement, expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_rna_fusion_legacy_1(self):
        statement = 'r(HGNC:BCR, fus(HGNC:JAK2, 1875, 2626))'
        result = self.parser.rna.parseString(statement)

        expected_dict = {
            FUNCTION: RNA,
            FUSION: {
                PARTNER_5P: {NAMESPACE: 'HGNC', NAME: 'BCR'},
                PARTNER_3P: {NAMESPACE: 'HGNC', NAME: 'JAK2'},
                RANGE_5P: {
                    FUSION_REFERENCE: 'r',
                    FUSION_START: '?',
                    FUSION_STOP: 1875

                },
                RANGE_3P: {
                    FUSION_REFERENCE: 'r',
                    FUSION_START: 2626,
                    FUSION_STOP: '?'
                }
            }

        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = rna_fusion(
            partner_5p=rna('HGNC', 'BCR'),
            partner_3p=rna('HGNC', 'JAK2'),
            range_5p=fusion_range('r', '?', 1875),
            range_3p=fusion_range('r', 2626, '?')
        )
        self.assertEqual('r(fus(HGNC:BCR, r.?_1875, HGNC:JAK2, r.2626_?))', expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_rna_fusion_legacy_2(self):
        statement = 'r(HGNC:CHCHD4, fusion(HGNC:AIFM1))'
        result = self.parser.rna.parseString(statement)

        expected_dict = {
            FUNCTION: RNA,
            FUSION: {
                PARTNER_5P: {NAMESPACE: 'HGNC', NAME: 'CHCHD4'},
                PARTNER_3P: {NAMESPACE: 'HGNC', NAME: 'AIFM1'},
                RANGE_5P: {FUSION_MISSING: '?'},
                RANGE_3P: {FUSION_MISSING: '?'}
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = rna_fusion(
            partner_5p=rna('HGNC', 'CHCHD4'),
            partner_3p=rna('HGNC', 'AIFM1'),
        )
        self.assertEqual('r(fus(HGNC:CHCHD4, ?, HGNC:AIFM1, ?))', expected_node.as_bel())
        self.assertHasNode(expected_node)

    def test_rna_variant_codingReference(self):
        """2.2.2 RNA coding reference sequence"""
        statement = 'r(HGNC:CFTR, var(r.1521_1523delcuu))'
        result = self.parser.rna.parseString(statement)

        expected_dict = {
            FUNCTION: RNA,
            NAMESPACE: 'HGNC',
            NAME: 'CFTR',
            VARIANTS: [hgvs('r.1521_1523delcuu')]
        }
        self.assertEqual(expected_dict, result.asDict())

        node = rna('HGNC', 'CFTR', variants=hgvs('r.1521_1523delcuu'))
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node)
        self.help_test_parent_in_graph(node)


class TestComplex(TestTokenParserBase):
    """2.1.2 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#XcomplexA"""

    def setUp(self):
        self.parser.clear()
        self.parser.complex_abundances.setParseAction(self.parser.handle_term)

    def test_complex_singleton(self):
        statement = 'complex(SCOMP:"AP-1 Complex")'
        result = self.parser.complex_abundances.parseString(statement)

        expected_result = [COMPLEX, 'SCOMP', 'AP-1 Complex']
        self.assertEqual(expected_result, result.asList())

        expected_dict = {
            FUNCTION: COMPLEX,
            NAMESPACE: 'SCOMP',
            NAME: 'AP-1 Complex'
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = named_complex_abundance('SCOMP', 'AP-1 Complex')
        self.assertEqual(statement, expected_node.as_bel())
        self.assertHasNode(expected_node, **expected_node)

    def test_complex_list_short(self):
        statement = 'complex(p(HGNC:FOS), p(HGNC:JUN))'
        result = self.parser.complex_abundances.parseString(statement)

        expected_result = [COMPLEX, [PROTEIN, 'HGNC', 'FOS'], [PROTEIN, 'HGNC', 'JUN']]
        self.assertEqual(expected_result, result.asList())

        expected_result = {
            FUNCTION: COMPLEX,
            MEMBERS: [
                {
                    FUNCTION: PROTEIN,
                    NAMESPACE: 'HGNC',
                    NAME: 'FOS'
                }, {
                    FUNCTION: PROTEIN,
                    NAMESPACE: 'HGNC',
                    NAME: 'JUN'
                }
            ]
        }
        self.assertEqual(expected_result, result.asDict())

        fos = protein('HGNC', 'FOS')
        jun = protein('HGNC', 'JUN')
        node = complex_abundance([fos, jun])
        self.assertEqual(statement, node.as_bel())
        self.assertHasNode(node, **node)

        self.assertHasNode(fos, **fos)
        self.assertHasEdge(node, fos, relation=HAS_COMPONENT)

        self.assertHasNode(jun, **jun)
        self.assertHasEdge(node, fos, relation=HAS_COMPONENT)

    def test_complex_list_long(self):
        statement = 'complexAbundance(proteinAbundance(HGNC:HBP1),geneAbundance(HGNC:NCF1))'
        self.parser.complex_abundances.parseString(statement)


class TestComposite(TestTokenParserBase):
    """Tests the parsing of the composite function

    .. seealso::

            `BEL 2.0 Specification 2.1.3 <http://openbel.org/language/version_2.0/bel_specification_version_2.0.html#XcompositeA>`_
    """

    def setUp(self):
        self.parser.clear()
        self.parser.composite_abundance.setParseAction(self.parser.handle_term)

    def test_213a(self):
        """Evidence: ``IL-6 and IL-23 synergistically induce Th17 differentiation"""
        statement = 'composite(p(HGNC:IL6), complex(GOCC:"interleukin-23 complex"))'
        result = self.parser.composite_abundance.parseString(statement)

        expected_result = [COMPOSITE, [PROTEIN, 'HGNC', 'IL6'], [COMPLEX, 'GOCC', 'interleukin-23 complex']]
        self.assertEqual(expected_result, result.asList())

        expected_dict = {
            FUNCTION: COMPOSITE,
            MEMBERS: [
                {
                    FUNCTION: PROTEIN,
                    NAMESPACE: 'HGNC',
                    NAME: 'IL6'
                }, {
                    FUNCTION: COMPLEX,
                    NAMESPACE: 'GOCC',
                    NAME: 'interleukin-23 complex'
                }
            ]
        }
        self.assertEqual(expected_dict, result.asDict())

        il23 = named_complex_abundance('GOCC', 'interleukin-23 complex')
        il6 = protein('HGNC', 'IL6')
        expected_node = composite_abundance([il23, il6])

        # needs sorting
        self.assertEqual('composite(complex(GOCC:"interleukin-23 complex"), p(HGNC:IL6))', expected_node.as_bel())

        self.assertEqual(3, self.parser.graph.number_of_nodes())
        self.assertEqual(2, self.parser.graph.number_of_edges())

        self.assertHasNode(il23, **il23)
        self.assertHasEdge(expected_node, il23, relation=HAS_COMPONENT)

        self.assertHasNode(il6, **il6)
        self.assertHasEdge(expected_node, il6, relation=HAS_COMPONENT)


class TestBiologicalProcess(TestTokenParserBase):
    def setUp(self):
        self.parser.clear()
        self.parser.biological_process.setParseAction(self.parser.handle_term)

    def test_231a(self):
        """"""
        statement = 'bp(GOBP:"cell cycle arrest")'
        result = self.parser.biological_process.parseString(statement)

        expected_result = [BIOPROCESS, 'GOBP', 'cell cycle arrest']
        self.assertEqual(expected_result, result.asList())

        expected_dict = {
            FUNCTION: BIOPROCESS,
            NAMESPACE: 'GOBP',
            NAME: 'cell cycle arrest'
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = bioprocess('GOBP', 'cell cycle arrest')
        self.assertEqual(statement, expected_node.as_bel())
        self.assertHasNode(expected_node, **expected_node)


class TestPathology(TestTokenParserBase):
    def setUp(self):
        self.parser.clear()
        self.parser.pathology.setParseAction(self.parser.handle_term)

    def test_232a(self):
        statement = 'pathology(MESHD:adenocarcinoma)'
        result = self.parser.pathology.parseString(statement)

        expected_result = [PATHOLOGY, 'MESHD', 'adenocarcinoma']
        self.assertEqual(expected_result, result.asList())

        expected_dict = {
            FUNCTION: PATHOLOGY,
            NAMESPACE: 'MESHD',
            NAME: 'adenocarcinoma'
        }
        self.assertEqual(expected_dict, result.asDict())

        expected_node = pathology('MESHD', 'adenocarcinoma')
        self.assertEqual('path(MESHD:adenocarcinoma)', expected_node.as_bel())
        self.assertHasNode(expected_node, **expected_node)


class TestActivity(TestTokenParserBase):
    def setUp(self):
        self.parser.clear()
        self.parser.activity.setParseAction(self.parser.handle_term)

    def test_activity_bare(self):
        """"""
        statement = 'act(p(HGNC:AKT1))'
        result = self.parser.activity.parseString(statement)

        expected_result = [ACTIVITY, [PROTEIN, 'HGNC', 'AKT1']]
        self.assertEqual(expected_result, result.asList())

        mod = modifier_po_to_dict(result)
        expected_mod = {
            MODIFIER: ACTIVITY
        }
        self.assertEqual(expected_mod, mod)

    def test_activity_withMolecularActivityDefault(self):
        """Tests activity modifier with molecular activity from default BEL namespace"""
        statement = 'act(p(HGNC:AKT1), ma(kin))'
        result = self.parser.activity.parseString(statement)

        expected_dict = {
            MODIFIER: ACTIVITY,
            EFFECT: {
                NAME: 'kin',
                NAMESPACE: BEL_DEFAULT_NAMESPACE
            },
            TARGET: {
                FUNCTION: PROTEIN,
                NAMESPACE: 'HGNC',
                NAME: 'AKT1'
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        mod = modifier_po_to_dict(result)
        expected_mod = {
            MODIFIER: ACTIVITY,
            EFFECT: {
                NAME: 'kin',
                NAMESPACE: BEL_DEFAULT_NAMESPACE
            }
        }
        self.assertEqual(expected_mod, mod)

    def test_activity_withMolecularActivityDefaultLong(self):
        """Tests activity modifier with molecular activity from custom namespaced"""
        statement = 'act(p(HGNC:AKT1), ma(catalyticActivity))'
        result = self.parser.activity.parseString(statement)

        expected_dict = {
            MODIFIER: ACTIVITY,
            EFFECT: {
                NAME: 'cat',
                NAMESPACE: BEL_DEFAULT_NAMESPACE
            },
            TARGET: {
                FUNCTION: PROTEIN,
                NAMESPACE: 'HGNC',
                NAME: 'AKT1'
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        mod = modifier_po_to_dict(result)
        expected_mod = {
            MODIFIER: ACTIVITY,
            EFFECT: {
                NAME: 'cat',
                NAMESPACE: BEL_DEFAULT_NAMESPACE
            }
        }
        self.assertEqual(expected_mod, mod)

    def test_activity_withMolecularActivityCustom(self):
        """Tests activity modifier with molecular activity from custom namespaced"""
        statement = 'act(p(HGNC:AKT1), ma(GOMF:"catalytic activity"))'
        result = self.parser.activity.parseString(statement)

        expected_dict = {
            MODIFIER: ACTIVITY,
            EFFECT: {
                NAMESPACE: 'GOMF',
                NAME: 'catalytic activity'
            },
            TARGET: {
                FUNCTION: PROTEIN,
                NAMESPACE: 'HGNC',
                NAME: 'AKT1'
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        mod = modifier_po_to_dict(result)
        expected_mod = {
            MODIFIER: ACTIVITY,
            EFFECT: {
                NAMESPACE: 'GOMF',
                NAME: 'catalytic activity'
            }
        }
        self.assertEqual(expected_mod, mod)

    def test_activity_legacy(self):
        """Test BEL 1.0 style molecular activity annotation"""
        statement = 'kin(p(HGNC:AKT1))'
        result = self.parser.activity.parseString(statement)

        expected_dict = {
            MODIFIER: ACTIVITY,
            EFFECT: {
                NAME: 'kin',
                NAMESPACE: BEL_DEFAULT_NAMESPACE
            },
            TARGET: {
                FUNCTION: PROTEIN,
                NAMESPACE: 'HGNC',
                NAME: 'AKT1'
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        mod = modifier_po_to_dict(result)
        expected_mod = {
            MODIFIER: ACTIVITY,
            EFFECT: {
                NAME: 'kin',
                NAMESPACE: BEL_DEFAULT_NAMESPACE
            }
        }
        self.assertEqual(expected_mod, mod)

        node = protein('HGNC', 'AKT1')
        self.assertHasNode(node)


class TestTranslocationPermissive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.graph = BELGraph()
        cls.parser = BelParser(cls.graph, allow_unqualified_translocations=True)

    def setUp(self):
        self.parser.clear()
        self.parser.transformation.setParseAction(self.parser.handle_term)

    def assertHasNode(self, member, **kwargs):
        assertHasNode(self, member, self.parser.graph, **kwargs)

    def assertHasEdge(self, u, v, **kwargs):
        assertHasEdge(self, u, v, self.parser.graph, **kwargs)

    def test_unqualified_translocation_single(self):
        """translocation example"""
        statement = 'tloc(p(HGNC:EGFR))'
        result = self.parser.transformation.parseString(statement)

        expected_dict = {
            MODIFIER: TRANSLOCATION,
            TARGET: {
                FUNCTION: PROTEIN,
                NAMESPACE: 'HGNC',
                NAME: 'EGFR'
            },
        }

        self.assertEqual(expected_dict, result.asDict())

        mod = modifier_po_to_dict(result)
        expected_mod = {
            MODIFIER: TRANSLOCATION,
        }
        self.assertEqual(expected_mod, mod)

        node = protein('HGNC', 'EGFR')
        self.assertHasNode(node)

    def test_unqualified_translocation_relation(self):
        """
        3.1.2 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#XdIncreases
        Test translocation in object
        """
        update_provenance(self.parser.control_parser)

        statement = 'a(ADO:"Abeta_42") => tloc(a(CHEBI:"calcium(2+)"))'
        result = self.parser.relation.parseString(statement)

        expected_dict = {
            SUBJECT: {
                FUNCTION: ABUNDANCE,
                NAMESPACE: 'ADO',
                NAME: 'Abeta_42'
            },
            RELATION: DIRECTLY_INCREASES,
            OBJECT: {
                TARGET: {
                    FUNCTION: ABUNDANCE,
                    NAMESPACE: 'CHEBI',
                    NAME: 'calcium(2+)'
                },
                MODIFIER: TRANSLOCATION,
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        sub = abundance('ADO', 'Abeta_42')
        self.assertHasNode(sub)

        obj = abundance('CHEBI', 'calcium(2+)')
        self.assertHasNode(obj)

        expected_annotations = {
            RELATION: DIRECTLY_INCREASES,
            OBJECT: {
                MODIFIER: TRANSLOCATION,
            }
        }

        self.assertHasEdge(sub, obj, **expected_annotations)


class TestTransformation(TestTokenParserBase):
    def setUp(self):
        self.parser.clear()
        self.parser.transformation.setParseAction(self.parser.handle_term)

    def test_degredation_short(self):
        """Test the short form of degradation works"""
        statement = 'deg(p(HGNC:AKT1))'
        result = self.parser.transformation.parseString(statement)

        expected_result = [DEGRADATION, [PROTEIN, 'HGNC', 'AKT1']]
        self.assertEqual(expected_result, result.asList())

        expected_dict = {
            MODIFIER: DEGRADATION,
            TARGET: {
                FUNCTION: PROTEIN,
                NAMESPACE: 'HGNC',
                NAME: 'AKT1'
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        mod = modifier_po_to_dict(result)
        expected_mod = {
            MODIFIER: DEGRADATION,
        }
        self.assertEqual(expected_mod, mod)

    def test_degradation_long(self):
        """Test the long form of degradation works"""
        statement = 'degradation(p(HGNC:EGFR))'
        result = self.parser.transformation.parseString(statement)

        expected_result = [DEGRADATION, [PROTEIN, 'HGNC', 'EGFR']]
        self.assertEqual(expected_result, result.asList())

        expected_dict = {
            MODIFIER: DEGRADATION,
            TARGET: {
                FUNCTION: PROTEIN,
                NAMESPACE: 'HGNC',
                NAME: 'EGFR'
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        mod = modifier_po_to_dict(result)
        expected_mod = {
            MODIFIER: DEGRADATION,
        }
        self.assertEqual(expected_mod, mod)

        node = protein('HGNC', 'EGFR')
        self.assertHasNode(node)

    def test_translocation_standard(self):
        """translocation example"""
        statement = 'tloc(p(HGNC:EGFR), fromLoc(GOCC:"cell surface"), toLoc(GOCC:endosome))'
        result = self.parser.transformation.parseString(statement)

        expected_dict = {
            MODIFIER: TRANSLOCATION,
            TARGET: {
                FUNCTION: PROTEIN,
                NAMESPACE: 'HGNC',
                NAME: 'EGFR'
            },
            EFFECT: {
                FROM_LOC: {NAMESPACE: 'GOCC', NAME: 'cell surface'},
                TO_LOC: {NAMESPACE: 'GOCC', NAME: 'endosome'}
            }
        }

        self.assertEqual(expected_dict, result.asDict())

        mod = modifier_po_to_dict(result)
        expected_mod = translocation(
            from_loc=entity(namespace='GOCC', name='cell surface'),
            to_loc=entity(namespace='GOCC', name='endosome'),
        )

        self.assertEqual(expected_mod, mod)

        node = protein('HGNC', 'EGFR')
        self.assertHasNode(node)

    def test_translocation_bare(self):
        """translocation example"""
        statement = 'tloc(p(HGNC:EGFR), GOCC:"cell surface", GOCC:endosome)'
        result = self.parser.transformation.parseString(statement)

        expected_dict = {
            MODIFIER: TRANSLOCATION,
            TARGET: {
                FUNCTION: PROTEIN,
                NAMESPACE: 'HGNC',
                NAME: 'EGFR'
            },
            EFFECT: {
                FROM_LOC: {NAMESPACE: 'GOCC', NAME: 'cell surface'},
                TO_LOC: {NAMESPACE: 'GOCC', NAME: 'endosome'}
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        mod = modifier_po_to_dict(result)
        expected_mod = {
            MODIFIER: TRANSLOCATION,
            EFFECT: {
                FROM_LOC: {NAMESPACE: 'GOCC', NAME: 'cell surface'},
                TO_LOC: {NAMESPACE: 'GOCC', NAME: 'endosome'}
            }
        }
        self.assertEqual(expected_mod, mod)

        node = protein('HGNC', 'EGFR')
        self.assertHasNode(node)

    def test_unqualified_translocation_strict(self):
        """Fail on an improperly written single argument translocation"""
        statement = 'tloc(a(NS:"T-Lymphocytes"))'
        with self.assertRaises(MalformedTranslocationWarning):
            self.parser.translocation.parseString(statement)

    def test_translocation_secretion(self):
        """cell secretion short form"""
        statement = 'sec(p(HGNC:EGFR))'
        result = self.parser.transformation.parseString(statement)

        expected_result = ['CellSecretion', [PROTEIN, 'HGNC', 'EGFR']]
        self.assertEqual(expected_result, result.asList())

        mod = modifier_po_to_dict(result)
        expected_mod = secretion()
        self.assertEqual(expected_mod, mod)

        node = protein('HGNC', 'EGFR')
        self.assertHasNode(node)

    def test_translocation_surface(self):
        """cell surface expression short form"""
        statement = 'surf(p(HGNC:EGFR))'
        result = self.parser.transformation.parseString(statement)

        expected_result = ['CellSurfaceExpression', [PROTEIN, 'HGNC', 'EGFR']]
        self.assertEqual(expected_result, result.asList())

        expected_mod = cell_surface_expression()
        self.assertEqual(expected_mod, modifier_po_to_dict(result))

        node = protein('HGNC', 'EGFR')
        self.assertHasNode(node)

    def test_reaction_1(self):
        """"""
        statement = 'rxn(reactants(a(CHEBI:superoxide)), products(a(CHEBI:"hydrogen peroxide"), a(CHEBI:oxygen)))'
        result = self.parser.transformation.parseString(statement)

        expected_result = [
            REACTION,
            [
                [ABUNDANCE, 'CHEBI', 'superoxide']
            ],
            [
                [ABUNDANCE, 'CHEBI', 'hydrogen peroxide'],
                [ABUNDANCE, 'CHEBI', 'oxygen']
            ]
        ]
        self.assertEqual(expected_result, result.asList())

        expected_dict = {
            FUNCTION: REACTION,
            REACTANTS: [
                {
                    FUNCTION: ABUNDANCE,
                    NAMESPACE: 'CHEBI',
                    NAME: 'superoxide'
                }
            ],
            PRODUCTS: [
                {
                    FUNCTION: ABUNDANCE,
                    NAMESPACE: 'CHEBI',
                    NAME: 'hydrogen peroxide'
                }, {

                    FUNCTION: ABUNDANCE,
                    NAMESPACE: 'CHEBI',
                    NAME: 'oxygen'
                }

            ]
        }
        self.assertEqual(expected_dict, result.asDict())

        superoxide_node = abundance('CHEBI', 'superoxide')
        hydrogen_peroxide = abundance('CHEBI', 'hydrogen peroxide')
        oxygen_node = abundance('CHEBI', 'oxygen')

        expected_node = reaction(reactants=superoxide_node, products=[hydrogen_peroxide, oxygen_node])
        self.assertEqual(statement, expected_node.as_bel())
        self.assertHasNode(expected_node)

        self.assertHasNode(superoxide_node)
        self.assertHasEdge(expected_node, superoxide_node, relation=HAS_REACTANT)

        self.assertHasNode(hydrogen_peroxide)
        self.assertHasEdge(expected_node, hydrogen_peroxide, relation=HAS_PRODUCT)

        self.assertHasNode(oxygen_node)
        self.assertHasEdge(expected_node, oxygen_node, relation=HAS_PRODUCT)

    def test_reaction_2(self):
        statement = 'rxn(reactants(p(HGNC:APP)), products(p(HGNC:APP, frag(672_713))))'
        result = self.parser.transformation.parseString(statement)

        app = protein('HGNC', 'APP')
        self.assertHasNode(app)

        amyloid_beta = protein('HGNC', 'APP', variants=fragment(672, 713))
        self.assertHasNode(amyloid_beta)

        expected_node = reaction(app, amyloid_beta)
        self.assertHasNode(expected_node)

        self.assertHasEdge(expected_node, app, relation=HAS_REACTANT)
        self.assertHasEdge(expected_node, amyloid_beta, relation=HAS_PRODUCT)

    def test_clearance(self):
        """Tests that after adding things, the graph and parser can be cleared properly"""
        s1 = 'surf(p(HGNC:EGFR))'
        s2 = 'rxn(reactants(a(CHEBI:superoxide)),products(a(CHEBI:"hydrogen peroxide"), a(CHEBI:"oxygen")))'

        self.parser.transformation.parseString(s1)
        self.parser.transformation.parseString(s2)
        self.assertGreater(self.parser.graph.number_of_nodes(), 0)
        self.assertGreater(self.parser.graph.number_of_edges(), 0)

        self.parser.clear()
        self.assertEqual(0, self.parser.graph.number_of_nodes())
        self.assertEqual(0, self.parser.graph.number_of_edges())
        self.assertEqual(0, len(self.parser.control_parser.annotations))
        self.assertEqual(0, len(self.parser.control_parser.citation))


class TestSemantics(unittest.TestCase):
    def test_lenient_semantic_no_failure(self):
        graph = BELGraph()
        parser = BelParser(graph, allow_naked_names=True)

        update_provenance(parser.control_parser)

        parser.bel_term.addParseAction(parser.handle_term)
        parser.bel_term.parseString('bp(ABASD)')

        node = bioprocess(namespace=DIRTY, name='ABASD')
        self.assertIn(node, graph)
