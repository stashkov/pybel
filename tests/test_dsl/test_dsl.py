# -*- coding: utf-8 -*-

"""This module tests the DSL"""

import itertools as itt

import unittest

from pybel.constants import (
    ABUNDANCE, BEL_DEFAULT_NAMESPACE, BIOPROCESS, COMPLEX, COMPOSITE, FRAGMENT, FUNCTION, GENE, IDENTIFIER, NAME,
    NAMESPACE, PATHOLOGY, PMOD, PROTEIN, REACTION, RNA,
)
from pybel.dsl import (
    abundance, bioprocess, complex_abundance, composite_abundance, fragment, fusion_range, gene, gene_fusion,
    gene_substitution, gmod, hgvs, hgvs_reference, hgvs_unspecified, mirna, missing_fusion_range,
    named_complex_abundance, pathology, pmod, protein, protein_deletion, protein_substitution, reaction, rna,
    rna_fusion,
)
from pybel.dsl.exc import InferCentralDogmaException, PyBELDSLException
from pybel.dsl.utils import entity
from pybel.utils import ensure_quotes
from tests.utils import n


def variant_to_bel(tokens):
    return tokens.as_bel()


def fusion_range_to_bel(tokens):
    return tokens.as_bel()


class TestDSL(unittest.TestCase):
    def test_str_has_name(self):
        namespace, name = n(), n()
        node = abundance(namespace=namespace, name=name)
        self.assertEqual('a({namespace}:{name})'.format(namespace=namespace, name=ensure_quotes(name)), node.as_bel())

    def test_str_has_identifier(self):
        namespace, identifier = n(), n()
        node = abundance(namespace=namespace, identifier=identifier)
        self.assertEqual(
            'a({namespace}:{identifier})'.format(namespace=namespace, identifier=ensure_quotes(identifier)),
            node.as_bel())

    def test_str_has_both(self):
        namespace, identifier = n(), n()
        node = abundance(namespace=namespace, identifier=identifier)
        self.assertEqual(
            'a({namespace}:{identifier})'.format(namespace=namespace, identifier=ensure_quotes(identifier)),
            node.as_bel())

    def test_as_tuple(self):
        namespace, name = n(), n()

        node_tuple = ABUNDANCE, namespace, name
        node = abundance(namespace=namespace, name=name)

        self.assertEqual(node_tuple, node.as_tuple())
        self.assertEqual(hash(node_tuple), hash(node))

    def test_complex_with_name(self):
        """Tests a what happens with a named complex

        .. code-block::

            complex(SCOMP:"9-1-1 Complex") hasComponent p(HGNC:HUS1)
            complex(SCOMP:"9-1-1 Complex") hasComponent p(HGNC:RAD1)
            complex(SCOMP:"9-1-1 Complex") hasComponent p(HGNC:RAD9A)

        """
        hus1 = protein(namespace='HGNC', name='HUS1')
        rad1 = protein(namespace='HGNC', name='RAD1')
        rad9a = protein(namespace='HGNC', name='RAD9A')
        members = [hus1, rad1, rad9a]

        nine_one_one = complex_abundance(members=members, namespace='SCOMP', name='9-1-1 Complex')

        node_tuple = (COMPLEX,) + tuple(member.as_tuple() for member in members)

        self.assertEqual(node_tuple, nine_one_one.as_tuple())
        self.assertEqual(hash(node_tuple), hash(nine_one_one))

    def test_missing_parent(self):
        app = protein(name='APP', namespace='HGNC')
        self.assertIsNone(app.get_parent(), msg='parent should be none for reference node')

    def test_get_parent(self):
        ab42 = protein(name='APP', namespace='HGNC', variants=[fragment(start=672, stop=713)])
        app = ab42.get_parent()
        self.assertEqual('p(HGNC:APP)', app.as_bel())


class TestInferCentralDogma(unittest.TestCase):
    def test_protein(self):
        namespace, name = n(), n()
        p = protein(namespace=namespace, name=name)
        r = p.get_rna()
        self.assertIsInstance(r, rna)
        self.assertIn(FUNCTION, r)
        self.assertEqual(RNA, r[FUNCTION])
        self.assertIn(NAMESPACE, r)
        self.assertEqual(namespace, r[NAMESPACE])
        self.assertIn(NAME, r)
        self.assertEqual(name, r[NAME])
        self.assertNotIn(IDENTIFIER, r)

    def test_protein_identifier(self):
        namespace, name, identifier = n(), n(), n()
        p = protein(namespace=namespace, name=name, identifier=identifier)
        r = p.get_rna()
        self.assertIsInstance(r, rna)
        self.assertIn(FUNCTION, r)
        self.assertEqual(RNA, r[FUNCTION])
        self.assertIn(NAMESPACE, r)
        self.assertEqual(namespace, r[NAMESPACE])
        self.assertIn(NAME, r)
        self.assertEqual(name, r[NAME])
        self.assertIn(IDENTIFIER, r)
        self.assertEqual(identifier, r[IDENTIFIER])

    def test_protein_variants(self):
        namespace, name, hgvs_variant = n(), n(), n()
        p = protein(namespace=namespace, name=name, variants=hgvs(hgvs_variant))
        with self.assertRaises(InferCentralDogmaException):
            p.get_rna()

    def test_rna(self):
        namespace, name = n(), n()
        p = rna(namespace=namespace, name=name)
        r = p.get_gene()
        self.assertIsInstance(r, gene)
        self.assertIn(FUNCTION, r)
        self.assertEqual(GENE, r[FUNCTION])
        self.assertIn(NAMESPACE, r)
        self.assertEqual(namespace, r[NAMESPACE])
        self.assertIn(NAME, r)
        self.assertEqual(name, r[NAME])
        self.assertNotIn(IDENTIFIER, r)

    def test_mirna(self):
        namespace, name = n(), n()
        p = mirna(namespace=namespace, name=name)
        r = p.get_gene()
        self.assertIsInstance(r, gene)
        self.assertIn(FUNCTION, r)
        self.assertEqual(GENE, r[FUNCTION])
        self.assertIn(NAMESPACE, r)
        self.assertEqual(namespace, r[NAMESPACE])
        self.assertIn(NAME, r)
        self.assertEqual(name, r[NAME])
        self.assertNotIn(IDENTIFIER, r)

    def test_rna_variants(self):
        namespace, name, hgvs_variant = n(), n(), n()
        p = rna(namespace=namespace, name=name, variants=hgvs(hgvs_variant))
        with self.assertRaises(InferCentralDogmaException):
            p.get_gene()


class TestSort(unittest.TestCase):
    def test_sort_variants(self):
        a, b, c = hgvs('c.1521_1523delCTT'), hgvs('c.308G>A'), hgvs('p.Phe508del')
        for l in itt.permutations([a, b, c]):
            self.assertEqual([a, b, c], sorted(l))


class TestCanonicalize(unittest.TestCase):
    def test_entity_dsl_error(self):
        """"""
        with self.assertRaises(PyBELDSLException):
            entity(namespace=n(), name=None, identifier=None)

    def test_canonicalize_variant(self):
        self.assertEqual('var(p.Val600Glu)', variant_to_bel(hgvs('p.Val600Glu')))
        self.assertEqual('var(p.Val600Glu)', variant_to_bel(protein_substitution('Val', 600, 'Glu')))
        self.assertEqual('var(c.600C>G)', variant_to_bel(gene_substitution('C', 600, 'G')))
        self.assertEqual('var(p.Val600del)', variant_to_bel(protein_deletion('Val', 600)))
        self.assertEqual('var(?)', variant_to_bel(hgvs_unspecified()))
        self.assertEqual('var(=)', variant_to_bel(hgvs_reference()))

        self.assertEqual('pmod(Ph)', variant_to_bel(pmod('Ph')))
        self.assertEqual('pmod(TEST:Ph)', variant_to_bel(pmod('Ph', namespace='TEST')))
        self.assertEqual('pmod(TEST:Ph, Ser)', variant_to_bel(pmod('Ph', namespace='TEST', code='Ser')))
        self.assertEqual('pmod(TEST:Ph, Ser, 5)', variant_to_bel(pmod('Ph', namespace='TEST', code='Ser', position=5)))
        self.assertEqual('pmod(GO:"protein phosphorylation", Thr, 308)',
                         variant_to_bel(pmod(name='protein phosphorylation', namespace='GO', code='Thr', position=308)))

        self.assertEqual('frag(?)', variant_to_bel(fragment()))
        self.assertEqual('frag(672_713)', variant_to_bel(fragment(start=672, stop=713)))
        self.assertEqual('frag(672_?)', variant_to_bel(fragment(start=672)))
        self.assertEqual('frag(?_713)', variant_to_bel(fragment(stop=713)))
        self.assertEqual('frag(?, "descr")', variant_to_bel(fragment(description='descr')))
        self.assertEqual('frag(672_713, "descr")', variant_to_bel(fragment(start=672, stop=713, description='descr')))

        self.assertEqual('gmod(Me)', variant_to_bel(gmod('Me')))
        self.assertEqual('gmod(TEST:Me)', variant_to_bel(gmod('Me', namespace='TEST')))
        self.assertEqual('gmod(GO:"DNA Methylation")', variant_to_bel(gmod('DNA Methylation', namespace='GO')))

    def test_canonicalize_variant_dsl(self):
        """Uses the __str__ functions in the DSL to create BEL instead of external pybel.canonicalize"""
        self.assertEqual('var(p.Val600Glu)', str(hgvs('p.Val600Glu')))
        self.assertEqual('var(p.Val600Glu)', str(protein_substitution('Val', 600, 'Glu')))

        self.assertEqual('pmod(Ph)', str(pmod('Ph')))
        self.assertEqual('pmod(TEST:Ph)', str(pmod('Ph', namespace='TEST')))
        self.assertEqual('pmod(TEST:Ph, Ser)', str(pmod('Ph', namespace='TEST', code='Ser')))
        self.assertEqual('pmod(TEST:Ph, Ser, 5)', str(pmod('Ph', namespace='TEST', code='Ser', position=5)))
        self.assertEqual('pmod(GO:"protein phosphorylation", Thr, 308)',
                         str(pmod(name='protein phosphorylation', namespace='GO', code='Thr', position=308)))

        self.assertEqual('frag(?)', str(fragment()))
        self.assertEqual('frag(672_713)', str(fragment(start=672, stop=713)))
        self.assertEqual('frag(?, "descr")', str(fragment(description='descr')))
        self.assertEqual('frag(672_713, "descr")', str(fragment(start=672, stop=713, description='descr')))

        self.assertEqual('gmod(Me)', str(gmod('Me')))
        self.assertEqual('gmod(TEST:Me)', str(gmod('Me', namespace='TEST')))
        self.assertEqual('gmod(GO:"DNA Methylation")', str(gmod('DNA Methylation', namespace='GO')))

    def test_canonicalize_fusion_range(self):
        self.assertEqual('p.1_15', fusion_range_to_bel(fusion_range('p', 1, 15)))
        self.assertEqual('p.*_15', fusion_range_to_bel(fusion_range('p', '*', 15)))
        self.assertEqual('?', fusion_range_to_bel(missing_fusion_range()))

    def test_canonicalize_fusion_range_dsl(self):
        self.assertEqual('p.1_15', str(fusion_range('p', 1, 15)))
        self.assertEqual('p.*_15', str(fusion_range('p', '*', 15)))

    def test_abundance(self):
        short = abundance(namespace='CHEBI', name='water')
        self.assertEqual('a(CHEBI:water)', str(short))
        self.assertEqual((ABUNDANCE, 'CHEBI', 'water'), short.as_tuple())

        long = abundance(namespace='CHEBI', name='test name')
        self.assertEqual('a(CHEBI:"test name")', str(long))
        self.assertEqual((ABUNDANCE, 'CHEBI', 'test name'), long.as_tuple())

    def test_gene_reference(self):
        node = gene(namespace='EGID', name='780')
        self.assertEqual('g(EGID:780)', str(node))
        self.assertEqual((GENE, 'EGID', '780'), node.as_tuple())

    def test_protein_reference(self):
        self.assertEqual('p(HGNC:AKT1)', str(protein(namespace='HGNC', name='AKT1')))

    def test_protein_pmod(self):
        node = protein(name='PLCG1', namespace='HGNC', variants=pmod(name='Ph', code='Tyr'))
        self.assertEqual('p(HGNC:PLCG1, pmod(Ph, Tyr))', str(node))
        self.assertEqual((PROTEIN, 'HGNC', 'PLCG1', (PMOD, (BEL_DEFAULT_NAMESPACE, 'Ph'), 'Tyr')), node.as_tuple())

    def test_protein_fragment(self):
        node = protein(name='APP', namespace='HGNC', variants=fragment(start=672, stop=713))
        self.assertEqual('p(HGNC:APP, frag(672_713))', str(node))
        self.assertEqual((PROTEIN, 'HGNC', 'APP', (FRAGMENT, (672, 713))), node.as_tuple())

    def test_mirna_reference(self):
        self.assertEqual('m(HGNC:MIR1)', str(mirna(namespace='HGNC', name='MIR1')))

    def test_rna_fusion_specified(self):
        node = rna_fusion(
            partner5p=rna(namespace='HGNC', name='TMPRSS2'),
            range5p=fusion_range('r', 1, 79),
            partner3p=rna(namespace='HGNC', name='ERG'),
            range3p=fusion_range('r', 312, 5034)
        )
        self.assertEqual('r(fus(HGNC:TMPRSS2, r.1_79, HGNC:ERG, r.312_5034))', str(node))

    def test_rna_fusion_unspecified(self):
        node = rna_fusion(
            partner5p=rna(namespace='HGNC', name='TMPRSS2'),
            partner3p=rna(namespace='HGNC', name='ERG'),
        )
        self.assertEqual('r(fus(HGNC:TMPRSS2, ?, HGNC:ERG, ?))', str(node))

        t = RNA, ('HGNC', 'TMPRSS2'), ('?',), ('HGNC', 'ERG'), ('?',)
        self.assertEqual(t, node.as_tuple())

    def test_gene_fusion_specified(self):
        node = gene_fusion(
            partner5p=gene(namespace='HGNC', name='TMPRSS2'),
            range5p=fusion_range('c', 1, 79),
            partner3p=gene(namespace='HGNC', name='ERG'),
            range3p=fusion_range('c', 312, 5034)
        )

        self.assertEqual('g(fus(HGNC:TMPRSS2, c.1_79, HGNC:ERG, c.312_5034))', str(node))
        t = GENE, ('HGNC', 'TMPRSS2'), ('c', 1, 79), ('HGNC', 'ERG'), ('c', 312, 5034)
        self.assertEqual(t, node.as_tuple())

    def test_pathology(self):
        node = pathology(namespace='DO', name='Alzheimer disease')
        self.assertEqual('path(DO:"Alzheimer disease")', str(node))
        self.assertEqual((PATHOLOGY, 'DO', 'Alzheimer disease'), node.as_tuple())

    def test_bioprocess(self):
        node = bioprocess(namespace='GO', name='apoptosis')
        self.assertEqual('bp(GO:apoptosis)', str(node))
        self.assertEqual((BIOPROCESS, 'GO', 'apoptosis'), node.as_tuple())

    def test_named_complex_abundance(self):
        node = named_complex_abundance(namespace='SCOMP', name='Calcineurin Complex')
        self.assertEqual('complex(SCOMP:"Calcineurin Complex")', str(node))
        self.assertEqual((COMPLEX, 'SCOMP', 'Calcineurin Complex'), node.as_tuple())

    def test_complex_abundance(self):
        node = complex_abundance(members=[protein(namespace='HGNC', name='FOS'), protein(namespace='HGNC', name='JUN')])
        t = COMPLEX, (PROTEIN, 'HGNC', 'FOS'), (PROTEIN, 'HGNC', 'JUN')
        self.assertEqual('complex(p(HGNC:FOS), p(HGNC:JUN))', str(node))
        self.assertEqual(t, node.as_tuple())

    def test_composite_abundance(self):
        node = composite_abundance(members=[
            protein(namespace='HGNC', name='FOS'),
            protein(namespace='HGNC', name='JUN')
        ])
        t = COMPOSITE, (PROTEIN, 'HGNC', 'FOS'), (PROTEIN, 'HGNC', 'JUN')
        self.assertEqual('composite(p(HGNC:FOS), p(HGNC:JUN))', str(node))
        self.assertEqual(t, node.as_tuple())

    def test_reaction_single(self):
        node = reaction(
            reactants=abundance(namespace='CHEBI', name='A'),
            products=abundance(namespace='CHEBI', name='B')
        )
        t = REACTION, ((ABUNDANCE, 'CHEBI', 'A'),), ((ABUNDANCE, 'CHEBI', 'B'),)
        self.assertEqual('rxn(reactants(a(CHEBI:A)), products(a(CHEBI:B)))', str(node))
        self.assertEqual(t, node.as_tuple())

    def test_reaction_multi(self):
        node = reaction(
            reactants=[abundance(namespace='CHEBI', name='A'), abundance(namespace='CHEBI', name='C')],
            products=[abundance(namespace='CHEBI', name='B'), abundance(namespace='CHEBI', name='D')]
        )
        t = REACTION, ((ABUNDANCE, 'CHEBI', 'A'), (ABUNDANCE, 'CHEBI', 'C')), (
            (ABUNDANCE, 'CHEBI', 'B'), (ABUNDANCE, 'CHEBI', 'D'))
        self.assertEqual('rxn(reactants(a(CHEBI:A), a(CHEBI:C)), products(a(CHEBI:B), a(CHEBI:D)))', str(node))
        self.assertEqual(t, node.as_tuple())
