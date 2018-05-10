# -*- coding: utf-8 -*-

"""Tests parsing of BEL relations"""

import logging
import unittest

from pyparsing import ParseException

from pybel import BELGraph
from pybel.constants import *
from pybel.dsl import (
    abundance, activity, bioprocess, complex_abundance, composite_abundance, entity, gene, gene_substitution, gmod,
    hgvs, named_complex_abundance, pathology, pmod, protein, reaction, rna,
)
from pybel.parser import BelParser
from pybel.parser.exc import (
    MissingNamespaceNameWarning, NestedRelationWarning, RelabelWarning, UndefinedNamespaceWarning,
)
from pybel.tokens import dict_to_entity
from tests.constants import TestTokenParserBase, test_citation_dict, test_evidence_text

log = logging.getLogger(__name__)

psoriasis = pathology(namespace='MESH', name='Psoriasis')
skin_diseases = pathology(namespace='MESH', name='Skin Diseases')


class TestRelations(TestTokenParserBase):
    @classmethod
    def setUpClass(cls):
        """Sets up the token parser and streamlines it"""
        super(TestRelations, cls).setUpClass()
        cls.parser.relation.streamline()

    def setUp(self):
        """Adds the default provenance at the beginning of each test"""
        super(TestRelations, self).setUp()
        self.add_default_provenance()

    def test_ensure_no_dup_nodes(self):
        """Ensure node isn't added twice, even if from different statements"""
        self.parser.gene.addParseAction(self.parser.handle_term)
        result = self.parser.bel_term.parseString('g(HGNC:AKT1)')

        expected_result_dict = {
            FUNCTION: GENE,
            NAMESPACE: 'HGNC',
            NAME: 'AKT1'
        }

        self.assertEqual(expected_result_dict, result.asDict())

        self.parser.degradation.addParseAction(self.parser.handle_term)
        self.parser.degradation.parseString('deg(g(HGNC:AKT1))')

        node = gene('HGNC', 'AKT1')

        self.assertNumberNodes(1)
        self.assertHasNode(node)

    def test_singleton(self):
        """Test singleton composite in subject."""
        statement = 'composite(p(HGNC:CASP8),p(HGNC:FADD),a(ADO:"Abeta_42"))'
        result = self.parser.statement.parseString(statement)

        expected = [
            COMPOSITE,
            [PROTEIN, 'HGNC', 'CASP8'],
            [PROTEIN, 'HGNC', 'FADD'],
            [ABUNDANCE, 'ADO', 'Abeta_42']
        ]
        self.assertEqual(expected, result.asList())

        amyloid_beta = abundance('ADO', 'Abeta_42')
        self.assertHasNode(amyloid_beta)

        casp8 = protein('HGNC', 'CASP8')
        self.assertHasNode(casp8)

        fadd = protein('HGNC', 'FADD')
        self.assertHasNode(fadd)

        node = composite_abundance([amyloid_beta, casp8, fadd])
        self.assertHasNode(node)

        self.assertHasEdge(node, amyloid_beta, relation=HAS_COMPONENT)
        self.assertHasEdge(node, casp8, relation=HAS_COMPONENT)
        self.assertHasEdge(node, fadd, relation=HAS_COMPONENT)

    def test_predicate_failure(self):
        """Checks that if there's a problem with the relation/object, that an error gets thrown"""
        statement = 'composite(p(HGNC:CASP8),p(HGNC:FADD),a(ADO:"Abeta_42")) -> nope(GOBP:"neuron apoptotic process")'

        with self.assertRaises(ParseException):
            self.parser.relation.parseString(statement)

    def test_increases(self):
        """Test composite in subject. See BEL 2.0 specification
        `3.1.1 <http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#Xincreases>`_
        """
        statement = 'composite(p(HGNC:CASP8),p(HGNC:FADD),a(ADO:"Abeta_42")) -> bp(GOBP:"neuron apoptotic process")'
        result = self.parser.relation.parseString(statement)

        expected = [
            [COMPOSITE, [PROTEIN, 'HGNC', 'CASP8'], [PROTEIN, 'HGNC', 'FADD'],
             [ABUNDANCE, 'ADO', 'Abeta_42']],
            INCREASES,
            [BIOPROCESS, 'GOBP', 'neuron apoptotic process']
        ]
        self.assertEqual(expected, result.asList())

        self.assertNumberNodes(5)
        self.assertNumberEdges(4)

        amyloid_beta = abundance('ADO', 'Abeta_42')
        self.assertHasNode(amyloid_beta)

        casp8 = protein('HGNC', 'CASP8')
        self.assertHasNode(casp8)

        fadd = protein('HGNC', 'FADD')
        self.assertHasNode(fadd)

        casp8_and_fadd_and_amyloid_beta = composite_abundance([casp8, fadd, amyloid_beta])
        self.assertHasNode(casp8_and_fadd_and_amyloid_beta)

        self.assertHasEdge(casp8_and_fadd_and_amyloid_beta, casp8, relation=HAS_COMPONENT)
        self.assertHasEdge(casp8_and_fadd_and_amyloid_beta, fadd, relation=HAS_COMPONENT)

        neuron_apoptotic_process = bioprocess('GOBP', 'neuron apoptotic process')
        self.assertHasNode(neuron_apoptotic_process)

        edge_data = {
            RELATION: INCREASES,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text
        }
        self.assertHasEdge(casp8_and_fadd_and_amyloid_beta, neuron_apoptotic_process, **edge_data)

    def test_directlyIncreases_withTlocObject(self):
        """Test translocation in object. See BEL 2.0 specification
        `3.1.2 <http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#XdIncreases>`_
        """
        statement = 'a(ADO:"Abeta_42") => tloc(a(CHEBI:"calcium(2+)"),fromLoc(MESHCS:"Cell Membrane"),' \
                    'toLoc(MESHCS:"Intracellular Space"))'
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
                EFFECT: {
                    FROM_LOC: {NAMESPACE: 'MESHCS', NAME: 'Cell Membrane'},
                    TO_LOC: {NAMESPACE: 'MESHCS', NAME: 'Intracellular Space'}
                }
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        amyloid_beta_42 = abundance('ADO', 'Abeta_42')
        self.assertHasNode(amyloid_beta_42)

        calcium_2_ion = abundance('CHEBI', 'calcium(2+)')
        self.assertHasNode(calcium_2_ion)

        expected_annotations = {
            RELATION: DIRECTLY_INCREASES,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text,
            OBJECT: {
                MODIFIER: TRANSLOCATION,
                EFFECT: {
                    FROM_LOC: {NAMESPACE: 'MESHCS', NAME: 'Cell Membrane'},
                    TO_LOC: {NAMESPACE: 'MESHCS', NAME: 'Intracellular Space'}
                }
            }
        }

        self.assertHasEdge(amyloid_beta_42, calcium_2_ion, **expected_annotations)

    def test_decreases(self):
        """
        3.1.3 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#Xdecreases
        Test decreases with reaction"""
        statement = ('pep(p(SFAM:"CAPN Family", location(GOCC:intracellular))) -| '
                     'reaction(reactants(p(HGNC:CDK5R1)),products(p(HGNC:CDK5)))')
        result = self.parser.relation.parseString(statement)

        expected_dict = {
            SUBJECT: {
                MODIFIER: ACTIVITY,
                TARGET: {
                    FUNCTION: PROTEIN,
                    NAMESPACE: 'SFAM',
                    NAME: 'CAPN Family',
                    LOCATION: {NAMESPACE: 'GOCC', NAME: 'intracellular'}
                },
                EFFECT: {
                    NAME: 'pep',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                },
            },
            RELATION: 'decreases',
            OBJECT: {
                FUNCTION: REACTION,
                REACTANTS: [
                    {FUNCTION: PROTEIN, NAMESPACE: 'HGNC', NAME: 'CDK5R1'}
                ],
                PRODUCTS: [
                    {FUNCTION: PROTEIN, NAMESPACE: 'HGNC', NAME: 'CDK5'}
                ]

            }
        }
        self.assertEqual(expected_dict, result.asDict())

        sub = protein('SFAM', 'CAPN Family')
        self.assertHasNode(sub)

        cdk5r1 = protein('HGNC', 'CDK5R1')
        cdk5 = protein('HGNC', 'CDK5')
        obj = reaction(cdk5r1, cdk5)

        self.assertHasNode(obj)
        self.assertHasNode(cdk5r1)
        self.assertHasNode(cdk5)
        self.assertHasEdge(obj, cdk5r1, relation=HAS_REACTANT)
        self.assertHasEdge(obj, cdk5, relation=HAS_PRODUCT)

        expected_edge_attributes = {
            RELATION: DECREASES,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text,
            SUBJECT: {
                MODIFIER: ACTIVITY,
                EFFECT: {
                    NAME: 'pep',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                },
                LOCATION: {NAMESPACE: 'GOCC', NAME: 'intracellular'}
            }
        }

        # check that the activity function works correctly
        self.assertEqual(expected_edge_attributes[SUBJECT],
                         activity(name='pep', location=entity(name='intracellular', namespace='GOCC')))

        self.assertHasEdge(sub, obj, **expected_edge_attributes)

    def test_directlyDecreases(self):
        """3.1.4 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#XdDecreases
        Tests simple triple"""
        statement = ('proteinAbundance(HGNC:CAT, location(GOCC:intracellular)) directlyDecreases '
                     'abundance(CHEBI:"hydrogen peroxide")')
        result = self.parser.relation.parseString(statement)

        expected_dict = {
            SUBJECT: {
                FUNCTION: PROTEIN,
                NAMESPACE: 'HGNC',
                NAME: 'CAT',
                LOCATION: {NAMESPACE: 'GOCC', NAME: 'intracellular'}
            },
            RELATION: 'directlyDecreases',
            OBJECT: {
                FUNCTION: ABUNDANCE,
                NAMESPACE: 'CHEBI',
                NAME: 'hydrogen peroxide'
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        sub = protein('HGNC', 'CAT')
        self.assertHasNode(sub)

        obj = abundance('CHEBI', 'hydrogen peroxide')
        self.assertHasNode(obj)

        expected_attrs = {
            RELATION: DIRECTLY_DECREASES,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text,
            SUBJECT: {
                LOCATION: {NAMESPACE: 'GOCC', NAME: 'intracellular'}
            },
        }
        self.assertHasEdge(sub, obj, **expected_attrs)

    def test_directlyDecreases_annotationExpansion(self):
        """
        3.1.4 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#XdDecreases
        Tests simple triple"""
        statement = 'g(HGNC:CAT, location(GOCC:intracellular)) directlyDecreases abundance(CHEBI:"hydrogen peroxide")'

        annotations = {
            'ListAnnotation': {'a', 'b'},
            'ScalarAnnotation': {'c'}
        }

        self.parser.control_parser.annotations.update(annotations)

        result = self.parser.relation.parseString(statement)

        expected_dict = {
            SUBJECT: {
                FUNCTION: GENE,
                NAMESPACE: 'HGNC',
                NAME: 'CAT',
                LOCATION: {
                    NAMESPACE: 'GOCC',
                    NAME: 'intracellular'
                }
            },
            RELATION: DIRECTLY_DECREASES,
            OBJECT: {
                FUNCTION: ABUNDANCE,
                NAMESPACE: 'CHEBI',
                NAME: 'hydrogen peroxide'
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        sub = gene('HGNC', 'CAT')
        self.assertHasNode(sub)

        obj = abundance('CHEBI', 'hydrogen peroxide')
        self.assertHasNode(obj)

        expected_attrs = {
            SUBJECT: {
                LOCATION: {
                    NAMESPACE: 'GOCC',
                    NAME: 'intracellular'
                }
            },
            RELATION: DIRECTLY_DECREASES,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text,
            ANNOTATIONS: {
                'ListAnnotation': {'a': True, 'b': True},
                'ScalarAnnotation': {'c': True}
            }
        }
        self.assertHasEdge(sub, obj, **expected_attrs)

    def test_rateLimitingStepOf_subjectActivity(self):
        """3.1.5 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#_ratelimitingstepof"""
        statement = 'act(p(HGNC:HMGCR), ma(cat)) rateLimitingStepOf bp(GOBP:"cholesterol biosynthetic process")'
        result = self.parser.relation.parseString(statement)

        expected_dict = {
            SUBJECT: {
                MODIFIER: ACTIVITY,
                TARGET: {
                    FUNCTION: PROTEIN,
                    NAMESPACE: 'HGNC',
                    NAME: 'HMGCR'
                },
                EFFECT: {
                    NAME: 'cat',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                },
            },
            RELATION: RATE_LIMITING_STEP_OF,
            OBJECT: {
                FUNCTION: BIOPROCESS,
                NAMESPACE: 'GOBP',
                NAME: 'cholesterol biosynthetic process'
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        sub = protein('HGNC', 'HMGCR')
        self.assertHasNode(sub)

        obj = bioprocess('GOBP', 'cholesterol biosynthetic process')
        self.assertHasNode(obj)

        edge_data = {
            RELATION: RATE_LIMITING_STEP_OF,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text,
            SUBJECT: {
                MODIFIER: ACTIVITY,
                EFFECT: {
                    NAME: 'cat',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                },
            },
        }
        self.assertHasEdge(sub, obj, **edge_data)

    def test_cnc_withSubjectVariant(self):
        """
        3.1.6 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#Xcnc
        Test SNP annotation
        """
        statement = 'g(HGNC:APP,sub(G,275341,C)) cnc path(MESHD:"Alzheimer Disease")'
        result = self.parser.relation.parseString(statement)

        expected_dict = {
            SUBJECT: {
                FUNCTION: GENE,
                NAMESPACE: 'HGNC',
                NAME: 'APP',
                VARIANTS: [
                    {
                        KIND: HGVS,
                        IDENTIFIER: 'c.275341G>C'
                    }
                ]
            },
            RELATION: CAUSES_NO_CHANGE,
            OBJECT: {
                FUNCTION: PATHOLOGY,
                NAMESPACE: 'MESHD',
                NAME: 'Alzheimer Disease'
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        sub = gene('HGNC', 'APP', variants=gene_substitution('G', 275341, 'C'))
        self.assertHasNode(sub)

        obj = pathology('MESHD', 'Alzheimer Disease')
        self.assertHasNode(obj)

        self.assertHasEdge(sub, obj, **{
            RELATION: CAUSES_NO_CHANGE,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text,
        })

    def test_regulates_multipleAnnotations(self):
        """
        3.1.7 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#_regulates_reg
        Test nested definitions"""
        statement = 'pep(complex(p(HGNC:F3),p(HGNC:F7))) regulates pep(p(HGNC:F9))'
        result = self.parser.relation.parseString(statement)

        expected_dict = {
            SUBJECT: {
                MODIFIER: ACTIVITY,
                EFFECT: {
                    NAME: 'pep',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                },
                TARGET: {
                    FUNCTION: COMPLEX,
                    MEMBERS: [
                        {FUNCTION: PROTEIN, NAMESPACE: 'HGNC', NAME: 'F3'},
                        {FUNCTION: PROTEIN, NAMESPACE: 'HGNC', NAME: 'F7'}
                    ]
                }
            },
            RELATION: REGULATES,
            OBJECT: {
                MODIFIER: ACTIVITY,
                EFFECT: {
                    NAME: 'pep',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                },
                TARGET: {
                    FUNCTION: PROTEIN,
                    NAMESPACE: 'HGNC',
                    NAME: 'F9'
                }

            }
        }
        self.assertEqual(expected_dict, result.asDict())

        f3 = protein('HGNC', 'F3')
        self.assertHasNode(f3)

        f7 = protein('HGNC', 'F7')
        self.assertHasNode(f7)

        sub = complex_abundance([f3, f7])

        self.assertHasEdge(sub, f3, relation=HAS_COMPONENT)
        self.assertHasEdge(sub, f7, relation=HAS_COMPONENT)

        f9 = protein('HGNC', 'F9')
        self.assertHasNode(f9)

        edge_data = {
            RELATION: REGULATES,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text,
            SUBJECT: {
                MODIFIER: ACTIVITY,
                EFFECT: {
                    NAME: 'pep',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                },
            },
            OBJECT: {
                MODIFIER: ACTIVITY,
                EFFECT: {
                    NAME: 'pep',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                },
            }
        }
        self.assertHasEdge(sub, f9, **edge_data)

    def test_nested_failure(self):
        """
        3.1 \
        Test nested statement"""
        statement = 'p(HGNC:CAT) -| (a(CHEBI:"hydrogen peroxide") -> bp(GO:"apoptotic process"))'
        with self.assertRaises(NestedRelationWarning):
            self.parser.relation.parseString(statement)

    def test_nested_lenient(self):
        """ 3.1 \ Test nested statement"""
        statement = 'p(HGNC:CAT) -| (a(CHEBI:"hydrogen peroxide") -> bp(GO:"apoptotic process"))'
        self.parser.allow_nested = True

        self.parser.relation.parseString(statement)

        cat = protein('HGNC', 'CAT')
        hydrogen_peroxide = abundance('CHEBI', "hydrogen peroxide")
        apoptotic_process = bioprocess('GO', "apoptotic process")

        self.assertHasEdge(cat, hydrogen_peroxide, **{
            RELATION: DECREASES,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text,
        })
        self.assertHasEdge(hydrogen_peroxide, apoptotic_process, **{
            RELATION: INCREASES,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text,
        })

        self.parser.lenient = False

    def test_negativeCorrelation_withObjectVariant(self):
        """Tests the phosphoralation tag

        3.2.1 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#XnegCor
        """
        statement = 'kin(p(SFAM:"GSK3 Family")) neg p(HGNC:MAPT,pmod(P))'
        result = self.parser.relation.parseString(statement)

        expected_dict = {
            SUBJECT: {
                MODIFIER: ACTIVITY,
                EFFECT: {
                    NAME: 'kin',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                },
                TARGET: {
                    FUNCTION: PROTEIN,
                    NAMESPACE: 'SFAM',
                    NAME: 'GSK3 Family'
                }
            },
            RELATION: NEGATIVE_CORRELATION,
            OBJECT: {
                FUNCTION: PROTEIN,
                NAMESPACE: 'HGNC',
                NAME: 'MAPT',
                VARIANTS: [pmod('Ph')]
            }
        }
        self.assertEqual(expected_dict, result.asDict())

        sub = protein('SFAM', 'GSK3 Family')
        self.assertHasNode(sub)

        obj = protein('HGNC', 'MAPT', variants=pmod('Ph'))
        self.assertHasNode(obj)

        edge_data = {
            RELATION: NEGATIVE_CORRELATION,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text,
            SUBJECT: {
                MODIFIER: ACTIVITY,
                EFFECT: {
                    NAME: 'kin',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                },
            },
        }
        self.assertHasEdge(sub, obj, **edge_data)

        edge_data_reverse = {
            RELATION: NEGATIVE_CORRELATION,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text,
            OBJECT: {
                MODIFIER: ACTIVITY,
                EFFECT: {
                    NAME: 'kin',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                },
            },
        }
        self.assertHasEdge(obj, sub, **edge_data_reverse)

    def test_positiveCorrelation_withSelfReferential(self):
        """
        3.2.2 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#XposCor
        Self-referential relationships"""
        statement = 'p(HGNC:GSK3B, pmod(P, S, 9)) pos act(p(HGNC:GSK3B), ma(kin))'
        result = self.parser.relation.parseString(statement)

        expected_dict = {
            SUBJECT: {
                FUNCTION: PROTEIN,
                NAMESPACE: 'HGNC',
                NAME: 'GSK3B',
                VARIANTS: [pmod('Ph', position=9, code='Ser')]
            },
            RELATION: POSITIVE_CORRELATION,
            OBJECT: {
                MODIFIER: ACTIVITY,
                TARGET: {
                    FUNCTION: PROTEIN,
                    NAMESPACE: 'HGNC',
                    NAME: 'GSK3B'
                },
                EFFECT: {
                    NAME: 'kin',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                }
            },
        }
        self.assertEqual(expected_dict, result.asDict())

        subject_node = protein('HGNC', 'GSK3B', variants=pmod('Ph', code='Ser', position=9))
        self.assertHasNode(subject_node)

        object_node = protein('HGNC', 'GSK3B')
        self.assertHasNode(object_node)

        edge_data = {
            RELATION: POSITIVE_CORRELATION,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text,
            OBJECT: {
                MODIFIER: ACTIVITY,
                EFFECT: {
                    NAME: 'kin',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                }
            },
        }
        self.assertHasEdge(subject_node, object_node, **edge_data)

        edge_data_reverse = {
            RELATION: POSITIVE_CORRELATION,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text,
            SUBJECT: {
                MODIFIER: ACTIVITY,
                EFFECT: {
                    NAME: 'kin',
                    NAMESPACE: BEL_DEFAULT_NAMESPACE
                }
            },
        }
        self.assertHasEdge(object_node, subject_node, **edge_data_reverse)

    def test_orthologous(self):
        """
        3.3.1 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#_orthologous
        """
        statement = 'g(HGNC:AKT1) orthologous g(MGI:AKT1)'
        result = self.parser.relation.parseString(statement)
        expected_result = [[GENE, 'HGNC', 'AKT1'], ORTHOLOGOUS, [GENE, 'MGI', 'AKT1']]
        self.assertEqual(expected_result, result.asList())

        sub = gene('HGNC', 'AKT1')
        self.assertHasNode(sub)

        obj = gene('MGI', 'AKT1')
        self.assertHasNode(obj)

        self.assertHasEdge(sub, obj, relation=ORTHOLOGOUS)
        self.assertHasEdge(obj, sub, relation=ORTHOLOGOUS)

    def test_transcription(self):
        """
        3.3.2 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#_transcribedto
        """
        statement = 'g(HGNC:AKT1) :> r(HGNC:AKT1)'
        result = self.parser.relation.parseString(statement)

        expected_result = [[GENE, 'HGNC', 'AKT1'], TRANSCRIBED_TO, [RNA, 'HGNC', 'AKT1']]
        self.assertEqual(expected_result, result.asList())

        sub = gene('HGNC', 'AKT1')
        self.assertHasNode(sub)

        obj = rna('HGNC', 'AKT1')
        self.assertHasNode(obj)

        self.assertHasEdge(sub, obj, relation=TRANSCRIBED_TO)

    def test_translation(self):
        """
        3.3.3 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#_translatedto
        """
        statement = 'r(HGNC:AKT1,loc(GOCC:intracellular)) >> p(HGNC:AKT1)'
        result = self.parser.relation.parseString(statement)

        # [[RNA, ['HGNC', 'AKT1']], TRANSLATED_TO, [PROTEIN, ['HGNC', 'AKT1']]]
        expected_result = {
            SUBJECT: {
                FUNCTION: RNA,
                NAMESPACE: 'HGNC',
                NAME: 'AKT1',
                LOCATION: {
                    NAMESPACE: 'GOCC',
                    NAME: 'intracellular'
                }
            },
            RELATION: TRANSLATED_TO,
            OBJECT: {
                FUNCTION: PROTEIN,
                NAMESPACE: 'HGNC',
                NAME: 'AKT1',
            }
        }
        self.assertEqual(expected_result, result.asDict())

        self.assertNumberNodes(2)

        source = rna('HGNC', 'AKT1')
        self.assertHasNode(source)

        target = protein('HGNC', 'AKT1')
        self.assertHasNode(target)

        self.assertNumberEdges(1)
        self.assertTrue(self.graph.has_edge(source, target))

        key_data = self.parser.graph.edge[source][target]
        self.assertEqual(1, len(key_data))

        key = list(key_data)[0]
        data = key_data[key]

        self.assertIn(RELATION, data)
        self.assertEqual(TRANSLATED_TO, data[RELATION])

        calculated_source_data = self.graph.node[source]
        self.assertTrue(calculated_source_data)

        calculated_target_data = self.graph.node[target]
        self.assertTrue(calculated_target_data)

    def test_component_list(self):
        s = 'complex(SCOMP:"C1 Complex") hasComponents list(p(HGNC:C1QB), p(HGNC:C1S))'
        result = self.parser.relation.parseString(s)

        expected_result_list = [
            [COMPLEX, 'SCOMP', 'C1 Complex'],
            'hasComponents',
            [
                [PROTEIN, 'HGNC', 'C1QB'],
                [PROTEIN, 'HGNC', 'C1S']
            ]
        ]
        self.assertEqual(expected_result_list, result.asList())

        sub = named_complex_abundance('SCOMP', 'C1 Complex')
        self.assertHasNode(sub)
        child_1 = protein('HGNC', 'C1QB')
        self.assertHasNode(child_1)
        self.assertHasEdge(sub, child_1, relation=HAS_COMPONENT)
        child_2 = protein('HGNC', 'C1S')
        self.assertHasNode(child_2)
        self.assertHasEdge(sub, child_2, relation=HAS_COMPONENT)

    def test_member_list(self):
        """
        3.4.2 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#_hasmembers
        """
        statement = 'p(PKC:a) hasMembers list(p(HGNC:PRKCA), p(HGNC:PRKCB), p(HGNC:PRKCD), p(HGNC:PRKCE))'
        result = self.parser.relation.parseString(statement)
        expected_result = [
            [PROTEIN, 'PKC', 'a'],
            'hasMembers',
            [
                [PROTEIN, 'HGNC', 'PRKCA'],
                [PROTEIN, 'HGNC', 'PRKCB'],
                [PROTEIN, 'HGNC', 'PRKCD'],
                [PROTEIN, 'HGNC', 'PRKCE']
            ]
        ]
        self.assertEqual(expected_result, result.asList())

        sub = protein('PKC', 'a')
        obj1 = protein('HGNC', 'PRKCA')
        obj2 = protein('HGNC', 'PRKCB')
        obj3 = protein('HGNC', 'PRKCD')
        obj4 = protein('HGNC', 'PRKCE')

        self.assertHasNode(sub)

        self.assertHasNode(obj1)
        self.assertHasEdge(sub, obj1, relation=HAS_MEMBER)

        self.assertHasNode(obj2)
        self.assertHasEdge(sub, obj2, relation=HAS_MEMBER)

        self.assertHasNode(obj3)
        self.assertHasEdge(sub, obj3, relation=HAS_MEMBER)

        self.assertHasNode(obj4)
        self.assertHasEdge(sub, obj4, relation=HAS_MEMBER)

    def test_is_a(self):
        """3.4.5 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#_isa"""
        statement = 'pathology(MESH:Psoriasis) isA pathology(MESH:"Skin Diseases")'
        result = self.parser.relation.parseString(statement)

        expected_result = [[PATHOLOGY, 'MESH', 'Psoriasis'], 'isA', [PATHOLOGY, 'MESH', 'Skin Diseases']]
        self.assertEqual(expected_result, result.asList())

        self.assertIn(psoriasis, self.graph)
        self.assertIn(skin_diseases, self.graph)

        self.assertIn(skin_diseases, self.graph[psoriasis])

        keys = list(self.graph[psoriasis][skin_diseases])
        self.assertEqual(1, len(keys))
        key = keys[0]
        self.assertIn(RELATION, self.graph[psoriasis][skin_diseases][key])
        self.assertEqual(IS_A, self.graph[psoriasis][skin_diseases][key][RELATION])

    def test_label_1(self):
        """Tests that a gene with variants can be relabeled. Note that the variant order has to be canonicalized"""
        statement = 'g(HGNC:APOE, var(c.526C>T), var(c.388T>C)) labeled "APOE E2"'
        result = self.parser.relation.parseString(statement)

        expected_dict = {
            SUBJECT: {
                FUNCTION: GENE,
                NAMESPACE: 'HGNC',
                NAME: 'APOE',
                VARIANTS: [
                    {
                        KIND: HGVS,
                        IDENTIFIER: 'c.526C>T'
                    }, {
                        KIND: HGVS,
                        IDENTIFIER: 'c.388T>C'
                    }
                ]
            },
            OBJECT: 'APOE E2'
        }
        self.assertEqual(expected_dict, result.asDict())

        self.assertNumberNodes(2)
        self.assertNumberEdges(1)

        # Variants will get canonicalized
        node = gene('HGNC', 'APOE', variants=[gene_substitution('C', 526, 'T'), gene_substitution('T', 388, 'C')])

        self.assertEqual('g(HGNC:APOE, var(c.388T>C), var(c.526C>T))', node.as_bel())
        self.assertHasNode(node)
        self.help_test_parent_in_graph(node)

        self.assertTrue(self.parser.graph.has_node_description(node))
        self.assertEqual('APOE E2', self.parser.graph.get_node_description(node))

    def test_raise_on_relabel(self):
        s1 = 'g(HGNC:APOE, var(c.526C>T), var(c.388T>C)) labeled "APOE E2"'
        s2 = 'g(HGNC:APOE, var(c.526C>T), var(c.388T>C)) labeled "APOE E2 Variant"'
        self.parser.relation.parseString(s1)
        with self.assertRaises(RelabelWarning):
            self.parser.relation.parseString(s2)

    def test_equivalentTo(self):
        statement = 'g(dbSNP:"rs123456") eq g(HGNC:YFG, var(c.123G>A))'
        result = self.parser.relation.parseString(statement)

        expected_result = {
            SUBJECT: {
                FUNCTION: GENE,
                NAMESPACE: 'dbSNP',
                NAME: 'rs123456',
            },
            RELATION: EQUIVALENT_TO,
            OBJECT: {
                FUNCTION: GENE,
                NAMESPACE: 'HGNC',
                NAME: 'YFG',
                VARIANTS: [
                    {
                        KIND: HGVS,
                        IDENTIFIER: 'c.123G>A'
                    }
                ]
            }
        }
        self.assertEqual(expected_result, result.asDict())

        sub = gene('dbSNP', 'rs123456')
        self.assertHasNode(sub)

        obj = gene('HGNC', 'YFG', variants=hgvs('c.123G>A'))
        self.assertHasNode(obj)

        self.assertTrue(self.graph.has_unqualified_edge(sub, obj, EQUIVALENT_TO))
        self.assertTrue(self.graph.has_unqualified_edge(obj, sub, EQUIVALENT_TO))

        self.assertHasEdge(sub, obj, **{RELATION: EQUIVALENT_TO})
        self.assertHasEdge(obj, sub, **{RELATION: EQUIVALENT_TO})

    def test_partOf(self):
        statement = 'a(UBERON:"corpus striatum") partOf a(UBERON:"basal ganglion")'
        self.parser.relation.parseString(statement)

        corpus_striatum = abundance(namespace='UBERON', name='corpus striatum')
        basal_ganglion = abundance(namespace='UBERON', name='basal ganglion')

        self.assertTrue(self.parser.graph.has_node(corpus_striatum))
        self.assertTrue(self.parser.graph.has_node(basal_ganglion))

        self.assertHasNode(corpus_striatum)
        self.assertHasNode(basal_ganglion)

        v = list(self.parser.graph.edge[corpus_striatum][basal_ganglion].values())
        self.assertEqual(1, len(v))

        v = v[0]
        self.assertIn(RELATION, v)
        self.assertEqual(PART_OF, v[RELATION])

    def test_subProcessOf(self):
        """
        3.4.6 http://openbel.org/language/web/version_2.0/bel_specification_version_2.0.html#_subprocessof
        """
        statement = 'rxn(reactants(a(CHEBI:"(S)-3-hydroxy-3-methylglutaryl-CoA"),a(CHEBI:NADPH), \
            a(CHEBI:hydron)),products(a(CHEBI:mevalonate), a(CHEBI:"CoA-SH"), a(CHEBI:"NADP(+)"))) \
            subProcessOf bp(GOBP:"cholesterol biosynthetic process")'
        result = self.parser.relation.parseString(statement)
        expected_result = [
            [
                REACTION,
                [
                    [ABUNDANCE, 'CHEBI', '(S)-3-hydroxy-3-methylglutaryl-CoA'],
                    [ABUNDANCE, 'CHEBI', 'NADPH'],
                    [ABUNDANCE, 'CHEBI', 'hydron'],
                ],
                [
                    [ABUNDANCE, 'CHEBI', 'mevalonate'],
                    [ABUNDANCE, 'CHEBI', 'CoA-SH'],
                    [ABUNDANCE, 'CHEBI', 'NADP(+)']
                ]
            ],
            SUBPROCESS_OF,
            [BIOPROCESS, 'GOBP', 'cholesterol biosynthetic process']]
        self.assertEqual(expected_result, result.asList())

        self.assertNumberNodes(8)
        self.assertNumberEdges(7)

        sub_reactant_1 = abundance('CHEBI', '(S)-3-hydroxy-3-methylglutaryl-CoA')
        sub_reactant_2 = abundance('CHEBI', 'NADPH')
        sub_reactant_3 = abundance('CHEBI', 'hydron')
        sub_product_1 = abundance('CHEBI', 'mevalonate')
        sub_product_2 = abundance('CHEBI', 'CoA-SH')
        sub_product_3 = abundance('CHEBI', 'NADP(+)')
        sub = reaction(reactants=[sub_reactant_1, sub_reactant_2, sub_reactant_3],
                       products=[sub_product_1, sub_product_2, sub_product_3])

        self.assertHasNode(sub)
        self.assertHasNode(sub_reactant_1)
        self.assertHasNode(sub_reactant_2)
        self.assertHasNode(sub_reactant_3)
        self.assertHasNode(sub_product_1)
        self.assertHasNode(sub_product_2)
        self.assertHasNode(sub_product_3)

        self.assertHasEdge(sub, sub_reactant_1, relation=HAS_REACTANT)
        self.assertHasEdge(sub, sub_reactant_2, relation=HAS_REACTANT)
        self.assertHasEdge(sub, sub_reactant_3, relation=HAS_REACTANT)
        self.assertHasEdge(sub, sub_product_1, relation=HAS_PRODUCT)
        self.assertHasEdge(sub, sub_product_2, relation=HAS_PRODUCT)
        self.assertHasEdge(sub, sub_product_3, relation=HAS_PRODUCT)

        cholesterol_biosynthetic_process = bioprocess('GOBP', 'cholesterol biosynthetic process')
        self.assertHasNode(cholesterol_biosynthetic_process)

        edge_data = {
            RELATION: SUBPROCESS_OF,
            CITATION: test_citation_dict,
            EVIDENCE: test_evidence_text,
        }
        self.assertHasEdge(sub, cholesterol_biosynthetic_process, **edge_data)

    def test_extra_1(self):
        statement = 'abundance(CHEBI:"nitric oxide") increases cellSurfaceExpression(complexAbundance(proteinAbundance(HGNC:ITGAV),proteinAbundance(HGNC:ITGB3)))'
        self.parser.relation.parseString(statement)

    def test_has_variant(self):
        statement = 'g(HGNC:AKT1) hasVariant g(HGNC:AKT1, gmod(M))'
        result = self.parser.relation.parseString(statement)

        expected_result = {
            SUBJECT: {
                FUNCTION: GENE,
                NAMESPACE: 'HGNC',
                NAME: 'AKT1',
            },
            RELATION: HAS_VARIANT,
            OBJECT: {
                FUNCTION: GENE,
                NAMESPACE: 'HGNC',
                NAME: 'AKT1',
                VARIANTS: [
                    {
                        KIND: GMOD,
                        IDENTIFIER: {
                            NAME: 'Me',
                            NAMESPACE: BEL_DEFAULT_NAMESPACE,
                        },
                    }
                ]
            }
        }
        self.assertEqual(expected_result, result.asDict(), msg='relation parsed wrong')

        sub = gene('HGNC', 'AKT1')
        obj = gene('HGNC', 'AKT1', variants=gmod('Me'))

        self.assertNumberNodes(2)
        self.assertNumberEdges(1)

        self.assertHasNode(sub)
        self.assertHasNode(obj)
        self.assertHasEdge(sub, obj, relation=HAS_VARIANT)

    def test_has_reaction_component(self):
        statement = 'rxn(reactants(a(CHEBI:"(S)-3-hydroxy-3-methylglutaryl-CoA"),a(CHEBI:NADPH), \
                    a(CHEBI:hydron)),products(a(CHEBI:mevalonate), a(CHEBI:"CoA-SH"), a(CHEBI:"NADP(+)"))) \
                    hasReactant a(CHEBI:"(S)-3-hydroxy-3-methylglutaryl-CoA")'
        result = self.parser.relation.parseString(statement)

        sub_reactant_1 = abundance('CHEBI', '(S)-3-hydroxy-3-methylglutaryl-CoA')
        sub_reactant_2 = abundance('CHEBI', 'NADPH')
        sub_reactant_3 = abundance('CHEBI', 'hydron')
        sub_product_1 = abundance('CHEBI', 'mevalonate')
        sub_product_2 = abundance('CHEBI', 'CoA-SH')
        sub_product_3 = abundance('CHEBI', 'NADP(+)')
        sub = reaction(
            reactants=[sub_reactant_1, sub_reactant_2, sub_reactant_3],
            products=[sub_product_1, sub_product_2, sub_product_3]
        )

        self.assertNumberNodes(7)
        self.assertNumberEdges(6)

        self.assertHasNode(sub)
        self.assertHasNode(sub_reactant_1)
        self.assertHasNode(sub_reactant_2)
        self.assertHasNode(sub_reactant_3)
        self.assertHasNode(sub_product_1)
        self.assertHasNode(sub_product_2)
        self.assertHasNode(sub_product_3)

        self.assertHasEdge(sub, sub_reactant_1, relation=HAS_REACTANT)
        self.assertHasEdge(sub, sub_reactant_2, relation=HAS_REACTANT)
        self.assertHasEdge(sub, sub_reactant_3, relation=HAS_REACTANT)
        self.assertHasEdge(sub, sub_product_1, relation=HAS_PRODUCT)
        self.assertHasEdge(sub, sub_product_2, relation=HAS_PRODUCT)
        self.assertHasEdge(sub, sub_product_3, relation=HAS_PRODUCT)


class TestCustom(unittest.TestCase):
    def setUp(self):
        graph = BELGraph()

        namespace_dict = {
            'HGNC': {
                'AKT1': 'GRP',
                'YFG': 'GRP'
            },
            'MESHCS': {
                'nucleus': 'A'
            }
        }

        self.parser = BelParser(graph, namespace_dict=namespace_dict, autostreamline=False)

    def test_tloc_undefined_namespace(self):
        s = 'tloc(p(HGNC:AKT1), fromLoc(MESHCS:nucleus), toLoc(MISSING:"undefined"))'

        with self.assertRaises(UndefinedNamespaceWarning):
            self.parser.translocation.parseString(s)

    def test_tloc_undefined_name(self):
        s = 'tloc(p(HGNC:AKT1), fromLoc(MESHCS:nucleus), toLoc(MESHCS:"undefined"))'

        with self.assertRaises(MissingNamespaceNameWarning):
            self.parser.translocation.parseString(s)

    def test_location_undefined_namespace(self):
        s = 'p(HGNC:AKT1, loc(MISSING:"nucleus")'

        with self.assertRaises(UndefinedNamespaceWarning):
            self.parser.protein.parseString(s)

    def test_location_undefined_name(self):
        s = 'p(HGNC:AKT1, loc(MESHCS:"undefined")'

        with self.assertRaises(MissingNamespaceNameWarning):
            self.parser.protein.parseString(s)


class TestWrite(TestTokenParserBase):
    def test_1(self):
        cases = [
            ('abundance(CHEBI:"superoxide")', 'a(CHEBI:superoxide)'),
            ('g(HGNC:AKT1,var(p.Phe508del))', 'g(HGNC:AKT1, var(p.Phe508del))'),
            ('geneAbundance(HGNC:AKT1, variant(p.Phe508del), sub(G,308,A), var(c.1521_1523delCTT))',
             'g(HGNC:AKT1, var(c.1521_1523delCTT), var(c.308G>A), var(p.Phe508del))'),
            ('p(HGNC:MAPT,proteinModification(P))', 'p(HGNC:MAPT, pmod(Ph))'),
            ('proteinAbundance(HGNC:SFN)', 'p(HGNC:SFN)'),
            ('complex(proteinAbundance(HGNC:SFN), p(HGNC:YWHAB))', 'complex(p(HGNC:SFN), p(HGNC:YWHAB))'),
            ('composite(proteinAbundance(HGNC:SFN), p(HGNC:YWHAB))', 'composite(p(HGNC:SFN), p(HGNC:YWHAB))'),
            ('reaction(reactants(a(CHEBI:superoxide)),products(a(CHEBI:"oxygen"),a(CHEBI:"hydrogen peroxide")))',
             'rxn(reactants(a(CHEBI:superoxide)), products(a(CHEBI:"hydrogen peroxide"), a(CHEBI:oxygen)))'),
            ('rxn(reactants(a(CHEBI:superoxide)),products(a(CHEBI:"hydrogen peroxide"), a(CHEBI:"oxygen")))',
             'rxn(reactants(a(CHEBI:superoxide)), products(a(CHEBI:"hydrogen peroxide"), a(CHEBI:oxygen)))'),
            ('g(HGNC:AKT1, geneModification(M))', 'g(HGNC:AKT1, gmod(Me))'),
            'g(fus(HGNC:TMPRSS2, p.1_79, HGNC:ERG, p.312_5034))',
            'g(fus(HGNC:TMPRSS2, r.1_?, HGNC:ERG, r.312_5034))',
            'g(fus(HGNC:TMPRSS2, r.1_79, HGNC:ERG, r.?_5034))',
            ('g(HGNC:CHCHD4, fusion(HGNC:AIFM1))', 'g(fus(HGNC:CHCHD4, ?, HGNC:AIFM1, ?))'),
            ('g(HGNC:CHCHD4, fusion(HGNC:AIFM1, ?, ?))', 'g(fus(HGNC:CHCHD4, ?, HGNC:AIFM1, ?))'),
            'g(fus(HGNC:TMPRSS2, ?, HGNC:ERG, ?))',
        ]

        self.parser.bel_term.addParseAction(self.parser.handle_term)

        for case in cases:
            source_bel, expected_bel = case if 2 == len(case) else (case, case)

            result = self.parser.bel_term.parseString(source_bel)
            tokens = result.asDict()
            entity = dict_to_entity(tokens)
            bel = entity.as_bel()
            self.assertEqual(expected_bel, bel)
