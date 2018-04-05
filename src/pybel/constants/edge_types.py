# -*- coding: utf-8 -*-

#: A BEL relationship
HAS_REACTANT = 'hasReactant'
#: A BEL relationship
HAS_PRODUCT = 'hasProduct'
#: A BEL relationship
HAS_COMPONENT = 'hasComponent'
#: A BEL relationship
HAS_VARIANT = 'hasVariant'
#: A BEL relationship
HAS_MEMBER = 'hasMember'
#: A BEL relationship
#: :data:`GENE` to :data:`RNA` is called transcription
TRANSCRIBED_TO = 'transcribedTo'
#: A BEL relationship
#: :data:`RNA` to :data:`PROTEIN` is called translation
TRANSLATED_TO = 'translatedTo'
#: A BEL relationship
INCREASES = 'increases'
#: A BEL relationship
DIRECTLY_INCREASES = 'directlyIncreases'
#: A BEL relationship
DECREASES = 'decreases'
#: A BEL relationship
DIRECTLY_DECREASES = 'directlyDecreases'
#: A BEL relationship
CAUSES_NO_CHANGE = 'causesNoChange'
#: A BEL relationship
REGULATES = 'regulates'
#: A BEL relationship
NEGATIVE_CORRELATION = 'negativeCorrelation'
#: A BEL relationship
POSITIVE_CORRELATION = 'positiveCorrelation'
#: A BEL relationship
ASSOCIATION = 'association'
#: A BEL relationship
ORTHOLOGOUS = 'orthologous'
#: A BEL relationship
ANALOGOUS_TO = 'analogousTo'
#: A BEL relationship
IS_A = 'isA'
#: A BEL relationship
RATE_LIMITING_STEP_OF = 'rateLimitingStepOf'
#: A BEL relationship
SUBPROCESS_OF = 'subProcessOf'
#: A BEL relationship
BIOMARKER_FOR = 'biomarkerFor'
#: A BEL relationship
PROGONSTIC_BIOMARKER_FOR = 'prognosticBiomarkerFor'
#: A BEL relationship, added by PyBEL
EQUIVALENT_TO = 'equivalentTo'
#: A BEL relationship, added by PyBEL
PART_OF = 'partOf'

#: A set of all causal relationships that have an increasing effect
CAUSAL_INCREASE_RELATIONS = {INCREASES, DIRECTLY_INCREASES}
#: A set of all causal relationships that have a decreasing effect
CAUSAL_DECREASE_RELATIONS = {DECREASES, DIRECTLY_DECREASES}
#: A set of direct causal relations
DIRECT_CAUSAL_RELATIONS = {DIRECTLY_DECREASES, DIRECTLY_INCREASES}
#: A set of direct causal relations
INDIRECT_CAUSAL_RELATIONS = {DECREASES, INCREASES}
#: A set of all causal relationships
CAUSAL_RELATIONS = {INCREASES, DIRECTLY_INCREASES, DECREASES, DIRECTLY_DECREASES}

#: A set of all relationships that are inherently directionless, and are therefore added to the graph twice
TWO_WAY_RELATIONS = {
    NEGATIVE_CORRELATION,
    POSITIVE_CORRELATION,
    ASSOCIATION,
    ORTHOLOGOUS,
    ANALOGOUS_TO,
    EQUIVALENT_TO,
}

#: A set of all correlative relationships
CORRELATIVE_RELATIONS = {
    POSITIVE_CORRELATION,
    NEGATIVE_CORRELATION
}

#: A set of polar relations
POLAR_RELATIONS = CAUSAL_RELATIONS | CORRELATIVE_RELATIONS

#: A list of relationship types that don't require annotations or evidence
#: This must be maintained as a list, since the :data:`unqualified_edge_code` is calculated based on the order
#: and needs to be consistent
unqualified_edges = [
    HAS_REACTANT,
    HAS_PRODUCT,
    HAS_COMPONENT,
    HAS_VARIANT,
    TRANSCRIBED_TO,
    TRANSLATED_TO,
    HAS_MEMBER,
    IS_A,
    EQUIVALENT_TO,
    PART_OF,
    ORTHOLOGOUS,
]

UNQUALIFIED_EDGES = set(unqualified_edges)

#: Unqualified edges are given negative keys since the standard NetworkX edge key factory starts at 0 and counts up
unqualified_edge_code = {relation: -i for i, relation in enumerate(unqualified_edges, start=1)}
