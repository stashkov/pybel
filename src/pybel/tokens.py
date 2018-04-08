# -*- coding: utf-8 -*-

"""This module helps handle node data dictionaries"""

from .constants import *
from .dsl import cell_surface_expression, secretion
from .dsl.nodes import (
    BaseAbundance, abundance, bioprocess, complex_abundance, composite_abundance, fragment, fusion_range, gene,
    gene_fusion, gmod, hgvs, mirna, missing_fusion_range, named_complex_abundance, pathology, pmod, protein,
    protein_fusion, reaction, rna, rna_fusion,
)

__all__ = [
    'dict_to_entity',
]


def safe_get_dict(tokens):  # FIXME get rid of this
    if hasattr(tokens, 'asDict'):
        return tokens.asDict()
    return dict(tokens)


def fusion_range_po_to_dict(tokens):
    """
    :type tokens: ParseResult
    :rtype: pybel.dsl.nodes.FusionRangeBase
    """
    if FUSION_MISSING in tokens:
        return missing_fusion_range()

    return fusion_range(
        reference=tokens[FUSION_REFERENCE],
        start=tokens[FUSION_START],
        stop=tokens[FUSION_STOP]
    )


_fusion_func_to_dsl = {
    GENE: gene_fusion,
    RNA: rna_fusion,
    PROTEIN: protein_fusion
}

_func_to_dsl = {
    GENE: gene,
    RNA: rna,
    PROTEIN: protein,
    MIRNA: mirna,
    ABUNDANCE: abundance,
    PATHOLOGY: pathology,
    BIOPROCESS: bioprocess,
    COMPLEX: named_complex_abundance
}


def fusion_po_to_dict(tokens):
    """Converts a PyParsing data dictionary to a PyBEL fusion data dictionary

    :param tokens: A PyParsing data dictionary representing a fusion
    :type tokens: ParseResult
    :rtype: pybel.dsl.nodes.FusionBase
    """
    tokens_func = tokens[FUNCTION]
    dsl_fusion_func = _fusion_func_to_dsl[tokens_func]
    dsl_func = _func_to_dsl[tokens_func]

    partner5p = dsl_func(
        namespace=tokens[FUSION][PARTNER_5P][NAMESPACE],
        name=tokens[FUSION][PARTNER_5P][NAME]
    )

    partner3p = dsl_func(
        namespace=tokens[FUSION][PARTNER_3P][NAMESPACE],
        name=tokens[FUSION][PARTNER_3P][NAME]
    )

    range5p = fusion_range_po_to_dict(tokens[FUSION][RANGE_5P])
    range3p = fusion_range_po_to_dict(tokens[FUSION][RANGE_3P])

    return dsl_fusion_func(
        partner5p=partner5p,
        range5p=range5p,
        partner3p=partner3p,
        range3p=range3p
    )


def _variant_dict_to_dsl(d):
    """

    :param dict d:
    :rtype: pybel.dsl.nodes.Variant
    """
    kind = d[KIND]

    if kind == FRAGMENT:
        return fragment_po_to_dsl(d)

    if kind == HGVS:
        return hgvs(d[IDENTIFIER])

    if kind == GMOD:
        return gmod(
            namespace=d[IDENTIFIER][NAMESPACE],
            name=d[IDENTIFIER][NAME],
        )

    if kind == PMOD:
        return pmod(
            namespace=d[IDENTIFIER][NAMESPACE],
            name=d[IDENTIFIER][NAME],
            position=d.get(PMOD_POSITION),
            code=d.get(PMOD_CODE)
        )

    raise ValueError


def variant_po_to_dict_helper(tokens):
    """Converts a PyParsing data dictionary to a PyBEL variant data dictionary

    :type tokens: ParseResult
    :rtype: list[pybel.dsl.nodes.Variant]
    """
    return [
        _variant_dict_to_dsl(safe_get_dict(variant))
        for variant in tokens[VARIANTS]
        # TODO check variant[KIND] ?
    ]


def variant_po_to_dict(tokens):
    """Converts a PyParsing data dictionary to a PyBEL variant data dictionary

    :type tokens: ParseResult
    :rtype: pybel.dsl.nodes.BaseAbundance
    """
    if tokens[FUNCTION] not in {PROTEIN, GENE, RNA, MIRNA}:
        raise TypeError

    parent = simple_dict_to_entity(tokens)
    variants = variant_po_to_dict_helper(tokens)
    return parent.with_variants(variants)


def reaction_dict_to_entity(tokens):
    """
    :type tokens: ParseResult
    :rtype: pybel.dsl.nodes.reaction
    """
    return reaction(
        reactants=[dict_to_entity(reactant) for reactant in tokens[REACTANTS]],
        products=[dict_to_entity(token) for token in tokens[PRODUCTS]],
    )


def simple_dict_to_entity(tokens):
    """
    :type tokens: ParseResult
    :rtype: BaseAbundance
    """
    return _func_to_dsl[tokens[FUNCTION]](
        namespace=tokens[NAMESPACE],
        name=tokens.get(NAME),
        identifier=tokens.get(IDENTIFIER),
    )


def fragment_po_to_dsl(tokens):
    """Converts the tokens for a fragment into the DSL object

    :param dict tokens:
    :rtype: fragment
    """
    if FRAGMENT_MISSING in tokens:
        return fragment(description=tokens.get(FRAGMENT_DESCRIPTION))

    return fragment(
        start=tokens[FRAGMENT_START],
        stop=tokens[FRAGMENT_STOP],
        description=tokens.get(FRAGMENT_DESCRIPTION)
    )


_list_function_to_dsl = {
    COMPOSITE: composite_abundance,
    COMPLEX: complex_abundance,
}


def list_po_to_dict(tokens):
    """
    :param tokens: PyParsing ParseObject
    :rtype: pybel.dsl.nodes.ListAbundance
    """
    members = [
        dict_to_entity(token)
        for token in tokens[MEMBERS]
    ]

    dsl_func = _list_function_to_dsl[tokens[FUNCTION]]

    return dsl_func(members=members)


def dict_to_entity(tokens):
    """Converts a dictionary to a BaseEntity

    :type tokens: ParseResult or dict
    :rtype: pybel.dsl.nodes.BaseEntity
    """
    if MODIFIER in tokens:
        return dict_to_entity(tokens[TARGET])

    elif REACTION == tokens[FUNCTION]:
        return reaction_dict_to_entity(tokens)

    elif VARIANTS in tokens:
        return variant_po_to_dict(tokens)

    elif MEMBERS in tokens:
        return list_po_to_dict(tokens)

    elif FUSION in tokens:
        return fusion_po_to_dict(tokens)

    return simple_dict_to_entity(tokens)


def modifier_po_to_dict(tokens):
    """Get activity, transformation, or transformation information as a dictionary

    :return: a dictionary describing the modifier
    :rtype: Optional[dict]
    """
    if LOCATION not in tokens and MODIFIER not in tokens:
        return  # rather return None than an empty dict

    attrs = {}

    if LOCATION in tokens:
        attrs[LOCATION] = tokens[LOCATION].asDict()

    if MODIFIER not in tokens:
        return attrs

    if LOCATION in tokens[TARGET]:
        attrs[LOCATION] = tokens[TARGET][LOCATION].asDict()

    if tokens[MODIFIER] == DEGRADATION:
        attrs[MODIFIER] = tokens[MODIFIER]

    elif tokens[MODIFIER] == ACTIVITY:
        attrs[MODIFIER] = tokens[MODIFIER]

        if EFFECT in tokens:
            attrs[EFFECT] = dict(tokens[EFFECT])

    elif tokens[MODIFIER] == TRANSLOCATION:
        attrs[MODIFIER] = tokens[MODIFIER]

        if EFFECT in tokens:
            attrs[EFFECT] = tokens[EFFECT].asDict()

    elif tokens[MODIFIER] == CELL_SECRETION:
        attrs.update(secretion())

    elif tokens[MODIFIER] == CELL_SURFACE_EXPRESSION:
        attrs.update(cell_surface_expression())

    else:
        raise ValueError('Invalid value for tokens[MODIFIER]: {}'.format(tokens[MODIFIER]))

    return attrs
