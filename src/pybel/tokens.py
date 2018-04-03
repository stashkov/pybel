# -*- coding: utf-8 -*-

"""This module helps handle node data dictionaries"""

from .constants import *
from .dsl import cell_surface_expression, secretion
from .dsl.nodes import (
    BaseAbundance, Variant, abundance, bioprocess, complex_abundance, composite_abundance, fragment, fusion_range, gene,
    gene_fusion, gmod, hgvs, mirna, missing_fusion_range, named_complex_abundance, pathology, pmod, protein,
    protein_fusion, reaction, rna, rna_fusion,
)

__all__ = [
    'sort_abundances',
]


def safe_get_dict(tokens):
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

    else:
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

    partner_5p = dsl_func(
        namespace=tokens[FUSION][PARTNER_5P][NAMESPACE],
        name=tokens[FUSION][PARTNER_5P][NAME]
    )

    partner_3p = dsl_func(
        namespace=tokens[FUSION][PARTNER_3P][NAMESPACE],
        name=tokens[FUSION][PARTNER_3P][NAME]
    )

    range_5p = fusion_range_po_to_dict(tokens[FUSION][RANGE_5P])
    range_3p = fusion_range_po_to_dict(tokens[FUSION][RANGE_3P])

    return dsl_fusion_func(
        partner_5p=partner_5p,
        range_5p=range_5p,
        partner_3p=partner_3p,
        range_3p=range_3p
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


def _as_tuple(v):
    return v.as_tuple()


def sort_variants(variants):
    """
    :param list[Variant] variants: A list of variants
    :rtype: list[Variant]
    """
    return sorted([variant for variant in variants if variant], key=_as_tuple)


def variant_po_to_dict_helper(tokens):
    """Converts a PyParsing data dictionary to a PyBEL variant data dictionary

    :type tokens: ParseResult
    :rtype: list[pybel.dsl.nodes.Variant]
    """
    variants = [
        _variant_dict_to_dsl(safe_get_dict(variant))
        for variant in tokens[VARIANTS]
        # TODO check variant[KIND] ?
    ]

    return sort_variants(variants)


def variant_po_to_dict(tokens):
    """Converts a PyParsing data dictionary to a PyBEL variant data dictionary

    :type tokens: ParseResult
    :rtype: pybel.dsl.nodes.BaseAbundance
    """
    attr_data = simple_po_to_dict(tokens)
    attr_data[VARIANTS] = variant_po_to_dict_helper(tokens)
    return attr_data


def sort_abundances(tokens):
    """Sorts a list of PyBEL data dictionaries to their canonical ordering

    :param iter[BaseAbundance] tokens:
    :rtype: list[BaseAbundance]
    """
    return sorted(tokens, key=_as_tuple)


def reaction_part_po_to_dict(tokens):
    """
    :type tokens: ParseResult
    :rtype: list[BaseAbundance]
    """
    return sort_abundances(
        po_to_dict(token)
        for token in tokens
    )


def reaction_po_to_dict(tokens):
    """
    :type tokens: ParseResult
    :rtype: pybel.dsl.nodes.reaction
    """
    return reaction(
        reactants=reaction_part_po_to_dict(tokens[REACTANTS]),
        products=reaction_part_po_to_dict(tokens[PRODUCTS]),
    )


def simple_po_to_dict(tokens):
    """
    :type tokens: ParseResult
    :rtype: BaseAbundance
    """
    return _func_to_dsl[tokens[FUNCTION]](
        namespace=tokens[NAMESPACE],
        name=tokens[NAME]
    )


def fragment_po_to_dsl(tokens):
    if FRAGMENT_MISSING in tokens:
        return fragment()

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
    list_entries = [
        po_to_dict(token)
        for token in tokens[MEMBERS]
    ]

    return _list_function_to_dsl[tokens[FUNCTION]](
        members=sort_abundances(list_entries)
    )


def po_to_dict(tokens):
    """Converts a dictionary to a BaseEntity

    :type tokens: ParseResult or dict
    :rtype: pybel.dsl.nodes.BaseEntity
    """
    if MODIFIER in tokens:
        return po_to_dict(tokens[TARGET])

    elif REACTION == tokens[FUNCTION]:
        return reaction_po_to_dict(tokens)

    elif VARIANTS in tokens:
        return variant_po_to_dict(tokens)

    elif MEMBERS in tokens:
        return list_po_to_dict(tokens)

    elif FUSION in tokens:
        return fusion_po_to_dict(tokens)

    return simple_po_to_dict(tokens)


def modifier_po_to_dict(tokens):
    """Get activity, transformation, or transformation information as a dictionary

    :return: a dictionary describing the modifier
    :rtype: dict
    """
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
