# -*- coding: utf-8 -*-

import six

from ...pipeline import in_place_mutator
from ....constants import FUNCTION, MIRNA, NAMESPACE, PROTEIN, RNA, VARIANTS

__all__ = [
    'infer_central_dogmatic_translations_by_namespace',
    'infer_central_dogmatic_transcriptions',
    'infer_central_dogmatic_translations',
    'infer_central_dogma',
]


def _infer_converter_helper(node, data, new_function):
    new_tup = list(node)
    new_tup[0] = new_function
    new_tup = tuple(new_tup)
    new_dict = data.copy()
    new_dict[FUNCTION] = new_function
    return new_tup, new_dict


@in_place_mutator
def infer_central_dogmatic_translations_by_namespace(graph, namespaces):
    """For all Protein entities in the given namespaces, adds the missing origin RNA and RNA-Protein translation edge

    :param pybel.BELGraph graph: A BEL graph
    :param str or iter[str] namespaces: The namespaces over which to do this
    """
    namespaces = {namespaces} if isinstance(namespaces, six.string_types) else set(namespaces)

    for node in list(graph):
        if node[FUNCTION] != PROTEIN:
            continue

        if NAMESPACE not in node:
            continue

        if VARIANTS in node:
            continue

        if node[NAMESPACE] not in namespaces:
            continue

        graph.add_translation(node.get_rna(), node)


@in_place_mutator
def infer_central_dogmatic_translations(graph):
    """For all HGNC Protein entities, adds the missing origin RNA and RNA-Protein translation edge

    :param pybel.BELGraph graph: A BEL graph
    """
    infer_central_dogmatic_translations_by_namespace(graph, 'HGNC')


@in_place_mutator
def infer_central_dogmatic_transcriptions(graph):
    """For all RNA entities, adds the missing origin Gene and Gene-RNA transcription edge

    :param pybel.BELGraph graph: A BEL graph
    """
    for node in list(graph):
        if node[FUNCTION] not in {MIRNA, RNA}:
            continue

        if NAMESPACE not in node:
            continue

        if VARIANTS in node:
            continue

        graph.add_transcription(node.get_gene(), node)


@in_place_mutator
def infer_central_dogma(graph):
    """Adds all RNA-Protein translations then all Gene-RNA transcriptions by applying
    :func:`infer_central_dogmatic_translations` then :func:`infer_central_dogmatic_transcriptions`

    :param pybel.BELGraph graph: A BEL graph
    """
    infer_central_dogmatic_translations(graph)
    infer_central_dogmatic_transcriptions(graph)
