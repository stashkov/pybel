# -*- coding: utf-8 -*-

import logging
from collections import defaultdict

from .utils import update_metadata
from ..graph import BELGraph
from ..pipeline import splitter
from ...constants import ANNOTATIONS

log = logging.getLogger(__name__)

__all__ = [
    'get_subgraphs_by_annotation',
    'get_subgraphs_by_annotation_filtered',
]


def _get_subgraphs_by_annotation_disregard_undefined(graph, annotation):
    result = defaultdict(BELGraph)

    for source, target, key, data in graph.edges_iter(keys=True, data=True):
        annotation_dict = data.get(ANNOTATIONS)

        if annotation_dict is None or annotation not in annotation_dict:
            continue

        for value in annotation_dict[annotation]:
            result[value].add_edge(source, target, key=key, **data)

    return dict(result)


def _get_subgraphs_by_annotation_keep_undefined(graph, annotation, sentinel):
    result = defaultdict(BELGraph)

    for source, target, key, data in graph.edges_iter(keys=True, data=True):
        annotation_dict = data.get(ANNOTATIONS)

        if annotation_dict is None or annotation not in annotation_dict:
            result[sentinel].add_edge(source, target, key=key, **data)
        else:
            for value in annotation_dict[annotation]:
                result[value].add_edge(source, target, key=key, **data)

    return dict(result)


@splitter
def get_subgraphs_by_annotation(graph, annotation, keep_undefined=True, sentinel='Undefined'):
    """Stratifies the given graph into subgraphs based on the values for edges' annotations

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation: The annotation to group by
    :param bool keep_undefined: If true, uses the sentinel value to store a subgraph of edges not matching the given
     annotation.
    :param str sentinel: The value to stick unannotated edges into
    :rtype: dict[str,pybel.BELGraph]
    """
    if keep_undefined:
        rv = _get_subgraphs_by_annotation_keep_undefined(graph, annotation, sentinel)
    else:
        rv = _get_subgraphs_by_annotation_disregard_undefined(graph, annotation)

    for value in rv.values():
        update_metadata(value, graph)

    return rv


@splitter
def get_subgraphs_by_annotation_filtered(graph, annotation, values):
    """Stratifies the given graph into subgraphs based on the values for edges' annotations, but filter by a set
    of given values

    :param pybel.BELGraph graph: A BEL graph
    :param str annotation: The annotation to group by
    :param iter[str] values: The values to keep
    :rtype: dict[str,pybel.BELGraph]
    """
    result = defaultdict(BELGraph)
    values = set(values)

    for source, target, key, data in graph.edges_iter(keys=True, data=True):
        annotation_dict = data.get(ANNOTATIONS)

        if annotation_dict is None or annotation not in annotation_dict:
            continue

        for value in annotation_dict[annotation]:
            if value in values:
                result[value].add_edge(source, target, key=key, **data)

    for value in result.values():
        update_metadata(value, graph)

    return dict(result)
