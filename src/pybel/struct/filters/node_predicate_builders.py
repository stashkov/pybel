# -*- coding: utf-8 -*-

from collections import Iterable

from ...constants import FUNCTION

__all__ = [
    'function_inclusion_filter_builder',
    'data_missing_key_builder',
]


def _single_function_inclusion_filter_builder(func):
    def function_inclusion_filter(graph, node):
        """Passes only for a node that has the enclosed function

        :param BELGraph graph: A BEL Graph
        :param tuple node: A BEL node
        :return: If the node doesn't have the enclosed function
        :rtype: bool
        """
        return graph.node[node][FUNCTION] == func

    return function_inclusion_filter


def _collection_function_inclusion_builder(funcs):
    funcs = set(funcs)

    def functions_inclusion_filter(graph, node):
        """Passes only for a node that is one of the enclosed functions

        :param BELGraph graph: A BEL Graph
        :param tuple node: A BEL node
        :return: If the node doesn't have the enclosed functions
        :rtype: bool
        """
        return graph.node[node][FUNCTION] in funcs

    return functions_inclusion_filter


def function_inclusion_filter_builder(func):
    """Builds a filter that only passes on nodes of the given function(s)

    :param func: A BEL Function or list/set/tuple of BEL functions
    :type func: str or iter[str]
    :return: A node filter (graph, node) -> bool
    :rtype: types.FunctionType
    """
    if isinstance(func, str):
        return _single_function_inclusion_filter_builder(func)

    elif isinstance(func, Iterable):
        return _collection_function_inclusion_builder(func)

    raise ValueError('Invalid type for argument: {}'.format(func))


def data_missing_key_builder(key):
    """Builds a filter that passes only on nodes that don't have the given key in their data dictionary

    :param str key: A key for the node's data dictionary
    :return: A node filter (graph, node) -> bool
    :rtype: types.FunctionType
    """

    def data_does_not_contain_key(graph, node):
        """Passes only for a node that doesn't contain the enclosed key in its data dictionary

        :param pybel.BELGraph graph: A BEL Graph
        :param BaseEntity node: A BEL node
        :return: If the node doesn't contain the enclosed key in its data dictionary
        :rtype: bool
        """
        return key not in graph.node[node]

    return data_does_not_contain_key
