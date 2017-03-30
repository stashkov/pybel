# -*- coding: utf-8 -*-

"""This module contains IO functions for BEL scripts"""

import codecs
import logging
import os

import requests
from requests_file import FileAdapter

from ..graph import BELGraph, parse_lines

__all__ = [
    'from_lines',
    'from_path',
    'from_url'
]

log = logging.getLogger(__name__)


def from_lines(lines, manager=None, allow_naked_names=False, allow_nested=False, citation_clearing=True, **kwargs):
    """Loads a BEL graph from an iterable over the lines of a BEL script. This can be a list of strings, file, or other.
    This function is a *very* thin wrapper around :class:`BELGraph`.

    :param lines: An iterable of strings (the lines in a BEL script)
    :type lines: iter
    :param manager: database connection string to cache, pre-built CacheManager, pre-built MetadataParser
                        or None to use default cache
    :type manager: str or :class:`pybel.manager.CacheManager` or :class:`pybel.parser.MetadataParser`
    :param allow_naked_names: if true, turn off naked namespace failures
    :type allow_naked_names: bool
    :param allow_nested: if true, turn off nested statement failures
    :type allow_nested: bool
    :param citation_clearing: Should :code:`SET Citation` statements clear evidence and all annotations?
                                Delegated to :class:`pybel.parser.ControlParser`
    :type citation_clearing: bool
    :param kwargs: keyword arguments to pass to :class:`networkx.MultiDiGraph`
    :type kwargs: dict
    :return: a parsed BEL graph
    :rtype: BELGraph
    """
    graph = BELGraph(**kwargs)
    parse_lines(graph, lines, manager=manager, allow_naked_names=allow_naked_names, allow_nested=allow_nested,
                citation_clearing=citation_clearing)
    return graph


def from_path(path, manager=None, allow_naked_names=False, allow_nested=False, citation_clearing=True,
              encoding='utf-8', **kwargs):
    """Loads a BEL graph from a file resource

    :param path: A file path
    :type path: str
    :param manager: database connection string to cache, pre-built CacheManager, pre-built MetadataParser
                        or None to use default cache
    :type manager: str or :class:`pybel.manager.CacheManager` or :class:`pybel.parser.MetadataParser`
    :param allow_naked_names: if true, turn off naked namespace failures
    :type allow_naked_names: bool
    :param allow_nested: if true, turn off nested statement failures
    :type allow_nested: bool
    :param citation_clearing: Should :code:`SET Citation` statements clear evidence and all annotations?
                                Delegated to :class:`pybel.parser.ControlParser`
    :type citation_clearing: bool
    :param encoding: the encoding to use when reading this file. Is passed to :code:`codecs.open`.
                     See the python `docs <https://docs.python.org/3/library/codecs.html#standard-encodings>`_ for a
                     list of standard encodings. For example, files starting with a UTF-8 BOM should use
                     :code:`utf_8_sig`
    :type encoding: str
    :param kwargs: Keyword arguments to pass to :class:`networkx.MultiDiGraph`
    :type kwargs: dict
    :return: A parsed BEL graph
    :rtype: BELGraph
    """
    log.info('Loading from path: %s', path)
    with codecs.open(os.path.expanduser(path), encoding=encoding) as file:
        return from_lines(file, manager=manager, allow_naked_names=allow_naked_names, allow_nested=allow_nested,
                          citation_clearing=citation_clearing, **kwargs)


def from_url(url, manager=None, allow_naked_names=False, allow_nested=False, citation_clearing=True, **kwargs):
    """Loads a BEL graph from a URL resource

    :param url: A valid URL pointing to a BEL resource
    :type url: str
    :param manager: database connection string to cache, pre-built CacheManager, pre-built MetadataParser
                        or None to use default cache
    :type manager: str or :class:`pybel.manager.CacheManager` or :class:`pybel.parser.MetadataParser`
    :param allow_naked_names: if true, turn off naked namespace failures
    :type allow_naked_names: bool
    :param allow_nested: if true, turn off nested statement failures
    :type allow_nested: bool
    :param citation_clearing: Should :code:`SET Citation` statements clear evidence and all annotations?
                                Delegated to :class:`pybel.parser.ControlParser`
    :type citation_clearing: bool
    :param kwargs: Keyword arguments to pass to :class:`networkx.MultiDiGraph`
    :type kwargs: dict
    :return: A parsed BEL graph
    :rtype: BELGraph
    """
    log.info('Loading from url: %s', url)

    session = requests.session()
    session.mount('file://', FileAdapter())

    response = session.get(url)
    response.raise_for_status()

    lines = (line.decode('utf-8') for line in response.iter_lines())

    return from_lines(lines, manager=manager, allow_naked_names=allow_naked_names, allow_nested=allow_nested,
                      citation_clearing=citation_clearing, **kwargs)