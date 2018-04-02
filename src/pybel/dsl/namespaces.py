# -*- coding: utf-8 -*-

"""This module contains simple wrappers around node DSL functions for common namespaces"""

from .nodes import abundance, protein


class chebi(abundance):
    def __init__(self, name=None, identifier=None):
        super(chebi, self).__init__(namespace='CHEBI', name=name, identifier=identifier)


class hgnc(protein):
    def __init__(self, name=None, identifier=None):
        super(hgnc, self).__init__(namespace='HGNC', name=name, identifier=identifier)
