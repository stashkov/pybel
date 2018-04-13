# -*- coding: utf-8 -*-

class QueryMissingNetworksError(KeyError):
    """Raised if a query is created from json but doesn't have a listing of network identifiers"""


class MissingPipelineFunctionError(KeyError):
    """Raised when trying to run the pipeline with a function that isn't registered"""
