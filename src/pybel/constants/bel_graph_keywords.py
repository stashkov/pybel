# -*- coding: utf-8 -*-

# Internal metadata representation. See BELGraph documentation, since these are shielded from the user by properties.

#: The key for the document metadata dictionary. Can be accessed by :code:`graph.graph[GRAPH_METADATA]`, or by using
#: the property built in to the :class:`pybel.BELGraph`, :func:`pybel.BELGraph.document`
GRAPH_METADATA = 'document_metadata'
GRAPH_NAMESPACE_URL = 'namespace_url'
GRAPH_NAMESPACE_PATTERN = 'namespace_pattern'
GRAPH_ANNOTATION_URL = 'annotation_url'
GRAPH_ANNOTATION_PATTERN = 'annotation_pattern'
GRAPH_ANNOTATION_LIST = 'annotation_list'
GRAPH_WARNINGS = 'warnings'
GRAPH_PYBEL_VERSION = 'pybel_version'
GRAPH_UNCACHED_NAMESPACES = 'namespaces_uncached'

#: The key for the document name. Can be accessed by :code:`graph.document[METADATA_NAME]` or by using the property
#: built into the :class:`pybel.BELGraph` class, :func:`pybel.BELGraph.name`
METADATA_NAME = 'name'
#: The key for the document version. Can be accessed by :code:`graph.document[METADATA_VERSION]`
METADATA_VERSION = 'version'
#: The key for the document description. Can be accessed by :code:`graph.document[METADATA_DESCRIPTION]`
METADATA_DESCRIPTION = 'description'
#: The key for the document authors. Can be accessed by :code:`graph.document[METADATA_NAME]`
METADATA_AUTHORS = 'authors'
#: The key for the document contact email. Can be accessed by :code:`graph.document[METADATA_CONTACT]`
METADATA_CONTACT = 'contact'
#: The key for the document licenses. Can be accessed by :code:`graph.document[METADATA_LICENSES]`
METADATA_LICENSES = 'licenses'
#: The key for the document copyright information. Can be accessed by :code:`graph.document[METADATA_COPYRIGHT]`
METADATA_COPYRIGHT = 'copyright'
#: The key for the document disclaimer. Can be accessed by :code:`graph.document[METADATA_DISCLAIMER]`
METADATA_DISCLAIMER = 'disclaimer'
#: The key for the document project. Can be accessed by :code:`graph.document[METADATA_PROJECT]`
METADATA_PROJECT = 'project'

#: The keys to use when inserting a graph to the cache
METADATA_INSERT_KEYS = {
    METADATA_NAME,
    METADATA_VERSION,
    METADATA_DESCRIPTION,
    METADATA_AUTHORS,
    METADATA_CONTACT,
    METADATA_LICENSES,
    METADATA_COPYRIGHT,
    METADATA_DISCLAIMER,
}

#: A set representing the required metadata during BEL document parsing
REQUIRED_METADATA = {
    METADATA_NAME,
    METADATA_VERSION,
    METADATA_DESCRIPTION,
    METADATA_AUTHORS,
    METADATA_CONTACT
}
