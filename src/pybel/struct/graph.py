# -*- coding: utf-8 -*-

"""Contains the main data structure for PyBEL"""

import logging
import warnings
from copy import deepcopy

import networkx
from six import string_types

from .operations import left_full_join, left_node_intersection_join, left_outer_join
from ..canonicalize import edge_to_bel
from ..constants import *
from ..dsl import activity
from ..dsl.nodes import BaseAbundance, BaseEntity
from ..utils import get_version, hash_edge

__all__ = [
    'BELGraph',
]

log = logging.getLogger(__name__)

RESOURCE_DICTIONARY_NAMES = (
    GRAPH_NAMESPACE_URL,
    GRAPH_NAMESPACE_OWL,
    GRAPH_NAMESPACE_PATTERN,
    GRAPH_ANNOTATION_URL,
    GRAPH_ANNOTATION_OWL,
    GRAPH_ANNOTATION_PATTERN,
    GRAPH_ANNOTATION_LIST,
)


def _clean_annotations(annotations_dict):
    """Fixes formatting of annotation dict

    :type annotations_dict: dict[str,str] or dict[str,set] or dict[str,dict[str,bool]]
    :rtype: dict[str,dict[str,bool]]
    """
    return {
        key: (
            values if isinstance(values, dict) else
            {v: True for v in values} if isinstance(values, set) else
            {values: True}
        )
        for key, values in annotations_dict.items()
    }


class BELGraph(networkx.MultiDiGraph):
    """This class represents biological knowledge assembled in BEL as a network by extending the
    :class:`networkx.MultiDiGraph`.
    """

    def __init__(self, name=None, version=None, description=None, authors=None, contact=None, license=None,
                 copyright=None, disclaimer=None, data=None, **kwargs):
        """The default constructor parses a BEL graph using the built-in :mod:`networkx` methods. For IO, see
        the :mod:`pybel.io` module

        :param str name: The graph's name
        :param str version: The graph's version. Recommended to use `semantic versioning <http://semver.org/>`_ or
                            ``YYYYMMDD`` format.
        :param str description: A description of the graph
        :param str authors: The authors of this graph
        :param str contact: The contact email for this graph
        :param str license: The license for this graph
        :param str copyright: The copyright for this graph
        :param str disclaimer: The disclaimer for this graph
        :param data: initial graph data to pass to :class:`networkx.MultiDiGraph`
        :param kwargs: keyword arguments to pass to :class:`networkx.MultiDiGraph`
        """
        super(BELGraph, self).__init__(data=data, **kwargs)

        self._warnings = []

        #: Maps node sha512 hashes to their actual nodes
        self.hash_to_node = {}  # TODO

        if GRAPH_METADATA not in self.graph:
            self.graph[GRAPH_METADATA] = {}

        if name:
            self.name = name

        if version:
            self.version = version

        if description:
            self.description = description

        if authors:
            self.authors = authors

        if contact:
            self.contact = contact

        if license:
            self.license = license

        if copyright:
            self.copyright = copyright

        if disclaimer:
            self.disclaimer = disclaimer

        if GRAPH_PYBEL_VERSION not in self.graph:
            self.graph[GRAPH_PYBEL_VERSION] = get_version()

        for resource_dict in RESOURCE_DICTIONARY_NAMES:
            if resource_dict not in self.graph:
                self.graph[resource_dict] = {}

        if GRAPH_UNCACHED_NAMESPACES not in self.graph:
            self.graph[GRAPH_UNCACHED_NAMESPACES] = set()

    @property
    def document(self):
        """A dictionary holding the metadata from the "Document" section of the BEL script. All keys are normalized
        according to :data:`pybel.constants.DOCUMENT_KEYS`

        :rtype: dict[str,str]
        """
        return self.graph[GRAPH_METADATA]

    @property
    def name(self, *attrs):  # Needs *attrs since it's an override
        """The graph's name, from the ``SET DOCUMENT Name = "..."`` entry in the source BEL script

        :rtype: str
        """
        return self.document.get(METADATA_NAME)

    @name.setter
    def name(self, *attrs, **kwargs):  # Needs *attrs and **kwargs since it's an override
        self.document[METADATA_NAME] = attrs[0]

    @property
    def version(self):
        """The graph's version, from the ``SET DOCUMENT Version = "..."`` entry in the source BEL script

        :rtype: str
        """
        return self.document.get(METADATA_VERSION)

    @version.setter
    def version(self, version):
        self.document[METADATA_VERSION] = version

    @property
    def description(self):
        """The graph's description, from the ``SET DOCUMENT Description = "..."`` entry in the source BEL Script

        :rtype: str
        """
        return self.document.get(METADATA_DESCRIPTION)

    @description.setter
    def description(self, description):
        self.document[METADATA_DESCRIPTION] = description

    @property
    def authors(self):
        """The graph's description, from the ``SET DOCUMENT Authors = "..."`` entry in the source BEL Script

        :rtype: str
        """
        return self.document[METADATA_AUTHORS]

    @authors.setter
    def authors(self, authors):
        self.document[METADATA_AUTHORS] = authors

    @property
    def contact(self):
        """The graph's description, from the ``SET DOCUMENT ContactInfo = "..."`` entry in the source BEL Script

        :rtype: str
        """
        return self.document.get(METADATA_CONTACT)

    @contact.setter
    def contact(self, contact):
        self.document[METADATA_CONTACT] = contact

    @property
    def license(self):
        """The graph's license, from the `SET DOCUMENT Licenses = "..."`` entry in the source BEL Script

        :rtype: Optional[str]
        """
        return self.document.get(METADATA_LICENSES)

    @license.setter
    def license(self, license):
        self.document[METADATA_LICENSES] = license

    @property
    def copyright(self):
        """The graph's copyright, from the `SET DOCUMENT Copyright = "..."`` entry in the source BEL Script

        :rtype: Optional[str]
        """
        return self.document.get(METADATA_COPYRIGHT)

    @copyright.setter
    def copyright(self, copyright):
        self.document[METADATA_COPYRIGHT] = copyright

    @property
    def disclaimer(self):
        """The graph's disclaimer, from the `SET DOCUMENT Disclaimer = "..."`` entry in the source BEL Script

        :rtype: Optional[str]
        """
        return self.document.get(METADATA_DISCLAIMER)

    @disclaimer.setter
    def disclaimer(self, disclaimer):
        self.document[METADATA_DISCLAIMER] = disclaimer

    @property
    def namespace_url(self):
        """A dictionary mapping the keywords used to create this graph to the URLs of the BELNS files from the
        ``DEFINE NAMESPACE [key] AS URL "[value]"`` entries in the definitions section.

        :rtype: dict[str,str]
        """
        return self.graph[GRAPH_NAMESPACE_URL]

    @property
    def namespace_owl(self):
        """A dictionary mapping the keywords used to create this graph to the URLs of the OWL files from the
        ``DEFINE NAMESPACE [key] AS OWL "[value]"`` entries in the definitions section

        :rtype: dict[str,str]
        """
        return self.graph[GRAPH_NAMESPACE_OWL]

    @property
    def defined_namespace_keywords(self):
        """Returns the set of all keywords defined as namespaces in this graph

        :rtype: set[str]
        """
        return (
                set(self.namespace_pattern) |
                set(self.namespace_url) |
                set(self.namespace_owl)
        )

    @property
    def uncached_namespaces(self):
        """Returns a list of namespaces's URLs that are present in the graph, but cannot be cached due to their
        corresponding resources' cachable flags being set to "no."

        :rtype: set[str]
        """
        return self.graph[GRAPH_UNCACHED_NAMESPACES]

    @property
    def namespace_pattern(self):
        """A dictionary mapping the namespace keywords used to create this graph to their regex patterns from the
        ``DEFINE NAMESPACE [key] AS PATTERN "[value]"`` entries in the definitions section

        :rtype: dict[str,str]
        """
        return self.graph[GRAPH_NAMESPACE_PATTERN]

    @property
    def annotation_url(self):
        """A dictionary mapping the annotation keywords used to create this graph to the URLs of the BELANNO files
        from the ``DEFINE ANNOTATION [key] AS URL "[value]"`` entries in the definitions section

        :rtype: dict[str,str]
        """
        return self.graph[GRAPH_ANNOTATION_URL]

    @property
    def annotation_owl(self):
        """A dictionary mapping the annotation keywords used to creat ethis graph to the URLs of the OWL files
        from the ``DEFINE ANNOTATION [key] AS OWL "[value]"`` entries in the definitions section

        :rtype: dict[str,str]
        """
        return self.graph[GRAPH_ANNOTATION_OWL]

    @property
    def annotation_pattern(self):
        """A dictionary mapping the annotation keywords used to create this graph to their regex patterns
        from the ``DEFINE ANNOTATION [key] AS PATTERN "[value]"`` entries in the definitions section

        :rtype: dict[str,str]
        """
        return self.graph[GRAPH_ANNOTATION_PATTERN]

    @property
    def annotation_list(self):
        """A dictionary mapping the keywords of locally defined annotations to a set of their values
        from the ``DEFINE ANNOTATION [key] AS LIST {"[value]", ...}`` entries in the definitions section

        :rtype: dict[str,set[str]]
        """
        return self.graph[GRAPH_ANNOTATION_LIST]

    @property
    def defined_annotation_keywords(self):
        """Returns the set of all keywords defined as annotations in this graph

        :rtype: set[str]
        """
        return (
                set(self.annotation_pattern) |
                set(self.annotation_url) |
                set(self.annotation_owl) |
                set(self.annotation_list)
        )

    @property
    def pybel_version(self):
        """Stores the version of PyBEL with which this graph was produced as a string

        :rtype: str
        """
        return self.graph[GRAPH_PYBEL_VERSION]

    @property
    def warnings(self):
        """Warnings are stored in a list of 4-tuples that is a property of the graph object.
        This tuple respectively contains the line number, the line text, the exception instance, and the context
        dictionary from the parser at the time of error.

        :rtype: list[tuple[int,str,Exception,dict[str,str]]]
        """
        return self._warnings

    def __str__(self):
        """Stringifies this graph as its name and version pair"""
        return '{} v{}'.format(self.name, self.version)

    def skip_storing_namespace(self, namespace):
        """Checks if the namespace should be skipped

        :param Optional[str] namespace:
        :rtype: bool
        """
        return (
                namespace is not None and
                namespace in self.namespace_url and
                self.namespace_url[namespace] in self.uncached_namespaces
        )

    def add_warning(self, line_number, line, exception, context=None):
        """Adds a warning to the internal warning log in the graph, with optional context information"""
        self.warnings.append((line_number, line, exception, {} if context is None else context))

    def add_unqualified_edge(self, u, v, relation):
        """Adds unique edge that has no annotations

        :param u: Either a PyBEL node tuple or PyBEL node data dictionary representing the source node
        :type u: tuple or dict
        :param v: Either a PyBEL node tuple or PyBEL node data dictionary representing the target node
        :type v: tuple or dict
        :param str relation: A relationship label from :mod:`pybel.constants`

        :return: The hash of this edge
        :rtype: str
        """
        key = unqualified_edge_code[relation]

        if not isinstance(u, BaseEntity):
            raise TypeError

        if not isinstance(v, BaseEntity):
            raise TypeError

        self.add_node_from_data(u)
        self.add_node_from_data(v)

        attr = {RELATION: relation}
        attr[HASH] = hash_edge(u, v, attr)

        if not self.has_edge(u, v, key):
            self.add_edge(u, v, key=key, **attr)

        return attr[HASH]

    def add_transcription(self, u, v):
        """Adds a transcription relation from a gene to an RNA or miRNA node

        :param u: Either a PyBEL node tuple or PyBEL node data dictionary representing the source node
        :type u: tuple or dict
        :param v: Either a PyBEL node tuple or PyBEL node data dictionary representing the target node
        :type v: tuple or dict
        """
        self.add_unqualified_edge(u, v, TRANSCRIBED_TO)

    def add_translation(self, u, v):
        """Adds a translation relation from a RNA to a protein

        :param u: Either a PyBEL node tuple or PyBEL node data dictionary representing the source node
        :type u: tuple or dict
        :param v: Either a PyBEL node tuple or PyBEL node data dictionary representing the target node
        :type v: tuple or dict
        """
        self.add_unqualified_edge(u, v, TRANSLATED_TO)

    def _add_two_way_unqualified_edge(self, u, v, relation):
        """Adds an unqualified edge both ways"""
        self.add_unqualified_edge(u, v, relation)
        self.add_unqualified_edge(v, u, relation)

    def add_equivalence(self, u, v):
        """Adds two equivalence relations for the nodes

        :param u: Either a PyBEL node tuple or PyBEL node data dictionary representing the source node
        :type u: tuple or dict
        :param v: Either a PyBEL node tuple or PyBEL node data dictionary representing the target node
        :type v: tuple or dict
        """
        self._add_two_way_unqualified_edge(u, v, EQUIVALENT_TO)

    def add_orthology(self, u, v):
        """Adds two orthology relations for the nodes

        :param u: Either a PyBEL node tuple or PyBEL node data dictionary representing the source node
        :type u: tuple or dict
        :param v: Either a PyBEL node tuple or PyBEL node data dictionary representing the target node
        :type v: tuple or dict
        """
        self._add_two_way_unqualified_edge(u, v, ORTHOLOGOUS)

    def add_is_a(self, u, v):
        """Adds an isA relationship such that u isA v.

        :param u: Either a PyBEL node tuple or PyBEL node data dictionary representing the source node
        :type u: tuple or dict
        :param v: Either a PyBEL node tuple or PyBEL node data dictionary representing the target node
        :type v: tuple or dict
        """
        self.add_unqualified_edge(u, v, IS_A)

    def iter_node_data_pairs(self):
        """Iterates over pairs of nodes and their data dictionaries

        :rtype: iter[tuple[tuple,dict]]
        """
        return self.nodes_iter(data=True)

    def iter_data(self):
        """Iterates over the node data dictionaries

        :rtype: iter[dict]
        """
        return (
            data
            for _, data in self.nodes_iter(data=True)
        )

    def add_node_from_data(self, node):
        """Converts a PyBEL node data dictionary to a canonical PyBEL node tuple and ensures it is in the graph.

        :param pybel.dsl.nodes.BaseEntity or dict node: A PyBEL node data dictionary
        :return: A PyBEL node tuple
        """
        if node is None:
            raise ValueError('can not add None')

        if not isinstance(node, BaseEntity):
            raise TypeError('should be using base entity {} {}'.format(node.__class__.__name__, node))

        if node in self:
            return

        super(BELGraph, self).add_node(node, **node)  # be a little redundant, for backwards compatability

        variants = node.get(VARIANTS)
        if variants is not None:
            parent_node = node.get_parent()
            self.add_unqualified_edge(parent_node, node, HAS_VARIANT)

        elif MEMBERS in node:
            for member in node[MEMBERS]:
                self.add_unqualified_edge(node, member, HAS_COMPONENT)

        elif PRODUCTS in node and REACTANTS in node:
            for reactant in node[REACTANTS]:
                self.add_unqualified_edge(node, reactant, HAS_REACTANT)

            for product in node[PRODUCTS]:
                self.add_unqualified_edge(node, product, HAS_PRODUCT)

    def has_node_with_data(self, attr_dict):
        """Checks if this graph has a node with the given data dictionary

        :param BaseEntity attr_dict: A PyBEL node data dictionary
        :rtype: bool
        """
        warnings.warn('just use function `has_node`', DeprecationWarning)

        return self.has_node(attr_dict)

    def add_qualified_edge(self, u, v, relation, evidence, citation, annotations=None, subject_modifier=None,
                           object_modifier=None, **attr):
        """Adds an edge, qualified with a relation, evidence, citation, and optional annotations, subject modifications,
        and object modifications

        :param BaseEntity u: Either a PyBEL node tuple or PyBEL node data dictionary representing the source node
        :param BaseEntity v: Either a PyBEL node tuple or PyBEL node data dictionary representing the target node
        :param str relation: The type of relation this edge represents
        :param str evidence: The evidence string from an article
        :param dict[str,str] or str citation: The citation data dictionary for this evidence. If a string is given,
                                                assumes it's a PubMed identifier and auto-fills the citation type.
        :param annotations: The annotations data dictionary
        :type annotations: Optional[dict[str,str] or dict[str,set] or dict[str,dict[str,bool]]]
        :param Optional[dict] subject_modifier: The modifiers (like activity) on the subject node. See data model documentation.
        :param Optional[dict] object_modifier: The modifiers (like activity) on the object node. See data model documentation.

        :return: The hash of the edge
        :rtype: str
        """
        if u is None or v is None:
            raise ValueError('arguments cant be none')

        attr.update({
            RELATION: relation,
            EVIDENCE: evidence,
        })

        if isinstance(citation, string_types):
            attr[CITATION] = {
                CITATION_TYPE: CITATION_TYPE_PUBMED,
                CITATION_REFERENCE: citation
            }
        elif isinstance(citation, dict):
            attr[CITATION] = citation
        else:
            raise TypeError

        if annotations:  # clean up annotations
            attr[ANNOTATIONS] = _clean_annotations(annotations)

        if subject_modifier:
            attr[SUBJECT] = subject_modifier

        if object_modifier:
            attr[OBJECT] = object_modifier

        self.add_node_from_data(u)
        self.add_node_from_data(v)

        attr[HASH] = hash_edge(u, v, attr)

        if not isinstance(u, BaseEntity):
            raise RuntimeError('u is not a bel thing! {} {}'.format(u.__class__.__name__, u))

        self.add_edge(u, v, **attr)

        return attr[HASH]

    def add_inhibits(self, u, v, evidence, citation, annotations=None, object_modifier=None):
        """A more specific version of add_qualified edge that automatically populates the relation and object
        modifier

        :param BaseEntity u: Either a PyBEL node tuple or PyBEL node data dictionary representing the source node
        :param BaseEntity or dict v: Either a PyBEL node tuple or PyBEL node data dictionary representing the target node
        :param str evidence: The evidence string from an article
        :param dict[str,str] or str citation: The citation data dictionary for this evidence. If a string is given,
                                                assumes it's a PubMed identifier and autofills the citation type.
        :param annotations: The annotations data dictionary
        :type annotations: Optional[dict[str,str] or dict[str,set] or dict[str,dict[str,bool]]]
        :param Optional[dict] object_modifier: A non-default activity.

        :return: The hash of the edge
        :rtype: str
        """
        return self.add_qualified_edge(
            u,
            v,
            relation=DECREASES,
            evidence=evidence,
            citation=citation,
            annotations=annotations,
            object_modifier=object_modifier or activity()
        )

    def has_edge_citation(self, u, v, key):
        """Does the given edge have a citation?

        :param BaseEntity u: A PyBEL node
        :param BaseEntity v: A PyBEL node
        :param key: The edge key
        :rtype: bool
        """
        return CITATION in self.edge[u][v][key]

    def get_edge_citation(self, u, v, key):
        """Gets the citation for a given edge

        :param BaseEntity u: A PyBEL node
        :param BaseEntity v: A PyBEL node
        :param key: The edge key
        :rtype: Optional[dict]
        """
        return self.edge[u][v][key].get(CITATION)

    def has_edge_evidence(self, u, v, key):
        """Does the given edge have evidence?

        :param BaseEntity u: A PyBEL node
        :param BaseEntity v: A PyBEL node
        :param key: The edge key
        :rtype: boolean
        """
        return EVIDENCE in self.edge[u][v][key]

    def get_edge_evidence(self, u, v, key):
        """Gets the evidence for a given edge

        :param BaseEntity u: A PyBEL node
        :param BaseEntity v: A PyBEL node
        :param key: The edge key
        :rtype: Optional[str]
        """
        return self.edge[u][v][key].get(EVIDENCE)

    def get_edge_annotations(self, u, v, key):
        """Gets the annotations for a given edge

        :param BaseEntity u: A PyBEL node
        :param BaseEntity v: A PyBEL node
        :param key: The edge key
        :rtype: Optional[dict]
        """
        return self.edge[u][v][key].get(ANNOTATIONS)

    def get_node_description(self, node):
        """Gets the description for a given node

        :rtype: Optional[str]
        """
        return self.node[node].get(DESCRIPTION)

    def has_node_description(self, node):
        """Returns if a node description is already present

        :param BaseEntity node: A PyBEL node
        :rtype: bool
        """
        return DESCRIPTION in self.node[node]

    def set_node_description(self, node, description):
        """Sets the description for a given node

        :param BaseEntity node: A PyBEL node
        :type description: str
        """
        self.node[node][DESCRIPTION] = description

    def __add__(self, other):
        """Creates a deep copy of this graph and full joins another graph with it using
        :func:`pybel.struct.left_full_join`.

        :param BELGraph other: Another BEL graph
        :rtype: BELGraph

        Example usage:

        >>> import pybel
        >>> g = pybel.from_path('...')
        >>> h = pybel.from_path('...')
        >>> k = g + h
        """
        if not isinstance(other, BELGraph):
            raise TypeError('{} is not a {}'.format(other, self.__class__.__name__))

        result = deepcopy(self)
        left_full_join(result, other)
        return result

    def __iadd__(self, other):
        """Full joins another graph into this one using :func:`pybel.struct.left_full_join`.

        :param BELGraph other: Another BEL graph
        :rtype: BELGraph

        Example usage:

        >>> import pybel
        >>> g = pybel.from_path('...')
        >>> h = pybel.from_path('...')
        >>> g += h
        """
        if not isinstance(other, BELGraph):
            raise TypeError('{} is not a {}'.format(other, self.__class__.__name__))

        left_full_join(self, other)
        return self

    def __and__(self, other):
        """Creates a deep copy of this graph and outer joins another graph with it using
        :func:`pybel.struct.left_outer_join`.

        :param BELGraph other: Another BEL graph
        :rtype: BELGraph

        Example usage:

        >>> import pybel
        >>> g = pybel.from_path('...')
        >>> h = pybel.from_path('...')
        >>> k = g & h
        """
        if not isinstance(other, BELGraph):
            raise TypeError('{} is not a {}'.format(other, self.__class__.__name__))

        result = deepcopy(self)
        left_outer_join(result, other)
        return result

    def __iand__(self, other):
        """Outer joins another graph into this one using :func:`pybel.struct.left_outer_join`.

        :param BELGraph other: Another BEL graph
        :rtype: BELGraph

        Example usage:

        >>> import pybel
        >>> g = pybel.from_path('...')
        >>> h = pybel.from_path('...')
        >>> g &= h
        """
        if not isinstance(other, BELGraph):
            raise TypeError('{} is not a {}'.format(other, self.__class__.__name__))

        left_outer_join(self, other)
        return self

    def __xor__(self, other):
        """Node intersection joins another graph using :func:`pybel.struct.left_node_intersection_join`

        :param BELGraph other: Another BEL graph

        Example usage:

        >>> import pybel
        >>> g = pybel.from_path('...')
        >>> h = pybel.from_path('...')
        >>> k = g ^ h
        """
        if not isinstance(other, BELGraph):
            raise TypeError('{} is not a {}'.format(other, self.__class__.__name__))

        return left_node_intersection_join(self, other)

    def node_to_bel(self, n):
        """Serializes a node as BEL

        :param tuple n: A PyBEL node tuple
        :rtype: str
        """
        if isinstance(n, BaseAbundance):
            return n.as_bel()

        raise RuntimeError('stop using tuples!')

    def edge_to_bel(self, u, v, data, sep=None):
        """Serializes a pair of nodes and related edge data as a BEL relation

        :param BaseEntity u: A PyBEL node tuple for the soure node
        :param BaseEntity v: A PyBEL node tuple for the target node
        :param dict data: A PyBEL edge data dictionary
        :param str sep: The separator between the source, relation, and target. Defaults to ' '
        :rtype: str
        """
        return edge_to_bel(u, v, data=data, sep=sep)

    def _equivalent_node_iterator_helper(self, node, visited):
        """Iterates over nodes and their data that are equal to the given node, starting with the original

        :param tuple node: A PyBEL node tuple
        :rtype: iter[tuple]
        """
        for v in self.edge[node]:
            if v in visited:
                continue

            if unqualified_edge_code[EQUIVALENT_TO] not in self.edge[node][v]:
                continue

            yield v
            visited.add(v)

            for w in self._equivalent_node_iterator_helper(v, visited):
                yield w

    def iter_equivalent_nodes(self, node):
        """Iterates over node tuples that are equivalent to the given node, including the original

        :param tuple node: A PyBEL node tuple
        :rtype: iter[tuple]
        """
        yield node

        for n in self._equivalent_node_iterator_helper(node, {node}):
            yield n

    def get_equivalent_nodes(self, node):
        """Gets a set of equivalent nodes to this node. Does not include the given node.

        :param tuple node: A PyBEL node tuple
        :rtype: set[tuple]
        """
        return set(self.iter_equivalent_nodes(node))

    def _node_has_namespace_helper(self, node, namespace):
        """Check that the node has namespace information. Might have cross references in future

        :param tuple node: A PyBEL node tuple
        :rtype: bool
        """
        if isinstance(node, dict):
            return namespace == node.get(NAMESPACE)

        return namespace == self.node[node].get(NAMESPACE)

    def node_has_namespace(self, node, namespace):
        """Does the node have the given namespace? This also should look in the equivalent nodes.

        :param tuple node: A PyBEL node tuple
        :param str namespace: A namespace
        :rtype: bool
        """
        for n in self.iter_equivalent_nodes(node):
            if self._node_has_namespace_helper(n, namespace):
                return True

        return False
