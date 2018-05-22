# -*- coding: utf-8 -*-

"""This module contains the SQLAlchemy database models that support the definition cache and graph cache."""

import datetime
import hashlib

from sqlalchemy import (
    Boolean, Column, Date, DateTime, ForeignKey, Integer, LargeBinary, String, Table, Text,
    UniqueConstraint,
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import backref, relationship

from .utils import int_or_str
from ..constants import *
from ..dsl import (
    abundance, bioprocess, complex_abundance, composite_abundance, fragment, fusion_range, gene, gene_fusion, gmod,
    hgvs, mirna, missing_fusion_range, named_complex_abundance, pathology, pmod, protein, protein_fusion, reaction, rna,
    rna_fusion,
)
from ..io.gpickle import from_bytes, to_bytes
from ..utils import hash_citation, hash_edge

__all__ = [
    'Base',
    'Namespace',
    'NamespaceEntry',
    'Annotation',
    'AnnotationEntry',
    'Network',
    'Node',
    'Modification',
    'Author',
    'Citation',
    'Evidence',
    'Edge',
    'Property',
]

BASE_TABLE_NAME = 'bel'
NAMESPACE_TABLE_NAME = f'{BASE_TABLE_NAME}_namespace'
NAMESPACE_ENTRY_TABLE_NAME = f'{BASE_TABLE_NAME}_namespaceEntry'
NAMESPACE_HIERARCHY_TABLE_NAME = f'{BASE_TABLE_NAME}_namespace_hierarchy'

ANNOTATION_TABLE_NAME = f'{BASE_TABLE_NAME}_annotation'
ANNOTATION_ENTRY_TABLE_NAME = f'{BASE_TABLE_NAME}_annotationEntry'
ANNOTATION_HIERARCHY_TABLE_NAME = f'{BASE_TABLE_NAME}_annotation_hierarchy'

NETWORK_TABLE_NAME = f'{BASE_TABLE_NAME}_network'
NETWORK_NODE_TABLE_NAME = f'{BASE_TABLE_NAME}_network_node'
NETWORK_EDGE_TABLE_NAME = f'{BASE_TABLE_NAME}_network_edge'
NETWORK_NAMESPACE_TABLE_NAME = f'{BASE_TABLE_NAME}_network_namespace'
NETWORK_ANNOTATION_TABLE_NAME = f'{BASE_TABLE_NAME}_network_annotation'
NETWORK_CITATION_TABLE_NAME = f'{BASE_TABLE_NAME}_network_citation'

NODE_TABLE_NAME = f'{BASE_TABLE_NAME}_node'
NODE_MODIFICATION_TABLE_NAME = f'{BASE_TABLE_NAME}_node_modification'

MODIFICATION_TABLE_NAME = f'{BASE_TABLE_NAME}_modification'

EDGE_TABLE_NAME = f'{BASE_TABLE_NAME}_edge'
EDGE_ANNOTATION_TABLE_NAME = f'{BASE_TABLE_NAME}_edge_annotationEntry'
EDGE_PROPERTY_TABLE_NAME = f'{BASE_TABLE_NAME}_edge_property'

AUTHOR_TABLE_NAME = f'{BASE_TABLE_NAME}_author'
AUTHOR_CITATION_TABLE_NAME = f'{BASE_TABLE_NAME}_author_citation'

CITATION_TABLE_NAME = f'{BASE_TABLE_NAME}_citation'
EVIDENCE_TABLE_NAME = f'{BASE_TABLE_NAME}_evidence'
PROPERTY_TABLE_NAME = f'{BASE_TABLE_NAME}_property'

LONGBLOB = 4294967295

Base = declarative_base()

_fn_dsl_map = {
    PROTEIN: protein,
    GENE: gene,
    RNA: rna,
    MIRNA: mirna,
    PATHOLOGY: pathology,
    BIOPROCESS: bioprocess,
    ABUNDANCE: abundance
}

namespace_hierarchy = Table(
    NAMESPACE_HIERARCHY_TABLE_NAME,
    Base.metadata,
    Column('left_id', Integer, ForeignKey(f'{NAMESPACE_ENTRY_TABLE_NAME}.id'), primary_key=True),
    Column('right_id', Integer, ForeignKey(f'{NAMESPACE_ENTRY_TABLE_NAME}.id'), primary_key=True)
)

annotation_hierarchy = Table(
    ANNOTATION_HIERARCHY_TABLE_NAME,
    Base.metadata,
    Column('left_id', Integer, ForeignKey(f'{ANNOTATION_ENTRY_TABLE_NAME}.id'), primary_key=True),
    Column('right_id', Integer, ForeignKey(f'{ANNOTATION_ENTRY_TABLE_NAME}.id'), primary_key=True)
)


class Namespace(Base):
    """Represents a BEL Namespace.

    Example:

    .. code-block:: python

        hgnc_namespace = Namespace(
            miriam_id='MIR:00000080',
            name='HGNC',
            keyword='hgnc',
            description='The HGNC (HUGO Gene Nomenclature Committee) provides an approved gene name and symbol (short-form abbreviation) for each known human gene. All approved symbols are stored in the HGNC database, and each symbol is unique. HGNC identifiers refer to records in the HGNC symbol database.',
            uri='http://identifiers.org/hgnc/',
        )
    """
    __tablename__ = NAMESPACE_TABLE_NAME

    id = Column(Integer, primary_key=True)

    uploaded = Column(DateTime, nullable=False, default=datetime.datetime.utcnow, doc='The date of upload to PyBEL')

    miriam_id = Column(String(32), nullable=True, index=True, doc='MIRIAM identifier of this namespace')
    keyword = Column(String(16), nullable=True, index=True,
                     doc='Keyword that is used in a BEL file to identify a specific namespace. Should correspond to Identifiers.org "Namespace" field.')
    name = Column(String(255), nullable=True,
                  doc='Name of the given namespace. Should correspond to Identifiers.org "Name" field')
    domain = Column(String(255), nullable=True,
                    doc='Domain for which this namespace is valid. Corresponds to identifiers.org "Definition" field.')


    # A namespace either needs a URL or a pattern
    pattern = Column(String(255), nullable=True, unique=True, index=True,
                     doc="Contains regex pattern for value identification.")

    url = Column(String(255), nullable=True, unique=True, index=True, doc='BELNS Resource location as URL')

    species = Column(String(255), nullable=True, doc='Taxonomy identifiers for which this namespace is valid')
    description = Column(Text, nullable=True, doc='Optional short description of the namespace')
    version = Column(String(255), nullable=True, doc='Version of the namespace')
    created = Column(DateTime, nullable=True, doc='DateTime of the creation of the namespace definition file')
    query_url = Column(Text, nullable=True, doc='URL that can be used to query the namespace (externally from PyBEL)')

    author = Column(String(255), doc='The author of the namespace')
    license = Column(String(255), nullable=True, doc='License information')
    contact = Column(String(255), nullable=True, doc='Contact information')

    citation = Column(String(255))
    citation_description = Column(Text, nullable=True)
    citation_version = Column(String(255), nullable=True)
    citation_published = Column(Date, nullable=True)
    citation_url = Column(String(255), nullable=True)

    entries = relationship('NamespaceEntry', lazy='dynamic', cascade='all, delete-orphan')

    def __str__(self):
        return '[{}]{}'.format(self.id, self.keyword)

    def to_values(self):
        """Returns this namespace as a dictionary of names to their encodings. Encodings are represented as a
        string, and lookup operations take constant time O(8).

        :rtype: dict[str,str]
        """
        return {
            entry.name: entry.encoding if entry.encoding else BELNS_ENCODING_STR
            for entry in self.entries
        }

    def to_tree_list(self):
        """Returns an edge set of the tree represented by this namespace's hierarchy

        :rtype: set[tuple[str,str]]
        """
        return {
            (parent.name, child.name)
            for parent in self.entries
            for child in parent.children
        }

    def to_json(self, include_id=False):
        """Returns the most useful entries as a dictionary

        :param bool include_id: If true, includes the model identifier
        :rtype: dict[str,str]
        """
        result = {
            'keyword': self.keyword,
            'name': self.name,
            'version': self.version,
        }

        if self.url:
            result['url'] = self.url
        else:
            result['pattern'] = self.pattern

        if include_id:
            result['id'] = self.id

        return result


class NamespaceEntry(Base):
    """Represents a name within a BEL namespace"""
    __tablename__ = NAMESPACE_ENTRY_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String(1023), index=True, nullable=False,
                  doc='Name that is defined in the corresponding namespace definition file')
    identifier = Column(String(255), index=True, nullable=True, doc='The database accession number')
    encoding = Column(String(8), nullable=True, doc='The biological entity types for which this name is valid')

    namespace_id = Column(Integer, ForeignKey(f'{NAMESPACE_TABLE_NAME}.id'), nullable=False, index=True)
    namespace = relationship('Namespace')

    children = relationship(
        'NamespaceEntry',
        secondary=namespace_hierarchy,
        primaryjoin=(id == namespace_hierarchy.c.left_id),
        secondaryjoin=(id == namespace_hierarchy.c.right_id),
    )

    __table_args__ = (
        UniqueConstraint(namespace_id, identifier),
    )

    def __str__(self):
        return '[{}]{}:[{}]{}'.format(self.namespace.id, self.namespace.keyword, self.identifier, self.name)

    def _as_dsl(self, fn, **kwargs):
        """

        :param type[BaseAbundance] fn:
        :rtype: BaseAbundance
        """
        return fn(namespace=self.namespace.keyword, name=self.name, identifier=self.identifier, **kwargs)

    def _as_central_dogma(self, fn, variants=None):
        return self._as_dsl(fn, variants=variants)

    def as_gene(self, variants=None):
        """
        :rtype: pybel.dsl.nodes.gene
        """
        return self._as_central_dogma(gene, variants=variants)

    def as_rna(self, variants=None):
        """
        :rtype: pybel.dsl.nodes.rna
        """
        return self._as_central_dogma(rna, variants=variants)

    def as_protein(self, variants=None):
        """
        :rtype: pybel.dsl.nodes.protein
        """
        return self._as_central_dogma(protein, variants=variants)

    def as_complex(self, members):
        return self._as_dsl(complex_abundance, members=members)

    def as_named_complex(self):
        """
        :rtype: pybel.dsl.nodes.named_complex_abundance
        """
        return self._as_dsl(named_complex_abundance)

    def to_json(self, include_id=False):
        """Describes the namespaceEntry as dictionary of Namespace-Keyword and Name.

        :param bool include_id: If true, includes the model identifier
        :rtype: dict[str,str]
        """
        result = {
            NAMESPACE: self.namespace.keyword,
            NAME: self.name
        }

        if self.identifier:
            result[IDENTIFIER] = self.identifier

        if include_id:
            result['id'] = self.id

        return result


class Annotation(Base):
    """Represents a BEL Annotation"""
    __tablename__ = ANNOTATION_TABLE_NAME

    id = Column(Integer, primary_key=True)
    uploaded = Column(DateTime, default=datetime.datetime.utcnow, doc='The date of upload')

    url = Column(String(255), nullable=False, unique=True, index=True,
                 doc='Source url of the given annotation definition file (.belanno)')
    keyword = Column(String(50), index=True, doc='Keyword that is used in a BEL file to identify a specific annotation')
    type = Column(String(255), doc='Annotation type')
    description = Column(Text, nullable=True, doc='Optional short description of the given annotation')
    usage = Column(Text, nullable=True)
    version = Column(String(255), nullable=True, doc='Version of the annotation')
    created = Column(DateTime, doc='DateTime of the creation of the given annotation definition')

    name = Column(String(255), doc='Name of the annotation definition')
    author = Column(String(255), doc='Author information')
    license = Column(String(255), nullable=True, doc='License information')
    contact = Column(String(255), nullable=True, doc='Contact information')

    citation = Column(String(255))
    citation_description = Column(Text, nullable=True)
    citation_version = Column(String(255), nullable=True)
    citation_published = Column(Date, nullable=True)
    citation_url = Column(String(255), nullable=True)

    def get_entries(self):
        """Gets a set of the names of all entries

        :rtype: set[str]
        """
        return {
            entry.name
            for entry in self.entries
        }

    def to_tree_list(self):
        """Returns an edge set of the tree represented by this namespace's hierarchy

        :rtype: set[tuple[str,str]]
        """
        return {
            (parent.name, child.name)
            for parent in self.entries
            for child in parent.children
        }

    def to_json(self, include_id=False):
        """Returns this annotation as a JSON dictionary

        :param bool include_id: If true, includes the model identifier
        :rtype: dict[str,str]
        """
        result = {
            'url': self.url,
            'keyword': self.keyword,
            'version': self.version,
            'name': self.name
        }

        if include_id:
            result['id'] = self.id

        return result

    def __str__(self):
        return self.keyword


class AnnotationEntry(Base):
    """Represents a value within a BEL Annotation"""
    __tablename__ = ANNOTATION_ENTRY_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String(255), nullable=False, index=True,
                  doc='Name that is defined in the corresponding annotation definition file')
    label = Column(Text, nullable=True)

    annotation_id = Column(Integer, ForeignKey(f'{ANNOTATION_TABLE_NAME}.id'), index=True)
    annotation = relationship('Annotation', backref=backref('entries', lazy='dynamic'))

    children = relationship(
        'AnnotationEntry',
        secondary=annotation_hierarchy,
        primaryjoin=(id == annotation_hierarchy.c.left_id),
        secondaryjoin=(id == annotation_hierarchy.c.right_id)
    )

    def to_json(self, include_id=False):
        """Describes the annotationEntry as dictionary of Annotation-Keyword and Annotation-Name.

        :param bool include_id: If true, includes the model identifier
        :rtype: dict[str,str]
        """
        result = {
            'annotation_keyword': self.annotation.keyword,
            'annotation': self.name
        }

        if include_id:
            result['id'] = self.id

        return result

    def __str__(self):
        return '{}:{}'.format(self.annotation, self.name)


network_edge = Table(
    NETWORK_EDGE_TABLE_NAME, Base.metadata,
    Column('network_id', Integer, ForeignKey(f'{NETWORK_TABLE_NAME}.id'), primary_key=True),
    Column('edge_id', Integer, ForeignKey(f'{EDGE_TABLE_NAME}.id'), primary_key=True)
)

network_node = Table(
    NETWORK_NODE_TABLE_NAME, Base.metadata,
    Column('network_id', Integer, ForeignKey(f'{NETWORK_TABLE_NAME}.id'), primary_key=True),
    Column('node_id', Integer, ForeignKey(f'{NODE_TABLE_NAME}.id'), primary_key=True)
)


class Network(Base):
    """Represents a collection of edges, specified by a BEL Script"""
    __tablename__ = NETWORK_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String(255), nullable=False, index=True, doc='Name of the given Network (from the BEL file)')
    version = Column(String(16), nullable=False, doc='Release version of the given Network (from the BEL file)')
    sha512 = Column(String(255), nullable=True, index=True)

    authors = Column(Text, nullable=True, doc='Authors of the underlying BEL file')
    contact = Column(String(255), nullable=True, doc='Contact email from the underlying BEL file')
    description = Column(Text, nullable=True, doc='Descriptive text from the underlying BEL file')
    copyright = Column(Text, nullable=True, doc='Copyright information')
    disclaimer = Column(Text, nullable=True, doc='Disclaimer information')
    licenses = Column(Text, nullable=True, doc='License information')

    created = Column(DateTime, nullable=False, default=datetime.datetime.utcnow)
    blob = Column(LargeBinary(LONGBLOB), doc='A pickled version of this network')

    nodes = relationship('Node', secondary=network_node, lazy='dynamic', backref=backref('networks', lazy='dynamic'))
    edges = relationship('Edge', secondary=network_edge, lazy='dynamic', backref=backref('networks', lazy='dynamic'))

    __table_args__ = (
        UniqueConstraint(name, version),
    )

    def to_json(self, include_id=False):
        """Returns this network as JSON

        :param bool include_id: If true, includes the model identifier
        :rtype: dict[str,str]
        """
        result = {
            METADATA_NAME: self.name,
            METADATA_VERSION: self.version,
        }

        if self.created:
            result['created'] = str(self.created)

        if include_id:
            result['id'] = self.id

        if self.authors:
            result[METADATA_AUTHORS] = self.authors

        if self.contact:
            result[METADATA_CONTACT] = self.contact

        if self.description:
            result[METADATA_DESCRIPTION] = self.description

        if self.copyright:
            result[METADATA_COPYRIGHT] = self.copyright

        if self.disclaimer:
            result[METADATA_DISCLAIMER] = self.disclaimer

        if self.licenses:
            result[METADATA_LICENSES] = self.licenses

        return result

    def __repr__(self):
        return '{} v{}'.format(self.name, self.version)

    def __str__(self):
        return repr(self)

    def as_bel(self):
        """Gets this network and loads it into a :class:`BELGraph`

        :rtype: pybel.BELGraph
        """
        return from_bytes(self.blob)

    def store_bel(self, graph):
        """Inserts a BEL graph and stores its hash

        :param pybel.BELGraph graph: A BEL Graph
        """
        self.blob = to_bytes(graph)
        self.sha512 = hashlib.sha512(self.blob).hexdigest()

    @staticmethod
    def filter_ids(network_ids):
        """

        :param network_ids:
        :return:
        """
        return Network.id.in_(network_ids)


node_modification = Table(
    NODE_MODIFICATION_TABLE_NAME, Base.metadata,
    Column('node_id', Integer, ForeignKey(f'{NODE_TABLE_NAME}.id'), primary_key=True),
    Column('modification_id', Integer, ForeignKey(f'{MODIFICATION_TABLE_NAME}.id'), primary_key=True)
)


class Node(Base):
    """Represents a BEL Term."""

    __tablename__ = NODE_TABLE_NAME

    id = Column(Integer, primary_key=True)

    type = Column(String(255), nullable=False, doc='The type of the represented biological entity e.g. Protein or Gene')
    is_variant = Column(Boolean, default=False, doc='Identifies weather or not the given node is a variant')
    has_fusion = Column(Boolean, default=False, doc='Identifies weather or not the given node is a fusion')
    bel = Column(String(255), nullable=False, doc='Canonical BEL term that represents the given node')
    sha512 = Column(String(255), nullable=True, index=True)

    namespace_entry_id = Column(Integer, ForeignKey(f'{NAMESPACE_ENTRY_TABLE_NAME}.id'), nullable=True)
    namespace_entry = relationship('NamespaceEntry', foreign_keys=[namespace_entry_id])

    modifications = relationship("Modification", secondary=node_modification, lazy='dynamic',
                                 backref=backref('nodes', lazy='dynamic'))

    def __str__(self):
        return self.bel

    def __repr__(self):
        return '<Node {}: {}>'.format(self.sha512[:10], self.bel)

    def _update_to_json_rv(self, result, include_id=False, include_hash=False):
        if include_id:
            result['id'] = self.id

        if include_hash:
            result[HASH] = self.sha512

        return result

    def as_bel(self, include_id=False, include_hash=False):
        """Serializes this node as a PyBEL node data dictionary

        :param bool include_id: Include the database identifier?
        :param bool include_hash: Include the node hash?
        :rtype: BaseEntity
        """
        if self.has_fusion:
            fusion_modification = self.modifications[0]  # should only have one modification
            result = fusion_modification.to_fusion_dsl(self.type)
            self._update_to_json_rv(result, include_id=include_id, include_hash=include_hash)
            return result

        if self.type == REACTION:
            reactants = [
                edge.target.as_bel()
                for edge in self.out_edges.filter(Edge.relation == HAS_REACTANT)
            ]
            products = [
                edge.target.as_bel()
                for edge in self.out_edges.filter(Edge.relation == HAS_PRODUCT)
            ]
            result = reaction(
                reactants=reactants,
                products=products,
            )
            self._update_to_json_rv(result, include_id=include_id, include_hash=include_hash)
            return result

        if self.type == COMPOSITE:
            result = composite_abundance(members=[
                edge.target.as_bel()
                for edge in self.out_edges.filter(Edge.relation == HAS_COMPONENT)
            ])
            self._update_to_json_rv(result, include_id=include_id, include_hash=include_hash)
            return result

        if self.type == COMPLEX or (self.type == COMPLEX and not self.namespace_entry):
            members = [
                edge.target.as_bel()
                for edge in self.out_edges.filter(Edge.relation == HAS_COMPONENT)
            ]

            if not members and self.namespace_entry:
                result = self.namespace_entry.as_named_complex()

            elif members and self.namespace_entry:
                result = self.namespace_entry.as_complex(members=members)

            elif members and not self.namespace_entry:  # means there are members, otherwise there's something wrong
                result = complex_abundance(members=members)

            else:
                raise ValueError

            self._update_to_json_rv(result, include_id=include_id, include_hash=include_hash)
            return result

        if self.is_variant:
            variants = [
                modification.to_variant()
                for modification in self.modifications
            ]

            if self.type == PROTEIN:
                result = self.namespace_entry.as_protein(variants=variants)
            elif self.type == RNA:
                result = self.namespace_entry.as_rna(variants=variants)
            elif self.type == GENE:
                result = self.namespace_entry.as_gene(variants=variants)
            else:
                raise ValueError

            self._update_to_json_rv(result, include_id=include_id, include_hash=include_hash)
            return result

        dsl_func = _fn_dsl_map[self.type]

        result = dsl_func(
            namespace=self.namespace_entry.namespace.keyword,
            name=self.namespace_entry.name,
            identifier=self.namespace_entry.identifier
        )
        self._update_to_json_rv(result, include_id=include_id, include_hash=include_hash)
        return result

    def to_tuple(self):
        """Converts this node to a PyBEL tuple

        :rtype: tuple
        """
        return self.to_json().as_tuple()


class Modification(Base):
    """The modifications that are present in the network are stored in this table."""
    __tablename__ = MODIFICATION_TABLE_NAME

    id = Column(Integer, primary_key=True)

    type = Column(String(255), nullable=False, doc='Type of the stored modification e.g. Fusion, gmod, pmod, etc')

    variantString = Column(String(255), nullable=True, doc='HGVS string if sequence modification')

    p3_partner_id = Column(Integer, ForeignKey(f'{NAMESPACE_ENTRY_TABLE_NAME}.id'), nullable=True)
    p3_partner = relationship("NamespaceEntry", foreign_keys=[p3_partner_id])

    p3_reference = Column(String(10), nullable=True)
    p3_start = Column(String(255), nullable=True)
    p3_stop = Column(String(255), nullable=True)

    p5_partner_id = Column(Integer, ForeignKey(f'{NAMESPACE_ENTRY_TABLE_NAME}.id'), nullable=True)
    p5_partner = relationship("NamespaceEntry", foreign_keys=[p5_partner_id])

    p5_reference = Column(String(10), nullable=True)
    p5_start = Column(String(255), nullable=True)
    p5_stop = Column(String(255), nullable=True)

    identifier_id = Column(Integer, ForeignKey(f'{NAMESPACE_ENTRY_TABLE_NAME}.id'), nullable=True)
    identifier = relationship("NamespaceEntry", foreign_keys=[identifier_id])

    residue = Column(String(3), nullable=True, doc='Three letter amino acid code if PMOD')
    position = Column(Integer, nullable=True, doc='Position of PMOD or GMOD')

    sha512 = Column(String(255), index=True)

    def _get_range3p(self):
        """
        :rtype: pybel.dsl.nodes.FusionRangeBase
        """
        if not self.p3_reference:
            return missing_fusion_range()

        return fusion_range(
            reference=str(self.p3_reference),
            start=int_or_str(self.p3_start),
            stop=int_or_str(self.p3_stop),
        )

    def _get_range5p(self):
        """
        :rtype: pybel.dsl.nodes.FusionRangeBase
        """
        if not self.p5_reference:
            return missing_fusion_range()

        return fusion_range(
            reference=str(self.p5_reference),
            start=int_or_str(self.p5_start),
            stop=int_or_str(self.p5_stop),
        )

    def to_fusion_dsl(self, fn):
        """
        :param str fn: Either GENE, RNA, or PROTEIN
        :rtype: FusionBase
        """
        range5p = self._get_range5p()
        range3p = self._get_range3p()

        if fn == PROTEIN:
            return protein_fusion(
                partner5p=self.p5_partner.as_protein(),
                partner3p=self.p3_partner.as_protein(),
                range5p=range5p,
                range3p=range3p,
            )
        elif fn == RNA:
            return rna_fusion(
                partner5p=self.p5_partner.as_rna(),
                partner3p=self.p3_partner.as_rna(),
                range5p=range5p,
                range3p=range3p,
            )
        elif fn == GENE:
            return gene_fusion(
                partner5p=self.p5_partner.as_gene(),
                partner3p=self.p3_partner.as_gene(),
                range5p=range5p,
                range3p=range3p,
            )
        else:
            raise ValueError

    def to_fragment(self):
        """
        :rtype: pybel.dsl.nodes.fragment
        """
        return fragment(
            start=int_or_str(self.p3_start),
            stop=int_or_str(self.p3_stop)
        )

    def to_variant(self):
        """Recreates a is_variant dictionary for :class:`BELGraph`

        :return: Dictionary that describes a variant or a fusion.
        :rtype: pybel.dsl.nodes.Variant
        """
        if self.type == FRAGMENT:
            return self.to_fragment()

        if self.type == HGVS:
            return hgvs(self.variantString)

        if self.type == GMOD:
            return gmod(
                namespace=self.identifier.namespace.keyword,
                name=self.identifier.name,
                identifier=self.identifier.identifier
            )

        if self.type == PMOD:
            return pmod(
                namespace=self.identifier.namespace.keyword,
                name=self.identifier.name,
                identifier=self.identifier.identifier,
                position=self.position,
                code=self.residue
            )

        raise ValueError


author_citation = Table(
    AUTHOR_CITATION_TABLE_NAME, Base.metadata,
    Column('author_id', Integer, ForeignKey(f'{AUTHOR_TABLE_NAME}.id'), primary_key=True),
    Column('citation_id', Integer, ForeignKey(f'{CITATION_TABLE_NAME}.id'), primary_key=True)
)


class Author(Base):
    """Contains all author names."""
    __tablename__ = AUTHOR_TABLE_NAME

    id = Column(Integer, primary_key=True)

    sha512 = Column(String(255), index=True, unique=True, nullable=False)
    name = Column(String(255), nullable=False, index=True, doc='The name of the author')

    @staticmethod
    def _hash_name(name):
        return hashlib.sha512(name.strip().casefold().encode('utf8')).hexdigest()

    @staticmethod
    def from_name(name):
        """Makes an author model from the name and wraps hashing the name

        :param str name: The name of the author
        :rtype: Author
        """
        return Author(name=name, sha512=Author._hash_name(name))

    @staticmethod
    def filter_name(name):
        return Author.sha512 == Author._hash_name(name)

    def __str__(self):
        return self.name


class Citation(Base):
    """The information about the citations that are used to prove a specific relation are stored in this table."""
    __tablename__ = CITATION_TABLE_NAME

    id = Column(Integer, primary_key=True)

    type = Column(String(16), nullable=False, doc='Type of the stored publication e.g. PubMed')
    reference = Column(String(255), nullable=False, doc='Reference identifier of the publication e.g. PubMed_ID')
    sha512 = Column(String(255), index=True)

    name = Column(String(255), nullable=True, doc='Journal name')
    title = Column(Text, nullable=True, doc='Title of the publication')
    volume = Column(Text, nullable=True, doc='Volume of the journal')
    issue = Column(Text, nullable=True, doc='Issue within the volume')
    pages = Column(Text, nullable=True, doc='Pages of the publication')
    date = Column(Date, nullable=True, doc='Publication date')

    first_id = Column(Integer, ForeignKey(f'{AUTHOR_TABLE_NAME}.id'), nullable=True, doc='First author')
    first = relationship("Author", foreign_keys=[first_id])

    last_id = Column(Integer, ForeignKey(f'{AUTHOR_TABLE_NAME}.id'), nullable=True, doc='Last author')
    last = relationship("Author", foreign_keys=[last_id])

    authors = relationship("Author", secondary=author_citation, backref='citations')

    __table_args__ = (
        UniqueConstraint(type, reference),
    )

    def __str__(self):
        return '{}:{}'.format(self.type, self.reference)

    @property
    def is_pubmed(self):
        """Returns if this is a PubMed citation

        :rtype:
        """
        return CITATION_TYPE_PUBMED == self.type

    @property
    def is_enriched(self):
        """Returns if this citation has been enriched for name, title, etc.

        :rtype:
        """
        return self.title is not None and self.name is not None

    @staticmethod
    def filter_reference(reference_type, reference):
        return Citation.sha512 == hash_citation(type=reference_type, reference=reference)

    def to_json(self, include_id=False):
        """Creates a citation dictionary that is used to recreate the edge data dictionary of a :class:`BELGraph`.

        :param bool include_id: If true, includes the model identifier
        :return: Citation dictionary for the recreation of a :class:`BELGraph`.
        :rtype: dict[str,str]
        """
        result = {
            CITATION_REFERENCE: self.reference,
            CITATION_TYPE: self.type
        }

        if include_id:
            result['id'] = self.id

        if self.name:
            result[CITATION_NAME] = self.name

        if self.title:
            result[CITATION_TITLE] = self.title

        if self.volume:
            result[CITATION_VOLUME] = self.volume

        if self.pages:
            result[CITATION_PAGES] = self.pages

        if self.date:
            result[CITATION_DATE] = self.date.strftime('%Y-%m-%d')

        if self.first:
            result[CITATION_FIRST_AUTHOR] = self.first.name

        if self.last:
            result[CITATION_LAST_AUTHOR] = self.last.name

        if self.authors:
            result[CITATION_AUTHORS] = sorted(
                author.name
                for author in self.authors
            )

        return result


class Evidence(Base):
    """This table contains the evidence text that proves a specific relationship and refers the source that is cited."""
    __tablename__ = EVIDENCE_TABLE_NAME

    id = Column(Integer, primary_key=True)

    text = Column(Text, nullable=False, doc='Supporting text from a given publication')
    sha512 = Column(String(255), index=True, nullable=False)

    citation_id = Column(Integer, ForeignKey('{}.id'.format(CITATION_TABLE_NAME)), nullable=False)
    citation = relationship('Citation', backref=backref('evidences'))

    def __str__(self):
        return '{}:{}'.format(self.citation, self.text)

    def to_json(self, include_id=False):
        """Creates a dictionary that is used to recreate the edge data dictionary for a :class:`BELGraph`.

        :param bool include_id: If true, includes the model identifier
        :return: Dictionary containing citation and evidence for a :class:`BELGraph` edge.
        :rtype: dict
        """
        result = {
            CITATION: self.citation.to_json(),
            EVIDENCE: self.text
        }

        if include_id:
            result['id'] = self.id

        return result


edge_annotation = Table(
    EDGE_ANNOTATION_TABLE_NAME, Base.metadata,
    Column('edge_id', Integer, ForeignKey('{}.id'.format(EDGE_TABLE_NAME)), primary_key=True),
    Column('annotationEntry_id', Integer, ForeignKey('{}.id'.format(ANNOTATION_ENTRY_TABLE_NAME)), primary_key=True)
)

edge_property = Table(
    EDGE_PROPERTY_TABLE_NAME, Base.metadata,
    Column('edge_id', Integer, ForeignKey('{}.id'.format(EDGE_TABLE_NAME)), primary_key=True),
    Column('property_id', Integer, ForeignKey('{}.id'.format(PROPERTY_TABLE_NAME)), primary_key=True)
)


class Edge(Base):
    """Relationships are represented in this table. It shows the nodes that are in a relation to eachother and provides
    information about the context of the relation by refaring to the annotation, property and evidence tables.
    """
    __tablename__ = EDGE_TABLE_NAME

    id = Column(Integer, primary_key=True)

    bel = Column(Text, nullable=False, doc='Valid BEL statement that represents the given edge')
    relation = Column(String(255), nullable=False, index=True)

    source_id = Column(Integer, ForeignKey(f'{NODE_TABLE_NAME}.id'), nullable=False)
    source = relationship('Node', foreign_keys=[source_id],
                          backref=backref('out_edges', lazy='dynamic', cascade='all, delete-orphan'))

    target_id = Column(Integer, ForeignKey(f'{NODE_TABLE_NAME}.id'), nullable=False)
    target = relationship('Node', foreign_keys=[target_id],
                          backref=backref('in_edges', lazy='dynamic', cascade='all, delete-orphan'))

    evidence_id = Column(Integer, ForeignKey('{}.id'.format(EVIDENCE_TABLE_NAME)), nullable=True)
    evidence = relationship("Evidence", backref=backref('edges', lazy='dynamic'))

    annotations = relationship('AnnotationEntry', secondary=edge_annotation, lazy="dynamic",
                               backref=backref('edges', lazy='dynamic'))
    properties = relationship('Property', secondary=edge_property, lazy="dynamic")  # , cascade='all, delete-orphan')

    sha512 = Column(String(255), index=True, doc='The hash of the source, target, and associated metadata')

    def __str__(self):
        return self.bel

    def __repr__(self):
        return '<Edge {}: {}>'.format(self.sha512[:10], self.bel)

    def get_annotations_json(self):
        """Formats the annotations properly

        :rtype: Optional[dict[str,dict[str,bool]]
        """
        annotations = {}

        for entry in self.annotations:
            if entry.annotation.keyword not in annotations:
                annotations[entry.annotation.keyword] = {entry.name: True}
            else:
                annotations[entry.annotation.keyword][entry.name] = True

        return annotations or None

    def get_data_json(self):
        """Gets the PyBEL edge data dictionary this edge represents

        :rtype: dict
        """
        data = {
            RELATION: self.relation,
        }

        annotations = self.get_annotations_json()
        if annotations:
            data[ANNOTATIONS] = annotations

        if self.evidence:
            data.update(self.evidence.to_json())

        for prop in self.properties:  # FIXME this is also probably broken for translocations or mixed activity/degrad
            if prop.side not in data:
                data[prop.side] = prop.to_json()
            else:
                data[prop.side].update(prop.to_json())

        return data

    def to_json(self, include_id=False, include_hash=False):
        """Creates a dictionary of one BEL Edge that can be used to create an edge in a :class:`BELGraph`.

        :param bool include_id: Include the database identifier?
        :param bool include_hash: Include the node hash?
        :return: Dictionary that contains information about an edge of a :class:`BELGraph`. Including participants
                 and edge data information.
        :rtype: dict
        """
        result = {
            'source': self.source.as_bel(),
            'target': self.target.as_bel(),
            'data': self.get_data_json(),
        }

        if include_id:
            result['id'] = self.id

        if include_hash:
            result[HASH] = self.sha512

        return result

    def insert_into_graph(self, graph):
        """Inserts this edge into a BEL Graph

        :param pybel.BELGraph graph: A BEL graph
        :return: The hash of the edge added
        :rtype: str
        """
        u = self.source.as_bel()
        v = self.target.as_bel()
        data = self.get_data_json()
        key = hash_edge(u, v, data)
        graph.add_edge(u, v, key=key, **data)
        return key


class Property(Base):
    """The property table contains additional information that is used to describe the context of a relation."""
    __tablename__ = PROPERTY_TABLE_NAME

    id = Column(Integer, primary_key=True)

    is_subject = Column(Boolean, doc='Identifies which participant of the edge if affected by the given property')
    modifier = Column(String(255), doc='The modifier: one of activity, degradation, location, or translocation')

    relative_key = Column(String(255), nullable=True, doc='Relative key of effect e.g. to_tloc or from_tloc')

    sha512 = Column(String(255), index=True)

    effect_id = Column(Integer, ForeignKey(f'{NAMESPACE_ENTRY_TABLE_NAME}.id'), nullable=True)
    effect = relationship('NamespaceEntry')

    @property
    def side(self):
        """Returns either :data:`pybel.constants.SUBJECT` or :data:`pybel.constants.OBJECT`

        :rtype: str
        """
        return SUBJECT if self.is_subject else OBJECT

    def to_json(self):
        """Creates a property dict that is used to recreate an edge dictionary for a :class:`BELGraph`.

        :return: Property dictionary of an edge that is participant (sub/obj) related.
        :rtype: dict
        """
        participant = self.side

        prop_dict = {
            participant: {
                MODIFIER: self.modifier  # FIXME this is probably wrong for location
            }
        }

        if self.modifier == LOCATION:
            prop_dict[participant] = {
                LOCATION: self.effect.to_json()
            }
        if self.relative_key:  # for translocations
            prop_dict[participant][EFFECT] = {
                self.relative_key: self.effect.to_json()
            }
        elif self.effect:  # for activities
            prop_dict[participant][EFFECT] = self.effect.to_json()

        # degradations don't have modifications

        return prop_dict
