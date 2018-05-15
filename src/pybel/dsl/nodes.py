# -*- coding: utf-8 -*-

import abc
import hashlib
import pickle
from functools import total_ordering

from .exc import InferCentralDogmaException
from .utils import entity
from ..constants import *
from ..utils import ensure_quotes

__all__ = [
    'BaseEntity',
    'BaseAbundance',
    'Variant',
    'abundance',
    'gene',
    'rna',
    'mirna',
    'protein',
    'complex_abundance',
    'composite_abundance',
    'bioprocess',
    'pathology',
    'named_complex_abundance',
    'reaction',
    'pmod',
    'gmod',
    'hgvs',
    'hgvs_reference',
    'hgvs_unspecified',
    'gene_substitution',
    'protein_substitution',
    'protein_deletion',
    'fragment',
    'fusion_range',
    'missing_fusion_range',
    'protein_fusion',
    'rna_fusion',
    'gene_fusion',
]


def _sort_by_bel(elements):
    return sorted(elements, key=lambda e: e._get_bel())


def _digest_bytes(b):
    tuple_sha512 = hashlib.sha512(b)
    return tuple_sha512.hexdigest()[:8]


def _hash_tuple(t):
    tuple_bytes = pickle.dumps(t)
    return _digest_bytes(tuple_bytes)


def _hash_node(node_tuple):
    """Converts a PyBEL node tuple to a hash

    :param tuple node_tuple: A BEL node
    :return: A hashed version of the node tuple using :func:`hashlib.sha512` hash of the binary pickle dump
    :rtype: str
    """
    return _hash_tuple(node_tuple)


class BelConvertable(dict, metaclass=abc.ABCMeta):
    def __init__(self, *args, **kwargs):
        super(BelConvertable, self).__init__(*args, **kwargs)
        self._bel = None

    @abc.abstractmethod
    def as_tuple(self):
        """Returns this entity as a canonical tuple

        :rtype: tuple
        """

    @abc.abstractmethod
    def as_bel(self):
        """Returns this entity as canonical BEL

        :rtype: str
        """

    def as_sha512(self):
        """Returns this entity as a hash

        :rtype: str
        """
        return _hash_node(self.as_tuple())

    def _get_bel(self):
        if self._bel is None:
            self._bel = self.as_bel()

        return self._bel

    def __hash__(self):
        """Use the tuple serialization of this node as the hash

        :rtype: int
        """
        return hash(self._get_bel())

    def __str__(self):
        return self._get_bel()


@total_ordering
class BaseEntity(BelConvertable):
    """This class represents all BEL nodes. It can be converted to a tuple and hashed."""

    def __init__(self, func, **kwargs):
        """
        :param str func: The PyBEL function
        """
        super(BaseEntity, self).__init__(**{FUNCTION: func})

    @property
    def function(self):
        """Returns the function

        :rtype: str
        """
        return self[FUNCTION]

    @property
    def short_bel_function(self):
        """Returns the short version of the function, appropriate for serialization to BEL

        :rtype: str
        """
        return rev_abundance_labels[self.function]

    def __hash__(self):
        """Use the tuple serialization of this node as the hash

        :rtype: int
        """
        return hash(self.as_tuple())

    def __eq__(self, other):
        if not isinstance(other, BelConvertable):
            raise TypeError('wrong type: {} {}'.format(other.__class__.__name__, other))

        return self._get_bel() == other._get_bel()

    def __lt__(self, other):
        """
        :param BaseEntity other:
        :rtype: bool
        """
        if not isinstance(other, BelConvertable):
            raise TypeError('wrong type: {} {}'.format(other.__class__.__name__, other))

        return self._get_bel() < other._get_bel()


class BaseAbundance(BaseEntity):
    """The superclass for building node data dictionaries"""

    def __init__(self, func, namespace, name=None, identifier=None):
        """Raises an exception if name and identifier are both None.

        :param str func: The PyBEL function
        :param str namespace: The name of the namespace
        :param Optional[str] name:
        :param Optional[str] identifier:
        :raises: PyBELDSLException
        """
        super(BaseAbundance, self).__init__(func=func)
        self.update(entity(namespace=namespace, name=name, identifier=identifier))

    @property
    def namespace(self):
        """Returns the namespace

        :rtype: str
        """
        return self[NAMESPACE]

    @property
    def name(self):
        """Returns the name, if defined

        :rtype: Optional[str]
        """
        return self.get(NAME)

    @property
    def identifier(self):
        """Returns the identifier, if defined

        :rtype: Optional[str]
        """
        return self.get(IDENTIFIER)

    @property
    def _prioritized_name(self):
        """
        :rtype: str
        """
        return self.name if self.name else self.identifier

    def _get_ns_arg(self):
        """Gets the qualified name for this entity

        Usually will be NS:NAME or NS:IDENTIFIER, but there will be a proposal to BEL to allow for the namespace and
        identifier to be disambiguated with a hashtag, so one or both can be written like NS[:NAME][#IDENTIFIER]
        """
        return '{}:{}'.format(self.namespace, ensure_quotes(self._prioritized_name))

    def as_tuple(self):
        """Returns this node as a PyBEL node tuple

        :rtype: tuple
        """
        return self.function, self.namespace, self._prioritized_name

    def as_bel(self):
        """Returns this node as BEL

        :rtype: str
        """
        return "{}({})".format(
            self.short_bel_function,
            self._get_ns_arg()
        )


class abundance(BaseAbundance):
    """Builds an abundance node data dictionary"""

    def __init__(self, namespace, name=None, identifier=None):
        """
        :param str namespace: The name of the database used to identify this entity
        :param str name: The database's preferred name or label for this entity
        :param str identifier: The database's identifier for this entity

        Example:

        >>> bioprocess(namespace='CHEBI', name='water')
        """
        super(abundance, self).__init__(ABUNDANCE, namespace=namespace, name=name, identifier=identifier)


class bioprocess(BaseAbundance):
    """Builds a biological process node data dictionary"""

    def __init__(self, namespace, name=None, identifier=None):
        """
        :param str namespace: The name of the database used to identify this entity
        :param str name: The database's preferred name or label for this entity
        :param str identifier: The database's identifier for this entity

        Example:

        >>> bioprocess(namespace='GO', name='apoptosis')
        """
        super(bioprocess, self).__init__(BIOPROCESS, namespace=namespace, name=name, identifier=identifier)


class pathology(BaseAbundance):
    """Builds a pathology node data dictionary"""

    def __init__(self, namespace, name=None, identifier=None):
        """
        :param str namespace: The name of the database used to identify this entity
        :param str name: The database's preferred name or label for this entity
        :param str identifier: The database's identifier for this entity

        Example:

        >>> pathology(namespace='DO', name='Alzheimer Disease')
        """
        super(pathology, self).__init__(PATHOLOGY, namespace=namespace, name=name, identifier=identifier)


class CentralDogma(BaseAbundance):
    """Builds a central dogma (gene, mirna, rna, protein) node data dictionary"""

    def __init__(self, func, namespace, name=None, identifier=None, variants=None):
        """
        :param str func: The PyBEL function to use
        :param str namespace: The name of the database used to identify this entity
        :param str name: The database's preferred name or label for this entity
        :param str identifier: The database's identifier for this entity
        :param Optional[Variant or list[Variant]] variants: An optional variant or list of variants
        """
        super(CentralDogma, self).__init__(func, namespace, name=name, identifier=identifier)

        if isinstance(variants, Variant):
            self[VARIANTS] = [variants]
        elif variants is not None:
            self[VARIANTS] = _sort_by_bel(variants)

    @property
    def variants(self):
        """Returns the variants, if they exist

        :rtype: Optional[list[Variant]]
        """
        return self.get(VARIANTS)

    def as_tuple(self):
        """Converts to a tuple"""
        t = super(CentralDogma, self).as_tuple()

        if self.variants:
            return t + _tuplable_list_as_tuple(self.variants)

        return t

    def as_bel(self):
        """Converts to BEL

        :rtype: str
        """
        if not self.variants:
            return super(CentralDogma, self).as_bel()

        variants_bel = (variant.as_bel() for variant in self.variants)

        return '{}({}, {})'.format(
            self.short_bel_function,
            self._get_ns_arg(),
            ', '.join(variants_bel)
        )

    def get_parent(self):
        """Gets the parent, or none if it's already a reference node

        :rtype: Optional[CentralDogma]

        Example usage:

        >>> ab42 = protein(name='APP', namespace='HGNC', variants=[fragment(start=672, stop=713)])
        >>> app = ab42.get_parent()
        >>> assert 'p(HGNC:APP)' == app.as_bel()
        """
        if VARIANTS not in self:
            return

        return self.__class__(namespace=self.namespace, name=self.name, identifier=self.identifier)

    def with_variants(self, variants):
        """Creates a new entity with the given variants

        :param Variant or list[Variant] variants: An optional variant or list of variants
        :rtype: CentralDogma

        Example Usage:

        >>> app = protein(name='APP', namespace='HGNC')
        >>> ab42 = app.with_variants([fragment(start=672, stop=713)])
        >>> assert 'p(HGNC:APP, frag(672_713)' == ab42.as_bel()
        """
        return self.__class__(
            namespace=self.namespace,
            name=self.name,
            identifier=self.identifier,
            variants=variants
        )


class Variant(BelConvertable):
    """The superclass for variant dictionaries"""

    def __init__(self, kind):
        """
        :param str kind: The kind of variant
        """
        super(Variant, self).__init__({KIND: kind})

    def __lt__(self, other):
        """
        :param BaseEntity other:
        :rtype: bool
        """
        if not isinstance(other, Variant):
            raise TypeError('wrong type: {} {}'.format(other.__class__.__name__, other))

        return self._get_bel() < other._get_bel()


class pmod(Variant):
    """Builds a protein modification variant dictionary"""

    def __init__(self, name, code=None, position=None, namespace=None, identifier=None):
        """
        :param str name: The name of the modification
        :param str code: The three letter amino acid code for the affected residue. Capital first letter.
        :param int position: The position of the affected residue
        :param str namespace: The namespace to which the name of this modification belongs
        :param str identifier: The identifier of the name of the modification

        Either the name or the identifier must be used. If the namespace is omitted, it is assumed that a name is
        specified from the BEL default namespace.

        Example from BEL default namespace:

        >>> pmod('Ph', code='Thr', position=308)

        Example from custom namespace:

        >>> pmod(name='protein phosphorylation', namespace='GO', code='Thr', position=308)

        Example from custom namespace additionally qualified with identifier:

        >>> pmod(name='protein phosphorylation', namespace='GO', identifier='GO:0006468', code='Thr', position=308)
        """
        super(pmod, self).__init__(PMOD)

        self[IDENTIFIER] = entity(
            namespace=(namespace or BEL_DEFAULT_NAMESPACE),
            name=name,
            identifier=identifier
        )

        if code:
            self[PMOD_CODE] = code

        if position:
            self[PMOD_POSITION] = position

    def as_tuple(self):
        """Convert a pmod to tuple

        :rtype: tuple
        """
        identifier = self[IDENTIFIER][NAMESPACE], self[IDENTIFIER][NAME]
        params = tuple(self[key] for key in PMOD_ORDER[2:] if key in self)
        return (PMOD,) + (identifier,) + params

    def as_bel(self):
        """Convert a pmod to BEL

        :rtype: str
        """
        return 'pmod({}{})'.format(
            str(self[IDENTIFIER]),
            ''.join(', {}'.format(self[x]) for x in PMOD_ORDER[2:] if x in self)
        )


class gmod(Variant):
    """Builds a gene modification variant dictionary"""

    def __init__(self, name, namespace=None, identifier=None):
        """
        :param str name: The name of the gene modification
        :param Optional[str] namespace: The namespace of the gene modification
        :param Optional[str] identifier: The identifier of the name in the database

        Either the name or the identifier must be used. If the namespace is omitted, it is assumed that a name is
        specified from the BEL default namespace.

        Example from BEL default namespace:

        >>> gmod(name='Me')

        Example from custom namespace:

        >>> gmod(name='DNA methylation', namespace='GO', identifier='GO:0006306',)
        """
        super(gmod, self).__init__(GMOD)

        self[IDENTIFIER] = entity(
            namespace=(namespace or BEL_DEFAULT_NAMESPACE),
            name=name,
            identifier=identifier
        )

    def as_tuple(self):
        """Converts a gmod to a tuple

        :rtype: tuple
        """
        identifier = self[IDENTIFIER][NAMESPACE], self[IDENTIFIER][NAME]
        params = tuple(self[key] for key in GMOD_ORDER[2:] if key in self)
        return (GMOD,) + (identifier,) + params

    def as_bel(self):
        """Returns the gene modification as BEL

        :rtype: str
        """
        return 'gmod({})'.format(str(self[IDENTIFIER]))


class hgvs(Variant):
    """Builds a HGVS variant dictionary"""

    def __init__(self, variant):
        """
        :param str variant: The HGVS variant string

        Example:

        >>> protein(namespace='HGNC', name='AKT1', variants=[hgvs('p.Ala127Tyr')])
        """
        super(hgvs, self).__init__(HGVS)

        self[IDENTIFIER] = variant

    def as_tuple(self):
        """Converts the HGVS variant to a tuple

        :rtype: tuple
        """
        return self[KIND], self[IDENTIFIER]

    def as_bel(self):
        """Returns the HGVS variant as BEL

        :rtype: str
        """
        return 'var({})'.format(self[IDENTIFIER])


class hgvs_reference(hgvs):
    """Represents the "reference" variant in HGVS"""

    def __init__(self):
        super(hgvs_reference, self).__init__('=')


class hgvs_unspecified(hgvs):
    """Represents an unspecified variant in HGVS"""

    def __init__(self):
        super(hgvs_unspecified, self).__init__('?')


class gene_substitution(hgvs):
    """Builds a HGVS variant dictionary for the given protein substitution"""

    def __init__(self, from_nucleotide, position, to_nucleotide):
        """
        :param str from_nucleotide: The old nucleotide letter (ACTG) (capitalized)
        :param int position: The position of the residue
        :param str to_nucleotide: The new nucleotide letter (ACTG) (capitalized)

        Example:
        >>> gene(namespace='HGNC', name='APOE', variants=[gene_substitution('C', 526, 'T')])
        """
        super(gene_substitution, self).__init__('c.{}{}>{}'.format(position, from_nucleotide, to_nucleotide))


class protein_substitution(hgvs):
    """Builds a HGVS variant dictionary for the given protein substitution"""

    def __init__(self, from_aa, position, to_aa):
        """
        :param str from_aa: The 3-letter amino acid code of the original residue
        :param int position: The position of the residue
        :param str to_aa: The 3-letter amino acid code of the new residue

        Example:

        >>> protein(namespace='HGNC', name='AKT1', variants=[protein_substitution('Ala', 127, 'Tyr')])
        """
        super(protein_substitution, self).__init__('p.{}{}{}'.format(from_aa, position, to_aa))


class protein_deletion(hgvs):
    """Builds a HGVS variant dictionary for the given protein deletion"""

    def __init__(self, aa, position):
        """
        :param str aa: The 3-letter amino acid code of the original residue
        :param int position: The position of the residue

        Example:

        >>> protein(namespace='HGNC', name='AKT1', variants=[protein_deletion('Ala', 127)])
        """
        super(protein_deletion, self).__init__('p.{}{}del'.format(aa, position))


class fragment(Variant):
    def __init__(self, start=None, stop=None, description=None):
        """Make a protein fragment dictionary

        :param Optional[int or str] start: The starting position
        :param Optional[int or str] stop: The stopping position
        :param Optional[str] description: An optional description

        Example of specified fragment:

        >>> protein(name='APP', namespace='HGNC', variants=[fragment(start=672, stop=713)])

        Example of unspecified fragment:

        >>> protein(name='APP', namespace='HGNC', variants=[fragment()])

        Example of a fragment specified by its description:

        >>> protein(name='APP', namespace='HGNC', variants=[fragment(description='55kD')])
        """
        super(fragment, self).__init__(FRAGMENT)

        if start is None and stop is None:
            self[FRAGMENT_MISSING] = '?'
        elif start and stop:
            self[FRAGMENT_START] = start
            self[FRAGMENT_STOP] = stop
        elif start:
            self[FRAGMENT_START] = start
            self[FRAGMENT_STOP] = '?'
        else:  # just stop
            self[FRAGMENT_START] = '?'
            self[FRAGMENT_STOP] = stop

        if description:
            self[FRAGMENT_DESCRIPTION] = description

    def as_tuple(self):
        """Converts the fragment to a tuple

        :rtype: tuple
        """
        if FRAGMENT_MISSING in self:
            result = FRAGMENT, '?'
        else:
            result = FRAGMENT, (self[FRAGMENT_START], self[FRAGMENT_STOP])

        if FRAGMENT_DESCRIPTION in self:
            return result + (self[FRAGMENT_DESCRIPTION],)

        return result

    def as_bel(self):
        """Returns the fragment as BEL

        :rtype: str
        """
        if FRAGMENT_MISSING in self:
            res = '?'
        else:
            res = '{}_{}'.format(self[FRAGMENT_START], self[FRAGMENT_STOP])

        if FRAGMENT_DESCRIPTION in self:
            res += ', "{}"'.format(self[FRAGMENT_DESCRIPTION])

        return 'frag({})'.format(res)


class gene(CentralDogma):
    """Builds a gene node data dictionary"""

    def __init__(self, namespace, name=None, identifier=None, variants=None):
        """
        :param str namespace: The name of the database used to identify this entity
        :param str name: The database's preferred name or label for this entity
        :param str identifier: The database's identifier for this entity
        :param Optional[Variant or list[Variant]] variants: An optional variant or list of variants
        """
        super(gene, self).__init__(GENE, namespace, name=name, identifier=identifier, variants=variants)


class _Transcribable(CentralDogma):
    """An intermediate class between the CentralDogma and rna/mirna because both of them share the ability to
    get their corresponding gene"""

    def get_gene(self):
        """Gets the corresponding gene. Raises an exception if it's not the reference node

        :rtype: gene
        :raises: InferCentralDogmaException
        """
        if self.variants:
            raise InferCentralDogmaException('can not get gene for variant')

        return gene(
            namespace=self.namespace,
            name=self.name,
            identifier=self.identifier
        )


class rna(_Transcribable):
    """Builds an RNA node data dictionary"""

    def __init__(self, namespace, name=None, identifier=None, variants=None):
        """
        :param str namespace: The name of the database used to identify this entity
        :param str name: The database's preferred name or label for this entity
        :param str identifier: The database's identifier for this entity
        :param Optional[Variant or list[Variant]] variants: An optional variant or list of variants


        Example: AKT1 protein coding gene's RNA:

        >>> rna(namespace='HGNC', name='AKT1', identifier='391')

        Non-coding RNA's can also be encoded such as `U85 <https://www-snorna.biotoul.fr/plus.php?id=U85>`_:

        >>> rna(namespace='SNORNABASE', identifer='SR0000073')
        """
        super(rna, self).__init__(RNA, namespace, name=name, identifier=identifier, variants=variants)


class mirna(_Transcribable):
    """Builds a miRNA node data dictionary"""

    def __init__(self, namespace, name=None, identifier=None, variants=None):
        """
        :param str namespace: The name of the database used to identify this entity
        :param str name: The database's preferred name or label for this entity
        :param str identifier: The database's identifier for this entity
        :param Optional[Variant or list[Variant]] variants: An optional variant or list of variants

        Human miRNA's are listed on HUGO's `MicroRNAs (MIR) <https://www.genenames.org/cgi-bin/genefamilies/set/476>`_
        gene family.

        MIR1-1 from `HGNC <https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=31499>`_:

        >>> mirna(namespace='HGNC', name='MIR1-1', identifier='31499')

        MIR1-1 from `miRBase <http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0000651>`_:

        >>> mirna(namespace='MIRBASE', identifier='MI0000651')

        MIR1-1 from `Entrez Gene <https://view.ncbi.nlm.nih.gov/gene/406904>`_

        >>> mirna(namespace='ENTREZ', identifier='406904')
        """
        super(mirna, self).__init__(MIRNA, namespace, name=name, identifier=identifier, variants=variants)


class protein(CentralDogma):
    """Builds a protein node data dictionary"""

    def __init__(self, namespace, name=None, identifier=None, variants=None):
        """Returns the node data dictionary for a protein

        :param str namespace: The name of the database used to identify this entity
        :param str name: The database's preferred name or label for this entity
        :param str identifier: The database's identifier for this entity
        :param Optional[Variant or list[Variant]] variants: An optional variant or list of variants

        Example: AKT

        >>> protein(namespace='HGNC', name='AKT1')

        Example: AKT with optionally included HGNC database identifier

        >>> protein(namespace='HGNC', name='AKT1', identifier='391')

        Example: AKT with phosphorylation

        >>> protein(namespace='HGNC', name='AKT', variants=[pmod('Ph', code='Thr', position=308)])
        """
        super(protein, self).__init__(PROTEIN, namespace, name=name, identifier=identifier, variants=variants)

    def get_rna(self):
        """Gets the corresponding RNA. Raises an exception if it's not the reference node

        :rtype: rna
        :raises: InferCentralDogmaException
        """
        if self.variants:
            raise InferCentralDogmaException('can not get rna for variant')

        return rna(
            namespace=self.namespace,
            name=self.name,
            identifier=self.identifier
        )


def _tuplable_list_as_tuple(entities):
    """A helper function for converting reaction list

    :type entities: iter[BaseEntity]
    :rtype: tuple[tuple]
    """
    return tuple(e.as_tuple() for e in entities)


def _entity_list_as_bel(entities):
    """A helper function for stringifying a list of BEL entities

    :type entities: iter[BaseEntity]
    :rtype: str
    """
    return ', '.join(e.as_bel() for e in entities)


class reaction(BaseEntity):
    """Builds a reaction node data dictionary"""

    def __init__(self, reactants, products):
        """
        :param BaseAbundance or iter[BaseAbundance] reactants: An entity or list of entities participating as reactants
        :param BaseAbundance or iter[BaseAbundance] products: An entity or list of entities participating as products

        Example:

        >>> reaction([protein(namespace='HGNC', name='KNG1')], [abundance(namespace='CHEBI', name='bradykinin')])
        """
        super(reaction, self).__init__(func=REACTION)

        if isinstance(reactants, BaseAbundance):
            self[REACTANTS] = [reactants]
        else:
            self[REACTANTS] = _sort_by_bel(reactants)

        if isinstance(products, BaseAbundance):
            self[PRODUCTS] = [products]
        else:
            self[PRODUCTS] = _sort_by_bel(products)

    def as_tuple(self):
        """Returns the reaction as a canonicalized tuple

        :rtype: tuple
        """
        return (
            self.function,
            _tuplable_list_as_tuple(self[REACTANTS]),
            _tuplable_list_as_tuple(self[PRODUCTS])
        )

    def as_bel(self):
        """Returns the reaction as canonicalized BEL

        :rtype: str
        """
        return 'rxn(reactants({}), products({}))'.format(
            _entity_list_as_bel(self[REACTANTS]),
            _entity_list_as_bel(self[PRODUCTS])
        )


class ListAbundance(BaseEntity):
    """The superclass for building list abundance (complex, abundance) node data dictionaries"""

    def __init__(self, func, members):
        """
        :param str func: The PyBEL function
        :param list[BaseAbundance] members: A list of PyBEL node data dictionaries
        """
        super(ListAbundance, self).__init__(func=func)
        self[MEMBERS] = _sort_by_bel(members)

    def as_tuple(self):
        """Returns the list abundance as a canonicalized tuple

        :rtype: tuple
        """
        return (self.function,) + _tuplable_list_as_tuple(self[MEMBERS])

    def as_bel(self):
        """Returns the list abundance as canonicalized BEL

        :rtype: str
        """
        return '{}({})'.format(
            self.short_bel_function,
            _entity_list_as_bel(self[MEMBERS])
        )


class complex_abundance(ListAbundance):
    """Builds a complex abundance node data dictionary with the optional ability to specificy a name"""

    def __init__(self, members, namespace=None, name=None, identifier=None):
        """
        :param list[BaseAbundance] members: A list of PyBEL node data dictionaries
        :param Optional[str] namespace: The namespace from which the name originates
        :param Optional[str] name: The name of the complex
        :param Optional[str] identifier: The identifier in the namespace in which the name originates
        """
        super(complex_abundance, self).__init__(func=COMPLEX, members=members)

        if namespace:
            self.update(entity(namespace=namespace, name=name, identifier=identifier))


class named_complex_abundance(BaseAbundance):
    """Builds a named complex abundance node data dictionary"""

    def __init__(self, namespace=None, name=None, identifier=None):
        """
        :param str namespace: The name of the database used to identify this entity
        :param str name: The database's preferred name or label for this entity
        :param str identifier: The database's identifier for this entity

        Example:

        >>> named_complex_abundance(namespace='SCOMP', name='Calcineurin Complex')
        """
        super(named_complex_abundance, self).__init__(
            func=COMPLEX,
            namespace=namespace,
            name=name,
            identifier=identifier
        )


class composite_abundance(ListAbundance):
    """Builds a composite abundance node data dictionary"""

    def __init__(self, members):
        """
        :param list[BaseAbundance] members: A list of PyBEL node data dictionaries
        """
        super(composite_abundance, self).__init__(func=COMPOSITE, members=members)


class FusionRangeBase(BelConvertable):
    """The superclass for fusion range data dictionaries"""


class missing_fusion_range(FusionRangeBase):
    """Builds a missing fusion range data dictionary"""

    def __init__(self):
        super(missing_fusion_range, self).__init__({
            FUSION_MISSING: '?'
        })

    def __str__(self):
        return self.as_bel()

    def as_tuple(self):
        """
        :rtype: tuple
        """
        return self[FUSION_MISSING],

    def as_bel(self):
        return '?'


class fusion_range(FusionRangeBase):
    """Creates a fusion range data dictionary"""

    def __init__(self, reference, start, stop):
        """
        :param str reference: The reference code
        :param int or str start: The start position, either specified by its integer position, or '?'
        :param int or str stop: The stop position, either specified by its integer position, '?', or '*

        Example fully specified RNA fusion range:

        >>> fusion_range('r', 1, 79)

        """
        super(fusion_range, self).__init__({
            FUSION_REFERENCE: reference,
            FUSION_START: start,
            FUSION_STOP: stop
        })

    def as_tuple(self):
        """Returns the fusion range as a tuple

        :rtype: tuple
        """
        return (
            self[FUSION_REFERENCE],
            self[FUSION_START],
            self[FUSION_STOP]
        )

    def as_bel(self):
        """Returns the fusion range as BEL

        :rtype: str
        """
        return '{reference}.{start}_{stop}'.format(
            reference=self[FUSION_REFERENCE],
            start=self[FUSION_START],
            stop=self[FUSION_STOP]
        )

    def __str__(self):
        return self.as_bel()


class FusionBase(BaseEntity):
    """The superclass for building fusion node data dictionaries"""

    def __init__(self, func, partner5p, partner3p, range5p=None, range3p=None):
        """
        :param str func: A PyBEL function
        :param CentralDogma partner5p: A PyBEL node data dictionary for the 5-prime partner
        :param CentralDogma partner3p: A PyBEL node data dictionary for the 3-prime partner
        :param Optional[FusionRangeBase] range5p: A fusion range for the 5-prime partner
        :param Optional[FusionRangeBase] range3p: A fusion range for the 3-prime partner
        """
        super(FusionBase, self).__init__(func=func)
        self[FUSION] = {
            PARTNER_5P: partner5p,
            PARTNER_3P: partner3p,
            RANGE_5P: range5p or missing_fusion_range(),
            RANGE_3P: range3p or missing_fusion_range()
        }

    @property
    def partner5p(self):
        """
        :rtype: CentralDogma
        """
        return self[FUSION][PARTNER_5P]

    @property
    def partner3p(self):
        """
        :rtype: CentralDogma
        """
        return self[FUSION][PARTNER_3P]

    @property
    def range5p(self):
        """
        :rtype: FusionRangeBase
        """
        return self[FUSION][RANGE_5P]

    @property
    def range3p(self):
        """
        :rtype: FusionRangeBase
        """
        return self[FUSION][RANGE_3P]

    def as_bel(self):
        """Returns the fusion as BEL

        :rtype: str
        """
        return "{}(fus({}:{}, {}, {}:{}, {}))".format(
            self.short_bel_function,
            self.partner5p[NAMESPACE],
            self.partner5p[NAME],
            str(self.range5p),
            self.partner3p[NAMESPACE],
            self.partner3p[NAME],
            str(self.range3p)
        )

    def as_tuple(self):
        """Returns the fusion as a canonicalized tuple

        :rtype: tuple
        """
        fusion = self[FUSION]

        partner5p = fusion[PARTNER_5P][NAMESPACE], fusion[PARTNER_5P][NAME]
        partner3p = fusion[PARTNER_3P][NAMESPACE], fusion[PARTNER_3P][NAME]
        range5p = self.range5p.as_tuple()
        range3p = self.range3p.as_tuple()

        return (
            self.function,
            partner5p,
            range5p,
            partner3p,
            range3p,
        )


class protein_fusion(FusionBase):
    """Builds a protein fusion data dictionary"""

    def __init__(self, partner5p, partner3p, range5p=None, range3p=None):
        """
        :param pybel.dsl.protein partner5p: A PyBEL node data dictionary for the 5-prime partner
        :param pybel.dsl.protein partner3p: A PyBEL node data dictionary for the 3-prime partner
        :param Optional[FusionRangeBase] range5p: A fusion range for the 5-prime partner
        :param Optional[FusionRangeBase] range3p: A fusion range for the 3-prime partner
        """
        super(protein_fusion, self).__init__(PROTEIN, partner5p=partner5p, range5p=range5p, partner3p=partner3p,
                                             range3p=range3p)


class rna_fusion(FusionBase):
    """Builds an RNA fusion data dictionary"""

    def __init__(self, partner5p, partner3p, range5p=None, range3p=None):
        """
        :param pybel.dsl.rna partner5p: A PyBEL node data dictionary for the 5-prime partner
        :param pybel.dsl.rna partner3p: A PyBEL node data dictionary for the 3-prime partner
        :param Optional[FusionRangeBase] range5p: A fusion range for the 5-prime partner
        :param Optional[FusionRangeBase] range3p: A fusion range for the 3-prime partner

        Example, with fusion ranges using the 'r' qualifier:

        >>> rna_fusion(
        >>> ... partner5p=rna(namespace='HGNC', name='TMPRSS2'),
        >>> ... range5p=fusion_range('r', 1, 79),
        >>> ... partner3p=rna(namespace='HGNC', name='ERG'),
        >>> ... range3p=fusion_range('r', 312, 5034)
        >>> )


        Example with missing fusion ranges:

        >>> rna_fusion(
        >>> ... partner5p=rna(namespace='HGNC', name='TMPRSS2'),
        >>> ... partner3p=rna(namespace='HGNC', name='ERG'),
        >>> )
        """
        super(rna_fusion, self).__init__(RNA, partner5p=partner5p, range5p=range5p, partner3p=partner3p,
                                         range3p=range3p)


class gene_fusion(FusionBase):
    """Builds a gene fusion data dictionary"""

    def __init__(self, partner5p, partner3p, range5p=None, range3p=None):
        """
        :param pybel.dsl.gene partner5p: A PyBEL node data dictionary for the 5-prime partner
        :param pybel.dsl.gene partner3p: A PyBEL node data dictionary for the 3-prime partner
        :param Optional[FusionRangeBase] range5p: A fusion range for the 5-prime partner
        :param Optional[FusionRangeBase] range3p: A fusion range for the 3-prime partner

        Example, using fusion ranges with the 'c' qualifier

        >>> gene_fusion(
        >>> ... partner5p=gene(namespace='HGNC', name='TMPRSS2'),
        >>> ... range5p=fusion_range('c', 1, 79),
        >>> ... partner3p=gene(namespace='HGNC', name='ERG'),
        >>> ... range3p=fusion_range('c', 312, 5034)
        >>> )


        Example with missing fusion ranges:

        >>> gene_fusion(
        >>> ... partner5p=gene(namespace='HGNC', name='TMPRSS2'),
        >>> ... partner3p=gene(namespace='HGNC', name='ERG'),
        >>> )
        """
        super(gene_fusion, self).__init__(GENE, partner5p=partner5p, range5p=range5p, partner3p=partner3p,
                                          range3p=range3p)
