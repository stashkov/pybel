# -*- coding: utf-8 -*-

# Internal edge data keys

#: The key for an internal edge data dictionary for the relation string
RELATION = 'relation'
#: The key for an internal edge data dictionary for the citation dictionary
CITATION = 'citation'
#: The key for an internal edge data dictionary for the evidence string
EVIDENCE = 'evidence'
#: The key for an internal edge data dictionary for the annotations dictionary
ANNOTATIONS = 'annotations'
#: The key for an internal edge data dictionary for the subject modifier dictionary
SUBJECT = 'subject'
#: The key for an internal edge data dictionary for the object modifier dictionary
OBJECT = 'object'
#: The key or an internal edge data dictionary for the line number
LINE = 'line'
#: The key representing the hash of the other
HASH = 'hash'

#: The group of all BEL-provided keys for edge data dictionaries, used for hashing.
PYBEL_EDGE_DATA_KEYS = {
    RELATION,
    CITATION,
    EVIDENCE,
    ANNOTATIONS,
    SUBJECT,
    OBJECT,
}

#: The group of all PyBEL-specific keys for edge data dictionaries, not used for hashing.
PYBEL_EDGE_METADATA_KEYS = {
    LINE,
    HASH,
}

#: The group of all PyBEL annotated keys for edge data dictionaries
PYBEL_EDGE_ALL_KEYS = PYBEL_EDGE_DATA_KEYS | PYBEL_EDGE_METADATA_KEYS
