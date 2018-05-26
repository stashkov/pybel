# -*- coding: utf-8 -*-

"""

CX JSON
~~~~~~~
CX is an aspect-oriented network interchange format encoded in JSON with a format inspired by the JSON-LD encoding of
Resource Description Framework (RDF). It is primarily used by the Network Data Exchange (NDEx) and more recent versions
of Cytoscape.

.. seealso::

    - The NDEx Data Model `Specification <http://www.home.ndexbio.org/data-model/>`_
    - `Cytoscape.js <http://js.cytoscape.org/>`_
    - CX Support for Cytoscape.js on the Cytoscape `App Store <http://apps.cytoscape.org/apps/cxsupport>`_

"""

import json
import logging
from collections import defaultdict

import time

from ..constants import *
from ..dsl.nodes import BaseEntity
from ..struct import BELGraph
from ..tokens import dict_to_entity
from ..utils import expand_dict, flatten_dict

CX_NODE_NAME = 'label'

__all__ = [
    'to_cx',
    'to_cx_jsons',
    'to_cx_file',
    'from_cx',
    'from_cx_jsons',
    'from_cx_file',
    'NDEX_SOURCE_FORMAT',
]

log = logging.getLogger(__name__)

NDEX_SOURCE_FORMAT = "ndex:sourceFormat"

_listy_keys = {PRODUCTS, REACTANTS, MEMBERS}


def _cx_to_dict(list_of_dicts, key_tag='k', value_tag='v'):
    """
    :param list[dict] list_of_dicts:
    :param str key_tag:
    :param str value_tag:
    :rtype: dict
    """
    return {
        d[key_tag]: d[value_tag]
        for d in list_of_dicts
    }


def _calculate_canonical_cx_identifier(node):
    """Calculates the canonical name for a given node. If it is a simple node, uses the namespace:name combination.
    Otherwise, it uses the BEL string.

    :param BaseEntity node: A node
    :return: Appropriate identifier for the node for CX indexing
    :rtype: str
    """
    if node[FUNCTION] == COMPLEX and NAMESPACE in node:
        return '{}:{}'.format(node[NAMESPACE], node[NAME])

    if VARIANTS in node or FUSION in node or node[FUNCTION] in {REACTION, COMPOSITE, COMPLEX}:
        return node.as_bel()

    if VARIANTS not in node and FUSION not in node:  # this is should be a simple node
        return '{}:{}'.format(node[NAMESPACE], node[NAME])

    raise ValueError('Unexpected node data: {}'.format(node))


def node_to_cx(nid, node):
    """Converts an entity to a list of CX attributes

    :param int nid:
    :param BaseEntity node:
    :rtype: list[dict]
    """
    rv = []

    for key, v in node.items():
        if key == VARIANTS:
            for i, el in enumerate(v):
                for a, b in flatten_dict(el).items():
                    rv.append({
                        'po': nid,
                        'n': '{}_{}_{}'.format(key, i, a),
                        'v': b
                    })
        elif key == FUSION:
            for a, b in flatten_dict(v).items():
                rv.append({
                    'po': nid,
                    'n': '{}_{}'.format(key, a),
                    'v': b
                })

        elif key == NAME:
            rv.append({
                'po': nid,
                'n': CX_NODE_NAME,
                'v': v
            })

        elif key in _listy_keys:
            rv.append({
                'po': nid,
                'n': key,
                'v': json.dumps(v)
            })

        else:
            rv.append({
                'po': nid,
                'n': key,
                'v': v
            })

    return rv


def edge_to_cx(key, data):
    """Converts the edge hash (key) and edge data to a list of CX entries

    :param int key: The edge key
    :param dict data: The edge data dictionary
    :rtype: list[dict]
    """
    rv = []

    if EVIDENCE in data:
        rv.append({
            'po': key,
            'n': EVIDENCE,
            'v': data[EVIDENCE]
        })

        for citation_key, value in data[CITATION].items():
            rv.append({
                'po': key,
                'n': '{}_{}'.format(CITATION, citation_key),
                'v': value
            })

    if ANNOTATIONS in data:
        for annotation, values in data[ANNOTATIONS].items():
            rv.append({
                'po': key,
                'n': annotation,
                'v': sorted(values),
                'd': 'list_of_string',
            })

    if SUBJECT in data:
        for subject_key, value in flatten_dict(data[SUBJECT]).items():
            rv.append({
                'po': key,
                'n': '{}_{}'.format(SUBJECT, subject_key),
                'v': value
            })

    if OBJECT in data:
        for object_key, value in flatten_dict(data[OBJECT]).items():
            rv.append({
                'po': key,
                'n': '{}_{}'.format(OBJECT, object_key),
                'v': value
            })

    return rv


def to_cx(graph):
    """Converts BEL Graph to CX data format (as in-memory JSON) for use with `NDEx <http://www.ndexbio.org/>`_

    :param BELGraph graph: A BEL Graph
    :return: The CX JSON for this graph
    :rtype: list

    .. seealso::

        - `NDEx Python Client <https://github.com/ndexbio/ndex-python>`_
        - `PyBEL / NDEx Python Client Wrapper <https://github.com/pybel/pybel2ndex>`_

    """
    nodes_aspect = []
    nodes_attributes_aspect = []

    node_to_id = {}

    for nid, node in enumerate(sorted(graph)):
        node_to_id[node] = nid
        nodes_aspect.append({
            '@id': nid,
            'n': _calculate_canonical_cx_identifier(node)
        })

        nodes_attributes_aspect.extend(node_to_cx(nid, node))

    edges_aspect = []
    edges_attributes_aspect = []

    for eid, (u, v, data) in enumerate(graph.edges(data=True)):
        edges_aspect.append({
            '@id': eid,
            's': node_to_id[u],
            't': node_to_id[v],
            'i': data[RELATION],
        })
        edges_attributes_aspect.extend(edge_to_cx(eid, data))

    context_legend = {}

    for key in graph.namespace_url:
        context_legend[key] = GRAPH_NAMESPACE_URL

    for key in graph.namespace_pattern:
        context_legend[key] = GRAPH_NAMESPACE_PATTERN

    for key in graph.annotation_url:
        context_legend[key] = GRAPH_ANNOTATION_URL

    for key in graph.annotation_pattern:
        context_legend[key] = GRAPH_ANNOTATION_PATTERN

    for key in graph.annotation_list:
        context_legend[key] = GRAPH_ANNOTATION_LIST

    context_legend_aspect = []
    for keyword, resource_type in context_legend.items():
        context_legend_aspect.append({
            'k': keyword,
            'v': resource_type
        })

    annotation_list_keys_lookup = {keyword: i for i, keyword in enumerate(sorted(graph.annotation_list))}
    annotation_lists_aspect = []
    for keyword, values in graph.annotation_list.items():
        for value in values:
            annotation_lists_aspect.append({
                'k': annotation_list_keys_lookup[keyword],
                'v': value
            })

    context_entry_dict = {}
    context_entry_dict.update(graph.namespace_url)
    context_entry_dict.update(graph.namespace_pattern)
    context_entry_dict.update(graph.annotation_url)
    context_entry_dict.update(graph.annotation_pattern)
    context_entry_dict.update(annotation_list_keys_lookup)
    context_entry_dict.update(graph.namespace_url)
    context_aspect = [context_entry_dict]

    network_attributes_aspect = [{
        "n": NDEX_SOURCE_FORMAT,
        "v": "PyBEL"
    }]
    for key, v in graph.document.items():
        network_attributes_aspect.append({
            'n': key,
            'v': v
        })

    # Coalesce to cx
    # cx = create_aspect.number_verification()
    cx = [{'numberVerification': [{'longNumber': 281474976710655}]}]

    cx_pairs = [
        ('@context', context_aspect),
        ('context_legend', context_legend_aspect),
        ('networkAttributes', network_attributes_aspect),
        ('nodes', nodes_aspect),
        ('nodeAttributes', nodes_attributes_aspect),
        ('edges', edges_aspect),
        ('edgeAttributes', edges_attributes_aspect),
    ]

    if annotation_lists_aspect:  # don't include an empty entry!
        cx_pairs.append(('annotation_lists', annotation_lists_aspect))

    cx_metadata = []

    for key, aspect in cx_pairs:
        aspect_dict = {
            "name": key,
            "elementCount": len(aspect),
            "lastUpdate": time.time(),
            "consistencyGroup": 1,
            "properties": [],
            "version": "1.0"
        }

        if key in {'citations', 'supports', 'nodes', 'edges'}:
            aspect_dict['idCounter'] = len(aspect)

        cx_metadata.append(aspect_dict)

    cx.append({
        'metaData': cx_metadata
    })

    for key, aspect in cx_pairs:
        cx.append({
            key: aspect
        })

    cx.append({"status": [{"error": "", "success": True}]})

    return cx


def to_cx_file(graph, file, indent=2, **kwargs):
    """Writes this graph to a JSON file in CX format

    :param BELGraph graph: A BEL graph
    :param file file: A writable file or file-like
    :param Optional[int] indent: How many spaces to use to pretty print. Change to None for no pretty printing

    Example:

    >>> from pybel import from_url, to_cx_file
    >>> from pybel.constants import SMALL_CORPUS_URL
    >>> graph = from_url(SMALL_CORPUS_URL)
    >>> with open('graph.cx', 'w') as f:
    >>> ... to_cx_file(graph, f)
    """
    graph_cx_json_dict = to_cx(graph)
    json.dump(graph_cx_json_dict, file, ensure_ascii=False, indent=indent, **kwargs)


def to_cx_jsons(graph, **kwargs):
    """Dumps a BEL graph as CX JSON to a string
    
    :param BELGraph graph: A BEL Graph
    :return: CX JSON string
    :rtype: str
    """
    graph_cx_json_str = to_cx(graph)
    return json.dumps(graph_cx_json_str, **kwargs)


def _iterate_list_of_dicts(lod):
    """
    :type lod: list[dict[A,B]]
    :rtype: iter[tuple[A,B]]
    """
    for l in lod:
        for k, v in l.items():
            yield k, v


def from_cx(cx):
    """Rebuilds a BELGraph from CX JSON output from PyBEL

    :param list[dict] cx: The CX JSON for this graph
    :rtype: BELGraph
    """
    graph = BELGraph()

    context_aspect = {}
    context_legend_aspect = []
    annotation_lists_aspect = []
    network_attributes_aspect = []
    nodes_aspect = []
    nodes_attributes_aspect = []
    edges_attributes_aspect = []
    edges_aspect = []
    meta_entries = defaultdict(list)

    for key, value in _iterate_list_of_dicts(cx):
        if key == 'context_legend':
            context_legend_aspect.extend(value)

        elif key == 'annotation_lists':
            annotation_lists_aspect.extend(value)

        elif key == '@context':
            context_aspect.update(value[0])

        elif key == 'networkAttributes':
            network_attributes_aspect.extend(value)

        elif key == 'nodes':
            nodes_aspect.extend(value)

        elif key == 'nodeAttributes':
            nodes_attributes_aspect.extend(value)

        elif key == 'edges':
            edges_aspect.extend(value)

        elif key == 'edgeAttributes':
            edges_attributes_aspect.extend(value)

        else:
            meta_entries[key].extend(value)

    context_legend = _cx_to_dict(context_legend_aspect)

    annotation_lists = defaultdict(set)
    for data in annotation_lists_aspect:
        annotation_lists[data['k']].add(data['v'])

    for keyword, entry in context_aspect.items():
        if context_legend[keyword] == GRAPH_NAMESPACE_URL:
            graph.namespace_url[keyword] = entry
        elif context_legend[keyword] == GRAPH_NAMESPACE_PATTERN:
            graph.namespace_pattern[keyword] = entry
        elif context_legend[keyword] == GRAPH_ANNOTATION_URL:
            graph.annotation_url[keyword] = entry
        elif context_legend[keyword] == GRAPH_ANNOTATION_PATTERN:
            graph.annotation_pattern[keyword] = entry
        elif context_legend[keyword] == GRAPH_ANNOTATION_LIST:
            graph.annotation_list[keyword] = annotation_lists[entry]

    for data in network_attributes_aspect:
        if data['n'] == NDEX_SOURCE_FORMAT:
            continue
        graph.graph[GRAPH_METADATA][data['n']] = data['v']

    nid_to_name = {}
    for data in nodes_aspect:
        nid = data['@id']
        node_name = data['n']
        nid_to_name[nid] = node_name

    nid_to_flattened_data = defaultdict(dict)
    for data in nodes_attributes_aspect:
        nid = data['po']
        nid_to_flattened_data[nid][data['n']] = data['v']

    # put all normal data here
    node_data_pp = defaultdict(dict)

    # Group all fusion-related data here
    node_data_fusion = defaultdict(dict)

    # Group all variant-related data
    node_data_variants = defaultdict(lambda: defaultdict(dict))

    for nid, flattened_data in nid_to_flattened_data.items():
        for key, value in flattened_data.items():
            if key.startswith(FUSION):
                node_data_fusion[nid][key] = value
            elif key.startswith(VARIANTS):
                _, i, vls = key.split('_', 2)
                node_data_variants[nid][i][vls] = value
            elif key in _listy_keys:
                node_data_pp[nid][key] = json.loads(value)
            else:
                node_data_pp[nid][key] = value

    for nid, data in node_data_fusion.items():
        node_data_pp[nid].update(expand_dict(data))

    for nid, data in node_data_variants.items():
        node_data_pp[nid][VARIANTS] = [
            expand_dict(value)
            for _, value in sorted(data.items())
        ]

    nid_node = {}
    for nid, data in node_data_pp.items():
        if CX_NODE_NAME in data:
            data[NAME] = data.pop(CX_NODE_NAME)

        node = dict_to_entity(data)
        nid_node[nid] = node
        graph.add_entity(node)

    edge_relation = {}
    edge_source = {}
    edge_target = {}
    for data in edges_aspect:
        edge_key = data['@id']  # SHA512 hash of the edge

        edge_relation[edge_key] = data['i']

        source_nid =data['s']
        source_node = nid_node[source_nid]
        edge_source[edge_key] = source_node

        target_nid = data['t']
        target_node = nid_node[target_nid]
        edge_target[edge_key] = target_node

    edge_data = defaultdict(dict)
    for data in edges_attributes_aspect:
        edge_key = data['po']  # SHA512 hash of the edge
        edge_data[edge_key][data['n']] = data['v']

    edge_citation = defaultdict(dict)
    edge_subject = defaultdict(dict)
    edge_object = defaultdict(dict)
    edge_annotations = defaultdict(lambda: defaultdict(dict))

    edge_data_pp = defaultdict(dict)

    for edge_key, data in edge_data.items():
        for key, value in data.items():

            if key.startswith(CITATION):
                _, vl = key.split('_', 1)
                edge_citation[edge_key][vl] = value

            elif key.startswith(SUBJECT):
                _, vl = key.split('_', 1)
                edge_subject[edge_key][vl] = value

            elif key.startswith(OBJECT):
                _, vl = key.split('_', 1)
                edge_object[edge_key][vl] = value

            elif key == EVIDENCE:
                edge_data_pp[edge_key][EVIDENCE] = value

            else:
                edge_annotations[edge_key][key] = value

    for key, data in edge_citation.items():
        edge_data_pp[key][CITATION] = data

    for key, data in edge_subject.items():
        edge_data_pp[key][SUBJECT] = expand_dict(data)

    for key, data in edge_object.items():
        edge_data_pp[key][OBJECT] = expand_dict(data)

    for key in edge_relation:
        if key in edge_annotations:  # FIXME stick this in edge_data.items() iteration
            edge_data_pp[key][ANNOTATIONS] = {
                key: {value: True for value in values}
                for key, values in edge_annotations[key].items()
            }

        u = edge_source[key]
        v = edge_target[key]
        relation = edge_relation[key]

        if key in edge_citation:  # being present in this dictionary this means that there is a citation
            graph.add_qualified_edge(
                u,
                v,
                relation=relation,
                citation=edge_data_pp[key][CITATION],
                evidence=edge_data_pp[key][EVIDENCE],
                subject_modifier=edge_data_pp[key].get(SUBJECT),
                object_modifier=edge_data_pp[key].get(OBJECT),
                annotations=edge_data_pp[key].get(ANNOTATIONS)
            )
        elif relation in unqualified_edges:
            graph.add_unqualified_edge(u, v, relation)
        else:
            raise ValueError('problem adding edge: {}'.format(key))

    return graph


def from_cx_jsons(graph_cx_json_str):
    """Reconstitutes a BEL graph from a CX JSON string
    
    :param str graph_cx_json_str: CX JSON string 
    :return: A BEL graph representing the CX graph contained in the string
    :rtype: BELGraph
    """
    graph_cx_json_dict = json.loads(graph_cx_json_str)
    return from_cx(graph_cx_json_dict)


def from_cx_file(file):
    """Reads a file containing CX JSON and converts to a BEL graph

    :param file file: A readable file or file-like containing the CX JSON for this graph
    :return: A BEL Graph representing the CX graph contained in the file
    :rtype: BELGraph
    """
    graph_cx_json_dict = json.load(file)
    return from_cx(graph_cx_json_dict)
