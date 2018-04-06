# -*- coding: utf-8 -*-

"""

Node-Link JSON
--------------
This module contains IO functions for interconversion between BEL graphs and Node-Link JSON

"""

import json
import os
from itertools import chain, count

from .utils import ensure_version
from ..constants import GRAPH_ANNOTATION_LIST, GRAPH_UNCACHED_NAMESPACES
from ..struct import BELGraph
from ..tokens import dict_to_entity
from ..utils import list2tuple

__all__ = [
    'to_json',
    'to_json_file',
    'to_json_path',
    'to_jsons',
    'from_json',
    'from_json_file',
    'from_json_path',
    'from_jsons',
]


def _augment_node_with_sha512(node):
    """

    :param BaseEntity node:
    :return:
    """
    v = node.copy()
    v['id'] = node.as_sha512()
    return v


def _node_link_data(graph):
    """Converts a BEL graph to a node-link format

    Adapted from :func:`networkx.readwrite.json_graph.node_link_data`

    :param pybel.BELGraph graph:
    :rtype: dict
    """
    nodes = sorted(graph)

    mapping = dict(zip(nodes, count()))

    return {
        'directed': True,
        'multigraph': True,
        'graph': graph.graph,
        'nodes': [
            _augment_node_with_sha512(node)
            for node in nodes
        ],
        'links': [
            dict(chain(
                data.items(),
                [('source', mapping[u]), ('target', mapping[v]), ('key', key)]
            ))
            for u, v, key, data in graph.edges(keys=True, data=True)
        ]
    }


def to_json(graph):
    """Converts this graph to a Node-Link JSON object

    :param BELGraph graph: A BEL graph
    :return: A Node-Link JSON object representing the given graph
    :rtype: dict
    """
    rv = _node_link_data(graph)

    # Convert annotation list definitions (which are sets) to canonicalized/sorted lists
    rv['graph'][GRAPH_ANNOTATION_LIST] = {
        key: list(sorted(values))
        for key, values in rv['graph'][GRAPH_ANNOTATION_LIST].items()
    }

    # Convert set to list
    rv['graph'][GRAPH_UNCACHED_NAMESPACES] = list(sorted(rv['graph'][GRAPH_UNCACHED_NAMESPACES]))

    return rv


def to_json_path(graph, path, **kwargs):
    """Writes this graph to the given path as a Node-Link JSON

    :param BELGraph graph: A BEL graph
    :param str path: A file path
    """
    with open(os.path.expanduser(path), 'w') as f:
        return to_json_file(graph, file=f, **kwargs)


def to_json_file(graph, file, **kwargs):
    """Writes this graph as Node-Link JSON to a file

    :param BELGraph graph: A BEL graph
    :param file file: A write-supporting file or file-like object
    """
    graph_json_dict = to_json(graph)
    json.dump(graph_json_dict, file, ensure_ascii=False, **kwargs)


def to_jsons(graph, **kwargs):
    """Dumps this graph as a Node-Link JSON object to a string

    :param BELGraph graph: A BEL graph
    :return: A string representation of the Node-Link JSON produced for this graph by :func:`pybel.to_json`
    :rtype: str
    """
    graph_json_str = to_json(graph)
    return json.dumps(graph_json_str, ensure_ascii=False, **kwargs)


def _node_link_graph(data):
    """Return graph from node-link data format

    Adapted from :func:`networkx.readwrite.json_graph.node_link_graph`

    :param dict data:
    :rtype: BELGraph
    """
    graph = BELGraph()
    graph.graph = data.get('graph', {})

    mapping = []

    for node_data in data['nodes']:
        node = dict_to_entity(node_data)
        graph.add_node_from_data(node)
        mapping.append(node)

    for data in data['links']:
        src = data['source']
        tgt = data['target']
        key = data['key']

        u = mapping[src]
        v = mapping[tgt]

        edgedata = {
            k: v
            for k, v in data.items()
            if k not in {'source', 'target', 'key'}
        }
        graph.add_edge(u, v, key=key, **edgedata)

    return graph


def from_json(graph_json_dict, check_version=True):
    """Builds a graph from Node-Link JSON Object

    :param dict graph_json_dict: A JSON dictionary representing a graph
    :param bool check_version: Checks if the graph was produced by this version of PyBEL
    :rtype: BELGraph
    """
    for i, node in enumerate(graph_json_dict['nodes']):
        graph_json_dict['nodes'][i]['id'] = list2tuple(graph_json_dict['nodes'][i]['id'])

    graph = _node_link_graph(graph_json_dict)

    # FIXME need to convert all nodes into BaseEntities

    graph = BELGraph(data=graph)
    return ensure_version(graph, check_version=check_version)


def from_json_path(path, check_version=True):
    """Builds a graph from a file containing Node-Link JSON

    :param str path: A file path. Expands user.
    :param bool check_version: Checks if the graph was produced by this version of PyBEL
    :rtype: BELGraph
    """
    with open(os.path.expanduser(path)) as f:
        return from_json_file(f, check_version=check_version)


def from_json_file(file, check_version=True):
    """Builds a graph from the Node-Link JSON contained in the given file

    :param file file: A readable file or file-like
    :param bool check_version: Checks if the graph was produced by this version of PyBEL
    :rtype: BELGraph
    """
    graph_json_dict = json.load(file)
    return from_json(graph_json_dict, check_version=check_version)


def from_jsons(graph_json_str, check_version=True):
    """Reads a BEL graph from a Node-Link JSON string

    :param str graph_json_str: A Node-Link JSON string produced by PyBEL
    :param bool check_version: Checks if the graph was produced by this version of PyBEL
    :rtype: BELGraph
    """
    graph_json_dict = json.loads(graph_json_str)
    return from_json(graph_json_dict, check_version=check_version)
