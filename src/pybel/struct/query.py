# -*- coding: utf-8 -*-

import json
import logging
from collections import Iterable

import numpy as np

from .exc import QueryMissingNetworksError
from .graph import union
from .mutation import get_subgraph
from .pipeline import Pipeline
from ..constants.induce_subgraph_keys import (
    NONNODE_SEED_TYPES, SEED_TYPE_ANNOTATION, SEED_TYPE_INDUCTION, SEED_TYPE_NEIGHBORS, SEED_TYPE_SAMPLE,
)
from ..dsl.nodes import BaseEntity
from ..manager.models import Node
from ..tokens import dict_to_entity

__all__ = [
    'Query',
]

log = logging.getLogger(__name__)

SEED_METHOD = 'type'
SEED_DATA = 'data'


def _get_seed():
    return np.random.randint(0, np.iinfo('i').max)


class Query(object):
    """Wraps a query over the network store"""

    def __init__(self, network_ids=None, seeding=None, pipeline=None):
        """
        :param int or iter[int] network_ids: Database network identifiers identifiers
        :param Optional[list[dict]] seeding:
        :param Optional[Pipeline] pipeline: Instance of a pipeline
        """
        if isinstance(network_ids, int):
            self.network_ids = [network_ids]

        elif isinstance(network_ids, Iterable):
            pool = tuple(network_ids)

            if any(not isinstance(entry, int) for entry in network_ids):
                raise TypeError('network identifiers entry is not int: {}'.format(network_ids))

            self.network_ids = list(pool)

        elif network_ids is None:
            self.network_ids = []

        else:
            raise TypeError('network identifiers is not list: {}'.format(network_ids))

        self.seeding = seeding or []
        self.pipeline = pipeline or Pipeline()

    def append_network(self, network_id):
        """Adds a network to this query

        :param int network_id: The database identifier of the network
        """
        self.network_ids.append(network_id)

    def _append_seed(self, seed_type, data):
        """Adds a seeding method

        :param str seed_type:
        :param data:
        """
        self.seeding.append({
            SEED_METHOD: seed_type,
            SEED_DATA: data
        })

    def append_seeding_induction(self, nodes):
        """Adds a seed induction method

        :param list[tuple or Node or BaseEntity] nodes: A list of PyBEL node tuples
        """
        self._append_seed(SEED_TYPE_INDUCTION, nodes)

    def append_seeding_neighbors(self, nodes):
        """Adds a seed by neighbors

        :param list[tuple or Node or BaseEntity] nodes: A list of PyBEL node tuples
        """
        self._append_seed(SEED_TYPE_NEIGHBORS, nodes)

    def append_seeding_annotation(self, annotation, values):
        """Adds a seed induction method for single annotation's values

        :param str annotation: The annotation to filter by
        :param set[str] values: The values of the annotation to keep
        """
        log.debug('appending seed by annotation=%s and values=%s', annotation, values)
        self._append_seed(SEED_TYPE_ANNOTATION, {
            'annotations': {
                annotation: values
            }
        })

    def append_seeding_sample(self, **kwargs):
        """Adds seed induction methods.

        Kwargs can have ``number_edges`` or ``number_seed_nodes``.
        """
        data = {
            'seed': _get_seed()
        }
        data.update(kwargs)

        self._append_seed(SEED_TYPE_SAMPLE, data)

    def append_pipeline(self, name, *args, **kwargs):
        """Adds an entry to the pipeline. Defers to :meth:`pybel_tools.pipeline.Pipeline.append`.

        :param name: The name of the function
        :type name: str or (pybel.BELGraph -> pybel.BELGraph)
        :return: This pipeline for fluid query building
        :rtype: Pipeline
        """
        return self.pipeline.append(name, *args, **kwargs)

    def __call__(self, manager, in_place=True):
        """Runs this query and returns the resulting BEL graph with :meth:`Query.run`

        :param pybel.manager.Manager manager: A cache manager
        :param bool in_place: Should the graph be copied before applying the algorithm?
        :rtype: Optional[pybel.BELGraph]
        """
        return self.run(manager, in_place=in_place)

    def run(self, manager, in_place=True):
        """Runs this query and returns the resulting BEL graph

        :param pybel.manager.Manager manager: A cache manager
        :param bool in_place: Should the graph be copied before applying the algorithm?
        :rtype: Optional[pybel.BELGraph]
        """
        log.debug('query universe consists of networks: %s', self.network_ids)

        if not self.network_ids:
            log.debug('can not run query without network identifiers')
            return

        universe = manager.get_graph_by_ids(self.network_ids)

        log.debug(
            'query universe has %d nodes/%d edges',
            universe.number_of_nodes(),
            universe.number_of_edges()
        )

        # parse seeding stuff

        if not self.seeding:
            return self.pipeline.run(universe, universe=universe, in_place=in_place)

        subgraphs = []

        for seed in self.seeding:
            seed_method, seed_data = seed[SEED_METHOD], seed[SEED_DATA]

            log.debug('seeding with %s: %s', seed_method, seed_data)
            subgraph = get_subgraph(universe, seed_method=seed_method, seed_data=seed_data)

            if subgraph is None:
                log.debug('seed returned empty graph: %s', seed)
                continue

            subgraphs.append(subgraph)

        if not subgraphs:
            log.debug('no subgraphs returned')
            return

        graph = union(subgraphs)

        return self.pipeline.run(graph, universe=universe, in_place=in_place)

    def seeding_to_jsons(self):
        """Returns seeding JSON as a string

        :rtype: str
        """
        return json.dumps(self.seeding)

    def to_json(self):
        """Returns this query as a JSON object

        :rtype: dict
        """
        rv = {
            'network_ids': list(self.network_ids),
        }

        if self.seeding:
            rv['seeding'] = self.seeding

        if self.pipeline:
            rv['pipeline'] = self.pipeline.to_json()

        return rv

    def to_jsons(self):
        """Returns this query as a stringified JSON object

        :rtype: str
        """
        return json.dumps(self.to_json())

    @staticmethod
    def from_jsons(s):
        """Loads a query from a stringified JSON object

        :param str s: A stringified JSON query
        :rtype: Query
        :raises: QueryMissingNetworks
        """
        return Query.from_json(json.loads(s))

    @staticmethod
    def from_json(data):
        """Loads a query from a JSON dictionary

        :param dict data: A JSON dictionary
        :rtype: Query
        :raises: QueryMissingNetworks
        """
        network_ids = data.get('network_ids')
        if network_ids is None:
            raise QueryMissingNetworksError('query JSON did not have key "network_ids"')

        if 'pipeline' in data:
            pipeline = Pipeline.from_json(data['pipeline'])
        else:
            pipeline = None

        if 'seeding' in data:
            seeding = Query._process_seeding(data['seeding'])
        else:
            seeding = None

        return Query(
            network_ids=network_ids,
            seeding=seeding,
            pipeline=pipeline
        )

    @staticmethod
    def _process_seeding(seeds):  # TODO write documentation
        """Makes sure nodes are tuples and not lists once back in

        :param seeds:
        :rtype: list[dict]
        """
        return [
            {
                SEED_METHOD: seed[SEED_METHOD],
                SEED_DATA: [
                    dict_to_entity(node)
                    for node in seed[SEED_DATA]
                ]
                if seed[SEED_METHOD] not in NONNODE_SEED_TYPES else seed[SEED_DATA]
            }
            for seed in seeds
        ]
