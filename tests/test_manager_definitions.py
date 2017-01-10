import os
import tempfile
import unittest

import tests.constants
from pybel.manager.cache import CacheManager
from tests.constants import HGNC_URL, help_check_hgnc, CELL_LINE_URL, HGNC_KEYWORD
from tests.constants import wine_iri, mock_bel_resources

test_ns1 = 'file:///' + tests.constants.test_ns_1
test_ns2 = 'file:///' + tests.constants.test_ns_2
test_an1 = 'file:///' + tests.constants.test_an_1


class TestCachePersistient(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.db_path = os.path.join(self.dir, 'test.db')
        self.connection = 'sqlite:///' + self.db_path

    def tearDown(self):
        os.remove(self.db_path)
        os.rmdir(self.dir)

    @mock_bel_resources
    def test_insert_namespace(self, mock_get):
        cm1 = CacheManager(connection=self.connection)
        cm1.ensure_namespace(HGNC_URL)
        help_check_hgnc(self, {HGNC_KEYWORD: cm1.namespace_cache[HGNC_URL]})

        cm2 = CacheManager(connection=self.connection)
        cm2.ensure_namespace(HGNC_URL)
        help_check_hgnc(self, {HGNC_KEYWORD: cm2.namespace_cache[HGNC_URL]})


class TestCache(unittest.TestCase):
    def setUp(self):
        self.connection = 'sqlite:///'
        self.cm = CacheManager(connection=self.connection)

    @mock_bel_resources
    def test_insert_namespace(self, mock_get):
        self.cm.ensure_namespace(HGNC_URL)
        help_check_hgnc(self, {HGNC_KEYWORD: self.cm.namespace_cache[HGNC_URL]})

    @mock_bel_resources
    def test_insert_annotation(self, mock_get):
        self.cm.ensure_annotation(CELL_LINE_URL)
        self.assertIn(CELL_LINE_URL, self.cm.annotation_cache)
        self.assertIn('1321N1 cell', self.cm.annotation_cache[CELL_LINE_URL])
        self.assertEqual('CLO_0001072', self.cm.annotation_cache[CELL_LINE_URL]['1321N1 cell'])

    def test_insert_owl(self):
        self.cm.ensure_owl(wine_iri)
        self.assertIn(wine_iri, self.cm.term_cache)
        self.assertIn('ChateauMorgon', self.cm.term_cache[wine_iri])
        self.assertIn('Winery', self.cm.term_cache[wine_iri])
