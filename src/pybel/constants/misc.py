# -*- coding: utf-8 -*-

from json import load
from logging import getLogger
from os import environ, makedirs, mkdir, path

from .bel_graph_keywords import *
from .bel_keywords import *

log = getLogger(__name__)

#: The last PyBEL version where the graph data definition changed
PYBEL_MINIMUM_IMPORT_VERSION = 0, 11, 0

BELFRAMEWORK_DOMAIN = 'http://resource.belframework.org'
OPENBEL_DOMAIN = 'http://resources.openbel.org'

SMALL_CORPUS_URL = OPENBEL_DOMAIN + '/belframework/20150611/knowledge/small_corpus.bel'
LARGE_CORPUS_URL = OPENBEL_DOMAIN + '/belframework/20150611/knowledge/large_corpus.bel'

FRAUNHOFER_RESOURCES = 'https://owncloud.scai.fraunhofer.de/index.php/s/JsfpQvkdx3Y5EMx/download?path='
OPENBEL_NAMESPACE_RESOURCES = OPENBEL_DOMAIN + '/belframework/20150611/namespace/'
OPENBEL_ANNOTATION_RESOURCES = OPENBEL_DOMAIN + '/belframework/20150611/annotation/'

#: GOCC is the only namespace that needs to be stored because translocations use some of its values by default
GOCC_LATEST = 'https://arty.scai.fraunhofer.de/artifactory/bel/namespace/go-cellular-component/go-cellular-component-20170511.belns'
GOCC_KEYWORD = 'GOCC'

#: The environment variable that contains the default SQL connection information for the PyBEL cache
PYBEL_CONNECTION = 'PYBEL_CONNECTION'

#: The default directory where PyBEL files, including logs and the  default cache, are stored. Created if not exists.
PYBEL_DIR = environ.get('PYBEL_RESOURCE_DIRECTORY', path.join(path.expanduser('~'), '.pybel'))
if not path.exists(PYBEL_DIR):
    try:
        mkdir(PYBEL_DIR)
    except FileExistsError:
        log.debug('pybel data directory was created already: %s', PYBEL_DIR)

DEFAULT_CACHE_NAME = 'pybel_{}.{}.{}_cache.db'.format(*PYBEL_MINIMUM_IMPORT_VERSION)
#: The default cache location is ``~/.pybel/data/pybel_cache.db``
DEFAULT_CACHE_LOCATION = path.join(PYBEL_DIR, DEFAULT_CACHE_NAME)
#: The default cache connection string uses sqlite.
DEFAULT_CACHE_CONNECTION = 'sqlite:///' + DEFAULT_CACHE_LOCATION


def get_config_dir():
    """Returns the path to the directory where configuration is stored for PyBEL. Can be overridden by setting the
    environment variable ``PYBEL_CONFIG_DIRECTORY``.

    :rtype: str
    """
    return environ.get('PYBEL_CONFIG_DIRECTORY', path.join(path.expanduser('~'), '.config', 'pybel'))


_config_dir = get_config_dir()
if not path.exists(_config_dir):
    try:
        makedirs(_config_dir)
    except FileExistsError:
        log.debug('config folder was already created: %s', _config_dir)

#: The global configuration for PyBEL is stored here. By default, it loads from ``~/.config/pybel/config.json``
config = {
    PYBEL_CONNECTION: DEFAULT_CACHE_CONNECTION
}


def get_config_path():
    """Returns the path of the configuration file, which should just be a file called ``config.json`` inside the
    directory returned by :func:`get_config_dir`

    :rtype: str
    """
    return path.join(_config_dir, 'config.json')


_config_path = get_config_path()
if path.exists(_config_path):
    with open(_config_path) as f:
        config.update(load(f))


def get_cache_connection(connection=None):
    """Returns the default cache connection string. If a connection string is explicitly given, passes it through

    :param str connection: RFC connection string
    :rtype: str
    """
    if connection is not None:
        log.debug('getting user-defined connection: %s', connection)
        return connection

    if PYBEL_CONNECTION in environ:
        log.debug('getting environment-defined connection: %s', environ[PYBEL_CONNECTION])
        return environ[PYBEL_CONNECTION]

    log.debug('getting default connection %s', config[PYBEL_CONNECTION])
    return config[PYBEL_CONNECTION]


PYBEL_CONTEXT_TAG = 'pybel_context'
PYBEL_AUTOEVIDENCE = 'Automatically added by PyBEL'

#: The default namespace given to entities in the BEL language
BEL_DEFAULT_NAMESPACE = 'bel'

#: Provides a mapping from BEL language keywords to internal PyBEL strings
DOCUMENT_KEYS = {
    BEL_KEYWORD_METADATA_AUTHORS: METADATA_AUTHORS,
    BEL_KEYWORD_METADATA_CONTACT: METADATA_CONTACT,
    BEL_KEYWORD_METADATA_COPYRIGHT: METADATA_COPYRIGHT,
    BEL_KEYWORD_METADATA_DESCRIPTION: METADATA_DESCRIPTION,
    BEL_KEYWORD_METADATA_DISCLAIMER: METADATA_DISCLAIMER,
    BEL_KEYWORD_METADATA_LICENSES: METADATA_LICENSES,
    BEL_KEYWORD_METADATA_NAME: METADATA_NAME,
    BEL_KEYWORD_METADATA_VERSION: METADATA_VERSION,
    BEL_KEYWORD_METADATA_PROJECT: METADATA_PROJECT,
}

#: Provides a mapping from internal PyBEL strings to BEL language keywords. Is the inverse of :data:`DOCUMENT_KEYS`
INVERSE_DOCUMENT_KEYS = {v: k for k, v in DOCUMENT_KEYS.items()}

#: The default location of BEL Commons
DEFAULT_SERVICE_URL = 'https://bel-commons.scai.fraunhofer.de'
