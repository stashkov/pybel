# -*- coding: utf-8 -*-

#: Represents the key for the citation type in a citation dictionary
CITATION_TYPE = 'type'
#: Represents the key for the citation name in a citation dictionary
CITATION_NAME = 'name'
#: Represents the key for the citation reference in a citation dictionary
CITATION_REFERENCE = 'reference'
#: Represents the key for the citation date in a citation dictionary
CITATION_DATE = 'date'
#: Represents the key for the citation authors in a citation dictionary
CITATION_AUTHORS = 'authors'
#: Represents the key for the citation comment in a citation dictionary
CITATION_COMMENTS = 'comments'

#: Represents the key for the optional PyBEL citation title entry in a citation dictionary
CITATION_TITLE = 'title'
#: Represents the key for the optional PyBEL citation volume entry in a citation dictionary
CITATION_VOLUME = 'volume'
#: Represents the key for the optional PyBEL citation issue entry in a citation dictionary
CITATION_ISSUE = 'issue'
#: Represents the key for the optional PyBEL citation pages entry in a citation dictionary
CITATION_PAGES = 'pages'
#: Represents the key for the optional PyBEL citation first author entry in a citation dictionary
CITATION_FIRST_AUTHOR = 'first'
#: Represents the key for the optional PyBEL citation last author entry in a citation dictionary
CITATION_LAST_AUTHOR = 'last'

#: Represents the ordering of the citation entries in a control statement (SET Citation = ...)
CITATION_ENTRIES = CITATION_TYPE, CITATION_NAME, CITATION_REFERENCE, CITATION_DATE, CITATION_AUTHORS, CITATION_COMMENTS
