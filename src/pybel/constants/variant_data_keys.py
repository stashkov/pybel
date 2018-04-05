# -*- coding: utf-8 -*-

from .node_data_keys import IDENTIFIER
from .variant_types import KIND

# Modifier parser constants

#: The key for the starting position of a fragment range
FRAGMENT_START = 'start'
#: The key for the stopping position of a fragment range
FRAGMENT_STOP = 'stop'
#: The key signifying that there is neither a start nor stop position defined
FRAGMENT_MISSING = 'missing'
#: The key for any additional descriptive data about a fragment
FRAGMENT_DESCRIPTION = 'description'

#: The order for serializing gene modification data
GMOD_ORDER = [KIND, IDENTIFIER]

#: The key for the reference nucleotide in a gene substitution.
#: Only used during parsing  since this is converted to HGVS.
GSUB_REFERENCE = 'reference'
#: The key for the position of a gene substitution.
#: Only used during parsing  since this is converted to HGVS
GSUB_POSITION = 'position'
#: The key for the effect of a gene substitution.
#: Only used during parsing since this is converted to HGVS
GSUB_VARIANT = 'variant'

#: The key for the protein modification code.
PMOD_CODE = 'code'
#: The key for the protein modification position.
PMOD_POSITION = 'pos'
#: The order for serializing information about a protein modification
PMOD_ORDER = [KIND, IDENTIFIER, PMOD_CODE, PMOD_POSITION]

#: The key for the reference amino acid in a protein substitution.
#: Only used during parsing since this is concerted to HGVS
PSUB_REFERENCE = 'reference'
#: The key for the position of a protein substitution. Only used during parsing since this is converted to HGVS.
PSUB_POSITION = 'position'
#: The key for the variant of a protein substitution.Only used during parsing since this is converted to HGVS.
PSUB_VARIANT = 'variant'

#: The key for the position at which a protein is truncated
TRUNCATION_POSITION = 'position'
