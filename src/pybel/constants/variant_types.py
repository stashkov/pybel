# -*- coding: utf-8 -*-

#: The key representing what kind of variation is being represented
KIND = 'kind'
#: The value for :data:`KIND` for an HGVS variant
HGVS = 'hgvs'
#: The value for :data:`KIND` for a protein modification
PMOD = 'pmod'
#: The value for :data:`KIND` for a gene modification
GMOD = 'gmod'
#: The value for :data:`KIND` for a fragment
FRAGMENT = 'frag'

#: The allowed values for :data:`KIND`
PYBEL_VARIANT_KINDS = {
    HGVS,
    PMOD,
    GMOD,
    FRAGMENT
}
