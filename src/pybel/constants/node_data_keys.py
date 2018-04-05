# -*- coding: utf-8 -*-

# Internal node data format keys
#: The node data key specifying the node's function (e.g. :data:`GENE`, :data:`MIRNA`, :data:`BIOPROCESS`, etc.)
FUNCTION = 'function'
#: The key specifying an identifier dictionary's namespace. Used for nodes, activities, and transformations.
NAMESPACE = 'namespace'
#: The key specifying an identifier dictionary's name. Used for nodes, activities, and transformations.
NAME = 'name'
#: The key specifying an identifier dictionary
IDENTIFIER = 'identifier'
#: The key specifying an optional label for the node
LABEL = 'label'
#: The key specifying an optional description for the node
DESCRIPTION = 'description'

#: They key representing the nodes that are a member of a composite or complex
MEMBERS = 'members'
#: The key representing the nodes appearing in the reactant side of a biochemical reaction
REACTANTS = 'reactants'
#: The key representing the nodes appearing in the product side of a biochemical reaction
PRODUCTS = 'products'

#: The node data key specifying a fusion dictionary, containing :data:`PARTNER_3P`, :data:`PARTNER_5P`,
# :data:`RANGE_3P`, and :data:`RANGE_5P`
FUSION = 'fusion'
#: The key specifying the identifier dictionary of the fusion's 3-Prime partner
PARTNER_3P = 'partner_3p'
#: The key specifying the identifier dictionary of the fusion's 5-Prime partner
PARTNER_5P = 'partner_5p'
#: The key specifying the range dictionary of the fusion's 3-Prime partner
RANGE_3P = 'range_3p'
#: The key specifying the range dictionary of the fusion's 5-Prime partner
RANGE_5P = 'range_5p'

FUSION_REFERENCE = 'reference'
FUSION_START = 'left'
FUSION_STOP = 'right'
FUSION_MISSING = 'missing'

#: The key specifying the node has a list of associated variants
VARIANTS = 'variants'

#: The group of all BEL-provided keys for node data dictionaries, used for hashing.
PYBEL_NODE_DATA_KEYS = {
    FUNCTION,
    NAMESPACE,
    NAME,
    IDENTIFIER,
    VARIANTS,
    FUSION,
    MEMBERS,
    REACTANTS,
    PRODUCTS,
}

#: Used as a namespace when none is given when lenient parsing mode is turned on. Not recommended!
DIRTY = 'dirty'

