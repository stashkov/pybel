# -*- coding: utf-8 -*-

#: Represents the BEL abundance, abundance()
ABUNDANCE = 'Abundance'

#: Represents the BEL abundance, geneAbundance()
#: .. seealso:: http://openbel.org/language/version_2.0/bel_specification_version_2.0.html#Xabundancea
GENE = 'Gene'

#: Represents the BEL abundance, rnaAbundance()
RNA = 'RNA'

#: Represents the BEL abundance, microRNAAbundance()
MIRNA = 'miRNA'

#: Represents the BEL abundance, proteinAbundance()
PROTEIN = 'Protein'

#: Represents the BEL function, biologicalProcess()
BIOPROCESS = 'BiologicalProcess'

#: Represents the BEL function, pathology()
PATHOLOGY = 'Pathology'

#: Represents the BEL abundance, compositeAbundance()
COMPOSITE = 'Composite'

#: Represents the BEL abundance, complexAbundance()
COMPLEX = 'Complex'

#: Represents the BEL transformation, reaction()
REACTION = 'Reaction'

#: A set of all of the valid PyBEL node functions
PYBEL_NODE_FUNCTIONS = {
    ABUNDANCE,
    GENE,
    RNA,
    MIRNA,
    PROTEIN,
    BIOPROCESS,
    PATHOLOGY,
    COMPOSITE,
    COMPLEX,
    REACTION
}

#: The mapping from PyBEL node functions to BEL strings
rev_abundance_labels = {
    ABUNDANCE: 'a',
    GENE: 'g',
    MIRNA: 'm',
    PROTEIN: 'p',
    RNA: 'r',
    BIOPROCESS: 'bp',
    PATHOLOGY: 'path',
    COMPLEX: 'complex',
    COMPOSITE: 'composite'
}
