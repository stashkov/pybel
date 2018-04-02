# -*- coding: utf-8 -*-

"""This example shows the downstream effects of carbamazepine on the GABAergic receptors in Epilepsy

.. code-block::

    SET Citation = {"PubMed", "22379998"}
    SET Evidence = "Carbamazepine (CBZ), lamotrigine (LTG) and phenytoin (PHT), among a class of antiepileptic drugs
    (AEDs) involved in sodium-channel modulation, are three of the leading first-line treatments for epilepsy"

    a(CHEBI:carbamazepine) -| path(MESHD:Epilepsy)
    a(CHEBI:carbamazepine) -- bp(GOBP:"regulation of voltage-gated sodium channel activity")

    UNSET ALL

    SET Citation = {"PubMed", "23535492"}
    SET Evidence = "The role of adenosine as an endogenous anticonvulsant is mediated via ADORA1"

    a(CHEBI:adenosine) -| path(MESHD:Seizures)
    a(CHEBI:adenosine) isA a(CHEBI:anticonvulsant)
    p(HGNC:ADORA1) -> a(CHEBI:adenosine)

    SET Citation = {"PubMed", "23535492"}
    SET Evidence = "In a rat model of seizure kindling, ADORA1 in the hippocampal CA1 region reduces seizures, whereas
    ADORA2A promotes them."

    SET MeSHAnatomy="CA1 Region, Hippocampal"
    SET Species = "10116"

    p(HGNC:ADORA1) -| path(MESHD:Seizures)
    p(HGNC:ADORA2A) -> path(MESHD:Seizures)

    UNSET ALL

    SET Citation = {"PubMed", "22389222"}
    SET Evidence = "The discovery that through a Ca2+-dependent glutamate release astrocytes can directly excite groups
    of neighbouring neurons (Parpura et al. 1994) and favour sySCOMPronised activities mediated by extrasynaptic
    N-methyl-D-aspartate (NMDA) receptor activation (Fellin et al. 2004) were the initial observations that hinted at
    a more direct role of Ca2+-dependent gliotransmission in the generation of epileptiform activities"

    SET MeSHAnatomy = "Astrocytes"
    SET Subgraph = {"Glutamatergic subgraph", "Neurotransmitter release subgraph"}

    a(CHEBI:"calcium(2+)") -> bp(GOBP:"glutamate secretion, neurotransmission")

    UNSET ALL

    SET Citation = {"PubMed", "2553432"}
    SET Evidence = "Potentialisation by muscimol the effect of CBZ may suggest that the reduction of immobility by CBZ
    in the behavioral despair test is mediated through the activation of GABA-a receptors, linked to the chloride
    ionophore"

    a(CHEBI:carbamazepine) -> p(HGNC:GABRA1)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRA2)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRA3)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRA4)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRA5)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRA6)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRB1)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRB2)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRB3)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRD)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRE)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRG1)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRG2)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRG3)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRP)
    a(CHEBI:carbamazepine) -> p(HGNC:GABRQ)

    UNSET ALL

    SET Citation = {"PubMed", "18093662"}
    SET Evidence = "Unlike conventional antidepressant drugs, of which mechanisms underlying the therapeutic actions
    are largely derived from enhancement of central catecholaminergic or indoleaminergic systems by inhibiting relevant
    metabolic enzymes or reuptake process, psychotrophic effects of anticonvulsants are thought to be mainly associated
    with the facilitation of central gamma-aminobutyric acid (GABA) neurotransmission or blockade of glutaminergic
    hyperexcitability"

    a(CHEBI:anticonvulsant) -> bp(GOBP:"synaptic transmission, GABAergic")
    a(CHEBI:anticonvulsant) -| bp(GOBP:"synaptic transmission, glutamatergic")

    UNSET ALL

    SET Citation = {"PubMed", "17054941"}
    SET Evidence = "It has been demonstrated that the initial depression of excitatory synaptic transmission and PS
    during energy deprivation is caused primarily by an increase of extracellular adenosine (Fowler, 1989, 1990; Cassar
    et al., 1998; Latini et al., 1998, 1999; Dale et al., 2000; Gervitz et al., 2001; Pearson et al., 2003). This
    increase in extracellular adenosine constitutes an endogenous neuroprotective process that is mediated by adenosine
    A1-receptors (Fowler, 1989, 1990; Latini et al., 1999), leading to reduced neuronal activity with lowered energy
    consumption, hampering the excitotoxic effects of glutamate"

    a(CHEBI:carbamazepine) -> a(CHEBI:adenosine)

    UNSET ALL
"""

from ..constants import *
from ..dsl.edges import translocation
from ..dsl.nodes import bioprocess, complex_abundance, protein, abundance, pathology
from ..dsl.utils import entity
from ..struct.graph import BELGraph

__all__ = [
    'gaba_graph'
]

gaba_graph = BELGraph(
    name='Carbamazepine Pathway',
    version='1.0.0',
    description="The downstream effects of carbamazepine in Epilepsy",
    authors='Charles Tapley Hoyt',
    contact='charles.hoyt@scai.fraunhofer.de',
)

gaba_graph.namespace_url.update({
    'HGNC': 'https://arty.scai.fraunhofer.de/artifactory/bel/namespace/hgnc-human-genes/hgnc-human-genes-20170725.belns',
    'CHEBI': 'https://arty.scai.fraunhofer.de/artifactory/bel/namespace/chebi/chebi-20170725.belns',
    'GOBP': 'https://arty.scai.fraunhofer.de/artifactory/bel/namespace/go-biological-process/go-biological-process-20170725.belns'
})

gaba_graph.annotation_url.update({
    'Confidence': 'https://arty.scai.fraunhofer.de/artifactory/bel/annotation/confidence/confidence-1.0.0.belanno',
    'Species': 'https://arty.scai.fraunhofer.de/artifactory/bel/annotation/species-taxonomy-id/species-taxonomy-id-20170511.belanno'
})


carbamazepine = abundance(namespace='CHEBI', name='carbamazeine')
epilepsy = pathology(namespace='MESH', name='Epilepsy')
sodium_channel_modulation = bioprocess(namespace='GO', name='regulation of voltage-gated sodium channel activity')

"""
SET Citation = {"PubMed", "22379998"}
    SET Evidence = "Carbamazepine (CBZ), lamotrigine (LTG) and phenytoin (PHT), among a class of antiepileptic drugs
    (AEDs) involved in sodium-channel modulation, are three of the leading first-line treatments for epilepsy"

    a(CHEBI:carbamazepine) -| path(MESHD:Epilepsy)
    a(CHEBI:carbamazepine) -- bp(GOBP:"regulation of voltage-gated sodium channel activity")

    UNSET ALL
"""


evidence = "Carbamazepine (CBZ), lamotrigine (LTG) and phenytoin (PHT), among a class of antiepileptic drugs (AEDs) " \
           "involved in sodium-channel modulation, are three of the leading first-line treatments for epilepsy"

gaba_graph.add_is_a(carbamazepine,)

gaba_graph.add_qualified_edge(
    carbamazepine,
    epilepsy,
    relation=DECREASES,
    citation="22379998",
    evidence=evidence
)

gaba_graph.add_qualified_edge(
    carbamazepine,
    sodium_channel_modulation,
    relation=ASSOCIATION,
    citation="22379998",
    evidence=evidence
)