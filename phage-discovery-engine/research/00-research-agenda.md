# Phase 0 — Self-Generated Research Agenda

**Status:** written before any literature work, so that it can be scored against
what I actually find. Questions marked **[SEED-RISK]** are ones where I expect the
answer may contradict or deflate the seed idea. Questions marked **[BLOCKING]**
change research direction and are being put to the client now.

The organizing principle: I am not allowed to propose an architecture until I can
answer, with citations, *what a machine can see in phage data that a human cannot,
that is not already seen by an existing tool, and that a bench experiment could
kill in under a month.*

---

## Thread A — How much of phage biology is actually unknown, and is "unknown" the same as "discoverable"?

A1. What is the real size and composition of the public phage sequence corpus in
2026 — isolate genomes (GenBank/INPHARED/PhagesDB) vs. metagenome-assembled
(IMG/VR, GPD, MGV)? How much is redundant?

A2. What fraction of phage ORFs are "hypothetical"? Is the widely repeated
"50–70%" figure measured per-gene, per-family, or per-genome, and against which
database? **[SEED-RISK]** — if most hypotheticals are singletons rather than large
conserved families, the "dark matter structure" mechanism collapses.

A3. How much of the hypothetical fraction is not dark matter at all but
*annotation failure*: missed small ORFs, stop-codon recoding, programmed
frameshifts, introns/inteins, wrong start sites? What is the published magnitude
of each? This is a competing explanation the seed does not consider.

A4. Has corpus-scale clustering of phage dark matter already been done (PHROGs,
efam, VOGdb, AFDB/ESM Atlas structural clusters, Foldseek-based)? If yes, what
did it *not* deliver, and why? **[SEED-RISK]** — "cluster the dark matter" may be
a solved and published step.

A5. What is the current state of structure-based phage annotation (phold,
ProstT5/3Di, Foldseek vs. PDB+AFDB+PHROG structures)? What is the marginal yield
of structure over profile HMMs, quantitatively?

A6. Where does taxonomy stand post-2021 (abolition of morphology-based families,
Caudoviricetes, genus/species thresholds)? Does taxonomic instability break any
proposed corpus-level statistic?

## Thread B — Host range: what is genuinely predictable today?

B1. What accuracy do host-prediction tools (iPHoP, CHERRY, RaFAH, PHIST, WIsH,
vHULK, PhageHostLearn) achieve, at what taxonomic rank, on what benchmark? What
is their performance at *strain* level — which is the only rank the lab cares
about? **[SEED-RISK]** — if genus-level is the ceiling, EOP prediction is out of
reach and the system must be about mechanism, not prediction.

B2. What is the mechanistic determinant stack — receptor identity, RBP/tailspike
architecture, depolymerases, capsule/LPS/O-antigen/pili loci, adsorption vs.
post-adsorption barriers, superinfection exclusion, abortive infection? Which
steps are individually predictable from sequence?

B3. How well can receptor specificity be predicted from RBP tip domains? What
happened when people tried (Boeckaerts, PhageRBPdetect, PhageHostLearn on
Klebsiella K-loci)? What are the reported ceilings?

B4. How reproducible is EOP itself? What are the known confounders (lysis from
without, MOI, agar overlay, temperature, plating host, phage stock heterogeneity,
lysogeny)? **[SEED-RISK]** — if EOP noise floor is ~1 log, half the "patterns" a
model finds in a lab matrix are measurement artifacts.

B5. Diversity-generating retroelements, phase variation, and phase-variable host
receptors: how much host-range variation is *stochastic* and therefore
unpredictable in principle?

## Thread C — The defense/counter-defense arms race

C1. How many anti-phage defense systems are now known, and what is the discovery
rate curve? Is the field saturating or still exponential? Who found them and by
what method (defense island guilt-by-association, Hachiman-style screens,
pan-immune analysis)?

C2. What do DefenseFinder and PADLOC actually cover, and what is their false
positive/negative behavior? What fraction of a typical genome's defense repertoire
do they miss?

C3. What is the state of anti-defense (anti-CRISPR, anti-restriction,
anti-Thoeris Tad1/Tad2, anti-CBASS Acb1/Acb2, nucleus-forming shielding, PARIS
antagonists)? Is there a database (dbAPIS, AntiDefenseFinder)? How were these
found — and could a machine have found them first? **[This is the retrospective
benchmark I most want to build.]**

C4. Has phage-gene × host-defense co-occurrence/mutual-exclusion already been
mined at corpus scale? By whom, with what statistics? **[SEED-RISK]** — the seed's
"COUPLING" mechanism may already be a 2023–2025 literature.

C5. What is the confounding structure of such analyses? Bacterial and phage
genomes are not independent samples — how badly does population structure inflate
co-occurrence signals, and what corrections exist (treeWAS, pyseer LMM, Coinfinder,
phylogenetic profiling with Pagel/Fisher-on-tree)? **[I expect this to be the
single most important technical point in the whole design.]**

## Thread D — What data actually exists, and can we use it?

D1. Genome corpora: scale, curation, update cadence, licensing, programmatic
access. INPHARED, PhagesDB, IMG/VR, GPD, MGV, RefSeq viral, PhageScope.

D2. Interaction data: is there any public phage–host *interaction matrix* at
scale, or only scattered supplementary tables? (This determines whether a learned
EOP model can be pre-trained at all.)

D3. Structure: AFDB coverage of phage proteins, ESM Atlas, PDB phage entries,
Foldseek cluster resources, BFVD (viral structure DB) if it exists.

D4. Literature: what is legally and technically retrievable at scale — PubMed
abstracts, PMC OA subset, bioRxiv, Europe PMC full text? What are the licensing
limits on full-text mining? What about figure/table extraction?

D5. Host genomes: how many sequenced strains of the relevant species exist, with
what metadata? Are defense repertoires already precomputed anywhere?

## Thread E — What already exists as a *tool*, and where does each stop?

E1. Annotation: Pharokka, phold, PHANOTATE, DRAMv, MultiPHATE. Where do they stop?
E2. Comparative/taxonomy: vConTACT3, VIRIDIC, VIRCLUST, GRAViTy, PhageClouds.
E3. Novelty detection: is there anything that flags "this is unlike everything"?
E4. Hypothesis generation: Google AI co-scientist, FutureHouse (Robin/Kosmos),
Coscientist, Sakana AI Scientist, SciAgents, BioDiscoveryAgent, Aviary/PaperQA.
What did they validate, and how honest were the validations? **[SEED-RISK]** — if
a general co-scientist already does this, the domain-specific system must justify
itself on evidence quality, not on generation.
E5. What is the *phage-specific* AI landscape in 2025–2026 — phage foundation
models, genomic LLMs (Evo/Evo 2, gLM), PLM-based phage tools?

## Thread F — Does machine-driven discovery actually work? (evidence, not vibes)

F1. Protein structure/design: what genuinely transferred (AF2/3, ESMFold,
RFdiffusion, ProteinMPNN) and what were the failure modes (orphan proteins,
conformational ensembles, disorder, complexes)?
F2. Drug repurposing: what fraction of ML-generated repurposing hypotheses ever
survived a trial? Name the successes and the base rate.
F3. Materials: GNoME's 2.2M "new materials" and the A-Lab — what did the critical
follow-ups (Cheetham & Seshadri; Leeman et al.) actually establish? **[This is the
cautionary tale most structurally similar to what is being proposed here.]**
F4. Literature-based discovery: Swanson/ARROWSMITH; the Tshitoyan et al. word2vec
"latent knowledge" claim and its replication status.
F5. ML-for-science reproducibility: what does the leakage literature (Kapoor &
Narayanan) say about pipelines exactly like the one being proposed?
F6. The single most on-point precedent: the Google AI co-scientist / Penadés
cf-PICI tail-piracy result. What was the claim, what was the control, and what did
skeptics say? Was the answer latent in prior literature?

## Thread G — Calibration, falsification, and not lying to the lab

G1. How do you measure the calibration of a hypothesis generator when the lab only
ever tests the top of the queue? (Selection bias makes the observed hit rate
uninterpretable.) What does the active-learning / bandit literature prescribe?
G2. What is the right retrospective benchmark? Can I freeze the corpus at a past
date and require the system to rediscover findings published after that date, with
a pre-registered scoring rubric?
G3. What formal machinery exists for "cheapest discriminating experiment" —
expected value of sample information, Bayesian experimental design, expected
information gain? What is actually computable here?
G4. How do you check novelty without an LLM hallucinating either the citation or
its absence? What is the failure mode of "I searched and found nothing"?
G5. Does adversarial self-critique (debate, critic agents, self-refine) measurably
improve factuality, or does it mostly produce confident-sounding revisions? What
is the actual evidence?
G6. What multiple-testing regime applies when N detectors × M genomes × K
hypotheses are generated? Hierarchical FDR? And how do you stop the LLM layer from
laundering a non-significant statistic into a plausible narrative?

## Thread H — The lab loop

H1. What experiments are cheap and fast for a phage lab (spot assay, EOP,
adsorption assay, capsule/LPS mutants, plasmid-borne defense in E. coli,
escape-mutant sequencing, cloning + heterologous expression, phage engineering)?
What do they cost, in days and dollars, and what do they discriminate?
H2. What is the natural cadence of hypothesis → experiment → data → model update,
and what latency makes the loop useless?
H3. What is the failure mode where the lab stops trusting the queue, and how many
bad hypotheses does it take?

## Thread I — Things the seed did not mention that might be better

I1. Nucleotide modification systems (5hmC, 7-deazaguanine/dpd, ADG) as a
restriction-escape predictor — underexploited?
I2. Stop-codon recoding and small-ORF blindness as a source of *systematic*,
correctable dark matter.
I3. Phage satellites, PICIs, and molecular piracy.
I4. Mosaicism/module-boundary detection → phage engineering targets.
I5. Selection signal (dN/dS, hypervariability) on hypothetical genes as a
"this is functional and under conflict" detector.
I6. The lab's own EOP residuals as a discovery instrument: pairs that are
*unexplained* by every known determinant are the most valuable rows in the matrix.
I7. Jumbo phages / nucleus-forming phages as a defense-agnostic chassis.
I8. Whether the highest-value output is actually a *ranked project queue* at all,
versus a smaller number of deeply-worked, pre-registered projects.

## Thread J — Questions for the client [BLOCKING]

J1. Which host species? (Determines defense landscape, receptor biology, public
corpus size, and whether capsule typing tools apply.)
J2. Bank size and matrix density/quantitativeness — how many phages, how many
strains, what fraction of cells filled, EOP values or binary?
J3. Wet-lab throughput per quarter and which assays are routine — this defines the
"cheapest discriminating experiment" library and therefore the ranking utility.
J4. Objective function: therapeutic candidate selection, basic-science
publications, or engineerable IP? These rank hypotheses differently.
