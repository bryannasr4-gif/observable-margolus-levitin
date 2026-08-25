# The unmet need, what not to build, and the discovery mechanisms

## 1. What is actually unmet

Strip away everything Phase 0 showed is solved and one gap remains:

> **Nobody converts the phage-side arms-race gene pool into target-assigned,
> calibrated, falsifiable experiments.**

Every adjacent piece exists. Hotspots are annotated (iLund4u). Defense systems are
found (DefenseFinder, DefensePredictor, Mordret). Anti-defense is catalogued at
small scale (dbAPIS, EVADES, AntiDefenseFinder). Structures are predicted (BFVD,
Viro3D). Function is transferred (phold, PLM annotation). What nobody does is take
an unannotated 140-amino-acid ORF sitting in a variable cassette of a phage genome
and say:

> *"This inhibits **Thoeris**, specifically. Here is the evidence chain. Here is
> the 8-day experiment. Here is the result that kills it. We think it's 35%
> likely, and here is why you should believe our 35%."*

**Target assignment is the whole thing.** "Is this an anti-defense gene?" needs a
panel of defense-carrying strains and weeks. "Does this inhibit Thoeris?" needs one
plate. The difference between those two questions is the difference between a
research programme and a wish list.

## 2. What NOT to build — explicitly

I would refuse to build each of these, and the charter asked me to say so.

| Don't build | Because | Use instead |
|---|---|---|
| **Bacterial defense-system discovery from PLMs** | Published twice in *Science* in the same month (DefensePredictor; Mordret et al.) | Consume their outputs as *host features* |
| **A phage annotation pipeline** | Pharokka + phold already annotate >50% of genes with 1.36M structures | Run phold; extend it only for recoding and small ORFs |
| **Genus-level host prediction** | iPHoP is integrated, benchmarked, <10% FDR | Call iPHoP; treat its output as a prior, not a fact |
| **Strain-level EOP prediction** | Ceiling is AUC 0.818 in *Klebsiella*, the easiest system in the field | Predict *mechanism*, not phenotype |
| **A protein structure database** | BFVD (351k) and Viro3D exist; ColabFold is ~1,000 structures/GPU-day | Predict structures only for the top ~10³ candidates |
| **Taxonomy / clustering tooling** | vConTACT3, VIRIDIC, taxmyPHAGE, 95%/70% ANI thresholds are settled | Import |
| **A knowledge graph as the product** | KGs are infrastructure that photograph well and decide nothing. The output is a queue | Build the provenance DAG the queue needs, nothing more |
| **A fine-tuned phage LLM** | No evidence it beats retrieval + a frontier model on this task; freezes to a corpus snapshot | Retrieval over a versioned corpus |
| **Generative phage genome design** | Evo 2 fails long-range organisation, k-mer and constraint structure; megaDNA is functionally unvalidated | Not yet. Revisit in 2 years |
| **A multi-agent debate society** | ~75% of multi-agent failures are silent; scaling agents often increases failure | 3 tool-bound roles with genuine information asymmetry |
| **Training on the lab's EOP matrix** | Too small; and spot/EOP disagreement means the labels are heterogeneous | Use it as held-out truth and for residual analysis |

**One more, and it is the seed's own headline idea:** do not build ABSENCE
detection as a primary mechanism. In a corpus whose novel fraction is *falling*
while its size doubles, "conspicuously missing" is mostly "not yet sequenced."
It survives only as a conditional mechanism behind a sampling-power model (B6).

## 3. The discovery mechanisms

Graded honestly. **A = build**, **B = conditional**, **C = rejected**. For each:
what it finds, why humans miss it, data, achievability, falsifier, and prior-art
risk — because a mechanism whose prior art I haven't checked is a liability.

---

### TIER A — build these

#### A1. Hotspot occupancy → target assignment (*the arms-race address book*)

**What it finds.** Phage genomes have conserved syntenic anchors with small
(~1–5 kb) variable cassettes between them. Catalogue every occupant of every
cassette across the whole corpus, then **assign each occupant family a target
defense system** by testing whether that family's presence in phages tracks the
presence of a specific defense system in those phages' hosts — corrected for
phylogeny and sampling.

**Why humans miss it.** It requires holding thousands of phage genomes ×
hundreds of host defense repertoires × positional context simultaneously. No
person does this; the published instances (P2/P4) each took a lab-years.

**Data.** INPHARED + PhageScope + IMG/VR + PhagesDB genomes; DefenseFinder/PADLOC
on host genomes; iLund4u-style hotspot annotation; iPHoP/CRISPR-spacer/MVP host
links; EVADES/dbAPIS as positive controls.

**Achievable.** Yes. Detection uses existing tooling; the novel work is the
coupling statistic and its correction.

**Falsifier.** Clone the occupant ORF onto a low-copy plasmid in a host carrying
the predicted target defense system; challenge with a phage that system restricts.
**No ≥10-fold EOP rescue ⇒ dead.** (This is precisely the assay that validated
42/94 DefensePredictor candidates.)

**Prior-art risk: HIGH — and this is the honest position.** Hotspot detection is
built (iLund4u). Family-specific mining is published (Rousset, PARIS). A Feb-2026
preprint did P2/P4 at scale. [*Cell Host & Microbe* 2025 reports phages encode
"broad and specific counter-defense repertoires"](https://www.cell.com/cell-host-microbe/abstract/S1931-3128(25)00239-2)
— close enough that it must be read in full before committing.
**The defensible novelty is corpus-wide, cross-family, *target-assigned*,
calibrated triage — not hotspot discovery.** If full-text reading shows target
assignment is also done, A1 collapses to A3 + A2 and the programme should be
re-scoped. I would want that checked in week one.

---

#### A2. The small-protein blind spot (*a deliberately anti-biased detector*)

**What it finds.** Anti-defense proteins that every existing method structurally
cannot see. 80% of known anti-defense proteins are <200 aa; gene callers impose a
~100-codon floor; PLM embeddings and structure prediction both degrade on short
sequences; λ alone hides 55 unannotated translated ORFs.

**Why machines miss it.** Not an oversight — a *convention* baked into every
annotation pipeline, inherited by every downstream tool. The blind spot is
systematic, which means it is enumerable.

**Data.** Raw nucleotide genomes — **not GenBank CDS records**. This mechanism
requires re-calling ORFs from sequence with no length floor, recoding-aware, and
scoring candidates by cross-corpus conservation, positional context (A1), and
selection signal (A5) rather than by length.

**Achievable.** Yes, but it is the highest-false-positive mechanism in the set:
short ORFs occur by chance constantly. It is only viable *in conjunction with*
A1 and A5 — position and conservation carry the signal that length cannot.

**Falsifier.** Same heterologous assay; plus translation evidence via Ribo-seq
for the strongest claims.

**Prior-art risk: LOW-MEDIUM.** Small-protein discovery is an active field in
bacteria; applying it corpus-wide to phages *as an anti-defense enrichment
strategy* appears open.

**This is the mechanism I am most confident is genuinely unexploited**, precisely
because it is unglamorous: the gap exists because of a 100-codon convention, not
because of a hard biological problem.

---

#### A3. Phage-gene × host-defense coupling, done properly

**What it finds.** The seed's COUPLING mechanism — but the value is almost
entirely in the correction, not the correlation.

Bacterial and phage genomes are **not independent samples**. Clonal expansion and
genome-wide LD interrupted by recombination are the two documented confounders in
bacterial GWAS, and both produce type-I errors by hitchhiking. Uncorrected
co-occurrence over 10⁴ genome pairs will produce thousands of beautiful,
significant, meaningless associations.

**Method.** Lineage-corrected association (phylogenetic profiling in the
PhyloCorrelate tradition; LMM/lineage-effect handling as in pyseer; permutation
*within* clades to preserve structure). Report corrected and uncorrected effect
sizes side by side; **if correction kills it, it dies.**

**A warning from the literature I am taking seriously:** in defense–defense
analyses, mutual exclusivity was *not* strict, and the authors attributed it to
"genetic drift due to functional redundancy or selection against redundancy" —
not to incompatibility. **Exclusion has boring explanations.** Any exclusion
finding must clear that alternative before it is called mechanistic.

**Prior-art risk: MEDIUM on the bacterial side (done), LOW on the phage side.**

---

#### A4. Trigger–antagonist pairing (*second-order arms-race grammar*)

**What it finds.** PARIS is the template: a defense system whose *trigger* is a
specific phage counter-defense protein (T7 Ocr). That is a second-order
relationship — defense₂ senses anti-defense₁ — and it is a *grammar*, not a
one-off. Search systematically for defense systems whose distribution tracks a
specific phage protein family's distribution, conditioned on the first-order
defense that family already defeats.

**Why humans miss it.** It requires reasoning about two coupled layers at once
across the corpus. The field has only just started: *Nature Microbiology* 2025,
["A phage protein screen identifies triggers of the bacterial innate immune
system"](https://www.nature.com/articles/s41564-025-02239-6).

**Achievable: yes, but this is my most speculative Tier-A entry and I label it
as such.** It depends on A1/A3 working first, and on enough validated
anti-defense families to define first-order pairs.

**Falsifier.** Express the candidate antagonist in a host carrying the candidate
sensor; a true trigger should *induce* the defense phenotype (abortive infection /
growth arrest) in the absence of phage. Clean, cheap, and unambiguous.

**Prior-art risk: LOW.** Highest scientific upside in the whole catalogue.

---

#### A5. Selection-signal detection on hypotheticals

**What it finds.** "This gene is under conflict." Elevated dN/dS,
hypervariability, and positive selection on unannotated families is strong,
homology-independent evidence of arms-race participation. It answers "is this
functional and contested?" without needing to know what it does.

**Critical confounder, stated because it would otherwise sink this:** phage
genomes are **pervasively mosaic**, and recombination badly inflates dN/dS
estimates. Recombination-aware methods are mandatory, and any finding that
disappears under recombination correction dies.

**Prior-art risk: LOW-MEDIUM.** Hypervariability annotation exists; using it as
an *evidence channel in a triage system* is open.

---

### TIER B — conditional, later, or cheap-and-corrective

- **B1. Recoding & alternative-code sweep.** Corpus-wide re-annotation under
  tables 4/11/15. Cheap, corrective, publishable as a *resource*. Not discovery —
  and I will not dress it up as discovery.
- **B2. Nucleotide-modification-based restriction escape.** 180 phages/archaeal
  viruses encode 7-deazaguanine pathway enzymes, **60% over-represented among
  viruses of pathogens**. Predicts which phages should plate on RM⁺ hosts.
  Directly testable against an EOP matrix — an unusually clean prediction.
- **B3. DGR-driven host-range plasticity.** 92 known DGRs, temperate-biased,
  mutation rates to 1.38×10⁻² per base per generation. This predicts *which
  phages can evolve* new host range — a different and arguably more valuable
  question than which they currently infect. Directly relevant to therapy.
- **B4. Module-boundary / engineerability map.** Mosaic boundaries fall at gene
  and sometimes domain boundaries. A corpus-wide map of *where swaps are
  tolerated* is a phage-engineering asset.
- **B5. EOP residual analysis (needs the lab's matrix).** For each measured pair,
  regress out every known determinant — receptor locus, defense repertoire,
  modification state, RBP class. **The residuals are the most information-dense
  cells in the lab's data**: pairs that no known mechanism explains. Small-n, but
  each residual is a targeted, pre-localised discovery opportunity. This is the
  one place the lab's private data beats the public corpus outright.
- **B6. ABSENCE, behind a sampling model.** Only permitted with a power
  calculation: given sampling density in this clade, how surprising is zero? A
  bare "we never see X with Y" is not admissible.

---

### TIER C — rejected, with reasons

| Rejected | Why |
|---|---|
| End-to-end EOP matrix prediction | AUC 0.818 ceiling in the easiest species; four host groups is far harder |
| Bacterial defense discovery | Two *Science* 2026 papers |
| Generative phage design | Evo 2 fails long-range structure; no functional validation path |
| LLM-scored "novelty" as a ranking signal | LLM judges overvalue novel-*sounding* ideas and undervalue ideas that anticipate real research |
| Metagenomic dark matter as primary substrate | 8.7M vOTUs but only 231k high-quality representatives; fragment artifacts would dominate. Use isolate genomes as the spine, metagenomes for family expansion only |
| Automated literature-only hypothesis generation | This is the Co-Scientist mode. It works, it's published, and it is not what a lab with *data* should build |

## 4. One honest structural observation

The four seed mechanisms are not four mechanisms. **ABSENCE, DARK MATTER, and
COUPLING are three views of one object**: a sparse binary matrix of gene families
× genomes, read for holes, for unlabelled rows, and for correlated columns
respectively. CROSS-TAXON ANALOGY is genuinely different — it operates on
*literature*, not on the matrix, and it is the one job LLMs demonstrably do well.

That collapse is useful: it means one well-built, well-corrected
family-by-genome substrate serves three mechanisms, and the LLM layer sits
orthogonally on top rather than being wired into each.
