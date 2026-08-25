# Phase 0 — Immersion: what I found, what surprised me, and what contradicts the seed

**Provenance note, stated up front because this document is about honesty.** This
session's egress policy blocked every scholarly host (nature.com, science.org,
academic.oup.com, ncbi/PMC, europepmc.org, biorxiv.org, arxiv.org) with a 403 at
CONNECT. Only web *search* was reachable. Every finding below therefore rests on
search-result synthesis, titles, and abstract-level content — **not on full-text
reading**. Numbers are attributed to the source that stated them. Where a number
would change a design decision if wrong, I say so. A system built to stop people
from over-claiming should not open by over-claiming its own sourcing.

---

## 1. The corpus is growing, and getting less novel

INPHARED held 14,244 complete phage genomes in 2021 and 28,777 by 2026 — but the
five-year retrospective reports that **the proportion representing novel
species-level diversity has *declined*: redundant sequencing is outpacing new
discovery** ([INPHARED 5-year retrospective, 2026](https://www.biorxiv.org/content/10.64898/2026.05.06.722914v1);
[Cook et al., PHAGE 2021](https://pubmed.ncbi.nlm.nih.gov/36159887/)).

Metagenomics is a different scale entirely: IMG/VR v4 holds **15.6M uncultivated
viral genomes clustering into 8.7M vOTUs, of which only 231,408 have a
high-quality representative**
([Camargo et al., NAR 2023](https://academic.oup.com/nar/article/51/D1/D733/6833254)).
PhageScope aggregates 873,718 phage sequences
([Wang et al., NAR 2024](https://academic.oup.com/nar/article/52/D1/D756/7334092)).
PhagesDB holds 6,011 sequenced actinobacteriophages under one uniform protocol
([PhagesDB](https://academic.oup.com/bioinformatics/article/33/5/784/2731030)).

**Why this matters, and it is the most consequential thing in this document:**
the seed's first mechanism is ABSENCE — "what is conspicuously missing relative
to a learned model of expectation." A corpus whose novel fraction is shrinking
while its size doubles is a corpus dominated by *sampling bias*. In such a
corpus, "absent" and "never sequenced" are nearly indistinguishable. **I am
demoting ABSENCE from the seed's headline mechanism to a conditional one** — it
is only usable behind an explicit sampling model, and I say how below.

## 2. Dark matter is real, deeper than the seed said, and partly not biology at all

The seed said 50–70% of phage genes are hypothetical. The right numbers are
stratified, and the per-family figure is far worse than the per-gene figure:

- **~40–45% of proteins and ~75–85% of protein *families*** in curated viral
  databases are hypothetical or unknown; 40–90% of viral genes in environmental
  data lack known homologues
  ([Viral Dark Matter review, *Biochemistry* 2025](https://pubs.acs.org/doi/10.1021/acs.biochem.5c00349)).
- PHROGs: 38,880 remote-homology families over 868,340 proteins, roughly half of
  families assignable to a function
  ([Terzian et al., NARGAB 2021](https://academic.oup.com/nargab/article/3/3/lqab067/6342220)).
- Structure helps but does not close it. Phold (ProstT5 → 3Di → Foldseek against
  1.36M predicted phage structures) annotates **>50% of genes on an average
  phage**, beating sequence-homology tools especially on metagenomic phages
  ([Bouras et al., NAR 2026](https://academic.oup.com/nar/article/54/1/gkaf1448/8415830)).
- BFVD holds 351,242 predicted viral structures, and **>62% show no or low
  structural similarity to existing repositories**
  ([Kim et al., NAR 2025](https://academic.oup.com/nar/article/53/D1/D340/7906834)).
- PLM-based annotation expanded the annotated fraction of viral protein families
  in ocean virome data by **29%**
  ([Flamholz et al., *Nature Microbiology* 2024](https://www.nature.com/articles/s41564-023-01584-8)).

**SURPRISE, and a competing explanation the seed does not consider: a
material share of "dark matter" is annotation failure, not unknown biology.**

- **Stop-codon recoding.** Some lineages reassign TAG→Gln (table 15) or TGA→Trp
  (table 4). Lak megaphages annotated under the standard code show ~70% coding
  density against ~90% for normal phages. 76 INPHARED genomes and 712 UHGV vOTUs
  were predicted to recode; Pharokka-gv fixes it
  ([Peters et al., *ISME Communications* 2024](https://academic.oup.com/ismecommun/article/4/1/ycae079/7696150)).
- **Small ORFs.** The 100-codon convention systematically erases them. A
  ribosome-profiling study of **phage λ — the most studied phage on Earth —
  found 55 previously unannotated ORFs with translation evidence**
  ([smORFer / small-protein literature](https://pmc.ncbi.nlm.nih.gov/articles/PMC8421149/);
  [Storz et al., *Mol Syst Biol*](https://link.springer.com/article/10.15252/msb.20188290)).

If λ is under-annotated by 55 ORFs, the corpus-wide figure is not a rounding
error. **This is cheaper to fix than the biology is to discover, and it changes
what the first build phase should be.**

## 3. Host-range prediction has a hard ceiling, well below what a lab needs

- iPHoP: **genus level, <10% FDR** — and genus is explicitly its resolution limit
  ([Roux et al., *PLOS Biology* 2023](https://journals.plos.org/plosbiology/article?id=10.1371%2Fjournal.pbio.3002083)).
- PhageHostLearn, on *Klebsiella* — capsule-driven specificity, the most
  favourable system in all of phage biology, with receptor-binding proteins and
  K-locus typing available — reaches **ROC AUC 81.8%** at strain level, held in
  lab validation
  ([Boeckaerts et al., *Nature Communications* 2024](https://www.nature.com/articles/s41467-024-48675-6)).

**That is the ceiling, in the easiest case.** Across four host groups with LPS,
teichoic acid and mycolic acid surfaces, strain-level EOP prediction is not
attainable today. Any part of the seed that implicitly hopes to *predict the
matrix* should be abandoned. The system's job is mechanism, not prediction.

## 4. The measurement itself is noisy

Spot tests **systematically overestimate both host range and virulence, and do
not correlate with EOP**
([Mirzaei & Nilsson, *PLOS ONE* 2015](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0118557)).
Much published "host range" data is therefore measuring a different quantity from
the one a model would be trained to predict. This is label noise, not model
weakness, and it caps achievable performance independent of method.

## 5. The defense field is being industrialised — and two Science 2026 papers just took the obvious lane

This is the finding that most changes the design.

- **DeWeirdt et al., "DefensePredictor", *Science* 2026** — ESM2 embeddings +
  gradient boosting. Applied to 69 diverse *E. coli* strains. **94 predicted
  systems cloned into a susceptible host; 42 protected against ≥1 of 24 phages
  → a 44.7% bench validation rate**, with 15 protein domains never previously
  validated as defensive; ~3,000 no-homology protein clusters across 1,000
  genomes ([Science](https://www.science.org/doi/10.1126/science.adv7924)).
- **Mordret et al., *Science* 2026** — fine-tuned ESM2 for distant homology plus
  an ALBERT-architecture *genomic* language model for genomic context; **>120M
  proteins**, up to 99% precision / 92% recall on held-out knowns, >32,000
  genomes, ~1.5% of a bacterial genome devoted to defense, **>85% of predicted
  families uncharacterised**; validated 12 systems in *Escherichia* and
  *Streptomyces* and 6 more in Actinomycetota
  ([Science](https://www.science.org/doi/10.1126/science.adv8275)).

Two conclusions, in tension, and both matter:

1. **Do not build this.** "PLM finds new bacterial defense systems" was published
   twice, in *Science*, in the same month. That lane is closed.
2. **44.7% is the number to calibrate against.** A state-of-the-art discovery
   queue validated by a cheap binary assay converts roughly *forty-five percent*
   of its top candidates into bench hits. Any system claiming much better should
   be disbelieved; any system delivering far worse is not competitive with
   published practice. This is the most useful single number I found.

## 6. The asymmetry that defines the opportunity

The bacterial (defense) side is saturated with machine-learning effort. The phage
(anti-defense) side is not:

| Side | Scale of catalogue |
|---|---|
| Bacterial defense | ~200 known systems; **hundreds of thousands** of candidate families (Mordret) |
| Phage anti-defense | ~180 proteins known to inhibit defenses; **EVADES: 268 proteins / 247 systems** targeting 45 defense families; dbAPIS 41 validated families / 4,428 homologs; AntiDefenseFinder 156 systems, 47,981 instances |

Sources: [AntiDefenseFinder / Tesson et al., NAR 2025](https://academic.oup.com/nar/article/53/1/gkae1171/7919512);
[dbAPIS, NAR 2024](https://academic.oup.com/nar/article/52/D1/D419/7331021).

**Anti-defense is roughly two orders of magnitude less explored than defense —
and there is a mechanistic reason for the gap, not just a sociological one:
80% of anti-defense proteins are shorter than 200 amino acids**, which the
AntiDefenseFinder authors name explicitly as a bottleneck for *both* sequence-
and structure-based detection.

Short proteins are exactly where homology search, PLM embeddings, and structure
prediction all degrade. **The blind spot and the opportunity are the same object.**

## 7. Both Science-2026 models are blind in a documented place

DefensePredictor's own authors note it "appears to have learned to leverage the
clustering of defense genes and the association between defense genes and mobile
elements" — that is shortcut learning on *genomic context*, not defense function.
It cannot, by construction, find systems outside defense islands and MGEs.

That blind spot is not hypothetical. A **context-agnostic functional selection
across 71 diverse *E. coli* strains recovered 21 conserved defense systems, none
of which had ever been detected as enriched in defense islands**
([Vassallo et al., *Nature Microbiology* 2022](https://www.nature.com/articles/s41564-022-01219-4)).

**A new engine earns its existence by looking where the context models cannot.**

## 8. The positional signal is real, mechanistic — and already half-mined

Arms-race genes occupy small, conserved genomic *addresses*:

- P2-like phages and P4-like satellites carry hotspots of **~1–5 kb between two
  conserved genes**, which are large reservoirs of anti-phage systems; mining
  them is how **PARIS** was found — a defense system triggered specifically by
  the T7 anti-restriction protein **Ocr**
  ([Rousset et al., *Cell Host & Microbe* 2022](https://www.cell.com/cell-host-microbe/fulltext/S1931-3128(22)00104-4)).
- *P. aeruginosa* genomic islands encode **11 conserved hotspots** carrying both
  defence and anti-defence
  ([NARGAB 2025](https://academic.oup.com/nargab/article/7/4/lqaf148/8328389)).
- **iLund4u** now annotates hyper-variability hotspots across *millions* of
  sequences via proteome communities
  ([bioRxiv 2024](https://www.biorxiv.org/content/10.1101/2024.10.15.618418v1)).
- Phage pangenome synteny is characterised: across 3,425 actinobacteriophage
  genomes, core gene order is well defined and **most accessory genes localise to
  a few positions** ([Synteny and linkage decay, 2025](https://www.biorxiv.org/content/10.1101/2025.08.12.669904v1.full)).
- A Feb-2026 preprint characterised defence hotspots across **thousands** of
  P2-like phages and P4-like satellites.

**I red-teamed my own best idea here and it survived only in modified form.**
Hotspot *detection* is built. Hotspot *mining for a specific phage family* is
published. What is not built is **target assignment** — not "this ORF is an
anti-defense gene" but "this ORF inhibits *Thoeris specifically*". That
distinction is the whole design, because it is what makes the experiment one
plate instead of a panel.

## 9. Machine-driven discovery in adjacent fields: one honest success, two cautionary tales, one measurement

**The closest precedent** — Google's Co-Scientist proposing cf-PICI tail piracy,
matching an unpublished experimental result
([*Cell* 2025](https://www.cell.com/cell/fulltext/S0092-8674(25)00973-0);
[*Nature* 2026](https://www.nature.com/articles/s41586-026-10644-y)). The honest
reading, which the reporting supports, is that **the model assembled hints
scattered across existing microbiology and virology literature into a coherent
new hypothesis.** That is literature-based discovery. It is genuinely valuable,
and it precisely delimits what LLMs contribute: **recombination of published
fragments, not perception of new signal in data.**

**The cautionary tales:**
- **GNoME.** Cheetham & Seshadri found "scant evidence for compounds that fulfil
  the trifecta of novelty, credibility, and utility"; many proposals were trivial
  dopant variants or chemically implausible
  ([The Register summary](https://www.theregister.com/2024/04/11/google_deepmind_material_study/)).
- **A-Lab.** Leeman, Palgrave et al. concluded **no new materials were
  discovered**, citing "systematic errors all the way through" and mishandled
  compositional disorder ([ChemRxiv 2024](https://chemrxiv.org/doi/full/10.26434/chemrxiv-2024-5p9j4);
  [Chemistry World](https://www.chemistryworld.com/news/new-analysis-raises-doubts-over-autonomous-labs-materials-discoveries/4018791.article)).
  The *Nature* paper was later corrected
  ([C&EN, 2026](https://cen.acs.org/research-integrity/Nature-robot-chemist-paper-corrected/104/web/2026/01)).

Both failed on **novelty verification**, not on generation. That is exactly the
failure mode this system must be built against.

**The measurements that constrain architecture:**
- **Kosmos**: independent scientists judged **79.4–80% of statements in its
  reports accurate** — about **one statement in five is wrong**
  ([Edison Scientific](https://edisonscientific.com/news/announcing-kosmos)).
- **Kapoor & Narayanan**: leakage documented across **17 fields, 294 papers**
  ([*Patterns* 2023](https://www.cell.com/patterns/fulltext/S2666-3899(23)00159-9)).
- **LLM ideation**: with 100+ NLP researchers, LLM ideas were rated *more novel*
  than expert-written ones — but LLM novelty judgments diverge substantially from
  expert gold, and **LLM judges overvalue novel-sounding ideas while undervaluing
  ideas that anticipate real research**. Novelty rating ≠ validity.
- **LLM self-correction**: without external ground truth, self-critique does not
  reliably improve reasoning and can degrade it.
- **Multi-agent systems**: 14 documented failure modes; **~75% of failures are
  silent**; scaling agents often *increases* failure rate
  ([Why Do Multi-Agent LLM Systems Fail?](https://openreview.net/forum?id=fAjbYBmonr)).
- **Drug repurposing**, the field most like "generate many plausible hypotheses":
  Phase I → approval base rate **7.9%** across 12,728 compounds (2011–2020).
- **Generative genomic LMs are not ready**: Evo 2 (40B params, 9.3T nucleotides)
  fails to preserve long-range genomic organisation, k-mer composition and
  evolutionary constraint, converging to an averaged k-mer landscape
  ([limitations preprint, 2026](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12871140/)).

**One encouraging counter-datapoint:** BioDiscoveryAgent, an LLM agent designing
genetic perturbation experiments, beat ML baselines by **21% average hit-ratio
improvement over five experimental rounds** (46% on the harder task)
([arXiv 2405.17631](https://arxiv.org/abs/2405.17631)). LLMs *are* good at
choosing what to test next. That is a narrower and more defensible role than
"co-scientist".

---

## Scorecard against my own agenda

| Seed mechanism | Verdict |
|---|---|
| **ABSENCE** | **Demoted.** Declining novel fraction in a doubling corpus makes "absent" and "unsampled" nearly inseparable. Usable only behind an explicit sampling/power model. |
| **DARK MATTER STRUCTURE** | **Partly solved, partly mis-specified.** Clustering is done (PHROGs, BFVD, phold). A material share of the residue is annotation failure — recoding and small ORFs — not unknown biology. The *remaining* structural dark matter is real (>62% of BFVD). |
| **CROSS-TAXON ANALOGY** | **Upheld, and it is the LLM's genuine job.** The cf-PICI result is exactly this, and it worked. |
| **COUPLING** | **Upheld on the phage side, largely done on the bacterial side.** Defense–defense co-occurrence across 42,925 genomes is published. Phage-gene × host-defense coupling with proper phylogenetic correction is open. |

**What surprised me most, in order:**
1. The corpus is getting *less* novel as it grows. That inverts the premise that more data means more discoverable structure.
2. λ has 55 unannotated translated ORFs. Our best-characterised object is not characterised.
3. Anti-defense is ~100× less catalogued than defense, and the reason is a protein-length bottleneck — a *methodological* gap, not a biological one.
4. The published bench hit-rate for a state-of-the-art discovery queue is 44.7%, which is both far better than I expected and a hard bar.
5. The two most-publicised machine-discovery successes in materials both failed their novelty checks, and neither failed at generation.
