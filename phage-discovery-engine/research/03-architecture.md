# Architecture

## 0. Six principles, each earned from Phase 0

| # | Principle | Evidence it comes from |
|---|---|---|
| **P1** | **Statistics decide; LLMs explain.** No claim enters the queue on LLM judgment alone. | LLM novelty judgments diverge from expert gold; LLM judges overvalue novel-*sounding* ideas |
| **P2** | **Critics have tools, not opinions.** Every adversarial critic is bound to an executable check. A critic that can only "think harder" is deleted. | LLMs cannot reliably self-correct without external ground truth; self-critique can degrade performance |
| **P3** | **Fewest agents with genuinely different information.** Roles are justified by information asymmetry, never by role-play. | ~75% of multi-agent failures are silent; scaling agents often increases failure |
| **P4** | **Re-call, don't inherit.** Derive ORFs from nucleotides; never trust upstream CDS records. | λ has 55 unannotated translated ORFs; recoding breaks standard-code annotation |
| **P5** | **Every number carries its null.** Each detector declares a null model and is certified against it before it may emit. | Leakage documented across 17 fields / 294 papers |
| **P6** | **Frozen, content-addressed corpus snapshots.** A novelty claim is meaningless without a stated snapshot. | GNoME and A-Lab both failed *novelty verification*, not generation |

**The single most important structural decision follows from P1 and P2:**

> **The confounder gauntlet runs BEFORE any LLM sees an observation.**

Because once a language model sees a spurious correlation, it will write a
compelling mechanism for it. Narrative capacity is the enemy of calibration.
Nothing reaches the narrative layer until it has survived the statistics.

---

## 1. Layer diagram

```
 [6] QUEUE ── pre-registered claims, calibrated p, falsifier, costed experiment
       ▲
 [5] RANK ── declared utility; expected information gain per dollar
       ▲
 [4] ADVERSARIAL VERIFICATION ── tool-bound critics (prior-art / confounder / design)
       ▲
 [3] SYNTHESIS (LLM) ── observation → mechanistic claim + explicit falsifier
       ▲
 ═════ CONFOUNDER GAUNTLET ═══ ← no LLM may see anything below this line
       ▲
 [2] SENSORS ── statistical detectors emitting OBSERVATIONS (effect, null, q, provenance)
       ▲
 [1] REPRESENTATION ── re-called ORFs · embeddings · 3Di · synteny graph · phylogenies
       ▲
 [0] SUBSTRATE ── versioned, content-addressed corpus snapshots
```

Data flows up. **Provenance flows down**: every queue item resolves to a DAG of
observations → representations → substrate objects, each content-addressed.

---

## 2. Layer 0 — Substrate

| Domain | Sources | Scale | Notes |
|---|---|---|---|
| Phage isolate genomes | INPHARED, PhageScope, PhagesDB | ~29k complete; 874k sequences; 6k actinophage | **The spine.** Curated, complete, low artifact rate |
| Metagenomic viral | IMG/VR v4 | 15.6M UViGs / 8.7M vOTUs | **Family expansion only** — 231k high-quality reps; fragments would dominate otherwise |
| Host genomes | RefSeq, GTDB | — | Needed for defense repertoires and phylogeny |
| Defense annotations | DefenseFinder, PADLOC, + DefensePredictor/Mordret outputs | ~200 systems; 10⁵ candidate families | Consumed, never rebuilt |
| Anti-defense | dbAPIS, EVADES, AntiDefenseFinder | 41 / 268 / 156 | **Positive-control set and benchmark truth** |
| Structures | BFVD, Viro3D, AFDB | 351k viral | >62% of BFVD has no close structural neighbour |
| Families | PHROGs | 38,880 | Remote-homology backbone |
| Interactions | PhagesDB host data, MVP, lab EOP | 26.5k inferred; lab matrix | **Mostly *inferred*, not experimental — see §6** |
| Literature | Europe PMC OA | ~9.3M full texts, 43M abstracts | Bulk download permitted for the OA subset |

Each snapshot is content-addressed and dated. **A claim's novelty is asserted
against a named snapshot or not at all** (P6).

---

## 3. Layer 1 — Representation

- **ORF re-calling** (P4): no length floor; genetic-code detection (tables 4/11/15)
  before calling; start-site correction. Output is a *superset* of GenBank CDS
  with a per-ORF support score. This is the substrate for A2.
- **Sequence embeddings**: ESM-2 650M / ESM-C 600M. Deliberately not the largest
  model — medium PLMs match large ones on realistic transfer tasks, and the
  8M/35M tiers are too weak. Cost is bounded (§7).
- **Structure proxy**: ProstT5 → 3Di tokens → Foldseek, the phold approach. Full
  ColabFold reserved for the top ~10³ candidates only (~1,000 structures/GPU-day).
- **Synteny graph**: conserved anchors and variable cassettes per phage family;
  the substrate for A1.
- **Phylogenies**: host and phage trees, required by the gauntlet — not optional.

---

## 4. Layer 2 — Sensors

Sensors emit **observations**, never hypotheses. An observation is:

```
Observation {
  detector_id, corpus_snapshot
  claim_frame          // e.g. "family F occupies cassette C in clade K"
  effect_size, ci
  null_model_id, p_raw, q_corrected
  n_effective          // after phylogenetic deflation, NOT raw n
  provenance []        // content-addressed
  known_confounds []   // gauntlet results, attached
}
```

Six sensors, mapping to the mechanisms: `hotspot-occupancy` (A1),
`coupling-lineage-corrected` (A3), `small-orf-conservation` (A2),
`selection-conflict` (A5), `trigger-pairing` (A4), `absence-power` (B6).

**No sensor may run on real data until it has been certified on its null** (P5).

---

## 5. The confounder gauntlet

Automated, non-negotiable, and it runs first. Six checks:

1. **Phylogenetic non-independence.** Lineage-corrected effect size and an
   *effective* n. A correlation across 400 genomes from three clonal expansions
   has n ≈ 3, not 400. Report corrected and uncorrected side by side; correction
   kills it or it survives.
2. **Sampling-bias / power.** Given sampling density in this clade, what is the
   probability of observing zero by chance? Gates every absence claim (B6).
3. **Database redundancy.** Near-identical genomes inflate every count. Dereplicate
   at declared ANI, and report both counts.
4. **Annotation-artifact check.** Is the "novel family" an artifact of a wrong
   genetic code, a split ORF, or a frameshift? Re-derive from nucleotides.
5. **Shortcut-learning check.** Can the association be reproduced from genomic
   context alone, with the protein sequence masked? If yes, it is DefensePredictor's
   documented failure mode and must be labelled as context-driven, not
   function-driven.
6. **Leakage check.** Do train and evaluation sets share homologues above
   threshold? Splits are by **family and by clade**, never at random.

---

## 6. Layers 3–4 — The LLM layer, and why exactly three roles

**Three roles, because there are exactly three genuine information asymmetries.**
Not five, not a debate society (P3).

| Role | Sees | Cannot see | Tools (mandatory — P2) |
|---|---|---|---|
| **Synthesiser** | Surviving observation + evidence chain + mechanism ontology | The literature; the prior-art verdict | Structure/domain lookup, defense-system ontology, homolog search |
| **Prior-art auditor** | Literature corpus + the *claim text only* | The claim's confidence, the evidence strength | Europe PMC retrieval, structured claim matching, EVADES/dbAPIS lookup |
| **Experiment designer** | Claim + falsifier + assay catalogue + costs | Novelty score, so it cannot be seduced by excitement | Assay catalogue, cost model, EOP power calculator |

The Synthesiser is blinded to the literature deliberately: if it knows the
prior-art verdict it will write toward it. The auditor is blinded to confidence so
it cannot rationalise a high-scoring claim. This is **blinding as an architectural
primitive**, borrowed from experimental design rather than from agent frameworks.

**Why not one LLM?** It would be Co-Scientist: literature recombination without
statistical grounding, and no separation between generating and criticising.

**Why not more?** Every additional agent adds a silent-failure surface, and 75% of
multi-agent failures are silent. Three roles is the minimum that achieves real
information asymmetry, and I would rather defend three than ten.

### The three-tier prior-art gate

| Tier | Finding | Action |
|---|---|---|
| **T1 — exact** | This claim is published | **Kill**, and *log it as a rediscovery* — rediscoveries are the calibration signal, not waste |
| **T2 — homologous** | Published for a homolog | Downgrade to confirmation; low priority unless the homolog is distant |
| **T3 — mechanistic analogy** | Mechanism known in another system | **Upgrade.** This is the cross-taxon analogy the seed wanted, and it is a *strength* |

**The failure mode to design against:** "I searched and found nothing" is weak
evidence of novelty and the auditor must never be allowed to return it bare. It is
required to state its search coverage and produce the **nearest published claim**
with a distance judgment. A verdict without a nearest neighbour is rejected by the
harness. This is the check that GNoME and A-Lab both lacked.

---

## 7. Cost

Compute is not the constraint; this surprised me and is worth stating plainly.

| Item | Estimate | Basis |
|---|---|---|
| ESM-2 650M embeddings, ~5–10M phage proteins | low single-digit GPU-weeks | ~5.6 GB GPU for 300–400 aa batches |
| ProstT5 3Di for the same set | Similar order; far cheaper than folding | phold does this at corpus scale |
| ColabFold on top 10³ candidates | ~1 GPU-day | ~1,000 structures/GPU-day |
| Foldseek / MMseqs2-GPU all-vs-all | Hours–days | GPU MMseqs2, *Nature Methods* 2025 |
| Literature auditing (LLM) | **The dominant recurring cost** | Every claim triggers retrieval + reading |
| Storage | Terabytes, not petabytes | Isolate spine, not full IMG/VR |

**The real costs are bench time and the calibration budget** (§ honesty doc, item
7), not GPUs. A design that spends its money on model size rather than on
falsification has misread this problem.

---

## 8. What the data model must carry that a normal one doesn't

Three unusual fields, each load-bearing:

1. **`n_effective`** alongside `n`. Phylogenetic deflation is not metadata — it is
   the number that decides whether a finding is real.
2. **`falsifier`** as a required, non-null field on every claim. A claim with no
   discriminating experiment cannot be stored, let alone ranked. Enforced by the
   schema, not by convention.
3. **`corpus_snapshot`** on every novelty assertion. Without it, "novel" is not a
   claim about the world; it is a claim about someone's cache.
