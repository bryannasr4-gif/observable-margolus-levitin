# The human–machine loop, roadmap, risks, and what I rejected

## 1. The experiment catalogue

You said you didn't know which assays are routine, so I am specifying the menu
rather than assuming it. **Costs and durations below are engineering estimates
from published protocols, not cited figures** — I could not find published
cost-per-assay data, and I would rather label the estimate than dress it up.

| Assay | Time | Est. cost | Discriminates | Verdict |
|---|---|---|---|---|
| **Spot test panel** | Hours | ~$1–3/pair | Almost nothing | **Pre-screen only.** Systematically overestimates host range and does not correlate with EOP |
| **Quantitative EOP** | 1–2 d | ~$10–15/pair | Productive infection vs. killing | The real phenotype measurement |
| **Heterologous expression + EOP shift** | ~1 wk + synthesis | ~$50–150/gene | **Does gene X confer/defeat defense Y** | **The workhorse.** Everything in Tier A depends on it |
| **Adsorption assay** | 1 d | ~$10 | Adsorption vs. post-adsorption barrier | Cheap disambiguation, under-used |
| **Escape-mutant isolation + WGS** | 1–2 wk | ~$50–100/genome | Receptor identity; defense target | Highest information per dollar in the whole menu |
| **Tn-seq / RB-TnSeq** | 3–6 wk | High | Receptor + host factors, genome-wide | Identified receptors for **5 of 6** previously uncharacterised phages |
| **Ribo-seq** | Weeks | High | Small-ORF translation | The only way to close out A2 claims definitively |

**The dependency I need to flag plainly:** the heterologous expression assay is the
falsifier for A1, A2, A3 and A4 — that is, for the entire Tier-A programme. It is
exactly how DefensePredictor validated 42/94 candidates (low-copy plasmid,
*E. coli* MG1655, panel of 24 phages, hit = ≥10-fold EOP drop or smaller plaques)
and how the functional selections that found the 21 non-island systems worked.

**If the lab cannot clone and express genes in a tractable host, Tier A is not
falsifiable in-house.** That is not a reason to abandon it — it is a reason to
partner, or to re-scope toward Tier B, which leans on EOP and sequencing you
already do. I would want to know this before Phase 3, not after.

## 2. The loop

```
   observation → claim (+falsifier, +predicted p) → PRE-REGISTER (hashed)
        ↑                                                    ↓
   recalibrate ← outcome joined on hash ← experiment ← rank by EIG/$
        │
        └── T1 rediscovery? → log to calibration set, not to waste
```

Two properties matter more than the diagram:

- **Latency.** If hypothesis→result exceeds roughly one quarter, the model
  updates too slowly to steer and the queue becomes a static report. Design target:
  **≤6 weeks** for the workhorse assay.
- **Negative results are first-class.** A killed hypothesis updates the detector
  that produced it. Under pre-registration, a negative is as informative as a
  positive, and the system must not be allowed to quietly drop them — the schema
  joins outcomes on the pre-registration hash precisely so that it can't.

**Where the LLM genuinely earns its place in this loop:** choosing what to test
next. BioDiscoveryAgent improved hit ratio by 21% over ML baselines across five
experimental rounds (46% on the harder task). That is a narrow, measured,
defensible claim — and a much better description of the role than "co-scientist".

## 3. Roadmap — each phase independently useful

**Phase 1 (0–3 mo) — Substrate and re-annotation.**
Corpus snapshots; recoding-aware, no-floor ORF re-calling; synteny/cassette atlas.
*Independently useful even if everything downstream fails:* a corrected,
small-ORF-inclusive re-annotation of the phage corpus is a publishable community
resource. Given λ's 55 missing ORFs, the corpus-wide correction is a result in
itself.
**Gate:** measurable increase in coding density and in recovered ORFs vs. GenBank.

**Phase 2 (3–6 mo) — Sensors, gauntlet, retrospective benchmark. This is the
go/no-go.**
Build detectors, certify each on structure-matched nulls, build the confounder
gauntlet, then run the 2021-frozen rediscovery benchmark **using only non-PLM
detectors** for the clean number.
**Gate:** clean recall of known-later discoveries meaningfully above chance, with
PARIS and the non-island systems ranked highly. *If this fails, stop.* No amount of
LLM layer rescues a sensor stack that cannot rediscover known biology.

**Phase 3 (6–9 mo) — Synthesis, critique, queue.**
Three LLM roles, three-tier prior-art gate, pre-registration schema.
**Gate:** null-corpus yield (§ honesty 2) is a small fraction of real-corpus yield.
If noise produces as many compelling hypotheses as data, stop.

**Phase 4 (9–15 mo) — Bench loop and calibration.**
First pre-registered experiments, including the 20–25% stratified calibration
sample.
**Gate:** realised top-bin hit rate within a defensible distance of the 44.7%
anchor, and a calibration curve that is not wildly off-diagonal.

**Phase 5 — Expansion.** Tier B mechanisms; the lab's EOP residuals (B5); target
assignment extended to new host groups.

## 4. Risks and abandon criteria

| Risk | Signal | Response |
|---|---|---|
| **Sensor stack has no real signal** | Rediscovery benchmark at chance (clean subset) | **Abandon.** This is the honest kill switch |
| **Pipeline manufactures narrative** | Null-corpus yield ≈ real-corpus yield | **Abandon** or rebuild the gauntlet |
| **Phylogeny eats everything** | Correction kills essentially all coupling | Drop A3; A1/A2/A5 survive independently |
| **A1 already published** | Full-text reading of the 2025 CHM counter-defense paper and the Feb-2026 P2/P4 preprint shows target assignment is done | Re-scope to A2/A4. **Check this in week one** |
| **No molecular capability in-house** | — | Partner, or restrict to Tier B |
| **Contamination invalidates the benchmark** | Clean and contaminated numbers diverge wildly | Report clean only; accept a weaker claim |
| **A2 drowns in false positives** | Short-ORF candidates dominate the queue | Require positional + selection evidence as a hard conjunction, not a score |
| **Prior-art gate keeps firing T1** | Rising rediscovery rate | The field is outrunning the loop — narrow the niche |
| **Lab trust collapses** | 2–3 confidently-ranked hypotheses fail publicly | Fixed by the calibration budget *before* it happens, not after. Trust is lost in bulk and regained slowly |

**The risk I rate highest is not technical.** It is that the system produces 200
well-formatted, well-cited, falsifiable hypotheses and the lab has bandwidth for
six. A queue longer than throughput is a queue nobody reads. **The system should
be tuned to emit roughly 2–3× the lab's actual experimental capacity, and no
more** — deliberately throwing away real findings to protect the credibility of the
ones it surfaces. That is an uncomfortable design choice and I think it is correct.

## 5. Red team: what I rejected, including my own ideas

**Rejected from the seed:**
- **ABSENCE as a headline mechanism.** Demoted to conditional (B6). A corpus with a
  falling novel fraction cannot distinguish absent from unsampled.
- **"Dark matter clustering" as novel.** Largely done — PHROGs, BFVD, phold. And a
  material share of the residue is annotation failure, not biology.
- **The four mechanisms as four things.** ABSENCE, DARK MATTER and COUPLING are
  three readings of one family×genome matrix.

**Rejected from my own first drafts:**
- **A phage foundation model.** Attractive, expensive, and Evo 2's failure on
  long-range genomic structure says the generative direction isn't ready. Killed.
- **A learned EOP predictor over the lab's matrix.** The matrix is too small, the
  labels are heterogeneous (spot vs. EOP), and the published ceiling is AUC 0.818
  in an easier setting. Killed and replaced with residual analysis (B5).
- **A debate/consensus multi-agent layer.** Killed on evidence: ~75% of multi-agent
  failures are silent and scaling agents often increases failure.
- **LLM-scored novelty as a ranking term.** Killed: LLM judges overvalue
  novel-*sounding* ideas and undervalue ideas that anticipate real research. Novelty
  enters only through the retrieval-grounded three-tier gate.
- **My own best mechanism, in its first form.** "Find arms-race hotspots" is built
  (iLund4u) and mined for P2/P4 (Rousset; Feb-2026). A1 survives *only* as target
  assignment plus calibrated triage. I would rather narrow it now than discover it
  in month five.

**What I would abandon the whole approach for:** a clean retrospective benchmark at
chance. Everything else here is recoverable; that isn't.

## 6. The one thing I would build first if forced to pick

**The recoding-aware, no-length-floor re-annotation of the phage corpus (Phase 1),
scored against λ's 55 known-missing ORFs as a positive control.**

It is cheap, it is verifiable against ground truth that already exists, it is
useful to the community whether or not the rest is built, and it produces exactly
the substrate that the anti-defense mechanisms need — because 80% of anti-defense
proteins are shorter than the floor that every current pipeline imposes.

If that single artifact is all this project ever ships, it will still have been
worth building.
