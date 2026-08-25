# Staying honest: calibration and falsification machinery

> A system that generates plausible hypotheses at scale is worse than useless.
> The difficulty is not generation — it is knowing which hypotheses deserve a
> bench experiment.

Everything here is designed against one observation from Phase 0: **GNoME and
A-Lab did not fail at generation. They failed at novelty verification and at
knowing what they had.** That is the failure mode to engineer against.

---

## 1. Null-first detector certification

No sensor may run on real data until it has been certified on data where the
signal it claims to detect has been destroyed. Permutations must be
**structure-matched**, or they certify nothing:

| Sensor | Correct permutation | Wrong permutation (would falsely pass) |
|---|---|---|
| Coupling (A3) | Permute host defense labels **within clade** | Global label shuffle — destroys phylogeny, inflates significance |
| Hotspot occupancy (A1) | Permute occupant identity **within cassette class** | Shuffle genome-wide — destroys positional structure |
| Small-ORF (A2) | Shuffle codons preserving GC and length distribution | Uniform random sequence — trivially distinguishable |
| Selection (A5) | Simulate under recombination at the observed rate | Assume no recombination — the classic dN/dS inflation |

A detector that cannot beat its own structure-matched null is deleted, not tuned.

---

## 2. End-to-end pipeline FDR (the null corpus) — *the central idea*

Per-test FDR control is necessary and badly insufficient here, because the LLM
layer is not a statistical test and does not carry a p-value. So:

> **Run the entire stack — sensors, gauntlet, synthesis, critique, ranking — on a
> corpus whose biological signal has been destroyed but whose structure has been
> preserved. Count how many compelling, well-argued, falsifiable hypotheses come
> out the other end. That number is the pipeline's false-discovery floor.**

If the null corpus yields 40 exciting hypotheses and the real corpus yields 60,
the system has produced roughly 20 findings and a great deal of prose.

This is cheap, it is run every release, and **it is the number I would put at the
top of the dashboard.** I have not seen it done in any AI-for-science system I
found in Phase 0, and I believe its absence is why the materials-discovery results
were able to fail so publicly: nobody ever asked their pipeline what it produced
from noise.

Construction: shuffle gene-family↔genome assignments within clade and within
cassette class, preserving phylogeny, GC, length distributions, and corpus
redundancy. The result should look exactly like real data to every summary
statistic except the ones carrying biological signal.

---

## 3. Phylogenetic non-independence as a first-class citizen

This is the load-bearing technical point in the whole design.

Bacterial genomes are related by descent, and the two documented confounders in
bacterial association studies — **genome-wide LD interrupted by recombination**,
and **strong population structure from clonal expansion** — both generate
type-I errors by hitchhiking. A correlation across 400 genomes drawn from three
clonal expansions has an effective n near 3.

Therefore: every association carries a **lineage-corrected effect size** and an
**`n_effective`**, and both the corrected and uncorrected values are reported.
Correction methods draw on the microbial-GWAS tradition (lineage-effect/LMM
handling as in pyseer; phylogeny-aware profiling as in PhyloCorrelate).

**The honest expectation:** correction will kill most of the coupling signal. That
is not a failure of the system, it is the system working. If it kills *all* of it,
mechanism A3 is dead and I say so in the risk register.

---

## 4. Absence requires a power calculation

Given INPHARED's *declining novel fraction*, an absence claim without a sampling
model is not a claim. Every B6 output must carry: given the sampling density of
this clade and the detection sensitivity, what is P(observe zero | the thing
exists at frequency f)? An absence claim with no power statement is rejected at
the schema level, like a claim with no falsifier.

---

## 5. The retrospective rediscovery benchmark

The most convincing validation available, fully achievable on public data, and it
should be built **before** the LLM layer, not after.

**Design.** Freeze the corpus at **2021-01-01**. Require the system to surface
discoveries published afterwards, with the rubric pre-registered before the run:

| Target | Published | Why it's a good test |
|---|---|---|
| PARIS | 2022 | Found by hotspot mining — the direct analogue of A1 |
| Tad1 (anti-Thoeris) | 2022 | Small protein; tests A2 |
| Acb2 (anti-CBASS) | 2023 | Small protein; tests A2 + A4 |
| The 21 non-island defense systems | 2022 | Tests the blind spot A1/A2 are built for |
| EVADES additions | 2024–25 | Volume test for recall@k |

Metrics: recall@k, and the **rank** of each known-later discovery. A system that
buries PARIS at rank 4,000 has not rediscovered it.

**The caveat that must not be buried — and it is serious.** ESM-2, ProstT5 and
every frontier LLM were trained on data that includes post-2021 sequences and
literature. **The benchmark is contaminated by construction.** Mitigations, in
descending order of rigour:

1. Run the strict version using **only non-PLM detectors** (positional,
   phylogenetic, selection) — these are computed from the frozen corpus alone and
   are genuinely clean.
2. Report PLM-assisted results separately and label them **contaminated**.
3. Never report a single headline number that mixes the two.

I would rather report a clean recall of 0.3 than a contaminated 0.8. **Anyone
reporting an uncontaminated-looking number on a PLM-based system over a temporal
holdout is either doing something subtle or is wrong**, and the distinction should
be visible in the methods.

---

## 6. Pre-registration, enforced by the schema

Every hypothesis is frozen and hash-committed *before* its experiment, carrying:
claim, falsifier, predicted probability of confirmation, the specific result that
kills it, and the assay. Post-hoc reinterpretation is structurally impossible
because the pre-registration hash is what the result is joined against.

The schema rejects a claim with a null `falsifier`. This single constraint does
more work than any amount of process discipline: **a hypothesis nobody can kill
cannot be stored.**

---

## 7. The calibration budget — the recommendation most likely to be resisted

**If the lab only ever tests the top of the queue, the calibration curve is
unidentifiable.** You will measure the hit rate of your top bin and learn nothing
about whether your 0.3 predictions are really 0.3 — which is precisely what you
need to know to trust the ranking.

> **Recommendation: allocate 20–25% of experimental capacity to a stratified
> random sample across confidence bins, including bins you expect to fail.**

This is a real cost and I will not disguise it. It is the price of knowing whether
the machine works. A lab that spends 100% of its bench time on the top of the
queue is running an unfalsifiable system, however good its statistics look.

The stratified sample also detects the specific failure where the model is
*ranking* well but *calibrated* terribly — common, and invisible to top-bin-only
testing.

---

## 8. Live calibration, and automatic quarantine

Track per-detector Brier score and calibration curves against realised outcomes.
**If a detector's realised hit rate falls outside its predicted band for n
consecutive experiments, quarantine it automatically** — remove it from the queue
and flag for audit. Not a review meeting; a code path.

The dashboard leads with three numbers, in this order:

1. **Null-corpus yield** (§2) — the false-discovery floor.
2. **Realised hit rate vs. predicted**, per confidence bin.
3. **Rediscovery count** — how often the prior-art gate fired T1. A rising T1 rate
   means the field is outrunning the loop, which is a strategic signal, not a bug.

---

## 9. The external anchor: 44.7%

DefensePredictor cloned 94 candidates and 42 worked. **That is the bar**, and it
is a good one because it comes from the same assay class this system would use.

- Substantially above ~45% in the top bin: **be suspicious**, and check for leakage
  and for a prior-art gate that is failing to fire.
- Far below: the system is not competitive with published practice and should not
  be consuming bench time.

Publishing this number as the target *before* the first experiment is itself a
falsification device — it makes the system's own success criterion pre-registered.

---

## 10. Conformal prediction, with the caveat that actually matters

Conformal prediction gives distribution-free finite-sample coverage **under
exchangeability**. Our data is phylogenetically structured and temporally shifting,
so **exchangeability fails exactly where we need the guarantee.**

Use it anyway — but validate coverage empirically on **held-out clades** and on
**held-out time periods**, and report the observed coverage rather than the
nominal one. A nominal 95% interval that achieves 60% coverage on a held-out clade
is worse than no interval, because it launders miscalibration as rigour.

---

## 11. The "so what" gate

A calibrated, true, novel finding that changes nothing is still waste. Every queue
item carries an impact statement: what becomes possible if this is true, and what
follow-on work it unlocks. Items that cannot answer are archived as *correct and
uninteresting* — a category the system should be willing to name out loud.

---

## 12. Where this machinery is still weak — stated because it is

Four honest gaps:

1. **The Synthesiser can still write a beautiful mechanism for a surviving-but-wrong
   observation.** The gauntlet reduces this; it does not eliminate it. Narrative
   fluency is not something you can fully engineer around.
2. **The prior-art gate depends on retrieval quality.** A claim novel only because
   the retriever missed a 2011 paper in a low-visibility journal will pass. The
   nearest-neighbour requirement (§ architecture 6) mitigates but does not solve it.
3. **Benchmark contamination is not fully solvable** with current PLMs. The clean
   subset is smaller and weaker than the contaminated one.
4. **Calibration on rare-but-important events is data-starved by definition.** The
   most valuable hypotheses are the ones we have the least evidence about, and no
   amount of machinery repairs that.
