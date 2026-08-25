# An Autonomous Discovery Engine for Bacteriophage Biology — research design

Research-grade design study for Pride Laboratories. Phase 0 immersion plus a
proposed architecture, honesty machinery, loop, and roadmap.

**Read in order:**

| Doc | What it answers |
|---|---|
| [00 — Research agenda](research/00-research-agenda.md) | The questions I generated before reading anything, so they could be scored against what I found |
| [01 — Phase 0 immersion](research/01-phase0-immersion.md) | What I learned, what surprised me, and what contradicts the seed |
| [02 — Unmet need & mechanisms](research/02-unmet-need-and-mechanisms.md) | What NOT to build; the expanded mechanism catalogue, graded A/B/C |
| [03 — Architecture](research/03-architecture.md) | Layers, sensors, the confounder gauntlet, and why exactly three LLM roles |
| [04 — Honesty machinery](research/04-honesty-machinery.md) | Null-corpus FDR, phylogenetic correction, retrospective benchmark, calibration budget |
| [05 — Loop, roadmap, risks](research/05-loop-roadmap-risks.md) | Assay catalogue, phased plan with kill gates, red team |
| [notes/RESEARCH-LOG.md](research/notes/RESEARCH-LOG.md) | Provenance and access limits |

## The argument in six sentences

1. The public corpus doubled to ~29k complete phage genomes while its *novel*
   fraction **declined** — so "conspicuously absent" is mostly "not yet sequenced,"
   and the seed's headline mechanism has to be demoted.
2. A material share of phage "dark matter" is **annotation failure** — stop-codon
   recoding and a 100-codon gene-calling floor — not unknown biology; phage λ alone
   hides 55 unannotated translated ORFs.
3. Two *Science* 2026 papers just did "protein language models discover bacterial
   defense systems," so that lane is closed — **but one of them validated 42 of 94
   candidates at the bench, giving a 44.7% calibration anchor** for what a good
   discovery queue looks like.
4. The **phage side is ~100× less explored** than the bacterial side, for a
   methodological reason: 80% of anti-defense proteins are under 200 amino acids,
   where homology search, PLM embeddings and structure prediction all degrade.
5. So the engine should be a **phage-side, small-protein-biased, target-assigning**
   machine — not "is this anti-defense?" but "does this inhibit *Thoeris*?", because
   that turns a panel into a single plate.
6. And it stays honest by running its **entire pipeline on a null corpus** to
   measure how many compelling hypotheses it manufactures from noise, by treating
   **phylogenetic non-independence** as a first-class quantity, and by spending
   **20–25% of bench capacity on experiments it expects to fail** — the only way the
   calibration curve is identifiable at all.

## Provenance limit

Session egress policy blocked every scholarly host (403 at CONNECT). Findings rest
on search synthesis and abstracts, **not full-text reading**. Stated up front
because a document about intellectual honesty should not open by overstating its
own sourcing. See the research log.
