# Research log — provenance and access limits

**Access constraint (must be disclosed in the final report):** this session's
organization egress policy returns 403 on CONNECT for essentially all scholarly
hosts (nature.com, academic.oup.com, ncbi.nlm.nih.gov, pmc, europepmc.org,
biorxiv.org, arxiv.org, science.org, wikipedia.org, phagesdb.org). Only the
WebSearch tool and raw.githubusercontent.com are reachable.

Consequence: most findings below rest on **search-result synthesis, titles, and
abstracts**, not on full-text reading. Numbers are recorded with the source that
stated them. Any claim I could not corroborate from at least the abstract level
is marked LOW-CONFIDENCE. This is exactly the kind of provenance limitation the
system being designed is meant to make explicit, so it is stated up front rather
than buried.

## FINDING 1 (major, deflates seed mechanism #2 for the *bacterial* side)

Two Science 2026 papers already did "PLM embeddings -> discover new bacterial
immune systems", published back to back:

- **DeWeirdt et al., "DefensePredictor: a machine learning model to discover
  prokaryotic immune systems," Science 2026** (doi 10.1126/science.adv7924).
  ESM2 embeddings + gradient-boosting classifier. Applied to 69 diverse E. coli
  strains. **Cloned 94 predicted systems into a susceptible E. coli; 42 protected
  against >=1 of 24 phages => 44.7% bench validation rate.** 15 protein domains
  not previously validated as defensive. ~3,000 protein clusters with no homology
  to known systems across 1,000 prokaryotic genomes.
- **Mordret et al., "Protein and genomic language models uncover the unexplored
  diversity of bacterial immunity," Science 2026** (doi 10.1126/science.adv8275).
  >120 million bacterial proteins; hundreds of thousands of candidate antiphage
  families, many unannotated.

Two implications, in tension:
1. **DO NOT BUILD** a "PLM finds new bacterial defense systems" engine. Two groups
   published it simultaneously in Science. That lane is closed.
2. **44.7% is the empirical anchor for a calibrated discovery queue.** A
   state-of-the-art ML discovery pipeline, validated with a cheap binary assay,
   converts ~45% of its *top* candidates into bench hits. Any queue claiming
   better should be disbelieved; any queue delivering far worse is not competitive
   with published practice. This is the number to calibrate against.

**Confounder worth flagging:** the authors themselves note DefensePredictor
"appears to have learned to leverage the clustering of defense genes and the
association between defense genes and mobile elements." That is shortcut learning
on *genomic context*, not defense *function*. It works, but it structurally cannot
find defense systems outside defense islands / MGEs -- the exact blind spot the
guilt-by-association literature already had. A genuinely new engine must be
designed to look where context-based models are blind.

**The asymmetry that defines the opportunity:** the *bacterial* (defense) side is
now saturated with ML effort. The *phage* (anti-defense) side is not. EVADES
catalogues only 268 anti-defense proteins / 247 systems targeting 45 defense
families; dbAPIS has 41 experimentally validated families. AntiDefenseFinder
covers 156 systems. Compare with ~200 defense systems and >100k candidate defense
families. **Anti-defense is roughly two orders of magnitude less explored than
defense, and a phage lab is the natural owner of that side.**
