# Reproducibility — numerics for *A Margolus–Levitin speed limit for observables*

Every quantitative claim in the paper is pinned to a seeded script here. Conventions throughout:
`ℏ = 1`, ground energy `E0 = 0` (so `H ⪰ 0`), spectral spread `σ_A = (λ_max − λ_min)/2`,
`Δ = |⟨A(T)⟩ − ⟨A(0)⟩|`, and the working ratio `r := T(⟨H⟩−E0) σ_A² / Δ²` (the bound is `r ≥ C⋆`).

**The named constant**

    C⋆ = 1/(8K),  K = sup_{x>0}(1−cos x)/x = sin x⋆,  x⋆ = smallest positive root of tan(x/2)=x
    x⋆ = 2.3311223704144226…,  K = 0.7246113537767085…,  C⋆ = 0.1725062674611626…  (pinned to 52 digits)

## Quick start

```bash
python -m pip install numpy scipy mpmath sympy
python -u numerics/make_results.py         # regenerates results.json (constants, no-go, saturation, mixed, structured, battery)
python -u numerics/verify_refinements.py   # trade-off curve, bandwidth-resolved constant, two-observable bound
python -u numerics/make_figures.py         # regenerates paper/figs/*.pdf
```

`results.json` is the consolidated, machine-readable record of the pinned numbers.

## Claim → script map

| Paper location | Claim / quoted number | Script | Reproduced value |
|---|---|---|---|
| Eq. (main), §V | `C⋆ = 1/(8K) = 1/(8 sin x⋆) = x⋆/(8(1−cos x⋆)) = (1+x⋆²)/(16x⋆) = 0.172506267461…` (four closed forms agree to 52 digits) | `theory_core_v2.py`, `make_results.py` | all four agree to 3.3e-52 |
| Lemma (tangent line), Fig. 2 | `1−cos x ≤ Kx` for all `x≥0`, equality only at `0, x⋆`; `K = sin x⋆ = 0.724611…` | `theory_core_v2.py` (part B) | max of `(1−cos x)−Kx` = −1.7e-10 (numerically ≤ 0), touching 0 at x⋆ |
| §III (no-go) | `P_min(Δ)/Δ → 0`; log–log exponent of `P_min(Δ)` is 2, not 1 | `theory_core.py`, `explore_observable_ML.py`, `make_results.py` | slope → 2.017 (and 1.9998 independent); `P_min/Δ → 0` |
| §V (saturation) | free-phase two-level infimum ratio `→ C⋆` from above | `theory_core_v2.py` (part C), `make_results.py` | 0.172506… ; β=0 (phase-restricted) → `C_old = 0.185551…` |
| §V, Fig. 4 | global-optimizer floor over `d = 2..7`, no config below `C⋆` | `optimize_constant.py` | min 0.17250656 ≥ C⋆ (approached from above) |
| §V / §VIII | apparent sub-`C⋆` dips are float64 catastrophic cancellation at excited pop. `~1e-14`; vanish at 60 dps | `audit_d2.py` | d=2 point recovers to `C⋆ + 2.66e-12` at 60 dps |
| §VI (mixed) | full `C⋆` for arbitrary `ρ`; violation search, 0 violations | `ext_mixed_check.py`, `ext_mixed_purification_check.py`, `make_results.py` | min sampled ratio 0.176–0.192 ≥ C⋆; convexity + purification chains audited |
| §VII (linear recovery) | eigenvector-`A` bound `→ 1/(4K) = 2C⋆ = 0.345013…` | `struct_verify_independent.py`, `struct_linear.py`, `make_results.py` | inf `r_lin` = 0.345071 ≈ 1/(4K) |
| §VII (turning point ≠ eigenvector) | zero initial slope alone gives **no** positive linear constant | `struct_verify_independent.py` (§3b) | `r_lin` falls to ≈ 0.034–0.057 ≪ 1/(4K) |
| §VII (structured) | rank-1 ground coupling → `C⋆` (no improvement); real couplings → `C_old`; banded → `C⋆` | `struct_rank1_check.py`, `struct_multifreq.py`, `struct_fast.py`, `struct_adversary.py`, `struct_verify_independent.py` | rank-1 inf → 0.172506; real → 0.185551; banded → C⋆ |
| §VII (ground-gapped) | static floor `⟨H⟩−E0 ≥ (Δ_gap/2)Δ/σ_A`; `r_quad → ∞` as `ε→0` | `struct_gap_law.py`, `struct_linear.py` | worst slack +3.0e-3 (≥0); `r_quad` grows ~1/ε |
| §V, Fig. 5 (trade-off curve) | `P⋆(δ)` [Eq. (curve)] is a strict lower boundary in `d=2..6`; `P⋆→0` as `δ→0` (slope `4C⋆`), `P⋆→π/2` as `δ→1` | `verify_refinements.py` (part 1) | min slack `P − P⋆(δ)` = +1.4e-5 over 1.56M points; slope 0.6905 ≈ 4C⋆; endpoint 1.5708 = π/2 |
| §V, Eq. (bandwidth) | `C⋆(L)=1/(8K(L))`; gains `1.58×` at `L=1`, `2.96×` at `L=0.5`; violation search 0 violations | `verify_refinements.py` (part 2) | 1.5763× / 2.9596×; 0 violations in 8.0e5 trials |
| §VI, Eq. (multi) | two-observable `T(⟨H⟩−E0) ≥ C⋆(Δ₁²+Δ₂²)/σ_A²`; `2×` gain at `Δ₁=Δ₂` | `verify_refinements.py` (part 3) | min ratio 0.172649 ≥ C⋆, 0 violations in 4.0e6 trials; at `Δ₁≈Δ₂` ratio → 0.3469 ≈ 2C⋆ |
| §IX (battery) | saturation `R → C⋆`; hostile min `R = 0.218`, 0 violations | `app_battery_check.py`, `make_results.py` | R → 0.17264 → C⋆; hostile min 0.2181 over 4000 trials |

## Notes

- `theory_core.py` is the original first-pass script and reports the **phase-restricted** constant
  `C_old = 0.18555147…`; the free relative phase lowers it to the sharp `C⋆ = 0.17250627…`
  (see `theory_core_v2.py` for the correction). Both values, and their `7.56%` ratio, are intentional
  and documented — `C_old` is the `β=0` sub-optimum quoted in §V.
- All scripts are self-contained and seeded; expected wall time is seconds to a couple of minutes each
  (`verify_refinements.py` part 3 is the longest at ~a few minutes due to ~4M sampled configurations).
- Random-search "min ratio ≥ C⋆" results are *stress tests*, not proofs; the bound itself is proved
  analytically in the manuscript (`paper/main.tex`).
