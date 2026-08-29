# A Margolus–Levitin speed limit for observables

**Mean energy bounds the change of an expectation value — quadratically.**

[![arXiv](https://img.shields.io/badge/arXiv-2608.22658-b31b1b.svg)](https://arxiv.org/abs/2608.22658)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22165921.svg)](https://doi.org/10.5281/zenodo.22165921)
[![build-paper](https://github.com/bryannasr4-gif/observable-margolus-levitin/actions/workflows/build-paper.yml/badge.svg)](https://github.com/bryannasr4-gif/observable-margolus-levitin/actions/workflows/build-paper.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Paper: RevTeX](https://img.shields.io/badge/paper-RevTeX%204--2-informational.svg)](paper/main.tex)
[![Status: preprint](https://img.shields.io/badge/status-preprint-orange.svg)](#status)

This repository contains the manuscript and the complete, seeded numerical verification suite for a
quantum speed limit on **observables** powered by the **mean energy**.

> **What is new, and what is inherited.** New here: the **no-go theorem** — the optimal
> state-independent mean-energy exponent of $\Delta$ is exactly **two**, so no linear mean-energy bound
> exists — and the **eigenvector linear-recovery theorem**. Inherited, and cited as such in the manuscript: the
> constant $C_\star$ (the $q\to0$ Giovannetti–Lloyd–Maccone member), the tangent-line lemma, the proof
> route (Jones–Kok, corrected by Zwierz; Bhattacharyya), and the saturating two-level family
> (Hörnedal–Sönnerborn). See [History](#history) for a withdrawn correction notice.

---

## The result

Let $H$ be a time-independent Hamiltonian with ground energy $E_0$, let $A$ be a bounded Hermitian
observable with spectral spread $\sigma_A = (\lambda_{\max}-\lambda_{\min})/2$, and let
$\Delta = |\langle A(T)\rangle - \langle A(0)\rangle|$ be the change of its expectation value under
$e^{-iHt}$ ($\hbar = 1$). Then for every pure **or mixed** state,

$$T\,(\langle H\rangle - E_0)\ \ge\ C_\star\,\frac{\Delta^2}{\sigma_A^2},
\qquad
C_\star = \frac{1}{8\sin x_\star} = 0.172506267461\ldots$$

where $x_\star = 2.331122370\ldots$ is the smallest positive root of $\tan(x/2)=x$.

**What is new here** (see the positioning note above for what is inherited):

- **A no-go theorem — the original kernel.** There is *no* state-independent **linear** mean-energy
  bound on $\Delta$; the optimal mean-energy exponent of $\Delta$ is exactly **two**. The bound above
  is therefore the best possible in its exponent.
- **The eigenvector linear-recovery theorem.** The quadratic law degrades to a linear one when the
  initial state is an eigenvector of $A$ — a *sufficient* condition; a zero initial slope alone is not
  enough, as the manuscript shows by explicit counterexample.

**What is inherited, and cited as such.** $C_\star$ itself is the $q\to0$ Giovannetti–Lloyd–Maccone
constant, the tangent-line lemma is theirs, and the proof route follows Jones–Kok / Zwierz; the
saturating two-level family is the one analyzed by Hörnedal–Sönnerborn. $C_\star$ *is* saturated (as an
infimum) by that family, and the paper gives the exact trade-off curve $P_\star(\delta)$ of which
$C_\star$ is the small-swing slope — the observable analog of the GLM curve for states. That analogy is
the point: this is the observable-target member of an existing family, not a new corner of the table.

The $2\times2$ speed-limit table (what moves: *state* vs *observable*) $\times$ (resource: *variance*
vs *mean energy*) organizes the setting. The (mean-energy $\times$ observable) entry is the **only
quadratic** one, and that quadratic character — not the discovery of an empty cell — is what this work
establishes and explains.

|                        | resource = variance $\Delta H$ (Mandelstam–Tamm) | resource = mean energy $\langle H\rangle-E_0$ (Margolus–Levitin) |
| ---------------------- | ------------------------------------------------ | --------------------------------------------------------------- |
| **state** $1-F$        | $L \le \Delta H\,T$                               | $1-F \le K\,T(\langle H\rangle-E_0)$  (linear)                  |
| **observable** $\Delta\langle A\rangle$ | $T \ge \Delta/(2\Delta H\,\sigma_A)$  (linear) | $T(\langle H\rangle-E_0) \ge C_\star\,\Delta^2/\sigma_A^2$  **(quadratic; this work)** |

The paper also proves the full constant survives for mixed states (joint convexity of the trace
distance), sharpens it for bandwidth-limited generators and for several observables at once, and
shows the quadratic law degrades to a linear one when the initial state is an eigenvector of
$A$ (a sufficient condition). Its most natural physical setting is the **autonomous quantum clock**, where it gives a coherent
mean-energy resolution floor complementary to the known entropy and rate bounds. **That application —
and the quantum-battery one — is illustrative, not binding:** neither currently constrains any existing
platform, because the bound is quadratic and therefore weakest exactly in the operationally interesting
large-swing regime. Demonstrating a single regime where it genuinely binds is the main open task.

---

## Repository layout

```
observable-margolus-levitin/
├── paper/
│   ├── main.tex          # the manuscript (RevTeX 4-2)
│   ├── refs.bib          # bibliography
│   ├── main.bbl          # pre-built bibliography (so it compiles without bibtex)
│   └── figs/             # figures (fig1–fig5, PDF)
├── numerics/
│   ├── README.md         # claim → script map (every quoted number → the script that emits it)
│   ├── results.json      # consolidated machine-readable record of the pinned numbers
│   ├── theory_core.py    # first-pass script: the phase-restricted constant C_old
│   ├── theory_core_v2.py # the constant C*, pinned to 52 digits
│   ├── explore_observable_ML.py # exploratory scan behind the no-go exponent
│   ├── audit_d2.py       # the d = 2 optimizer minimum, re-evaluated four ways
│   ├── make_results.py   # regenerates results.json
│   ├── verify_refinements.py  # trade-off curve, bandwidth constant, two-observable bound
│   ├── optimize_constant.py   # global optimizer floor (d = 2..10)
│   ├── ext_mixed_check.py / ext_mixed_purification_check.py  # mixed-state bound
│   ├── struct_*.py            # structured observables & linear recovery
│   ├── app_battery_check.py   # quantum-battery application
│   └── make_figures.py        # regenerates paper/figs/*.pdf
├── requirements.txt
├── CITATION.cff
└── LICENSE
```

---

## Reproducing the numerics

Every quantitative claim in the paper is pinned to a seeded script. Conventions throughout:
$\hbar = 1$, $E_0 = 0$, $\sigma_A = (\lambda_{\max}-\lambda_{\min})/2$, and the working ratio
$r := T(\langle H\rangle-E_0)\,\sigma_A^2/\Delta^2$ (the bound is $r \ge C_\star$).

```bash
python -m pip install -r requirements.txt

python -u numerics/make_results.py         # regenerates results.json
python -u numerics/verify_refinements.py   # trade-off curve, bandwidth-resolved constant, two-observable bound
python -u numerics/theory_core_v2.py       # the constant C*, to 52 digits
python -u numerics/make_figures.py          # regenerates paper/figs/*.pdf
```

See [`numerics/README.md`](numerics/README.md) for the full **claim → script → reproduced value**
table. Random-search "min ratio $\ge C_\star$" results are stress tests, not proofs; the bound
itself is proved analytically in the manuscript.

Verified fresh on Python 3.12 with numpy 2.5.0 / scipy 1.18.0 / mpmath 1.3.0: `make_results.py` and
`verify_refinements.py` both reproduce the paper's numbers with zero bound violations across ~6 million
sampled configurations.

---

## Building the paper

The latest compiled manuscript is committed at [`paper/main.pdf`](paper/main.pdf). The GitHub Actions
workflow ([`.github/workflows/build-paper.yml`](.github/workflows/build-paper.yml)) recompiles
`paper/main.tex` on every push that touches `paper/` and uploads the result as a build artifact, so a
source change that would break the build cannot pass unnoticed. It does not overwrite the committed
PDF; that is regenerated and committed by hand.

Locally, with a TeX distribution that includes RevTeX 4-2:

```bash
cd paper
pdflatex main
bibtex   main
pdflatex main
pdflatex main
```

---

## Status

This is a **preprint** — a self-contained, independently cross-checked manuscript that has **not** yet
been peer-reviewed or submitted to a journal. Feedback and corrections are welcome via the issue
tracker.

The committed `paper/main.pdf` is compiled from the committed `paper/main.tex` and may be cited as-is.

## History

An August 2026 revision of this README carried a "correction notice" retracting several attribution and
formula claims. **That notice was itself in error and was withdrawn.** Checked against the first commit:

- The manuscript credits Giovannetti–Lloyd–Maccone for the tangent-line constant, in the abstract and
  again at the point of derivation; $C_\star$ was never claimed as new.
- It cites Jones–Kok, Zwierz and Bhattacharyya as the proof route it follows, and names its own new step
  (a bounded observable on the left via Hölder). All three were in `refs.bib` and cited from the start,
  as was Hörnedal–Sönnerborn, which it engages on the time-dependent-$H$ point.
- Sec. V prints $C_{\beta=0} = 1/(16\,\phi_\star\sin^2\phi_\star) = \phi_\star/(4(1-\cos\phi_\star)^2)
  = 0.185551\ldots$, with $\phi_\star = 2.78650\ldots$ defined in the same sentence as the smallest
  positive root of $1-\cos x = 2x\sin x$. The formula the notice called erroneous appears nowhere in the
  manuscript.

The one substantive point was kept: describing the (mean-energy $\times$ observable) corner as "empty"
overstates it, since Hörnedal–Sönnerborn analyse an adjacent near-ground two-level family in the
time-dependent-$H$ setting. The bound, the constant, and every number in the paper are unaffected.

## Citation

If you use this work, please cite the manuscript (see [`CITATION.cff`](CITATION.cff)):

> B. Nasr, *A Margolus–Levitin speed limit for observables: mean energy bounds expectation-value
> change quadratically*, arXiv:2608.22658 [quant-ph] (2026).

The code and data in this repository are archived independently and have their own DOI:
[10.5281/zenodo.22165921](https://doi.org/10.5281/zenodo.22165921) (concept DOI; always resolves to
the latest release).

## License

Code (everything under `numerics/`) is released under the [MIT License](LICENSE); this is the licence
recorded in [`CITATION.cff`](CITATION.cff). The manuscript text and figures under `paper/` are
© 2026 Bryan Nasr, released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Contact

Bryan Nasr — bryannasr4@gmail.com
