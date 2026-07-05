# A Margolus–Levitin speed limit for observables

**Mean energy bounds the change of an expectation value — quadratically.**

[![build-paper](https://github.com/bryannasr4-gif/observable-margolus-levitin/actions/workflows/build-paper.yml/badge.svg)](https://github.com/bryannasr4-gif/observable-margolus-levitin/actions/workflows/build-paper.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Paper: RevTeX](https://img.shields.io/badge/paper-RevTeX%204--2-informational.svg)](paper/main.tex)
[![Status: preprint](https://img.shields.io/badge/status-preprint-orange.svg)](#status)

This repository contains the manuscript and the complete, seeded numerical verification suite for a
quantum speed limit on **observables** powered by the **mean energy** — the missing member of the
Mandelstam–Tamm / Margolus–Levitin family.

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

Two facts make this new:

- **A no-go theorem.** There is *no* state-independent **linear** mean-energy bound on $\Delta$; the
  optimal mean-energy exponent of $\Delta$ is exactly **two**. The bound above is the best possible.
- **It is tight.** $C_\star$ is saturated (as an infimum) by a near-ground two-level family, and the
  paper gives the full **exact trade-off curve** $P_\star(\delta)$ of which $C_\star$ is the
  small-swing slope — the observable analog of the Giovannetti–Lloyd–Maccone curve for states.

This completes the $2\times2$ speed-limit table (what moves: *state* vs *observable*) $\times$
(resource: *variance* vs *mean energy*). Three corners were known; the (mean-energy $\times$
observable) corner — the **only quadratic** one — was empty.

|                        | resource = variance $\Delta H$ (Mandelstam–Tamm) | resource = mean energy $\langle H\rangle-E_0$ (Margolus–Levitin) |
| ---------------------- | ------------------------------------------------ | --------------------------------------------------------------- |
| **state** $1-F$        | $L \le \Delta H\,T$                               | $1-F \le K\,T(\langle H\rangle-E_0)$  (linear)                  |
| **observable** $\Delta\langle A\rangle$ | $T \ge \Delta/(2\Delta H\,\sigma_A)$  (linear) | $T(\langle H\rangle-E_0) \ge C_\star\,\Delta^2/\sigma_A^2$  **(quadratic; this work)** |

The paper also proves the full constant survives for mixed states (joint convexity of the trace
distance), sharpens it for bandwidth-limited generators and for several observables at once, and
shows the quadratic law degrades to a linear one exactly when the initial state is an eigenvector of
$A$. Its cleanest physical home is the **autonomous quantum clock**, where it gives a coherent
mean-energy resolution floor complementary to the known entropy and rate bounds.

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
│   ├── theory_core_v2.py # the constant C*, pinned to 52 digits
│   ├── make_results.py   # regenerates results.json
│   ├── verify_refinements.py  # trade-off curve, bandwidth constant, two-observable bound
│   ├── optimize_constant.py   # global optimizer floor (d = 2..7)
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
workflow ([`.github/workflows/build-paper.yml`](.github/workflows/build-paper.yml)) also recompiles
`paper/main.tex` on every push and uploads the PDF as a build artifact, so the source and the PDF stay
in sync.

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

## Citation

If you use this work, please cite the manuscript (see [`CITATION.cff`](CITATION.cff)):

> B. Nasr, *A Margolus–Levitin speed limit for observables: mean energy bounds expectation-value
> change quadratically* (2026).

## License

Code (everything under `numerics/`) is released under the [MIT License](LICENSE). The manuscript text
and figures under `paper/` are © 2026 Bryan Nasr, released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Contact

Bryan Nasr — bryannasr4@gmail.com
