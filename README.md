# A Margolus–Levitin speed limit for observables

**Mean energy bounds the change of an expectation value — quadratically.**

[![build-paper](https://github.com/bryannasr4-gif/observable-margolus-levitin/actions/workflows/build-paper.yml/badge.svg)](https://github.com/bryannasr4-gif/observable-margolus-levitin/actions/workflows/build-paper.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Paper: RevTeX](https://img.shields.io/badge/paper-RevTeX%204--2-informational.svg)](paper/main.tex)
[![Status: preprint](https://img.shields.io/badge/status-preprint-orange.svg)](#status)

This repository contains the manuscript and the complete, seeded numerical verification suite for a
quantum speed limit on **observables** powered by the **mean energy**.

> ### Correction notice (August 2026) — novelty re-scoped, one formula corrected
>
> An external prior-art and correctness audit found that **this repository previously over-scoped its
> novelty**. The mathematics is sound, but two of the claims made for it were not. Both are corrected
> here, and the manuscript is being revised to match. Nothing about the *validity* of the bound has
> changed.
>
> **1. The constant $C_\star$ is not new.** $x_\star = 2.33112237\ldots$ is *bit-identical* to the
> $q\to 0$ (orthogonality) member of the tangent-line optimization of **Giovannetti, Lloyd & Maccone,
> PRA 67, 052109 (2003)** — it is the root of $1-\cos y = y\sin y$, and $K=\sin x_\star = 0.724611\ldots$
> is the standard tangent-to-cosine constant of the generalized-ML literature (used off-the-shelf as
> $\cos x \ge 1-\alpha x$, $\alpha \approx 0.724$, in Carabba–Hörnedal–del Campo, Quantum 6, 884 (2022)).
> Quoting it to 52 digits as a fresh discovery was wrong. *(Verify:
> `mpmath.findroot(lambda y: 1-cos(y)-y*sin(y), 2.33)` → $x_\star$.)*
>
> **2. The proof route is inherited.** The two-step derivation and the purification-to-mixed argument
> follow the statistical-distance route to Margolus–Levitin of **Jones & Kok, PRA 82, 022107 (2010)**,
> corrected by **Zwierz, PRA 86, 016101 (2012)**; the underlying inequality traces to **Bhattacharyya,
> J. Phys. A 16, 2993 (1983)**. Placing a bounded observable on the left via Hölder — converting a state
> distinguishability measure into an expectation-value swing — is a genuine repackaging, but it is a
> repackaging, and these references were absent from the bibliography.
>
> **3. Adjacent geometry already in the literature.** **Hörnedal & Sönnerborn, PRA 108, 052421 (2023)**
> analyze the same near-ground two-level family and prove that no ML-type bound survives for *closed*
> (time-dependent-$H$) systems. That is a different setting and does **not** contradict the result here
> — it is why the time-independent-$H$ restriction is essential — but it must be engaged rather than
> omitted.
>
> **4. What is actually new.** The genuinely original kernel is the **no-go theorem** (the optimal
> mean-energy exponent of $\Delta$ is exactly **two**, so no state-independent *linear* mean-energy
> bound exists) together with the **eigenvector dichotomy** (the quadratic law degrades to linear
> exactly when the initial state is an eigenvector of $A$). That is a real contribution. It is
> *solid-incremental*, not a landmark, and the earlier framing of an "empty corner" oversold it.
>
> **5. Erratum, Sec. V.** The phase-restricted constant was printed as
> $C_{\beta=0} = s_\star^2/(8\sqrt{2s_\star-1}) = 0.185551$. **The value 0.185551 is correct; the
> formula is not** — it evaluates to $0.453865\ldots$, and $s_\star$ was never defined. The correct
> closed form is
> $$C_{\beta=0} \;=\; \frac{1}{16\,s_\star \sin^2 s_\star} \;=\; \frac{s_\star}{4\,(1-\cos s_\star)^2} \;=\; 0.18555147\ldots,$$
> where $s_\star = 2.78649815\ldots$ is the root of $1-\cos x = 2x\sin x$ — a *different* transcendental
> equation from $x_\star$'s. Verified three ways (both closed forms, and direct minimization of
> $x/[4(1-\cos x)^2]$, which attains $0.18555147$ at $x = 2.786496$).
>
> **6. The applications do not currently bind.** The autonomous-clock and quantum-battery sections are
> illustrative, not binding constraints on any existing platform. They should be read as such.
>
> The audit that produced this notice is an automated multi-agent prior-art review; every item above was
> independently re-verified before being published here, and its numerical claims are reproducible with
> the one-line recipes given.

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

**What is new here** (see the correction notice above for what is not):

- **A no-go theorem — the original kernel.** There is *no* state-independent **linear** mean-energy
  bound on $\Delta$; the optimal mean-energy exponent of $\Delta$ is exactly **two**. The bound above
  is therefore the best possible in its exponent.
- **The eigenvector dichotomy.** The quadratic law degrades to a linear one *exactly* when the initial
  state is an eigenvector of $A$ — a sharp characterization of when the quadratic penalty disappears.

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
shows the quadratic law degrades to a linear one exactly when the initial state is an eigenvector of
$A$. Its most natural physical setting is the **autonomous quantum clock**, where it gives a coherent
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

**Revision in progress (August 2026).** The README above has been corrected; `paper/main.tex` has
**not yet** been updated to match, so the manuscript still carries the over-scoped novelty framing, the
missing citations (GLM 2003, Jones–Kok 2010, Zwierz 2012, Bhattacharyya 1983, Hörnedal–Sönnerborn 2023)
and the Sec. V formula erratum described in the correction notice. **Where the manuscript and this
README disagree, the README is current.** Do not cite the Sec. V closed form from the PDF.

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
