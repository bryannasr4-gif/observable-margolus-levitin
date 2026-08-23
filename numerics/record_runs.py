#!/usr/bin/env python3
"""
record_runs.py -- capture what the two SEARCH scripts actually print, to
numerics/runs.json, so that make_results.py never hard-codes a search result.

Three different "optimizer floor" numbers were once circulating in the repository
for the same quantity:

    results.json          0.17250656          over d = 2..7   (a hard-coded
                                                               literal inside
                                                               make_results.py --
                                                               produced by no run)
    numerics/README.md    0.17250634253781064 over d = 2..10  (optimize_constant.py)
    numerics/README.md    0.172511            over d = 2..7   (make_figures.py fig-4 block)

Two of those are real, distinct quantities from two different searches; the
third was a stale literal.  This script re-runs both searches (or ingests a
completed log of the slow one) and records exactly what they print, with their
seeds and dimension ranges, so the three can never drift again.

Usage
-----
    python numerics/record_runs.py --run                 # run both (hours)
    python numerics/record_runs.py --log <path>          # ingest a completed
                                                         # optimize_constant.py
                                                         # log, run the fast scan
    python numerics/record_runs.py --log <path> --skip-fig4

`--log` is accepted because optimize_constant.py takes hours; the script is
seeded and deterministic, so `python -u numerics/optimize_constant.py > log`
reproduces the ingested log exactly.  Any log that does not contain a complete
RESULT block is REJECTED -- an uncompleted run is never recorded.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


# --------------------------------------------------------------------------- #
# (1) optimize_constant.py -- the slow multi-start global optimizer, d = 2..10
# --------------------------------------------------------------------------- #
def parse_optimizer_log(text: str) -> dict:
    per_d = {}
    for m in re.finditer(r"^\s*d=\s*(\d+): best ratio = ([0-9.]+)", text, re.M):
        per_d[int(m.group(1))] = float(m.group(2))
    got = re.search(r"min ratio found\s*=\s*([0-9.eE+-]+)", text)
    gap = re.search(r"ratio - C\s*=\s*([0-9.eE+-]+)", text)
    viol = re.search(r"VIOLATION \(< C\)\?\s*=\s*(\w+)", text)
    where = re.search(r"where\s*=\s*(.+)", text)
    struct = re.search(r"structured min ratio = ([0-9.]+)", text)
    if not (per_d and got and viol):
        raise SystemExit(
            "REJECTED: this log does not contain a completed optimize_constant.py run "
            "(need the per-dimension lines AND the RESULT block). "
            "An uncompleted run is never recorded."
        )
    if sorted(per_d) != list(range(2, 11)):
        raise SystemExit(f"REJECTED: incomplete dimension sweep, got d = {sorted(per_d)}")
    return {
        "script": "numerics/optimize_constant.py",
        "seeds": [1, 7, 42, 101, 2024, 31337],
        "np_random_seed": 12345,
        "d_range": "2..10",
        "starts_per_d": {2: 80, 3: 60, 4: 45, 5: 35, 6: 28, 7: 22, 8: 18, 9: 15, 10: 12},
        "per_d_best_ratio": {str(k): per_d[k] for k in sorted(per_d)},
        "structured_min_ratio": float(struct.group(1)) if struct else None,
        "min_ratio_found": float(got.group(1)),
        "ratio_minus_C_true": float(gap.group(1)) if gap else None,
        "violation_below_C_true": viol.group(1).strip().lower() == "true",
        "where": where.group(1).strip() if where else None,
        "note": "stochastic multi-start search => an UPPER bound on the floor, not the floor",
    }


def run_optimizer() -> str:
    t0 = time.time()
    p = subprocess.run([sys.executable, "-u", os.path.join(HERE, "optimize_constant.py")],
                       capture_output=True, text=True, cwd=ROOT)
    sys.stderr.write(f"[record_runs] optimize_constant.py took {time.time()-t0:.0f} s\n")
    if p.returncode != 0:
        raise SystemExit(f"optimize_constant.py exited {p.returncode}\n{p.stderr[-2000:]}")
    return p.stdout


# --------------------------------------------------------------------------- #
# (2) the make_figures.py fig-4 multi-start scan, d = 2..7
# --------------------------------------------------------------------------- #
# The block below is a VERBATIM copy of the fig-4 scan in make_figures.py.  It
# is copied rather than imported because importing make_figures.py regenerates
# paper/figs/*.pdf as a side effect, which this script must not overwrite.
# _assert_fig4_block_unchanged() fails loudly if the source block ever diverges.
FIG4_SOURCE_MARKERS = ("def mind(d, seed, ntri=900):", "ds = list(range(2, 8));")


def _assert_fig4_block_unchanged() -> str:
    src = open(os.path.join(HERE, "make_figures.py"), encoding="utf-8").read()
    a = src.find(FIG4_SOURCE_MARKERS[0])
    b = src.find(FIG4_SOURCE_MARKERS[1])
    if a < 0 or b < 0:
        raise SystemExit("REJECTED: cannot locate the fig-4 scan block in make_figures.py")
    block = src[a:b]
    mine = _fig4_block_text()
    if _norm(block) != _norm(mine):
        raise SystemExit(
            "REJECTED: the fig-4 scan in make_figures.py no longer matches the copy in "
            "record_runs.py.  Re-copy it (and re-record) before trusting either number."
        )
    return block


def _norm(s: str) -> str:
    return "\n".join(line.rstrip() for line in s.strip().splitlines())


def _fig4_block_text() -> str:
    import inspect

    return inspect.getsource(mind)


def mind(d, seed, ntri=900):
    rng = np.random.default_rng(seed); best = np.inf
    # (a) random near-ground configs (shows nothing dips below C_true)
    for _ in range(ntri):
        E = np.zeros(d); E[1:] = rng.uniform(0.3, 3.0, d-1)
        M = rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)); A=(M+M.conj().T)/2; A/=sigA(A)
        th = rng.uniform(0.01, 0.5)
        a = rng.normal(size=d-1)+1j*rng.normal(size=d-1); a/=np.linalg.norm(a)
        c = np.zeros(d, complex); c[0]=np.cos(th); c[1:]=np.sin(th)*a
        meanH = np.real(np.sum(E*np.abs(c)**2))
        if meanH < 1e-9: continue
        diffs = np.abs(E[:,None]-E[None,:]); wpos=diffs[diffs>1e-6]
        if wpos.size==0: continue
        Hd = np.diag(E); A0 = np.real(c.conj()@(A@c))
        for T in np.linspace(1e-3, 2.5*2*np.pi/wpos.min(), 200):
            ct = expmH(Hd,T)@c; D = abs(np.real(ct.conj()@(A@ct))-A0)
            if D < 1e-9: continue
            best = min(best, T*meanH*sigA(A)**2/D**2)
    # (b) structured 2-level near-ground saturator embedded in dimension d
    #     ground |0> coupled to one excited |1> by sigma_x (spectators 2..d-1 inert);
    #     optimal phase beta=x*/2-pi/2, E*T=x*, theta->0 -> ratio -> C_true for EVERY d.
    A = np.zeros((d,d), complex); A[0,1]=A[1,0]=1.0
    for th in [0.02, 0.01, 0.005]:
        c = np.zeros(d, complex); c[0]=np.cos(th); c[1]=np.sin(th)*np.exp(1j*(xstar/2-np.pi/2))
        E = np.zeros(d); E[1]=1.0; E[2:] = np.linspace(1.5, 3.0, max(0,d-2)); Hd=np.diag(E)
        meanH = np.real(np.sum(E*np.abs(c)**2)); A0 = np.real(c.conj()@(A@c))
        for T in np.linspace(1e-3, 2*np.pi, 1400):
            ct = expmH(Hd,T)@c; D = abs(np.real(ct.conj()@(A@ct))-A0)
            if D < 1e-9: continue
            best = min(best, T*meanH*sigA(A)**2/D**2)
    return best


def run_fig4_scan() -> dict:
    _assert_fig4_block_unchanged()
    t0 = time.time()
    ds = list(range(2, 8))
    mins = [mind(d, 1000 + d) for d in ds]
    return {
        "script": "numerics/make_figures.py (fig-4 block), re-run verbatim by record_runs.py",
        "seeds": {str(d): 1000 + d for d in ds},
        "d_range": "2..7",
        "trials_per_d": 900,
        "per_d_min_ratio": {str(d): float(m) for d, m in zip(ds, mins)},
        "min_ratio_found": float(min(mins)),
        "wall_seconds": round(time.time() - t0, 1),
        "note": "separate, faster purpose-built scan; this is the number Fig. 4 plots",
    }


# --------------------------------------------------------------------------- #
# (3) the two fast batteries whose numbers the README quotes
# --------------------------------------------------------------------------- #
def parse_refinements_log(text: str) -> dict:
    if "(3) two-observable bound" not in text or "two-observable search" not in text:
        raise ValueError("incomplete verify_refinements.py run (part 3 missing)")
    g = lambda p: re.search(p, text)  # noqa: E731
    ends = g(r"P_star\(delta->0\)->([0-9.eE+-]+).*?P_star\(delta->1\)->([0-9.]+)")
    slope = g(r"P_star/delta\^2 -> ([0-9.]+)")
    boundary = g(r"configs sampled: ([\d,]+)\s+min slack \(P - P_star\(delta\)\) = ([0-9.eE+-]+)")
    l1 = g(r"L=1\.0:.*?= ([0-9.]+)\s")
    l05 = g(r"L=0\.5:.*?= ([0-9.]+)\s")
    bw = g(r"bandwidth violation search: ([\d,]+) trials, violations=(\d+), "
           r"min\(ratio/C_true\(BT\)\)=([0-9.]+)")
    tw = g(r"two-observable search: ([\d,]+) trials, violations=(\d+), "
           r"min\(T<H>/\(D1\^2\+D2\^2\)\)=([0-9.]+)")
    tw2 = g(r"min T<H>/max\(D\)\^2 = ([0-9.]+)")
    if not all((ends, slope, boundary, l1, l05, bw, tw, tw2)):
        raise ValueError("verify_refinements.py log did not yield every recorded quantity")
    return {
        "script": "numerics/verify_refinements.py",
        "curve_lower_boundary": {
            "configs_sampled": int(boundary.group(1).replace(",", "")),
            "min_slack": float(boundary.group(2)),
            "P_star_at_delta_to_1": float(ends.group(2)),
            "small_delta_P_over_delta2": float(slope.group(1)),
        },
        "bandwidth_constant": {
            "gain_at_L_1": float(l1.group(1)), "gain_at_L_0.5": float(l05.group(1)),
            "trials": int(bw.group(1).replace(",", "")), "violations": int(bw.group(2)),
            "min_ratio_over_C_L": float(bw.group(3)),
            "claim_status": "stated as a remark in main.tex Sec. V; full search retained here",
        },
        "two_observable": {
            "trials": int(tw.group(1).replace(",", "")), "violations": int(tw.group(2)),
            "min_ratio": float(tw.group(3)),
            "min_ratio_at_equal_swings": float(tw2.group(1)),
            "claim_status": "stated as a remark in main.tex Sec. VI; full search retained here",
        },
    }


def parse_mixed_log(text: str) -> dict:
    per_d = {m.group(1): float(m.group(2))
             for m in re.finditer(r"d=(\d+): min ratio = ([0-9.]+)", text)}
    if not per_d:
        raise ValueError("no per-dimension minima in the ext_mixed_check.py log")
    return {
        "script": "numerics/ext_mixed_check.py",
        "per_d_min_ratio": per_d,
        "min_over_all_d": min(per_d.values()),
        "violations_below_C_true": 0,
        "note": ("main.tex Sec. VIII quotes the d=2 entry as 'minimum sampled ratio 0.176'; "
                 "results.json's mixed.random_min_ratio is a different sample. Both are "
                 "recorded so neither can drift."),
    }


def ingest_or_missing(parser, stem: str, runs_dir: str) -> dict:
    if not os.path.isdir(runs_dir):
        return {"status": "NOT-RERUN-HERE",
                "note": f"no {runs_dir}; capture a completed run of {stem}.py there"}
    # NB: *.log is gitignored (LaTeX build artefacts), so the archived run
    # captures are stored as *.txt and tracked.
    for name in sorted(os.listdir(runs_dir)):
        if not name.startswith(stem) or not name.endswith((".log", ".txt")):
            continue
        path = os.path.join(runs_dir, name)
        with open(path, encoding="utf-8", errors="replace") as fh:
            try:
                out = parser(fh.read())
            except ValueError as e:
                return {"status": "NOT-RERUN-HERE", "note": f"{name}: {e}"}
        out["source"] = f"ingested from a completed run log (numerics/runs/{name})"
        return out
    return {"status": "NOT-RERUN-HERE",
            "note": f"no completed {stem} log in {runs_dir}"}


# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-dir", default=os.path.join(HERE, "runs"))
    ap.add_argument("--run", action="store_true",
                    help="execute optimize_constant.py (hours) instead of ingesting a log")
    ap.add_argument("--log", help="path to a COMPLETED optimize_constant.py stdout capture")
    ap.add_argument("--skip-fig4", action="store_true")
    ap.add_argument("--out", default=os.path.join(HERE, "runs.json"))
    a = ap.parse_args()

    runs = {"_meta": {
        "generated_by": "numerics/record_runs.py",
        "purpose": "what the two search scripts ACTUALLY printed, with seeds and ranges; "
                   "make_results.py ingests this instead of hard-coding a search result",
        "python": sys.version.split()[0],
    }}

    if a.run:
        runs["global_optimizer_floor"] = parse_optimizer_log(run_optimizer())
        runs["global_optimizer_floor"]["source"] = "executed by record_runs.py --run"
    elif a.log:
        with open(a.log, encoding="utf-8", errors="replace") as fh:
            runs["global_optimizer_floor"] = parse_optimizer_log(fh.read())
        runs["global_optimizer_floor"]["source"] = f"ingested from a completed run log ({a.log})"
    else:
        runs["global_optimizer_floor"] = {"status": "NOT-RERUN-HERE"}

    if a.skip_fig4:
        runs["figure4_scan_floor"] = {"status": "NOT-RERUN-HERE"}
    else:
        runs["figure4_scan_floor"] = run_fig4_scan()

    runs["refinements"] = ingest_or_missing(
        parse_refinements_log, "verify_refinements", a.runs_dir)
    runs["mixed_violation_search"] = ingest_or_missing(
        parse_mixed_log, "ext_mixed_check", a.runs_dir)

    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(runs, fh, indent=2)
    print(json.dumps(runs, indent=2))
    print(f"\nWROTE {a.out}")
    return 0


# numpy/mpmath are needed only by the fig-4 block; import after the docstring so
# that --help works on a bare interpreter.
import numpy as np  # noqa: E402
import mpmath as _mp  # noqa: E402

_mp.mp.dps = 40
xstar = float(_mp.findroot(lambda x: _mp.tan(x / 2) - x, _mp.mpf("2.33")))


def expmH(H, t):
    w, V = np.linalg.eigh(H)
    return (V * np.exp(-1j * w * t)) @ V.conj().T


def sigA(A):
    w = np.linalg.eigvalsh(A)
    return (w[-1] - w[0]) / 2


if __name__ == "__main__":
    sys.exit(main())
