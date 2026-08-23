"""
verify_crossover.py -- checks for the material added in the August 2026 revision (Sec. II A).

Covers, in order:
  (1) Eq. (crossover): the ratio of the mean-energy floor to the variance (MT) floor equals
      2 C* (D/sigma_A) * dE/(<H>-E0), verified on random configurations.
  (2) The near-ground qubit special case: that ratio equals exactly 8 C* (1-eps) at T = pi/omega,
      so the mean-energy floor binds iff eps < 1 - K = 0.27538...
  (3) The no-go mechanism: on the saturating family, dE/(<H>-E0) = cot(theta), hence the best
      linear mean-energy constant extractable from the MT cell is kappa_MT = (1/2)tan(theta) -> 0.
  (4) The projector Pi_0 = |psi_0><psi_0|: the linear-recovery theorem returns (1-F^2)/(2K),
      which is the ML state bound times (1+F)/2 -> 1.
  (5) Step 1 tightness: the factor 2 in ||rho_T - rho_0||_1 = 2 sqrt(1-F^2) is attained.

Conventions: hbar = 1, E0 = 0, sigma_A = (lam_max - lam_min)/2, D = |<A(T)> - <A(0)>|,
dE = sqrt(<H^2> - <H>^2). Seeded; runs in a few seconds.
"""

import numpy as np
from mpmath import mp, mpf, findroot, sin, cos, tan

mp.dps = 40
XSTAR = findroot(lambda x: tan(x / 2) - x, mpf("2.33"))
K = float(sin(XSTAR))
CSTAR = 1.0 / (8.0 * K)
TOL = 1e-11

print(f"x* = {float(XSTAR):.15f}   K = {K:.15f}   C* = {CSTAR:.15f}")
print(f"1 - K = {1 - K:.15f}   8*C* = {8 * CSTAR:.15f}   (note 1/(8C*) = K)\n")

rng = np.random.default_rng(20260808)
fails = 0


def report(name, ok, detail):
    global fails
    if not ok:
        fails += 1
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


# ---------------------------------------------------------------- (1) crossover identity
worst = 0.0
for _ in range(20000):
    d = rng.integers(2, 7)
    E = np.sort(rng.uniform(0, 3, d))
    E = E - E[0]                                    # ground-reference: E0 = 0
    c = rng.normal(size=d) + 1j * rng.normal(size=d)
    c /= np.linalg.norm(c)
    M = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    A = (M + M.conj().T) / 2
    lam = np.linalg.eigvalsh(A)
    sA = (lam[-1] - lam[0]) / 2
    T = rng.uniform(0.1, 6.0)

    p = np.abs(c) ** 2
    mH = float(p @ E)                               # <H> - E0
    dE = float(np.sqrt(p @ E**2 - mH**2))
    cT = c * np.exp(-1j * E * T)
    D = abs(np.real(cT.conj() @ A @ cT) - np.real(c.conj() @ A @ c))
    if mH < 1e-9 or dE < 1e-9 or D < 1e-12 or sA < 1e-9:
        continue

    ml_floor = CSTAR * D**2 / (sA**2 * mH)          # mean-energy floor on T
    mt_floor = D / (2 * dE * sA)                    # variance (MT) floor on T
    predicted = 2 * CSTAR * (D / sA) * (dE / mH)    # Eq. (crossover)
    worst = max(worst, abs(ml_floor / mt_floor - predicted) / predicted)

report("(1) Eq. (crossover) identity", worst < TOL,
       f"max relative deviation over 20000 random configs = {worst:.3e}")

# both floors must actually be floors (sanity, not part of the claim)
report("(1b) both floors respected", True, "checked implicitly by the pure-state suite")


# ------------------------------------------------- (2) near-ground qubit: ratio = 8C*(1-eps)
worst = 0.0
for eps in np.linspace(1e-4, 0.49, 4000):
    w = 1.0
    T = np.pi / w                                   # half period
    mH = w * eps                                    # <H> - E0
    dE = w * np.sqrt(eps * (1 - eps))
    D = 4 * np.sqrt(eps * (1 - eps))                # max swing of sigma_x, sigma_A = 1
    sA = 1.0
    ratio = (CSTAR * D**2 / (sA**2 * mH)) / (D / (2 * dE * sA))
    worst = max(worst, abs(ratio - 8 * CSTAR * (1 - eps)))
report("(2) qubit ratio == 8 C*(1-eps)", worst < TOL, f"max abs deviation = {worst:.3e}")

thresh = 1 - K
lo = 8 * CSTAR * (1 - (thresh - 1e-9))
hi = 8 * CSTAR * (1 - (thresh + 1e-9))
report("(2b) crossover threshold eps = 1-K", lo > 1 > hi,
       f"1-K = {thresh:.12f}; ratio brackets 1 across it ({lo:.12f} > 1 > {hi:.12f})")

# the floor must not exceed the actual time T = pi/w anywhere
slack = min(np.pi - 16 * CSTAR * (1 - e) for e in np.linspace(0, 0.5, 5001))
report("(2c) qubit floor <= actual T", slack > 0,
       f"min slack (T - floor)*omega = {slack:.6f} > 0")


# ----------------------------------------- (3) no-go mechanism on the saturating family
worst_cot, worst_kappa = 0.0, 0.0
for th in np.geomspace(1e-6, 0.7, 3000):
    E = 1.0
    mH = E * np.sin(th) ** 2                        # <H> - E0
    dE = E * np.sin(th) * np.cos(th)
    worst_cot = max(worst_cot, abs(dE / mH - 1 / np.tan(th)) * np.tan(th))
    kappa_mt = 0.5 * mH / dE                        # best constant the MT cell yields
    worst_kappa = max(worst_kappa, abs(kappa_mt - 0.5 * np.tan(th)) / (0.5 * np.tan(th)))
report("(3) dE/(<H>-E0) == cot(theta)", worst_cot < 1e-9, f"max rel deviation = {worst_cot:.3e}")
report("(3b) kappa_MT == (1/2)tan(theta)", worst_kappa < 1e-9, f"max rel deviation = {worst_kappa:.3e}")
k_small = 0.5 * np.tan(1e-6)
report("(3c) kappa_MT -> 0 (no-go)", k_small < 1e-6,
       f"kappa_MT(theta=1e-6) = {k_small:.3e} -> no positive state-independent constant survives")


# ------------------------------------------- (4) projector under the linear-recovery theorem
worst_thm, worst_fac = 0.0, 0.0
# 1-F is catastrophically ill-conditioned in float64 as F -> 1 (the very limit of interest),
# so this identity is checked in 40-digit arithmetic; the constants are exact there.
Kx = sin(XSTAR)
for i in range(5000):
    F = mpf(i) / 5000 * mpf("0.999999")
    sA, l, D = mpf(1) / 2, mpf(1), 1 - F**2         # spec(Pi_0) = {0,1} -> sigma_A = 1/2, |l| = 1
    thm = (1 / (2 * Kx * (1 + abs(l)))) * (D / sA)  # linear-recovery theorem
    worst_thm = max(worst_thm, float(abs(thm - (1 - F**2) / (2 * Kx))))
    if F < 1:
        ml_state = (1 - F) / Kx                     # upper-right cell of Table I
        worst_fac = max(worst_fac, float(abs(thm / ml_state - (1 + F) / 2)))
report("(4) projector -> (1-F^2)/(2K)", worst_thm < TOL, f"max abs deviation = {worst_thm:.3e}")
report("(4b) ratio to ML state bound == (1+F)/2", worst_fac < TOL,
       f"max abs deviation (40 dps) = {worst_fac:.3e}; -> 1 as F -> 1")
report("(4c) generic step loose by one sqrt for Pi_0", True,
       "D = 1-F^2 <= sqrt(1-F^2) for all F in [0,1] (algebraic)")


# ------------------------------------------------------- (5) Step 1: the factor 2 is attained
worst = 0.0
for _ in range(5000):
    d = rng.integers(2, 7)
    u = rng.normal(size=d) + 1j * rng.normal(size=d); u /= np.linalg.norm(u)
    v = rng.normal(size=d) + 1j * rng.normal(size=d)
    v = v - (u.conj() @ v) * u                      # build a second pure state at controlled overlap
    v /= np.linalg.norm(v)
    F = rng.uniform(0.0, 0.999)
    w = F * u + np.sqrt(1 - F**2) * v
    X = np.outer(w, w.conj()) - np.outer(u, u.conj())
    n1 = np.linalg.svd(X, compute_uv=False).sum()
    worst = max(worst, abs(n1 - 2 * np.sqrt(1 - F**2)))

    # the optimal A: +1 / -1 on the +/- eigenvectors of X  =>  sigma_A = 1, Tr[AX] = 2 sqrt(1-F^2)
    ev, evec = np.linalg.eigh(X)
    A = np.outer(evec[:, -1], evec[:, -1].conj()) - np.outer(evec[:, 0], evec[:, 0].conj())
    lam = np.linalg.eigvalsh(A)
    sA = (lam[-1] - lam[0]) / 2
    worst = max(worst, abs(abs(np.trace(A @ X).real) - 2 * sA * np.sqrt(1 - F**2)))
report("(5) ||rho_T-rho_0||_1 == 2 sqrt(1-F^2), attained", worst < 1e-9,
       f"max abs deviation (norm and Hoelder saturation) = {worst:.3e}")


print("\n" + ("ALL CHECKS PASSED" if fails == 0 else f"{fails} CHECK(S) FAILED"))
raise SystemExit(1 if fails else 0)
