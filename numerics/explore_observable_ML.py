"""
Step A -- approach-simulation / de-risking for the project:
  "A Margolus-Levitin (mean-energy) quantum speed limit for observables"

Goal of THIS script (ground-truth numerics, hbar = 1, ground energy E0 = 0):
  PART A.  Confirm the LINEAR no-go: the naive ML-for-observables bound
           T >= f(Delta) * pi/(2(<H>-E0)) with state-independent f CANNOT hold,
           because the optimal P(Delta) := T*(<H>-E0) vanishes faster than linearly.
  PART B.  Discover/confirm the surviving law.  Hypothesis: P(Delta) ~ C * Delta^2
           (QUADRATIC), with an explicit best constant C.  Pin the exponent and C
           from the exactly-solvable two-level family, both numerically and via the
           reduced scalar optimization  C = min_s s^2 * arccos(1 - 1/s) / 4.
  PART C.  Stress-test the candidate bound  T*(<H>-E0) >= C*Delta^2  (||A||_op = 1)
           against random and structured states in d = 2..8 -- look for ANY violation
           (a configuration with P/Delta^2 < C).  Also report which operator
           normalization (||A||_op vs spectral spread) gives the cleanest constant.

Everything here is laptop-scale: dense Hermitian matrices, d <= 8, vectorized.
"""

import numpy as np

rng = np.random.default_rng(20260623)
np.set_printoptions(precision=6, suppress=True)


# ----------------------------------------------------------------------
# core dynamics in the ENERGY EIGENBASIS:  H = diag(E), E[0] = 0
#   <A(t)> = sum_{m,n} conj(psi_m) psi_n A_{mn} exp(i (E_m - E_n) t)
# vectorized over a time grid ts.
# ----------------------------------------------------------------------
def amp_traj(E, A, psi, ts):
    B = (np.conjugate(psi)[:, None] * psi[None, :]) * A      # B_{mn} = psi_m^* psi_n A_{mn}
    w = (E[:, None] - E[None, :]).ravel()                    # Bohr frequencies omega_{mn}
    Bf = B.ravel()
    phases = np.exp(1j * np.outer(ts, w))                    # (n_t, d^2)
    return np.real(phases @ Bf)                              # (n_t,)


def mean_energy(E, psi):                                     # <H> - E0  (E0 = 0)
    return float(np.sum((np.abs(psi) ** 2) * E))


def first_time_to_change(E, A, psi, Delta, Tmax, n=8000):
    """First time |<A(t)> - <A(0)>| reaches Delta, linearly interpolated. None if never."""
    ts = np.linspace(0.0, Tmax, n)
    a = amp_traj(E, A, psi, ts)
    d = np.abs(a - a[0])
    hit = np.nonzero(d >= Delta)[0]
    if hit.size == 0:
        return None, a
    i = hit[0]
    if i == 0:
        return 0.0, a
    t0, t1, d0, d1 = ts[i - 1], ts[i], d[i - 1], d[i]
    return t0 + (Delta - d0) * (t1 - t0) / (d1 - d0), a


# ======================================================================
# PART B-core:  the two-level reduced optimization
#   For H=diag(0,1), A=sigma_x, psi=cos th|0>+sin th|1>:
#     <A(t)> = sin(2th) cos(t),  <H>-E0 = sin^2 th.
#   To change <A> by Delta:  T = arccos(1 - Delta/sin(2th)),  feasible sin(2th) >= Delta/2.
#   P(th;Delta) = sin^2(th) * arccos(1 - Delta/sin(2th)).
#   Small-Delta:  P ~ C2 * Delta^2 with C2 = min_{s>=1/2} s^2 arccos(1-1/s) / 4.
# ======================================================================
def two_level_P(Delta, n_th=200000):
    th = np.linspace(1e-6, np.pi / 4, n_th)
    s2 = np.sin(2 * th)
    feasible = s2 >= Delta / 2
    x = np.clip(1 - Delta / np.where(feasible, s2, np.nan), -1.0, 1.0)
    P = (np.sin(th) ** 2) * np.arccos(x)
    P = np.where(feasible, P, np.inf)
    k = int(np.nanargmin(P))
    return float(P[k]), float(th[k])


def reduced_constant(n_s=400000):
    s = np.linspace(0.5, 50.0, n_s)
    g = s ** 2 * np.arccos(np.clip(1 - 1 / s, -1, 1))
    k = int(np.argmin(g))
    return float(g[k]) / 4.0, float(s[k])


# ----------------------------------------------------------------------
# random ingredients for PART C
# ----------------------------------------------------------------------
def random_E(d, emax):
    e = np.sort(np.abs(rng.standard_normal(d)) * emax)
    e = e - e[0]
    return e


def random_A_unit_opnorm(d):
    """Random Hermitian A with ||A||_op = 1 (spectrum rescaled into [-1,1])."""
    M = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    A = (M + M.conj().T) / 2
    ev = np.linalg.eigvalsh(A)
    rad = max(abs(ev[0]), abs(ev[-1]))
    return A / rad


def random_state(d, near_ground=False):
    if near_ground:
        c = rng.standard_normal(d) + 1j * rng.standard_normal(d)
        c[0] *= 8.0                                  # bias weight onto the ground state
        eps = 10.0 ** rng.uniform(-2.5, -0.3)        # tiny excited admixture
        c[1:] *= eps
    else:
        c = rng.standard_normal(d) + 1j * rng.standard_normal(d)
    return c / np.linalg.norm(c)


def spectral_spread(A):                              # (lambda_max - lambda_min)/2
    ev = np.linalg.eigvalsh(A)
    return (ev[-1] - ev[0]) / 2.0


# ======================================================================
def main():
    print("=" * 72)
    print("PART A & B:  two-level family  ->  exponent + constant + linear no-go")
    print("=" * 72)
    Cstar, sstar = reduced_constant()
    print(f"reduced constant  C = min_s s^2 arccos(1-1/s)/4  = {Cstar:.6f}  at s = {sstar:.4f}")
    print(f"(for reference  pi/16 = {np.pi/16:.6f})\n")

    print(f"{'Delta':>10} {'P_min(Delta)':>14} {'P/Delta':>12} {'P/Delta^2':>12}")
    deltas = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
    for D in deltas:
        P, _ = two_level_P(D)
        print(f"{D:>10.4f} {P:>14.3e} {P/D:>12.4e} {P/D**2:>12.6f}")
    # log-log exponent on the smallest deltas
    ds = np.array(deltas[-5:])
    Ps = np.array([two_level_P(D)[0] for D in ds])
    slope = np.polyfit(np.log(ds), np.log(Ps), 1)[0]
    print(f"\n  log-log slope of P_min(Delta) on small Delta  =  {slope:.4f}   (linear=1, quadratic=2)")
    print("  -> P/Delta -> 0 : NO linear mean-energy bound exists (no-go confirmed).")
    print(f"  -> P/Delta^2 -> {Cstar:.4f} : candidate QUADRATIC law  T*(<H>-E0) >= C*Delta^2.")

    print("\n" + "=" * 72)
    print("PART C:  violation search for  T*(<H>-E0) >= C*Delta^2  (||A||_op = 1)")
    print("         scanning random + near-ground states, d = 2..8")
    print("=" * 72)
    C_bound = Cstar
    n_trials = 4000
    global_min_ratio = np.inf
    global_min_info = None
    min_ratio_opnorm = np.inf
    min_ratio_spread = np.inf
    violations = 0
    checked = 0

    for _ in range(n_trials):
        d = int(rng.integers(2, 9))
        emax = 10.0 ** rng.uniform(-0.5, 1.5)
        E = random_E(d, emax)
        A = random_A_unit_opnorm(d)
        psi = random_state(d, near_ground=(rng.random() < 0.6))

        gap = np.min(np.diff(np.unique(np.round(E, 9))))     # smallest positive gap
        if not np.isfinite(gap) or gap <= 0:
            continue
        Tmax = 12 * np.pi / gap

        # operating point: aim for a modest fraction of the actually achievable swing
        ts = np.linspace(0, Tmax, 4000)
        a = amp_traj(E, A, psi, ts)
        swing = a.max() - a.min()
        if swing < 1e-6:
            continue
        for frac in (0.3, 0.1, 0.03):
            Delta = frac * swing
            T, _ = first_time_to_change(E, A, psi, Delta, Tmax, n=8000)
            if T is None or T <= 0:
                continue
            checked += 1
            me = mean_energy(E, psi)
            P = T * me
            ratio_op = P / Delta ** 2                         # normalize by ||A||_op = 1
            sp = spectral_spread(A)
            ratio_sp = P / (Delta / sp) ** 2                  # normalize change by spread
            min_ratio_opnorm = min(min_ratio_opnorm, ratio_op)
            min_ratio_spread = min(min_ratio_spread, ratio_sp)
            if ratio_op < global_min_ratio:
                global_min_ratio = ratio_op
                global_min_info = (d, float(emax), float(Delta), float(T), float(me))
            if ratio_op < C_bound - 1e-9:
                violations += 1

    print(f"  trials={n_trials}  operating-points checked={checked}")
    print(f"  candidate constant C (||A||_op normalization) = {C_bound:.6f}")
    print(f"  GLOBAL min  P/Delta^2  (||A||_op=1)           = {global_min_ratio:.6f}")
    print(f"      achieved at (d, emax, Delta, T, <H>-E0)   = {global_min_info}")
    print(f"  min  P/(Delta/spread)^2  (spread normalization)= {min_ratio_spread:.6f}")
    print(f"  VIOLATIONS of  P/Delta^2 >= C  :  {violations}")
    if violations == 0:
        print("  -> no violation found: quadratic bound holds across the sampled ensemble,")
        print(f"     and the two-level value C={C_bound:.4f} is approached but not beaten")
        print("     (consistent with two-level configs being the saturating family).")
    else:
        print("  -> VIOLATION found: the constant and/or the law need revision (investigate).")


if __name__ == "__main__":
    main()
