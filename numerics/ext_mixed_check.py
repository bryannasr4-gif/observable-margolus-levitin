"""
P2.2 — Mixed-state constant: numerical stress test.

CLAIM to test: for ANY density matrix rho on finite dim, time-independent H
(E0=0), bounded Hermitian A with spectral spread sigma_A=(lmax-lmin)/2,
    T * <H> * sigma_A^2 / Delta^2  >=  C_true = 1/(8 sin x*) = 0.172506267...
where Delta = |Tr[A (rho_T - rho_0)]|, rho_T = U rho_0 U^dag, U=exp(-i H T).

This SHARPENS the previously-proven mixed bound (C_true/2, via Uhlmann fidelity)
to the full pure-state constant.  Proof route under test (joint convexity):
    Delta <= 2 sigma_A * D(rho_0,rho_T)                       (Holder)
    D(rho_0,rho_T) <= sum_k p_k sqrt(1-F_k^2)                 (joint convexity, p_k=eigs of rho_0)
    sqrt(1-F_k^2) <= sqrt(2(1-F_k)) <= sqrt(2 K T <k|H|k>)    (per-eigenstate ML)
    sum_k p_k sqrt(2KT<k|H|k>) <= sqrt(2KT<H>)               (Jensen, concavity of sqrt)
  => Delta^2 <= 8 K sigma_A^2 T <H>  => ratio >= 1/(8K) = C_true.

We (1) search hard for a VIOLATION of ratio >= C_true over random mixed states,
and (2) directly verify each link of the convexity chain on random instances.
hbar=1.  Fixed seeds for reproducibility.
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 40
xstar = mp.findroot(lambda x: mp.tan(x/2) - x, mp.mpf('2.33'))
K = float((1 - mp.cos(xstar)) / xstar)
C_true = 1.0 / (8 * K)
print(f"C_true = {C_true:.15f}   K = {K:.15f}\n")


def expm_herm(H, t):
    w, V = np.linalg.eigh(H)
    return (V * np.exp(-1j * w * t)) @ V.conj().T


def sigma_A(A):
    w = np.linalg.eigvalsh(A)
    return (w[-1] - w[0]) / 2.0


def trace_distance(rho, sigma):
    return 0.5 * np.sum(np.abs(np.linalg.eigvalsh(rho - sigma)))


def random_state(d, rng, pure=False):
    if pure:
        c = rng.normal(size=d) + 1j * rng.normal(size=d)
        c /= np.linalg.norm(c)
        return np.outer(c, c.conj())
    # random density matrix: random eigenvalues (Dirichlet) in random basis
    p = rng.dirichlet(np.ones(d) * rng.uniform(0.3, 2.0))
    Z = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    Q, _ = np.linalg.qr(Z)
    return (Q * p) @ Q.conj().T


def random_H(d, rng):
    E = np.abs(rng.normal(scale=rng.uniform(0.5, 3.0), size=d))
    E[0] = 0.0
    E -= E.min()
    return np.diag(E), E


def random_A(d, rng):
    M = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    A = (M + M.conj().T) / 2.0
    s = sigma_A(A)
    return A / s  # normalize spectral spread to 1


# ----------------------------------------------------------------------
# (1) Violation search: min ratio over random mixed configs and time grid
# ----------------------------------------------------------------------
def min_ratio_config(H, A, rho, Tgrid):
    Evec = np.linalg.eigvalsh(H)
    meanH = np.real(np.trace(rho @ H))
    if meanH < 1e-12:
        return np.inf
    sA = sigma_A(A)
    A0 = np.real(np.trace(A @ rho))
    best = np.inf
    for T in Tgrid:
        U = expm_herm(H, T)
        rhoT = U @ rho @ U.conj().T
        Delta = abs(np.real(np.trace(A @ rhoT)) - A0)
        if Delta < 1e-9:
            continue
        ratio = T * meanH * sA**2 / Delta**2
        if ratio < best:
            best = ratio
    return best


print("=" * 70)
print("(1) VIOLATION SEARCH over random mixed states (looking for ratio < C_true)")
print("=" * 70)
global_min = np.inf
viol = 0
n_per_d = 4000
for d in range(2, 6):
    rng = np.random.default_rng(7000 + d)
    dmin = np.inf
    for _ in range(n_per_d):
        H, E = random_H(d, rng)
        A = random_A(d, rng)
        rho = random_state(d, rng, pure=(rng.random() < 0.25))
        # time grid: cover slowest Bohr oscillation a couple times
        diffs = np.abs(E[:, None] - E[None, :])
        wpos = diffs[diffs > 1e-6]
        if wpos.size == 0:
            continue
        Tmax = 2.2 * 2 * np.pi / wpos.min()
        Tgrid = np.linspace(1e-4, Tmax, 400)
        r = min_ratio_config(H, A, rho, Tgrid)
        if r < dmin:
            dmin = r
        if r < C_true - 1e-7:
            viol += 1
    print(f"  d={d}: min ratio = {dmin:.10f}   (>= C_true? {dmin >= C_true - 1e-7})")
    global_min = min(global_min, dmin)
print(f"\n  GLOBAL min ratio (random mixed) = {global_min:.10f}")
print(f"  C_true                          = {C_true:.10f}")
print(f"  violations (< C_true - 1e-7)    = {viol}")
print(f"  => {'CONSISTENT with full C_true' if viol == 0 else 'VIOLATION FOUND'}\n")

# ----------------------------------------------------------------------
# (2) Direct check of the convexity chain on random instances
# ----------------------------------------------------------------------
print("=" * 70)
print("(2) DIRECT CHECK of the convexity proof chain (random instances)")
print("=" * 70)
rng = np.random.default_rng(99)
worst_slack = np.inf
for trial in range(20000):
    d = rng.integers(2, 6)
    H, E = random_H(d, rng)
    A = random_A(d, rng)
    rho = random_state(d, rng, pure=(rng.random() < 0.3))
    diffs = np.abs(E[:, None] - E[None, :])
    wpos = diffs[diffs > 1e-6]
    if wpos.size == 0:
        continue
    T = rng.uniform(0.05, 2.0) * 2 * np.pi / wpos.min()
    U = expm_herm(H, T)
    rhoT = U @ rho @ U.conj().T
    # links
    Delta = abs(np.real(np.trace(A @ (rhoT - rho))))
    sA = sigma_A(A)
    D = trace_distance(rho, rhoT)
    # link a: Delta <= 2 sA D
    la = 2 * sA * D - Delta
    # link b: D <= sum_k p_k sqrt(1-F_k^2)
    p, V = np.linalg.eigh(rho)
    rhs_b = 0.0
    for k in range(d):
        ket = V[:, k]
        Fk = abs(ket.conj() @ (U @ ket))
        rhs_b += p[k] * np.sqrt(max(0.0, 1 - Fk**2))
    lb = rhs_b - D
    # link c (per eigenstate): 1-Fk <= K T <k|H|k>
    lc = np.inf
    for k in range(d):
        ket = V[:, k]
        Fk = abs(ket.conj() @ (U @ ket))
        ekH = np.real(ket.conj() @ (H @ ket))
        lc = min(lc, K * T * ekH - (1 - Fk))
    # overall: Delta^2 <= 8 K sA^2 T <H>
    meanH = np.real(np.trace(rho @ H))
    overall = 8 * K * sA**2 * T * meanH - Delta**2
    for name, slack in [("a", la), ("b", lb), ("c", lc), ("overall", overall)]:
        if slack < worst_slack:
            worst_slack = slack
        if slack < -1e-9:
            print(f"  !! link {name} VIOLATED: slack={slack:.3e} (trial {trial}, d={d})")
print(f"  worst slack across all links/trials = {worst_slack:.3e}  (>=0 means all hold)")
print(f"  => convexity chain {'HOLDS on all instances' if worst_slack >= -1e-9 else 'FAILS'}")
