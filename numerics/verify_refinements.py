"""
Verification of the three Sec. V/VI refinements that were stated as numerically
checked in the paper but previously lacked a persisted reproducing script:

  (1) Exact trade-off curve  P_star(delta)  [Eq. (curve)] is a STRICT LOWER
      BOUNDARY: over random (H, A, pure rho) in d=2..6, no config's operating
      point (delta, P) with P = T(<H>-E0), delta = Delta/(2 sigma_A) falls below
      P_star(delta).  Also reproduces the endpoints P_star->0 (delta->0),
      P_star->pi/2 (delta->1) and the small-delta slope 4*C_true.

  (2) Bandwidth-resolved constant  C_true(L) = 1/(8 K(L)),
      K(L) = sup_{0<x<=L}(1-cos x)/x = (1-cos L)/L for L<x*, else K.
      Violation search over spectrally bounded H (spec in [0,B]); reproduces the
      quoted gains  C_true(L)/C_true = 1.58x at L=1,  2.96x at L=0.5.

  (3) Two-observable bound  T(<H>-E0) >= C_true (Delta_1^2 + Delta_2^2)/sigma_A^2
      for orthonormal A1=sigma_x, A2=sigma_y; no violation, and the 2x gain over
      the single-observable floor is saturated at Delta_1 = Delta_2.

hbar = 1, E0 = 0.  Seeded.  Run: python -u numerics/verify_refinements.py
"""
import numpy as np, mpmath as mp
mp.mp.dps = 40

xstar = float(mp.findroot(lambda x: mp.tan(x/2) - x, mp.mpf('2.33')))
K = (1 - np.cos(xstar)) / xstar
C_true = 1.0 / (8 * K)
print(f"x*={xstar:.12f}  K={K:.12f}  C_true={C_true:.12f}\n")

def expmH(H, t):
    w, V = np.linalg.eigh(H)
    return (V * np.exp(-1j * w * t)) @ V.conj().T
def sigA(A):
    w = np.linalg.eigvalsh(A); return (w[-1] - w[0]) / 2

# ======================================================================
# (1) exact trade-off curve P_star(delta) as a lower boundary
# ======================================================================
print("=" * 70)
print("(1) exact trade-off curve P_star(delta) is a strict lower boundary")
print("=" * 70)
# parametric curve, tau = E*T in [x*, pi):  sin^2 th = (tau - tan(tau/2))/(tau - 2 tan(tau/2))
taus = np.linspace(xstar + 1e-9, np.pi - 1e-6, 200000)
tt = np.tan(taus / 2)
s2 = (taus - tt) / (taus - 2 * tt)         # sin^2(theta)
s2 = np.clip(s2, 0.0, 1.0)
Pstar = taus * s2
delta = np.sqrt(4 * s2 * (1 - s2)) * np.sin(taus / 2)   # sin(2th)=2 sinth costh -> sin^2(2th)=4 s2(1-s2)
order = np.argsort(delta)
d_tab, P_tab = delta[order], Pstar[order]
def Pstar_of(d):                            # monotone-interp lower envelope
    return np.interp(d, d_tab, P_tab)
print(f"  endpoints: P_star(delta->0)->{P_tab[0]:.3e} (=0),  "
      f"P_star(delta->1)->{P_tab[-1]:.6f}  (pi/2={np.pi/2:.6f})")
# small-delta slope: P_star ~ 4 C_true delta^2
small = d_tab < 0.05
slope = np.polyfit(d_tab[small] ** 2, P_tab[small], 1)[0]
print(f"  small-delta:  P_star/delta^2 -> {slope:.6f}   (4*C_true={4*C_true:.6f})")

rng = np.random.default_rng(4242)
min_slack = np.inf; n = 0
for _ in range(6000):
    d = rng.integers(2, 7)
    E = np.abs(rng.normal(scale=rng.uniform(0.4, 2.5), size=d)); E[0] = 0; E -= E.min()
    H = np.diag(E)
    M = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)); A = (M + M.conj().T) / 2
    A = A / sigA(A)                          # sigma_A = 1
    v = rng.normal(size=d) + 1j * rng.normal(size=d); psi = v / np.linalg.norm(v)
    meanH = np.real(psi.conj() @ (H @ psi))
    if meanH < 1e-9: continue
    A0 = np.real(psi.conj() @ (A @ psi))
    gaps = np.abs(E[:, None] - E[None, :]); gpos = gaps[gaps > 1e-6]
    if gpos.size == 0: continue
    for T in np.linspace(1e-3, 1.5 * 2 * np.pi / gpos.min(), 260):
        pt = expmH(H, T) @ psi
        Delta = abs(np.real(pt.conj() @ (A @ pt)) - A0)
        dd = Delta / 2.0                     # sigma_A = 1 -> delta = Delta/2
        if dd < 1e-4 or dd >= 0.999: continue
        P = T * meanH
        slack = P - Pstar_of(dd)             # must be >= 0
        n += 1
        if slack < min_slack: min_slack = slack
print(f"  configs sampled: {n:,}   min slack (P - P_star(delta)) = {min_slack:.3e}")
print(f"  => {'PASS (curve is a lower boundary)' if min_slack > -1e-6 else 'VIOLATION'}\n")

# ======================================================================
# (2) bandwidth-resolved constant
# ======================================================================
print("=" * 70)
print("(2) bandwidth-resolved constant  C_true(L)=1/(8 K(L))")
print("=" * 70)
def KL(L):
    return (1 - np.cos(L)) / L if L < xstar else K
for L in (1.0, 0.5):
    gain = K / KL(L)                         # C_true(L)/C_true
    print(f"  L={L}:  C_true(L)/C_true = K/K(L) = {gain:.4f}   (paper: {'1.58' if L==1 else '2.96'}x)")
# violation search over spectrally bounded H: spec in [0, B]
rng = np.random.default_rng(909)
viol = 0; trials = 0; min_ratio_over_CL = np.inf
for _ in range(4000):
    d = rng.integers(2, 6)
    B = rng.uniform(0.3, 3.0)
    E = np.sort(rng.uniform(0, B, size=d)); E[0] = 0                    # spec in [0,B]
    H = np.diag(E)
    M = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)); A = (M + M.conj().T) / 2
    A = A / sigA(A)
    v = rng.normal(size=d) + 1j * rng.normal(size=d); psi = v / np.linalg.norm(v)
    meanH = np.real(psi.conj() @ (H @ psi))
    if meanH < 1e-9: continue
    A0 = np.real(psi.conj() @ (A @ psi))
    for T in np.linspace(1e-3, 3.0, 200):
        pt = expmH(H, T) @ psi
        Delta = abs(np.real(pt.conj() @ (A @ pt)) - A0)
        if Delta < 1e-6: continue
        L = B * T
        ratio = T * meanH / Delta ** 2       # sigma_A=1; must be >= C_true(L)=1/(8 K(L))
        CL = 1.0 / (8 * KL(L))
        trials += 1
        min_ratio_over_CL = min(min_ratio_over_CL, ratio / CL)
        if ratio < CL - 1e-9: viol += 1
print(f"  bandwidth violation search: {trials:,} trials, violations={viol}, "
      f"min(ratio/C_true(BT))={min_ratio_over_CL:.5f}")
print(f"  => {'PASS' if viol == 0 else 'VIOLATION'}\n")

# ======================================================================
# (3) two-observable bound
# ======================================================================
print("=" * 70)
print("(3) two-observable bound  T<H> >= C_true (Delta_1^2+Delta_2^2)/sigma_A^2")
print("=" * 70)
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]], complex)
sz = np.array([[1, 0], [0, -1]], complex)
# sigma_A(sx)=sigma_A(sy)=1
rng = np.random.default_rng(31337)
viol = 0; trials = 0; min_ratio = np.inf; best_2x = np.inf
for _ in range(20000):
    h = rng.normal(size=3); H = h[0]*sx + h[1]*sy + h[2]*sz
    H = H - np.linalg.eigvalsh(H)[0] * np.eye(2)
    v = rng.normal(size=2) + 1j * rng.normal(size=2); psi = v / np.linalg.norm(v)
    meanH = np.real(psi.conj() @ (H @ psi))
    if meanH < 1e-9: continue
    r0 = np.outer(psi, psi.conj())
    A1_0 = np.real(np.trace(sx @ r0)); A2_0 = np.real(np.trace(sy @ r0))
    dd = np.diff(np.linalg.eigvalsh(H))[0]
    if dd < 1e-6: continue
    for T in np.linspace(1e-3, 2 * 2 * np.pi / dd, 200):
        U = expmH(H, T); rT = U @ r0 @ U.conj().T
        D1 = abs(np.real(np.trace(sx @ rT)) - A1_0)
        D2 = abs(np.real(np.trace(sy @ rT)) - A2_0)
        s = D1 ** 2 + D2 ** 2
        if s < 1e-10: continue
        ratio = T * meanH / s                # must be >= C_true
        trials += 1
        min_ratio = min(min_ratio, ratio)
        if ratio < C_true - 1e-9: viol += 1
        if abs(D1 - D2) < 0.02 * max(D1, D2):   # near Delta_1=Delta_2: track 2x saturation
            best_2x = min(best_2x, T * meanH / max(D1, D2) ** 2)  # -> 2*C_true if saturated
print(f"  two-observable search: {trials:,} trials, violations={viol}, "
      f"min(T<H>/(D1^2+D2^2))={min_ratio:.6f}  (C_true={C_true:.6f})")
print(f"  at Delta_1~=Delta_2:  min T<H>/max(D)^2 = {best_2x:.6f}   "
      f"(2*C_true={2*C_true:.6f}; the 2x single-observable floor)")
print(f"  => {'PASS' if viol == 0 else 'VIOLATION'}")
