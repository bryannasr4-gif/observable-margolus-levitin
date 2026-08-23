"""
Audit the suspected sub-C d=2 optimizer result.

For d=2 the dynamics is EXACTLY the two-level problem and the analytic floor C
applies rigorously.  Any reported ratio < C at d=2 must be a numerical artifact
of the trajectory/first-passage evaluator (under-resolved Delta or T).

Here we:
  1. Re-run the SAME d=2 multistart, but SAVE the winning parameter vector.
  2. Reduce the winning A,E,psi to its canonical two-level form and compute the
     ratio THREE ways:
       (i)  the search evaluator (coarse grid),
       (ii) a very-dense-grid first-passage evaluator,
       (iii) the EXACT closed-form two-level ratio, in float64:
              <A(t)> = A0 + R cos(omega t + phi); first passage to target Delta.
             For a 2-level system with energies {0,E}, state (c0,c1), Hermitian A,
             <A(t)> = a + 2|c0||c1| |A01| cos(E t + delta).  The change from t=0,
             |<A(t)>-<A(0)>|, as a function of t has analytic turning structure;
             ratio = T*<H>*sigmaA^2/Delta^2 with <H>=E|c1|^2, sigmaA=(lmax-lmin)/2.
       (iv) the SAME closed form re-evaluated in 60-digit arithmetic (mpmath),
            with the first-passage time solved analytically rather than scanned,
            so that no grid resolution enters at all.
We print all four.  The winning d=2 configuration is a near-ground state whose
excited population is ~1e-11 and whose swing Delta is ~1e-5, so the coarse
first-passage evaluator of the search is the least reliable of the four; (ii),
(iii) and (iv) agreeing above C is what settles the d=2 floor.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import mpmath as mp
from scipy.optimize import minimize
from optimize_constant import (unpack, pack_dim, objective, C,
                               mean_A_trajectory, sigma_A)

d = 2
rng = np.random.default_rng(1002)
nE, nA, nc = pack_dim(d)
n = nE + nA + nc
best = np.inf; bestx = None
for s in range(160):
    x0 = rng.normal(size=n)
    x0[:nE] = np.abs(rng.normal(scale=1.5, size=nE)); x0[0] = 0.0
    x0[nE+nA:nE+nA+nc] *= 0.3
    x0[nE+nA] = 1.0
    res = minimize(objective, x0, args=(d,), method='Nelder-Mead',
                   options=dict(maxiter=3000, xatol=1e-7, fatol=1e-9))
    if res.fun < best:
        best = res.fun; bestx = res.x
print(f"search-evaluator best d=2 ratio = {best:.10f}  (C={C:.10f})", flush=True)

E, A, c = unpack(bestx, d)
E = E - E.min()
A = (A + A.conj().T)/2
c = c/np.linalg.norm(c)
print("E =", E)
print("A =\n", A)
print("c =", c)

sA = sigma_A(A)
meanH = np.real(np.sum(E*np.abs(c)**2))
print(f"sigma_A={sA:.8f}  <H>={meanH:.8f}")

def dense_eval(E, A, c, Npts=400000):
    Ediff = abs(E[1]-E[0])
    Tmax = 6*2*np.pi/Ediff
    Tg = np.linspace(1e-12, Tmax, Npts)
    traj = mean_A_trajectory(E, A, c, Tg)
    A0 = mean_A_trajectory(E, A, c, np.array([0.0]))[0]
    ch = np.abs(traj - A0)
    cmax = ch.max()
    best = np.inf; binfo=None
    for frac in np.linspace(0.01, 0.999, 400):
        lvl = frac*cmax
        idx = np.argmax(ch >= lvl)
        if ch[idx] >= lvl and idx>0:
            T=Tg[idx]; Delta=ch[idx]
            r = T*meanH*sA**2/Delta**2
            if r < best: best=r; binfo=(T,Delta,frac)
    return best, binfo

dr, dinfo = dense_eval(E, A, c)
print(f"dense-grid first-passage ratio   = {dr:.10f}   info(T,Delta,frac)={dinfo}")

# Exact closed form: <A(t)> = a + R cos(E t + delta)
# a = |c0|^2 A00 + |c1|^2 A11 ; R cos delta + ... from cross term 2Re(conj(c0)c1 A01 e^{-iEt})
c0, c1 = c[0], c[1]
a_const = (abs(c0)**2*A[0,0] + abs(c1)**2*A[1,1]).real
z = np.conj(c0)*c1*A[0,1]      # cross amplitude: term = 2 Re(z e^{i(E0-E1)t}) = 2 Re(z e^{-iE t})
R = 2*abs(z)
# <A(t)> = a_const + 2 Re(z e^{-iE t}) = a_const + R cos(E t - arg(z))
A0_exact = a_const + 2*np.real(z)
Ew = E[1]-E[0]
# change(t) = |a_const + R cos(Ew t - arg z) - A0_exact|
# scan exact change for min ratio over first-passage levels
Tg = np.linspace(1e-12, 6*2*np.pi/abs(Ew), 2000000)
val = a_const + R*np.cos(Ew*Tg - np.angle(z))
ch = np.abs(val - A0_exact)
cmax = ch.max()
bestE=np.inf; bestinfo=None
for frac in np.linspace(0.01,0.999,500):
    lvl=frac*cmax
    idx=np.argmax(ch>=lvl)
    if ch[idx]>=lvl and idx>0:
        T=Tg[idx]; Delta=ch[idx]
        r=T*meanH*sA**2/Delta**2
        if r<bestE: bestE=r; bestinfo=(T,Delta,frac)
print(f"EXACT closed-form 2-level ratio  = {bestE:.10f}   info={bestinfo}")

# ----------------------------------------------------------------------
# (iv) the same closed form at 60 digits, first passage solved analytically.
#      <A(t)> - <A(0)> = R[cos(w t - phi) - cos phi],  phi = arg z, w = E1 - E0.
#      For a level L the first passage solves cos(w t - phi) = cos phi +- L/R
#      exactly, so no time grid is involved.
# ----------------------------------------------------------------------
mp.mp.dps = 60
c0m, c1m = mp.mpc(c[0]), mp.mpc(c[1])
A00, A11, A01 = mp.mpf(A[0, 0].real), mp.mpf(A[1, 1].real), mp.mpc(A[0, 1])
a_m = mp.fabs(c0m)**2 * A00 + mp.fabs(c1m)**2 * A11
z_m = mp.conj(c0m) * c1m * A01
R_m = 2 * mp.fabs(z_m)
phi_m = mp.arg(z_m)
w_m = mp.mpf(E[1]) - mp.mpf(E[0])
meanH_m = w_m * mp.fabs(c1m)**2
half_tr, half_df = (A00 + A11) / 2, (A00 - A11) / 2
rad = mp.sqrt(half_df**2 + mp.fabs(A01)**2)
sA_m = ((half_tr + rad) - (half_tr - rad)) / 2          # (lmax - lmin)/2, exact for 2x2
excited_pop = mp.fabs(c1m)**2

def first_passage(level):
    """smallest t > 0 with |R(cos(w t - phi) - cos phi)| = level, or None."""
    best_t = None
    for sgn in (1, -1):
        target = mp.cos(phi_m) + sgn * level / R_m
        if abs(target) > 1:
            continue
        base = mp.acos(target)
        for root in (base, -base):
            for k in (0, 1, 2):
                t = (phi_m + root + 2 * mp.pi * k) / w_m
                if t > 0 and (best_t is None or t < best_t):
                    best_t = t
    return best_t

cmax_m = min(2 * R_m, R_m * (1 + mp.fabs(mp.cos(phi_m))))
best60, info60 = None, None
for i in range(1, 4001):
    lvl = cmax_m * mp.mpf(i) / 4001
    t = first_passage(lvl)
    if t is None:
        continue
    r = t * meanH_m * sA_m**2 / lvl**2
    if best60 is None or r < best60:
        best60, info60 = r, (t, lvl)
print(f"60-digit exact 2-level ratio      = {mp.nstr(best60, 12)}"
      f"   (ratio - C = {mp.nstr(best60 - mp.mpf(C), 6)})")
print(f"    excited population |c1|^2 = {mp.nstr(excited_pop, 6)}, "
      f"swing at the minimum = {mp.nstr(info60[1], 6)}")
print(f"C = {C:.10f}")
print()
print("Interpretation:")
print(f"  search(coarse) = {best:.6f}   dense = {dr:.6f}   exact = {bestE:.6f}"
      f"   exact@60dps = {mp.nstr(best60, 8)}")
ok60 = best60 >= mp.mpf(C) - mp.mpf('1e-9')
if bestE >= C - 1e-6 and dr >= C - 1e-6 and ok60:
    print("  => every resolution-independent evaluation is >= C.  The d=2 floor holds;")
    print("     any sub-C value reported by the coarse search is a GRID ARTIFACT of its")
    print("     first-passage evaluator on this near-ground configuration, not a violation.")
else:
    print("  => a resolution-independent evaluation is ALSO below C: investigate as a")
    print("     genuine candidate violation.")
