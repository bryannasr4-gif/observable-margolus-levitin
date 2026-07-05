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
       (iii) the EXACT closed-form two-level ratio:
              <A(t)> = A0 + R cos(omega t + phi); first passage to target Delta.
             For a 2-level system with energies {0,E}, state (c0,c1), Hermitian A,
             <A(t)> = a + 2|c0||c1| |A01| cos(E t + delta).  The change from t=0,
             |<A(t)>-<A(0)>|, as a function of t has analytic turning structure;
             ratio = T*<H>*sigmaA^2/Delta^2 with <H>=E|c1|^2, sigmaA=(lmax-lmin)/2.
We print all three; agreement of (ii),(iii) and a value >= C exposes (i) as artifact.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize
from optimize_constant import (best_ratio_for_config, unpack, pack_dim,
                               objective, C, mean_A_trajectory, sigma_A)
np.random.seed(12345)

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
delta = -np.angle(z)
A0 = a_const + R*np.cos(0 - 0*delta)  # at t=0: a_const + 2 Re(z) = a_const + R cos(delta)... careful
A0_exact = a_const + 2*np.real(z)
Ew = E[1]-E[0]
# change(t) = |a_const + R cos(Ew t - arg z) - A0_exact|
def change_exact(t):
    val = a_const + R*np.cos(Ew*t - np.angle(z))
    return abs(val - A0_exact)
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
print(f"C = {C:.10f}")
print()
print("Interpretation:")
print(f"  search(coarse) = {best:.6f}   dense = {dr:.6f}   exact = {bestE:.6f}")
if bestE >= C - 1e-6 and dr >= C - 1e-6:
    print("  => dense+exact are >= C: the sub-C search value is a GRID ARTIFACT (under-resolved).")
else:
    print("  => dense/exact ALSO below C: investigate as a genuine candidate violation.")
