"""
P2.1 structured-A investigation.

Question: for A that couples the ground state to the excited space in a rank-1 way,
  A = |0><v| + |v><0|,  v unit, v _|_ |0>  (so spec(A) = {+1, -1, 0,...,0}, sigma_A=1),
or energy-banded A (A connects only levels within an energy window), does the
mean-energy speed-limit constant improve below

  C_true = 1/(8K) = 0.17250626746116262 ?

ratio r := T*(<H>-E0)*sigma_A^2 / Delta^2  ;  proven floor r >= C_true (general A).

We compute the *infimum* of r over each structured class numerically and compare to
C_true, and we identify the saturating configuration.

Conventions: hbar=1, E0=0 (energies shifted so min=0). sigma_A=(lmax-lmin)/2.
We use FIRST-PASSAGE / turning-point ratio (same as optimize_constant.py) but here,
because we control the structure analytically, we ALSO use the exact reductions.
"""
import numpy as np
import mpmath as mp
mp.mp.dps = 40

# ------------------------------------------------------------------ constants
xstar = mp.findroot(lambda x: mp.tan(x/2) - x, mp.mpf('2.33'))
K     = (1 - mp.cos(xstar))/xstar
C_true = float(1/(8*K))
print(f"C_true = {C_true:.15f}   K={float(K):.12f}  x*={float(xstar):.10f}")
print("="*72)

# ------------------------------------------------------------------ evaluator
def best_ratio(E, A, c, Npts=200000, ncycles=3.0):
    """min over first-passage levels of T*<H>*sigmaA^2/Delta^2."""
    E = np.asarray(E, float); E = E - E.min()
    A = np.asarray(A, complex); A = (A + A.conj().T)/2
    c = np.asarray(c, complex); c = c/np.linalg.norm(c)
    w = np.linalg.eigvalsh(A); sA = (w[-1]-w[0])/2
    meanH = float(np.sum(E*np.abs(c)**2))
    if sA < 1e-12 or meanH < 1e-14: return np.inf, None
    diffs = np.abs(E[:,None]-E[None,:])
    wmax = diffs.max(); wmin = diffs[diffs>1e-9].min() if np.any(diffs>1e-9) else wmax
    Tmax = ncycles*2*np.pi/wmin
    Tg = np.linspace(1e-12, Tmax, Npts)
    M = np.outer(np.conj(c), c)*A
    Wm = E[:,None]-E[None,:]
    phase = np.exp(1j*Wm[None,:,:]*Tg[:,None,None])
    traj = np.real(np.einsum('mn,tmn->t', M, phase))
    A0 = traj[0]
    ch = np.abs(traj - A0)
    cmax = ch.max(); best=np.inf; binfo=None
    for frac in np.linspace(0.01,0.999,300):
        lvl = frac*cmax
        idx = np.argmax(ch>=lvl)
        if ch[idx]>=lvl and idx>0:
            T=Tg[idx]; Delta=ch[idx]
            r = T*meanH*sA**2/Delta**2
            if r<best: best=r; binfo=dict(T=T,Delta=Delta,meanH=meanH,sigmaA=sA,frac=frac)
    return best, binfo

# ================================================================= CLASS 1
# rank-1 ground coupling, single excited level reached: A = |0><1|+|1><0| in 2-level
# but EMBEDDED so the ground couples to a normalized v that is a SINGLE eigenstate
# of H with energy E1.  This is exactly the 2-level saturating family.
print("\n[1] rank-1 ground coupling, v = single H-eigenstate (energy E1)")
print("    A = |0><1| + |1><0|, state ground-dominated; sweep theta, E, relative phase")
best1 = np.inf; where1=None
for E1 in [0.5,1.0,2.0,5.0]:
    for theta in np.linspace(0.0008, 0.5, 200):
        for beta in np.linspace(0, 2*np.pi, 60, endpoint=False):
            Ev = np.array([0.0, E1])
            A  = np.array([[0,1],[1,0]], complex)
            c  = np.array([np.cos(theta), np.sin(theta)*np.exp(1j*beta)], complex)
            r,info = best_ratio(Ev,A,c,Npts=4000)
            if r<best1: best1=r; where1=(E1,theta,beta,info)
print(f"    min ratio (class 1) = {best1:.10f}   (C_true={C_true:.10f})")
print(f"    at E1={where1[0]}, theta={where1[1]:.4f}, beta={where1[2]:.4f}")

# ================================================================= CLASS 2
# rank-1 ground coupling where v is a SUPERPOSITION of several excited levels.
# A = |0><v| + |v><0|, v = sum_{k>=1} v_k |k>, H=diag(0,E1,...,E_{d-1}).
# sigma_A = |v| = 1 still (eigs +-|v|).  state: ground + small bit in v-direction etc.
print("\n[2] rank-1 ground coupling, v = superposition of MANY excited levels")
print("    A = |0><v|+|v><0|; H spreads excited energies; full state search")
from scipy.optimize import minimize
import warnings; warnings.filterwarnings("ignore")
def class2_obj(x, d, Evec):
    # params: v (d-1 reals, >=1 component), state c (2d reals)
    v = x[:d-1]
    nv = np.linalg.norm(v)
    if nv < 1e-9: return 10.0
    v = v/nv
    A = np.zeros((d,d), complex)
    A[0,1:] = v; A[1:,0] = v
    cr = x[d-1:d-1+d]; ci = x[d-1+d:d-1+2*d]
    c = cr + 1j*ci
    r,_ = best_ratio(Evec, A, c, Npts=3000)
    return r if np.isfinite(r) else 10.0
best2=np.inf; where2=None
rng = np.random.default_rng(7)
for d in [3,4,5]:
    for trial in range(40):
        Evec = np.zeros(d); Evec[1:] = np.sort(np.abs(rng.normal(scale=2,size=d-1)))+0.1
        n = (d-1) + 2*d
        x0 = rng.normal(size=n)
        x0[d-1] = 1.0  # ground amp real large
        x0[d-1+1:d-1+d] *= 0.2; x0[d-1+d:] *= 0.2
        res = minimize(class2_obj, x0, args=(d,Evec), method='Nelder-Mead',
                       options=dict(maxiter=4000, xatol=1e-6, fatol=1e-9))
        if res.fun < best2:
            best2 = res.fun; where2=(d, Evec.copy(), res.x.copy())
print(f"    min ratio (class 2) = {best2:.10f}   (C_true={C_true:.10f})")
print(f"    d={where2[0]}, E={np.round(where2[1],3)}")

# ================================================================= CLASS 3
# energy-banded A: A connects only levels within an energy window [Elo,Ehi].
# Worst-case intuition: still reduces to a 2-level inside the band where the
# band touches the GROUND, so the smallest Bohr gap available is the band width.
print("\n[3] energy-banded A (A_mn != 0 only if |E_m - E_n| <= W)")
print("    does NOT touch ground unless band includes E0; test both")
def class3_obj(x, d, Evec, W):
    # build Hermitian A with zeros outside band, then full state
    A = np.zeros((d,d), complex)
    idx=0
    for i in range(d):
        for j in range(i+1,d):
            if abs(Evec[i]-Evec[j]) <= W + 1e-12:
                A[i,j] = x[idx]+1j*x[idx+1]; A[j,i]=np.conj(A[i,j]); idx+=2
    # diag
    # (diag doesn't move <A> and is allowed; skip)
    if np.allclose(A,0): return 10.0
    nfree = idx
    c = x[nfree:nfree+d] + 1j*x[nfree+d:nfree+2*d]
    r,_ = best_ratio(Evec, A, c, Npts=3000)
    return r if np.isfinite(r) else 10.0
best3=np.inf; where3=None
for d in [3,4,5]:
    for W in [0.3, 0.8, 1.5]:
        for trial in range(25):
            Evec = np.zeros(d); Evec[1:] = np.sort(np.abs(rng.normal(scale=2,size=d-1)))+0.1
            # count band pairs
            npairs = sum(1 for i in range(d) for j in range(i+1,d) if abs(Evec[i]-Evec[j])<=W+1e-12)
            if npairs==0: continue
            n = 2*npairs + 2*d
            x0 = rng.normal(size=n); x0[2*npairs]=1.0
            x0[2*npairs+1:2*npairs+d]*=0.2; x0[2*npairs+d:]*=0.2
            res = minimize(class3_obj, x0, args=(d,Evec,W), method='Nelder-Mead',
                           options=dict(maxiter=4000, xatol=1e-6, fatol=1e-9))
            if res.fun<best3: best3=res.fun; where3=(d,W,Evec.copy())
print(f"    min ratio (class 3, banded) = {best3:.10f}   (C_true={C_true:.10f})")
print(f"    d={where3[0]}, W={where3[1]}, E={np.round(where3[2],3)}")

print("\n" + "="*72)
print("SUMMARY")
print(f"  C_true                       = {C_true:.12f}")
print(f"  class1 (rank-1, single lvl)  = {best1:.12f}   ratio/C_true = {best1/C_true:.6f}")
print(f"  class2 (rank-1, multi lvl)   = {best2:.12f}   ratio/C_true = {best2/C_true:.6f}")
print(f"  class3 (energy-banded)       = {best3:.12f}   ratio/C_true = {best3/C_true:.6f}")
