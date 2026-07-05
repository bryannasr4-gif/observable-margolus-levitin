"""
P2.1 CRUX: ground-dominated infimum for rank-1 ground coupling, MULTI-FREQUENCY.

Derivation (see ext_structured.md):  state |psi> = c0|0> + sqrt(eps) sum_k a_k|k>,
A = |0><v|+|v><0|, v = sum_k v_k|k>, ||a||=||v||=1, H=diag(0,E1,..),  eps->0.

  ratio r(T) = T * Ebar / ( 4 * | Re sum_k a_k v_k^* (e^{-i E_k T} - 1) |^2 ),
  Ebar = sum_k |a_k|^2 E_k.

We MINIMIZE r over a_k, v_k (complex, unit norm each), energies E_k>0, time T>0,
for K (= number of excited levels) = 1,2,3,4.  If the multi-frequency infimum is
strictly below C_true, the structured class would IMPROVE the bound is FALSE
(it would VIOLATE the proven floor -> would be a bug).  If it equals C_true, the
constant does NOT improve for rank-1 ground coupling.  If it is > C_true for a
RESTRICTED structure (e.g. real couplings, or banded energies), that restricted
class gets a BETTER (larger) constant.
"""
import numpy as np
import mpmath as mp
from scipy.optimize import minimize
import warnings; warnings.filterwarnings("ignore")
mp.mp.dps = 40
xstar = mp.findroot(lambda x: mp.tan(x/2)-x, mp.mpf('2.33'))
K = (1-mp.cos(xstar))/xstar
C_true = float(1/(8*K))
C_old  = 0.18555147172183604   # phase-restricted (real) single-level value
print(f"C_true (free phase) = {C_true:.15f}")
print(f"C_old  (real/beta=0)= {C_old:.15f}")
print("="*72)

def r_of(params, Kn, complex_couple=True, free_T=True):
    """params packs: a (Kn complex=2Kn reals), v (Kn complex=2Kn reals),
       E (Kn reals, >0), T (1 real). Returns ratio r."""
    i=0
    ar = params[i:i+Kn]; i+=Kn
    ai = params[i:i+Kn]; i+=Kn
    vr = params[i:i+Kn]; i+=Kn
    vi = params[i:i+Kn]; i+=Kn
    E  = np.abs(params[i:i+Kn]); i+=Kn
    T  = abs(params[i])
    a = ar + 1j*ai; v = vr + 1j*vi
    na = np.linalg.norm(a); nv = np.linalg.norm(v)
    if na<1e-9 or nv<1e-9 or T<1e-9: return 10.0
    a/=na; v/=nv
    Ebar = float(np.sum(np.abs(a)**2 * E))
    if Ebar<1e-12: return 10.0
    w = a*np.conj(v)
    denom = np.real(np.sum(w*(np.exp(-1j*E*T)-1.0)))
    if abs(denom)<1e-12: return 10.0
    return T*Ebar/(4*denom**2)

def search(Kn, n_starts=400, complex_couple=True, real_only=False, seed=0):
    rng = np.random.default_rng(seed)
    best=np.inf; bx=None
    npar = 5*Kn+1
    for s in range(n_starts):
        x0 = rng.normal(size=npar)
        x0[4*Kn:5*Kn] = np.abs(rng.normal(scale=1.5,size=Kn))+0.2  # E>0
        x0[-1] = rng.uniform(1.0,4.0)/ (x0[4*Kn]+0.1)              # T ~ scale
        if real_only:
            # zero out imaginary parts of a and v
            x0[Kn:2*Kn]=0; x0[3*Kn:4*Kn]=0
        def obj(x):
            if real_only:
                x=x.copy(); x[Kn:2*Kn]=0; x[3*Kn:4*Kn]=0
            return r_of(x,Kn)
        res = minimize(obj, x0, method='Nelder-Mead',
                       options=dict(maxiter=8000,xatol=1e-9,fatol=1e-12))
        if res.fun<best: best=res.fun; bx=res.x
    return best,bx

print("\n[FREE-PHASE complex couplings]  infimum of r vs #excited-levels Kn")
for Kn in [1,2,3,4]:
    b,_ = search(Kn, n_starts=300 if Kn<=2 else 200, seed=10+Kn)
    flag = "  <-- BELOW C_true (would be VIOLATION/bug)" if b<C_true-1e-7 else ("  == C_true" if abs(b-C_true)<1e-5 else "")
    print(f"   Kn={Kn}: inf r = {b:.12f}   r/C_true={b/C_true:.6f}{flag}")

print("\n[REAL couplings & real amplitudes]  (a_k, v_k real)  infimum of r")
for Kn in [1,2,3,4]:
    b,_ = search(Kn, n_starts=300 if Kn<=2 else 200, real_only=True, seed=50+Kn)
    rel = b/C_old
    note = "  == C_old(0.18555)" if abs(b-C_old)<1e-4 else ("  BELOW C_old" if b<C_old-1e-5 else "")
    print(f"   Kn={Kn}: inf r = {b:.12f}   r/C_true={b/C_true:.6f}  r/C_old={rel:.6f}{note}")
