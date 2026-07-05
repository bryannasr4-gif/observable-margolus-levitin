"""
P2.1 FAST targeted verification of the analytic claims.

Claim A (rank-1 ground coupling does NOT improve): inf over multi-frequency
   r = T*Ebar / (4 (Re sum_k a_k v_k^* (e^{-iE_k T}-1))^2)  equals C_true,
   and the saturating config is single-frequency (all active E_k equal x*/T),
   a ∝ v, optimal phase.  PROOF: |D| <= sqrt(2 K T Ebar)  => r >= 1/(8K).

Claim B (real-coupling restriction): if a_k,v_k REAL (no phase freedom) the
   infimum rises to C_old = 0.18555147... (single-level beta=0).

Claim C (ground-GAPPED A, i.e. A does NOT couple to the ground state):
   in the ground-dominated limit the ratio DIVERGES (no first-order motion),
   so the constant is unbounded-improved; quantify via the EXACT min-gap law.

Claim D (energy-banded, ground touching, allowed transition energies in [Elo,Ehi]):
   a single level inside the band still saturates -> no improvement (band that
   includes a ground-coupling does not help).

We verify A,B by random sampling of the analytic ratio (fast), and the tightness
of |D| <= sqrt(2 K T Ebar).  We verify C with an exact 3-level computation.
"""
import numpy as np
import mpmath as mp
mp.mp.dps = 40
xstar = mp.findroot(lambda x: mp.tan(x/2)-x, mp.mpf('2.33'))
K = (1-mp.cos(xstar))/xstar
C_true = float(1/(8*K)); C_old = 0.18555147172183604
xstar_f = float(xstar)
print(f"C_true={C_true:.15f}  C_old={C_old:.15f}  x*={xstar_f:.12f}")
print("="*72)
rng = np.random.default_rng(2025)

def ratio_and_bound(a, v, E, T):
    a = a/np.linalg.norm(a); v=v/np.linalg.norm(v)
    Ebar = np.sum(np.abs(a)**2*E)
    w = a*np.conj(v)
    D = np.real(np.sum(w*(np.exp(-1j*E*T)-1.0)))
    r = T*Ebar/(4*D**2) if abs(D)>1e-14 else np.inf
    # analytic upper bound on |D|: sqrt(2 K T Ebar)
    Dbound = np.sqrt(2*float(K)*T*Ebar)
    return r, D, Dbound, Ebar

# ---- (A1) verify |D| <= sqrt(2 K T Ebar) over random configs => r >= C_true
print("\n[A1] verify  |D| <= sqrt(2 K T Ebar)  (=> r >= C_true) over 200k random configs")
worst_violation = -np.inf; min_r = np.inf
for _ in range(200000):
    Kn = rng.integers(1,5)
    a = rng.normal(size=Kn)+1j*rng.normal(size=Kn)
    v = rng.normal(size=Kn)+1j*rng.normal(size=Kn)
    E = np.abs(rng.normal(scale=2,size=Kn))+0.05
    T = rng.uniform(0.05, 10.0)
    r,D,Db,Eb = ratio_and_bound(a,v,E,T)
    worst_violation = max(worst_violation, abs(D)-Db)   # must be <= 0
    if np.isfinite(r): min_r = min(min_r, r)
print(f"   max(|D| - sqrt(2KT Ebar)) = {worst_violation:.3e}  (must be <= 0)")
print(f"   min ratio sampled         = {min_r:.10f}   (>= C_true={C_true:.10f}? {min_r>=C_true-1e-9})")

# ---- (A2) saturating single-frequency free-phase config hits C_true
print("\n[A2] single-frequency, a=v=1, optimal phase, E*T=x*  ->  r = C_true")
# put the optimal phase into v so that Re[w(e^{-i x*}-1)] = |e^{-ix*}-1| = 2 sin(x*/2)
E=np.array([1.0]); T=xstar_f
# choose w phase to maximize Re[w(e^{-iET}-1)]: w = conj(e^{-iET}-1)/|...|
z = np.exp(-1j*xstar_f)-1.0
# want D = Re[ (a v*) z ] maximal = |z|, i.e. a v* = conj(z)/|z|.  a=1 => v = z/|z|.
a=np.array([1.0+0j]); v=np.array([z/abs(z)])
r,D,Db,Eb = ratio_and_bound(a,v,E,T)
print(f"   r = {r:.12f}   (C_true={C_true:.12f})   match: {abs(r-C_true)<1e-6}")

# ---- (B) real couplings -> infimum = C_old (scan single + 2-level real)
print("\n[B] real a,v: single-level scan of r(phi)=phi/(4(cos phi -1)^2) -> C_old")
phis = np.linspace(0.01, 2*np.pi-0.01, 2_000_00)
rphi = phis/(4*(np.cos(phis)-1.0)**2)
print(f"   single-level real inf = {rphi.min():.12f}  at phi={phis[np.argmin(rphi)]:.6f}  (C_old={C_old:.12f})")
# 2-level real multi-freq random search for anything below C_old
bb=np.inf
for _ in range(300000):
    a=rng.normal(size=2); v=rng.normal(size=2); E=np.abs(rng.normal(scale=2,size=2))+0.05
    T=rng.uniform(0.05,12)
    r,_,_,_=ratio_and_bound(a.astype(complex),v.astype(complex),E,T)
    if np.isfinite(r): bb=min(bb,r)
print(f"   real 2-level random inf = {bb:.10f}  (>= C_old? {bb>=C_old-1e-4};  >= C_true? {bb>=C_true})")

# ---- (C) ground-GAPPED A (A couples excited<->excited only): exact 3-level
print("\n[C] ground-GAPPED A (no ground coupling): ratio diverges as eps->0")
# H=diag(0,E1,E2); A couples |1><2| only.  state = sqrt(1-eps)|0> + sqrt(eps)(cos g|1>+ e^{ib} sin g|2>)
def gapped_ratio(eps, E1=1.0, E2=2.5, g=np.pi/4, b=0.0, Npts=200000):
    Ev=np.array([0.0,E1,E2]);
    A=np.zeros((3,3),complex); A[1,2]=1; A[2,1]=1
    c=np.array([np.sqrt(1-eps), np.sqrt(eps)*np.cos(g), np.sqrt(eps)*np.sin(g)*np.exp(1j*b)])
    sA=1.0  # eigs of A block are +-1
    meanH=np.sum(Ev*np.abs(c)**2)
    Tg=np.linspace(1e-9, 6*2*np.pi/abs(E2-E1), Npts)
    M=np.outer(np.conj(c),c)*A; W=Ev[:,None]-Ev[None,:]
    traj=np.real(np.einsum('mn,tmn->t', M, np.exp(1j*W[None]*Tg[:,None,None])))
    ch=np.abs(traj-traj[0]); cmax=ch.max(); best=np.inf
    for frac in np.linspace(0.02,0.999,300):
        idx=np.argmax(ch>=frac*cmax)
        if ch[idx]>=frac*cmax and idx>0:
            r=Tg[idx]*meanH*sA**2/ch[idx]**2; best=min(best,r)
    return best
for eps in [0.3,0.1,0.03,0.01]:
    print(f"   eps={eps:5.3f}: min ratio = {gapped_ratio(eps):.6f}   (C_true={C_true:.4f}; grows ~ 1/eps)")

# ---- (D) energy-banded ground-touching: single level in band saturates -> no improvement
print("\n[D] ground-touching banded A: single level in band -> still C_true (no improvement)")
# trivially: pick E in band, T=x*/E, optimal phase -> r=C_true (covered by A2). Confirm:
for Eband in [0.5, 1.7, 4.0]:
    T=xstar_f/Eband; E=np.array([Eband])
    z=np.exp(-1j*Eband*T)-1; a=np.array([1+0j]); v=np.array([z/abs(z)])
    r,_,_,_=ratio_and_bound(a,v,E,T)
    print(f"   band level E={Eband}: r={r:.10f}  (=C_true)")
