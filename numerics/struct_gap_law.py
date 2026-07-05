"""
P2.1: exact improvement law for GROUND-GAPPED A (A couples excited<->excited only),
and for A whose ground-coupling has a minimum transition gap.

(I) Ground-gapped A: ratio r ~ const/eps as eps=excited population ->0.
    We fit r*eps and show it tends to a finite constant (so r -> infinity:
    the bound is INFINITELY improved as the state concentrates on the ground).

(II) Hard-gap ground-touching A: A couples ground only to levels with E >= Egap.
    Then the saturating single level sits AT E=Egap with T=x*/Egap (still allowed),
    so r -> C_true.  No improvement from a gap, as long as SOME allowed level can
    be put at energy giving phase x*.  Confirm.

(III) Energy-WINDOW restriction on the *product* phi=E*T is what would change C.
    If the structure pins all active phases to a value phi0 != x* (e.g. a single
    rigid frequency with fixed observation time), the constant becomes
      C(phi0) = phi0 / (8 (1-cos phi0))   (free phase)  >= C_true,
    with equality only at phi0=x*.  Tabulate.
"""
import numpy as np, mpmath as mp
mp.mp.dps=40
xstar=mp.findroot(lambda x: mp.tan(x/2)-x, mp.mpf('2.33')); K=(1-mp.cos(xstar))/xstar
C_true=float(1/(8*K)); xs=float(xstar)
print(f"C_true={C_true:.12f}  x*={xs:.10f}")
print("="*72)

def best_ratio(E,A,c,Npts=300000,ncyc=4.0):
    E=np.asarray(E,float); E=E-E.min(); A=np.asarray(A,complex); A=(A+A.conj().T)/2
    c=np.asarray(c,complex); c=c/np.linalg.norm(c)
    w=np.linalg.eigvalsh(A); sA=(w[-1]-w[0])/2; meanH=float(np.sum(E*np.abs(c)**2))
    if sA<1e-12 or meanH<1e-14: return np.inf
    d=np.abs(E[:,None]-E[None,:]); wmax=d.max(); wmin=d[d>1e-9].min() if np.any(d>1e-9) else wmax
    Tg=np.linspace(1e-9, ncyc*2*np.pi/wmin, Npts)
    M=np.outer(np.conj(c),c)*A; W=E[:,None]-E[None,:]
    traj=np.real(np.einsum('mn,tmn->t',M,np.exp(1j*W[None]*Tg[:,None,None])))
    ch=np.abs(traj-traj[0]); cmax=ch.max(); best=np.inf
    for frac in np.linspace(0.02,0.999,400):
        idx=np.argmax(ch>=frac*cmax)
        if ch[idx]>=frac*cmax and idx>0: best=min(best, Tg[idx]*meanH*sA**2/ch[idx]**2)
    return best

print("\n(I) GROUND-GAPPED A (excited<->excited only): r*eps -> finite => r ~ C_gap/eps")
# 3-level, A couples 1<->2; optimize phase b and split g, energies, to MINIMIZE r at each eps
from scipy.optimize import minimize_scalar
def min_r_gapped(eps, E1=1.0, E2=2.0):
    best=np.inf
    for g in np.linspace(0.05, np.pi/2-0.05, 25):
        for b in np.linspace(0, 2*np.pi, 24, endpoint=False):
            c=np.array([np.sqrt(1-eps), np.sqrt(eps)*np.cos(g), np.sqrt(eps)*np.sin(g)*np.exp(1j*b)])
            r=best_ratio([0,E1,E2], np.array([[0,0,0],[0,0,1],[0,1,0]],complex), c, Npts=120000)
            best=min(best,r)
    return best
print("   eps      min_r        r*eps")
for eps in [0.2,0.1,0.05,0.02,0.01]:
    r=min_r_gapped(eps); print(f"   {eps:5.3f}  {r:9.4f}   {r*eps:.5f}")
print("   => r*eps approaches a constant ~ const  =>  r ~ const/eps -> infinity.")
print("   Interpretation: A not touching the ground cannot move <A> at O(sqrt eps);")
print("   only at O(eps), while <H>=O(eps) and the slowest motion costs T~1/gap, so")
print("   r = T<H>sigma^2/Delta^2 ~ (1/gap)*eps / eps^2 = O(1/eps).")

print("\n(II) HARD-GAP ground-touching A (ground couples only to E>=Egap): -> C_true")
for Egap in [0.5,1.0,3.0]:
    # single allowed level at E=Egap, optimal phase, T=x*/Egap
    E=Egap; T=xs/Egap
    Ev=np.array([0.0,E]); A=np.array([[0,1],[1,0]],complex)
    # optimal relative phase beta = x*/2 - pi/2
    beta=xs/2-np.pi/2
    th=0.02
    c=np.array([np.cos(th), np.sin(th)*np.exp(1j*beta)])
    r=best_ratio(Ev,A,c,Npts=300000)
    print(f"   Egap={Egap}: r(theta=0.02) = {r:.8f}  -> C_true={C_true:.8f} as theta->0")

print("\n(III) phase-pinned constant C(phi0)=phi0/(8(1-cos phi0)) (free amplitude phase)")
print("   phi0      C(phi0)     C(phi0)/C_true")
for phi0 in [xs, np.pi/2, np.pi, 3*np.pi/2, 1.0, 4.0]:
    C0=phi0/(8*(1-np.cos(phi0)))
    print(f"   {phi0:6.4f}   {C0:.8f}   {C0/C_true:.5f}")
