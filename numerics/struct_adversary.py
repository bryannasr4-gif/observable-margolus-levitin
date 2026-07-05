"""
P2.1 ADVERSARIAL: full real-time dynamics (no perturbative reduction).
Try hard to push the ratio r = T*<H>*sigma_A^2/Delta^2 BELOW C_true using
rank-1 ground-coupling A = |0><v|+|v><0| with v a superposition of MANY excited
levels at GENERIC energies, full (not ground-dominated) states, optimal phases.
If the floor C_true is real, nothing should dip below it. Flushed output.
"""
import sys, numpy as np
from scipy.optimize import minimize
import warnings; warnings.filterwarnings("ignore")
import mpmath as mp
mp.mp.dps=30
xstar=mp.findroot(lambda x: mp.tan(x/2)-x, mp.mpf('2.33')); K=(1-mp.cos(xstar))/xstar
C_true=float(1/(8*K))
print(f"C_true={C_true:.12f}", flush=True)

def best_ratio(E,A,c,Npts=60000,ncyc=3.0):
    E=np.asarray(E,float); E=E-E.min(); A=np.asarray(A,complex); A=(A+A.conj().T)/2
    c=np.asarray(c,complex); c=c/np.linalg.norm(c)
    w=np.linalg.eigvalsh(A); sA=(w[-1]-w[0])/2; mH=float(np.sum(E*np.abs(c)**2))
    if sA<1e-12 or mH<1e-14: return np.inf
    d=np.abs(E[:,None]-E[None,:]); wmax=d.max(); wmin=d[d>1e-9].min() if np.any(d>1e-9) else wmax
    Tg=np.linspace(1e-9, ncyc*2*np.pi/wmin, Npts)
    M=np.outer(np.conj(c),c)*A; W=E[:,None]-E[None,:]
    traj=np.real(np.einsum('mn,tmn->t',M,np.exp(1j*W[None]*Tg[:,None,None])))
    ch=np.abs(traj-traj[0]); cmax=ch.max(); best=np.inf
    if cmax<1e-12: return np.inf
    for frac in np.linspace(0.02,0.999,250):
        idx=np.argmax(ch>=frac*cmax)
        if ch[idx]>=frac*cmax and idx>0: best=min(best, Tg[idx]*mH*sA**2/ch[idx]**2)
    return best

def obj(x,d,Evec):
    # v: d-1 complex (excited coupling, unit), c: d complex (full state, any pop)
    vr=x[:d-1]; vi=x[d-1:2*(d-1)]
    v=vr+1j*vi; nv=np.linalg.norm(v)
    if nv<1e-9: return 10.0
    v=v/nv
    A=np.zeros((d,d),complex); A[0,1:]=v; A[1:,0]=np.conj(v)
    cr=x[2*(d-1):2*(d-1)+d]; ci=x[2*(d-1)+d:2*(d-1)+2*d]
    c=cr+1j*ci
    r=best_ratio(Evec,A,c,Npts=40000)
    return r if np.isfinite(r) else 10.0

gmin=np.inf; gwhere=None
for d in [2,3,4,5]:
    rng=np.random.default_rng(900+d)
    nstart={2:120,3:90,4:60,5:45}[d]
    bestd=np.inf; bx=None; bE=None
    for s in range(nstart):
        Evec=np.zeros(d); Evec[1:]=np.sort(np.abs(rng.normal(scale=2,size=d-1)))+0.1
        n=2*(d-1)+2*d
        x0=rng.normal(size=n)
        # bias toward ground-dominated (where the floor is approached)
        sc=rng.choice([0.05,0.15,0.4,1.0])
        x0[2*(d-1)]=1.0  # Re c0
        x0[2*(d-1)+1:2*(d-1)+d]*=sc; x0[2*(d-1)+d:]*=sc
        res=minimize(obj,x0,args=(d,Evec),method='Nelder-Mead',
                     options=dict(maxiter=6000,xatol=1e-8,fatol=1e-11))
        if res.fun<bestd: bestd=res.fun; bx=res.x; bE=Evec.copy()
    # high-res re-eval
    if bx is not None:
        vr=bx[:d-1]; vi=bx[d-1:2*(d-1)]; v=(vr+1j*vi); v/=np.linalg.norm(v)
        A=np.zeros((d,d),complex); A[0,1:]=v; A[1:,0]=np.conj(v)
        c=bx[2*(d-1):2*(d-1)+d]+1j*bx[2*(d-1)+d:2*(d-1)+2*d]
        rhi=best_ratio(bE,A,c,Npts=300000)
        bestd=min(bestd,rhi)
    flag="  <-- BELOW C_true !!" if bestd<C_true-1e-6 else ""
    print(f"  d={d}: min r = {bestd:.10f}   r/C_true={bestd/C_true:.6f}{flag}", flush=True)
    if bestd<gmin: gmin=bestd; gwhere=d
print(f"\nGLOBAL min r (rank-1 ground coupling, full dynamics) = {gmin:.10f}", flush=True)
print(f"  C_true = {C_true:.10f}   floor respected: {gmin>=C_true-1e-6}", flush=True)
