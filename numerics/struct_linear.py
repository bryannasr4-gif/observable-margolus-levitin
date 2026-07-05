"""
P2.1 LINEAR-RECOVERY test.

For GROUND-GAPPED A (A_0n = 0, A acts within the excited block), claim:
  (i)  instantaneous LINEAR mean-energy bound:  <H>-E0 >= (Delta_gap/2) * Delta/sigma_A,
       i.e. <H> >= kappa * Delta/sigma_A with kappa = Delta_gap/2  (Delta_gap = lowest excited gap).
       Mechanism: ||P_e psi||^2 = eps conserved => Delta <= 2 sigma_A eps, and <H> >= Delta_gap eps.
  (ii) the SPEED-LIMIT scaling P=T(<H>-E0) vs Delta is LINEAR (slope 1), not quadratic,
       in the ground-dominated regime -> the no-go (quadratic) law is BROKEN by the gap.

Contrast: ground-TOUCHING A gives slope 2 (quadratic, C_true).

We verify (i) as an inequality over random ground-gapped configs, and (ii) by the
log-log slope of P_min(Delta).
"""
import numpy as np, mpmath as mp
mp.mp.dps=30
xstar=mp.findroot(lambda x: mp.tan(x/2)-x, mp.mpf('2.33')); K=(1-mp.cos(xstar))/xstar
C_true=float(1/(8*K))
print(f"C_true={C_true:.10f}")
print("="*72)
rng=np.random.default_rng(11)

def evolve_change(E,A,c,Tg):
    E=np.asarray(E,float); E=E-E.min(); A=np.asarray(A,complex); A=(A+A.conj().T)/2
    c=np.asarray(c,complex); c=c/np.linalg.norm(c)
    M=np.outer(np.conj(c),c)*A; W=E[:,None]-E[None,:]
    traj=np.real(np.einsum('mn,tmn->t',M,np.exp(1j*W[None]*Tg[:,None,None])))
    return traj, traj[0]

# ---- (i) instantaneous linear bound for ground-gapped A
print("\n(i) ground-gapped A:  <H> >= (Delta_gap/2)*Delta/sigma_A  over random configs")
worst=np.inf  # min of (<H> / ((Dgap/2)*Delta/sA)) ; must be >= 1
for _ in range(40000):
    d=rng.integers(2,6)
    E=np.zeros(d); E[1:]=np.sort(np.abs(rng.normal(scale=2,size=d-1)))+rng.uniform(0.2,1.5)
    Dgap=E[1]-E[0]
    # A in excited block only
    A=np.zeros((d,d),complex)
    for i in range(1,d):
        for j in range(i+1,d):
            A[i,j]=rng.normal()+1j*rng.normal(); A[j,i]=np.conj(A[i,j])
    if np.allclose(A,0): continue
    w=np.linalg.eigvalsh(A); sA=(w[-1]-w[0])/2
    if sA<1e-9: continue
    c=rng.normal(size=d)+1j*rng.normal(size=d); c[0]=abs(c[0])+rng.uniform(0,3)  # ground weight
    c=c/np.linalg.norm(c)
    mH=float(np.sum(E*np.abs(c)**2))
    # max |Delta| over time
    Tg=np.linspace(0, 8*2*np.pi/max(Dgap,1e-3), 4000)
    traj,A0=evolve_change(E,A,c,Tg); Delta=np.abs(traj-A0).max()
    if Delta<1e-9: continue
    lhs=mH; rhs=(Dgap/2)*Delta/sA
    worst=min(worst, lhs/rhs)
print(f"   min( <H> / [(Dgap/2)Delta/sigma_A] ) = {worst:.6f}   (LINEAR bound holds iff >= 1: {worst>=1-1e-6})")

# ---- (ii) slope of P_min vs Delta: ground-gapped (linear) vs ground-touching (quadratic)
print("\n(ii) log-log slope of P=T<H> vs Delta near Delta->0")
def Pmin_vs_Delta(ground_touch, E1=1.0, E2=2.3):
    # 3-level. ground-touch: A couples 0<->1. gapped: A couples 1<->2.
    Ev=np.array([0.0,E1,E2])
    if ground_touch:
        A=np.array([[0,1,0],[1,0,0],[0,0,0]],complex)
        coupled=(0,1)
    else:
        A=np.array([[0,0,0],[0,0,1],[0,1,0]],complex)
        coupled=(1,2)
    sA=1.0
    Ds=[]; Ps=[]
    for eps in np.geomspace(1e-4,1e-2,12):
        # put the excited weight into the coupled pair, optimal phase, scan
        best=(np.inf,None)
        for b in np.linspace(0,2*np.pi,40,endpoint=False):
            c=np.zeros(3,complex); c[0]=np.sqrt(1-eps)
            if ground_touch:
                c[1]=np.sqrt(eps)*np.exp(1j*b)
            else:
                c[1]=np.sqrt(eps/2); c[2]=np.sqrt(eps/2)*np.exp(1j*b)
            mH=float(np.sum(Ev*np.abs(c)**2))
            wgap=abs(Ev[coupled[0]]-Ev[coupled[1]])
            Tg=np.linspace(1e-9, 4*2*np.pi/wgap, 80000)
            traj,A0=evolve_change(Ev,A,c,Tg); ch=np.abs(traj-A0); cmax=ch.max()
            if cmax<1e-12: continue
            # first passage to a fixed fraction; take the min-P first turning
            idx=np.argmax(ch>=0.5*cmax)
            if ch[idx]<0.5*cmax or idx==0: continue
            T=Tg[idx]; Delta=ch[idx]; P=T*mH
            if P<best[0]: best=(P,(Delta,P))
        if best[1] is not None:
            Ds.append(best[1][0]); Ps.append(best[1][1])
    Ds=np.array(Ds); Ps=np.array(Ps)
    slope=np.polyfit(np.log(Ds),np.log(Ps),1)[0]
    return slope, Ds, Ps
s_gt,_,_=Pmin_vs_Delta(True)
s_gp,_,_=Pmin_vs_Delta(False)
print(f"   ground-TOUCHING A: slope(log P vs log Delta) = {s_gt:.4f}   (expect ~2, quadratic)")
print(f"   ground-GAPPED   A: slope(log P vs log Delta) = {s_gp:.4f}   (expect ~1, LINEAR)")
print("\n   => the GAP turns the quadratic no-go into a LINEAR mean-energy bound.")
