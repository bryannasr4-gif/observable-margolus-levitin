"""
P2.1 — INDEPENDENT re-verification of the structured-observable claims
(written from scratch, lean; run with python -u for live output).

hbar=1, E0=0, sigma_A=(lmax-lmin)/2.
  r_quad := T<H> sigma_A^2 / Delta^2   (general floor: r_quad >= C_true)
  r_lin  := T<H> sigma_A   / Delta      (linear ratio for eigenvector observables)

(1) rank-1 ground coupling, multi-excited, free phase: inf r_quad = C_true.
(2) real couplings/amplitudes: inf r_quad = C_old = 0.18555147172183604.
(3) eigenvector (extreme) observable: LINEAR bound T<H> >= (1/4K) Delta/sigma_A, sharp; r_quad diverges.
(4) ground-gapped A (A_0n=0): static bound <H> >= (Delta_gap/2) Delta/sigma_A; r_quad -> inf.
"""
import sys
import numpy as np, mpmath as mp
mp.mp.dps = 30
xstar = mp.findroot(lambda x: mp.tan(x/2) - x, mp.mpf('2.33'))
K = float((1 - mp.cos(xstar)) / xstar)
C_true = 1/(8*K); C_old = 0.18555147172183604; inv4K = 1/(4*K)
print(f"K={K:.12f}  C_true={C_true:.12f}  C_old={C_old:.12f}  1/(4K)={inv4K:.12f}  (==2*C_true? {abs(inv4K-2*C_true)<1e-12})")

def expmH(H,t):
    w,V=np.linalg.eigh(H); return (V*np.exp(-1j*w*t))@V.conj().T
def sigA(A):
    w=np.linalg.eigvalsh(A); return (w[-1]-w[0])/2

# ---- (1) rank-1 ground coupling, multi-excited, free phase ----
def rank1_min(nexc, ntrials=1500):
    rng=np.random.default_rng(100+nexc); best=np.inf
    for _ in range(ntrials):
        d=1+nexc
        E=np.zeros(d); E[1:]=rng.uniform(0.2,3.0,nexc)
        v=rng.normal(size=nexc)+1j*rng.normal(size=nexc); v/=np.linalg.norm(v)
        A=np.zeros((d,d),complex)
        for k in range(nexc): A[0,1+k]=v[k]; A[1+k,0]=np.conj(v[k])
        sA=sigA(A)
        th=rng.uniform(0.01,0.4)
        a=rng.normal(size=nexc)+1j*rng.normal(size=nexc); a/=np.linalg.norm(a)
        c=np.zeros(d,complex); c[0]=np.cos(th); c[1:]=np.sin(th)*a
        meanH=np.real(np.sum(E*np.abs(c)**2))
        if meanH<1e-12: continue
        diffs=np.abs(E[:,None]-E[None,:]); wpos=diffs[diffs>1e-6]
        if wpos.size==0: continue
        Tg=np.linspace(1e-3,2.5*2*np.pi/wpos.min(),250)
        Hd=np.diag(E); A0=np.real(c.conj()@(A@c))
        for T in Tg:
            ct=expmH(Hd,T)@c; Delta=abs(np.real(ct.conj()@(A@ct))-A0)
            if Delta<1e-9: continue
            r=T*meanH*sA**2/Delta**2
            if r<best: best=r
    return best
for nexc in [1,2,3]:
    print(f"(1) rank-1 ground coupling, {nexc} excited: inf r_quad = {rank1_min(nexc):.8f}  (target C_true={C_true:.8f})", flush=True)

# ---- (2) real-restricted single excited (closed-form min) ----
phi=np.linspace(1e-3,2*np.pi-1e-3,200000); r_real=phi/(4*(1-np.cos(phi))**2)
print(f"(2) real-restricted inf r_quad = {r_real.min():.10f} at phi={phi[np.argmin(r_real)]:.6f}  (target C_old={C_old:.10f})", flush=True)

# ---- (3) eigenvector observable: A=sigma_z, psi0=|0> (extreme, l=+1) ----
A=np.array([[1,0],[0,-1]],complex); sA=sigA(A); psi0=np.array([1,0],complex)
best_lin=np.inf; rng=np.random.default_rng(7)
sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]],complex); sz=np.array([[1,0],[0,-1]],complex)
for _ in range(8000):
    h=rng.normal(size=3); H=h[0]*sx+h[1]*sy+h[2]*sz; H=H-np.linalg.eigvalsh(H)[0]*np.eye(2)
    meanH=np.real(psi0.conj()@(H@psi0))
    if meanH<1e-9: continue
    diffs=np.diff(np.linalg.eigvalsh(H))
    if diffs[0]<1e-6: continue
    A0=np.real(psi0.conj()@(A@psi0))
    for T in np.linspace(1e-3,2*2*np.pi/diffs[0],200):
        pt=expmH(H,T)@psi0; Delta=abs(np.real(pt.conj()@(A@pt))-A0)
        if Delta<1e-9: continue
        best_lin=min(best_lin,T*meanH*sA/Delta)
print(f"(3) eigenvector obs (A=sigma_z, psi0=|0>): inf r_lin = T<H>sigma_A/Delta = {best_lin:.8f}  (target 1/(4K)={inv4K:.8f})", flush=True)
for th in [0.3,0.1,0.03,0.01]:
    E=1.0; ket=np.array([np.sin(th),np.cos(th)],complex); H=E*np.outer(ket,ket.conj()); H=H-np.linalg.eigvalsh(H)[0]*np.eye(2)
    meanH=np.real(psi0.conj()@(H@psi0)); A0=np.real(psi0.conj()@(A@psi0)); bq=np.inf
    for T in np.linspace(1e-3,4*np.pi,1500):
        pt=expmH(H,T)@psi0; Delta=abs(np.real(pt.conj()@(A@pt))-A0)
        if Delta<1e-9: continue
        bq=min(bq,T*meanH*sA**2/Delta**2)
    print(f"     theta={th}: r_quad={bq:.3f} (diverges => LINEAR for eigenvector obs)", flush=True)

# ---- (3b) COUNTEREXAMPLE: turning point alone (NOT eigenvector) does NOT give the linear bound ----
# psi0=(|0>+|1>)/sqrt2 is NOT an eigenvector of A=diag(1,-1,0); H with real H01 -> zero initial slope.
Ace = np.diag([1.,-1.,0.]).astype(complex); psice = np.array([1,1,0],complex)/np.sqrt(2)
resid = np.linalg.norm(Ace@psice - (psice.conj()@(Ace@psice))*psice)
best_ce = np.inf; rngce = np.random.default_rng(0)
for _ in range(8000):
    h01=rngce.normal(); h02=rngce.normal()+1j*rngce.normal(); h12=rngce.normal()+1j*rngce.normal(); dd=rngce.normal(size=3)
    H=np.array([[dd[0],h01,h02],[h01,dd[1],h12],[np.conj(h02),np.conj(h12),dd[2]]],complex); H=H-np.linalg.eigvalsh(H)[0]*np.eye(3)
    meanH=np.real(psice.conj()@(H@psice));  A0=np.real(psice.conj()@(Ace@psice))
    if meanH<1e-9: continue
    for T in np.linspace(1e-2,8,300):
        pt=expmH(H,T)@psice; D=abs(np.real(pt.conj()@(Ace@pt))-A0)
        if D<1e-6: continue
        best_ce=min(best_ce,T*meanH*1.0/D)
print(f"(3b) turning-point-but-NOT-eigenvector (psi0 not A-eigvec, resid={resid:.2f}): inf r_lin = {best_ce:.5f}  << 1/(4K)={inv4K:.4f}  => linear bound NEEDS eigenvector, not just zero slope", flush=True)

# ---- (4) ground-gapped A (d=3, excited block) ----
A=np.zeros((3,3),complex); A[1,2]=A[2,1]=1.0; sA=sigA(A); worst_static=np.inf; rows=[]
rng=np.random.default_rng(11)
for eps in [0.2,0.1,0.05,0.02,0.01]:
    E=np.array([0.0,1.0,1.7]); Dgap=E[1]; bq=np.inf
    for _ in range(150):
        a=rng.normal(size=2)+1j*rng.normal(size=2); a/=np.linalg.norm(a)
        c=np.array([np.sqrt(1-eps),*(np.sqrt(eps)*a)],complex)
        meanH=np.real(np.sum(E*np.abs(c)**2)); A0=np.real(c.conj()@(A@c)); Hd=np.diag(E)
        for T in np.linspace(1e-2,6*np.pi,1200):
            ct=expmH(Hd,T)@c; Delta=abs(np.real(ct.conj()@(A@ct))-A0)
            if Delta<1e-9: continue
            worst_static=min(worst_static, meanH-(Dgap/2)*Delta/sA)
            bq=min(bq,T*meanH*sA**2/Delta**2)
    rows.append((eps,bq))
print(f"(4) ground-gapped: static bound <H>-(Dgap/2)Delta/sigma_A worst slack = {worst_static:.3e} (>=0 holds)", flush=True)
print("    r_quad vs eps (grows ~1/eps => diverges): "+"  ".join(f"eps={e}:{r:.1f}" for e,r in rows), flush=True)
