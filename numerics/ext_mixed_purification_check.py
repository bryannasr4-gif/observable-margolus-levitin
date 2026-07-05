"""
P2.2 — Purification / Uhlmann route: persistent numeric cross-check
(replaces the workflow's scratch/ scripts so notes/ext_mixed.md is reproducible).

Verifies, for random mixed states (hbar=1, E0=0, rho_T=U rho U^dag, U=exp(-iHT)):
  (A) Purification preservation: <H(x)I>_Psi = <H>, sigma_{A(x)I} = sigma_A, Delta' = Delta.
  (B) LEMMA U1: Uhlmann fidelity F(rho, U rho U^dag) >= |Tr[rho U]|.
  (C) No-factor-2 deficit: 1 - F <= K T <H>  (the pure-state ML inequality, mixed).
  (D) Full Uhlmann chain ratio  T<H> sigma_A^2 / Delta^2 >= C_true.
"""
import numpy as np, mpmath as mp
from scipy.linalg import sqrtm
mp.mp.dps = 30
xstar = mp.findroot(lambda x: mp.tan(x/2) - x, mp.mpf('2.33'))
K = float((1 - mp.cos(xstar)) / xstar); C_true = 1/(8*K)
print(f"K={K:.12f}  C_true={C_true:.12f}\n")

def expmH(H,t):
    w,V=np.linalg.eigh(H); return (V*np.exp(-1j*w*t))@V.conj().T
def sigA(A):
    w=np.linalg.eigvalsh(A); return (w[-1]-w[0])/2
def uhlmann_F(rho,sig):
    s=sqrtm(rho); M=s@sig@s
    w=np.linalg.eigvalsh((M+M.conj().T)/2); w=np.clip(w,0,None)
    return float(np.sum(np.sqrt(w)))

rng=np.random.default_rng(2024)
wa=wb=wc=wd=np.inf; minratio=np.inf
NT=4000
for _ in range(NT):
    d=rng.integers(2,5)
    E=np.abs(rng.normal(scale=rng.uniform(0.5,2.5),size=d)); E[0]=0; E-=E.min(); H=np.diag(E)
    M=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)); A=(M+M.conj().T)/2; A/=sigA(A); sAv=1.0
    p=rng.dirichlet(np.ones(d)*rng.uniform(0.3,2.0))
    Z=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)); Q,_=np.linalg.qr(Z)
    rho=(Q*p)@Q.conj().T
    diffs=np.abs(E[:,None]-E[None,:]); wpos=diffs[diffs>1e-6]
    if wpos.size==0: continue
    T=rng.uniform(0.05,1.5)*2*np.pi/wpos.min()
    U=expmH(H,T); rhoT=U@rho@U.conj().T
    meanH=np.real(np.trace(rho@H)); Delta=abs(np.real(np.trace(A@(rhoT-rho))))
    # (A) purification preservation
    vals,vecs=np.linalg.eigh(rho)
    Psi=np.zeros(d*d,complex)
    for n in range(d):
        Psi += np.sqrt(max(vals[n],0))*np.kron(vecs[:,n], np.eye(d)[:,n])
    Hp=np.kron(H,np.eye(d)); Ap=np.kron(A,np.eye(d)); Up=np.kron(U,np.eye(d))
    meanHp=np.real(Psi.conj()@(Hp@Psi)); sAp=sigA(Ap)
    PsiT=Up@Psi; A0p=np.real(Psi.conj()@(Ap@Psi)); DeltaP=abs(np.real(PsiT.conj()@(Ap@PsiT))-A0p)
    wa=min(wa, -max(abs(meanHp-meanH),abs(sAp-sAv),abs(DeltaP-Delta)))  # want ~0; store as nonneg slack
    # (B) LEMMA U1
    F=uhlmann_F(rho,rhoT); trU=abs(np.trace(rho@U))
    wb=min(wb, F-trU)
    # (C) deficit 1-F <= K T <H>
    wc=min(wc, K*T*meanH-(1-F))
    # (D) ratio
    if Delta>1e-7 and meanH>1e-9:
        ratio=T*meanH*sAv**2/Delta**2; minratio=min(minratio,ratio)
        wd=min(wd, 8*K*sAv**2*T*meanH - Delta**2)
print(f"(A) purification preservation: worst |mismatch| = {-wa:.2e}  (==0 => exact)")
print(f"(B) LEMMA U1  F - |Tr[rho U]|: worst = {wb:.2e}  (>=0 => holds)")
print(f"(C) deficit  K T<H> - (1-F):  worst = {wc:.2e}  (>=0 => no factor 2)")
print(f"(D) chain 8K sigma^2 T<H> - Delta^2: worst = {wd:.2e}  (>=0 => holds)")
print(f"    min Uhlmann-route ratio T<H>sigma^2/Delta^2 = {minratio:.8f}  (>= C_true={C_true:.8f}? {minratio>=C_true-1e-7})")
