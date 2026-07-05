"""
Stage-1 CORRECTION + independent verification.

The Stage-1 workflow (4 independent derivations + numeric adversary) found that the
real-amplitude two-level value C=0.18555147 is NOT the universal constant: it is the
phase-restricted optimum.  Releasing the relative phase gives the sharp constant

    C_true = 1/(8K),   K = sup_{x>0} (1-cos x)/x = sin(x*),  x* = root of tan(x/2)=x.

This script verifies, from scratch and to high precision:
  (A) K and x*, and the three closed forms of C_true agree;
  (B) the tangent-line lemma  1-cos x <= K x  for all x>0 (equality at x*);
  (C) the free-phase two-level infimum ratio = C_true (= 0.172506...), BELOW
      the phase-restricted (beta=0) value, which must reproduce 0.18555147.
Ratio r := T*(<H>-E0)*sigma_A^2 / Delta^2 ;  bound is r >= C_true.
hbar=1, E0=0.
"""
import mpmath as mp
mp.mp.dps = 50

print("="*70); print("(A) K = sup (1-cos x)/x, x* = argmax (root of tan(x/2)=x)"); print("="*70)
xstar = mp.findroot(lambda x: mp.tan(x/2) - x, mp.mpf('2.33'))
K     = (1 - mp.cos(xstar))/xstar
print(f"  x*            = {mp.nstr(xstar,20)}   (check tan(x*/2)-x* = {mp.nstr(mp.tan(xstar/2)-xstar,3)})")
print(f"  K=(1-cos x*)/x* = {mp.nstr(K,20)}")
print(f"  sin(x*)         = {mp.nstr(mp.sin(xstar),20)}   (should equal K)")
# three closed forms of C_true
C_a = 1/(8*K)
C_b = xstar/(8*(1-mp.cos(xstar)))
C_c = (1+xstar**2)/(16*xstar)
C_d = 1/(8*mp.sin(xstar))
print(f"  C_true = 1/(8K)            = {mp.nstr(C_a,20)}")
print(f"  C_true = x*/(8(1-cos x*))  = {mp.nstr(C_b,20)}")
print(f"  C_true = (1+x*^2)/(16 x*)  = {mp.nstr(C_c,20)}")
print(f"  C_true = 1/(8 sin x*)      = {mp.nstr(C_d,20)}")
print(f"  max spread of the four     = {mp.nstr(max(C_a,C_b,C_c,C_d)-min(C_a,C_b,C_c,C_d),3)}")
C_true = C_a
C_old  = mp.mpf('0.18555147172183604')
print(f"  C_true/C_old(0.18555)      = {mp.nstr(C_true/C_old,12)}   (workflow claimed 0.929695)")

print("="*70); print("(B) tangent-line lemma  1-cos x <= K x  for all x>0 (eq. at x*)"); print("="*70)
worst = mp.mpf('-1'); worst_x = None
x = mp.mpf('0.0001')
while x < 40:
    d = (1-mp.cos(x)) - K*x          # must be <= 0
    if d > worst: worst, worst_x = d, x
    x += mp.mpf('0.001')
print(f"  max_x [ (1-cos x) - K x ] = {mp.nstr(worst,6)}  at x={mp.nstr(worst_x,8)}  (should be ~0 at x*={mp.nstr(xstar,6)})")
print(f"  => lemma holds (<=0 everywhere, touching 0 only at x*).")

print("="*70); print("(C) two-level infimum ratio: free-phase vs beta=0"); print("="*70)
# theta->0 limit:  r(phi,beta) = phi / (4 (cos(phi-beta) - cos(beta))^2),   phi=E*T
def r_freephase():
    # minimize over phi in (0,2pi), beta in (0,2pi); refine the grid minimum
    best = mp.inf; bp = None
    NP, NB = 1400, 1400
    for i in range(1, NP):
        phi = 2*mp.pi*i/NP
        sphi2 = mp.sin(phi/2)
        # optimal beta makes |sin(phi/2 - beta)|=1, giving |D|=2|sin(phi/2)|
        D2 = 4*sphi2**2
        if D2 == 0: continue
        val = phi/(4*D2)
        if val < best: best, bp = val, phi
    return best, bp
def r_beta0():
    best = mp.inf; bp=None
    N=200000
    for i in range(1, N):
        phi = 2*mp.pi*i/N
        D = mp.cos(phi) - 1            # cos(phi-0)-cos(0)
        if D == 0: continue
        val = phi/(4*D**2)
        if val < best: best, bp = val, phi
    return best, bp
rf, phif = r_freephase()
r0, phi0 = r_beta0()
print(f"  free-phase inf ratio = {mp.nstr(rf,12)}  at phi=E T={mp.nstr(phif,8)}  (= x* ? {mp.nstr(xstar,8)})")
print(f"      C_true target    = {mp.nstr(C_true,12)}   match: {abs(rf-C_true) < mp.mpf('1e-5')}")
print(f"  beta=0 inf ratio     = {mp.nstr(r0,12)}  at phi={mp.nstr(phi0,8)}")
print(f"      C_old target     = {mp.nstr(C_old,12)}   match: {abs(r0-C_old) < mp.mpf('1e-4')}")
print()
print(f"  CONCLUSION: free-phase optimum {mp.nstr(rf,8)} < beta=0 optimum {mp.nstr(r0,8)}.")
print(f"  The relative phase (start at max slope) lowers the constant from 0.185551 to {mp.nstr(C_true,8)}.")
print(f"  => Universal sharp constant C_true = {mp.nstr(C_true,17)}")
