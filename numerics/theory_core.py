"""
Stage 1 analytic core for:
  A Margolus-Levitin (mean-energy) speed limit for observables.

  !!! CORRECTION (Stage 1 result): the constant 0.18555147 derived below is the
  PHASE-RESTRICTED optimum (real amplitudes, A=sigma_x, zero relative phase ->
  start at a turning point). It is NOT the universal constant. The sharp universal
  constant is  C_true = 1/(8 sin x*) = 0.17250626746116262  (x* = root of
  tan(x/2)=x). See numerics/theory_core_v2.py and notes/thm_constant.md.
  The NO-GO content below (P/Delta -> 0) is unaffected.

Establishes rigorously, for the exactly-solvable TWO-LEVEL family
(H=diag(0,E), A=sigma_x, psi=cos(th)|0>+sin(th)|1>):
  (1) <A(t)> = sin(2 th) cos(E t)   [symbolic],
  (2) P(th;Delta) := T*(<H>-E0) = sin^2(th) * arccos(1 - Delta/sin(2th)),
  (3) NO-GO: P_min(Delta)/Delta -> 0  (no linear mean-energy bound),
  (4) the optimal constant C and its CLOSED-FORM characterization:
        C = min_{s>=1/2} s^2 arccos(1-1/s) / 4,
      with optimum s* solving the transcendental equation
        2 arccos(1-1/s*) * sqrt(2 s*-1) = 1,
      and the closed form  C = s*^2 / (8 sqrt(2 s*-1)).
hbar = 1, E0 = 0.
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 40

# ---------------------------------------------------------------- PART 1
print("="*70); print("PART 1  symbolic <A(t)> for the two-level family"); print("="*70)
t, E, th = sp.symbols('t E theta', real=True, positive=True)
c0, c1 = sp.cos(th), sp.sin(th)
# psi(t) = c0|0> + c1 e^{-iEt}|1>;  A=sigma_x => <A>=2 Re(c0 c1 e^{-iEt})
At = sp.simplify(2*c0*c1*sp.cos(E*t))
print("  <A(t)> =", At, "   (= sin(2 theta) cos(E t))")
meanH = sp.simplify(c1**2 * E)               # <H> - E0,  E0 = 0
print("  <H>-E0 =", meanH, "   (= E sin^2 theta)")

# ---------------------------------------------------------------- PART 2
print("="*70); print("PART 2  g(s), its derivative, and the transcendental optimum"); print("="*70)
s = sp.symbols('s', positive=True)
g = s**2 * sp.acos(1 - 1/s)
gp = sp.simplify(sp.diff(g, s))
print("  g(s)   =", g)
print("  g'(s)  =", gp)
# Show g'(s)=0  <=>  2 acos(1-1/s) = 1/sqrt(2s-1)
lhs_minus_rhs = sp.simplify(gp/s - 0)        # divide out s>0
print("  g'(s)/s =", sp.simplify(gp/s), "  (set = 0)")

# ---------------------------------------------------------------- PART 3
print("="*70); print("PART 3  high-precision constant (mpmath)"); print("="*70)
gm   = lambda x: x**2 * mp.acos(1 - 1/x)
trans = lambda x: 2*mp.acos(1 - 1/x)*mp.sqrt(2*x - 1) - 1     # =0 at optimum
sstar = mp.findroot(trans, mp.mpf('0.516'))
C_min   = gm(sstar)/4
C_closed = sstar**2/(8*mp.sqrt(2*sstar - 1))
# independent brute-force minimization as a cross-check
xs = [mp.mpf('0.5') + mp.mpf(i)/200000 for i in range(0, 200000)]
gbrute = min(gm(x) for x in xs)
print(f"  s*                      = {mp.nstr(sstar, 20)}")
print(f"  transcendental residual = {mp.nstr(trans(sstar), 5)}  (should be ~0)")
print(f"  C = g(s*)/4             = {mp.nstr(C_min, 20)}")
print(f"  C = s*^2/(8 sqrt(2s*-1))= {mp.nstr(C_closed, 20)}   (closed form)")
print(f"  C  (brute min g/4)      = {mp.nstr(gbrute/4, 20)}   (independent check)")
print(f"  agreement |min - closed|= {mp.nstr(abs(C_min - C_closed), 5)}")
print(f"  matches explore script 0.185551 ? -> {mp.nstr(C_min,6)}")

# ---------------------------------------------------------------- PART 4
print("="*70); print("PART 4  NO-GO: P_min(Delta)/Delta -> 0  (no linear bound)"); print("="*70)
def Pmin(Delta):
    Delta = mp.mpf(Delta)
    # minimize sin^2(th) arccos(1 - Delta/sin(2th)) over feasible th (sin2th>=Delta/2)
    def P(th):
        s2 = mp.sin(2*th)
        if s2 < Delta/2:
            return mp.inf
        x = 1 - Delta/s2
        if x < -1: x = mp.mpf(-1)
        return mp.sin(th)**2 * mp.acos(x)
    best = mp.inf
    N = 4000
    for i in range(1, N):
        th = mp.pi/4 * i/N
        v = P(th)
        if v < best: best = v
    return best
print(f"  {'Delta':>10}{'P_min':>16}{'P_min/Delta':>16}{'P_min/Delta^2':>16}")
for D in ['0.2','0.1','0.05','0.02','0.01','0.005']:
    P = Pmin(D); Df = mp.mpf(D)
    print(f"  {D:>10}{mp.nstr(P,8):>16}{mp.nstr(P/Df,8):>16}{mp.nstr(P/Df**2,8):>16}")
print("  -> P_min/Delta -> 0 : NO linear mean-energy observable bound (no-go).")
print(f"  -> P_min/Delta^2 -> C = {mp.nstr(C_min,8)} : quadratic law with the closed-form constant.")
