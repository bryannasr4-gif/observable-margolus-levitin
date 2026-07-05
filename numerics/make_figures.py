"""
P3.3 — Figures for the paper -> paper/figs/*.pdf  (matplotlib Agg, seeded).

fig1_convergence.pdf : 2-level ratio T<H>sigma_A^2/Delta^2 vs theta, free-phase -> C_true,
                       phase-restricted (beta=0) -> C_old, both from above.
fig2_tangent.pdf     : g(x)=(1-cos x)/x with maximum K=sin x* at x* (the ML tangent-line lemma).
fig3_trajectory.pdf  : saturating-state <A(t)> trajectory (free-phase: starts at max slope;
                       beta=0: starts at a turning point) with Delta and T=x* marked.
fig4_multid.pdf      : per-dimension min ratio (random + near-ground search) sitting above C_true.

Run: .venv/bin/python -u numerics/make_figures.py
"""
import numpy as np, mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mp.mp.dps = 30
xstar = float(mp.findroot(lambda x: mp.tan(x/2) - x, mp.mpf('2.33')))
K = (1 - np.cos(xstar)) / xstar
C_true = 1/(8*K); C_old = 0.18555147172183604
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 140})

def expmH(H, t):
    w, V = np.linalg.eigh(H); return (V*np.exp(-1j*w*t)) @ V.conj().T
sx = np.array([[0,1],[1,0]], complex); H2 = np.array([[0,0],[0,1.0]])

def twolevel_ratio(th, freephase):
    betas = np.linspace(0, 2*np.pi, 160) if freephase else np.array([0.0])
    best = np.inf
    for b in betas:
        c = np.array([np.cos(th), np.sin(th)*np.exp(1j*b)])
        A0 = np.real(c.conj()@(sx@c)); meanH = np.sin(th)**2
        for T in np.linspace(1e-3, 2*np.pi, 1400):
            ct = expmH(H2, T)@c; D = abs(np.real(ct.conj()@(sx@ct))-A0)
            if D < 1e-9: continue
            best = min(best, T*meanH/D**2)
    return best

# ---- fig 1: convergence ----
thetas = np.linspace(0.6, 0.02, 26)
rf = [twolevel_ratio(t, True) for t in thetas]
r0 = [twolevel_ratio(t, False) for t in thetas]
fig, ax = plt.subplots(figsize=(5.2, 3.8))
ax.plot(thetas, rf, 'o-', ms=4, label=r'free phase $\to C_{\rm true}$', color='C0')
ax.plot(thetas, r0, 's-', ms=4, label=r'$\beta=0$ (real) $\to C_{\rm old}$', color='C3')
ax.axhline(C_true, ls='--', color='C0', lw=1); ax.axhline(C_old, ls='--', color='C3', lw=1)
ax.text(0.45, C_true-0.004, rf'$C_{{\rm true}}={C_true:.5f}$', color='C0', fontsize=9)
ax.text(0.45, C_old+0.002, rf'$C_{{\rm old}}={C_old:.5f}$', color='C3', fontsize=9)
ax.set_xlabel(r'ground-dominance angle $\theta$'); ax.set_ylabel(r'$T\langle H\rangle\,\sigma_A^2/\Delta^2$')
ax.set_title('Two-level ratio approaches the sharp constant'); ax.legend(loc='upper left', fontsize=9)
ax.set_ylim(0.16, 0.30); fig.tight_layout(); fig.savefig("paper/figs/fig1_convergence.pdf"); plt.close(fig)

# ---- fig 2: tangent-line lemma ----
x = np.linspace(0.01, 8, 1000); g = (1-np.cos(x))/x
fig, ax = plt.subplots(figsize=(5.2, 3.8))
ax.plot(x, g, color='C0', lw=1.8, label=r'$(1-\cos x)/x$')
ax.axhline(K, ls='--', color='k', lw=1, label=rf'$K=\sin x^*={K:.5f}$')
ax.plot([xstar],[K],'o',color='C3',ms=7,zorder=5)
ax.axvline(xstar, ls=':', color='C3', lw=1)
ax.annotate(rf'$x^*={xstar:.4f}$', xy=(xstar, K), xytext=(xstar+1.0, K-0.13),
            arrowprops=dict(arrowstyle='->', color='C3'), color='C3', fontsize=9)
ax.set_xlabel(r'$x=E\,T$'); ax.set_ylabel(r'$(1-\cos x)/x$')
ax.set_title(r'ML tangent-line lemma: $\max_x (1-\cos x)/x = K$'); ax.legend(loc='upper right', fontsize=9)
fig.tight_layout(); fig.savefig("paper/figs/fig2_tangent.pdf"); plt.close(fig)

# ---- fig 3: saturating trajectory ----
th = 0.35; tt = np.linspace(0, 2*xstar, 600)
def traj(b):
    c = np.array([np.cos(th), np.sin(th)*np.exp(1j*b)]); A0 = np.real(c.conj()@(sx@c))
    return np.array([np.real((expmH(H2,t)@c).conj()@(sx@(expmH(H2,t)@c))) for t in tt]), A0
b_free = xstar/2 - np.pi/2
yf, A0f = traj(b_free); y0, A00 = traj(0.0)
fig, ax = plt.subplots(figsize=(5.6, 3.8))
ax.plot(tt, yf, color='C0', lw=1.8, label=r'free phase $\beta=x^*/2-\pi/2$ (max slope at $t=0$)')
ax.plot(tt, y0, color='C3', lw=1.8, label=r'$\beta=0$ (turning point at $t=0$)')
ax.axvline(xstar, ls=':', color='gray'); ax.text(xstar+0.05, ax.get_ylim()[0]+0.05, r'$T=x^*$', fontsize=9)
# mark Delta for free-phase at first max
iT = np.argmin(np.abs(tt-xstar))
ax.annotate('', xy=(xstar, yf[iT]), xytext=(xstar, A0f),
            arrowprops=dict(arrowstyle='<->', color='C0'))
ax.text(xstar-1.15, (yf[iT]+A0f)/2, r'$\Delta$', color='C0', fontsize=11)
ax.set_xlabel(r'$t$'); ax.set_ylabel(r'$\langle A(t)\rangle$ ($A=\sigma_x$)')
ax.set_title(rf'Saturating trajectory ($\theta={th}$)'); ax.legend(loc='lower right', fontsize=8)
fig.tight_layout(); fig.savefig("paper/figs/fig3_trajectory.pdf"); plt.close(fig)

# ---- fig 4: multi-d min ratio above C_true ----
def sigA(A):
    w = np.linalg.eigvalsh(A); return (w[-1]-w[0])/2
def mind(d, seed, ntri=900):
    rng = np.random.default_rng(seed); best = np.inf
    # (a) random near-ground configs (shows nothing dips below C_true)
    for _ in range(ntri):
        E = np.zeros(d); E[1:] = rng.uniform(0.3, 3.0, d-1)
        M = rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)); A=(M+M.conj().T)/2; A/=sigA(A)
        th = rng.uniform(0.01, 0.5)
        a = rng.normal(size=d-1)+1j*rng.normal(size=d-1); a/=np.linalg.norm(a)
        c = np.zeros(d, complex); c[0]=np.cos(th); c[1:]=np.sin(th)*a
        meanH = np.real(np.sum(E*np.abs(c)**2))
        if meanH < 1e-9: continue
        diffs = np.abs(E[:,None]-E[None,:]); wpos=diffs[diffs>1e-6]
        if wpos.size==0: continue
        Hd = np.diag(E); A0 = np.real(c.conj()@(A@c))
        for T in np.linspace(1e-3, 2.5*2*np.pi/wpos.min(), 200):
            ct = expmH(Hd,T)@c; D = abs(np.real(ct.conj()@(A@ct))-A0)
            if D < 1e-9: continue
            best = min(best, T*meanH*sigA(A)**2/D**2)
    # (b) structured 2-level near-ground saturator embedded in dimension d
    #     ground |0> coupled to one excited |1> by sigma_x (spectators 2..d-1 inert);
    #     optimal phase beta=x*/2-pi/2, E*T=x*, theta->0 -> ratio -> C_true for EVERY d.
    A = np.zeros((d,d), complex); A[0,1]=A[1,0]=1.0
    for th in [0.02, 0.01, 0.005]:
        c = np.zeros(d, complex); c[0]=np.cos(th); c[1]=np.sin(th)*np.exp(1j*(xstar/2-np.pi/2))
        E = np.zeros(d); E[1]=1.0; E[2:] = np.linspace(1.5, 3.0, max(0,d-2)); Hd=np.diag(E)
        meanH = np.real(np.sum(E*np.abs(c)**2)); A0 = np.real(c.conj()@(A@c))
        for T in np.linspace(1e-3, 2*np.pi, 1400):
            ct = expmH(Hd,T)@c; D = abs(np.real(ct.conj()@(A@ct))-A0)
            if D < 1e-9: continue
            best = min(best, T*meanH*sigA(A)**2/D**2)
    return best
ds = list(range(2, 8)); mins = [mind(d, 1000+d) for d in ds]
fig, ax = plt.subplots(figsize=(5.2, 3.8))
ax.axhline(C_true, ls='--', color='C0', lw=1.3, label=rf'$C_{{\rm true}}={C_true:.6f}$')
ax.scatter(ds, mins, color='C3', zorder=5, label='min ratio found (random search)')
for d, m in zip(ds, mins): ax.annotate(f'{m:.4f}', (d, m), textcoords='offset points', xytext=(0,7), fontsize=8, ha='center')
ax.set_xlabel('Hilbert-space dimension $d$'); ax.set_ylabel(r'min $T\langle H\rangle\,\sigma_A^2/\Delta^2$')
ax.set_title('No configuration beats $C_{\\rm true}$'); ax.legend(loc='upper left', fontsize=9)
ax.set_ylim(C_true-0.01, max(mins)+0.03); fig.tight_layout(); fig.savefig("paper/figs/fig4_multid.pdf"); plt.close(fig)

print("min ratios per d:", {d: round(m,6) for d,m in zip(ds,mins)})
print("WROTE paper/figs/fig1_convergence.pdf, fig2_tangent.pdf, fig3_trajectory.pdf, fig4_multid.pdf")
