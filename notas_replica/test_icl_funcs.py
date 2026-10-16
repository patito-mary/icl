# %% Imports y configuración
"""Pruebas sintéticas de icl_funcs.py (réplica Mayes+2026). Uso: python test_icl_funcs.py ../mayes_2026_methodology"""
import sys, types
import numpy as np
from scipy.integrate import quad
from scipy.special import gammaincinv, gammainc, gamma
sys.path.insert(0, sys.argv[1])
import icl_funcs as F
import params_icl as P
rng = np.random.default_rng(42)
OK = []

def check(name, cond, info=''):
    OK.append(bool(cond))
    print(f"[{'OK ' if cond else 'FALLA'}] {name}  {info}")

# %% 1. Área del anillo gaussiano: fórmula exacta vs integral numérica y límite (2π)^1.5 r0 σ
for r0, s in [(0.1, 1.0), (1.0, 1.0), (10.0, 1.5), (300.0, 45.0)]:
    num = quad(lambda R: 2*np.pi*R*np.exp(-(R-r0)**2/(2*s**2)), 0, r0 + 40*s)[0]
    ex = F.gaussian_ring_area_kpc2(r0, s)
    check(f"área kernel r0={r0} σ={s}", abs(ex/num - 1) < 1e-8, f"exacta/num={ex/num:.10f}  lim/num={(2*np.pi)**1.5*r0*s/num:.4f}")

# %% 2. Σ constante (disco uniforme R=600 kpc): anillos y kernel deben dar Δμ ≈ 0
N = 2_000_000; Rd = 600.0; Ltot = 1e12
R = Rd*np.sqrt(rng.random(N)); lum = np.full(N, Ltot/N)
mu_true = F.mu_from_sigma(Ltot/(np.pi*(Rd*1e3)**2), P.SB_CONST)
for meth in ['anillos', 'kernel']:
    for binning in ['lineal', 'log']:
        rm, mu = F.sb_profile(R, lum, 500.0, P.SB_CONST, method=meth, binning=binning)
        sel = (rm > 20) & (rm < 450)
        d = np.median(mu[sel] - mu_true)
        check(f"Σ cte {meth}/{binning}", abs(d) < 0.01, f"Δμ mediana={d:+.4f} mag")

# %% 3. Σ ∝ R^-α (μ exactamente lineal en log R): el ajuste debe recuperar a, b y R_H
alpha = 1.6; Rin, Rout = 1.0, 1000.0
u = rng.random(N)   # p(R) ∝ R^{1-α} entre Rin y Rout
R = (Rin**(2-alpha) + u*(Rout**(2-alpha) - Rin**(2-alpha)))**(1/(2-alpha))
Ltot = 3e12; lum = np.full(N, Ltot/N)
norm = Ltot/(2*np.pi*(Rout**(2-alpha)-Rin**(2-alpha))/(2-alpha))   # Σ(R)=norm R^-α [L☉/kpc²]
mu_R = lambda r: F.mu_from_sigma(norm*r**-alpha/1e6, P.SB_CONST)
b_true = 2.5*alpha; a_true = mu_R(1.0)
rh_true = 10**((P.MU_HOLMBERG - a_true)/b_true)
rm, mu = F.sb_profile(R, lum, 900.0, P.SB_CONST, method='anillos', binning='lineal', bin_kpc=5.0)
hol = F.holmberg_radius(rm, mu, P.MU_HOLMBERG, method='ajuste', r_fit_min=5.0, r_fit_max=900.0)
check("ajuste log: pendiente", abs(hol['b']/b_true - 1) < 0.02, f"b={hol['b']:.3f} (verdad {b_true:.3f})")
check("ajuste log: R_H", abs(hol['R_H']/rh_true - 1) < 0.03, f"R_H={hol['R_H']:.1f} kpc (verdad {rh_true:.1f}); interp={hol['R_H_interp']:.1f}")

# %% 4. Sérsic n=4 (BCG): R_H del ajuste log vs interpolación vs verdad, y f_ICL
n, Re = 4.0, 40.0; bn = gammaincinv(2*n, 0.5)
R = Re*(gammaincinv(2*n, rng.random(N))/bn)**n
for Ltot in [3e11, 1e12, 3e12]:
    lum = np.full(N, Ltot/N); mass = np.ones(N)
    Ie = Ltot/(2*np.pi*n*(Re*1e3)**2*np.exp(bn)*bn**(-2*n)*gamma(2*n))
    Scut = 10**(-0.4*(P.MU_HOLMBERG - P.SB_CONST))
    rh_true = Re*((np.log(Ie/Scut) + bn)/bn)**n
    fi_true = 1 - gammainc(2*n, bn*(rh_true/Re)**(1/n))
    rmax = min(np.percentile(R, 99.5), 1500)
    for meth in ['anillos', 'kernel']:
        rm, mu = F.sb_profile(R, lum, rmax, P.SB_CONST, method=meth, binning='lineal')
        hol = F.holmberg_radius(rm, mu, P.MU_HOLMBERG, r_fit_max=rmax)
        mb, mi = F.split_bcg_icl(R, hol['R_H_interp'])
        fr = F.icl_fractions(mass, mb, mi, mass.sum())
        print(f"  Sérsic L={Ltot:.0e} {meth:7s}: verdad R_H={rh_true:6.1f} f={fi_true:.3f} | "
              f"interp R_H={hol['R_H_interp']:6.1f} f={fr['f_icl']:.3f} | ajuste log R_H={hol['R_H_fit']:6.1f}")
    check(f"Sérsic interp R_H (L={Ltot:.0e})", abs(hol['R_H_interp']/rh_true-1) < 0.05)
    hl = F.holmberg_radius(rm, mu, P.MU_HOLMBERG, method='ajuste_local', r_fit_max=rmax, mu_window=1.0)
    check(f"Sérsic ajuste_local R_H (L={Ltot:.0e})", abs(hl['R_H']/rh_true-1) < 0.05,
          f"local={hl['R_H']:.1f}  global={hl['R_H_fit']:.1f}  verdad={rh_true:.1f}")

# %% 5. f_ICL con denominador M*200 vs BCG+ICL
mass = np.array([1., 1., 1., 1.]); mb = np.array([1, 1, 0, 0], bool); mi = ~mb
fr = F.icl_fractions(mass, mb, mi, mstar200=8.0)
check("f_icl = M_ICL/M*200", np.isclose(fr['f_icl'], 0.25) and np.isclose(fr['f_icl_bcgicl'], 0.5), str(fr))

# %% 6. Rotación por tensor de inercia: recupera ejes de un elipsoide rotado
pts = rng.normal(size=(200000, 3))*np.array([30., 15., 5.])
th = np.radians(35); Rz = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
pr, Rm = F.rotate_by_inertia_tensor(pts @ Rz.T, np.ones(len(pts)), reduced=False)
sd = pr.std(axis=0)
check("rotación: ejes ordenados", sd[0] > sd[1] > sd[2] and abs(sd[0]/30-1) < 0.02, f"σ={sd.round(2)} det={np.linalg.det(Rm):+.0f}")

# %% 7. Cuña de 15° sobre el eje mayor
p = np.array([[1, 0, 0], [-1, 0.1, 0], [1, 0.2, 0], [0, 1, 0]], float)
check("cuña ±7.5°", list(F.wedge_mask(p, 15.0)) == [True, True, False, False])

# %% 8. Estado dinámico
st = F.dynamical_state([0.05, 0.1, 0.25, 0.4, 0.5, np.nan])
check("dyn_state", list(st) == [0, 1, 1, 1, 2, -1], str(st))

# %% 9. NFW: recuperar c de un halo sintético c=6
c_true, r200 = 6.0, 1000.0; rs = r200/c_true
mu_ = lambda x: np.log(1+x) - x/(1+x)
x = np.linspace(1e-4, c_true, 200000); cdf = mu_(x)/mu_(c_true)
r_dm = rs*np.interp(rng.random(1_000_000), cdf, x)
c_fit, rs_fit, _, _ = F.fit_nfw_concentration(r_dm, 1.0, r200)
check("NFW c", abs(c_fit/c_true - 1) < 0.05, f"c_fit={c_fit:.2f}")

# %% 10. Último major merger en un árbol SubLink sintético
# MPB: 0(s99) → 1(s98) → 2(s97) → 3(s96). Progenitores secundarios:
#   nodo 4 (s97, rama 4→5) se une a 1 (ratio con max pasado 3/10 → major) ;
#   nodo 6 (s98, rama 6) se une a 0 con ratio 1/10 (menor)
tree = {'SubhaloID': np.arange(7), 'SnapNum': np.array([99, 98, 97, 96, 97, 96, 98]),
        'FirstProgenitorID': np.array([1, 2, 3, -1, 5, -1, -1]),
        'NextProgenitorID': np.array([-1, 6, 4, -1, -1, -1, -1]),
        'MainLeafProgenitorID': np.array([3, 3, 3, 3, 5, 5, 6]),
        'SubhaloMassType': np.zeros((7, 6))}
tree['SubhaloMassType'][:, 4] = [12, 11, 10, 9, 2, 3, 1.1]
snap, ratio = F.last_major_merger_snap(tree, 0.2)
check("último major merger", snap == 98 and np.isclose(ratio, 0.3), f"snap={snap} ratio={ratio}")
tree['SubhaloMassType'][4:6, 4] = [1, 1.5]
check("sin major merger → -1", F.last_major_merger_snap(tree, 0.2)[0] == -1)

# %% 11. Gradiente con bins de igual número: Z ∝ R^-0.3 → pendiente −0.3 en log
Rg = 10**rng.uniform(0, 3, 300000); Zg = 0.02*Rg**-0.3
sl, er, _, xg, yg = F.radial_gradient(Rg, 30, lambda idx: F.mean_log_metallicity(Zg[idx], np.ones(len(idx)), 0.0127))
check("gradiente log", abs(sl + 0.3) < 0.01, f"pendiente={sl:.3f}±{er:.3f}")

# %% 12. Componentes SA, progenitores compartidos y significativos
comp = F.classify_components(np.array([1, 0, 0, 0, -1, 0]), np.array([-1, 0, 2, 1, -1, -1]))
check("componentes", list(comp) == [0, 1, 2, 2, -1, -1], str(comp))
cf = F.component_fractions(np.ones(6), comp)
check("fracciones componentes", np.isclose(cf[0], 0.25) and np.isclose(cf[2], 0.5) and np.isclose(cf['f_sin_clasificar'], 2/6))
check("shared_fraction", np.isclose(F.shared_fraction([1, 1, 2, 3, -1], [2, 3, 5]), 0.5))
top = F.significant_progenitors(np.array([7, 7, 7, 8, 9, -1]), np.array([1, 1, 1, 2, 1, 5.]))
check("progenitor más significativo", top[0] == 7 and np.isclose(top[1], 0.5) and top[2] == 3, str(top))

# %% 13. Subhalo por índice local del FoF y despegue del central
lens = np.array([5, 3, 2]); loc = np.array([0, 4, 5, 7, 8, 9, 10, 12])
check("subhalo_of_local_index", list(F.subhalo_of_local_index(loc, lens, 100)) == [100, 100, 101, 101, 102, 102, -1, -1])
ids = np.array([10, 11, 12, 13]); done = np.zeros(4, bool); out = np.full(4, -1)
F.detach_update(ids, done, np.array([10, 12, 13]), 98, out)
F.detach_update(ids, done, np.array([10]), 97, out)
check("detach_update", list(out) == [-1, 98, 97, 97], str(out))

# %% 14. Color integrado y radio de semi-masa
check("color integrado", np.isclose(F.integrated_color(np.array([1., 1.]), np.array([0., 0.])), 1.0))
check("half mass radius", abs(F.half_mass_radius(np.linspace(0, 1, 100001), np.ones(100001)) - 0.5) < 1e-4)

print(f"\n{sum(OK)}/{len(OK)} pruebas OK")
