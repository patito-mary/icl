from common import *

IMPORTS = r'''
from astropy.cosmology import FlatLambdaCDM

import numpy as np
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import binned_statistic_2d

import illustris_python as il
import params_icl as P
import icl_funcs as F
'''

LOAD_CAT = r'''
with h5py.File(P.CATALOG_OUT, 'r') as f:
    group_idx     = f['group_idx'][:]
    M200c         = f['M200c_Msun'][:]
    R200c         = f['R200c_kpc'][:]
    GroupPos      = f['GroupPos_kpc'][:]
    bcg_sub_idx   = f['bcg_sub_idx'][:]
    Mstar200      = f['Mstar200_Msun'][:]
    dyn_state     = f['dyn_state'][:]
    concentration = f['concentration'][:]
    t_last_merger = f['t_last_merger_Gyr'][:]

n_cl = len(group_idx)
'''

cells = [
md(r'''
# 01 · Separación BCG/ICL y fracción de masa del ICL
**Metodología:** §2.2 y §3 de Mayes+2026 · Figs. 1–5

1. Estrellas del **subhalo central** (BCG+ICL según Subfind; sin satélites), sin viento, con r₃D ≤ R₂₀₀c.
2. Rotación por el tensor de inercia y proyección sobre el plano (eje mayor, eje intermedio).
3. Perfil 1D de brillo superficial en **banda r** (SDSS, AB) hasta R₂₀₀c y **ajuste logarítmico** μ(R) = a + b·log₁₀R.
4. Radio de Holmberg R_H: donde el ajuste cruza μ_r = 26.5 mag/arcsec². BCG: R ≤ R_H; ICL: R_H < R (y r₃D ≤ R₂₀₀c).
5. **f_ICL = M_ICL / M*₂₀₀** (definición del paper); se guarda también M_ICL/(M_BCG+M_ICL).

Los supuestos (binning, rango del ajuste, proyección, R₂₀₀ en 3D) están en `params_icl.py` marcados `# SUPUESTO:`.
Corre sólo en el TNG Lab.
'''),
code(REPO_CELL),
code(IMPORTS),
code(STYLE_CELL),
code(fig_dirs_cell('01_separacion_BCG_ICL')),
code(COSMO_CELL),

md('## 1 · Parámetros del perfil (definidos en `params_icl.py`)'),
code(r'''
print(f"Banda          = {P.BAND_NAME} (índice {P.BAND_IDX}, AB)   M_sun = {P.M_SUN_R_AB}")
print(f"mu_Holmberg    = {P.MU_HOLMBERG} mag/arcsec2")
print(f"Perfil         = {P.PROFILE_METHOD}, binning {P.PROFILE_BINNING}, ancho {P.SB_BIN_KPC} kpc")
print(f"R_H            = {P.RH_METHOD}  (ajuste mu = a + b log10 R desde {P.R_FIT_MIN_KPC} kpc hasta R200)")
print(f"Límite externo = {P.R_OUTER_FRAC} R200 (3D)")
'''),

md('## 2 · Cargar catálogo (salida de `00_catalogo`)'),
code(LOAD_CAT),
code(r'''
Header  = il.groupcat.loadHeader(P.basePath, P.SNAP)
box_kpc = Header['BoxSize'] * P.UL
print(f"Catálogo: {n_cl} halos   (paper: {P.PAPER_VALUES['N_halos']})")
'''),
code(MASS_LABEL_CELL),

md(r'''
## 3 · Demostración con un cúmulo (Figs. 1 y 2)
`# SUPUESTO:` el cúmulo de las Figs. 1–2 del paper parece ser el tercero más masivo (su mapa llega a ±1465 kpc ≈ R₂₀₀ de ese halo).
'''),
code(r'''
i_demo = np.argsort(M200c)[::-1][2]
print(f"Cúmulo demo: i = {i_demo}   FoF = {group_idx[i_demo]}   sub = {bcg_sub_idx[i_demo]}")
print(f"log M200c = {lM200[i_demo]:.2f} [{MASS_UNIT}]   R200c = {R200c[i_demo]:.0f} kpc")
'''),
code(r'''
stars_demo = F.load_central_stars(P.basePath, P.SNAP, bcg_sub_idx[i_demo], GroupPos[i_demo],
                                  R200c[i_demo], box_kpc, P.UL, P.UM, r_outer_frac=P.R_OUTER_FRAC)
print(f"Estrellas del central (r3D <= R200, sin viento): {len(stars_demo['mass']):,} de {stars_demo['n_total']:,}")
'''),
code(r'''
res_demo = F.bcg_icl_from_stars(stars_demo, Mstar200[i_demo], R200c[i_demo], P)
'''),
code(r'''
print(f"Ajuste: mu = {res_demo['a']:.2f} + {res_demo['b']:.2f} log10(R)   ({res_demo['n_fit']} anillos)")
print(f"R_H (ajuste)    = {res_demo['R_H_fit']:.1f} kpc")
print(f"R_H (interp)    = {res_demo['R_H_interp']:.1f} kpc")
print(f"R_H (local)     = {res_demo['R_H_local']:.1f} kpc")
print(f"M_BCG           = {res_demo['M_bcg']:.3e} Msun")
print(f"M_ICL           = {res_demo['M_icl']:.3e} Msun")
print(f"f_ICL (M*200)   = {res_demo['f_icl']:.3f}")
print(f"f_ICL (BCG+ICL) = {res_demo['f_icl_bcgicl']:.3f}")
'''),

md('### Figura 1 – Perfil 1D de brillo superficial en banda r (= Fig. 1 del paper)'),
code(r'''
rr = np.linspace(1, R200c[i_demo], 500)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(res_demo['r_mid'], res_demo['mu'], 'o', ms=4, color=my_colors[0], label='Perfil 1D')
ax.plot(rr, res_demo['a'] + res_demo['b'] * np.log10(rr), '-', color='k', lw=1.5,
        label='Ajuste logarítmico')
ax.axhline(P.MU_HOLMBERG, color='gray', ls='--', lw=1)
ax.axvline(res_demo['R_H'], color='r', lw=1.5, label=f"R_H = {res_demo['R_H']:.0f} kpc")
ax.set_ylim(32, 15)
ax.set_xlabel('Radio [kpc]')
ax.set_ylabel(r'$\mu_r$ [mag arcsec$^{-2}$]')
ax.set_title(f'FoF {group_idx[i_demo]}   snap {P.SNAP}')
ax.legend()

plt.savefig(f'{FIG_PDF}/fig01_perfil_sb_fof{group_idx[i_demo]}.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig01_perfil_sb_fof{group_idx[i_demo]}.png', dpi=150, bbox_inches='tight')
plt.show()
'''),

md('### Figura 2 – Mapa 2D de brillo superficial hasta R₂₀₀ (= Fig. 2 del paper)'),
code(r'''
r_plot  = R200c[i_demo]
n_pix   = 512
edges   = np.linspace(-r_plot, r_plot, n_pix + 1)
pix_pc2 = ((2 * r_plot / n_pix) * 1e3)**2

L_map, _, _, _ = binned_statistic_2d(res_demo['pos_rot'][:, 0], res_demo['pos_rot'][:, 1],
                                     res_demo['lum'], statistic='sum', bins=[edges, edges])
mu_map = F.mu_from_sigma(L_map / pix_pc2, P.SB_CONST)
'''),
code(r'''
fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(mu_map.T, origin='lower', cmap='ocean', vmin=16, vmax=32,
               extent=[-r_plot, r_plot, -r_plot, r_plot])
plt.colorbar(im, ax=ax, label=r'$\mu_r$ [mag arcsec$^{-2}$]')
ax.add_patch(plt.Circle((0, 0), res_demo['R_H'], fill=False, color='k', lw=1.5, label='R_H'))
ax.add_patch(plt.Circle((0, 0), r_plot, fill=False, color='gray', lw=1, ls='--', label='R$_{200}$'))
ax.set_xlabel('x [kpc]')
ax.set_ylabel('y [kpc]')
ax.grid(False)
ax.legend(loc='upper right')

plt.savefig(f'{FIG_PDF}/fig02_mapa_sb_fof{group_idx[i_demo]}.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig02_mapa_sb_fof{group_idx[i_demo]}.png', dpi=150, bbox_inches='tight')
plt.show()
'''),

md('## 4 · R_H y f_ICL para toda la muestra'),
code(r'''
R_H        = np.full(n_cl, np.nan)
R_H_fit    = np.full(n_cl, np.nan)
R_H_interp = np.full(n_cl, np.nan)
R_H_local  = np.full(n_cl, np.nan)      # diagnóstico: ajuste sólo en 26.5 ± 1 mag
sb_a       = np.full(n_cl, np.nan)
sb_b       = np.full(n_cl, np.nan)
M_bcg      = np.full(n_cl, np.nan)
M_icl      = np.full(n_cl, np.nan)
f_icl      = np.full(n_cl, np.nan)       # M_ICL / M*200
f_icl_bcgicl = np.full(n_cl, np.nan)     # M_ICL / (M_BCG + M_ICL)
R_rot      = np.full((n_cl, 3, 3), np.nan)
'''),
code(r'''
for i in range(n_cl):
    try:
        st  = F.load_central_stars(P.basePath, P.SNAP, bcg_sub_idx[i], GroupPos[i], R200c[i],
                                   box_kpc, P.UL, P.UM, r_outer_frac=P.R_OUTER_FRAC)
        res = F.bcg_icl_from_stars(st, Mstar200[i], R200c[i], P)
    except Exception as e:
        print(f"  [WARN] i={i} sub={bcg_sub_idx[i]}: {e}")
        continue
    R_H[i], R_H_fit[i], R_H_interp[i] = res['R_H'], res['R_H_fit'], res['R_H_interp']
    R_H_local[i] = res['R_H_local']
    sb_a[i], sb_b[i] = res['a'], res['b']
    M_bcg[i], M_icl[i] = res['M_bcg'], res['M_icl']
    f_icl[i], f_icl_bcgicl[i] = res['f_icl'], res['f_icl_bcgicl']
    R_rot[i] = res['R_mat']
    if (i + 1) % 10 == 0:
        print(f"  {i+1}/{n_cl}  f_ICL_med = {np.nanmedian(f_icl[:i+1]):.3f}", end="\r")
'''),
code(r'''
ok = np.isfinite(f_icl)
print(f"Con R_H válido           = {ok.sum()}/{n_cl}")
print(f"f_ICL (M*200) mediana    = {np.nanmedian(f_icl):.3f}   (paper: ~{P.PAPER_VALUES['f_icl_mediana_aprox']})")
print(f"f_ICL 5–95 %             = {np.nanpercentile(f_icl, 5):.3f} – {np.nanpercentile(f_icl, 95):.3f}   (paper: 0.05–0.25)")
print(f"f_ICL (BCG+ICL) mediana  = {np.nanmedian(f_icl_bcgicl):.3f}")
print(f"R_H mediana              = {np.nanmedian(R_H):.0f} kpc")
print(f"R_H ajuste / interp      = {np.nanmedian(R_H_fit / R_H_interp):.2f} (mediana)")
print(f"R_H ajuste / local       = {np.nanmedian(R_H_fit / R_H_local):.2f} (mediana; en la Fig. 1 del paper ≈ 0.87)")
'''),

md('## 5 · Guardar en el catálogo'),
code(r'''
nuevos = {'R_H_kpc': R_H, 'R_H_fit_kpc': R_H_fit, 'R_H_interp_kpc': R_H_interp, 'R_H_local_kpc': R_H_local,
          'sb_fit_a': sb_a, 'sb_fit_b': sb_b,
          'M_bcg_Msun': M_bcg, 'M_icl_Msun': M_icl,
          'f_icl': f_icl, 'f_icl_bcgicl': f_icl_bcgicl, 'R_rot': R_rot}

with h5py.File(P.CATALOG_OUT, 'a') as f:
    for k, v in nuevos.items():
        if k in f:
            del f[k]
        f.create_dataset(k, data=v)
    f['f_icl'].attrs['definicion']        = 'M_ICL / M*200 (Mayes+2026 §3)'
    f['f_icl_bcgicl'].attrs['definicion'] = 'M_ICL / (M_BCG + M_ICL)'
    f['R_rot'].attrs['uso']               = 'pos_rot = pos_rel @ R_rot.T (eje mayor = x)'
    f.attrs['perfil'] = (f"{P.PROFILE_METHOD}/{P.PROFILE_BINNING} {P.SB_BIN_KPC} kpc, banda {P.BAND_NAME}, "
                         f"R_H {P.RH_METHOD}, mu {P.MU_HOLMBERG}")
print(f"Guardado en {P.CATALOG_OUT}")
'''),

md('## 6 · Figura 3 – f_ICL vs M₂₀₀ por estado dinámico (= Fig. 3)'),
code(r'''
beta, beta_err, icpt = F.linfit_beta(lM200, f_icl)
xx = np.linspace(np.nanmin(lM200), np.nanmax(lM200), 100)
'''),
code(r'''
fig = plt.figure(figsize=(8, 7))
gs  = fig.add_gridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4], wspace=0.05, hspace=0.05)
ax    = fig.add_subplot(gs[1, 0])
ax_mx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_fy = fig.add_subplot(gs[1, 1], sharey=ax)

for s in [0, 1, 2]:
    m = (dyn_state == s) & ok
    ax.scatter(lM200[m], f_icl[m], s=30, color=COLORS_STATE[s], edgecolor='k', linewidth=0.3,
               label=P.DYN_LABELS[s])
    ax_mx.hist(lM200[m], bins=12, histtype='step', lw=2, color=COLORS_STATE[s])
    ax_mx.axvline(np.median(lM200[m]), color=COLORS_STATE[s], ls='--', lw=1)
    ax_fy.hist(f_icl[m], bins=12, histtype='step', lw=2, color=COLORS_STATE[s], orientation='horizontal')
    ax_fy.axhline(np.median(f_icl[m]), color=COLORS_STATE[s], ls='--', lw=1)

ax.plot(xx, beta * xx + icpt, 'k-', lw=1.5)
ax.text(0.5, 0.04, rf'$\beta$ = {beta:.3f} ± {beta_err:.3f}', transform=ax.transAxes, ha='center')
ax.set_xlabel(rf'$\log_{{10}}(M_{{200}}\,/\,{MASS_UNIT})$')
ax.set_ylabel(r'$f_{\rm ICL} = M_{\rm ICL}/M^{*}_{200}$')
ax.legend(loc='upper left')
ax_mx.tick_params(labelbottom=False)
ax_fy.tick_params(labelleft=False)

plt.savefig(f'{FIG_PDF}/fig03_ficl_vs_masa.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig03_ficl_vs_masa.png', dpi=150, bbox_inches='tight')
plt.show()
'''),

md('## 7 · Figura 4 – f_ICL vs concentración (= Fig. 4)'),
code(r'''
beta_c, beta_c_err, icpt_c = F.linfit_beta(concentration, f_icl)
xc = np.linspace(np.nanmin(concentration), np.nanmax(concentration), 100)

fig, ax = plt.subplots(figsize=(7, 5))
sc = ax.scatter(concentration, f_icl, c=lM200, cmap='viridis', s=30, edgecolor='k', linewidth=0.3)
plt.colorbar(sc, ax=ax, label=rf'$\log_{{10}}(M_{{200}}\,/\,{MASS_UNIT})$')
ax.plot(xc, beta_c * xc + icpt_c, 'k-', lw=1.5)
ax.text(0.05, 0.04, rf'$\beta$ = {beta_c:.3f} ± {beta_c_err:.3f}', transform=ax.transAxes)
ax.set_xlabel(r'Concentración $c_{200}$')
ax.set_ylabel(r'$f_{\rm ICL}$')

plt.savefig(f'{FIG_PDF}/fig04_ficl_vs_concentracion.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig04_ficl_vs_concentracion.png', dpi=150, bbox_inches='tight')
plt.show()
'''),

md('## 8 · Figura 5 – f_ICL vs lookback time del último major merger (= Fig. 5)'),
code(r'''
beta_t, beta_t_err, icpt_t = F.linfit_beta(t_last_merger, f_icl)
xt = np.linspace(np.nanmin(t_last_merger), np.nanmax(t_last_merger), 100)

fig, ax = plt.subplots(figsize=(7, 5))
sc = ax.scatter(t_last_merger, f_icl, c=lM200, cmap='viridis', s=30, edgecolor='k', linewidth=0.3)
plt.colorbar(sc, ax=ax, label=rf'$\log_{{10}}(M_{{200}}\,/\,{MASS_UNIT})$')
ax.plot(xt, beta_t * xt + icpt_t, 'k-', lw=1.5)
ax.text(0.95, 0.04, rf'$\beta$ = {beta_t:.4f} ± {beta_t_err:.4f}', transform=ax.transAxes, ha='right')
ax.set_xlabel('Lookback time último major merger [Gyr]')
ax.set_ylabel(r'$f_{\rm ICL}$')

plt.savefig(f'{FIG_PDF}/fig05_ficl_vs_merger.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig05_ficl_vs_merger.png', dpi=150, bbox_inches='tight')
plt.show()
'''),
code(r'''
print(f"beta f_ICL–M200    = {beta:+.3f} ± {beta_err:.3f}   (paper: positiva, con mucha dispersión)")
print(f"beta f_ICL–c200    = {beta_c:+.4f} ± {beta_c_err:.4f}   (paper: declive débil)")
print(f"beta f_ICL–t_merger= {beta_t:+.4f} ± {beta_t_err:.4f}   (paper: positiva)")
'''),
]
if __name__ == "__main__":
    write_nb('01_separacion_BCG_ICL.ipynb', cells)
