from common import *

cells = [
md(r'''
# 02 · Componentes de ensamblaje y propiedades estelares del BCG y el ICL
**Metodología:** §2.3, §4, §6 de Mayes+2026 · Figs. 6–12

- Cada partícula estelar del central (r₃D ≤ R₂₀₀c, sin viento) se etiqueta con el catálogo de ensamblaje
  **por partícula** de Rodriguez-Gomez+2016: **in situ**, **mergers completados** o **sobrevivientes**
  (acretadas de galaxias que siguen existiendo a z = 0).
- BCG / ICL con el **mismo R_H y la misma rotación** guardados por `01_` (no se recalculan).
- Metalicidad, color B−V (Buser, Vega) y edad del ICL, del BCG y de cada componente.
- Gradientes: cuña de 15° sobre el eje mayor, 30 bins de igual número de partículas hasta R₂₀₀c,
  ajuste lineal en log R (BCG+ICL y sólo ICL).

`# SUPUESTO:` ruta y nombres de campos del catálogo por partícula (`params_icl.py`); HMR en 3D;
medias pesadas en masa; color integrado (suma de flujos) por bin/componente.
Corre sólo en el TNG Lab. Requiere `00` y `01`.
'''),
code(REPO_CELL),
code(r'''
from astropy.cosmology import FlatLambdaCDM

import numpy as np
import h5py
import matplotlib.pyplot as plt
import seaborn as sns

import illustris_python as il
import params_icl as P
import icl_funcs as F
'''),
code(STYLE_CELL),
code(fig_dirs_cell('02_propiedades_estelares')),
code(COSMO_CELL),

md('## 1 · Parámetros'),
code(r'''
SA_PATH = P.SA_STARS_FILE.format(snap=P.SNAP)
print(f"Catálogo SA por partícula : {SA_PATH}")
print(f"Campos                    : {P.SA_FIELD_INSITU}, {P.SA_FIELD_ORIGIN}")
print(f"Cuña                      : {P.WEDGE_DEG}°   bins de igual número: {P.N_BINS_GRAD}")
print(f"Z_sun                     : {P.Z_SUN}")
'''),
code(r'''
# Chequeo temprano del catálogo por partícula: falla aquí con un mensaje claro si no está
if not os.path.exists(SA_PATH):
    raise FileNotFoundError(f"No está {SA_PATH}. Ajusta P.SA_STARS_FILE o TNG_SA_STARS_FILE "
                            "(ls /home/tnguser/sims.TNG/TNG100-1/postprocessing/StellarAssembly/).")
with h5py.File(SA_PATH, 'r') as f:
    sa_keys = list(f.keys())
    sa_len  = {k: f[k].shape for k in sa_keys if isinstance(f[k], h5py.Dataset)}
print(f"Campos en el archivo: {sa_len}")
'''),
code(r'''
with h5py.File(il.snapshot.snapPath(P.basePath, P.SNAP), 'r') as f:
    npart = f['Header'].attrs['NumPart_Total']
    hword = f['Header'].attrs['NumPart_Total_HighWord']
n_star_snap = int(npart[4]) + (int(hword[4]) << 32)
print(f"PartType4 en el snapshot : {n_star_snap:,}")
print(f"Largo del catálogo SA    : {sa_len.get(P.SA_FIELD_INSITU)}")
print("(deben coincidir: el catálogo se asume en el orden de PartType4)")
'''),

md('## 2 · Cargar catálogo (salidas de `00` y `01`)'),
code(r'''
with h5py.File(P.CATALOG_OUT, 'r') as f:
    group_idx     = f['group_idx'][:]
    M200c         = f['M200c_Msun'][:]
    R200c         = f['R200c_kpc'][:]
    GroupPos      = f['GroupPos_kpc'][:]
    bcg_sub_idx   = f['bcg_sub_idx'][:]
    dyn_state     = f['dyn_state'][:]
    R_H           = f['R_H_kpc'][:]
    R_rot         = f['R_rot'][:]

n_cl = len(group_idx)
Header  = il.groupcat.loadHeader(P.basePath, P.SNAP)
box_kpc = Header['BoxSize'] * P.UL
print(f"Catálogo: {n_cl} halos   con R_H: {np.isfinite(R_H).sum()}")
'''),
code(MASS_LABEL_CELL),

md('## 3 · Funciones auxiliares'),
code(r'''
def load_cluster(i):
    """Estrellas del central (≤ R200, sin viento) con la rotación y R_H de 01_, componente SA, Z y edad."""
    st = F.load_central_stars(P.basePath, P.SNAP, bcg_sub_idx[i], GroupPos[i], R200c[i], box_kpc,
                              P.UL, P.UM, extra_fields=['GFM_Metallicity'], r_outer_frac=P.R_OUTER_FRAC)
    insitu, origin = F.load_sa_stars(P.SA_STARS_FILE, P.SNAP, st['offset'], st['length'],
                                     P.SA_FIELD_INSITU, P.SA_FIELD_ORIGIN)
    st['comp']    = F.classify_components(insitu[st['idx_sub']], origin[st['idx_sub']])
    st['pos_rot'] = st['pos'] @ R_rot[i].T
    st['R']       = F.projected_radius(st['pos_rot'])
    st['mask_bcg'], st['mask_icl'] = F.split_bcg_icl(st['R'], R_H[i])
    st['age']     = F.lookback_from_scale_factor(st['sft'], cosmo)       # Gyr
    return st


def mean_Z(st, m):
    return F.mean_log_metallicity(st['GFM_Metallicity'][m], st['mass'][m], P.Z_SUN)


def color_BV(st, m):
    return F.integrated_color(st['phot'][m, P.BV_IDX[0]], st['phot'][m, P.BV_IDX[1]])


def mean_age(st, m):
    return np.average(st['age'][m], weights=st['mass'][m]) if m.sum() > 0 else np.nan


def scatter_beta(ax, x, y, color, label, fmt='.2f'):
    """Puntos + ajuste lineal con β (± error) en la leyenda."""
    b, eb, c = F.linfit_beta(x, y)
    ax.scatter(x, y, s=18, color=color, edgecolor='k', linewidth=0.2, alpha=0.8,
               label=rf'{label}: $\beta$ = {b:{fmt}} ± {eb:{fmt}}')
    xx = np.linspace(np.nanmin(x), np.nanmax(x), 50)
    ax.plot(xx, b * xx + c, color=color, lw=1.5)
'''),

md('## 4 · Figura 6 – distribución radial de las tres componentes en un cúmulo masivo'),
code(r'''
i_demo  = np.argsort(M200c)[::-1][2]     # mismo cúmulo que las Figs. 1–2 (ver 01_)
st_demo = load_cluster(i_demo)
print(f"FoF {group_idx[i_demo]}: {len(st_demo['mass']):,} estrellas   sin clasificar: "
      f"{(st_demo['comp'] < 0).mean():.3%}")
'''),
code(r'''
bins_r = np.linspace(0, 1, 41)
x_demo = st_demo['r3d'] / R200c[i_demo]

fig, ax = plt.subplots(figsize=(7, 5))
for k in [0, 1, 2]:
    m = st_demo['comp'] == k
    ax.hist(x_demo[m], bins=bins_r, weights=st_demo['mass'][m], histtype='step', lw=2,
            color=COLORS_COMP[k], label=F.COMP_LABELS[k])
ax.axvline(R_H[i_demo] / R200c[i_demo], color='k', ls='--', lw=1.2, label='R_H (proyectado)')
ax.set_yscale('log')
ax.set_xlabel(r'$r_{3D}\,/\,R_{200}$')
ax.set_ylabel(r'$M_*$ por bin [$M_\odot$]')
ax.legend()

plt.savefig(f'{FIG_PDF}/fig06_componentes_radial_fof{group_idx[i_demo]}.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig06_componentes_radial_fof{group_idx[i_demo]}.png', dpi=150, bbox_inches='tight')
plt.show()
'''),

md('## 5 · Propiedades por cúmulo'),
code(r'''
f_comp_bcgicl = np.full((n_cl, 3), np.nan)    # fracciones de masa (in situ, mergers, sobrev.)
f_comp_icl    = np.full((n_cl, 3), np.nan)
f_comp_bcg    = np.full((n_cl, 3), np.nan)
f_unclass     = np.full(n_cl, np.nan)
hmr_bcgicl    = np.full((n_cl, 3), np.nan)    # kpc, 3D
hmr_icl       = np.full((n_cl, 3), np.nan)
Z_icl, Z_bcg  = np.full(n_cl, np.nan), np.full(n_cl, np.nan)
BV_icl, BV_bcg = np.full(n_cl, np.nan), np.full(n_cl, np.nan)
age_icl, age_bcg = np.full(n_cl, np.nan), np.full(n_cl, np.nan)
Z_comp, BV_comp, age_comp = (np.full((n_cl, 3), np.nan) for _ in range(3))
grad_Z_bcgicl, grad_Z_icl   = np.full(n_cl, np.nan), np.full(n_cl, np.nan)
grad_BV_bcgicl, grad_BV_icl = np.full(n_cl, np.nan), np.full(n_cl, np.nan)
'''),
code(r'''
for i in range(n_cl):
    if not np.isfinite(R_H[i]):
        continue
    try:
        st = load_cluster(i)
    except (FileNotFoundError, KeyError, ValueError):
        raise          # problemas con el catálogo SA: detener con su mensaje
    except Exception as e:
        print(f"  [WARN] i={i}: {e}")
        continue
    comp, mb, mi = st['comp'], st['mask_bcg'], st['mask_icl']
    allm = np.ones(len(comp), dtype=bool)

    for arr, msk in [(f_comp_bcgicl, allm), (f_comp_icl, mi), (f_comp_bcg, mb)]:
        fr = F.component_fractions(st['mass'], comp, msk)
        arr[i] = [fr[0], fr[1], fr[2]]
    f_unclass[i] = F.component_fractions(st['mass'], comp)['f_sin_clasificar']

    for k in [0, 1, 2]:
        mk = comp == k
        hmr_bcgicl[i, k] = F.half_mass_radius(st['r3d'][mk], st['mass'][mk])
        hmr_icl[i, k]    = F.half_mass_radius(st['r3d'][mk & mi], st['mass'][mk & mi])
        Z_comp[i, k], BV_comp[i, k], age_comp[i, k] = mean_Z(st, mk), color_BV(st, mk), mean_age(st, mk)

    Z_icl[i], Z_bcg[i]     = mean_Z(st, mi), mean_Z(st, mb)
    BV_icl[i], BV_bcg[i]   = color_BV(st, mi), color_BV(st, mb)
    age_icl[i], age_bcg[i] = mean_age(st, mi), mean_age(st, mb)

    # gradientes en la cuña del eje mayor (BCG+ICL y sólo ICL)
    w = F.wedge_mask(st['pos_rot'], P.WEDGE_DEG)
    for sel, gZ, gBV in [(w, grad_Z_bcgicl, grad_BV_bcgicl), (w & mi, grad_Z_icl, grad_BV_icl)]:
        idx_sel = np.nonzero(sel)[0]
        Rw = st['R'][idx_sel]
        gZ[i]  = F.radial_gradient(Rw, P.N_BINS_GRAD,
                                   lambda j: F.mean_log_metallicity(st['GFM_Metallicity'][idx_sel[j]],
                                                                    st['mass'][idx_sel[j]], P.Z_SUN))[0]
        gBV[i] = F.radial_gradient(Rw, P.N_BINS_GRAD,
                                   lambda j: F.integrated_color(st['phot'][idx_sel[j], P.BV_IDX[0]],
                                                                st['phot'][idx_sel[j], P.BV_IDX[1]]))[0]
    if (i + 1) % 10 == 0:
        print(f"  {i+1}/{n_cl}", end="\r")
print("\nListo.")
'''),
code(r'''
print(f"Masa sin clasificar (mediana) = {np.nanmedian(f_unclass):.2%}")
print(f"Todos los gradientes de Z negativos (BCG+ICL): {(grad_Z_bcgicl[np.isfinite(grad_Z_bcgicl)] < 0).all()}")
print(f"Todos los gradientes B-V negativos (BCG+ICL) : {(grad_BV_bcgicl[np.isfinite(grad_BV_bcgicl)] < 0).all()}")
'''),

md('## 6 · Comparación con los valores publicados (media ± 1σ)'),
code(r'''
pv = P.PAPER_VALUES
filas = [('BCG+ICL', f_comp_bcgicl, ['bcgicl_insitu', 'bcgicl_mergers', 'bcgicl_surviving']),
         ('ICL',     f_comp_icl,    ['icl_insitu', 'icl_mergers', 'icl_surviving']),
         ('BCG',     f_comp_bcg,    [None, None, 'bcg_surviving'])]
print(f"{'región':<8} {'componente':<20} {'este análisis':>16} {'Mayes+2026':>14}")
for region, arr, keys in filas:
    for k in [0, 1, 2]:
        ref = f"{pv[keys[k]][0]:.2f} ± {pv[keys[k]][1]:.2f}" if keys[k] else '—'
        print(f"{region:<8} {F.COMP_LABELS[k]:<20} {np.nanmean(arr[:, k]):>8.2f} ± {np.nanstd(arr[:, k]):.2f} {ref:>14}")
'''),

md('## 7 · Guardar en el catálogo'),
code(r'''
nuevos = {'f_comp_bcgicl': f_comp_bcgicl, 'f_comp_icl': f_comp_icl, 'f_comp_bcg': f_comp_bcg,
          'f_sa_unclassified': f_unclass, 'hmr_comp_bcgicl_kpc': hmr_bcgicl, 'hmr_comp_icl_kpc': hmr_icl,
          'logZ_icl': Z_icl, 'logZ_bcg': Z_bcg, 'logZ_comp': Z_comp,
          'BV_icl': BV_icl, 'BV_bcg': BV_bcg, 'BV_comp': BV_comp,
          'age_icl_Gyr': age_icl, 'age_bcg_Gyr': age_bcg, 'age_comp_Gyr': age_comp,
          'grad_logZ_bcgicl': grad_Z_bcgicl, 'grad_logZ_icl': grad_Z_icl,
          'grad_BV_bcgicl': grad_BV_bcgicl, 'grad_BV_icl': grad_BV_icl}
with h5py.File(P.CATALOG_OUT, 'a') as f:
    for k, v in nuevos.items():
        if k in f:
            del f[k]
        f.create_dataset(k, data=v)
    f['f_comp_bcgicl'].attrs['columnas'] = 'in situ, mergers completados, sobrevivientes'
print(f"Guardado en {P.CATALOG_OUT}")
'''),

md('## 8 · Figura 7 – radio de semi-masa (3D) de cada componente'),
code(r'''
fig, ax = plt.subplots(2, 1, figsize=(7, 9), sharex=True)
for k in [0, 1, 2]:
    scatter_beta(ax[0], lM200, np.log10(hmr_bcgicl[:, k]), COLORS_COMP[k], F.COMP_LABELS[k])
    scatter_beta(ax[1], lM200, np.log10(hmr_icl[:, k]), COLORS_COMP[k], F.COMP_LABELS[k])
ax[0].set_ylabel(r'$\log_{10}(r_{1/2}$ BCG+ICL / kpc)')
ax[1].set_ylabel(r'$\log_{10}(r_{1/2}$ ICL / kpc)')
ax[1].set_xlabel(rf'$\log_{{10}}(M_{{200}}\,/\,{MASS_UNIT})$')
ax[0].legend(fontsize=9)
ax[1].legend(fontsize=9)

plt.savefig(f'{FIG_PDF}/fig07_hmr_componentes.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig07_hmr_componentes.png', dpi=150, bbox_inches='tight')
plt.show()
'''),

md('## 9 · Figura 8 – fracciones de masa de las componentes del BCG+ICL'),
code(r'''
fig, ax = plt.subplots(figsize=(7, 5))
for k in [0, 1, 2]:
    scatter_beta(ax, lM200, f_comp_bcgicl[:, k], COLORS_COMP[k], F.COMP_LABELS[k])
ax.set_xlabel(rf'$\log_{{10}}(M_{{200}}\,/\,{MASS_UNIT})$')
ax.set_ylabel('Fracción de masa del BCG+ICL')
ax.legend(fontsize=9)

plt.savefig(f'{FIG_PDF}/fig08_componentes_bcgicl.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig08_componentes_bcgicl.png', dpi=150, bbox_inches='tight')
plt.show()
'''),

md('## 10 · Figura 9 – componentes del ICL (arriba) y del BCG (abajo)'),
code(r'''
fig, ax = plt.subplots(2, 1, figsize=(7, 9), sharex=True)
for k in [0, 1, 2]:
    scatter_beta(ax[0], lM200, f_comp_icl[:, k], COLORS_COMP[k], F.COMP_LABELS[k])
    scatter_beta(ax[1], lM200, f_comp_bcg[:, k], COLORS_COMP[k], F.COMP_LABELS[k])
ax[0].set_ylabel('Fracción de masa del ICL')
ax[1].set_ylabel('Fracción de masa del BCG')
ax[1].set_xlabel(rf'$\log_{{10}}(M_{{200}}\,/\,{MASS_UNIT})$')
ax[0].legend(fontsize=9)
ax[1].legend(fontsize=9)

plt.savefig(f'{FIG_PDF}/fig09_componentes_icl_bcg.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig09_componentes_icl_bcg.png', dpi=150, bbox_inches='tight')
plt.show()
'''),

md('## 11 · Figura 10 – metalicidad'),
code(r'''
fig, ax = plt.subplots(3, 1, figsize=(7, 13), sharex=True)
scatter_beta(ax[0], lM200, Z_icl, my_colors[0], 'ICL')
scatter_beta(ax[0], lM200, Z_bcg, my_colors[4], 'BCG')
for k in [0, 1, 2]:
    scatter_beta(ax[1], lM200, Z_comp[:, k], COLORS_COMP[k], F.COMP_LABELS[k])
scatter_beta(ax[2], lM200, grad_Z_bcgicl, my_colors[2], 'BCG+ICL', fmt='.3f')
scatter_beta(ax[2], lM200, grad_Z_icl, my_colors[0], 'ICL', fmt='.3f')
ax[0].set_ylabel(r'$\langle \log Z_*/Z_\odot \rangle$')
ax[1].set_ylabel(r'$\langle \log Z_*/Z_\odot \rangle$')
ax[2].set_ylabel(r'$d \log Z\, /\, d \log R$')
ax[2].set_xlabel(rf'$\log_{{10}}(M_{{200}}\,/\,{MASS_UNIT})$')
for a in ax:
    a.legend(fontsize=8)

plt.savefig(f'{FIG_PDF}/fig10_metalicidad.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig10_metalicidad.png', dpi=150, bbox_inches='tight')
plt.show()
'''),

md('## 12 · Figura 11 – color B−V'),
code(r'''
fig, ax = plt.subplots(3, 1, figsize=(7, 13), sharex=True)
scatter_beta(ax[0], lM200, BV_icl, my_colors[0], 'ICL', fmt='.3f')
scatter_beta(ax[0], lM200, BV_bcg, my_colors[4], 'BCG', fmt='.3f')
for k in [0, 1, 2]:
    scatter_beta(ax[1], lM200, BV_comp[:, k], COLORS_COMP[k], F.COMP_LABELS[k], fmt='.3f')
scatter_beta(ax[2], lM200, grad_BV_bcgicl, my_colors[2], 'BCG+ICL', fmt='.3f')
scatter_beta(ax[2], lM200, grad_BV_icl, my_colors[0], 'ICL', fmt='.3f')
ax[0].set_ylabel(r'$B-V$')
ax[1].set_ylabel(r'$B-V$')
ax[2].set_ylabel(r'$d(B-V)\, /\, d \log R$')
ax[2].set_xlabel(rf'$\log_{{10}}(M_{{200}}\,/\,{MASS_UNIT})$')
for a in ax:
    a.legend(fontsize=8)

plt.savefig(f'{FIG_PDF}/fig11_color_BV.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig11_color_BV.png', dpi=150, bbox_inches='tight')
plt.show()
'''),

md('## 13 · Figura 12 – edades (lookback de formación)'),
code(r'''
fig, ax = plt.subplots(2, 1, figsize=(7, 9), sharex=True)
scatter_beta(ax[0], lM200, age_icl, my_colors[0], 'ICL')
scatter_beta(ax[0], lM200, age_bcg, my_colors[4], 'BCG')
for k in [0, 1, 2]:
    scatter_beta(ax[1], lM200, age_comp[:, k], COLORS_COMP[k], F.COMP_LABELS[k])
ax[0].set_ylabel('Edad media [Gyr]')
ax[1].set_ylabel('Edad media [Gyr]')
ax[1].set_xlabel(rf'$\log_{{10}}(M_{{200}}\,/\,{MASS_UNIT})$')
ax[0].legend(fontsize=9)
ax[1].legend(fontsize=9)

plt.savefig(f'{FIG_PDF}/fig12_edades.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig12_edades.png', dpi=150, bbox_inches='tight')
plt.show()
'''),
]
if __name__ == "__main__":
    write_nb('02_propiedades_estelares.ipynb', cells)
