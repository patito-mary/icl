from common import *
from build_01 import IMPORTS, LOAD_CAT

cells = [
md(r'''
# 01e · Chequeo: BCG/ICL con el FoF local (central + fuzz) recortado a R₂₀₀c
**No es la definición del paper.** Mayes+2026 usa sólo el subhalo central (`01_separacion_BCG_ICL`).
Aquí se repite exactamente el mismo pipeline (mismas funciones de `icl_funcs.py`) agregando las
estrellas del FoF que no están ligadas a ningún subhalo (*fuzz*), para medir cuánto cambian R_H y f_ICL.

| Notebook | Partículas | Límite externo |
|---|---|---|
| `01_` | subhalo central | r₃D ≤ R₂₀₀c |
| **`01e_`** | **central + fuzz del FoF (satélites excluidos)** | **r₃D ≤ R₂₀₀c** |
| `01f_` | cutout de la API (mismo criterio que `01_` por defecto) | r₃D ≤ R₂₀₀c |

Corre sólo en el TNG Lab. Requiere haber corrido `00` y `01`.
'''),
code(REPO_CELL),
code(IMPORTS),
code(STYLE_CELL),
code(fig_dirs_cell('01e_BCG_ICL_R200_local')),
code(COSMO_CELL),

md('## 1 · Parámetros de la máscara (los usan la demo y el loop)'),
code(r'''
INCLUDE_SATELLITES = False   # satélites siempre fuera (como en el paper)
INCLUDE_FUZZ       = True    # True = central + fuzz (chequeo); False = sólo central (= 01_)
SUFIJO             = '_fof'  # sufijo de las columnas en el catálogo
'''),

md('## 2 · Cargar catálogo y longitudes de subhalos'),
code(LOAD_CAT),
code(r'''
with h5py.File(P.CATALOG_OUT, 'r') as f:
    R_H_01   = f['R_H_kpc'][:]
    f_icl_01 = f['f_icl'][:]
'''),
code(r'''
Header  = il.groupcat.loadHeader(P.basePath, P.SNAP)
box_kpc = Header['BoxSize'] * P.UL

halo_cat       = il.groupcat.loadHalos(P.basePath, P.SNAP, fields=['GroupNsubs', 'GroupFirstSub'])
SubhaloLenType = il.groupcat.loadSubhalos(P.basePath, P.SNAP, fields=['SubhaloLenType'])
'''),
code(r'''
def lens_star_of(fof_idx):
    """Número de estrellas de cada subhalo del FoF (sub 0 = central)."""
    first = halo_cat['GroupFirstSub'][fof_idx]
    nsub  = halo_cat['GroupNsubs'][fof_idx]
    return SubhaloLenType[first:first + nsub, 4].astype(int)
'''),
code(MASS_LABEL_CELL),

md('## 3 · Demostración con el mismo cúmulo que `01_`'),
code(r'''
i_demo = np.argsort(M200c)[::-1][2]
lens_demo = lens_star_of(group_idx[i_demo])
print(f"FoF {group_idx[i_demo]}: central {lens_demo[0]:,}  satélites {lens_demo[1:].sum():,}")
'''),
code(r'''
stars_fof = F.load_fof_stars(P.basePath, P.SNAP, group_idx[i_demo], lens_demo, GroupPos[i_demo],
                             R200c[i_demo], box_kpc, P.UL, P.UM,
                             include_satellites=INCLUDE_SATELLITES, include_fuzz=INCLUDE_FUZZ,
                             r_outer_frac=P.R_OUTER_FRAC)
stars_cen = F.load_central_stars(P.basePath, P.SNAP, bcg_sub_idx[i_demo], GroupPos[i_demo],
                                 R200c[i_demo], box_kpc, P.UL, P.UM, r_outer_frac=P.R_OUTER_FRAC)
print(f"Estrellas usadas: central {len(stars_cen['mass']):,}   central+fuzz {len(stars_fof['mass']):,}")
'''),
code(r'''
res_cen = F.bcg_icl_from_stars(stars_cen, Mstar200[i_demo], R200c[i_demo], P)
res_fof = F.bcg_icl_from_stars(stars_fof, Mstar200[i_demo], R200c[i_demo], P)
'''),
code(r'''
for name, r in [('central', res_cen), ('central+fuzz', res_fof)]:
    print(f"{name:<13} R_H = {r['R_H']:6.1f} kpc   f_ICL = {r['f_icl']:.3f}   f_ICL(BCG+ICL) = {r['f_icl_bcgicl']:.3f}")
print(f"catálogo 01_  R_H = {R_H_01[i_demo]:6.1f} kpc   f_ICL = {f_icl_01[i_demo]:.3f}")
'''),
code(r'''
fig, ax = plt.subplots(figsize=(7, 5))
for name, r, c in [('central (01_)', res_cen, my_colors[0]), ('central + fuzz', res_fof, my_colors[2])]:
    ax.plot(r['r_mid'], r['mu'], 'o', ms=3, color=c, label=name)
    ax.axvline(r['R_H'], color=c, lw=1.5)
ax.axhline(P.MU_HOLMBERG, color='gray', ls='--', lw=1)
ax.set_ylim(33, 15)
ax.set_xlabel('Radio [kpc]')
ax.set_ylabel(r'$\mu_r$ [mag arcsec$^{-2}$]')
ax.set_title(f'FoF {group_idx[i_demo]}')
ax.legend()

plt.savefig(f'{FIG_PDF}/fig01e_perfil_central_vs_fof{group_idx[i_demo]}.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig01e_perfil_central_vs_fof{group_idx[i_demo]}.png', dpi=150, bbox_inches='tight')
plt.show()
'''),

md('## 4 · R_H y f_ICL para toda la muestra (misma máscara que la demo)'),
code(r'''
R_H_fof   = np.full(n_cl, np.nan)
f_icl_fof = np.full(n_cl, np.nan)
f_bi_fof  = np.full(n_cl, np.nan)
M_bcg_fof = np.full(n_cl, np.nan)
M_icl_fof = np.full(n_cl, np.nan)

print(f"include_satellites = {INCLUDE_SATELLITES}   include_fuzz = {INCLUDE_FUZZ}")
for i in range(n_cl):
    try:
        st  = F.load_fof_stars(P.basePath, P.SNAP, group_idx[i], lens_star_of(group_idx[i]), GroupPos[i],
                               R200c[i], box_kpc, P.UL, P.UM, include_satellites=INCLUDE_SATELLITES,
                               include_fuzz=INCLUDE_FUZZ, r_outer_frac=P.R_OUTER_FRAC)
        res = F.bcg_icl_from_stars(st, Mstar200[i], R200c[i], P)
    except Exception as e:
        print(f"  [WARN] i={i}: {e}")
        continue
    R_H_fof[i], f_icl_fof[i], f_bi_fof[i] = res['R_H'], res['f_icl'], res['f_icl_bcgicl']
    M_bcg_fof[i], M_icl_fof[i] = res['M_bcg'], res['M_icl']
    if (i + 1) % 10 == 0:
        print(f"  {i+1}/{n_cl}", end="\r")
'''),
code(r'''
print(f"f_ICL mediana: central {np.nanmedian(f_icl_01):.3f}   central+fuzz {np.nanmedian(f_icl_fof):.3f}")
print(f"R_H(fuzz)/R_H(central) mediana = {np.nanmedian(R_H_fof / R_H_01):.3f}")
'''),
code(r'''
with h5py.File(P.CATALOG_OUT, 'a') as f:
    for k, v in [('R_H_kpc', R_H_fof), ('f_icl', f_icl_fof), ('f_icl_bcgicl', f_bi_fof),
                 ('M_bcg_Msun', M_bcg_fof), ('M_icl_Msun', M_icl_fof)]:
        key = k + SUFIJO
        if key in f:
            del f[key]
        f.create_dataset(key, data=v)
print(f"Guardado con sufijo '{SUFIJO}' en {P.CATALOG_OUT}")
'''),

md('## 5 · Comparación con `01_` (subhalo central)'),
code(r'''
fig, ax = plt.subplots(1, 2, figsize=(12, 5))

ax[0].scatter(f_icl_01, f_icl_fof, c=lM200, cmap='viridis', s=25, edgecolor='k', linewidth=0.3)
lim = [0, np.nanmax([f_icl_01, f_icl_fof]) * 1.05]
ax[0].plot(lim, lim, 'k--', lw=1)
ax[0].set_xlabel(r'$f_{\rm ICL}$ central (01_)')
ax[0].set_ylabel(r'$f_{\rm ICL}$ central + fuzz (01e_)')

ax[1].scatter(R_H_01, R_H_fof, c=lM200, cmap='viridis', s=25, edgecolor='k', linewidth=0.3)
lim = [0, np.nanmax([R_H_01, R_H_fof]) * 1.05]
ax[1].plot(lim, lim, 'k--', lw=1)
ax[1].set_xlabel(r'$R_H$ central [kpc]')
ax[1].set_ylabel(r'$R_H$ central + fuzz [kpc]')

plt.savefig(f'{FIG_PDF}/fig01e_comparacion_central_vs_fof.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig01e_comparacion_central_vs_fof.png', dpi=150, bbox_inches='tight')
plt.show()
'''),
]
if __name__ == "__main__":
    write_nb('01e_BCG_ICL_R200_local.ipynb', cells)
