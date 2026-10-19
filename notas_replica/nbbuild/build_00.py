from common import *

cells = [
md(r'''
# 00 · Catálogo de la muestra
**Metodología:** §2.1 y §3 de Mayes, Gómez & Monachesi 2026 (A&A 708, A124)

- Muestra: todos los halos FoF de TNG100-1 (snap 99) con M₂₀₀c ≥ 10¹³ M☉/h → **127 halos** (10 con M₂₀₀c > 10¹⁴).
  El paper rotula M☉, pero sólo el umbral en M☉/h reproduce sus números (ver `params_icl.py`).
- Sin cortes adicionales (no se excluyen aislados ni hay clasificación fósil en este notebook).
- M*₂₀₀: masa estelar dentro de R₂₀₀c (3D), todas las estrellas del FoF.
- Estado dinámico: offset 3D BCG – centro de masa del gas (relajado < 0.1 R₂₀₀, perturbado > 0.4 R₂₀₀; Aldás+2023).
- Concentración: c = R₂₀₀c/r_s con un ajuste NFW al perfil de densidad de la DM.
- Último major merger: razón de masa estelar > 1/5 en el árbol SubLink del BCG.

**Dónde corre:** la sección 2 (chequeo con el pickle local) corre también en la laptop; desde la sección 3 todo necesita los datos del TNG Lab.
'''),
code(REPO_CELL),
code(r'''
from astropy.cosmology import FlatLambdaCDM

import pickle
import numpy as np
import h5py
import matplotlib.pyplot as plt
import seaborn as sns

import params_icl as P
import icl_funcs as F
'''),
code(STYLE_CELL),
code(fig_dirs_cell('00_catalogo')),
code(COSMO_CELL),

md('## 1 · Parámetros'),
code(r'''
# Campos del catálogo de grupos/subhalos que usa este notebook
HALO_FIELDS = ['Group_M_Crit200', 'Group_R_Crit200', 'GroupPos', 'GroupFirstSub', 'GroupNsubs']
SUB_FIELDS  = ['SubhaloPos', 'SubhaloMassType', 'SubhaloFlag']
TREE_FIELDS = ['SubhaloID', 'NextProgenitorID', 'FirstProgenitorID', 'MainLeafProgenitorID',
               'SnapNum', P.MERGER_MASS_FIELD]

print(f"Umbral de masa : Group_M_Crit200 >= {P.M200_MIN_INTERNAL:g} [1e10 Msun/h]  "
      f"(= {P.M200_MIN_INTERNAL * P.UM_H:.2e} Msun/h = {P.M200_MIN_INTERNAL * P.UM:.3e} Msun)")
print(f"Estado dinámico: relajado < {P.DYN_RELAXED} R200, perturbado > {P.DYN_DISTURBED} R200")
print(f"Major merger   : razón estelar > {P.MAJOR_MERGER_RATIO:.2f}")
'''),

md(r'''
## 2 · Chequeo local de la muestra (pickle del código de shift)
Corre en la laptop. El pickle contiene 182 halos (umbral más bajo); contamos los que pasan el corte de la réplica.
'''),
code(r'''
if os.path.exists(P.CATALOG_PKL):
    with open(P.CATALOG_PKL, 'rb') as f:
        halos_shift = pickle.load(f)
    M_int = np.asarray(halos_shift['Group_M_Crit200'])          # 1e10 Msun/h
    n_ok  = (M_int >= P.M200_MIN_INTERNAL).sum()
    n_14  = (M_int * P.UM_H > 1e14).sum()
    print(f"Halos en el pickle          : {len(M_int)}")
    print(f"M200c >= 1e13 Msun/h        : {n_ok}   (paper: {P.PAPER_VALUES['N_halos']})")
    print(f"  de ellos > 1e14 Msun/h    : {n_14}   (paper: {P.PAPER_VALUES['N_halos_1e14']})")
    print(f"M200c máx                   : {M_int.max() * P.UM_H:.3e} Msun/h  (paper: 2.6e14)")
    print(f"Con el umbral en Msun físico: {(M_int * P.UM >= 1e13).sum()} halos")
else:
    print(f"No está {P.CATALOG_PKL}: se omite el chequeo local")
'''),

md('## 3 · Selección de la muestra desde el catálogo de grupos (Lab)'),
code(r'''
Header  = il.groupcat.loadHeader(P.basePath, P.SNAP)
box_kpc = Header['BoxSize'] * P.UL
print(f"Redshift = {Header['Redshift']:.3g}   BoxSize = {Header['BoxSize']:.0f} ckpc/h")
'''),
code(r'''
halos = il.groupcat.loadHalos(P.basePath, P.SNAP, fields=HALO_FIELDS)
'''),
code(r'''
sel       = (halos['Group_M_Crit200'] >= P.M200_MIN_INTERNAL) & (halos['GroupFirstSub'] >= 0)
group_idx = np.nonzero(sel)[0]          # índice FoF absoluto de TNG
n_cl      = len(group_idx)
'''),
code(r'''
M200c_h     = halos['Group_M_Crit200'][group_idx] * P.UM_H     # Msun/h
M200c       = halos['Group_M_Crit200'][group_idx] * P.UM       # Msun
R200c       = halos['Group_R_Crit200'][group_idx] * P.UL       # kpc físicos
GroupPos    = halos['GroupPos'][group_idx] * P.UL              # kpc (partícula más ligada)
bcg_sub_idx = halos['GroupFirstSub'][group_idx].astype(np.int64)
'''),
code(r'''
print(f"N halos              = {n_cl}   (paper: {P.PAPER_VALUES['N_halos']})")
print(f"N > 1e14 Msun/h      = {(M200c_h > 1e14).sum()}   (paper: {P.PAPER_VALUES['N_halos_1e14']})")
print(f"log M200c [Msun/h]   = {np.log10(M200c_h.min()):.2f} – {np.log10(M200c_h.max()):.2f}")
print(f"R200c [kpc]          = {R200c.min():.0f} – {R200c.max():.0f}")
'''),
code(r'''
subs = il.groupcat.loadSubhalos(P.basePath, P.SNAP, fields=SUB_FIELDS)

bcg_pos   = subs['SubhaloPos'][bcg_sub_idx] * P.UL              # kpc
bcg_Mstar = subs['SubhaloMassType'][bcg_sub_idx, 4] * P.UM       # Msun (BCG+ICL de Subfind)
bcg_flag  = subs['SubhaloFlag'][bcg_sub_idx]
'''),
code(r'''
# BCG = partícula más ligada del central: debe coincidir con GroupPos
d_pos = np.linalg.norm(F.periodic_delta(bcg_pos, GroupPos, box_kpc), axis=1)
print(f"max |SubhaloPos - GroupPos| = {d_pos.max():.2e} kpc")
print(f"BCGs con SubhaloFlag == 0   = {(bcg_flag == 0).sum()}")
'''),

md(r'''
## 4 · Masa estelar dentro de R₂₀₀ (M*₂₀₀)
`# SUPUESTO:` todas las estrellas del FoF (central, satélites y fuzz) con r₃D ≤ R₂₀₀c, sin partículas de viento.
Las estrellas de otros FoF que caen dentro de R₂₀₀ no se incluyen (efecto despreciable).
'''),
code(r'''
Mstar200 = np.full(n_cl, np.nan)

for i in range(n_cl):
    st = il.snapshot.loadHalo(P.basePath, P.SNAP, int(group_idx[i]), 'stars',
                              fields=['Coordinates', 'Masses', 'GFM_StellarFormationTime'])
    r3d  = np.linalg.norm(F.periodic_delta(st['Coordinates'] * P.UL, GroupPos[i], box_kpc), axis=1)
    keep = (r3d <= R200c[i]) & (st['GFM_StellarFormationTime'] > 0)
    Mstar200[i] = st['Masses'][keep].sum() * P.UM
    if (i + 1) % 10 == 0:
        print(f"  {i+1}/{n_cl}", end="\r")
print(f"\nlog M*200 [Msun] = {np.log10(np.nanmin(Mstar200)):.2f} – {np.log10(np.nanmax(Mstar200)):.2f}")
'''),

md(r'''
## 5 · Estado dinámico: offset 3D entre el BCG y el centro de masa del gas
`# SUPUESTO:` CM de las celdas de gas del FoF con r₃D ≤ `P.GAS_CM_R_FRAC`·R₂₀₀c (el paper no dice qué región usa).
'''),
code(r'''
gas_offset_kpc = np.full(n_cl, np.nan)

for i in range(n_cl):
    gas = il.snapshot.loadHalo(P.basePath, P.SNAP, int(group_idx[i]), 'gas',
                               fields=['Coordinates', 'Masses'])
    d   = F.periodic_delta(gas['Coordinates'] * P.UL, GroupPos[i], box_kpc)
    sel_gas = np.linalg.norm(d, axis=1) <= P.GAS_CM_R_FRAC * R200c[i]
    cm_rel  = F.center_of_mass(d[sel_gas], gas['Masses'][sel_gas])
    gas_offset_kpc[i] = np.linalg.norm(cm_rel)
    if (i + 1) % 10 == 0:
        print(f"  {i+1}/{n_cl}", end="\r")

dyn_offset = gas_offset_kpc / R200c
dyn_state  = F.dynamical_state(dyn_offset, P.DYN_RELAXED, P.DYN_DISTURBED)
'''),
code(r'''
for s in [0, 1, 2]:
    print(f"{P.DYN_LABELS[s]:<11} = {(dyn_state == s).sum():3d}")
print(f"D_BCG-CM/R200 mediana = {np.nanmedian(dyn_offset):.3f}")
'''),

md(r'''
## 6 · Concentración NFW (c = R₂₀₀c / r_s)
`# SUPUESTO:` partículas DM del FoF, `P.NFW_N_BINS` cáscaras logarítmicas entre `P.NFW_R_MIN_FRAC` y 1 R₂₀₀c, ajuste en log ρ.
'''),
code(r'''
# MassTable está en el Header del SNAPSHOT, no en el del groupcat (il.groupcat.loadHeader)
m_dm = F.dm_particle_mass_msun(P.basePath, P.SNAP, P.UM, P.M_DM_MSUN_FALLBACK)   # Msun
print(f"m_DM = {m_dm:.3e} Msun")

concentration = np.full(n_cl, np.nan)
r_s_kpc       = np.full(n_cl, np.nan)

for i in range(n_cl):
    dm  = il.snapshot.loadHalo(P.basePath, P.SNAP, int(group_idx[i]), 'dm', fields=['Coordinates'])
    rdm = np.linalg.norm(F.periodic_delta(dm * P.UL, GroupPos[i], box_kpc), axis=1)
    concentration[i], r_s_kpc[i], _, _ = F.fit_nfw_concentration(
        rdm, m_dm, R200c[i], P.NFW_R_MIN_FRAC, P.NFW_R_MAX_FRAC, P.NFW_N_BINS)
    if (i + 1) % 10 == 0:
        print(f"  {i+1}/{n_cl}  c_med = {np.nanmedian(concentration[:i+1]):.2f}", end="\r")
print(f"\nc mediana = {np.nanmedian(concentration):.2f}   (ajustes fallidos: {np.isnan(concentration).sum()})")
'''),

md(r'''
## 7 · Lookback time del último major merger del BCG
`# SUPUESTO:` razón con la masa estelar máxima pasada de cada progenitor (como `il.sublink.numMergers`);
el tiempo es el del snapshot descendiente en que ambos ya están fusionados.
'''),
code(r'''
z_snap = {}
for s in range(P.SNAP + 1):
    try:
        z_snap[s] = il.groupcat.loadHeader(P.basePath, s)['Redshift']
    except (OSError, KeyError):
        z_snap[s] = np.nan
'''),
code(r'''
snap_last_merger  = np.full(n_cl, -1, dtype=int)
ratio_last_merger = np.full(n_cl, np.nan)
t_last_merger     = np.full(n_cl, np.nan)      # Gyr (lookback)

for i in range(n_cl):
    tree = il.sublink.loadTree(P.basePath, P.SNAP, int(bcg_sub_idx[i]), fields=TREE_FIELDS)
    if tree is None:
        continue
    snap_last_merger[i], ratio_last_merger[i] = F.last_major_merger_snap(
        tree, P.MAJOR_MERGER_RATIO, mass_field=P.MERGER_MASS_FIELD)
    if snap_last_merger[i] >= 0:
        t_last_merger[i] = cosmo.lookback_time(z_snap[snap_last_merger[i]]).value
    if (i + 1) % 10 == 0:
        print(f"  {i+1}/{n_cl}", end="\r")
print(f"\nCon major merger: {np.isfinite(t_last_merger).sum()}/{n_cl}   "
      f"t_med = {np.nanmedian(t_last_merger):.2f} Gyr")
'''),

md('## 8 · Guardar el catálogo'),
code(r'''
with h5py.File(P.CATALOG_OUT, 'w') as f:
    f.attrs['basePath']   = P.basePath
    f.attrs['snap']       = P.SNAP
    f.attrs['n_clusters'] = n_cl
    f.attrs['seleccion']  = f"Group_M_Crit200 >= {P.M200_MIN_INTERNAL:g} [1e10 Msun/h] (M200c >= 1e13 Msun/h)"

    f.create_dataset('group_idx',          data=group_idx)          # índice FoF TNG
    f.create_dataset('M200c_Msun',         data=M200c)              # Msun
    f.create_dataset('M200c_Msun_h',       data=M200c_h)            # Msun/h
    f.create_dataset('R200c_kpc',          data=R200c)              # kpc físicos
    f.create_dataset('GroupPos_kpc',       data=GroupPos)           # kpc
    f.create_dataset('bcg_sub_idx',        data=bcg_sub_idx)        # subfind ID del central
    f.create_dataset('bcg_pos_kpc',        data=bcg_pos)            # kpc
    f.create_dataset('bcg_Mstar_Msun',     data=bcg_Mstar)          # Msun (todo el central)
    f.create_dataset('Mstar200_Msun',      data=Mstar200)           # Msun, r3D <= R200
    f.create_dataset('gas_offset_kpc',     data=gas_offset_kpc)     # kpc
    f.create_dataset('dyn_offset_R200',    data=dyn_offset)         # D_BCG-CM / R200
    f.create_dataset('dyn_state',          data=dyn_state)          # 0 rel, 1 int, 2 pert
    f.create_dataset('concentration',      data=concentration)      # R200c / r_s
    f.create_dataset('r_s_kpc',            data=r_s_kpc)            # kpc
    f.create_dataset('snap_last_merger',   data=snap_last_merger)   # -1 = sin major merger
    f.create_dataset('ratio_last_merger',  data=ratio_last_merger)
    f.create_dataset('t_last_merger_Gyr',  data=t_last_merger)      # Gyr lookback

print(f"Catálogo guardado en: {P.CATALOG_OUT}  ({n_cl} halos)")
'''),

md('## 9 · Figuras de control'),
code(MASS_LABEL_CELL),
code(r'''
fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))

ax[0].hist(lM200, bins=20, color=my_colors[0], edgecolor='white')
ax[0].axvline(14, color='gray', ls='--', lw=1)
ax[0].set_xlabel(rf'$\log_{{10}}(M_{{200c}}\,/\,{MASS_UNIT})$')
ax[0].set_ylabel('N halos')
ax[0].set_title(f'Muestra: {n_cl} halos')

ax[1].scatter(lM200, np.log10(Mstar200), s=20, color=my_colors[2], edgecolor='k', linewidth=0.3)
ax[1].set_xlabel(rf'$\log_{{10}}(M_{{200c}}\,/\,{MASS_UNIT})$')
ax[1].set_ylabel(r'$\log_{10}(M^{*}_{200}\,/\,M_\odot)$')

plt.savefig(f'{FIG_PDF}/fig00_muestra.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig00_muestra.png', dpi=150, bbox_inches='tight')
plt.show()
'''),
code(r'''
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))

ax[0].hist(dyn_offset, bins=25, color=my_colors[0], edgecolor='white')
for thr in [P.DYN_RELAXED, P.DYN_DISTURBED]:
    ax[0].axvline(thr, color='gray', ls='--', lw=1)
ax[0].set_xlabel(r'$D_{\rm BCG-CM}\,/\,R_{200}$')
ax[0].set_ylabel('N halos')

ax[1].scatter(lM200, concentration, s=20, color=my_colors[2], edgecolor='k', linewidth=0.3)
ax[1].set_xlabel(rf'$\log_{{10}}(M_{{200c}}\,/\,{MASS_UNIT})$')
ax[1].set_ylabel(r'$c_{200}$')

ax[2].hist(t_last_merger[np.isfinite(t_last_merger)], bins=20, color=my_colors[1], edgecolor='white')
ax[2].set_xlabel('Lookback time último major merger [Gyr]')

plt.savefig(f'{FIG_PDF}/fig00_propiedades_halo.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig00_propiedades_halo.png', dpi=150, bbox_inches='tight')
plt.show()
'''),
]

# illustris_python se importa en la celda de imports
cells[2] = code(cells[2]['source'].replace("import params_icl as P", "import illustris_python as il\nimport params_icl as P"))
if __name__ == "__main__":
    write_nb('00_catalogo.ipynb', cells)
