"""Edición parcial de 01f: conserva intactas las celdas de acceso a la API (3, 4, 5, 6) y metadatos (9, 10)."""
import json, hashlib
import nbformat
from common import *

PATH = os.path.join(REPO, '01f_api_tng_BCG_ICL_R200.ipynb')
nb = nbformat.read(PATH, as_version=4)
cells = nb.cells
by_id = {c['id']: c for c in cells}
KEEP = ['dfc6a91f-1028-4674-bd2a-7fba25275deb', '1a6de62f-ad61-4a56-8098-75ace3ce9ec0',
        '047aa74e-a3e1-4b8c-9950-9765d72bc720', '109faf5f-4fab-4279-a6b0-7aa30e0ea480',
        'cbf2101a-4087-4abd-9f72-860bbbe34d11', '37cd25fc-6e5e-4f47-8db6-5cedbd04b773',
        '58e54322-c585-4a09-af0c-14bd4158960a', '3f58f8d6-5b6d-4cbc-bc92-aa77309faf50',
        '1f2d85dd-cdc3-4340-b653-7ea0cacd85b7', 'aa7b155b-b46e-4593-bec0-25c072d6d1d8',
        '749db988-fab8-417a-ae2f-0d726cdf5c44']
digest_before = {k: hashlib.md5(json.dumps(by_id[k], sort_keys=True).encode()).hexdigest() for k in KEEP}


def replace(cid, new_src):
    c = by_id[cid]
    c['source'] = new_src.strip('\n')
    if c['cell_type'] == 'code':
        c['outputs'] = []
        c['execution_count'] = None


replace('6253175f-a64e-481c-a0fe-0c888abc640b', r'''
# 01f · Separación BCG/ICL con un cutout de la API de TNG (un cúmulo)
**Objetivo:** verificar con datos descargados de la API que el pipeline de `01_` (mismas funciones
de `icl_funcs.py`: banda r, perfil hasta R₂₀₀c, ajuste logarítmico, R_H) da el mismo R_H y f_ICL
para un cúmulo.

- Se descarga el cutout del **halo FoF** (`halos/{FOF_IDX}/cutout.hdf5`) y se seleccionan partículas con
  la máscara por longitudes de subhalo: por defecto **sólo el subhalo central** (definición del paper,
  igual que `01_`); con `INCLUDE_FUZZ = True` se reproduce el chequeo de `01e_`.
- Corte esférico r₃D ≤ R₂₀₀c y sin partículas de viento.
- Las longitudes de subhalo (`SubhaloLenType`) se leen de los archivos locales del Lab.

Requiere haber corrido `00` y `01`.
''')

replace('4ac62eca-96ec-44db-abcd-9075f06faceb', REPO_CELL)
imports = code(r'''
from astropy.cosmology import FlatLambdaCDM

import time
import numpy as np
import h5py
import requests
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import binned_statistic_2d

import illustris_python as il
import Catalogue
import params_icl as P
import icl_funcs as F
''')
style = code(STYLE_CELL + '\n' + fig_dirs_cell('01f_api_tng_BCG_ICL_R200') + '\n' + COSMO_CELL)
replace('acaf302a-ec1f-4f37-8ef2-59d1a6d73a6a', style['source'])

replace('21ef3692-4600-4879-9bac-3c8d8b3cf646', r'''
with h5py.File(P.CATALOG_OUT, 'r') as f:
    group_idx_cat = f['group_idx'][:]
    M200c_all     = f['M200c_Msun'][:]
    R200c_all     = f['R200c_kpc'][:]
    GroupPos_all  = f['GroupPos_kpc'][:]
    bcg_sub_all   = f['bcg_sub_idx'][:]
    Mstar200_all  = f['Mstar200_Msun'][:]
    R_H_01_all    = f['R_H_kpc'][:] if 'R_H_kpc' in f else None
    f_icl_01_all  = f['f_icl'][:] if 'f_icl' in f else None
    R_H_fof_all   = f['R_H_kpc_fof'][:] if 'R_H_kpc_fof' in f else None
    f_icl_fof_all = f['f_icl_fof'][:] if 'f_icl_fof' in f else None

FOF_IDX  = int(group_idx_cat[GROUP_IDX])
M200c    = M200c_all[GROUP_IDX]
R200c    = R200c_all[GROUP_IDX]
GroupPos = GroupPos_all[GROUP_IDX]
BCG_sub  = int(bcg_sub_all[GROUP_IDX])
Mstar200 = Mstar200_all[GROUP_IDX]

print(f"GROUP_IDX = {GROUP_IDX}   FOF_IDX = {FOF_IDX}   BCG = {BCG_sub}")
print(f"log M200c = {np.log10(M200c):.2f} Msun   R200c = {R200c:.1f} kpc   log M*200 = {np.log10(Mstar200):.2f}")
''')

# Cutout: se agrega GFM_StellarFormationTime para quitar viento (resto igual)
c12 = by_id['15b90ded-9a00-4f85-9195-a2a7ad81850e']
src = c12['source']
src = src.replace("cutout_fields = 'Coordinates,Masses,GFM_StellarPhotometrics'",
                  "cutout_fields = 'Coordinates,Masses,GFM_StellarPhotometrics,GFM_StellarFormationTime'")
src = src.replace("    phot_raw     = f['PartType4']['GFM_StellarPhotometrics'][:]  # mag Vega/AB\n",
                  "    phot_raw     = f['PartType4']['GFM_StellarPhotometrics'][:]  # mag Vega/AB\n"
                  "    sft_raw      = (f['PartType4']['GFM_StellarFormationTime'][:]\n"
                  "                    if 'GFM_StellarFormationTime' in f['PartType4'] else None)\n")
assert 'GFM_StellarFormationTime' in src and 'sft_raw' in src, 'no se pudo editar la celda del cutout'
replace(c12['id'], src)

replace('ae0d3aa9-f77b-484d-aa81-f63a99577a28', r'''
## Máscara por longitudes de subhalo (`F.mask_from_lens`, la misma de `01e_`)
Las partículas del cutout FoF vienen ordenadas `[subhalo 0] [subhalo 1] ... [subhalo N-1] [no ligadas]`.
''')
replace('f5c14da9-9c9b-4158-b814-8141cce9f6f6', r'''
INCLUDE_SATELLITES = False   # satélites fuera (paper)
INCLUDE_FUZZ       = False   # False = sólo subhalo central (= 01_);  True = central + fuzz (= 01e_)

mask_sub = F.mask_from_lens(lens_star, n_fof_total, INCLUDE_SATELLITES, INCLUDE_FUZZ)

print(f"Total FoF       : {n_fof_total:>9,}")
print(f"Subhalo central : {int(lens_star[0]):>9,}")
print(f"Satélites       : {int(lens_star[1:].sum()):>9,}  ({'incluidos' if INCLUDE_SATELLITES else 'excluidos'})")
print(f"No ligadas      : {n_fof_total - int(lens_star.sum()):>9,}  ({'incluidas' if INCLUDE_FUZZ else 'excluidas'})")
print(f"Seleccionadas   : {mask_sub.sum():>9,}")
''')
replace('557b3353-0c52-41d8-a5e1-1989952846ce', '## Corte esférico r₃D ≤ R₂₀₀c y partículas de viento')
replace('9709fc07-a289-43e1-ab5b-414527aa8b46', r'''
r_3d      = np.linalg.norm(pos_c, axis=1)
mask_r200 = r_3d <= P.R_OUTER_FRAC * R200c
if sft_raw is not None:
    mask_star = sft_raw > 0
else:
    mask_star = np.ones(n_fof_total, dtype=bool)
    print("[AVISO] el cutout no trae GFM_StellarFormationTime: vuelve a descargarlo con overwrite=True")

mask_final = mask_sub & mask_r200 & mask_star
stars_api  = {'pos': pos_c[mask_final], 'mass': mass_raw[mask_final], 'phot': phot_raw[mask_final]}

print(f"Fuera de R200    : {(mask_sub & ~mask_r200).sum():>9,}")
print(f"Viento           : {(mask_sub & mask_r200 & ~mask_star).sum():>9,}")
print(f"Usadas (01f_)    : {mask_final.sum():>9,}")
''')
replace('7b51a3b3-4a90-4f36-94ff-86cdee088b6f', r'''
## Funciones
Rotación, perfil μ_r(R), ajuste logarítmico y R_H vienen de `icl_funcs.py` (las mismas de `01_`, `01e_` y `03_`).
''')

replace('93702a05-ad31-4651-aa16-9a304eace306', r'''
res_api = F.bcg_icl_from_stars(stars_api, Mstar200, R200c, P)
pos_rot = res_api['pos_rot']

print(f"R_H (ajuste)    : {res_api['R_H_fit']:.1f} kpc    R_H (interp): {res_api['R_H_interp']:.1f} kpc")
print(f"M_BCG           : {res_api['M_bcg']:.3e} Msun")
print(f"M_ICL           : {res_api['M_icl']:.3e} Msun")
print(f"f_ICL (M*200)   : {res_api['f_icl']:.4f}")
print(f"f_ICL (BCG+ICL) : {res_api['f_icl_bcgicl']:.4f}")

ref = (R_H_fof_all, f_icl_fof_all, '01e_') if INCLUDE_FUZZ else (R_H_01_all, f_icl_01_all, '01_')
if ref[0] is not None:
    print(f"\nComparación con {ref[2]} (mismo cúmulo, archivos locales):")
    print(f"  R_H   : {ref[0][GROUP_IDX]:.1f} kpc  vs API {res_api['R_H']:.1f} kpc")
    print(f"  f_ICL : {ref[1][GROUP_IDX]:.4f}      vs API {res_api['f_icl']:.4f}")
    print("  (deben coincidir si el cutout de la API = loadHalo/loadSubhalo local)")
else:
    print(f"\nNo hay resultados de {ref[2]} en el catálogo: córrelo primero")
''')

replace('d03344d7-6b4a-42a7-ab12-d66f523df573', r'''
rr = np.linspace(1, R200c, 500)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(res_api['r_mid'], res_api['mu'], 'o', ms=4, color=my_colors[0], label='Perfil 1D (API)')
ax.plot(rr, res_api['a'] + res_api['b'] * np.log10(rr), 'k-', lw=1.5, label='Ajuste logarítmico')
ax.axhline(P.MU_HOLMBERG, color='gray', ls='--', lw=1)
ax.axvline(res_api['R_H'], color='r', lw=1.5, label=f"R_H = {res_api['R_H']:.0f} kpc")
ax.set_ylim(33, 15)
ax.set_xlabel('Radio [kpc]')
ax.set_ylabel(r'$\mu_r$ [mag arcsec$^{-2}$]')
ax.set_title(f'FoF {FOF_IDX}  (cutout API, snap {SNAP})')
ax.legend()

plt.savefig(f'{FIG_PDF}/fig01f_perfil_sb_fof{FOF_IDX}_snap{SNAP}.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig01f_perfil_sb_fof{FOF_IDX}_snap{SNAP}.png', dpi=150, bbox_inches='tight')
plt.show()
''')

replace('161fa16d-dac3-41d4-b272-c5c2aa31d452', r'''
r_plot  = R200c
n_pix   = 512
edges   = np.linspace(-r_plot, r_plot, n_pix + 1)
pix_pc2 = ((2 * r_plot / n_pix) * 1e3)**2

L_map, _, _, _ = binned_statistic_2d(pos_rot[:, 0], pos_rot[:, 1], res_api['lum'],
                                     statistic='sum', bins=[edges, edges])
mu_map = F.mu_from_sigma(L_map / pix_pc2, P.SB_CONST)

fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(mu_map.T, origin='lower', cmap='ocean', vmin=16, vmax=32,
               extent=[-r_plot, r_plot, -r_plot, r_plot])
plt.colorbar(im, ax=ax, label=r'$\mu_r$ [mag arcsec$^{-2}$]')
ax.add_patch(plt.Circle((0, 0), res_api['R_H'], fill=False, color='k', lw=1.5, label='R_H'))
ax.add_patch(plt.Circle((0, 0), R200c, fill=False, color='gray', lw=1, ls='--', label='R$_{200}$'))
ax.set_xlabel('x [kpc]')
ax.set_ylabel('y [kpc]')
ax.grid(False)
ax.legend(loc='upper right')

plt.savefig(f'{FIG_PDF}/fig01f_mapa_sb_fof{FOF_IDX}_snap{SNAP}.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig01f_mapa_sb_fof{FOF_IDX}_snap{SNAP}.png', dpi=150, bbox_inches='tight')
plt.show()
''')

# ensamblar: insertar imports tras la celda de repo y quitar la celda de funciones duplicadas
new_cells = []
for c in cells:
    if c['id'] == '25c615d7-e828-42be-9065-c11a1295e953':
        continue
    new_cells.append(c)
    if c['id'] == '4ac62eca-96ec-44db-abcd-9075f06faceb':
        new_cells.append(imports)
nb.cells = new_cells

digest_after = {k: hashlib.md5(json.dumps(by_id[k], sort_keys=True).encode()).hexdigest() for k in KEEP}
assert digest_before == digest_after, 'se modificó una celda que debía quedar intacta'
nbformat.validate(nb)
nbformat.write(nb, PATH)
print(f'escrito {PATH}  ({len(nb.cells)} celdas; {len(KEEP)} celdas intactas verificadas por hash)')
