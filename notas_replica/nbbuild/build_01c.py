"""Edición mínima de 01c (exploratorio): funciones desde icl_funcs, banda r, rutas robustas."""
import nbformat
from common import *

PATH = os.path.join(REPO, '01c_perfil_subhalos_grupo.ipynb')
nb = nbformat.read(PATH, as_version=4)
cells = nb.cells


def set_src(c, src):
    c['source'] = src.strip('\n')
    if c['cell_type'] == 'code':
        c['outputs'] = []
        c['execution_count'] = None


# 0 · aviso en el título
c0 = cells[0]
assert c0['source'].startswith('# 01c')
set_src(c0, c0['source'].replace(
    '# 01c · Perfiles de Holmberg por subhalo — análisis de un grupo',
    '# 01c · Perfiles de Holmberg por subhalo — análisis de un grupo\n\n'
    '> **Exploratorio, fuera de la réplica de Mayes+2026.** Usa las mismas funciones de perfil que la réplica '
    '(`icl_funcs.py`: banda r, kernel corregido o anillos según `params_icl.py`), pero el centrado, el radio '
    'máximo (percentil 99) y la selección de subhalos son propios de este notebook.'))

# 1 · rutas + imports
c1 = cells[1]
assert 'import params_icl as P' in c1['source']
set_src(c1, REPO_CELL + r'''

import numpy as np
import h5py
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from astropy.cosmology import FlatLambdaCDM
from IPython.display import display

import illustris_python as il
import Catalogue
import params_icl as P
import icl_funcs as F

FIG_PDF = './figuras/01c_perfil_subhalos_grupo/pdf'
FIG_PNG = './figuras/01c_perfil_subhalos_grupo/png'
os.makedirs(FIG_PDF, exist_ok=True)
os.makedirs(FIG_PNG, exist_ok=True)

%matplotlib inline
plt.rcParams.update({'figure.dpi': 200, 'font.size': 12,
                     'axes.spines.top': False, 'axes.spines.right': False})

cosmo = FlatLambdaCDM(H0=P.H0, Om0=P.Om0, Ob0=P.Ob0)
''')

# 6/7 · funciones: envoltorios finos sobre icl_funcs (mismos nombres que usaba el notebook)
c6, c7 = cells[6], cells[7]
assert c7['source'].startswith('def rotate_by_inertia_tensor')
set_src(c6, '## Funciones auxiliares (envoltorios de `icl_funcs.py`, las mismas de la réplica)')
set_src(c7, r'''
rotate_by_inertia_tensor = F.rotate_by_inertia_tensor


def sb_profile_smooth(r_2d, lum_r, r_max_kpc, n_bins=60):
    """Perfil μ_r(R) [mag/arcsec²] hasta r_max_kpc con el método de params_icl (anillos o kernel corregido)."""
    return F.sb_profile(r_2d, lum_r, max(r_max_kpc, 1.0), P.SB_CONST, method=P.PROFILE_METHOD,
                        binning='log', n_bins=n_bins, r_min=P.SB_R_MIN_KPC, frac_sigma=P.KERNEL_FRAC_SIGMA)


def holmberg_radius(r_mid, mu_r, mu_cut=P.MU_HOLMBERG):
    """R_H [kpc] con la misma función de la réplica (método P.RH_METHOD)."""
    return F.holmberg_radius(r_mid, mu_r, mu_cut, method=P.RH_METHOD, r_fit_max=r_mid.max())['R_H']
''')

# 9 · banda r en vez de B
c9 = cells[9]
s = c9['source']
old = ("        # índice 1 = B-band (Johnson, Vega) en GFM_StellarPhotometrics (U,B,V,K,g,r,i,z)\n"
       "        lum_b = 10**((P.M_SUN_B_VEGA - phot[:, 1]) / 2.5)\n")
assert old in s
s = s.replace(old, "        # banda r (SDSS, AB): índice P.BAND_IDX en GFM_StellarPhotometrics (U,B,V,K,g,r,i,z)\n"
                   "        lum_b = F.luminosity_from_mag(phot[:, P.BAND_IDX], P.M_SUN_R_AB)\n")
s = s.replace("fields_star = ['Coordinates', 'Masses', 'GFM_StellarPhotometrics']",
              "fields_star = ['Coordinates', 'Masses', 'GFM_StellarPhotometrics']   # (viento no excluido: exploratorio)")
s = s.replace("        # r_mid, mu_B = sb_profile_r(r_2d, lum_b, r_max)\n", "")
set_src(c9, s)

# 13, 15, 17 · etiquetas μ_B → μ_r ; quitar sb_profile_r comentado
for k in (13, 15, 17):
    s = cells[k]['source'].replace(r"r'$\mu_B$", r"r'$\mu_r$").replace(
        "    # r_mid_icl, mu_icl = sb_profile_r(r_2d_out, outer_lum, r_max_out, n_bins=50)\n", "")
    set_src(cells[k], s)

nbformat.validate(nb)
nbformat.write(nb, PATH)
print(f'escrito {PATH}  ({len(nb.cells)} celdas)')
