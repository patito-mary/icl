"""
Genera los 4 Jupyter Notebooks que replican la metodología de Mayes+2026:
"Coevolution of intracluster light and brightest cluster galaxies"

Versión adaptada para:
  - Datos TNG-100 en el cluster (sin API web)
  - Catálogo de grupos en formato pickle
  - illustris_python (il) para acceso a snapshots
  - Reutiliza Catalogue.py (sin modificarlo)
  - Sin cálculo de estado dinámico (se usa código externo)
"""
import nbformat as nbf

def code(src): return nbf.v4.new_code_cell(src)
def md(src):   return nbf.v4.new_markdown_cell(src)

def save(nb, name):
    with open(name, 'w') as f:
        nbf.write(nb, f)
    print(f"  → {name}")

# ═════════════════════════════════════════════════════════════════════════════
# NOTEBOOK 0 – Construcción del catálogo ICL (§2.1)
# ═════════════════════════════════════════════════════════════════════════════
nb0 = nbf.v4.new_notebook()
nb0.cells = [

md("""# 00 · Construcción del Catálogo de Cúmulos
**Metodología:** §2.1 de Mayes+2026

A partir del catálogo de group IDs (generado por tu código de estado dinámico),
cargamos las propiedades necesarias de TNG-100 usando `illustris_python`
y construimos el catálogo base para el análisis ICL/BCG.

**Campos que se extraen:**
- M₂₀₀c, R₂₀₀c, posición del halo
- ID del BCG (GroupFirstSub)
- Tiempo desde el último merger mayor del BCG (árbol SubLink)

**Lo que NO se hace aquí:** estado dinámico (lo calcula tu código existente).
El estado dinámico puede añadirse al catálogo como columna adicional.
"""),

code("""import sys, os, pickle
import numpy as np
import h5py
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

# Agregar el código original al path (solo lectura, no se modifica)
sys.path.insert(0, './original_shift_code')
import illustris_python as il
import params_icl as P          # copia local de params para este análisis

%matplotlib inline
plt.rcParams.update({'figure.dpi': 110, 'font.size': 12,
                     'axes.spines.top': False, 'axes.spines.right': False})

cosmo = FlatLambdaCDM(H0=67.74, Om0=0.3089)
print(f"basePath : {P.basePath}")
print(f"Snapshot : {P.SNAP}")
print(f"h = {P.h},  UL = {P.UL:.4f},  UM = {P.UM:.3e}")
"""),

md("## 1 · Leer el catálogo de group IDs"),

code("""# ── Cargar el pickle con los grupos a analizar ────────────────────────────
with open(P.CATALOG_PKL, 'rb') as f:
    cat_raw = pickle.load(f)

# Inspeccionar la estructura
print("Tipo:", type(cat_raw))
if isinstance(cat_raw, dict):
    print("Claves:", list(cat_raw.keys()))
    # Muestra las primeras entradas de cada campo
    for k, v in cat_raw.items():
        try: print(f"  {k}: shape={np.shape(v)}, primeros valores={np.asarray(v).flat[:3]}")
        except: print(f"  {k}: {v}")
else:
    print("Shape:", np.shape(cat_raw))
    print("Primeros valores:", np.asarray(cat_raw)[:5])
"""),

code("""# ── Extraer los IDs de grupos (ajustar según la estructura del pkl) ────────
#
# CASO A: el pkl es un dict con 'GroupNumber' (salida de TNG2PROPnum)
#         Los GroupNumber tienen la forma snap*1e10 + grupo_idx
#         → grupo_idx = GroupNumber - 99*1e10
#
# CASO B: el pkl es simplemente un array de índices de grupo (enteros)
#
# Descomenta el caso que corresponda:

# ── CASO A ────────────────────────────────────────────────────────────────
# group_numbers = np.array(cat_raw['GroupNumber'])
# group_idx     = (group_numbers - P.SNAP * 1e10).astype(int)

# ── CASO B ────────────────────────────────────────────────────────────────
group_idx = np.array(cat_raw, dtype=int)   # ajustar key/field si es dict

n_cl = len(group_idx)
print(f"N cúmulos a analizar : {n_cl}")
print(f"Rango de índices     : {group_idx.min()} – {group_idx.max()}")
"""),

md("## 2 · Cargar propiedades del catálogo de halos TNG"),

code("""# Cargar header y catálogos completos de halos y subhalos
Header = il.groupcat.loadHeader(P.basePath, P.SNAP)
print("BoxSize (ckpc/h):", Header['BoxSize'])
print("Redshift        :", Header['Redshift'])
print("Time (a)        :", Header['Time'])

# Catálogo de halos (todos los grupos de TNG100-1)
Halos_all = il.groupcat.loadHalos(P.basePath, P.SNAP, fields=P.HALOS_FIELDS)
print(f"\\nTotal de halos en TNG100-1 snap {P.SNAP}: {Halos_all['count']}")
"""),

code("""# Extraer sólo los halos que están en nuestro catálogo
M200c_raw = Halos_all['Group_M_Crit200'][group_idx]   # 1e10 M☉/h
R200c_raw = Halos_all['Group_R_Crit200'][group_idx]   # ckpc/h
GroupPos  = Halos_all['GroupPos'][group_idx] * P.UL   # kpc físicos
GroupCM   = Halos_all['GroupCM'][group_idx]  * P.UL   # kpc físicos
first_sub = Halos_all['GroupFirstSub'][group_idx]     # índice del BCG

# Convertir a unidades físicas
M200c = M200c_raw * P.UM     # M☉
R200c = R200c_raw * P.UL     # kpc físicos

print(f"N cúmulos         : {n_cl}")
print(f"log M200c (M☉)    : {np.log10(M200c.min()):.2f} – {np.log10(M200c.max()):.2f}")
print(f"R200c [kpc]       : {R200c.min():.0f} – {R200c.max():.0f}")
print(f"Primeros GroupFirstSub: {first_sub[:5]}")
"""),

md("## 3 · Posición del BCG para cada cúmulo"),

code("""# Cargar posiciones de subhalos para identificar el BCG
Subhalos_all = il.groupcat.loadSubhalos(P.basePath, P.SNAP, fields=P.SUBHALOS_FIELDS)

bcg_pos    = Subhalos_all['SubhaloPos'][first_sub] * P.UL   # kpc físicos
bcg_M_star = Subhalos_all['SubhaloMassType'][first_sub, 4] * P.UM  # masa estelar BCG [M☉]

print("Posición BCG (primeros 3):")
for i in range(min(3, n_cl)):
    print(f"  Cúmulo {group_idx[i]}: pos = {bcg_pos[i]}")
print(f"\\nMasa estelar BCG: {np.log10(bcg_M_star.min()):.2f} – {np.log10(bcg_M_star.max()):.2f} log M☉")
"""),

md("""## 4 · Tiempo desde el último merger mayor del BCG

Usamos el árbol SubLink de TNG para trazar el Main Progenitor Branch (MPB) del BCG
y encontrar el último evento de fusión con razón de masa estelar ≥ 1/5.

Referencia: Mayes+2026 §2 + Rodríguez-Gómez+2015
"""),

code("""def time_since_last_major_merger(bcg_id, basePath, snap, cosmo,
                                    mass_ratio_min=0.2):
    \"\"\"
    Traza el árbol MPB del BCG y encuentra el último merger mayor.

    Parámetros
    ----------
    bcg_id        : int, índice del subhalo BCG en snap
    mass_ratio_min: float, razón mínima de masa estelar (defecto 1/5)

    Devuelve
    --------
    t_lookback : float, tiempo de lookback en Gyr desde z=0 hasta el merger.
                 NaN si no hay merger mayor registrado.
    \"\"\"
    fields = ['SnapNum', 'SubhaloMassType', 'SubhaloID',
              'NextProgenitorID', 'FirstProgenitorID', 'MainLeafProgenitorID']
    try:
        tree = il.sublink.loadTree(basePath, snap, bcg_id, fields=fields,
                                    onlyMPB=False)
    except Exception as e:
        return np.nan

    if tree is None or 'SnapNum' not in tree:
        return np.nan

    snap_arr   = tree['SnapNum']
    mass_arr   = tree['SubhaloMassType'][:, 4]   # masa estelar (tipo 4)

    # Masa del progenitor principal en cada snapshot
    # Los snapshots están ordenados de z=0 hacia atrás en el MPB
    # Para detectar mergers: buscar el snap donde la masa del 2do progenitor
    # es significativa respecto a la del progenitor principal.
    # Usando la masa del PROG2: M_main[snap] - M_main[snap+1] (masa que se sumó)
    # Aproximación: Δm = masa sumada en cada paso del árbol
    # La masa del 2do progenitor = masa del subtree que se une = Δm

    # Identificamos mergers como: M_main(t-1) salta más que factor respecto a M_main(t)
    # (el MPB tiene masa creciente hacia z=0, decrece hacia atrás)
    # Un merger mayor ocurre cuando la masa cae significativamente hacia atrás:
    # ratio = M_main[i+1] / M_main[i] < (1 - mass_ratio_min / (1 + mass_ratio_min))
    t_last = np.nan
    for i in range(len(snap_arr) - 1):
        m_now  = mass_arr[i]
        m_prev = mass_arr[i + 1]
        if m_now <= 0 or m_prev <= 0:
            continue
        # Masa del 2do progenitor ≈ masa ganada en este paso
        m2 = m_now - m_prev
        if m2 <= 0:
            continue
        ratio = m2 / m_prev   # ≈ M2/M1 en el momento del merger
        if ratio >= mass_ratio_min:
            snap_merger = snap_arr[i]
            z_merger    = il.groupcat.loadHeader(basePath, int(snap_merger))['Redshift']
            t_last      = cosmo.lookback_time(z_merger).value   # Gyr
            break   # primer (más reciente) merger mayor encontrado

    return t_last

print("Calculando tiempo desde el último merger mayor...")
t_last_merger = np.full(n_cl, np.nan)
for i, sub_id in enumerate(first_sub):
    t_last_merger[i] = time_since_last_major_merger(
        int(sub_id), P.basePath, P.SNAP, cosmo)
    if (i+1) % 10 == 0:
        print(f"  {i+1}/{n_cl}  t_med = {np.nanmedian(t_last_merger[:i+1]):.2f} Gyr", end="\\r")
print(f"\\nMerger mayor encontrado en {np.isfinite(t_last_merger).sum()}/{n_cl} cúmulos")
"""),

md("## 5 · Guardar el catálogo"),

code("""# Guardar como HDF5 para uso en los siguientes notebooks
with h5py.File(P.CATALOG_OUT, 'w') as f:
    f.attrs['basePath']    = P.basePath
    f.attrs['snap']        = P.SNAP
    f.attrs['n_clusters']  = n_cl

    f.create_dataset('group_idx',     data=group_idx)    # índice FoF en TNG
    f.create_dataset('M200c_Msun',    data=M200c)        # M☉
    f.create_dataset('R200c_kpc',     data=R200c)        # kpc físicos
    f.create_dataset('GroupPos_kpc',  data=GroupPos)     # kpc físicos (3D)
    f.create_dataset('GroupCM_kpc',   data=GroupCM)      # kpc físicos (3D)
    f.create_dataset('bcg_sub_idx',   data=first_sub)    # índice subhalo BCG
    f.create_dataset('bcg_pos_kpc',   data=bcg_pos)      # kpc físicos (3D)
    f.create_dataset('bcg_Mstar_Msun',data=bcg_M_star)   # M☉
    f.create_dataset('t_last_merger_Gyr', data=t_last_merger)  # Gyr

print(f"Catálogo guardado en: {P.CATALOG_OUT}")
print(f"  N cúmulos          : {n_cl}")
print(f"  log M200c [M☉]     : {np.log10(M200c.min()):.2f} – {np.log10(M200c.max()):.2f}")
"""),

md("## 6 · Visualización rápida del catálogo"),

code("""fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.hist(np.log10(M200c), bins=20, color='steelblue', edgecolor='white', lw=0.5)
ax.axvline(np.log10(M200c.min()), ls='--', color='k', lw=1.2, label='Límite mínimo')
ax.set_xlabel(r'$\\log_{10}(M_{200c}\\,/\\,M_\\odot)$')
ax.set_ylabel('N cúmulos')
ax.set_title('Distribución de masas')
ax.legend(fontsize=9)

ax = axes[1]
valid = np.isfinite(t_last_merger)
ax.scatter(np.log10(M200c[valid]), t_last_merger[valid],
           s=20, alpha=0.8, color='tomato')
ax.set_xlabel(r'$\\log_{10}(M_{200c}\\,/\\,M_\\odot)$')
ax.set_ylabel('Lookback time último merger mayor [Gyr]')
ax.set_title('Tiempo desde merger (Fig. 5 input)')

plt.tight_layout()
plt.savefig('fig_catalogo.pdf', bbox_inches='tight')
plt.show()

# Añadir estado dinámico al catálogo si ya lo tienes calculado:
# with h5py.File(P.CATALOG_OUT, 'a') as f:
#     f.create_dataset('dyn_state', data=mi_array_dyn_state)
#     # 0=relajado, 1=intermedio, 2=perturbado
"""),

]

# ═════════════════════════════════════════════════════════════════════════════
# NOTEBOOK 1 – Separación BCG/ICL y fracción de masa (§2.2, §3)
# ═════════════════════════════════════════════════════════════════════════════
nb1 = nbf.v4.new_notebook()
nb1.cells = [

md("""# 01 · Separación BCG/ICL y Fracción de Masa del ICL
**Metodología:** §2.2 y §3 de Mayes+2026

1. **Rotación por tensor de inercia reducido** — usando la misma construcción que `Catalogue.comp_Reduced_Inertia_tensor()`
2. **Perfil de brillo superficial 1D** en banda r a lo largo del eje semi-mayor
3. **Radio de Holmberg** — corte en μ_r = 26.5 mag/arcsec² (separación BCG/ICL)
4. **Fracción de masa del ICL** vs masa, concentración y tiempo del último merger

**Equivalente a:** Figs. 1, 2, 3, 4, 5 de Mayes+2026
"""),

code("""import sys, os, pickle
import numpy as np
import h5py
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d
from scipy.stats import linregress
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

# ── Importar código del proyecto (solo lectura) ──────────────────────────
sys.path.insert(0, './original_shift_code')
import illustris_python as il
import Catalogue                 # para Distance_3D, CM_3D, etc.
import params_icl as P

%matplotlib inline
plt.rcParams.update({'figure.dpi': 110, 'font.size': 12,
                     'axes.spines.top': False, 'axes.spines.right': False})

cosmo = FlatLambdaCDM(H0=67.74, Om0=0.3089)
"""),

md("## Cargar catálogo"),

code("""with h5py.File(P.CATALOG_OUT, 'r') as f:
    group_idx     = f['group_idx'][:]
    M200c         = f['M200c_Msun'][:]
    R200c         = f['R200c_kpc'][:]
    GroupPos      = f['GroupPos_kpc'][:]
    bcg_sub_idx   = f['bcg_sub_idx'][:]
    t_last_merger = f['t_last_merger_Gyr'][:]
    # Estado dinámico (si ya fue añadido)
    dyn_state = f['dyn_state'][:] if 'dyn_state' in f else None

n_cl = len(group_idx)
Header = il.groupcat.loadHeader(P.basePath, P.SNAP)

print(f"Catálogo cargado: {n_cl} cúmulos")
print(f"Estado dinámico disponible: {dyn_state is not None}")

# Colores por estado dinámico
COLORS_STATE  = {0: '#2196F3', 1: '#FF9800', 2: '#F44336'}
LABELS_STATE  = {0: 'Relajado', 1: 'Intermedio', 2: 'Perturbado'}
COLORS_COMP   = {0: '#E91E63', 1: '#9C27B0', 2: '#00BCD4'}
LABELS_COMP   = {0: 'In situ', 1: 'Mergers', 2: 'Stripped'}
"""),

md("## §2.2 – Rotación por tensor de inercia\n\nUsamos la misma construcción que `Catalogue.comp_Reduced_Inertia_tensor()` pero pedimos eigenvectores para rotar las partículas."),

code("""def rotate_by_inertia_tensor(pos_rel, mass, r_lim=np.inf):
    \"\"\"
    Alinea las partículas con los ejes principales del tensor de inercia reducido.
    Usa la misma construcción que Catalogue.comp_Reduced_Inertia_tensor().

    pos_rel : (N,3) posiciones relativas al BCG [kpc]
    mass    : (N,)  masas de las partículas
    r_lim   : float, radio límite para el cálculo del tensor [kpc]

    Devuelve
    --------
    pos_rot : (N,3) posiciones rotadas (eje mayor → x̂, menor → ẑ)
    R_mat   : (3,3) matriz de rotación
    \"\"\"
    dist = np.linalg.norm(pos_rel, axis=1)
    ok   = (dist > 0) & (dist <= r_lim) & np.isfinite(mass)
    p, m = pos_rel[ok], mass[ok]

    if m.sum() == 0 or len(m) < 4:
        return pos_rel, np.eye(3)

    # Tensor de inercia reducido (ponderado por 1/r²)
    w    = 1.0 / dist[ok]**2
    mtot = np.sum(m)
    Ixx  = np.sum(m * p[:,0]**2 * w) / mtot
    Iyy  = np.sum(m * p[:,1]**2 * w) / mtot
    Izz  = np.sum(m * p[:,2]**2 * w) / mtot
    Ixy  = np.sum(m * p[:,0] * p[:,1] * w) / mtot
    Ixz  = np.sum(m * p[:,0] * p[:,2] * w) / mtot
    Iyz  = np.sum(m * p[:,1] * p[:,2] * w) / mtot
    I = np.array([[Ixx, Ixy, Ixz],
                  [Ixy, Iyy, Iyz],
                  [Ixz, Iyz, Izz]])

    eigvals, eigvecs = np.linalg.eigh(I)
    # eigh devuelve valores en orden ascendente
    # eigenvalor mayor → eje de mayor extensión (semi-eje mayor) → queremos que sea x̂
    idx   = np.argsort(eigvals)[::-1]   # descendente: mayor primero
    R_mat = eigvecs[:, idx].T           # filas = ejes principales
    pos_rot = pos_rel @ R_mat.T
    return pos_rot, R_mat
"""),

md("## §2.2 – Perfil de brillo superficial en banda r y radio de Holmberg"),

code("""def sb_profile_r(r_2d, lum_r, r_max_kpc, n_bins=60):
    \"\"\"
    Perfil 1D de brillo superficial en banda r.

    r_2d    : distancias 2D proyectadas [kpc]
    lum_r   : luminosidad r de cada partícula [L☉]
                GFM_StellarPhotometrics[:, 5] = magnitud abs r-band (AB)
                L = 10^((M_sun,r - mag) / 2.5)
    r_max_kpc: radio máximo [kpc]

    Conversión:
        μ_r [mag/arcsec²] = M_☉,r + 21.572 − 2.5 log₁₀(Σ_L [L☉/pc²])
    \"\"\"
    r_bins = np.logspace(np.log10(0.5), np.log10(r_max_kpc), n_bins + 1)
    r_mid  = np.sqrt(r_bins[:-1] * r_bins[1:])   # media geométrica
    mu_r   = np.full(n_bins, np.nan)

    for k, (r1, r2) in enumerate(zip(r_bins[:-1], r_bins[1:])):
        mk = (r_2d >= r1) & (r_2d < r2)
        if mk.sum() == 0:
            continue
        area_pc2 = np.pi * ((r2 * 1e3)**2 - (r1 * 1e3)**2)   # kpc→pc, anillo
        sigma_L  = lum_r[mk].sum() / area_pc2
        if sigma_L > 0:
            mu_r[k] = P.SB_CONST - 2.5 * np.log10(sigma_L)

    return r_mid, mu_r


def holmberg_radius(r_mid, mu_r, mu_cut=P.MU_HOLMBERG):
    \"\"\"Interpola el radio donde μ_r = mu_cut.\"\"\"
    valid = np.isfinite(mu_r) & (r_mid > 0)
    if valid.sum() < 3:
        return np.nan
    # μ crece con r; necesitamos que sea monótonamente creciente para invertir
    r_v, m_v = r_mid[valid], mu_r[valid]
    # Ordenar por radio por si hay ruido
    idx_s = np.argsort(r_v)
    r_v, m_v = r_v[idx_s], m_v[idx_s]
    if m_v[0] > mu_cut or m_v[-1] < mu_cut:
        return np.nan
    try:
        f = interp1d(m_v, r_v, kind='linear', fill_value='extrapolate')
        r_h = float(f(mu_cut))
        return r_h if 0 < r_h <= r_v[-1] * 1.2 else np.nan
    except Exception:
        return np.nan
"""),

md("## §2.2 – Separación BCG/ICL: demostración con un cúmulo"),

code("""# Elegir el cúmulo más masivo como ejemplo
i_demo = np.argmax(M200c)
sub_id = int(bcg_sub_idx[i_demo])
cen    = GroupPos[i_demo]   # kpc físicos

print(f"Cúmulo demo: group_idx={group_idx[i_demo]}, "
      f"sub_id={sub_id}, log M200c={np.log10(M200c[i_demo]):.2f}")

# Cargar partículas estelares del subhalo central (BCG+ICL)
# Todas las estrellas no ligadas a satélites están en el subhalo central
fields = ['Coordinates', 'Masses', 'GFM_StellarPhotometrics',
          'GFM_Metallicity', 'GFM_StellarFormationTime']
stars = il.snapshot.loadSubhalo(P.basePath, P.SNAP, sub_id, 'stars', fields=fields)

# Convertir unidades
pos   = stars['Coordinates'] * P.UL               # kpc físicos
mass  = stars['Masses'] * P.UM                    # M☉
phot  = stars['GFM_StellarPhotometrics']           # mag abs AB (U,B,V,K,g,r,i,z)
metal = stars['GFM_Metallicity']                   # fracción metálica Z
aform = stars['GFM_StellarFormationTime']          # factor de escala de formación

print(f"N partículas estelares: {len(mass):,}")

# Centrar respecto al BCG usando Distance_3D (corrige periodicidad)
pos_c = Catalogue.Distance_3D(pos, cen, Header['BoxSize'] * P.UL)

# Rotar al plano del eje mayor (face-on)
pos_rot, R_mat = rotate_by_inertia_tensor(pos_c, mass)

# Proyección 2D (plano xy)
r_2d = np.sqrt(pos_rot[:, 0]**2 + pos_rot[:, 1]**2)

# Luminosidad r-band: índice 5 en GFM_StellarPhotometrics (U,B,V,K,g,r,i,z)
lum_r = 10**((P.M_SUN_R_AB - phot[:, 5]) / 2.5)   # L☉

# Perfil de brillo superficial
r_max = np.percentile(r_2d, 99)
r_mid, mu_r = sb_profile_r(r_2d, lum_r, r_max)

# Radio de Holmberg
r_h = holmberg_radius(r_mid, mu_r)
print(f"Radio de Holmberg: {r_h:.1f} kpc")

# Separación BCG/ICL
mask_bcg = r_2d <= r_h
mask_icl = r_2d >  r_h
M_bcg_demo = mass[mask_bcg].sum()
M_icl_demo = mass[mask_icl].sum()
print(f"f_ICL = {M_icl_demo/(M_bcg_demo+M_icl_demo):.3f}")
"""),

md("### Figura 1 – Perfil de brillo superficial 1D (= Fig. 1 del paper)"),

code("""fig, ax = plt.subplots(figsize=(7, 5))
valid = np.isfinite(mu_r)

ax.plot(r_mid[valid], mu_r[valid], 'k-', lw=2, label='BCG+ICL')
ax.axhline(P.MU_HOLMBERG, color='r', ls='--', lw=1.5,
           label=f'μ_r = {P.MU_HOLMBERG} mag/arcsec²')
if np.isfinite(r_h):
    ax.axvline(r_h, color='r', ls=':', lw=1.5,
               label=f'R_Holmberg = {r_h:.0f} kpc')
    ax.fill_betweenx([mu_r[valid].min()-1, P.MU_HOLMBERG+2],
                      0, r_h, alpha=0.10, color='royalblue', label='BCG')
    ax.fill_betweenx([mu_r[valid].min()-1, P.MU_HOLMBERG+2],
                      r_h, r_mid[valid][-1], alpha=0.10, color='tomato', label='ICL')

ax.set_xscale('log')
ax.invert_yaxis()
ax.set_xlabel('Radio semi-mayor [kpc]')
ax.set_ylabel(r'$\mu_r$ [mag arcsec$^{-2}$]')
ax.set_title(f'TNG100-1  snap={P.SNAP}  grupo {group_idx[i_demo]}')
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig('fig01_perfil_sb.pdf', bbox_inches='tight')
plt.show()
"""),

md("### Figura 2 – Mapa 2D de brillo superficial (= Fig. 2 del paper)"),

code("""from scipy.stats import binned_statistic_2d

r_plot  = R200c[i_demo]
n_pix   = 400
edges   = np.linspace(-r_plot, r_plot, n_pix + 1)
pix_pc2 = ((2 * r_plot / n_pix) * 1e3)**2

H, _, _, _ = binned_statistic_2d(pos_rot[:, 0], pos_rot[:, 1],
                                   lum_r, statistic='sum',
                                   bins=[edges, edges])
with np.errstate(divide='ignore', invalid='ignore'):
    sigma = np.where(H > 0, H / pix_pc2, np.nan)
    mu_map = np.where(sigma > 0, P.SB_CONST - 2.5 * np.log10(sigma), np.nan)

fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(mu_map.T, origin='lower', cmap='magma_r',
               extent=[-r_plot, r_plot, -r_plot, r_plot],
               vmin=20, vmax=30)
plt.colorbar(im, ax=ax, label=r'$\mu_r$ [mag arcsec$^{-2}$]')
if np.isfinite(r_h):
    ax.add_patch(plt.Circle((0, 0), r_h,   fill=False, color='black',
                              lw=1.5, label=f'R_H={r_h:.0f} kpc'))
ax.add_patch(plt.Circle((0, 0), r_plot, fill=False, color='white',
                          lw=1, ls='--', label='R₂₀₀'))
ax.set_xlabel('x [kpc]'); ax.set_ylabel('y [kpc]')
ax.set_title(f'Mapa 2D – grupo {group_idx[i_demo]}')
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig('fig02_mapa_sb.pdf', bbox_inches='tight')
plt.show()
"""),

md("## §3 – Fracción de masa ICL para todos los cúmulos"),

code("""def bcg_icl_mass_fraction(sub_id, cen_pos, header=Header):
    \"\"\"
    Carga partículas del subhalo central, aplica el corte de Holmberg
    y devuelve (f_ICL, M_BCG, M_ICL) en M☉.
    \"\"\"
    try:
        fields = ['Coordinates', 'Masses', 'GFM_StellarPhotometrics']
        st = il.snapshot.loadSubhalo(P.basePath, P.SNAP, int(sub_id), 'stars', fields=fields)
        pos  = Catalogue.Distance_3D(st['Coordinates'] * P.UL, cen_pos, header['BoxSize'] * P.UL)
        mass = st['Masses'] * P.UM
        lum  = 10**((P.M_SUN_R_AB - st['GFM_StellarPhotometrics'][:, 5]) / 2.5)
        pos_rot, _ = rotate_by_inertia_tensor(pos, mass)
        r2   = np.sqrt(pos_rot[:,0]**2 + pos_rot[:,1]**2)
        r_m, mu = sb_profile_r(r2, lum, np.percentile(r2, 99))
        r_h  = holmberg_radius(r_m, mu)
        if not np.isfinite(r_h): return np.nan, np.nan, np.nan
        m_bcg = mass[r2 <= r_h].sum()
        m_icl = mass[r2 >  r_h].sum()
        m_tot = m_bcg + m_icl
        return (m_icl / m_tot, m_bcg, m_icl) if m_tot > 0 else (np.nan, np.nan, np.nan)
    except Exception as e:
        return np.nan, np.nan, np.nan

icl_frac = np.full(n_cl, np.nan)
M_bcg    = np.full(n_cl, np.nan)
M_icl    = np.full(n_cl, np.nan)

print("Calculando fracciones de masa ICL...")
for i in range(n_cl):
    icl_frac[i], M_bcg[i], M_icl[i] = bcg_icl_mass_fraction(
        bcg_sub_idx[i], GroupPos[i])
    if (i+1) % 10 == 0:
        print(f"  {i+1}/{n_cl}  f_ICL_med={np.nanmedian(icl_frac[:i+1]):.3f}", end="\\r")

# Guardar en el catálogo
with h5py.File(P.CATALOG_OUT, 'a') as f:
    for k, v in [('icl_frac', icl_frac), ('M_bcg_Msun', M_bcg), ('M_icl_Msun', M_icl)]:
        if k in f: del f[k]
        f.create_dataset(k, data=v)
print(f"\\nf_ICL = {np.nanmean(icl_frac):.3f} ± {np.nanstd(icl_frac):.3f}")
"""),

md("### Figura 3 – Fracción de masa ICL vs M₂₀₀ (= Fig. 3 del paper)"),

code("""def linfit(x, y, log_x=True):
    xx = np.log10(x) if log_x else np.asarray(x)
    yy = np.asarray(y)
    ok = np.isfinite(xx) & np.isfinite(yy)
    sl, ic, *_ = linregress(xx[ok], yy[ok])
    return sl, ic, xx[ok], yy[ok]

lM = np.log10(M200c)
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
if dyn_state is not None:
    for s in [0, 1, 2]:
        m = dyn_state == s
        ax.scatter(lM[m], icl_frac[m], color=COLORS_STATE[s],
                   label=LABELS_STATE[s], s=20, alpha=0.75)
else:
    ax.scatter(lM, icl_frac, color='steelblue', s=20, alpha=0.75)

sl, ic, lx, ly = linfit(M200c, icl_frac)
xx = np.linspace(lx.min(), lx.max(), 100)
ax.plot(xx, sl*xx+ic, 'k-', lw=1.8, label=f'β = {sl:.3f}')
ax.set_xlabel(r'$\\log_{10}(M_{200c}\\,/\\,M_\\odot)$')
ax.set_ylabel(r'$f_{\\rm ICL}$')
ax.set_title('Fracción de masa ICL (Fig. 3)')
ax.legend(fontsize=9)

ax = axes[1]
ax.hist(icl_frac[np.isfinite(icl_frac)], bins=20,
        color='steelblue', edgecolor='white', lw=0.5)
ax.axvline(np.nanmedian(icl_frac), ls='--', color='k', lw=1.5,
           label=f'Mediana = {np.nanmedian(icl_frac):.3f}')
ax.set_xlabel(r'$f_{\\rm ICL}$')
ax.set_ylabel('N cúmulos')
ax.legend()

plt.tight_layout()
plt.savefig('fig03_icl_vs_masa.pdf', bbox_inches='tight')
plt.show()
"""),

md("### Figura 4 – Fracción ICL vs concentración NFW (= Fig. 4 del paper)\n\nCalculamos c₂₀₀ ajustando un perfil NFW al perfil de masa encerrada de la materia oscura."),

code("""def nfw_enclosed_norm(r, r_s):
    \"\"\"Masa encerrada NFW normalizada a r200.\"\"\"
    x200 = R200c[_i_nfw] / r_s
    x    = r / r_s
    return (np.log(1 + x) - x/(1+x)) / (np.log(1+x200) - x200/(1+x200))

def compute_concentration(sub_id_or_halo_idx, r200, cen_pos, header=Header):
    \"\"\"Ajusta perfil NFW a las partículas DM del halo.\"\"\"
    try:
        dm = il.snapshot.loadHalo(P.basePath, P.SNAP, int(sub_id_or_halo_idx),
                                   'dm', fields=['Coordinates'])
        if not isinstance(dm, np.ndarray): dm = dm
        pos_dm = Catalogue.Distance_3D(dm * P.UL, cen_pos, header['BoxSize'] * P.UL)
        r_dm   = np.linalg.norm(pos_dm, axis=1)
        r_bins = np.logspace(np.log10(1.0), np.log10(r200), 31)
        r_mid  = np.sqrt(r_bins[:-1] * r_bins[1:])
        M_enc  = np.array([(r_dm < rb).sum() for rb in r_bins[1:]], dtype=float)
        if M_enc[-1] == 0: return np.nan
        M_enc /= M_enc[-1]
        popt, _ = curve_fit(lambda r, rs: (np.log(1+r/rs) - (r/rs)/(1+r/rs)) /
                                           (np.log(1+r200/rs) - (r200/rs)/(1+r200/rs)),
                             r_mid, M_enc, p0=[r200/5], bounds=(0.1, r200))
        return r200 / popt[0]
    except Exception:
        return np.nan

print("Calculando concentración NFW para cada cúmulo...")
concentration = np.full(n_cl, np.nan)
for i in range(n_cl):
    _i_nfw = i
    concentration[i] = compute_concentration(group_idx[i], R200c[i], GroupPos[i])
    if (i+1) % 10 == 0:
        print(f"  {i+1}/{n_cl}  c_med={np.nanmedian(concentration[:i+1]):.2f}", end="\\r")

# Guardar en catálogo
with h5py.File(P.CATALOG_OUT, 'a') as f:
    if 'concentration' in f: del f['concentration']
    f.create_dataset('concentration', data=concentration)

fig, ax = plt.subplots(figsize=(7, 5))
valid = np.isfinite(concentration) & np.isfinite(icl_frac)
sc = ax.scatter(concentration[valid], icl_frac[valid],
                c=lM[valid], cmap='viridis', s=25, alpha=0.8)
plt.colorbar(sc, ax=ax, label=r'$\\log M_{200c}$')
sl, ic, cx, cy = linfit(concentration[valid], icl_frac[valid], log_x=False)
xx = np.linspace(cx.min(), cx.max(), 100)
ax.plot(xx, sl*xx+ic, 'k-', lw=1.8, label=f'β = {sl:.3f}')
ax.set_xlabel('Concentración NFW $c_{200}$')
ax.set_ylabel(r'$f_{\\rm ICL}$')
ax.set_title('Fracción ICL vs Concentración (Fig. 4)')
ax.legend()
plt.tight_layout()
plt.savefig('fig04_icl_vs_concentracion.pdf', bbox_inches='tight')
plt.show()
"""),

md("### Figura 5 – Fracción ICL vs tiempo desde el último merger mayor (= Fig. 5)"),

code("""fig, ax = plt.subplots(figsize=(7, 5))
valid = np.isfinite(t_last_merger) & np.isfinite(icl_frac)
sc = ax.scatter(t_last_merger[valid], icl_frac[valid],
                c=lM[valid], cmap='viridis', s=25, alpha=0.8)
plt.colorbar(sc, ax=ax, label=r'$\\log M_{200c}$')
sl, ic, tx, ty = linfit(t_last_merger[valid], icl_frac[valid], log_x=False)
xx = np.linspace(tx.min(), tx.max(), 100)
ax.plot(xx, sl*xx+ic, 'k-', lw=1.8, label=f'β = {sl:.3f}')
ax.set_xlabel('Lookback time desde último merger mayor [Gyr]')
ax.set_ylabel(r'$f_{\\rm ICL}$')
ax.set_title('Fracción ICL vs Tiempo merger (Fig. 5)')
ax.legend()
plt.tight_layout()
plt.savefig('fig05_icl_vs_merger.pdf', bbox_inches='tight')
plt.show()
print(f"β = {sl:.3f}  (positivo → más tiempo desde merger = mayor f_ICL ✓)")
"""),

]

# ═════════════════════════════════════════════════════════════════════════════
# NOTEBOOK 2 – Propiedades estelares y componentes (§4, §5, §6)
# ═════════════════════════════════════════════════════════════════════════════
nb2 = nbf.v4.new_notebook()
nb2.cells = [

md("""# 02 · Propiedades Estelares y Componentes del BCG+ICL
**Metodología:** §4, §5, §6 de Mayes+2026

Clasificamos las partículas estelares usando los catálogos de ensamblaje estelar
(Rodriguez-Gomez+2016, disponibles en el cluster):
- **In situ** (canal 0): formadas en el subhalo central
- **Mergers completados** (canal 1): de galaxias totalmente disruptadas
- **Stripped / Sobrevivientes** (canal 2): de galaxias aún existentes

Luego estudiamos **metalicidad**, **color B−V** y **edad** del ICL y BCG
con sus gradientes radiales. Equivalente a **Figs. 7–12** de Mayes+2026.
"""),

code("""import sys, os, pickle
import numpy as np
import h5py
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import linregress
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, './original_shift_code')
import illustris_python as il
import Catalogue
import params_icl as P

%matplotlib inline
plt.rcParams.update({'figure.dpi': 110, 'font.size': 12,
                     'axes.spines.top': False, 'axes.spines.right': False})

cosmo = FlatLambdaCDM(H0=67.74, Om0=0.3089)
Header = il.groupcat.loadHeader(P.basePath, P.SNAP)

with h5py.File(P.CATALOG_OUT, 'r') as f:
    group_idx   = f['group_idx'][:]
    M200c       = f['M200c_Msun'][:]
    R200c       = f['R200c_kpc'][:]
    GroupPos    = f['GroupPos_kpc'][:]
    bcg_sub_idx = f['bcg_sub_idx'][:]
    icl_frac    = f['icl_frac'][:]
    dyn_state   = f['dyn_state'][:] if 'dyn_state' in f else None

n_cl = len(group_idx)
lM   = np.log10(M200c)
COLORS_COMP = {0: '#E91E63', 1: '#9C27B0', 2: '#00BCD4'}
LABELS_COMP = {0: 'In situ', 1: 'Mergers completados', 2: 'Stripped (sobrev.)'}
print(f"Catálogo cargado: {n_cl} cúmulos")
"""),

md("""## §2.3 – Catálogos de ensamblaje estelar

Ruta en el cluster: `P.SA_FILE`
(Archivo: `stellar_assembly_TNG100-1_099.hdf5`)

Estructura esperada:
```
ParticleID       : ID único de cada partícula estelar (para cruzar con el snapshot)
AssemblyChannel  : 0=in situ, 1=merger completado, 2=stripped
ProgGalaxyMass   : masa máxima histórica del progenitor [1e10 M☉/h]
ProgGalaxyID     : ID del subhalo progenitor
```
"""),

code("""with h5py.File(P.SA_FILE, 'r') as f:
    print("Campos disponibles:", list(f.keys()))
    sa_pid      = f['ParticleID'][:]
    sa_channel  = f['AssemblyChannel'][:]
    sa_progmass = f['ProgGalaxyMass'][:] * P.UM   # M☉
    sa_progid   = f['ProgGalaxyID'][:]

# Índice de búsqueda rápida id → posición en el catálogo SA
sa_map = dict(zip(sa_pid, np.arange(len(sa_pid))))

print(f"\\nTotal partículas SA : {len(sa_pid):,}")
print(f"  In situ (0) : {(sa_channel==0).sum():,}")
print(f"  Mergers (1) : {(sa_channel==1).sum():,}")
print(f"  Stripped(2) : {(sa_channel==2).sum():,}")
"""),

md("## Función principal: cargar y clasificar partículas de un cúmulo"),

code("""# Reutilizamos rotate_by_inertia_tensor del notebook 01
# (reproducida aquí para que este notebook sea autosuficiente)

from scipy.interpolate import interp1d

def rotate_by_inertia_tensor(pos_rel, mass, r_lim=np.inf):
    dist = np.linalg.norm(pos_rel, axis=1)
    ok   = (dist > 0) & (dist <= r_lim) & np.isfinite(mass)
    p, m = pos_rel[ok], mass[ok]
    if m.sum() == 0 or len(m) < 4:
        return pos_rel, np.eye(3)
    w    = 1.0 / dist[ok]**2
    mtot = np.sum(m)
    Ixx  = np.sum(m * p[:,0]**2 * w) / mtot
    Iyy  = np.sum(m * p[:,1]**2 * w) / mtot
    Izz  = np.sum(m * p[:,2]**2 * w) / mtot
    Ixy  = np.sum(m * p[:,0] * p[:,1] * w) / mtot
    Ixz  = np.sum(m * p[:,0] * p[:,2] * w) / mtot
    Iyz  = np.sum(m * p[:,1] * p[:,2] * w) / mtot
    I = np.array([[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]])
    ev, evec = np.linalg.eigh(I)
    R_mat = evec[:, np.argsort(ev)[::-1]].T
    return pos_rel @ R_mat.T, R_mat

def sb_profile_r(r_2d, lum_r, r_max_kpc, n_bins=60):
    r_bins = np.logspace(np.log10(0.5), np.log10(r_max_kpc), n_bins+1)
    r_mid  = np.sqrt(r_bins[:-1]*r_bins[1:])
    mu_r   = np.full(n_bins, np.nan)
    for k,(r1,r2) in enumerate(zip(r_bins[:-1],r_bins[1:])):
        mk = (r_2d>=r1)&(r_2d<r2)
        if mk.sum()==0: continue
        sl = lum_r[mk].sum() / (np.pi*((r2*1e3)**2-(r1*1e3)**2))
        if sl>0: mu_r[k] = P.SB_CONST - 2.5*np.log10(sl)
    return r_mid, mu_r

def holmberg_radius(r_mid, mu_r, mu_cut=P.MU_HOLMBERG):
    valid = np.isfinite(mu_r)&(r_mid>0)
    if valid.sum()<3: return np.nan
    rv,mv = r_mid[valid],mu_r[valid]
    idx   = np.argsort(rv)
    rv,mv = rv[idx],mv[idx]
    if mv[0]>mu_cut or mv[-1]<mu_cut: return np.nan
    try:
        r_h = float(interp1d(mv,rv,fill_value='extrapolate')(mu_cut))
        return r_h if 0<r_h<=rv[-1]*1.2 else np.nan
    except: return np.nan

def hmr(r, m):
    \"\"\"Radio de semi-masa (Half-Mass Radius).\"\"\"
    if len(r)==0 or m.sum()==0: return np.nan
    idx = np.argsort(r)
    cum = np.cumsum(m[idx])
    i50 = np.searchsorted(cum, cum[-1]/2)
    return r[idx[min(i50, len(r)-1)]]

def load_cluster_stars(sub_id, cen_pos, header=Header):
    \"\"\"
    Carga todas las partículas estelares del subhalo central,
    las rota, separa BCG/ICL y las clasifica en los 3 componentes.
    Devuelve un diccionario con todos los arrays necesarios.
    \"\"\"
    fields = ['Coordinates','Masses','ParticleIDs',
              'GFM_StellarPhotometrics','GFM_Metallicity','GFM_StellarFormationTime']
    st = il.snapshot.loadSubhalo(P.basePath, P.SNAP, int(sub_id), 'stars', fields=fields)

    pos   = Catalogue.Distance_3D(st['Coordinates']*P.UL, cen_pos, header['BoxSize']*P.UL)
    mass  = st['Masses'] * P.UM
    pids  = st['ParticleIDs']
    phot  = st['GFM_StellarPhotometrics']
    metal = st['GFM_Metallicity']
    aform = st['GFM_StellarFormationTime']

    # Clasificación SA
    channel   = np.full(len(pids), -1, dtype=int)
    progmass  = np.full(len(pids), np.nan)
    progid    = np.full(len(pids), -1, dtype=np.int64)
    for j, pid in enumerate(pids):
        if pid in sa_map:
            idx = sa_map[pid]
            channel[j]  = int(sa_channel[idx])
            progmass[j] = float(sa_progmass[idx])
            progid[j]   = int(sa_progid[idx])

    # Rotación
    pos_rot, _ = rotate_by_inertia_tensor(pos, mass)
    r_2d = np.sqrt(pos_rot[:,0]**2 + pos_rot[:,1]**2)
    lum_r = 10**((P.M_SUN_R_AB - phot[:,5]) / 2.5)

    # Radio de Holmberg
    r_m, mu = sb_profile_r(r_2d, lum_r, np.percentile(r_2d, 99))
    r_h = holmberg_radius(r_m, mu)
    mask_bcg = r_2d <= r_h if np.isfinite(r_h) else np.zeros(len(r_2d), bool)
    mask_icl = ~mask_bcg

    return {'mass':mass, 'phot':phot, 'metal':metal, 'aform':aform,
            'channel':channel, 'progmass':progmass, 'progid':progid,
            'r_2d':r_2d, 'pos_rot':pos_rot, 'lum_r':lum_r,
            'r_holmberg':r_h, 'mask_bcg':mask_bcg, 'mask_icl':mask_icl}
"""),

md("## §4 – Fracciones de masa por componente y HMR (Figs. 7, 8, 9)"),

code("""print("Calculando componentes y HMR para todos los cúmulos...")

# Resultados: fracciones y HMR para cada región × canal
rows = []
for i in range(n_cl):
    try:
        d = load_cluster_stars(bcg_sub_idx[i], GroupPos[i])
    except Exception as e:
        print(f"  Error cúmulo {i}: {e}")
        rows.append(None); continue

    row = {'i': i, 'M200c': M200c[i], 'r_h': d['r_holmberg']}
    for region, mask in [('bcgicl', np.ones(len(d['mass']),bool)),
                          ('icl',   d['mask_icl']),
                          ('bcg',   d['mask_bcg'])]:
        m_tot = d['mass'][mask].sum()
        for ch in [0,1,2]:
            mk = mask & (d['channel']==ch)
            row[f'{region}_f{ch}'] = d['mass'][mk].sum()/m_tot if m_tot>0 else np.nan
            row[f'{region}_hmr{ch}'] = hmr(d['r_2d'][mk], d['mass'][mk])
    rows.append(row)
    if (i+1)%10==0: print(f"  {i+1}/{n_cl}", end="\\r")

print("\\nListo.")
valid_rows = [r for r in rows if r is not None]
lM_v = np.array([r['M200c'] for r in valid_rows])
lM_v = np.log10(lM_v)
"""),

code("""def linfit(x, y, log_x=False):
    ok = np.isfinite(x)&np.isfinite(y)
    xx = np.log10(x[ok]) if log_x else x[ok]
    sl,ic,*_ = linregress(xx, y[ok])
    return sl, ic

# ── Figura 8: Fracciones BCG+ICL ─────────────────────────────────────────
fig, axes = plt.subplots(1,2,figsize=(13,5))
for ax, region, title in zip(axes,['bcgicl','icl'],['BCG+ICL','ICL']):
    for ch in [0,1,2]:
        fv = np.array([r[f'{region}_f{ch}'] for r in valid_rows])
        ok = np.isfinite(fv)
        ax.scatter(lM_v[ok], fv[ok], color=COLORS_COMP[ch],
                   label=LABELS_COMP[ch], s=15, alpha=0.7)
        sl,ic = linfit(lM_v[ok], fv[ok])
        xx = np.linspace(lM_v.min(), lM_v.max(), 100)
        ax.plot(xx, sl*xx+ic, color=COLORS_COMP[ch], lw=1.8, ls='--',
                label=f'β={sl:.3f}')
    ax.set_xlabel(r'$\\log M_{{200c}}$'); ax.set_ylabel('Fracción de masa')
    ax.set_title(f'Componentes del {title} (Fig. {"8" if title=="BCG+ICL" else "9"})')
    ax.legend(fontsize=8, ncol=2)
plt.tight_layout()
plt.savefig('fig08_09_componentes.pdf', bbox_inches='tight')
plt.show()
"""),

code("""# ── Figura 7: HMR de los componentes ─────────────────────────────────────
fig, axes = plt.subplots(1,2,figsize=(13,5))
for ax, region, title in zip(axes,['bcgicl','icl'],['BCG+ICL','ICL']):
    for ch in [0,1,2]:
        hv = np.array([r[f'{region}_hmr{ch}'] for r in valid_rows])
        ok = np.isfinite(hv) & (hv>0)
        ax.scatter(lM_v[ok], hv[ok], color=COLORS_COMP[ch],
                   label=LABELS_COMP[ch], s=15, alpha=0.7)
        if ok.sum()>3:
            sl,ic = linfit(lM_v[ok], np.log10(hv[ok]))
            xx = np.linspace(lM_v.min(), lM_v.max(), 100)
            ax.plot(xx, 10**(sl*xx+ic), color=COLORS_COMP[ch], lw=1.5, ls='--')
    ax.set_xlabel(r'$\\log M_{{200c}}$'); ax.set_ylabel('HMR [kpc]')
    ax.set_yscale('log'); ax.set_title(f'Radio de semi-masa – {title} (Fig. 7)')
    ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig('fig07_hmr_componentes.pdf', bbox_inches='tight')
plt.show()
"""),

md("## §6.1 – Metalicidad (equivalente a Fig. 10 del paper)"),

code("""print("Calculando propiedades estelares para todos los cúmulos...")

def lookback(aform):
    z = np.clip(1/np.where(aform>0,aform,1e-4) - 1, 0, 30)
    return cosmo.lookback_time(z).value

mean_Z_icl   = np.full(n_cl, np.nan)
mean_Z_bcg   = np.full(n_cl, np.nan)
mean_BV_icl  = np.full(n_cl, np.nan)
mean_BV_bcg  = np.full(n_cl, np.nan)
mean_age_icl = np.full(n_cl, np.nan)
mean_age_bcg = np.full(n_cl, np.nan)
grad_Z       = np.full(n_cl, np.nan)   # gradiente metalicidad BCG+ICL
grad_BV      = np.full(n_cl, np.nan)   # gradiente color BCG+ICL

for i in range(n_cl):
    try: d = load_cluster_stars(bcg_sub_idx[i], GroupPos[i])
    except: continue

    for mask, Z_arr, BV_arr, age_arr in [
        (d['mask_icl'], mean_Z_icl,  mean_BV_icl,  mean_age_icl),
        (d['mask_bcg'], mean_Z_bcg,  mean_BV_bcg,  mean_age_bcg),
    ]:
        if mask.sum()==0: continue
        Z_arr[i]   = np.log10(np.average(d['metal'][mask], weights=d['mass'][mask])
                               / P.Z_SUN)
        lr_k       = 10**((P.M_SUN_R_AB - d['phot'][mask,5]) / 2.5)
        BV_arr[i]  = np.average(d['phot'][mask,1] - d['phot'][mask,2], weights=lr_k)
        age_arr[i] = np.average(lookback(d['aform'][mask]), weights=d['mass'][mask])

    # Gradiente radial de metalicidad (BCG+ICL completo)
    r2 = d['r_2d']
    rg = np.logspace(np.log10(max(r2.min(),0.5)), np.log10(np.percentile(r2,99)), 21)
    rc = np.sqrt(rg[:-1]*rg[1:])
    Z_r = np.array([np.log10(np.average(d['metal'][(r2>=r1)&(r2<r2_)],
                                          weights=d['mass'][(r2>=r1)&(r2<r2_)])/P.Z_SUN)
                     if ((r2>=r1)&(r2<r2_)).sum()>0 else np.nan
                     for r1,r2_ in zip(rg[:-1],rg[1:])])
    BV_r= np.array([np.average(d['phot'][(r2>=r1)&(r2<r2_),1]-d['phot'][(r2>=r1)&(r2<r2_),2],
                                 weights=10**((P.M_SUN_R_AB-d['phot'][(r2>=r1)&(r2<r2_),5])/2.5))
                     if ((r2>=r1)&(r2<r2_)).sum()>0 else np.nan
                     for r1,r2_ in zip(rg[:-1],rg[1:])])
    ok_Z  = np.isfinite(Z_r) & np.isfinite(np.log10(rc))
    ok_BV = np.isfinite(BV_r)& np.isfinite(np.log10(rc))
    if ok_Z.sum()>=3:
        sl,*_ = linregress(np.log10(rc[ok_Z]), Z_r[ok_Z]); grad_Z[i]=sl
    if ok_BV.sum()>=3:
        sl,*_ = linregress(np.log10(rc[ok_BV]), BV_r[ok_BV]); grad_BV[i]=sl
    if (i+1)%10==0: print(f"  {i+1}/{n_cl}", end="\\r")
print("\\nListo.")
"""),

code("""fig, axes = plt.subplots(1,2,figsize=(13,5))

# Panel izq: metalicidad media
ax = axes[0]
for arr, col, lbl in [(mean_Z_icl,'royalblue','ICL'),(mean_Z_bcg,'tomato','BCG')]:
    ok = np.isfinite(arr)
    ax.scatter(lM[ok], arr[ok], color=col, s=20, alpha=0.7, label=lbl)
    sl,ic = linfit(lM[ok], arr[ok])
    xx = np.linspace(lM.min(), lM.max(), 100)
    ax.plot(xx, sl*xx+ic, color=col, lw=1.8, ls='--')
ax.set_xlabel(r'$\\log M_{{200c}}$')
ax.set_ylabel(r'$\\langle\\log Z_*/Z_\\odot\\rangle$')
ax.set_title('Metalicidad media ICL vs BCG (Fig. 10 top)'); ax.legend()

# Panel der: gradiente de metalicidad
ax = axes[1]
ok = np.isfinite(grad_Z)
ax.scatter(lM[ok], grad_Z[ok], color='slategray', s=20, alpha=0.7)
ax.axhline(0, color='k', lw=1, ls='--')
sl,ic = linfit(lM[ok], grad_Z[ok])
xx = np.linspace(lM.min(), lM.max(), 100)
ax.plot(xx, sl*xx+ic, 'k-', lw=1.8, label=f'β={sl:.3f}')
ax.set_xlabel(r'$\\log M_{{200c}}$')
ax.set_ylabel(r'Gradiente $d\\log Z / d\\log r$')
ax.set_title('Gradiente metalicidad (Fig. 10 bottom)'); ax.legend()

plt.tight_layout()
plt.savefig('fig10_metalicidad.pdf', bbox_inches='tight')
plt.show()
print(f"Todos los gradientes negativos: {(grad_Z[ok]<0).all()}")
"""),

md("### §6.2 – Color B−V (Fig. 11) · §6.3 – Edades (Fig. 12)"),

code("""fig, axes = plt.subplots(1,3,figsize=(16,5))

# ── B−V medio ──────────────────────────────────────────────────────────────
ax = axes[0]
for arr,col,lbl in [(mean_BV_icl,'royalblue','ICL'),(mean_BV_bcg,'tomato','BCG')]:
    ok = np.isfinite(arr)
    ax.scatter(lM[ok],arr[ok],color=col,s=20,alpha=0.7,label=lbl)
    sl,ic = linfit(lM[ok],arr[ok])
    ax.plot(np.linspace(lM.min(),lM.max(),100),
            sl*np.linspace(lM.min(),lM.max(),100)+ic,color=col,lw=1.8,ls='--')
ax.set_xlabel(r'$\\log M_{{200c}}$'); ax.set_ylabel(r'$\\langle B-V\\rangle$')
ax.set_title('Color B−V medio (Fig. 11 top)'); ax.legend()

# ── Gradiente de color ─────────────────────────────────────────────────────
ax = axes[1]
ok = np.isfinite(grad_BV)
ax.scatter(lM[ok],grad_BV[ok],color='slategray',s=20,alpha=0.7)
ax.axhline(0,color='k',lw=1,ls='--')
sl,ic = linfit(lM[ok],grad_BV[ok])
ax.plot(np.linspace(lM.min(),lM.max(),100),
        sl*np.linspace(lM.min(),lM.max(),100)+ic,'k-',lw=1.8,label=f'β={sl:.3f}')
ax.set_xlabel(r'$\\log M_{{200c}}$'); ax.set_ylabel(r'Gradiente $d(B-V)/d\\log r$')
ax.set_title('Gradiente color (Fig. 11 bottom)'); ax.legend()

# ── Edad estelar ──────────────────────────────────────────────────────────
ax = axes[2]
for arr,col,lbl in [(mean_age_icl,'royalblue','ICL'),(mean_age_bcg,'tomato','BCG')]:
    ok = np.isfinite(arr)
    ax.scatter(lM[ok],arr[ok],color=col,s=20,alpha=0.7,label=lbl)
    sl,ic = linfit(lM[ok],arr[ok])
    ax.plot(np.linspace(lM.min(),lM.max(),100),
            sl*np.linspace(lM.min(),lM.max(),100)+ic,color=col,lw=1.8,ls='--')
ax.set_xlabel(r'$\\log M_{{200c}}$'); ax.set_ylabel('Edad [Gyr lookback]')
ax.set_title('Edad estelar ICL vs BCG (Fig. 12)'); ax.legend()

plt.tight_layout()
plt.savefig('fig11_12_color_edad.pdf', bbox_inches='tight')
plt.show()
print(f"ICL más azul que BCG: {np.nanmean(mean_BV_icl)<np.nanmean(mean_BV_bcg)}")
print(f"Δ(B-V) = {np.nanmean(mean_BV_bcg)-np.nanmean(mean_BV_icl):.3f}")
"""),

]

# ═════════════════════════════════════════════════════════════════════════════
# NOTEBOOK 3 – Análisis de progenitores (§7)
# ═════════════════════════════════════════════════════════════════════════════
nb3 = nbf.v4.new_notebook()
nb3.cells = [

md("""# 03 · Análisis de Progenitores del ICL y BCG
**Metodología:** §7 de Mayes+2026

1. **§7.1** – Masa media de progenitores ICL vs BCG (Fig. 13)
2. **§7.2** – Fracción de progenitores compartidos ICL↔BCG (Fig. 14)
3. **§7.3** – Localización de material de progenitores alta/baja masa (Fig. 15)
4. **§7.4** – Progenitor más significativo por componente
"""),

code("""import sys, os, pickle
import numpy as np
import h5py
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import linregress
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, './original_shift_code')
import illustris_python as il
import Catalogue
import params_icl as P

%matplotlib inline
plt.rcParams.update({'figure.dpi': 110, 'font.size': 12,
                     'axes.spines.top': False, 'axes.spines.right': False})

cosmo = FlatLambdaCDM(H0=67.74, Om0=0.3089)
Header = il.groupcat.loadHeader(P.basePath, P.SNAP)

with h5py.File(P.CATALOG_OUT, 'r') as f:
    group_idx   = f['group_idx'][:]
    M200c       = f['M200c_Msun'][:]
    GroupPos    = f['GroupPos_kpc'][:]
    bcg_sub_idx = f['bcg_sub_idx'][:]
    dyn_state   = f['dyn_state'][:] if 'dyn_state' in f else None

n_cl = len(group_idx)
lM   = np.log10(M200c)
COLORS_STATE = {0:'#2196F3', 1:'#FF9800', 2:'#F44336'}
LABELS_STATE = {0:'Relajado', 1:'Intermedio', 2:'Perturbado'}
print(f"Catálogo cargado: {n_cl} cúmulos")
"""),

code("""# Cargar catálogo de ensamblaje estelar
with h5py.File(P.SA_FILE, 'r') as f:
    sa_pid      = f['ParticleID'][:]
    sa_channel  = f['AssemblyChannel'][:]
    sa_progmass = f['ProgGalaxyMass'][:] * P.UM
    sa_progid   = f['ProgGalaxyID'][:]

sa_map = dict(zip(sa_pid, np.arange(len(sa_pid))))
print(f"Catálogo SA: {len(sa_pid):,} partículas")
"""),

code("""# Funciones auxiliares (mismas que notebook 02)
def rotate_by_inertia_tensor(pos_rel, mass, r_lim=np.inf):
    dist = np.linalg.norm(pos_rel, axis=1)
    ok   = (dist>0)&(dist<=r_lim)&np.isfinite(mass)
    p,m  = pos_rel[ok], mass[ok]
    if m.sum()==0 or len(m)<4: return pos_rel, np.eye(3)
    w = 1/dist[ok]**2; mt = m.sum()
    Ixx=np.sum(m*p[:,0]**2*w)/mt; Iyy=np.sum(m*p[:,1]**2*w)/mt
    Izz=np.sum(m*p[:,2]**2*w)/mt; Ixy=np.sum(m*p[:,0]*p[:,1]*w)/mt
    Ixz=np.sum(m*p[:,0]*p[:,2]*w)/mt; Iyz=np.sum(m*p[:,1]*p[:,2]*w)/mt
    I=np.array([[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]])
    ev,evec=np.linalg.eigh(I)
    R=evec[:,np.argsort(ev)[::-1]].T
    return pos_rel@R.T, R

def sb_and_holmberg(r_2d, lum_r, n_bins=60):
    r_bins=np.logspace(np.log10(0.5),np.log10(np.percentile(r_2d,99)),n_bins+1)
    r_mid=np.sqrt(r_bins[:-1]*r_bins[1:]); mu_r=np.full(n_bins,np.nan)
    for k,(r1,r2) in enumerate(zip(r_bins[:-1],r_bins[1:])):
        mk=(r_2d>=r1)&(r_2d<r2)
        if mk.sum()==0: continue
        sl=lum_r[mk].sum()/(np.pi*((r2*1e3)**2-(r1*1e3)**2))
        if sl>0: mu_r[k]=P.SB_CONST-2.5*np.log10(sl)
    valid=np.isfinite(mu_r)&(r_mid>0)
    r_h=np.nan
    if valid.sum()>=3:
        rv,mv=r_mid[valid],mu_r[valid]
        idx=np.argsort(rv); rv,mv=rv[idx],mv[idx]
        if mv[0]<=P.MU_HOLMBERG<=mv[-1]:
            try: r_h=float(interp1d(mv,rv,fill_value='extrapolate')(P.MU_HOLMBERG))
            except: pass
    return r_h

def hmr(r, m):
    if len(r)==0 or m.sum()==0: return np.nan
    idx=np.argsort(r); cum=np.cumsum(m[idx])
    i50=np.searchsorted(cum,cum[-1]/2)
    return r[idx[min(i50,len(r)-1)]]

def linfit(x, y):
    ok=np.isfinite(x)&np.isfinite(y)
    if ok.sum()<3: return np.nan,np.nan
    sl,ic,*_=linregress(x[ok],y[ok]); return sl,ic

def load_prog_data(sub_id, cen_pos, header=Header):
    \"\"\"Carga partículas + clasificación SA + separa BCG/ICL.\"\"\"
    fields=['Coordinates','Masses','ParticleIDs','GFM_StellarPhotometrics']
    st=il.snapshot.loadSubhalo(P.basePath,P.SNAP,int(sub_id),'stars',fields=fields)
    pos  =Catalogue.Distance_3D(st['Coordinates']*P.UL, cen_pos, header['BoxSize']*P.UL)
    mass =st['Masses']*P.UM
    pids =st['ParticleIDs']
    phot =st['GFM_StellarPhotometrics']
    channel=np.full(len(pids),-1,dtype=int)
    progm  =np.full(len(pids),np.nan)
    progid =np.full(len(pids),-1,dtype=np.int64)
    for j,pid in enumerate(pids):
        if pid in sa_map:
            k=sa_map[pid]; channel[j]=int(sa_channel[k])
            progm[j]=float(sa_progmass[k]); progid[j]=int(sa_progid[k])
    pos_rot,_=rotate_by_inertia_tensor(pos,mass)
    r_2d=np.sqrt(pos_rot[:,0]**2+pos_rot[:,1]**2)
    lum_r=10**((P.M_SUN_R_AB-phot[:,5])/2.5)
    r_h=sb_and_holmberg(r_2d,lum_r)
    mask_bcg = r_2d<=r_h if np.isfinite(r_h) else np.zeros(len(r_2d),bool)
    return {'mass':mass,'r_2d':r_2d,'channel':channel,
            'progmass':progm,'progid':progid,'lum_r':lum_r,
            'mask_bcg':mask_bcg,'mask_icl':~mask_bcg,'r_h':r_h}
"""),

md("## §7.1 – Masa media de progenitores ICL vs BCG (Fig. 13)"),

code("""print("§7.1 – Masa media de progenitores...")
mean_pm_icl = np.full(n_cl, np.nan)
mean_pm_bcg = np.full(n_cl, np.nan)
frac_hi_icl = np.full(n_cl, np.nan)   # fracción ICL de progenitores M>thresh

for i in range(n_cl):
    try: d = load_prog_data(bcg_sub_idx[i], GroupPos[i])
    except: continue
    for mask, arr in [(d['mask_icl'], mean_pm_icl), (d['mask_bcg'], mean_pm_bcg)]:
        pm = d['progmass'][mask]; m = d['mass'][mask]
        ex = d['channel'][mask] != 0
        ok = ex & np.isfinite(pm)
        if ok.sum()>0: arr[i] = np.average(pm[ok], weights=m[ok])
    # fracción de alta masa en ICL
    ex_icl = d['mask_icl'] & (d['channel']!=0) & np.isfinite(d['progmass'])
    if ex_icl.sum()>0:
        hi = d['mass'][ex_icl & (d['progmass']>=P.M_PROG_THRESH)].sum()
        frac_hi_icl[i] = hi / d['mass'][ex_icl].sum()
    if (i+1)%10==0: print(f"  {i+1}/{n_cl}",end="\\r")
print("\\nListo.")
print(f"Fracción media ICL de prog. M>1e10 M☉: {np.nanmean(frac_hi_icl):.2f} ± {np.nanstd(frac_hi_icl):.2f}")
"""),

code("""fig, ax = plt.subplots(figsize=(7,6))
valid = np.isfinite(mean_pm_icl)&np.isfinite(mean_pm_bcg)
sc = ax.scatter(np.log10(mean_pm_bcg[valid]), np.log10(mean_pm_icl[valid]),
                c=lM[valid], cmap='viridis', s=25, alpha=0.8)
plt.colorbar(sc, ax=ax, label=r'$\\log M_{{200c}}$')
x = np.log10(mean_pm_bcg[valid]); y = np.log10(mean_pm_icl[valid])
sl,ic = linfit(x,y)
xx = np.linspace(x.min(),x.max(),100)
ax.plot(xx,sl*xx+ic,'orange',lw=2,label=f'β={sl:.3f}')
lim=[min(x.min(),y.min()),max(x.max(),y.max())]
ax.plot(lim,lim,'k--',lw=1.2,label='1:1')
ax.set_xlabel(r'$\\log\\langle M_{{prog,BCG}}\\rangle$')
ax.set_ylabel(r'$\\log\\langle M_{{prog,ICL}}\\rangle$')
ax.set_title('Masa media de progenitores (Fig. 13)'); ax.legend()
pct_below = (mean_pm_icl[valid]<mean_pm_bcg[valid]).mean()*100
ax.text(0.03,0.97,f'ICL<BCG en {pct_below:.0f}% de los cúmulos',
        transform=ax.transAxes,va='top',fontsize=10)
plt.tight_layout()
plt.savefig('fig13_masa_progenitores.pdf', bbox_inches='tight')
plt.show()
"""),

md("## §7.2 – Progenitores compartidos ICL ↔ BCG (Fig. 14)"),

code("""print("§7.2 – Progenitores compartidos...")
shared_icl = np.full(n_cl, np.nan)   # fracción ICL cuyos prog también están en BCG
shared_bcg = np.full(n_cl, np.nan)   # fracción BCG cuyos prog también están en ICL

for i in range(n_cl):
    try: d = load_prog_data(bcg_sub_idx[i], GroupPos[i])
    except: continue
    ex = d['channel'] != 0
    prog_icl = set(d['progid'][d['mask_icl']&ex&(d['progid']>=0)].tolist())
    prog_bcg = set(d['progid'][d['mask_bcg']&ex&(d['progid']>=0)].tolist())
    if len(prog_icl)>0: shared_icl[i] = len(prog_icl&prog_bcg)/len(prog_icl)
    if len(prog_bcg)>0: shared_bcg[i] = len(prog_icl&prog_bcg)/len(prog_bcg)
    if (i+1)%10==0: print(f"  {i+1}/{n_cl}",end="\\r")
print("\\nListo.")
print(f"Progenitores compartidos ICL→BCG: {np.nanmean(shared_icl)*100:.1f}%  (paper: 93%)")
print(f"Progenitores compartidos BCG→ICL: {np.nanmean(shared_bcg)*100:.1f}%  (paper: 91%)")
"""),

code("""fig, axes = plt.subplots(1,2,figsize=(13,5))
bins = np.linspace(0,1,25)
ax = axes[0]
ax.hist(shared_icl[np.isfinite(shared_icl)],bins=bins,color='royalblue',alpha=0.7,label='ICL→BCG')
ax.hist(shared_bcg[np.isfinite(shared_bcg)],bins=bins,color='tomato',alpha=0.7,label='BCG→ICL')
for x,col in [(shared_icl,'royalblue'),(shared_bcg,'tomato')]:
    ax.axvline(np.nanmedian(x),color=col,lw=2,ls='--')
ax.set_xlabel('Fracción de progenitores compartidos'); ax.set_ylabel('N cúmulos')
ax.set_title('Progenitores compartidos (Fig. 14 top)'); ax.legend()

ax = axes[1]
if dyn_state is not None:
    for s in [0,1,2]:
        m=dyn_state==s
        ax.scatter(shared_icl[m],shared_bcg[m],color=COLORS_STATE[s],
                   label=LABELS_STATE[s],s=20,alpha=0.8)
else:
    ax.scatter(shared_icl,shared_bcg,color='steelblue',s=20,alpha=0.7)
ax.plot([0,1],[0,1],'k--',lw=1)
ax.set_xlabel('Fracción ICL→BCG'); ax.set_ylabel('Fracción BCG→ICL')
ax.set_title('Por estado dinámico (Fig. 14 bottom)'); ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig('fig14_progenitores_compartidos.pdf', bbox_inches='tight')
plt.show()
"""),

md("## §7.3 – HMR del material por masa de progenitor (Fig. 15)"),

code("""print("§7.3 – HMR por masa de progenitor...")
# Para cada cúmulo: HMR del material de progenitores alta (>thresh) y baja masa
# separado por canal (todo, mergers, stripped)
results_hmr = {k: {'hi': np.full(n_cl,np.nan), 'lo': np.full(n_cl,np.nan)}
               for k in ['all','mer','str']}

for i in range(n_cl):
    try: d = load_prog_data(bcg_sub_idx[i], GroupPos[i])
    except: continue
    pm = d['progmass']; r = d['r_2d']; m = d['mass']; ch = d['channel']
    ex = ch != 0
    for key, filt in [('all', ex), ('mer', ex&(ch==1)), ('str', ex&(ch==2))]:
        for hi_lo, cond in [('hi', pm>=P.M_PROG_THRESH), ('lo', pm<P.M_PROG_THRESH)]:
            mk = filt & cond & np.isfinite(pm)
            if mk.sum()>0:
                results_hmr[key][hi_lo][i] = hmr(r[mk], m[mk])
    if (i+1)%10==0: print(f"  {i+1}/{n_cl}",end="\\r")
print("\\nListo.")
"""),

code("""fig, axes = plt.subplots(1,3,figsize=(15,5),sharey=True)
titles = ['Todo ex-situ','Mergers completados','Stripped (sobrev.)']
keys   = ['all','mer','str']

for ax,title,key in zip(axes,titles,keys):
    hi = results_hmr[key]['hi']
    lo = results_hmr[key]['lo']
    for arr,col,lbl in [(hi,'tomato',f'M>10¹⁰ M☉'),(lo,'royalblue',f'M<10¹⁰ M☉')]:
        ok = np.isfinite(arr)&(arr>0)
        ax.scatter(lM[ok],arr[ok],color=col,s=20,alpha=0.7,label=lbl)
        if ok.sum()>3:
            sl,ic=linfit(lM[ok],np.log10(arr[ok]))
            xx=np.linspace(lM.min(),lM.max(),100)
            ax.plot(xx,10**(sl*xx+ic),color=col,lw=1.8,ls='--',label=f'β={sl:.3f}')
    ax.set_xlabel(r'$\\log M_{{200c}}$')
    ax.set_ylabel('HMR [kpc]' if ax==axes[0] else '')
    ax.set_yscale('log'); ax.set_title(f'Fig. 15 – {title}'); ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig('fig15_hmr_por_masa_progenitor.pdf', bbox_inches='tight')
plt.show()
"""),

md("## §7.4 – Progenitor más significativo"),

code("""print("§7.4 – Progenitor más significativo...")
frac_top1_icl = np.full(n_cl, np.nan)
frac_top1_bcg = np.full(n_cl, np.nan)
n_for_90_icl  = np.full(n_cl, np.nan)

for i in range(n_cl):
    try: d = load_prog_data(bcg_sub_idx[i], GroupPos[i])
    except: continue
    ex = d['channel'] != 0
    for mask, arr1, arr90 in [(d['mask_icl'], frac_top1_icl, n_for_90_icl),
                               (d['mask_bcg'], frac_top1_bcg, None)]:
        mk = mask & ex & (d['progid']>=0)
        pids_k = d['progid'][mk]; ms_k = d['mass'][mk]
        if len(pids_k)==0: continue
        m_tot = ms_k.sum()
        contrib = {p: ms_k[pids_k==p].sum() for p in np.unique(pids_k)}
        sorted_c = sorted(contrib.values(), reverse=True)
        arr1[i] = sorted_c[0] / m_tot
        if arr90 is not None:
            cumf = np.cumsum(sorted_c)/m_tot
            arr90[i] = np.searchsorted(cumf, 0.9) + 1
    if (i+1)%10==0: print(f"  {i+1}/{n_cl}",end="\\r")
print("\\nListo.")
print(f"Contribución top progenitor ICL: {np.nanmean(frac_top1_icl):.3f}±{np.nanstd(frac_top1_icl):.3f}  (paper: 0.27±0.12)")
print(f"Contribución top progenitor BCG: {np.nanmean(frac_top1_bcg):.3f}±{np.nanstd(frac_top1_bcg):.3f}")
"""),

code("""fig, axes = plt.subplots(1,2,figsize=(13,5))

ax = axes[0]
bins = np.linspace(0,1,25)
ax.hist(frac_top1_icl[np.isfinite(frac_top1_icl)],bins=bins,color='royalblue',alpha=0.7,label='ICL')
ax.hist(frac_top1_bcg[np.isfinite(frac_top1_bcg)],bins=bins,color='tomato',alpha=0.7,label='BCG')
for x,col in [(frac_top1_icl,'royalblue'),(frac_top1_bcg,'tomato')]:
    ax.axvline(np.nanmedian(x),color=col,lw=2,ls='--')
ax.set_xlabel('Fracción del progenitor más significativo')
ax.set_ylabel('N cúmulos'); ax.set_title('Progenitor más significativo'); ax.legend()

ax = axes[1]
ok = np.isfinite(n_for_90_icl)
sc = ax.scatter(lM[ok], n_for_90_icl[ok], c=lM[ok], cmap='viridis', s=25, alpha=0.8)
plt.colorbar(sc, ax=ax, label=r'$\\log M_{{200c}}$')
ax.set_xlabel(r'$\\log M_{{200c}}$')
ax.set_ylabel('N progenitores para el 90% del ICL')
ax.set_title('Progenitores necesarios para 90% ICL')

plt.tight_layout()
plt.savefig('fig_progenitor_significativo.pdf', bbox_inches='tight')
plt.show()
"""),

md("""## Comparación con Mayes+2026

| Resultado | Mayes+2026 | Este análisis |
|-----------|-----------|---------------|
| % ICL de progenitores M>10¹⁰ M☉ | 65 ± 15% | *(ver arriba)* |
| Progenitores compartidos ICL→BCG | **93%** | *(ver arriba)* |
| Progenitores compartidos BCG→ICL | **91%** | *(ver arriba)* |
| Contribución top progenitor ICL | **0.27 ± 0.12** | *(ver arriba)* |
"""),

]

# ─────────────────────────────────────────────────────────────────────────────
# Guardar
# ─────────────────────────────────────────────────────────────────────────────
print("Generando notebooks...")
save(nb0, "00_catalogo.ipynb")
save(nb1, "01_separacion_BCG_ICL.ipynb")
save(nb2, "02_propiedades_estelares.ipynb")
save(nb3, "03_analisis_progenitores.ipynb")
print("¡Listo!")
