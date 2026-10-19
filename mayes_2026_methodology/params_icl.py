"""
Parámetros del análisis ICL/BCG: réplica de Mayes, Gómez & Monachesi 2026
(A&A 708, A124; arXiv:2506.16645) sobre TNG100-1, snapshot 99 (z = 0).

Todo número que use más de un notebook vive aquí. Lo que el paper NO
especifica está marcado con `# SUPUESTO:` (lista para confirmar con Mary).

Rutas: se pueden sobreescribir con variables de entorno; los valores por
defecto son los del TNG Lab (/home/tnguser/...).
"""
import os
from pathlib import Path

# ── Carpeta del repo (donde vive este archivo) ───────────────────────────
REPO_DIR = Path(__file__).resolve().parent

# ── Rutas a los datos de TNG ──────────────────────────────────────────────
basePath      = os.environ.get('TNG_BASEPATH', '/home/tnguser/sims.TNG/TNG100-1/output/')
POSTPROC_PATH = os.environ.get('TNG_POSTPROC', os.path.join(os.path.dirname(basePath.rstrip('/')), 'postprocessing'))

# Catálogo de ensamblaje estelar AGREGADO por subhalo (Rodriguez-Gomez+2016),
# el que se descarga desde la API (/api/TNG100-1/files/stellar_assembly.hdf5).
# Lo usa read_hdf5_SA.ipynb; la réplica usa el catálogo POR PARTÍCULA (abajo).
SA_FILE = os.environ.get('TNG_SA_FILE', str(REPO_DIR.parent / 'stellar_assembly.hdf5'))

# Catálogo de ensamblaje estelar POR PARTÍCULA (§2.3 del paper).
# SUPUESTO (no verificado): archivo stars_XXX.hdf5 en postprocessing/StellarAssembly,
# con un valor por partícula PartType4 en el MISMO orden que el snapshot.
# Si el nombre o los campos difieren en el Lab, cambiar aquí o con variables de entorno.
SA_STARS_FILE      = os.environ.get('TNG_SA_STARS_FILE',
                                    os.path.join(POSTPROC_PATH, 'StellarAssembly', 'stars_{snap:03d}.hdf5'))
SA_FIELD_INSITU    = os.environ.get('TNG_SA_FIELD_INSITU', 'InSitu')            # 1 in situ, 0 ex situ
SA_FIELD_ORIGIN    = os.environ.get('TNG_SA_FIELD_ORIGIN', 'AccretionOrigin')   # 0 merger completo, 1 ongoing, 2 flyby, -1 n/a

# ── Catálogos locales (código de shift, TNG2PROPnum) ──────────────────────
# Se conservan para inspeccion_catalogos*.ipynb y para el chequeo local de la muestra en 00.
CATALOG_PKL              = str(REPO_DIR / 'mary_haloshift_z0.pkl')
CATALOG_SUBHALOSHIFT_PKL = str(REPO_DIR / 'mary_subhaloshift_z0.pkl')

# ── Catálogo de salida de la réplica ──────────────────────────────────────
# Nombre nuevo para no pisar el catalogo_icl.hdf5 anterior (muestra fósil/no fósil).
CATALOG_OUT = os.environ.get('ICL_CATALOG_OUT', str(REPO_DIR / 'catalogo_icl_mayes2026.hdf5'))

# Carpeta local para cutouts de la API (01f)
CUTOUT_DIR = str(REPO_DIR / 'tng_data')

# ── Snapshot ──────────────────────────────────────────────────────────────
SNAP = 99            # z = 0

# ── Cosmología y unidades de TNG (Planck 2015 = Planck Collaboration XIII 2016) ──
# El paper cita Ωm=0.31, ΩΛ=0.69, Ωb=0.0486, h=0.677 (valores redondeados de TNG).
h     = 0.6774
H0    = 67.74
Om0   = 0.3089
Ob0   = 0.0486
L_BOX = 75000.0      # ckpc/h

# Posición: ckpc/h → kpc físicos = valor * a / h (a = 1 en snap 99)
# Masa:     1e10 M☉/h → M☉ = valor * 1e10 / h
UL   = 1.0 / h       # factor longitud (a = 1)
UM   = 1e10 / h      # factor masa → M☉
UM_H = 1e10          # factor masa → M☉/h

# ── Selección de la muestra (§2.1) ────────────────────────────────────────
# "127 clusters and groups above 1e13 M☉, 10 > 1e14, 1e13–2.6e14 M☉".
# Esos números se reproducen SOLO si el umbral está en M☉/h:
#   Group_M_Crit200 (unidades internas 1e10 M☉/h) ≥ 1e3  ⇔  M200c ≥ 1e13 M☉/h (= 1.476e13 M☉).
# SUPUESTO: el paper rotula M☉ pero usa M☉/h (inferencia fuerte: N=127, 10 > 1e14, máx 2.58e14).
M200_MIN_INTERNAL = 1e3      # Group_M_Crit200 ≥ 1e3 [1e10 M☉/h]
N_HALOS_PAPER     = 127

# SUPUESTO: las masas de halo y de progenitores se grafican/comparan en M☉/h
# (misma convención que reproduce la muestra). Cambiar a UM para M☉ físicas.
UM_PLOT = UM_H

# ── Estado dinámico (§2.1, §3; Aldás et al. 2023) ─────────────────────────
# offset 3D entre el BCG (partícula más ligada = SubhaloPos del central) y el CM del gas.
DYN_RELAXED   = 0.1      # D/R200 < 0.1 → relajado
DYN_DISTURBED = 0.4      # D/R200 > 0.4 → perturbado
# SUPUESTO: CM del gas calculado con las celdas de gas del FoF dentro de R200 (3D).
GAS_CM_R_FRAC = 1.0
DYN_LABELS = {0: 'Relajado', 1: 'Intermedio', 2: 'Perturbado'}

# ── Concentración NFW (§3; c = R200c/Rs, ajuste al perfil de densidad DM) ─
# SUPUESTO: partículas DM del FoF, 20 bins logarítmicos entre 0.05 y 1 R200,
# ajuste por mínimos cuadrados en log ρ (el paper sólo cita Anbajagane+2021).
NFW_R_MIN_FRAC = 0.05
NFW_R_MAX_FRAC = 1.0
NFW_N_BINS     = 20
# Masa de la partícula DM (constante en TNG; PartType1 no tiene 'Masses'):
# 1) MassTable[1] del Header del SNAPSHOT (el Header del catálogo de grupos no la trae);
# 2) si falta o es 0, el valor DOCUMENTADO de TNG100-1 (tng-project.org/data/downloads/TNG100-1/:
#    "dark matter particle mass [10^10 msun/h] 0.000505574296436975" ≈ 7.46e6 Msun).
# La concentración no depende de esta masa (sólo escala ρ).
M_DM_INTERNAL_DOC = 0.000505574296436975     # [1e10 Msun/h]

# ── Último major merger (§3) ──────────────────────────────────────────────
MAJOR_MERGER_RATIO = 1/5     # razón de masa estelar > 1/5
# SUPUESTO: razón con la masa estelar máxima pasada de cada rama (como
# illustris_python.sublink.numMergers / Rodriguez-Gomez+2015), usando
# SubhaloMassType[:,4]; el tiempo es el lookback del snapshot descendiente.
MERGER_MASS_FIELD = 'SubhaloMassType'

# ── Fotometría (§2.2): banda r ────────────────────────────────────────────
# GFM_StellarPhotometrics: U, B, V, K (Buser, Vega) ; g, r, i, z (SDSS, AB)
PHOT_BANDS  = ['U', 'B', 'V', 'K', 'g', 'r', 'i', 'z']
BAND_IDX    = 5          # r (SDSS, AB)
BAND_NAME   = 'r'
M_SUN_R_AB  = 4.65       # M☉ en r SDSS (AB; Willmer 2018)
# μ [mag/arcsec²] = M☉ + 21.572 − 2.5 log10(Σ [L☉/pc²])  (M☉ se cancela: μ queda en AB)
SB_CONST    = M_SUN_R_AB + 21.572
# SUPUESTO: μ intrínseco en reposo a z=0, sin atenuación cosmológica (1+z)^4,
# sin corrección K ni extinción por polvo (el paper no menciona ninguna).
MU_HOLMBERG = 26.5       # mag/arcsec² (Holmberg 1958), banda r según el paper
BV_IDX      = (1, 2)     # B y V (Buser, Vega) para el color B−V (§6.1)

# ── Perfil de brillo superficial y radio de Holmberg (§2.2, Fig. 1) ───────
# SUPUESTO: proyección sobre el plano de los ejes mayor e intermedio del tensor
# de inercia reducido (línea de visión = eje menor) y radio circular R = sqrt(x²+y²).
INERTIA_REDUCED = True
# SUPUESTO: anillos duros lineales de 5 kpc hasta R200 (la Fig. 1 tiene eje lineal
# 0–1400 kpc y puntos densos y equiespaciados). 'kernel' = kernel gaussiano corregido.
PROFILE_METHOD  = 'anillos'      # 'anillos' | 'kernel'
PROFILE_BINNING = 'lineal'       # 'lineal' | 'log'
SB_BIN_KPC      = 5.0            # ancho de anillo (binning lineal) [kpc]
SB_N_BINS_LOG   = 60             # número de anillos (binning log)
SB_R_MIN_KPC    = 0.5            # radio interno (binning log) [kpc]
KERNEL_FRAC_SIGMA = 0.15         # σ del kernel como fracción del radio (sólo 'kernel')
# SUPUESTO: R_H = radio donde el AJUSTE μ(R) = a + b·log10(R) cruza 26.5
# (leyenda de la Fig. 1: "Logarithmic fit to 1D surface brightness curve"),
# ajuste sin pesos a todos los anillos con μ finito entre R_FIT_MIN y R200.
# Evidencia (Fig. 1 digitalizada, 2026-10-16): la línea roja (200 kpc) coincide con el cruce
# del AJUSTE (μ_fit ≈ 26.55 a 200 kpc); los puntos medidos están en μ ≈ 26.15 ahí y cruzan
# 26.5 recién hacia ~230 kpc. Por eso el defecto es 'ajuste' (global). Ojo: en perfiles con
# curvatura tipo Sérsic este R_H queda ~15 % por dentro del cruce real del perfil.
RH_METHOD       = 'ajuste'       # 'ajuste' (global, defecto) | 'interp' (cruce del perfil) | 'ajuste_local'
RH_LOCAL_MU_WINDOW = 1.0         # sólo 'ajuste_local': anillos con |μ − 26.5| ≤ 1 mag
R_FIT_MIN_KPC   = 0.0
# Límite externo de la selección BCG+ICL: R200 en 3D (SUPUESTO: 3D, no proyectado).
R_OUTER_FRAC    = 1.0

# ── Gradientes (§6.1) ─────────────────────────────────────────────────────
# SUPUESTO: "15-degree selection along the major axis" = cuña de ±7.5° alrededor
# del eje mayor proyectado, a ambos lados del centro.
WEDGE_DEG     = 15.0
N_BINS_GRAD   = 30           # bins de igual número de partículas hasta R200
Z_SUN         = 0.0127       # fracción metálica solar (Asplund+2009), como en TNG

# ── Progenitores (§2.3, §7) ───────────────────────────────────────────────
M_PROG_THRESH = 1e10         # M☉ (SUPUESTO: misma convención de h que UM_PLOT)
# SUPUESTO: subhalo "anómalo" = SubhaloFlag == 0 en el snapshot de origen.
# SUPUESTO: búsqueda del subhalo de origen dentro del FoF anfitrión del progenitor
# principal del BCG; fuera del FoF sólo si TRACE_GLOBAL_LOOKUP=True (costoso).
TRACE_GLOBAL_LOOKUP = False
TRACE_SNAP_MIN      = 0      # snapshot más antiguo a rastrear
SIGNIFICANT_FRAC    = 0.90   # "significant progenitors": 90 % del material
# Archivo con el rastreo por partícula (03), para no repetir el cálculo
PROG_OUT = os.environ.get('ICL_PROG_OUT', str(REPO_DIR / 'progenitores_mayes2026.hdf5'))

# ── Valores publicados (para tablas de comparación) ───────────────────────
PAPER_VALUES = {
    'N_halos'                    : 127,
    'N_halos_1e14'               : 10,
    'f_icl_rango'                : (0.05, 0.25),
    'f_icl_mediana_aprox'        : 0.10,
    'bcgicl_mergers'             : (0.59, 0.13),
    'bcgicl_insitu'              : (0.28, 0.09),
    'bcgicl_surviving'           : (0.12, 0.10),
    'icl_mergers'                : (0.44, 0.18),
    'icl_insitu'                 : (0.15, 0.06),
    'icl_surviving'              : (0.41, 0.20),
    'bcg_surviving'              : (0.06, 0.08),
    'icl_prog_alta_masa'         : (0.654, 0.149),
    'icl_prog_baja_masa'         : (0.346, 0.149),
    'icl_comparte_con_bcg'       : 0.93,
    'bcg_comparte_con_icl'       : 0.91,
    'icl_top1_progenitor'        : (0.27, 0.12),
}
