"""
Copia local de params.py adaptada para el análisis ICL/BCG (Mayes+2026).
Modifica los valores según tu configuración del cluster.
"""
import numpy as np
import sys, os

# ── Ruta a los datos de TNG en el cluster ─────────────────────────────────
basePath = '/home/tnguser/sims.TNG/TNG100-1/output/'

# ── Snapshot a analizar (99 = z=0) ────────────────────────────────────────
SNAP = 99

# ── Parámetros cosmo/simulación ───────────────────────────────────────────
h       = 0.6774            # parámetro de Hubble (TNG)
L_BOX   = 75000.0           # kpc/h (= 110.7 Mpc * h)

# ── Factores de conversión de unidades internas → físicas ─────────────────
# Posición: ckpc/h  →  kpc físicos = value * a / h  (a=1 en snap 99)
# Masa:     1e10 M☉/h →  M☉ = value * 1e10 / h
UL = 1.0 / h      # factor longitud (a=1)
UM = 1e10 / h     # factor masa

# ── Archivo del catálogo con los group IDs a analizar ─────────────────────
# (Generado por tu código de estado dinámico; contiene los IDs de grupos TNG)
CATALOG_PKL = './catalogo_grupos.pkl'   # ← ajustar path según corresponda

# ── Archivo de salida del catálogo ICL ────────────────────────────────────
CATALOG_OUT = './catalogo_icl.hdf5'

# ── Campos del catálogo de halos que queremos extraer ─────────────────────
HALOS_FIELDS = [
    'GroupFirstSub',
    'Group_M_Crit200',
    'Group_R_Crit200',
    'GroupPos',
    'GroupCM',
    'GroupNsubs',
    'GroupMassType'
]

SUBHALOS_FIELDS = [
    'SubhaloPos',
    'SubhaloMassType',
    'SubhaloGrNr',
    'SubhaloFlag',
    'SubhaloStellarPhotometrics'
]

# ── Parámetros de la separación BCG/ICL ───────────────────────────────────
MU_HOLMBERG  = 26.5    # mag/arcsec² en banda r (corte de Holmberg)
M_SUN_R_AB   = 4.65    # Magnitud absoluta del Sol en r-band (AB)
SB_CONST     = M_SUN_R_AB + 21.572   # constante μ = SB_CONST - 2.5 log10(Σ_L)

# ── Parámetros del análisis de progenitores ───────────────────────────────
M_PROG_THRESH = 1e10   # M☉ — umbral alta/baja masa progenitor
Z_SUN         = 0.0127  # fracción metálica solar (Asplund+2009)

# ── Ruta a los catálogos de ensamblaje estelar (Rodriguez-Gomez+2016) ─────
SA_FILE = '../stellar_assembly.hdf5'   # ← ajustar
CAT_FOSSIL_50 = '../fossil_catalog_TNG50.hdf5'        # ← ajustar
