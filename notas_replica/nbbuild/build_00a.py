from common import *

cells = [
md(r'''
# 00a · Chequeos baratos en el Lab antes de correr la réplica
Sólo lee metadatos (headers, catálogo de grupos con un campo, claves del catálogo de ensamblaje,
un tramo pequeño de él y una rama principal). **No carga partículas.** Si algo falla, la última celda
lo lista con qué ajustar en `params_icl.py` o con qué variable de entorno.
'''),
code(REPO_CELL),
code(r'''
import importlib.metadata
import numpy as np
import h5py

import illustris_python as il
import params_icl as P
import icl_funcs as F
'''),
code(r'''
fallas = []

def chequeo(ok, texto, ayuda=''):
    """Imprime OK/FALLA y acumula las fallas con su ayuda."""
    print(f"[{'OK   ' if ok else 'FALLA'}] {texto}")
    if not ok:
        fallas.append(f"{texto}  →  {ayuda}")
    return ok
'''),

md('## 1 · illustris_python y rutas'),
code(r'''
try:
    il_version = importlib.metadata.version('illustris_python')
except importlib.metadata.PackageNotFoundError:
    il_version = getattr(il, '__version__', 'sin metadatos')
print(f"illustris_python = {il_version}   ({il.__file__})")
print(f"basePath         = {P.basePath}")
print(f"SA por partícula = {P.SA_STARS_FILE.format(snap=P.SNAP)}")
print(f"SA agregado      = {P.SA_FILE}")
print(f"Catálogo salida  = {P.CATALOG_OUT}")
'''),
code(r'''
chequeo(os.path.isdir(P.basePath), "basePath existe", "define TNG_BASEPATH o corrige P.basePath")
snap_file = il.snapshot.snapPath(P.basePath, P.SNAP)
chequeo(os.path.exists(snap_file), f"snapshot {P.SNAP} (chunk 0) existe", snap_file)
off_file = il.groupcat.offsetPath(P.basePath, P.SNAP)
chequeo(os.path.exists(off_file), f"offsets_{P.SNAP:03d}.hdf5 existe", off_file)
'''),

md(r'''
## 2 · Headers: catálogo de grupos y snapshot
`il.groupcat.loadHeader` lee el Header del catálogo de grupos (sin `MassTable`); `MassTable` y
`NumPart_Total` se leen del Header del snapshot con `F.snapshot_header`.
'''),
code(r'''
Header = il.groupcat.loadHeader(P.basePath, P.SNAP)
CLAVES_GC = ['BoxSize', 'Redshift', 'Time', 'HubbleParam']
faltan_gc = [k for k in CLAVES_GC if k not in Header]
chequeo(not faltan_gc, f"Header del groupcat trae {CLAVES_GC}", f"faltan {faltan_gc}; claves: {sorted(Header)}")
print(f"  (MassTable en el Header del groupcat: {'MassTable' in Header} — la réplica no la usa de ahí)")
chequeo(abs(Header['Redshift']) < 1e-6, f"z = {Header['Redshift']:.3g} en el snap {P.SNAP}", "revisa P.SNAP")
chequeo(abs(Header['BoxSize'] - P.L_BOX) < 1, f"BoxSize = {Header['BoxSize']:.0f} ckpc/h", "¿es TNG100-1?")
chequeo(abs(Header['HubbleParam'] - P.h) < 1e-4, f"h = {Header['HubbleParam']}", "revisa P.h")
'''),
code(r'''
HeaderSnap = F.snapshot_header(P.basePath, P.SNAP)
CLAVES_SNAP = ['MassTable', 'NumPart_Total', 'NumPart_Total_HighWord']
faltan_sn = [k for k in CLAVES_SNAP if k not in HeaderSnap]
chequeo(not faltan_sn, f"Header del snapshot trae {CLAVES_SNAP}", f"faltan {faltan_sn}; claves: {sorted(HeaderSnap)}")
n_star_snap = F.n_stars_snapshot(HeaderSnap)
m_dm = HeaderSnap['MassTable'][1] * P.UM
print(f"PartType4 en el snapshot = {n_star_snap:,}")
chequeo(5e6 < m_dm < 1e7, f"m_DM = {m_dm:.3e} Msun (TNG100-1: 7.5e6)", "MassTable[1] inesperado")
'''),

md('## 3 · Muestra: corte Group_M_Crit200 ≥ 1e3 (1e10 M☉/h)'),
code(r'''
M200_int = il.groupcat.loadHalos(P.basePath, P.SNAP, fields=['Group_M_Crit200'])
n_sel = int((M200_int >= P.M200_MIN_INTERNAL).sum())
n_14  = int((M200_int * P.UM_H > 1e14).sum())
chequeo(n_sel == P.PAPER_VALUES['N_halos'], f"N halos = {n_sel}  (paper {P.PAPER_VALUES['N_halos']})",
        "si no da 127, revisar la convención de h (params_icl.M200_MIN_INTERNAL)")
chequeo(n_14 == P.PAPER_VALUES['N_halos_1e14'], f"N > 1e14 Msun/h = {n_14}  (paper {P.PAPER_VALUES['N_halos_1e14']})")
print(f"M200c máx = {M200_int.max() * P.UM_H:.3e} Msun/h   (paper 2.6e14)")
'''),

md('## 4 · Catálogo de ensamblaje estelar por partícula (Rodriguez-Gomez+2016)'),
code(r'''
SA_PATH = P.SA_STARS_FILE.format(snap=P.SNAP)
sa_ok = chequeo(os.path.exists(SA_PATH), "archivo SA por partícula existe",
                "busca con: ls /home/tnguser/sims.TNG/TNG100-1/postprocessing/StellarAssembly/ "
                "y ajusta P.SA_STARS_FILE / TNG_SA_STARS_FILE")
'''),
code(r'''
if sa_ok:
    with h5py.File(SA_PATH, 'r') as f:
        sa_shapes = {k: (f[k].shape, str(f[k].dtype)) for k in f if isinstance(f[k], h5py.Dataset)}
        for k, (shape, dt) in sa_shapes.items():
            print(f"  {k:<35s} {str(shape):<14s} {dt}")
        for campo in (P.SA_FIELD_INSITU, P.SA_FIELD_ORIGIN):
            if chequeo(campo in f, f"campo '{campo}' presente", f"campos disponibles: {list(sa_shapes)}"):
                chequeo(f[campo].shape[0] == n_star_snap, f"largo de '{campo}' = PartType4 del snapshot",
                        "el catálogo no está en el orden de PartType4: hay que alinear por ParticleIDs")
                muestra = f[campo][:1_000_000]
                print(f"    valores en el primer millón: { {int(v): int(n) for v, n in zip(*np.unique(muestra, return_counts=True))} }")
'''),
code(r'''
# opcional: sólo lo usa read_hdf5_SA.ipynb (no detiene la réplica)
print(f"[{'OK   ' if os.path.exists(P.SA_FILE) else 'aviso'}] catálogo SA agregado: {P.SA_FILE}")
'''),

md('## 5 · Árboles SubLink (rama principal del BCG más masivo)'),
code(r'''
i_max   = int(np.argmax(M200_int))
first   = il.groupcat.loadSingle(P.basePath, P.SNAP, haloID=i_max)['GroupFirstSub']
try:
    mpb = il.sublink.loadTree(P.basePath, P.SNAP, int(first), onlyMPB=True,
                              fields=['SnapNum', 'SubfindID', 'SubhaloGrNr', 'SubhaloPos',
                                      'SubhaloMassInRadType', 'MainLeafProgenitorID'])
    chequeo(mpb is not None, f"MPB del BCG {first}: {len(mpb['SnapNum'])} snapshots, campos de 03 presentes")
except (OSError, KeyError) as e:
    chequeo(False, "árbol SubLink con los campos de 03", str(e))
'''),

md('## 6 · Resumen'),
code(r'''
if fallas:
    print("Hay que resolver antes de correr la réplica:")
    for x in fallas:
        print(f"  - {x}")
    raise RuntimeError(f"{len(fallas)} chequeo(s) fallaron")
print("Todo OK: se puede correr 00_catalogo.ipynb")
'''),
]
if __name__ == "__main__":
    write_nb('00a_chequeos_lab.ipynb', cells)
