from common import *

cells = [
md(r'''
# 03 · Progenitores del ICL y del BCG
**Metodología:** §2.3 y §7 de Mayes+2026 · Figs. 13–15 y tabla comparativa

Rastreo por partícula de las estrellas **ex situ** del BCG+ICL (subhalo central, r₃D ≤ R₂₀₀c, z = 0):

1. Hacia atrás en la rama principal (MPB) del BCG: snapshot en que cada partícula **deja de estar en el central**.
2. Subhalo que contiene a la partícula en ese snapshot (muestra 1).
3. Muestra 2: último snapshot en que la partícula está **fuera de R₂₀₀** del anfitrión, y su subhalo.
   Si el subhalo de la muestra 1 es anómalo, se reemplaza por el de la muestra 2; si ambos lo son, se excluye.
4. Masa del progenitor = máximo de la masa estelar dentro de 2×HMR (`SubhaloMassInRadType`) a lo largo de su MPB.
   Identidad del progenitor = `MainLeafProgenitorID` de su rama principal.

BCG / ICL con el mismo R_H y la misma rotación de `01_`. Las componentes vienen del mismo catálogo SA por partícula que `02_`.

`# SUPUESTO:` anómalo = `SubhaloFlag == 0`; búsqueda del subhalo dentro del FoF anfitrión del MPB
(fuera de él sólo con `P.TRACE_GLOBAL_LOOKUP = True`); una partícula que no está en el FoF anfitrión
cuenta como "fuera de R₂₀₀"; partículas que eran gas en el snapshot de despegue quedan sin progenitor.

**Costo:** recorre hasta 99 snapshots × 127 cúmulos en el Lab (horas). El resultado se guarda en
`P.PROG_OUT` y se puede recargar sin repetir el rastreo. Probar primero con `N_PRUEBA = 2`.
'''),
code(REPO_CELL),
code(r'''
from astropy.cosmology import FlatLambdaCDM

import time
import numpy as np
import h5py
import matplotlib.pyplot as plt
import seaborn as sns

import illustris_python as il
import params_icl as P
import icl_funcs as F
'''),
code(STYLE_CELL),
code(fig_dirs_cell('03_analisis_progenitores')),
code(COSMO_CELL),

md('## 1 · Parámetros'),
code(r'''
N_PRUEBA   = None    # p. ej. 2 para una corrida corta de validación; None = toda la muestra
RECALCULAR = False   # True = rehacer el rastreo aunque exista P.PROG_OUT

print(f"Archivo de rastreo : {P.PROG_OUT}  (existe: {os.path.exists(P.PROG_OUT)})")
print(f"Snapshots          : {P.SNAP - 1} → {P.TRACE_SNAP_MIN}")
print(f"Búsqueda global    : {P.TRACE_GLOBAL_LOOKUP}")
print(f"Umbral progenitor  : {P.M_PROG_THRESH:.0e} (unidades de P.UM_PLOT)")
'''),

md('## 2 · Cargar catálogo (salidas de `00`, `01` y `02`)'),
code(r'''
with h5py.File(P.CATALOG_OUT, 'r') as f:
    group_idx   = f['group_idx'][:]
    M200c       = f['M200c_Msun'][:]
    R200c       = f['R200c_kpc'][:]
    GroupPos    = f['GroupPos_kpc'][:]
    bcg_sub_idx = f['bcg_sub_idx'][:]
    dyn_state   = f['dyn_state'][:]
    R_H         = f['R_H_kpc'][:]
    R_rot       = f['R_rot'][:]
    f_comp_bcgicl = f['f_comp_bcgicl'][:] if 'f_comp_bcgicl' in f else None
    f_comp_icl    = f['f_comp_icl'][:] if 'f_comp_icl' in f else None
    f_comp_bcg    = f['f_comp_bcg'][:] if 'f_comp_bcg' in f else None

n_cl     = len(group_idx)
clusters = [i for i in range(n_cl) if np.isfinite(R_H[i])][:N_PRUEBA]
Header   = il.groupcat.loadHeader(P.basePath, P.SNAP)
box_kpc  = Header['BoxSize'] * P.UL
print(f"Cúmulos a procesar: {len(clusters)} de {n_cl}")
'''),
code(MASS_LABEL_CELL),
code(r'''
HACER_RASTREO = RECALCULAR or not os.path.exists(P.PROG_OUT)
print(f"Rastreo nuevo: {HACER_RASTREO}")
'''),

md('## 3 · Estrellas ex situ del BCG+ICL a z = 0'),
code(r'''
def ex_situ_z0(i):
    """IDs, masa [Msun], componente (1 mergers, 2 sobrev.), r3D [kpc] y región ICL de las estrellas ex situ."""
    st = F.load_central_stars(P.basePath, P.SNAP, bcg_sub_idx[i], GroupPos[i], R200c[i], box_kpc,
                              P.UL, P.UM, extra_fields=['ParticleIDs'], r_outer_frac=P.R_OUTER_FRAC)
    insitu, origin = F.load_sa_stars(P.SA_STARS_FILE, P.SNAP, st['offset'], st['length'],
                                     P.SA_FIELD_INSITU, P.SA_FIELD_ORIGIN)
    comp = F.classify_components(insitu[st['idx_sub']], origin[st['idx_sub']])
    R    = F.projected_radius(st['pos'] @ R_rot[i].T)
    ex   = (comp == 1) | (comp == 2)
    return {'ids': st['ParticleIDs'][ex].astype(np.int64), 'mass': st['mass'][ex], 'comp': comp[ex],
            'r3d': st['r3d'][ex], 'is_icl': R[ex] > R_H[i]}
'''),
code(r'''
z0 = {}
if HACER_RASTREO:
    for i in clusters:
        z0[i] = ex_situ_z0(i)
    print(f"Partículas ex situ a rastrear: {sum(len(v['ids']) for v in z0.values()):,}")
'''),

md('## 4 · Rama principal (MPB) de cada BCG y redshift de cada snapshot'),
code(r'''
z_snap = {}
for s in range(P.SNAP + 1):
    try:
        z_snap[s] = il.groupcat.loadHeader(P.basePath, s)['Redshift']
    except (OSError, KeyError):
        z_snap[s] = np.nan
'''),
code(r'''
mpb = {}
if HACER_RASTREO:
    for i in clusters:
        t = il.sublink.loadTree(P.basePath, P.SNAP, int(bcg_sub_idx[i]),
                                fields=['SnapNum', 'SubfindID', 'SubhaloGrNr', 'SubhaloPos'], onlyMPB=True)
        mpb[i] = {int(s): (int(sf), int(g), p) for s, sf, g, p in
                  zip(t['SnapNum'], t['SubfindID'], t['SubhaloGrNr'], t['SubhaloPos'])}
    print(f"Largo medio del MPB: {np.mean([len(v) for v in mpb.values()]):.0f} snapshots")
'''),

md(r'''
## 5 · Rastreo hacia atrás (pasos 1–3 de §2.3)
Un loop externo por snapshot (se carga una sola vez `SubhaloLenType` y `SubhaloFlag`) e interno por cúmulo.
'''),
code(r'''
def fof_star_table(s, grnr, lens_all):
    """IDs ordenados, subhalo (−1 = fuzz) y posición [ckpc/h] de las estrellas del FoF grnr en el snapshot s."""
    g  = il.groupcat.loadSingle(P.basePath, s, haloID=int(grnr))
    st = il.snapshot.loadHalo(P.basePath, s, int(grnr), 'stars', fields=['ParticleIDs', 'Coordinates'])
    if not isinstance(st, dict) or st.get('count', 0) == 0:
        return g, None
    first = int(g['GroupFirstSub'])
    lens  = lens_all[first:first + int(g['GroupNsubs']), 4] if first >= 0 else np.array([], int)
    sub   = F.subhalo_of_local_index(np.arange(st['count']), lens, max(first, 0))
    order = np.argsort(st['ParticleIDs'])
    return g, {'ids': st['ParticleIDs'][order].astype(np.int64), 'sub': sub[order],
               'pos': st['Coordinates'][order]}


def lookup(table, ids_q):
    """Posición en la tabla ordenada de cada ID pedido y si fue encontrado."""
    if table is None or len(table['ids']) == 0:
        return np.zeros(len(ids_q), int), np.zeros(len(ids_q), bool)
    k = np.clip(np.searchsorted(table['ids'], ids_q), 0, len(table['ids']) - 1)
    return k, table['ids'][k] == ids_q


_global_cache = {}
def global_subhalo(s, ids_q, lens_all):
    """Subhalo (−1 = no ligado / no encontrado) buscando en TODAS las estrellas del snapshot (costoso)."""
    if s not in _global_cache:
        _global_cache.clear()
        ids_all = il.snapshot.loadSubset(P.basePath, s, 'stars', fields=['ParticleIDs']).astype(np.int64)
        with h5py.File(il.groupcat.offsetPath(P.basePath, s), 'r') as f:
            off = f['Subhalo/SnapByType'][:, 4].astype(np.int64)
        order = np.argsort(ids_all)
        _global_cache[s] = (ids_all[order], order, off)
    ids_sorted, order, off = _global_cache[s]
    k = np.clip(np.searchsorted(ids_sorted, ids_q), 0, len(ids_sorted) - 1)
    found = ids_sorted[k] == ids_q
    gidx = order[k]
    sub = np.searchsorted(off, gidx, side='right') - 1
    inside = found & (sub >= 0) & (gidx < off[np.maximum(sub, 0)] + lens_all[np.maximum(sub, 0), 4])
    return np.where(inside, sub, -1)
'''),
code(r'''
trk = {}
if HACER_RASTREO:
    for i in clusters:
        n = len(z0[i]['ids'])
        trk[i] = {k: np.full(n, -1, dtype=np.int64) for k in
                  ['det_snap', 'det_sub', 'det_flag', 'out_snap', 'out_sub', 'out_flag']}
        trk[i]['det_done'] = np.zeros(n, dtype=bool)
        trk[i]['out_done'] = np.zeros(n, dtype=bool)
'''),
code(r'''
if HACER_RASTREO:
    t0 = time.time()
    for s in range(P.SNAP - 1, P.TRACE_SNAP_MIN - 1, -1):
        activos = [i for i in clusters if s in mpb[i]
                   and not (trk[i]['det_done'].all() and trk[i]['out_done'].all())]
        if not activos:
            continue
        sub_cat  = il.groupcat.loadSubhalos(P.basePath, s, fields=['SubhaloLenType', 'SubhaloFlag'])
        lens_all = sub_cat['SubhaloLenType']
        flag_all = sub_cat['SubhaloFlag'].astype(int)
        a_s      = 1.0 / (1.0 + z_snap[s])
        box_s    = Header['BoxSize'] * a_s / P.h

        for i in activos:
            sf, grnr, pos_cen = mpb[i][s]
            S, ids = trk[i], z0[i]['ids']
            g, tab = fof_star_table(s, grnr, lens_all)

            # paso 1: despegue del central
            ids_cen = il.snapshot.loadSubhalo(P.basePath, s, sf, 'stars', fields=['ParticleIDs'])
            ids_cen = ids_cen.astype(np.int64) if isinstance(ids_cen, np.ndarray) else np.array([], np.int64)
            new = F.detach_update(ids, S['det_done'], ids_cen, s, S['det_snap'])
            if new.any():
                k, found = lookup(tab, ids[new])
                sub = np.where(found, tab['sub'][k] if tab is not None else -1, -1)
                if P.TRACE_GLOBAL_LOOKUP and (~found).any():
                    sub[~found] = global_subhalo(s, ids[new][~found], lens_all)
                S['det_sub'][new]  = sub
                S['det_flag'][new] = np.where(sub >= 0, flag_all[np.maximum(sub, 0)], -1)

            # muestra 2: último snapshot fuera de R200 del anfitrión
            pend = np.nonzero(~S['out_done'])[0]
            if len(pend):
                k, found = lookup(tab, ids[pend])
                r = np.full(len(pend), np.inf)       # SUPUESTO: fuera del FoF ⇒ fuera de R200
                if found.any():
                    d = F.periodic_delta(tab['pos'][k[found]] * a_s / P.h, pos_cen * a_s / P.h, box_s)
                    r[found] = np.linalg.norm(d, axis=1)
                out = r > g['Group_R_Crit200'] * a_s / P.h
                idx = pend[out]
                sub = np.where(found[out], tab['sub'][k[out]] if tab is not None else -1, -1)
                if P.TRACE_GLOBAL_LOOKUP and (~found[out]).any():
                    sub[~found[out]] = global_subhalo(s, ids[idx][~found[out]], lens_all)
                S['out_done'][idx] = True
                S['out_snap'][idx] = s
                S['out_sub'][idx]  = sub
                S['out_flag'][idx] = np.where(sub >= 0, flag_all[np.maximum(sub, 0)], -1)

        n_pend = sum((~trk[i]['det_done']).sum() + (~trk[i]['out_done']).sum() for i in clusters)
        print(f"  snap {s:2d}  activos {len(activos):3d}  pendientes {n_pend:,}  ({(time.time()-t0)/60:.1f} min)", end="\r")
    print("\nRastreo terminado.")
'''),

md('## 6 · Elección del subhalo de origen, masa e identidad del progenitor (paso 3)'),
code(r'''
_prog_cache = {}
def progenitor_info(s, sub):
    """(MainLeafProgenitorID, masa estelar máx. en 2×HMR [Msun]) a lo largo del MPB del subhalo (s, sub)."""
    key = (int(s), int(sub))
    if key not in _prog_cache:
        t = il.sublink.loadTree(P.basePath, key[0], key[1], onlyMPB=True,
                                fields=['SubhaloMassInRadType', 'MainLeafProgenitorID'])
        _prog_cache[key] = ((-1, np.nan) if t is None else
                            (int(t['MainLeafProgenitorID'][0]), float(t['SubhaloMassInRadType'][:, 4].max() * P.UM)))
    return _prog_cache[key]
'''),
code(r'''
if HACER_RASTREO:
    for i in clusters:
        S = trk[i]
        use1 = (S['det_sub'] >= 0) & (S['det_flag'] == 1)
        use2 = ~use1 & (S['out_sub'] >= 0) & (S['out_flag'] == 1)
        S['origen']    = np.where(use1, 1, np.where(use2, 2, 0))      # 0 = excluida
        S['prog_id']   = np.full(len(use1), -1, dtype=np.int64)
        S['prog_mass'] = np.full(len(use1), np.nan)
        snaps = np.where(use1, S['det_snap'], S['out_snap'])
        subs  = np.where(use1, S['det_sub'], S['out_sub'])
        sel   = S['origen'] > 0
        if not sel.any():
            continue
        pares, inv = np.unique(np.stack([snaps[sel], subs[sel]], axis=1), axis=0, return_inverse=True)
        info = [progenitor_info(s, sub) for s, sub in pares]
        S['prog_id'][sel]   = np.array([x[0] for x in info], dtype=np.int64)[inv.ravel()]
        S['prog_mass'][sel] = np.array([x[1] for x in info])[inv.ravel()]
    print(f"Progenitores distintos consultados: {len(_prog_cache):,}")
'''),
code(r'''
if HACER_RASTREO:
    with h5py.File(P.PROG_OUT, 'w') as f:
        f.attrs['muestra'] = P.CATALOG_OUT
        for i in clusters:
            g = f.create_group(f'c{i:03d}')
            for k, v in z0[i].items():
                g.create_dataset(k, data=v)
            for k in ['det_snap', 'det_sub', 'det_flag', 'out_snap', 'out_sub', 'out_flag',
                      'origen', 'prog_id', 'prog_mass']:
                g.create_dataset(k, data=trk[i][k])
    print(f"Rastreo guardado en {P.PROG_OUT}")
'''),
code(r'''
# Recarga (también sirve para re-analizar sin rastrear de nuevo)
prog = {}
with h5py.File(P.PROG_OUT, 'r') as f:
    for name in f:
        prog[int(name[1:])] = {k: f[name][k][:] for k in f[name]}
clusters = sorted(prog)
print(f"Cúmulos con rastreo: {len(clusters)}")
'''),

md('## 7 · Chequeos del rastreo'),
code(r'''
n_tot   = sum(len(prog[i]['ids']) for i in clusters)
n_anom1 = sum(((prog[i]['det_sub'] >= 0) & (prog[i]['det_flag'] == 0)).sum() for i in clusters)
n_sub1  = sum((prog[i]['det_sub'] >= 0).sum() for i in clusters)
n_o1    = sum((prog[i]['origen'] == 1).sum() for i in clusters)
n_o2    = sum((prog[i]['origen'] == 2).sum() for i in clusters)
n_exc   = sum((prog[i]['origen'] == 0).sum() for i in clusters)

print(f"Partículas ex situ             : {n_tot:,}")
print(f"Muestra 1 con subhalo          : {n_sub1 / n_tot:.1%}")
print(f"  anómalos en muestra 1        : {n_anom1 / max(n_sub1, 1):.2%}   (paper: 1–2 %)")
print(f"Origen muestra 1 / 2 / excluida: {n_o1 / n_tot:.1%} / {n_o2 / n_tot:.1%} / {n_exc / n_tot:.1%}")
print(f"Progenitores con M < 1e7 Msun  : "
      f"{np.mean(np.concatenate([prog[i]['prog_mass'][prog[i]['origen'] > 0] for i in clusters]) < 1e7):.2%}"
      f"   (paper: ~0.1 %)")
'''),

md('## 8 · §7.1 – masa media de los progenitores (Fig. 13)'),
code(r'''
def m_conv(m_msun):
    """Masa de progenitor en la convención de P.UM_PLOT (para el umbral P.M_PROG_THRESH)."""
    return m_msun / P.UM * P.UM_PLOT

n_c = len(clusters)
lmean_icl, lmean_bcg = np.full(n_c, np.nan), np.full(n_c, np.nan)
f_hi_icl,  f_hi_bcg  = np.full(n_c, np.nan), np.full(n_c, np.nan)

for j, i in enumerate(clusters):
    d  = prog[i]
    ok = d['origen'] > 0
    mp = m_conv(d['prog_mass'])
    for reg, lm, fh in [(d['is_icl'], lmean_icl, f_hi_icl), (~d['is_icl'], lmean_bcg, f_hi_bcg)]:
        m = ok & reg
        if m.sum() == 0:
            continue
        lm[j] = np.log10(np.mean(mp[m]))                       # media pesada por número de partículas
        fh[j] = d['mass'][m & (mp > P.M_PROG_THRESH)].sum() / d['mass'][m].sum()

lM_c = lM200[clusters]
print(f"ICL: fracción de masa de progenitores de alta masa = {np.nanmean(f_hi_icl):.3f} ± {np.nanstd(f_hi_icl):.3f}"
      f"   (paper: {P.PAPER_VALUES['icl_prog_alta_masa'][0]} ± {P.PAPER_VALUES['icl_prog_alta_masa'][1]})")
'''),
code(r'''
b13, eb13, c13 = F.linfit_beta(lmean_bcg, lmean_icl)
xx = np.linspace(np.nanmin(lmean_bcg), np.nanmax(lmean_bcg), 50)

fig, ax = plt.subplots(figsize=(7, 6))
sc = ax.scatter(lmean_bcg, lmean_icl, c=lM_c, cmap='viridis', s=30, edgecolor='k', linewidth=0.3)
plt.colorbar(sc, ax=ax, label=rf'$\log_{{10}}(M_{{200}}\,/\,{MASS_UNIT})$')
ax.plot(xx, xx, 'k--', lw=1)
ax.plot(xx, b13 * xx + c13, color=my_colors[3], lw=2)
ax.text(0.95, 0.04, rf'$\beta$ = {b13:.2f} ± {eb13:.2f}', transform=ax.transAxes, ha='right')
ax.set_xlabel(r'$\log_{10}\langle M_{*,\rm prog}\rangle$ BCG')
ax.set_ylabel(r'$\log_{10}\langle M_{*,\rm prog}\rangle$ ICL')

plt.savefig(f'{FIG_PDF}/fig13_masa_progenitores.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig13_masa_progenitores.png', dpi=150, bbox_inches='tight')
plt.show()
'''),

md('## 9 · §7.2 – progenitores compartidos entre ICL y BCG (Fig. 14)'),
code(r'''
shared = {k: np.full((n_c, 2), np.nan) for k in ['todos', 'mergers', 'sobrevivientes']}
for j, i in enumerate(clusters):
    d  = prog[i]
    pid = np.where(d['origen'] > 0, d['prog_id'], -1)
    for key, cm in [('todos', d['comp'] > 0), ('mergers', d['comp'] == 1), ('sobrevivientes', d['comp'] == 2)]:
        icl, bcg = pid[cm & d['is_icl']], pid[cm & ~d['is_icl']]
        shared[key][j] = [F.shared_fraction(icl, bcg), F.shared_fraction(bcg, icl)]

print(f"ICL comparte con BCG: {np.nanmean(shared['todos'][:, 0]):.2f}   (paper: {P.PAPER_VALUES['icl_comparte_con_bcg']})")
print(f"BCG comparte con ICL: {np.nanmean(shared['todos'][:, 1]):.2f}   (paper: {P.PAPER_VALUES['bcg_comparte_con_icl']})")
print(f"Cúmulos con > 90 % compartido: ICL {(shared['todos'][:, 0] > 0.9).sum()}/{n_c}   "
      f"BCG {(shared['todos'][:, 1] > 0.9).sum()}/{n_c}   (paper: 108/127 y 91/127)")
'''),
code(r'''
bins = np.linspace(0, 1, 26)
fig, ax = plt.subplots(3, 1, figsize=(7, 10), sharex=True)
for a, key in zip(ax, ['todos', 'mergers', 'sobrevivientes']):
    a.hist(shared[key][:, 0], bins=bins, histtype='step', lw=2, color=my_colors[0], label='ICL comparte con BCG')
    a.hist(shared[key][:, 1], bins=bins, histtype='step', lw=2, color=my_colors[3], label='BCG comparte con ICL')
    a.set_ylabel('N cúmulos')
    a.set_title(key)
ax[0].legend(loc='upper left')
ax[2].set_xlabel('Fracción de partículas con progenitor compartido')

plt.savefig(f'{FIG_PDF}/fig14_progenitores_compartidos.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig14_progenitores_compartidos.png', dpi=150, bbox_inches='tight')
plt.show()
'''),

md('## 10 · §7.3 – dónde queda el material de progenitores de baja y alta masa (Fig. 15)'),
code(r'''
hmr_prog = {k: np.full((n_c, 2), np.nan) for k in ['todos', 'mergers', 'sobrevivientes']}   # [baja, alta] kpc
for j, i in enumerate(clusters):
    d  = prog[i]
    ok = d['origen'] > 0
    hi = m_conv(d['prog_mass']) > P.M_PROG_THRESH
    for key, cm in [('todos', d['comp'] > 0), ('mergers', d['comp'] == 1), ('sobrevivientes', d['comp'] == 2)]:
        for col, msk in [(0, ok & cm & ~hi), (1, ok & cm & hi)]:
            hmr_prog[key][j, col] = F.half_mass_radius(d['r3d'][msk], d['mass'][msk])
'''),
code(r'''
def scatter_beta(ax, x, y, color, label):
    b, eb, c = F.linfit_beta(x, y)
    ax.scatter(x, y, s=18, color=color, edgecolor='k', linewidth=0.2,
               label=rf'{label}: $\beta$ = {b:.2f} ± {eb:.2f}')
    xx = np.linspace(np.nanmin(x), np.nanmax(x), 50)
    ax.plot(xx, b * xx + c, color=color, lw=1.5)

fig, ax = plt.subplots(3, 1, figsize=(7, 12), sharex=True)
for a, key in zip(ax, ['todos', 'mergers', 'sobrevivientes']):
    scatter_beta(a, lM_c, np.log10(hmr_prog[key][:, 0]), my_colors[1], 'baja masa')
    scatter_beta(a, lM_c, np.log10(hmr_prog[key][:, 1]), my_colors[2], 'alta masa')
    a.set_ylabel(r'$\log_{10}(r_{1/2}$ / kpc)')
    a.set_title(key)
    a.legend(fontsize=9)
ax[2].set_xlabel(rf'$\log_{{10}}(M_{{200}}\,/\,{MASS_UNIT})$')

plt.savefig(f'{FIG_PDF}/fig15_hmr_progenitores.pdf', bbox_inches='tight')
plt.savefig(f'{FIG_PNG}/fig15_hmr_progenitores.png', dpi=150, bbox_inches='tight')
plt.show()
'''),

md('## 11 · §7.4 – progenitor más significativo'),
code(r'''
top_icl, top_bcg = np.full(n_c, np.nan), np.full(n_c, np.nan)
n90_icl, n90_bcg = np.full(n_c, np.nan), np.full(n_c, np.nan)
same_top = np.zeros(n_c, dtype=bool)
for j, i in enumerate(clusters):
    d   = prog[i]
    pid = np.where(d['origen'] > 0, d['prog_id'], -1)
    id_i, top_icl[j], n90_icl[j] = F.significant_progenitors(pid[d['is_icl']], d['mass'][d['is_icl']], P.SIGNIFICANT_FRAC)
    id_b, top_bcg[j], n90_bcg[j] = F.significant_progenitors(pid[~d['is_icl']], d['mass'][~d['is_icl']], P.SIGNIFICANT_FRAC)
    same_top[j] = (id_i >= 0) and (id_i == id_b)

print(f"Contribución del progenitor más significativo al ICL acretado: "
      f"{np.nanmean(top_icl):.2f} ± {np.nanstd(top_icl):.2f}   (paper: 0.27 ± 0.12)")
print(f"ICL y BCG comparten el progenitor más significativo: {same_top.sum()}/{n_c}")
print(f"N progenitores para el 90 %: ICL {np.nanmedian(n90_icl):.0f}   BCG {np.nanmedian(n90_bcg):.0f}")
'''),

md('## 12 · Tabla comparativa con Mayes+2026 (media ± 1σ)'),
code(r'''
pv = P.PAPER_VALUES
filas = []
if f_comp_bcgicl is not None:
    for reg, arr, keys in [('BCG+ICL', f_comp_bcgicl, ['bcgicl_insitu', 'bcgicl_mergers', 'bcgicl_surviving']),
                           ('ICL',     f_comp_icl,    ['icl_insitu', 'icl_mergers', 'icl_surviving'])]:
        for k in [0, 1, 2]:
            filas.append((f"{reg} {F.COMP_LABELS[k]}", np.nanmean(arr[:, k]), np.nanstd(arr[:, k]), pv[keys[k]]))
    filas.append(('BCG Sobrevivientes', np.nanmean(f_comp_bcg[:, 2]), np.nanstd(f_comp_bcg[:, 2]), pv['bcg_surviving']))
else:
    print("Faltan las fracciones de componentes: corre 02_ primero")
filas += [('ICL prog. alta masa', np.nanmean(f_hi_icl), np.nanstd(f_hi_icl), pv['icl_prog_alta_masa']),
          ('ICL prog. baja masa', np.nanmean(1 - f_hi_icl), np.nanstd(1 - f_hi_icl), pv['icl_prog_baja_masa']),
          ('ICL comparte con BCG', np.nanmean(shared['todos'][:, 0]), np.nanstd(shared['todos'][:, 0]),
           (pv['icl_comparte_con_bcg'], np.nan)),
          ('BCG comparte con ICL', np.nanmean(shared['todos'][:, 1]), np.nanstd(shared['todos'][:, 1]),
           (pv['bcg_comparte_con_icl'], np.nan)),
          ('ICL progenitor top', np.nanmean(top_icl), np.nanstd(top_icl), pv['icl_top1_progenitor'])]

print(f"{'cantidad':<32} {'este análisis':>16} {'Mayes+2026':>16}")
for name, m, s, ref in filas:
    ref_txt = f"{ref[0]:.2f}" + (f" ± {ref[1]:.2f}" if np.isfinite(ref[1]) else '')
    print(f"{name:<32} {m:>8.2f} ± {s:<5.2f} {ref_txt:>16}")
'''),
]
if __name__ == "__main__":
    write_nb('03_analisis_progenitores.ipynb', cells)
