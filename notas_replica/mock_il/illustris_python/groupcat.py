"""groupcat simulado: misma firma que illustris_python.groupcat para las funciones que usa la réplica."""
import os
import numpy as np
import h5py
from . import _world as W


def _pick(d, fields):
    if fields is None:
        return d
    if isinstance(fields, str):
        fields = [fields]
    missing = [f for f in fields if f not in d]
    if missing:
        raise KeyError(f"[mock] campos no simulados: {missing}")
    if len(fields) == 1:
        return d[fields[0]]
    out = {f: d[f] for f in fields}
    out['count'] = d['count']
    return out


def loadHeader(basePath, snapNum):
    z = W.z_of(snapNum)
    # Fiel a la doc de TNG (data/docs/specifications, Header de fof_subhalo_tab): Ngroups_*, Nsubgroups_*,
    # NumFiles, Time, Redshift. BoxSize no está en esa tabla pero SÍ en los archivos reales (outputs del
    # Lab de Mary: BoxSize = 75000). Sin MassTable, NumPart_* ni HubbleParam.
    return {'Redshift': z, 'Time': 1 / (1 + z), 'BoxSize': W.BOX,
            'Ngroups_ThisFile': W.NG, 'Ngroups_Total': W.NG,
            'Nsubgroups_ThisFile': W.NG * W.NSUB, 'Nsubgroups_Total': W.NG * W.NSUB, 'NumFiles': 1}


def loadHalos(basePath, snapNum, fields=None):
    d = {'count': W.NG,
         'Group_M_Crit200': W.M_INT.astype(np.float32),
         'Group_R_Crit200': W.R_INT.astype(np.float32),
         'GroupPos': W.POS.astype(np.float32),
         'GroupFirstSub': np.array([W.first_sub(g) for g in range(W.NG)], dtype=np.int32),
         'GroupNsubs': np.full(W.NG, W.NSUB, dtype=np.int32)}
    return _pick(d, fields)


def loadSubhalos(basePath, snapNum, fields=None):
    n = W.NG * W.NSUB
    pos = np.zeros((n, 3)); lent = np.zeros((n, 6), dtype=np.int32); grnr = np.zeros(n, dtype=np.int32)
    for g in range(W.NG):
        st = W.stars_at(g, snapNum)
        for k in range(W.NSUB):
            i = W.first_sub(g) + k
            pos[i] = W.POS[g] if k == 0 else W.sat_pos(g, k, snapNum)
            lent[i, 4] = st['lens'][k]
            grnr[i] = g
    mtype = np.zeros((n, 6), dtype=np.float32); mtype[:, 4] = lent[:, 4] * W.M_STAR
    hmr = np.zeros((n, 6), dtype=np.float32); hmr[:, 4] = 10.0
    d = {'count': n, 'SubhaloPos': pos.astype(np.float32), 'SubhaloLenType': lent,
         'SubhaloMassType': mtype, 'SubhaloFlag': W.sub_flag(snapNum), 'SubhaloGrNr': grnr,
         'SubhaloMass': mtype.sum(1) * 10, 'SubhaloHalfmassRadType': hmr}
    return _pick(d, fields)


def loadSingle(basePath, snapNum, haloID=-1, subhaloID=-1):
    if haloID >= 0:
        h = loadHalos(basePath, snapNum)
        return {k: v[haloID] for k, v in h.items() if k != 'count'}
    s = loadSubhalos(basePath, snapNum)
    return {k: v[subhaloID] for k, v in s.items() if k != 'count'}


def offsetPath(basePath, snapNum):
    W.ensure_files()
    path = os.path.join(W.MDIR, f'offsets_{snapNum:03d}.hdf5')
    if not os.path.exists(path):
        sub_off, _, grp_off, _ = W.global_offsets(snapNum)
        with h5py.File(path, 'w') as f:
            sb = np.zeros((len(sub_off), 6), dtype=np.int64); sb[:, 4] = sub_off
            gb = np.zeros((len(grp_off), 6), dtype=np.int64); gb[:, 4] = grp_off
            f.create_dataset('Subhalo/SnapByType', data=sb)
            f.create_dataset('Group/SnapByType', data=gb)
    return path
