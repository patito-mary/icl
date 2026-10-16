"""snapshot simulado: loadSubhalo / loadHalo / loadSubset / getSnapOffsets / snapPath."""
import os
import numpy as np
from . import _world as W
from .groupcat import _pick

STAR_MAP = {'Coordinates': 'pos', 'Masses': 'mass', 'GFM_StellarPhotometrics': 'phot',
            'GFM_StellarFormationTime': 'sft', 'GFM_Metallicity': 'Z', 'ParticleIDs': 'ids'}


def snapPath(basePath, snapNum, chunkNum=0):
    W.ensure_files()
    return os.path.join(W.MDIR, f'snap_{snapNum:03d}.{chunkNum}.hdf5')


def _stars_dict(st, sel):
    d = {k: st[v][sel] for k, v in STAR_MAP.items()}
    d['count'] = int(sel.sum())
    return d


def _gas(g):
    rng = np.random.default_rng(500 + g)
    R200 = W.R_INT[g]
    v = np.array([1.0, 0.4, -0.2]); v /= np.linalg.norm(v)
    cen = W.POS[g] + W.DYN_OFF[g] * R200 * v
    pos = (cen + rng.normal(scale=0.12 * R200, size=(6000, 3))) % W.BOX
    return {'Coordinates': pos, 'Masses': np.full(6000, 1e-3), 'count': 6000}


def _dm(g):
    rng = np.random.default_rng(700 + g)
    c, R200 = W.C_NFW[g], W.R_INT[g]
    x = np.linspace(1e-4, c, 200000)
    mu = np.log(1 + x) - x / (1 + x)
    r = R200 / c * np.interp(rng.random(60000), mu / mu[-1], x)
    u = rng.normal(size=(len(r), 3)); u /= np.linalg.norm(u, axis=1)[:, None]
    return {'Coordinates': (W.POS[g] + u * r[:, None]) % W.BOX, 'count': len(r)}


def loadHalo(basePath, snapNum, id, partType, fields=None):
    g = int(id)
    if partType in ('stars', 4):
        st = W.stars_at(g, snapNum)
        return _pick(_stars_dict(st, np.ones(len(st['ids']), bool)), fields)
    if partType in ('gas', 0):
        return _pick(_gas(g), fields)
    if partType in ('dm', 1):
        return _pick(_dm(g), fields)
    raise KeyError(partType)


def loadSubhalo(basePath, snapNum, id, partType, fields=None):
    assert partType in ('stars', 4), '[mock] sólo estrellas por subhalo'
    g, k = divmod(int(id), W.NSUB)
    st = W.stars_at(g, snapNum)
    sel = st['sub'] == k
    if sel.sum() == 0:
        return {'count': 0}
    return _pick(_stars_dict(st, sel), fields)


def getSnapOffsets(basePath, snapNum, id, type):
    assert type == 'Subhalo'
    sub_off, sub_len, _, _ = W.global_offsets(snapNum)
    off = np.zeros(6, np.int64); ln = np.zeros(6, np.int64)
    off[4], ln[4] = sub_off[int(id)], sub_len[int(id)]
    return {'offsetType': off, 'lenType': ln}


def loadSubset(basePath, snapNum, partType, fields=None, subset=None, mdi=None, sq=True, float32=False):
    assert partType in ('stars', 4)
    parts = [_stars_dict(W.stars_at(g, snapNum), np.ones(len(W.stars_at(g, snapNum)['ids']), bool))
             for g in range(W.NG)]
    d = {k: np.concatenate([p[k] for p in parts]) for k in STAR_MAP}
    d['count'] = len(d['ParticleIDs'])
    return _pick(d, fields)
