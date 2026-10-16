"""sublink simulado: árbol completo del central (merger inyectado) y MPBs de centrales y satélites."""
import numpy as np
from . import _world as W

# (nodo k de la rama principal que recibe al secundario, razón de masas) por FoF
MERGERS = {0: (10, 0.3), 1: (5, 0.1), 2: (30, 0.25), 3: (3, 0.3)}
L_MAIN, L_SEC = 60, 5


def _full_tree(g):
    k, ratio = MERGERS[g]
    n = L_MAIN + L_SEC
    t = {'SubhaloID': np.arange(n, dtype=np.int64),
         'SnapNum': np.zeros(n, dtype=np.int16),
         'FirstProgenitorID': np.full(n, -1, dtype=np.int64),
         'NextProgenitorID': np.full(n, -1, dtype=np.int64),
         'MainLeafProgenitorID': np.zeros(n, dtype=np.int64),
         'SubhaloMassType': np.zeros((n, 6), dtype=np.float32)}
    main = 10 * (1 - 0.005 * np.arange(L_MAIN))
    for i in range(L_MAIN):
        t['SnapNum'][i] = 99 - i
        t['FirstProgenitorID'][i] = i + 1 if i < L_MAIN - 1 else -1
        t['MainLeafProgenitorID'][i] = L_MAIN - 1
        t['SubhaloMassType'][i, 4] = main[i]
    for j in range(L_SEC):
        i = L_MAIN + j
        t['SnapNum'][i] = 99 - k - j
        t['FirstProgenitorID'][i] = i + 1 if j < L_SEC - 1 else -1
        t['MainLeafProgenitorID'][i] = n - 1
        t['SubhaloMassType'][i, 4] = ratio * main[k] * (1 - 0.05 * j)
    t['NextProgenitorID'][k] = L_MAIN
    return t


def loadTree(basePath, snapNum, id, fields=None, onlyMPB=False, onlyMDB=False, treeName="SubLink", cache=True):
    g, k = divmod(int(id), W.NSUB)
    if not onlyMPB:
        assert snapNum == 99 and k == 0, '[mock] árbol completo sólo para centrales a z=0'
        t = _full_tree(g)
    else:
        snaps = np.arange(snapNum, -1, -1)
        n = len(snaps)
        mir = np.zeros((n, 6), dtype=np.float32)
        mir[:, 4] = (5.0 if k == 0 else [None, 2.0, 0.3][k]) * (1 - 0.003 * (snapNum - snaps))
        pos = np.array([W.POS[g] if k == 0 else W.sat_pos(g, k, s) for s in snaps], dtype=np.float32)
        t = {'SnapNum': snaps, 'SubfindID': np.full(n, int(id)), 'SubhaloGrNr': np.full(n, g),
             'SubhaloPos': pos, 'SubhaloMassInRadType': mir,
             'MainLeafProgenitorID': np.full(n, 1000 + int(id), dtype=np.int64)}
    if fields is None:
        return t
    missing = [f for f in fields if f not in t]
    if missing:
        raise KeyError(f"[mock] campos de árbol no simulados: {missing}")
    out = {f: t[f] for f in fields}
    out['count'] = len(t['SnapNum'])
    return out
