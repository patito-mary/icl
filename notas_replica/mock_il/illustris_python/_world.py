"""
Universo sintético mínimo con la interfaz de illustris_python, para probar los notebooks de la réplica
sin datos de TNG. Unidades internas de TNG: ckpc/h, 1e10 Msun/h.

4 FoF (3 pasan el corte Group_M_Crit200 >= 1000), 3 subhalos por FoF (central + 2 satélites) + fuzz.
El FoF 0 cruza el borde de la caja (prueba de periodicidad).
Verdades inyectadas (para comparar con lo que miden los notebooks):
  - estado dinámico: offset del CM del gas = [0.05, 0.25, 0.55] R200   → [relajado, intermedio, perturbado]
  - concentración NFW de la DM = [5, 7, 9]
  - major merger (razón 0.3) del FoF 0 registrado en el snapshot 90; FoF 1 sólo minor (0.1); FoF 2 razón 0.25 en snap 70
  - cada estrella ex situ del central tiene un snapshot de despegue d ∈ [60, 98] y un satélite de origen
  - satélite 1 del FoF 1 tiene SubhaloFlag = 0 en snapshots 74–90 (dentro de R200) → reemplazo por la muestra 2 (snap 73)
"""
import os
import numpy as np
import h5py
from functools import lru_cache
from scipy.special import gammaincinv

H = 0.6774
BOX = 75000.0
NG = 4
M_INT = np.array([2500., 1400., 1050., 950.])
R_INT = np.array([900., 750., 680., 660.])          # ckpc/h
POS = np.array([[74950., 100., 37000.], [20000., 30000., 40000.],
                [50000., 10000., 60000.], [10000., 60000., 20000.]])
NSUB = 3
N_CEN = [24000, 16000, 12000, 10000]
N_SAT = 1500
N_FUZZ = 800
M_STAR = 1e-4                                         # 1e10 Msun/h por partícula
M_DM_INT = 0.000505574296436975                       # masa DM [1e10 Msun/h] (valor documentado TNG100-1)
DYN_OFF = [0.05, 0.25, 0.55, 0.1]
C_NFW = [5., 7., 9., 6.]
MDIR = os.environ.get('MOCK_TNG_DIR', '/tmp/mock_tng')


def z_of(s):
    return (99 - s) * 0.03


def sat_dist0(g, j):
    """Distancia de los satélites a z=0 (fracción de R200)."""
    return [0.5, 0.8][j - 1]


def sat_pos(g, j, s):
    """Posición comóvil del satélite j (1 o 2) del FoF g en el snapshot s: se aleja hacia atrás en el tiempo."""
    ang = 0.7 * j + g
    r = sat_dist0(g, j) * R_INT[g] * (1 + (99 - s) * 0.04)
    v = np.array([np.cos(ang), np.sin(ang), 0.3])
    return (POS[g] + r * v / np.linalg.norm(v)) % BOX


def rand_rot(rng):
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return q


@lru_cache(None)
def base_stars(g):
    """Estrellas del FoF g a z=0 en orden [sub0][sub1][sub2][fuzz] y metadatos por partícula."""
    rng = np.random.default_rng(100 + g)
    R200 = R_INT[g]
    # central: BCG (Sérsic n=4 en 3D aprox.) + halo ICL ρ ∝ r^-3 hasta 1.3 R200
    nc = N_CEN[g]; nb = int(0.6 * nc); ni = nc - nb
    bn = gammaincinv(8, 0.5)
    r_b = 25 * H * (gammaincinv(8, rng.random(nb)) / bn)**4
    r_b = np.minimum(r_b, 0.9 * R200)
    lo, hi = 15 * H, 1.3 * R200
    r_i = lo * (hi / lo)**rng.random(ni)                  # p(r) ∝ 1/r  ⇔ ρ ∝ r^-3
    r_c = np.concatenate([r_b, r_i])
    u = rng.normal(size=(nc, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
    pc = u * r_c[:, None] * np.array([1.0, 0.7, 0.5])
    pc = pc @ rand_rot(rng).T + POS[g]
    pos = [pc]; sub_local = [np.zeros(nc, int)]
    for j in (1, 2):
        ps = rng.normal(scale=4 * H, size=(N_SAT, 3)) + sat_pos(g, j, 99)
        pos.append(ps); sub_local.append(np.full(N_SAT, j))
    pf = POS[g] + rng.uniform(-1.2, 1.2, size=(N_FUZZ, 3)) * R200 / np.sqrt(3)
    pos.append(pf); sub_local.append(np.full(N_FUZZ, -1))
    pos = np.concatenate(pos) % BOX
    sub_local = np.concatenate(sub_local)
    n = len(pos)
    r_phys = np.linalg.norm(((pos - POS[g] + BOX / 2) % BOX - BOX / 2), axis=1) / H
    mass = np.full(n, M_STAR) * rng.uniform(0.8, 1.2, n)
    L_r = mass * 1e10 / H / 3.0                           # M/L_r = 3
    Mr = 4.65 - 2.5 * np.log10(L_r) + rng.normal(0, 0.1, n)
    phot = np.stack([Mr + 1.6, Mr + 1.1, Mr + 0.35, Mr - 2.4, Mr + 0.5, Mr, Mr - 0.2, Mr - 0.3], axis=1)
    phot[r_phys > 100, 1] -= 0.15                          # ICL más azul en B
    sft = rng.uniform(0.05, 1.0, n); sft[rng.random(n) < 0.01] = -0.3   # 1 % viento
    Z = np.clip(0.02 * (np.maximum(r_phys, 1) / 10)**-0.3, 1e-4, 0.06)
    ids = (g + 1) * 10_000_000 + rng.permutation(n)
    # ensamblaje: in situ concentrado; ex situ = mergers completos (0) o flybys/sobrevivientes (2)
    insitu = np.where((sub_local == 0) & (r_phys < 30) & (rng.random(n) < 0.6), 1, 0)
    insitu[sub_local > 0] = 1
    origin = np.where(insitu == 1, -1, np.where(r_phys < 150, 0, 2))
    origin[(insitu == 0) & (rng.random(n) < 0.02)] = -1   # ex situ formado fuera de galaxias
    insitu[(insitu == 0) & (origin == -1)] = 0
    # despegue y satélite de origen para estrellas ex situ del central
    det = np.full(n, -1)
    ex_c = (sub_local == 0) & (insitu == 0)
    det[ex_c] = rng.integers(60, 99, ex_c.sum())
    orig_sat = np.zeros(n, int)
    orig_sat[ex_c] = rng.integers(1, 3, ex_c.sum())
    return dict(pos=pos, sub=sub_local, mass=mass, phot=phot, sft=sft, Z=Z, ids=ids,
                insitu=insitu, origin=origin, det=det, orig_sat=orig_sat)


def first_sub(g):
    return NSUB * g


@lru_cache(None)
def stars_at(g, s):
    """Estrellas del FoF g en el snapshot s, ordenadas [sub0][sub1][sub2][fuzz]; posiciones comóviles."""
    b = base_stars(g)
    sub = b['sub'].copy(); pos = b['pos'].copy()
    if s < 99:
        moved = (b['det'] >= 0) & (s <= b['det'])        # aún en su satélite de origen
        sub[moved] = b['orig_sat'][moved]
        for j in (1, 2):
            m = moved & (b['orig_sat'] == j)
            pos[m] = (sat_pos(g, j, s) + (pos[m] - POS[g]) * 0.01) % BOX
            msat = (b['sub'] == j)
            pos[msat] = (sat_pos(g, j, s) + (b['pos'][msat] - sat_pos(g, j, 99))) % BOX
    key = np.where(sub < 0, NSUB, sub)
    order = np.argsort(key, kind='stable')
    lens = np.array([(sub == k).sum() for k in range(NSUB)])
    out = {k: (v[order] if isinstance(v, np.ndarray) and len(v) == len(sub) else v) for k, v in b.items()}
    out['pos'] = pos[order]; out['sub'] = sub[order]; out['lens'] = lens
    return out


def global_offsets(s):
    """Offsets de cada subhalo y FoF en el PartType4 global (orden FoF por FoF)."""
    sub_off = np.zeros(NG * NSUB, np.int64); sub_len = np.zeros(NG * NSUB, np.int64); grp_off = []
    o = 0
    for g in range(NG):
        st = stars_at(g, s); grp_off.append(o)
        for k in range(NSUB):
            sub_off[first_sub(g) + k] = o; sub_len[first_sub(g) + k] = st['lens'][k]; o += st['lens'][k]
        o += len(st['ids']) - st['lens'].sum()
    return sub_off, sub_len, np.array(grp_off), o


def sub_flag(s):
    f = np.ones(NG * NSUB, dtype=np.uint8)
    if 74 <= s <= 90:
        f[first_sub(1) + 1] = 0      # satélite 1 del FoF 1 anómalo mientras está dentro de R200
    return f


def ensure_files():
    """Archivos HDF5 que los notebooks abren directamente: header del snapshot, offsets y SA por partícula."""
    os.makedirs(os.path.join(MDIR, 'output'), exist_ok=True)     # basePath simulado
    sa = os.path.join(MDIR, 'stars_099.hdf5')
    if os.path.exists(sa):
        return
    _, _, _, ntot = global_offsets(99)
    with h5py.File(os.path.join(MDIR, 'snap_099.0.hdf5'), 'w') as f:
        hd = f.create_group('Header')
        hd.attrs['NumPart_Total'] = np.array([0, 0, 0, 0, ntot, 0], dtype=np.uint32)
        hd.attrs['NumPart_Total_HighWord'] = np.zeros(6, dtype=np.uint32)
        # claves de la doc de TNG para el Header del snapshot; MOCK_SIN_MASSTABLE=1 simula que falta
        if os.environ.get('MOCK_SIN_MASSTABLE') != '1':
            hd.attrs['MassTable'] = np.array([0, M_DM_INT, 0, 0, 0, 0])
        hd.attrs['NumPart_ThisFile'] = hd.attrs['NumPart_Total']
        hd.attrs['BoxSize'] = BOX; hd.attrs['Redshift'] = 0.0; hd.attrs['Time'] = 1.0
        hd.attrs['Omega0'] = 0.3089; hd.attrs['OmegaLambda'] = 0.6911; hd.attrs['NumFilesPerSnapshot'] = 1
    with h5py.File(sa, 'w') as f:
        ins = np.concatenate([base_stars(g)['insitu'][np.argsort(np.where(base_stars(g)['sub'] < 0, NSUB, base_stars(g)['sub']), kind='stable')] for g in range(NG)])
        ori = np.concatenate([base_stars(g)['origin'][np.argsort(np.where(base_stars(g)['sub'] < 0, NSUB, base_stars(g)['sub']), kind='stable')] for g in range(NG)])
        f.create_dataset('InSitu', data=ins.astype(np.int8))
        f.create_dataset('AccretionOrigin', data=ori.astype(np.int8))
