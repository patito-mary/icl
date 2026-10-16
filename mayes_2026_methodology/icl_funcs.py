"""
Funciones comunes del análisis ICL/BCG (réplica de Mayes, Gómez & Monachesi 2026).

Convención de unidades (salvo que se indique otra cosa):
    distancias en kpc físicos, masas en M☉, luminosidades en L☉ (banda r),
    brillo superficial en mag/arcsec² (AB), tiempos en Gyr.

Las funciones "puras" (perfil, kernel, R_H, ajuste log, f_ICL, merger,
gradientes) no dependen de TNG y se pueden probar con datos sintéticos.
Las funciones de carga usan illustris_python y sólo corren donde están los
datos (TNG Lab).
"""
import os
import numpy as np
from scipy.special import erf
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

try:
    import illustris_python as il
except ImportError:          # permite importar el módulo sin TNG (pruebas locales)
    il = None


# ═════════════════════════════════════════════════════════════════════════
# Geometría
# ═════════════════════════════════════════════════════════════════════════
def periodic_delta(pos, center, box):
    """Vector pos − center con condiciones periódicas (mismas unidades que box)."""
    d = np.asarray(pos, dtype=float) - np.asarray(center, dtype=float)
    d = np.where(d >  box / 2, d - box, d)
    d = np.where(d < -box / 2, d + box, d)
    return d


def rotate_by_inertia_tensor(pos_rel, mass, r_lim=np.inf, reduced=True):
    """
    Alinea las partículas con los ejes principales del tensor de inercia
    (reducido: pesos 1/r², como Catalogue.comp_Reduced_Inertia_tensor).

    pos_rel : (N,3) posiciones relativas al centro [kpc]
    mass    : (N,)  masas
    r_lim   : radio máximo de las partículas usadas para el tensor [kpc]

    Devuelve pos_rot (N,3) con x̂ = eje mayor, ẑ = eje menor, y R_mat (3,3)
    tal que pos_rot = pos_rel @ R_mat.T
    """
    dist = np.linalg.norm(pos_rel, axis=1)
    ok   = (dist > 0) & (dist <= r_lim) & np.isfinite(mass)
    p, m = pos_rel[ok], mass[ok]
    if len(m) < 4 or m.sum() == 0:
        return pos_rel.copy(), np.eye(3)

    w = 1.0 / dist[ok]**2 if reduced else np.ones(len(m))
    I = np.einsum('n,ni,nj->ij', m * w, p, p) / m.sum()

    eigvals, eigvecs = np.linalg.eigh(I)
    idx   = np.argsort(eigvals)[::-1]          # mayor extensión primero
    R_mat = eigvecs[:, idx].T
    # sistema derecho (det = +1) para que la rotación sea propia
    if np.linalg.det(R_mat) < 0:
        R_mat[2] *= -1
    return pos_rel @ R_mat.T, R_mat


def projected_radius(pos_rot):
    """Radio proyectado en el plano (eje mayor, eje intermedio) [kpc]."""
    return np.hypot(pos_rot[:, 0], pos_rot[:, 1])


def wedge_mask(pos_rot, wedge_deg=15.0):
    """
    Partículas dentro de una cuña de apertura total `wedge_deg` centrada en el
    eje mayor proyectado (x̂), a ambos lados del centro.
    """
    ang = np.degrees(np.arctan2(np.abs(pos_rot[:, 1]), np.abs(pos_rot[:, 0])))
    return ang <= wedge_deg / 2


# ═════════════════════════════════════════════════════════════════════════
# Fotometría y brillo superficial
# ═════════════════════════════════════════════════════════════════════════
def luminosity_from_mag(mag_abs, m_sun):
    """Luminosidad [L☉] a partir de la magnitud absoluta y M☉ en la misma banda."""
    return 10**(-0.4 * (np.asarray(mag_abs) - m_sun))


def mu_from_sigma(sigma_L_pc2, sb_const):
    """μ [mag/arcsec²] desde Σ [L☉/pc²]: μ = M☉ + 21.572 − 2.5 log10 Σ."""
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(sigma_L_pc2 > 0, sb_const - 2.5 * np.log10(sigma_L_pc2), np.nan)


def radial_edges(r_max, binning='lineal', bin_kpc=5.0, n_bins=60, r_min=0.5):
    """Bordes de anillos [kpc]: lineales de ancho bin_kpc desde 0, o logarítmicos."""
    if binning == 'lineal':
        n = max(int(np.ceil(r_max / bin_kpc)), 1)
        return np.linspace(0.0, n * bin_kpc, n + 1)
    if binning == 'log':
        return np.logspace(np.log10(r_min), np.log10(r_max), n_bins + 1)
    raise ValueError(f"binning desconocido: {binning}")


def sb_profile_rings(R, lum, edges, sb_const):
    """
    Perfil 1D con anillos duros disjuntos.

    R : radios proyectados [kpc] ; lum : luminosidades [L☉] ; edges : bordes [kpc]
    Devuelve r_mid [kpc] (centro del anillo en área), mu [mag/arcsec²], n_part.
    """
    L_bin, _ = np.histogram(R, bins=edges, weights=lum)
    n_bin, _ = np.histogram(R, bins=edges)
    area_pc2 = np.pi * ((edges[1:] * 1e3)**2 - (edges[:-1] * 1e3)**2)
    r_mid    = np.sqrt(0.5 * (edges[1:]**2 + edges[:-1]**2))
    return r_mid, mu_from_sigma(L_bin / area_pc2, sb_const), n_bin


def gaussian_ring_area_kpc2(r0, sigma):
    """
    Área efectiva [kpc²] de un anillo con peso gaussiano en R:
        A = ∫_0^∞ 2πR exp[−(R−r0)²/(2σ²)] dR
          = 2π [σ² exp(−r0²/2σ²) + r0 σ sqrt(π/2) (1 + erf(r0/(σ√2)))]
    Para r0 ≫ σ → (2π)^{3/2} r0 σ.
    """
    r0 = np.asarray(r0, dtype=float); sigma = np.asarray(sigma, dtype=float)
    return 2 * np.pi * (sigma**2 * np.exp(-r0**2 / (2 * sigma**2))
                        + r0 * sigma * np.sqrt(np.pi / 2) * (1 + erf(r0 / (sigma * np.sqrt(2)))))


def sb_profile_kernel(R, lum, r_centers, sb_const, frac_sigma=0.15):
    """
    Perfil 1D suavizado con kernel gaussiano en R de ancho σ = frac_sigma·r0.
    Σ(r0) = Σ_i L_i w_i / A_eff(r0), con el área efectiva exacta del anillo
    gaussiano (corrige el factor √(2π) que faltaba en la versión anterior).
    Los puntos vecinos quedan correlacionados.
    """
    r_centers = np.asarray(r_centers, dtype=float)
    mu = np.full(len(r_centers), np.nan)
    for k, r0 in enumerate(r_centers):
        sig = frac_sigma * r0
        if sig <= 0:
            continue
        w = np.exp(-0.5 * ((R - r0) / sig)**2)
        sigma_L = np.sum(lum * w) / (gaussian_ring_area_kpc2(r0, sig) * 1e6)
        mu[k] = mu_from_sigma(sigma_L, sb_const)
    return r_centers, mu


def sb_profile(R, lum, r_max, sb_const, method='anillos', binning='lineal',
               bin_kpc=5.0, n_bins=60, r_min=0.5, frac_sigma=0.15):
    """
    Perfil μ(R) hasta r_max con el método elegido ('anillos' o 'kernel').
    Devuelve r_mid [kpc] y mu [mag/arcsec²].
    """
    edges = radial_edges(r_max, binning, bin_kpc, n_bins, r_min)
    r_mid, mu, _ = sb_profile_rings(R, lum, edges, sb_const)
    if method == 'anillos':
        return r_mid, mu
    if method == 'kernel':
        return sb_profile_kernel(R, lum, r_mid, sb_const, frac_sigma)
    raise ValueError(f"método de perfil desconocido: {method}")


# ═════════════════════════════════════════════════════════════════════════
# Radio de Holmberg (una sola función para 01, 01e, 01f, 02 y 03)
# ═════════════════════════════════════════════════════════════════════════
def fit_log_profile(r_mid, mu, r_fit_min=0.0, r_fit_max=np.inf, mu_range=None):
    """
    Ajuste μ(R) = a + b·log10(R) (sin pesos) a los anillos con μ finito.
    mu_range = (μ_min, μ_max) restringe el ajuste a anillos en ese rango de brillo (ajuste local).
    Devuelve a, b, err_a, err_b, n_puntos.
    """
    ok = np.isfinite(mu) & (r_mid > 0) & (r_mid >= r_fit_min) & (r_mid <= r_fit_max)
    if mu_range is not None:
        ok &= (mu >= mu_range[0]) & (mu <= mu_range[1])
    if ok.sum() < 3:
        return np.nan, np.nan, np.nan, np.nan, int(ok.sum())
    x, y = np.log10(r_mid[ok]), mu[ok]
    (b, a), cov = np.polyfit(x, y, 1, cov=True)
    return a, b, np.sqrt(cov[1, 1]), np.sqrt(cov[0, 0]), int(ok.sum())


def holmberg_from_fit(a, b, mu_cut):
    """R_H [kpc] donde a + b·log10(R) = mu_cut (NaN si b ≤ 0)."""
    if not (np.isfinite(a) and np.isfinite(b)) or b <= 0:
        return np.nan
    return 10**((mu_cut - a) / b)


def holmberg_from_interp(r_mid, mu, mu_cut):
    """R_H [kpc] por interpolación lineal en el PRIMER cruce de mu_cut hacia afuera."""
    ok = np.isfinite(mu) & (r_mid > 0)
    r_v, m_v = r_mid[ok], mu[ok]
    if len(r_v) < 2 or m_v[0] > mu_cut:
        return np.nan
    above = np.where(m_v >= mu_cut)[0]
    if len(above) == 0:
        return np.nan
    k = above[0]
    if k == 0:
        return r_v[0]
    return float(interp1d(m_v[k-1:k+1], r_v[k-1:k+1])(mu_cut))


def holmberg_radius(r_mid, mu, mu_cut=26.5, method='ajuste', r_fit_min=0.0, r_fit_max=np.inf,
                    mu_window=1.0):
    """
    Radio de Holmberg [kpc] y diagnósticos. method:
      'ajuste'       R_H donde el ajuste logarítmico GLOBAL (todos los anillos hasta r_fit_max)
                     cruza mu_cut. Es lo que muestra la Fig. 1 del paper (la línea roja a 200 kpc
                     coincide con el cruce del ajuste, no con el de los datos). Defecto.
      'interp'       primer cruce de mu_cut sobre el perfil medido.
      'ajuste_local' ajuste logarítmico sólo con anillos en mu_cut ± mu_window.
    Devuelve dict con R_H (el del método elegido), R_H_fit, R_H_interp, R_H_local,
    a, b, err_a, err_b, n_fit (los coeficientes son siempre los del ajuste global).
    Los R_H de ajustes se descartan (NaN) si caen fuera de (0, r_fit_max].
    """
    def _fit_rh(mu_range):
        a_, b_, ea_, eb_, n_ = fit_log_profile(r_mid, mu, r_fit_min, r_fit_max, mu_range)
        r_ = holmberg_from_fit(a_, b_, mu_cut)
        if np.isfinite(r_) and not (0 < r_ <= r_fit_max):
            r_ = np.nan
        return r_, (a_, b_, ea_, eb_, n_)
    r_fit, (a, b, ea, eb, n) = _fit_rh(None)
    r_loc, _ = _fit_rh((mu_cut - mu_window, mu_cut + mu_window))
    r_int = holmberg_from_interp(r_mid, mu, mu_cut)
    choice = {'ajuste': r_fit, 'interp': r_int, 'ajuste_local': r_loc}
    if method not in choice:
        raise ValueError(f"método de R_H desconocido: {method}")
    return {'R_H': choice[method], 'R_H_fit': r_fit, 'R_H_interp': r_int, 'R_H_local': r_loc,
            'a': a, 'b': b, 'err_a': ea, 'err_b': eb, 'n_fit': n}


# ═════════════════════════════════════════════════════════════════════════
# Separación BCG / ICL y fracciones de masa
# ═════════════════════════════════════════════════════════════════════════
def split_bcg_icl(R, R_H):
    """Máscaras BCG (R ≤ R_H) e ICL (R > R_H) con el corte circular proyectado."""
    if not np.isfinite(R_H):
        z = np.zeros(len(R), dtype=bool)
        return z, z.copy()
    return R <= R_H, R > R_H


def icl_fractions(mass, mask_bcg, mask_icl, mstar200):
    """
    M_BCG, M_ICL [M☉] y fracciones:
        f_icl        = M_ICL / M*200            (definición del paper, §3)
        f_icl_bcgicl = M_ICL / (M_BCG + M_ICL)  (secundaria)
    """
    m_bcg = float(mass[mask_bcg].sum())
    m_icl = float(mass[mask_icl].sum())
    tot   = m_bcg + m_icl
    return {'M_bcg': m_bcg, 'M_icl': m_icl,
            'f_icl': m_icl / mstar200 if mstar200 > 0 else np.nan,
            'f_icl_bcgicl': m_icl / tot if tot > 0 else np.nan}


# ═════════════════════════════════════════════════════════════════════════
# Propiedades del halo: estado dinámico, concentración, último major merger
# ═════════════════════════════════════════════════════════════════════════
def center_of_mass(pos_rel, mass):
    """Centro de masa de posiciones ya relativas (periódicas resueltas) [kpc]."""
    return np.sum(pos_rel * mass[:, None], axis=0) / mass.sum()


def dynamical_state(offset_frac, thr_relaxed=0.1, thr_disturbed=0.4):
    """0 = relajado (< thr_relaxed), 2 = perturbado (> thr_disturbed), 1 = intermedio."""
    off = np.asarray(offset_frac, dtype=float)
    state = np.ones(off.shape, dtype=int)
    state[off < thr_relaxed]   = 0
    state[off > thr_disturbed] = 2
    state[~np.isfinite(off)]   = -1
    return state


def nfw_log_density(r, log_rho_s, r_s):
    """log10 ρ_NFW(r) = log10 ρ_s − log10[(r/r_s)(1 + r/r_s)²]."""
    x = r / r_s
    return log_rho_s - np.log10(x * (1 + x)**2)


def fit_nfw_concentration(r_dm, m_part, r200, r_min_frac=0.05, r_max_frac=1.0, n_bins=20):
    """
    Concentración c = R200/r_s ajustando NFW al perfil de densidad de DM en
    cáscaras logarítmicas entre r_min_frac y r_max_frac de R200.
    r_dm [kpc], m_part [M☉] (escalar o array), r200 [kpc].
    Devuelve c, r_s [kpc], r_mid, rho [M☉/kpc³].
    """
    m_part = np.broadcast_to(np.asarray(m_part, dtype=float), np.shape(r_dm))
    edges  = np.logspace(np.log10(r_min_frac * r200), np.log10(r_max_frac * r200), n_bins + 1)
    m_bin, _ = np.histogram(r_dm, bins=edges, weights=m_part)
    vol   = 4/3 * np.pi * (edges[1:]**3 - edges[:-1]**3)
    r_mid = np.sqrt(edges[1:] * edges[:-1])
    rho   = m_bin / vol
    ok = rho > 0
    if ok.sum() < 4:
        return np.nan, np.nan, r_mid, rho
    try:
        p0 = [np.log10(rho[ok][len(rho[ok]) // 2]) + 1, r200 / 5]
        popt, _ = curve_fit(nfw_log_density, r_mid[ok], np.log10(rho[ok]), p0=p0,
                            bounds=([-np.inf, 1e-3 * r200], [np.inf, 2 * r200]))
        return r200 / popt[1], popt[1], r_mid, rho
    except (RuntimeError, ValueError):
        return np.nan, np.nan, r_mid, rho


def max_past_mass(tree, index, mass_field='SubhaloMassType', part_idx=4):
    """Masa máxima a lo largo de la rama principal que empieza en `index` (unidades del árbol)."""
    n = tree['MainLeafProgenitorID'][index] - tree['SubhaloID'][index] + 1
    m = tree[mass_field][index:index + n]
    m = m[:, part_idx] if m.ndim == 2 else m
    return float(np.max(m)) if len(m) else 0.0


def last_major_merger_snap(tree, ratio_min=0.2, mass_field='SubhaloMassType', part_idx=4):
    """
    Último (más reciente) merger con razón de masa estelar > ratio_min en la
    rama principal de un árbol SubLink completo (raíz en el índice 0).
    La razón usa la masa máxima pasada de cada progenitor (Rodriguez-Gomez+2015).
    Devuelve (SnapNum del descendiente donde ya están fusionados, razón) o (-1, nan).
    """
    root = tree['SubhaloID'][0]
    idx  = 0
    fp   = tree['FirstProgenitorID'][0]
    while fp != -1:
        fpi = int(fp - root)
        m1  = max_past_mass(tree, fpi, mass_field, part_idx)
        nxt = tree['NextProgenitorID'][fpi]
        while nxt != -1:
            npi = int(nxt - root)
            m2  = max_past_mass(tree, npi, mass_field, part_idx)
            if m1 > 0 and m2 > 0:
                ratio = min(m1, m2) / max(m1, m2)
                if ratio > ratio_min:
                    return int(tree['SnapNum'][idx]), ratio
            nxt = tree['NextProgenitorID'][npi]
        idx = fpi
        fp  = tree['FirstProgenitorID'][fpi]
    return -1, np.nan


# ═════════════════════════════════════════════════════════════════════════
# Propiedades estelares
# ═════════════════════════════════════════════════════════════════════════
def half_mass_radius(r, m):
    """Radio que encierra la mitad de la masa [mismas unidades que r]."""
    if len(r) == 0 or np.sum(m) <= 0:
        return np.nan
    idx = np.argsort(r)
    cum = np.cumsum(m[idx])
    return float(np.interp(0.5 * cum[-1], cum, r[idx]))


def mean_log_metallicity(Z, m, z_sun):
    """log10(<Z>_masa / Z☉)."""
    if len(Z) == 0 or np.sum(m) <= 0:
        return np.nan
    return float(np.log10(np.average(Z, weights=m) / z_sun))


def integrated_color(mag1, mag2):
    """Color integrado m1 − m2 = −2.5 log10(Σ 10^{−0.4 m1} / Σ 10^{−0.4 m2})."""
    if len(mag1) == 0:
        return np.nan
    return float(-2.5 * np.log10(np.sum(10**(-0.4 * mag1)) / np.sum(10**(-0.4 * mag2))))


def lookback_from_scale_factor(a_form, cosmo, n_grid=2000):
    """Edad (lookback de formación) [Gyr] desde GFM_StellarFormationTime (a), interpolando."""
    a_grid = np.linspace(1e-3, 1.0, n_grid)
    t_grid = cosmo.lookback_time(1 / a_grid - 1).value
    return np.interp(np.asarray(a_form), a_grid, t_grid)


def equal_number_edges(x, n_bins):
    """Bordes con (aprox.) el mismo número de partículas por bin."""
    return np.quantile(x, np.linspace(0, 1, n_bins + 1))


def radial_gradient(R, n_bins, stat_func, r_floor=1e-3):
    """
    Gradiente d(stat)/d log10(R) con bins de igual número de partículas.
    stat_func(idx) → valor de la propiedad en el bin (idx = índices en R).
    Devuelve pendiente, error, intercepto, log10(R_med_bin), valores por bin.
    """
    if len(R) < 3 * n_bins:
        return np.nan, np.nan, np.nan, np.array([]), np.array([])
    order = np.argsort(R)
    chunks = np.array_split(order, n_bins)
    x = np.array([np.log10(max(np.median(R[c]), r_floor)) for c in chunks])
    y = np.array([stat_func(c) for c in chunks])
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return np.nan, np.nan, np.nan, x, y
    res = linregress(x[ok], y[ok])
    return res.slope, res.stderr, res.intercept, x, y


def linfit_beta(x, y):
    """Pendiente β (± error estándar) e intercepto de un ajuste lineal (ignora NaN)."""
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return np.nan, np.nan, np.nan
    res = linregress(x[ok], y[ok])
    return res.slope, res.stderr, res.intercept


# ═════════════════════════════════════════════════════════════════════════
# Componentes de ensamblaje (catálogo por partícula, Rodriguez-Gomez+2016)
# ═════════════════════════════════════════════════════════════════════════
COMP_LABELS = {0: 'In situ', 1: 'Mergers completados', 2: 'Sobrevivientes'}


def classify_components(insitu, origin):
    """
    Componente por partícula:
        0 = in situ (InSitu == 1)
        1 = mergers completados (ex situ, AccretionOrigin == 0)
        2 = sobrevivientes (ex situ, AccretionOrigin ∈ {1, 2}: ongoing + flybys;
            a z=0 'ongoing' es vacío por definición)
       −1 = sin clasificar (p. ej. ex situ formado fuera de galaxias, AccretionOrigin = −1)
    """
    insitu = np.asarray(insitu); origin = np.asarray(origin)
    comp = np.full(len(insitu), -1, dtype=int)
    comp[insitu == 1] = 0
    ex = insitu == 0
    comp[ex & (origin == 0)] = 1
    comp[ex & ((origin == 1) | (origin == 2))] = 2
    return comp


def component_fractions(mass, comp, mask=None):
    """Fracción de masa en cada componente (0,1,2) sobre la masa clasificada (comp ≥ 0)."""
    if mask is None:
        mask = np.ones(len(mass), dtype=bool)
    sel = mask & (comp >= 0)
    tot = mass[sel].sum()
    out = {k: (mass[sel & (comp == k)].sum() / tot if tot > 0 else np.nan) for k in (0, 1, 2)}
    out['f_sin_clasificar'] = (mass[mask & (comp < 0)].sum() / mass[mask].sum()
                               if mass[mask].sum() > 0 else np.nan)
    return out


def load_sa_stars(path_pattern, snap, offset, length, field_insitu, field_origin):
    """
    Lee el catálogo de ensamblaje POR PARTÍCULA para un tramo contiguo
    [offset, offset+length) de PartType4 del snapshot (orden del snapshot).
    Falla con un mensaje claro si el archivo o los campos no existen.
    """
    import h5py
    path = path_pattern.format(snap=snap)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"No se encontró el catálogo de ensamblaje por partícula: {path}\n"
            "  Ajusta P.SA_STARS_FILE (o la variable de entorno TNG_SA_STARS_FILE).\n"
            "  En el TNG Lab búscalo con: ls /home/tnguser/sims.TNG/TNG100-1/postprocessing/")
    with h5py.File(path, 'r') as f:
        keys = list(f.keys())
        missing = [k for k in (field_insitu, field_origin) if k not in f]
        if missing:
            raise KeyError(f"Campos {missing} no están en {path}. Campos disponibles: {keys}. "
                           "Ajusta P.SA_FIELD_INSITU / P.SA_FIELD_ORIGIN.")
        insitu = f[field_insitu][offset:offset + length]
        origin = f[field_origin][offset:offset + length]
    if len(insitu) != length:
        raise ValueError(f"El catálogo SA devolvió {len(insitu)} partículas y se esperaban {length}: "
                         "¿está alineado con el orden de PartType4 del snapshot?")
    return insitu, origin


# ═════════════════════════════════════════════════════════════════════════
# Progenitores (§2.3, §7)
# ═════════════════════════════════════════════════════════════════════════
def subhalo_of_local_index(local_idx, lens_star, first_sub):
    """
    Subfind ID de partículas de un FoF dado su índice local en loadHalo
    (orden [sub0][sub1]...[subN−1][fuzz]). −1 = fuzz (no ligado a subhalo).
    """
    cum = np.concatenate([[0], np.cumsum(lens_star)])
    k = np.searchsorted(cum, local_idx, side='right') - 1
    out = first_sub + k
    out[(local_idx >= cum[-1]) | (k < 0)] = -1
    return out


def detach_update(tracked_ids, done, ids_in_central, snap, out_snap):
    """
    Marca, en el snapshot `snap`, las partículas rastreadas que YA NO están en
    el progenitor principal del BCG (primera vez hacia atrás). Modifica done y
    out_snap in situ y devuelve la máscara de las recién marcadas.
    """
    new = (~done) & (~np.isin(tracked_ids, ids_in_central))
    done[new] = True
    out_snap[new] = snap
    return new


def shared_fraction(prog_a, prog_b):
    """Fracción de partículas de A cuyo progenitor también aporta a B (progenitores válidos ≥ 0)."""
    prog_a = np.asarray(prog_a); prog_b = np.asarray(prog_b)
    va = prog_a[prog_a >= 0]
    if len(va) == 0:
        return np.nan
    return float(np.isin(va, np.unique(prog_b[prog_b >= 0])).mean())


def significant_progenitors(prog_ids, mass, frac=0.90):
    """
    Contribución del progenitor más significativo y número de progenitores que
    suman `frac` de la masa acretada. Devuelve (id_top, f_top, n_frac).
    """
    prog_ids = np.asarray(prog_ids); mass = np.asarray(mass)
    ok = prog_ids >= 0
    if ok.sum() == 0:
        return -1, np.nan, 0
    u, inv = np.unique(prog_ids[ok], return_inverse=True)
    m_u = np.bincount(inv, weights=mass[ok])
    order = np.argsort(m_u)[::-1]
    cum = np.cumsum(m_u[order]) / m_u.sum()
    return int(u[order[0]]), float(cum[0]), int(np.searchsorted(cum, frac) + 1)


# ═════════════════════════════════════════════════════════════════════════
# Carga de partículas de TNG (sólo en el Lab)
# ═════════════════════════════════════════════════════════════════════════
STAR_FIELDS = ['Coordinates', 'Masses', 'GFM_StellarPhotometrics', 'GFM_StellarFormationTime']


def _require_il():
    if il is None:
        raise ImportError("illustris_python no está disponible: esta función sólo corre en el TNG Lab.")


def load_central_stars(basePath, snap, sub_id, center_kpc, r200_kpc, box_kpc,
                       UL, UM, extra_fields=(), r_outer_frac=1.0):
    """
    Estrellas del subhalo central (BCG+ICL según Subfind), sin partículas de
    viento (GFM_StellarFormationTime ≤ 0), centradas en el BCG y con r_3D ≤ R200.

    Devuelve dict con pos [kpc], r3d, mass [M☉], phot (N,8), sft y extras,
    más 'idx_sub' (índice de cada partícula dentro del subhalo, para alinear
    con catálogos por partícula) y 'offset'/'length' del subhalo en PartType4.
    """
    _require_il()
    fields = list(dict.fromkeys(STAR_FIELDS + list(extra_fields)))
    st = il.snapshot.loadSubhalo(basePath, snap, int(sub_id), 'stars', fields=fields)
    if not isinstance(st, dict) or st.get('count', 0) == 0:
        return None
    return _prepare_stars(st, fields, center_kpc, r200_kpc, box_kpc, UL, UM, r_outer_frac,
                          offsets=il.snapshot.getSnapOffsets(basePath, snap, int(sub_id), 'Subhalo'))


def mask_from_lens(lens_star, n_total_fof, include_satellites=False, include_fuzz=True):
    """
    Máscara sobre las estrellas de un FoF ordenadas [sub0][sub1]...[subN−1][fuzz]:
    siempre incluye el subhalo central (sub0); opcionalmente satélites y fuzz
    (partículas del FoF no ligadas a ningún subhalo).
    """
    lens_star = np.asarray(lens_star, dtype=int)
    n_central = int(lens_star[0]) if len(lens_star) else 0
    n_in_subs = int(lens_star.sum())
    mask = np.zeros(n_total_fof, dtype=bool)
    mask[:n_central] = True
    if include_satellites:
        mask[n_central:n_in_subs] = True
    if include_fuzz:
        mask[n_in_subs:] = True
    return mask


def load_fof_stars(basePath, snap, fof_idx, lens_star, center_kpc, r200_kpc, box_kpc, UL, UM,
                   include_satellites=False, include_fuzz=True, r_outer_frac=1.0):
    """
    Estrellas del FoF con la máscara de mask_from_lens, sin viento y con r_3D ≤ R200
    (chequeo 01e: central + fuzz). Mismo formato de salida que load_central_stars.
    """
    _require_il()
    st = il.snapshot.loadHalo(basePath, snap, int(fof_idx), 'stars', fields=STAR_FIELDS)
    if not isinstance(st, dict) or st.get('count', 0) == 0:
        return None
    sel = mask_from_lens(lens_star, st['count'], include_satellites, include_fuzz)
    st = {k: (v[sel] if k != 'count' else int(sel.sum())) for k, v in st.items()}
    return _prepare_stars(st, STAR_FIELDS, center_kpc, r200_kpc, box_kpc, UL, UM, r_outer_frac)


def _prepare_stars(st, fields, center_kpc, r200_kpc, box_kpc, UL, UM, r_outer_frac, offsets=None):
    """Conversión de unidades, máscara de viento y corte 3D en R200 (uso interno)."""
    phot = st['GFM_StellarPhotometrics']
    if phot.ndim == 1:
        phot = phot.reshape(-1, 8)
    pos  = periodic_delta(st['Coordinates'] * UL, center_kpc, box_kpc)
    r3d  = np.linalg.norm(pos, axis=1)
    keep = (st['GFM_StellarFormationTime'] > 0) & (r3d <= r_outer_frac * r200_kpc)
    out = {'pos': pos[keep], 'r3d': r3d[keep], 'mass': st['Masses'][keep] * UM,
           'phot': phot[keep], 'sft': st['GFM_StellarFormationTime'][keep],
           'idx_sub': np.nonzero(keep)[0], 'n_total': len(keep)}
    for k in fields:
        if k not in ('Coordinates', 'Masses', 'GFM_StellarPhotometrics', 'GFM_StellarFormationTime'):
            out[k] = st[k][keep]
    if offsets is not None:
        out['offset'] = int(offsets['offsetType'][4])
        out['length'] = int(offsets['lenType'][4])
    return out


def bcg_icl_from_stars(stars, mstar200, r200_kpc, P, R_mat=None):
    """
    Pipeline §2.2 sobre estrellas ya preparadas: rotación por tensor de inercia,
    perfil μ_r(R) hasta R200, R_H (función única) y masas BCG/ICL.
    Si se da R_mat se reutiliza esa rotación (para 02/03).
    Devuelve dict con R (radio proyectado), pos_rot, R_mat, perfil, R_H y fracciones.
    """
    if R_mat is None:
        pos_rot, R_mat = rotate_by_inertia_tensor(stars['pos'], stars['mass'], reduced=P.INERTIA_REDUCED)
    else:
        pos_rot = stars['pos'] @ R_mat.T
    R   = projected_radius(pos_rot)
    lum = luminosity_from_mag(stars['phot'][:, P.BAND_IDX], P.M_SUN_R_AB)
    r_max = P.R_OUTER_FRAC * r200_kpc
    r_mid, mu = sb_profile(R, lum, r_max, P.SB_CONST, method=P.PROFILE_METHOD,
                           binning=P.PROFILE_BINNING, bin_kpc=P.SB_BIN_KPC,
                           n_bins=P.SB_N_BINS_LOG, r_min=P.SB_R_MIN_KPC,
                           frac_sigma=P.KERNEL_FRAC_SIGMA)
    hol = holmberg_radius(r_mid, mu, P.MU_HOLMBERG, method=P.RH_METHOD,
                          r_fit_min=P.R_FIT_MIN_KPC, r_fit_max=r_max,
                          mu_window=getattr(P, 'RH_LOCAL_MU_WINDOW', 1.0))
    m_bcg, m_icl = split_bcg_icl(R, hol['R_H'])
    out = {'R': R, 'pos_rot': pos_rot, 'R_mat': R_mat, 'lum': lum,
           'r_mid': r_mid, 'mu': mu, 'mask_bcg': m_bcg, 'mask_icl': m_icl}
    out.update(hol)
    out.update(icl_fractions(stars['mass'], m_bcg, m_icl, mstar200))
    return out
