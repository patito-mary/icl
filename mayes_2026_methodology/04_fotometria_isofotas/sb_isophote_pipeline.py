"""
Fotometría de superficie sintética "como si fuese observacional".

A diferencia de 01_-01f_ (perfil 1D directo sobre el radio 2D de cada
partícula), este pipeline:
  1. Deposita la luminosidad de las partículas en una imagen 2D (grilla +
     suavizado gaussiano), imitando una imagen observacional real.
  2. Ajusta isofotas elípticas sobre esa imagen (photutils.isophote,
     Jedrzejewski 1987).
  3. Ajusta un perfil de Sérsic a mu_B(a).
  4. Estima el radio de Holmberg (mu_B = 26.5) de dos formas independientes:
     interpolación empírica del perfil medido, y a partir del modelo de Sérsic.

Nota sobre M_sun,B: aquí se usa 5.44 (valor pedido para este pipeline),
distinto del M_SUN_B_VEGA = 5.36 (Willmer 2018) usado en params_icl.py para
01_-01f_. Se mantiene independiente para no alterar esos notebooks.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse as MplEllipse
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from photutils.isophote import EllipseGeometry, Ellipse

# ── Constantes fotométricas ────────────────────────────────────────────────
M_SUN_B     = 5.44                    # magnitud absoluta del Sol en B (Vega)
MU_HOLMBERG = 26.5                    # mag/arcsec^2 (Holmberg 1958)


# ─────────────────────────────────────────────────────────────────────────
# 1. Carga, proyección e imagen sintética
# ─────────────────────────────────────────────────────────────────────────

def load_subhalo_stars(il, basePath, snap, sub_id, cen_pos, box_size_kpc,
                        UL, UM, Distance_3D):
    """
    Carga posiciones/masa/fotometría de las partículas estelares de un
    subhalo, convierte a unidades físicas y centra en `cen_pos` (kpc).

    `il` (illustris_python), `Distance_3D` (Catalogue.Distance_3D) se pasan
    como dependencias para no fijar rutas de import dentro del módulo.

    Retorna: pos_c (N,3) kpc, mass (N,) Msun, phot (N,8) mag abs.
    """
    fields = ['Coordinates', 'Masses', 'GFM_StellarPhotometrics']
    stars  = il.snapshot.loadSubhalo(basePath, snap, sub_id, 'stars', fields=fields)

    pos  = stars['Coordinates'] * UL
    mass = stars['Masses'] * UM
    phot = stars['GFM_StellarPhotometrics']
    if phot.ndim == 1:
        phot = phot.reshape(-1, 8)

    pos_c = Distance_3D(pos, cen_pos, box_size_kpc)
    return pos_c, mass, phot


def rotate_by_inertia_tensor(pos_rel, mass, r_lim=np.inf):
    """Alinea las partículas con los ejes principales del tensor de inercia reducido (idéntica a 01b_/01e_)."""
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
    I = np.array([[Ixx, Ixy, Ixz],
                  [Ixy, Iyy, Iyz],
                  [Ixz, Iyz, Izz]])

    eigvals, eigvecs = np.linalg.eigh(I)
    idx   = np.argsort(eigvals)[::-1]
    R_mat = eigvecs[:, idx].T
    return pos_rel @ R_mat.T, R_mat


_PROJ_AXES = {'xy': (0, 1), 'xz': (0, 2), 'yz': (1, 2)}


def project_particles(pos_c, mass, projection='xy', align_principal=False, r_lim=np.inf):
    """
    Proyecta las partículas centradas (pos_c, kpc) sobre un plano.

    projection      : 'xy' | 'xz' | 'yz' — plano de proyección.
    align_principal : si True, rota primero al sistema de ejes principales
                      del tensor de inercia (face-on respecto al eje menor)
                      antes de seleccionar el plano.

    Retorna: x, y (kpc).
    """
    if projection not in _PROJ_AXES:
        raise ValueError(f"projection debe ser una de {list(_PROJ_AXES)}")
    pos = pos_c
    if align_principal:
        pos, _ = rotate_by_inertia_tensor(pos_c, mass, r_lim=r_lim)
    i, j = _PROJ_AXES[projection]
    return pos[:, i], pos[:, j]


def luminosity_from_mag(phot, band_idx=1, M_sun_band=M_SUN_B):
    """Convierte magnitud absoluta (banda `band_idx` de GFM_StellarPhotometrics) a luminosidad L/Lsun."""
    return 10 ** ((M_sun_band - phot[:, band_idx]) / 2.5)


def make_synthetic_image(x, y, lum, pixel_scale_kpc=0.2, half_size_kpc=None, smooth_px=1.5):
    """
    Deposita la luminosidad en una grilla 2D y suaviza con un kernel
    gaussiano (imita el ruido de discretización por número finito de
    partículas / PSF observacional).

    Retorna: image (ny,nx, Lsun/pc^2), x0, y0 (centro en píxeles),
             pixel_scale_kpc, extent_kpc (para imshow).
    """
    if half_size_kpc is None:
        half_size_kpc = np.percentile(np.hypot(x, y), 99)
    n_pix = int(np.ceil(2 * half_size_kpc / pixel_scale_kpc))
    edges = np.linspace(-half_size_kpc, half_size_kpc, n_pix + 1)

    H, _, _ = np.histogram2d(x, y, bins=[edges, edges], weights=lum)
    pixel_area_pc2 = (pixel_scale_kpc * 1e3) ** 2
    image = H.T / pixel_area_pc2   # (ny, nx), Lsun/pc^2 — transpuesta para indexado [fila=y, col=x]

    if smooth_px > 0:
        image = gaussian_filter(image, sigma=smooth_px)

    x0 = y0 = n_pix / 2.0
    extent_kpc = [-half_size_kpc, half_size_kpc, -half_size_kpc, half_size_kpc]
    return image, x0, y0, pixel_scale_kpc, extent_kpc


def sb_from_sigma(sigma_Lsun_pc2, M_sun_band=M_SUN_B, z=0.0):
    """
    mu = M_sun_band + 21.572 - 2.5*log10(Sigma[Lsun/pc^2]) + 2.5*log10((1+z)^4)

    El término de dimming cosmológico es 0 para z=0 (snap 99, este proyecto);
    se deja como parámetro documentando cómo se sumaría a z>0 real.
    """
    sb_const = M_sun_band + 21.572 + 2.5 * np.log10((1 + z) ** 4)
    with np.errstate(divide='ignore', invalid='ignore'):
        mu = np.where(sigma_Lsun_pc2 > 0, sb_const - 2.5 * np.log10(sigma_Lsun_pc2), np.nan)
    return mu


# ─────────────────────────────────────────────────────────────────────────
# 2. Ajuste de isofotas elípticas
# ─────────────────────────────────────────────────────────────────────────

def fit_isophotes(image, x0, y0, sma0_px=10.0, eps0=0.3, pa0=0.0, step=0.1,
                   linear=False, minsma=0.0, maxsma=None, fix_center=True):
    """
    Ajusta isofotas elípticas (Jedrzejewski 1987) sobre `image` con
    photutils.isophote. `image` debe ser una intensidad (Sigma en Lsun/pc^2,
    no en magnitudes) para que el ajuste converja correctamente.

    Retorna un `IsophoteList` de photutils (sma en píxeles).
    """
    geometry = EllipseGeometry(x0=x0, y0=y0, sma=sma0_px, eps=eps0, pa=pa0,
                                fix_center=fix_center)
    ellipse = Ellipse(image, geometry=geometry)
    return ellipse.fit_image(step=step, linear=linear, minsma=minsma, maxsma=maxsma)


def _n_particles_in_annulus(x_part, y_part, a_kpc, eps, pa_rad, da_kpc):
    """Cuenta partículas dentro de un anillo elíptico (centrado en a_kpc, ancho da_kpc, orientado con eps/pa)."""
    cospa, sinpa = np.cos(pa_rad), np.sin(pa_rad)
    xr =  x_part * cospa + y_part * sinpa
    yr = -x_part * sinpa + y_part * cospa
    b_over_a = max(1 - eps, 1e-3)
    r_ell = np.sqrt(xr**2 + (yr / b_over_a)**2)
    return int(np.sum(np.abs(r_ell - a_kpc) < da_kpc / 2))


def isophote_table(isolist, x_part, y_part, pixel_scale_kpc, M_sun_band=M_SUN_B,
                    z=0.0, min_particles=50):
    """
    Convierte una IsophoteList de photutils a un DataFrame con columnas
    físicas: sma_kpc, a_circ_kpc (=sma*sqrt(1-eps)), eps, pa_rad, mu, mu_err,
    n_part (partículas dentro de la isofota — chequeo de ruido de Poisson) y
    `reliable` (n_part >= min_particles y mu finito).
    """
    sma_kpc = np.asarray(isolist.sma) * pixel_scale_kpc
    eps     = np.asarray(isolist.eps)
    pa      = np.asarray(isolist.pa)          # radianes
    a_circ  = sma_kpc * np.sqrt(np.clip(1 - eps, 0, 1))

    sigma     = np.asarray(isolist.intens)     # Lsun/pc^2 (mismas unidades que `image`)
    sigma_err = np.asarray(isolist.int_err)
    mu = sb_from_sigma(sigma, M_sun_band=M_sun_band, z=z)
    with np.errstate(divide='ignore', invalid='ignore'):
        mu_err = np.where(sigma > 0, (2.5 / np.log(10)) * np.abs(sigma_err / sigma), np.nan)

    n_part = np.array([
        _n_particles_in_annulus(x_part, y_part, a, e, p, pixel_scale_kpc)
        for a, e, p in zip(sma_kpc, eps, pa)
    ])
    reliable = (n_part >= min_particles) & np.isfinite(mu)

    return pd.DataFrame({
        'sma_kpc': sma_kpc, 'a_circ_kpc': a_circ, 'eps': eps, 'pa_rad': pa,
        'mu': mu, 'mu_err': mu_err, 'sigma_Lsun_pc2': sigma,
        'n_part': n_part, 'reliable': reliable,
    })


# ─────────────────────────────────────────────────────────────────────────
# 3. Ajuste de Sérsic y radio de Holmberg
# ─────────────────────────────────────────────────────────────────────────

def sersic_b_n(n):
    """Aproximación de Ciotti & Bertin (1999) para b_n del perfil de Sérsic."""
    return 1.9992 * n - 0.3271


def sersic_mu(a, mu_e, Re, n):
    """mu(a) = mu_e + (2.5 b_n / ln10) * [(a/Re)^(1/n) - 1]."""
    b_n = sersic_b_n(n)
    return mu_e + (2.5 * b_n / np.log(10)) * (np.abs(a / Re) ** (1.0 / n) - 1.0)


def fit_sersic(a_kpc, mu, mu_err=None, p0=None, bounds=((15, 0.2, 0.3), (32, 300, 10))):
    """
    Ajusta mu(a) al modelo de Sérsic con scipy.optimize.curve_fit.

    Retorna: popt=(mu_e, Re, n), perr=(mu_e_err, Re_err, n_err), pcov.
    """
    a_kpc, mu = np.asarray(a_kpc), np.asarray(mu)
    ok = np.isfinite(a_kpc) & np.isfinite(mu) & (a_kpc > 0)
    if mu_err is not None:
        mu_err = np.asarray(mu_err)
        ok &= np.isfinite(mu_err) & (mu_err > 0)
    a_fit, mu_fit = a_kpc[ok], mu[ok]
    sigma_fit = mu_err[ok] if mu_err is not None else None

    if len(a_fit) < 3:
        raise ValueError(
            f"Solo {len(a_fit)} isofotas confiables — se necesitan >= 3 para "
            "ajustar Sérsic. Revisa min_particles / pixel_scale_kpc.")

    lo, hi = np.asarray(bounds[0], dtype=float), np.asarray(bounds[1], dtype=float)
    if p0 is None:
        p0 = [np.median(mu_fit), np.median(a_fit), 2.0]
    p0 = np.clip(p0, lo + 1e-3 * (hi - lo), hi - 1e-3 * (hi - lo))

    popt, pcov = curve_fit(sersic_mu, a_fit, mu_fit, p0=p0, sigma=sigma_fit,
                            absolute_sigma=sigma_fit is not None,
                            bounds=bounds, maxfev=20000)
    perr = np.sqrt(np.diag(pcov))
    return popt, perr, pcov


def holmberg_radius_empirical(a_kpc, mu, mu_cut=MU_HOLMBERG):
    """Interpola el perfil medido mu(a) para hallar a(mu=mu_cut)."""
    a_kpc, mu = np.asarray(a_kpc), np.asarray(mu)
    valid = np.isfinite(mu) & (a_kpc > 0)
    if valid.sum() < 3:
        return np.nan
    a_v, m_v = a_kpc[valid], mu[valid]
    order = np.argsort(a_v)
    a_v, m_v = a_v[order], m_v[order]
    if m_v[0] > mu_cut or m_v[-1] < mu_cut:
        return np.nan
    f = interp1d(m_v, a_v, kind='linear', fill_value='extrapolate')
    r_h = float(f(mu_cut))
    return r_h if 0 < r_h <= a_v[-1] * 1.2 else np.nan


def holmberg_radius_model(mu_e, Re, n, mu_cut=MU_HOLMBERG):
    """R_H = Re * [1 + ln10/(2.5 b_n) * (mu_cut - mu_e)]^n."""
    b_n  = sersic_b_n(n)
    base = 1.0 + (np.log(10) / (2.5 * b_n)) * (mu_cut - mu_e)
    return Re * base ** n if base > 0 else np.nan


# ─────────────────────────────────────────────────────────────────────────
# 4. Orquestador de una proyección + robustez
# ─────────────────────────────────────────────────────────────────────────

def run_projection(pos_c, mass, phot, projection='xy', align_principal=False,
                    pixel_scale_kpc=0.2, half_size_kpc=None, smooth_px=1.5,
                    band_idx=1, M_sun_band=M_SUN_B, min_particles=50,
                    sma0_kpc=10.0, eps0=0.3, pa0=0.0, step=0.1, mu_cut=MU_HOLMBERG):
    """
    Corre el pipeline completo (imagen -> isofotas -> Sérsic -> R_H) para una
    proyección dada. Retorna un dict con todos los productos intermedios.
    """
    x, y = project_particles(pos_c, mass, projection=projection,
                              align_principal=align_principal)
    lum = luminosity_from_mag(phot, band_idx=band_idx, M_sun_band=M_sun_band)

    image, x0, y0, pxscale, extent = make_synthetic_image(
        x, y, lum, pixel_scale_kpc=pixel_scale_kpc,
        half_size_kpc=half_size_kpc, smooth_px=smooth_px)

    isolist = fit_isophotes(image, x0, y0, sma0_px=sma0_kpc / pxscale,
                             eps0=eps0, pa0=pa0, step=step)
    table = isophote_table(isolist, x, y, pxscale, M_sun_band=M_sun_band,
                            min_particles=min_particles)

    rel = table[table['reliable']]
    popt, perr, pcov = fit_sersic(rel['a_circ_kpc'].to_numpy(), rel['mu'].to_numpy(),
                                   rel['mu_err'].to_numpy())
    mu_e, Re, n = popt

    r_h_emp   = holmberg_radius_empirical(rel['a_circ_kpc'].to_numpy(), rel['mu'].to_numpy(), mu_cut=mu_cut)
    r_h_model = holmberg_radius_model(mu_e, Re, n, mu_cut=mu_cut)

    return {
        'projection': projection, 'align_principal': align_principal,
        'image': image, 'extent_kpc': extent, 'x0': x0, 'y0': y0,
        'pixel_scale_kpc': pxscale, 'smooth_px': smooth_px,
        'isolist': isolist, 'table': table,
        'mu_e': mu_e, 'Re': Re, 'n': n,
        'mu_e_err': perr[0], 'Re_err': perr[1], 'n_err': perr[2],
        'r_h_empirical': r_h_emp, 'r_h_model': r_h_model,
    }


def convergence_grid(pos_c, mass, phot, projection='xy',
                      pixel_scales=(0.1, 0.2, 0.3, 0.5), smooth_pxs=(1.0, 1.5, 2.0),
                      **kwargs):
    """
    Corre run_projection variando pixel_scale_kpc y smooth_px para chequear
    que R_H no dependa fuertemente de estos parámetros de discretización.
    Retorna un DataFrame con una fila por combinación (NaN si no converge).
    """
    rows = []
    for ps in pixel_scales:
        for sm in smooth_pxs:
            try:
                res = run_projection(pos_c, mass, phot, projection=projection,
                                      pixel_scale_kpc=ps, smooth_px=sm, **kwargs)
                rows.append({'pixel_scale_kpc': ps, 'smooth_px': sm,
                             'r_h_empirical': res['r_h_empirical'],
                             'r_h_model': res['r_h_model'],
                             'n': res['n'], 'Re_kpc': res['Re'], 'mu_e': res['mu_e']})
            except Exception as e:
                rows.append({'pixel_scale_kpc': ps, 'smooth_px': sm,
                             'r_h_empirical': np.nan, 'r_h_model': np.nan,
                             'n': np.nan, 'Re_kpc': np.nan, 'mu_e': np.nan,
                             'error': str(e)})
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────
# 5. Salidas: CSV, tabla resumen, figuras
# ─────────────────────────────────────────────────────────────────────────

def save_isophote_csv(table, path):
    table.to_csv(path, index=False)


def summary_row(result, sub_id, extra=None):
    """Fila de resumen (dict) a partir del resultado de run_projection."""
    r_h_e, r_h_m = result['r_h_empirical'], result['r_h_model']
    row = {
        'sub_id': sub_id, 'projection': result['projection'],
        'align_principal': result['align_principal'],
        'pixel_scale_kpc': result['pixel_scale_kpc'], 'smooth_px': result['smooth_px'],
        'mu_e': result['mu_e'], 'mu_e_err': result['mu_e_err'],
        'Re_kpc': result['Re'], 'Re_err_kpc': result['Re_err'],
        'n': result['n'], 'n_err': result['n_err'],
        'R_H_empirical_kpc': r_h_e, 'R_H_model_kpc': r_h_m,
        'delta_R_H_kpc': abs(r_h_e - r_h_m) if np.isfinite(r_h_e) and np.isfinite(r_h_m) else np.nan,
    }
    if extra:
        row.update(extra)
    return row


def plot_galaxy_isophotes(result, sub_id, fig_pdf_dir, fig_png_dir, mu_cut=MU_HOLMBERG, n_ellipses=6):
    """
    Figura de dos paneles: (izq) imagen sintética (mu_B) + isofotas
    superpuestas; (der) perfil mu(a) medido + modelo de Sérsic + línea en
    mu_cut, marcando ambos radios de Holmberg.
    """
    image  = result['image']
    extent = result['extent_kpc']
    table  = result['table']
    rel    = table[table['reliable']].reset_index(drop=True)

    mu_img = sb_from_sigma(image)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    ax = axes[0]
    im = ax.imshow(mu_img, origin='lower', cmap='ocean', extent=extent, vmin=18, vmax=30)
    plt.colorbar(im, ax=ax, label=r'$\mu_B$ [mag arcsec$^{-2}$]')
    idx = np.linspace(0, len(rel) - 1, min(n_ellipses, len(rel))).astype(int) if len(rel) else []
    for k in idx:
        a, eps, pa = rel.loc[k, ['sma_kpc', 'eps', 'pa_rad']]
        b = a * (1 - eps)
        patch = MplEllipse((0, 0), 2 * a, 2 * b, angle=np.degrees(pa),
                            fill=False, edgecolor='white', lw=1.0, alpha=0.85)
        ax.add_patch(patch)
    ax.set_xlabel('x [kpc]'); ax.set_ylabel('y [kpc]')
    ax.set_title(f"Imagen sintética + isofotas — subhalo {sub_id} ({result['projection']})")

    ax = axes[1]
    ax.errorbar(rel['a_circ_kpc'], rel['mu'], yerr=rel['mu_err'], fmt='o',
                ms=4, color='steelblue', ecolor='lightsteelblue', label='Isofotas (medido)')
    if len(rel):
        a_model  = np.linspace(rel['a_circ_kpc'].min(), rel['a_circ_kpc'].max(), 200)
        mu_model = sersic_mu(a_model, result['mu_e'], result['Re'], result['n'])
        ax.plot(a_model, mu_model, 'k-', lw=1.8,
                label=f"Sérsic: n={result['n']:.2f}, Re={result['Re']:.1f} kpc")
    ax.axhline(mu_cut, color='r', ls='--', lw=1.2, label=f'$\\mu_B$={mu_cut}')
    if np.isfinite(result['r_h_empirical']):
        ax.axvline(result['r_h_empirical'], color='tomato', ls=':',
                   label=f"R_H emp={result['r_h_empirical']:.1f} kpc")
    if np.isfinite(result['r_h_model']):
        ax.axvline(result['r_h_model'], color='seagreen', ls='-.',
                   label=f"R_H mod={result['r_h_model']:.1f} kpc")
    ax.invert_yaxis()
    ax.set_xlabel('Semieje mayor circularizado [kpc]')
    ax.set_ylabel(r'$\mu_B$ [mag arcsec$^{-2}$]')
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(f"{fig_pdf_dir}/fig04_isofotas_sub{sub_id}_{result['projection']}.pdf", bbox_inches='tight')
    plt.savefig(f"{fig_png_dir}/fig04_isofotas_sub{sub_id}_{result['projection']}.png", bbox_inches='tight', dpi=150)
    return fig
