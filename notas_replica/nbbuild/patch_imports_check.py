"""
Parche mínimo (conserva ids): (1) en la celda de rutas, el repo queda PRIMERO en sys.path
(antes original_shift_code se anteponía al repo); (2) tras la celda de imports se inserta una
celda que recarga P y F, imprime F.__file__ y verifica las funciones que necesita cada notebook.
Idempotente: si la celda de verificación ya existe, la reemplaza.
"""
import os, uuid
import nbformat

REPO = '/Users/mvs/surface_brightness_profile/mayes_2026_methodology'

OLD_PATH = """for p in [REPO_DIR, REPO_DIR / 'original_shift_code']:
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))"""
NEW_PATH = """# el repo va PRIMERO en sys.path (gana sobre cualquier otro icl_funcs/params_icl)
for p in [REPO_DIR / 'original_shift_code', REPO_DIR]:
    if str(p) in sys.path:
        sys.path.remove(str(p))
    sys.path.insert(0, str(p))"""

NECESITA = {
    '00a_chequeos_lab': (['snapshot_header', 'n_stars_snapshot', 'dm_particle_mass_msun'], ['M_DM_INTERNAL_DOC']),
    '00_catalogo': (['periodic_delta', 'center_of_mass', 'dynamical_state', 'fit_nfw_concentration',
                     'last_major_merger_snap', 'dm_particle_mass_msun'], ['M_DM_INTERNAL_DOC']),
    '01_separacion_BCG_ICL': (['load_central_stars', 'bcg_icl_from_stars', 'mu_from_sigma', 'linfit_beta'],
                              ['RH_METHOD', 'RH_LOCAL_MU_WINDOW']),
    '01e_BCG_ICL_R200_local': (['load_central_stars', 'load_fof_stars', 'bcg_icl_from_stars'], ['RH_METHOD']),
    '02_propiedades_estelares': (['load_central_stars', 'load_sa_stars', 'classify_components',
                                  'component_fractions', 'snapshot_header', 'n_stars_snapshot',
                                  'radial_gradient', 'wedge_mask'], ['SA_STARS_FILE']),
    '03_analisis_progenitores': (['load_central_stars', 'load_sa_stars', 'classify_components',
                                  'subhalo_of_local_index', 'detach_update', 'shared_fraction',
                                  'significant_progenitors', 'half_mass_radius'], ['PROG_OUT']),
}

MARCA = '# ── Verificación de módulos'


def check_src(funcs, pvars):
    return f"""{MARCA} (tras un git pull basta re-ejecutar desde aquí) ──
import importlib
P, F = importlib.reload(P), importlib.reload(F)
print(f"icl_funcs  : {{F.__file__}}\\nparams_icl : {{P.__file__}}")
faltan = [f for f in {funcs!r} if not hasattr(F, f)] + [f'P.{{v}}' for v in {pvars!r} if not hasattr(P, v)]
if faltan or os.path.realpath(os.path.dirname(F.__file__)) != os.path.realpath(REPO_DIR):
    raise ImportError(f"Módulo desactualizado o de otra carpeta (faltan {{faltan}}; F = {{F.__file__}}). "
                      "Ejecuta `git pull` en ~/icl_FS, confirma `git log --oneline -3` y REINICIA el kernel.")"""


def patch(name):
    path = os.path.join(REPO, name + '.ipynb')
    nb = nbformat.read(path, as_version=4)
    i_repo = next(i for i, c in enumerate(nb.cells) if c.cell_type == 'code' and 'REPO_CANDIDATOS' in c.source)
    c = nb.cells[i_repo]
    if OLD_PATH in c.source:
        c.source = c.source.replace(OLD_PATH, NEW_PATH)
        c.outputs, c.execution_count = [], None
    assert NEW_PATH in c.source, f'{name}: celda de rutas inesperada'
    if name in NECESITA:
        i_imp = next(i for i, c in enumerate(nb.cells)
                     if c.cell_type == 'code' and 'import icl_funcs as F' in c.source)
        src = check_src(*NECESITA[name])
        existing = [i for i, c in enumerate(nb.cells) if c.cell_type == 'code' and c.source.startswith(MARCA)]
        if existing:
            cell = nb.cells[existing[0]]
            cell.source, cell.outputs, cell.execution_count = src, [], None
        else:
            cell = nbformat.v4.new_code_cell(src)
            cell['id'] = uuid.uuid4().hex[:8]
            nb.cells.insert(i_imp + 1, cell)
    nbformat.validate(nb)
    nbformat.write(nb, path)
    print(f'parcheado {name}')


if __name__ == '__main__':
    for n in ['00a_chequeos_lab', '00_catalogo', '01_separacion_BCG_ICL', '01e_BCG_ICL_R200_local',
              '02_propiedades_estelares', '03_analisis_progenitores', '01f_api_tng_BCG_ICL_R200',
              '01c_perfil_subhalos_grupo']:
        patch(n)
