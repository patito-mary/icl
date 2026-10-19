"""Helpers para construir los notebooks de la réplica con nbformat (celdas sin outputs)."""
import json, os, uuid
import nbformat
from nbformat.v4 import new_notebook, new_code_cell, new_markdown_cell

REPO = '/Users/mvs/surface_brightness_profile/mayes_2026_methodology'

KERNEL_META = {
    'kernelspec': {'display_name': 'Python 3 (ipykernel)', 'language': 'python', 'name': 'python3'},
    'language_info': {'name': 'python'},
}


def _id():
    return uuid.uuid4().hex[:8]


def md(text, cid=None):
    c = new_markdown_cell(text.strip('\n'))
    c['id'] = cid or _id()
    return c


def code(text, cid=None):
    c = new_code_cell(text.strip('\n'))
    c['id'] = cid or _id()
    return c


def write_nb(name, cells, metadata=None):
    path = os.path.join(REPO, name)
    # reutiliza los ids del notebook existente si tiene el mismo número de celdas (diffs limpios)
    if os.path.exists(path):
        old = nbformat.read(path, as_version=4).cells
        if len(old) == len(cells):
            for c, o in zip(cells, old):
                if c['cell_type'] == o['cell_type']:
                    c['id'] = o['id']
    nb = new_notebook(cells=cells, metadata=metadata or KERNEL_META)
    nb.nbformat, nb.nbformat_minor = 4, 5
    nbformat.validate(nb)
    with open(path, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)
    print(f'escrito {path}  ({len(cells)} celdas)')


# ── Celdas comunes de cabecera ───────────────────────────────────────────
REPO_CELL = r'''
# ── Carpeta del repo: el kernel del Lab arranca en /home/tnguser ─────────
import os, sys
from pathlib import Path

REPO_CANDIDATOS = [Path(os.environ.get('ICL_REPO_DIR', Path.cwd())),
                   Path.home() / 'icl_FS' / 'mayes_2026_methodology',
                   Path('/Users/mvs/surface_brightness_profile/mayes_2026_methodology')]
REPO_DIR = next((p.resolve() for p in REPO_CANDIDATOS if (p / 'params_icl.py').is_file()), None)
if REPO_DIR is None:
    raise FileNotFoundError(f"No encuentro params_icl.py en {REPO_CANDIDATOS}; define ICL_REPO_DIR")

os.chdir(REPO_DIR)
for p in [REPO_DIR, REPO_DIR / 'original_shift_code']:
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
print(f"Carpeta de trabajo: {os.getcwd()}")
'''

STYLE_CELL = r'''
%matplotlib inline
my_colors = ["#4587f1", "#5fe0ca", "#aa89e9", "#f2cc8f", "#f1b4b4"]
sns.set_theme(style="whitegrid", palette=my_colors)

# estado dinámico: 0 relajado, 1 intermedio, 2 perturbado
COLORS_STATE = {0: my_colors[0], 1: my_colors[3], 2: my_colors[4]}
# componentes: 0 in situ, 1 mergers completados, 2 sobrevivientes
COLORS_COMP  = {0: my_colors[4], 1: my_colors[2], 2: my_colors[1]}
'''

COSMO_CELL = r'''
cosmo = FlatLambdaCDM(H0=P.H0, Om0=P.Om0, Ob0=P.Ob0)
print(f"basePath = {P.basePath}")
print(f"SNAP     = {P.SNAP}")
print(f"h        = {P.h}   UL = {P.UL:.4f}   UM = {P.UM:.3e}")
'''


def fig_dirs_cell(nbname):
    return f'''
FIG_PDF = './figuras/{nbname}/pdf'
FIG_PNG = './figuras/{nbname}/png'
os.makedirs(FIG_PDF, exist_ok=True)
os.makedirs(FIG_PNG, exist_ok=True)
'''


MASS_LABEL_CELL = r'''
# Masas de halo en la convención que reproduce la muestra del paper (ver params_icl.py)
MASS_UNIT  = r'M_\odot\,h^{-1}' if P.UM_PLOT == P.UM_H else r'M_\odot'
M200_plot  = M200c / P.UM * P.UM_PLOT
lM200      = np.log10(M200_plot)
'''
