"""
Prueba de integración: ejecuta COPIAS de los notebooks de la réplica contra el illustris_python simulado.
Uso:  /opt/anaconda3/bin/python run_mock.py <dir_de_trabajo> [notebooks...]
No toca el repo: copia params_icl.py, icl_funcs.py, original_shift_code/ y los notebooks a <dir_de_trabajo>.
"""
import os, sys, shutil, time
import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(os.path.dirname(HERE), 'mayes_2026_methodology')
WORK = os.path.abspath(sys.argv[1])
NBS = sys.argv[2:] or ['00a_chequeos_lab', '00_catalogo', '01_separacion_BCG_ICL', '01e_BCG_ICL_R200_local',
                       '02_propiedades_estelares', '03_analisis_progenitores']

os.makedirs(WORK, exist_ok=True)
for f in ['params_icl.py', 'icl_funcs.py', 'mary_haloshift_z0.pkl']:
    shutil.copy(os.path.join(REPO, f), WORK)
shutil.copytree(os.path.join(REPO, 'original_shift_code'), os.path.join(WORK, 'original_shift_code'), dirs_exist_ok=True)

MDIR = os.path.join(WORK, 'mock_tng')
os.environ.update({
    'PYTHONPATH': os.path.join(HERE, 'mock_il') + os.pathsep + os.environ.get('PYTHONPATH', ''),
    'PYTHONDONTWRITEBYTECODE': '1',
    'MOCK_TNG_DIR': MDIR,
    'ICL_REPO_DIR': WORK,
    'TNG_BASEPATH': os.path.join(MDIR, 'output') + '/',
    'TNG_SA_STARS_FILE': os.path.join(MDIR, 'stars_{snap:03d}.hdf5'),
    'MPLBACKEND': 'Agg',
})

for name in NBS:
    src = os.path.join(REPO, name + '.ipynb')
    nb = nbformat.read(src, as_version=4)
    t0 = time.time()
    client = NotebookClient(nb, timeout=1800, kernel_name='python3', resources={'metadata': {'path': WORK}})
    try:
        client.execute()
        status = 'OK'
    except CellExecutionError as e:
        status = 'ERROR\n' + str(e)[-2500:]
    nbformat.write(nb, os.path.join(WORK, name + '.ejecutado.ipynb'))
    print(f"==== {name}: {status}  ({time.time() - t0:.0f} s)")
    # prints de texto para revisar resultados
    for c in nb.cells:
        if c.cell_type != 'code':
            continue
        for o in c.get('outputs', []):
            if o.get('output_type') == 'stream':
                txt = o['text'].replace('\r', '\n')
                lines = [l for l in txt.splitlines() if l.strip()]
                for l in lines[-12:]:
                    print('   |', l[:160])
