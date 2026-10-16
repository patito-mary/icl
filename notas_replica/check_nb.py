import sys, ast, json, nbformat, re
from pyflakes.api import check as pfcheck
from pyflakes.reporter import Reporter
import io
for path in sys.argv[1:]:
    nb = nbformat.read(path, as_version=4)
    nbformat.validate(nb)
    src_all = []
    bad = 0
    for i, c in enumerate(nb.cells):
        if c.cell_type != 'code': continue
        s = '\n'.join(l for l in c.source.splitlines() if not l.lstrip().startswith(('%', '!')))
        try: ast.parse(s)
        except SyntaxError as e: bad += 1; print(f'  SINTAXIS {path} celda {i}: {e}')
        src_all.append(s)
    out, err = io.StringIO(), io.StringIO()
    pfcheck('\n\n'.join(src_all), path, Reporter(out, err))
    msgs = [m for m in out.getvalue().splitlines() if 'imported but unused' not in m]
    n_out = sum(len(c.get('outputs', [])) for c in nb.cells if c.cell_type == 'code')
    print(f'{path.split("/")[-1]:38s} celdas={len(nb.cells):3d} sintaxis_mal={bad} outputs={n_out} pyflakes={len(msgs)}')
    for m in msgs: print('   ', m.split(':',1)[1] if ':' in m else m)
