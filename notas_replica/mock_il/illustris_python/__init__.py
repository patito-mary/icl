"""illustris_python SIMULADO (sólo pruebas de integración de la réplica Mayes+2026; ver _world.py)."""
from . import groupcat, snapshot, sublink

__version__ = 'mock-0.1'

from . import _world as _W
_W.ensure_files()      # en el Lab estos archivos ya existen; en el mock se crean al importar
