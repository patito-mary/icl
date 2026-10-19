# Avance de la réplica de Mayes+2026 (rama replica-mayes2026, sin commits)

Actualizado: 2026-10-16, segunda sesión. Nada se ha corrido todavía con datos de TNG.

## Estado de los archivos en `mayes_2026_methodology/`

**Fuente de verdad.** Los scripts de `nbbuild/build_*.py` generan los notebooks:
- `build_00`, `build_00a`, `build_01`, `build_01e`, `build_02`, `build_03` reescriben el notebook completo.
- `build_01c` y `build_01f` editan in situ: NO hay que re-ejecutarlos (aplicarían la edición dos veces; el de 01c fallaría en un assert).

**Validación.**
- `nbformat.validate`, `ast` y pyflakes OK en todo; los notebooks no tienen outputs.
- Las celdas de API de 01f siguen idénticas a HEAD.
- El script de chequeo está en `notas_replica/check_nb.py`.

**Archivos:**
- `params_icl.py`, `icl_funcs.py`: terminados. Nuevo método de R_H `'ajuste_local'` y parámetro `RH_LOCAL_MU_WINDOW`.
- `00a_chequeos_lab.ipynb` (NUEVO): sólo lee metadatos y acumula fallas. Probado con el mock.
- `00`, `01`, `01e`, `02`, `03`: terminados y probados de punta a punta con el mock.
- `01c`: probado con el mock.
- `01f`: no se puede probar con el mock porque depende de la API web.
- `01b` y `01d`: eliminados (los reemplazan 01e y 01f).

## Prueba de integración: `notas_replica/mock_il` + `run_mock.py`

Uso:

```
/opt/anaconda3/bin/python notas_replica/run_mock.py <dir_temporal> [notebooks]
```

El script ejecuta COPIAS de los notebooks contra un `illustris_python` simulado. El universo sintético está descrito en `mock_il/illustris_python/_world.py`.

**Resultados contra las verdades inyectadas:**
- Muestra: N = 3 con el corte 1e3. La periodicidad funciona (el FoF 0 cruza el borde de la caja).
- Estado dinámico: offsets 0.050 / 0.251 / 0.552 → estados 0/1/2.
- Concentración c: 5.00 / 6.87 / 8.85 (verdad 5 / 7 / 9).
- Major merger: snap 90 (razón 0.30), ninguno, snap 70 (razón 0.25), como se inyectó.
- 03:
  - snapshot de despegue y subhalo de origen correctos en el 100 % de las partículas;
  - los subhalos anómalos se reemplazan por la muestra 2 (snap 73);
  - las anómalas en ambas muestras se excluyen.
- `global_subhalo` de 03 se probó aparte: coincide con la verdad, incluido el fuzz.
- Gradientes de 02: negativos, como se inyectó.

**No probado:**
- 01f (API).
- Rendimiento y memoria con TNG real.
- Un MPB con snapshots saltados (el código los salta, pero el mock no los tiene).
- `loadTree` devolviendo `None`.
- La estructura real del catálogo SA por partícula.
- Matplotlib: `plt.cm.get_cmap` (en la celda original de Mary en 01c) avisa deprecación con mpl 3.10.

**Bugs de la réplica encontrados con el mock:** ninguno. Sólo hubo que arreglar el propio mock.

## Decisión sobre R_H (sección B)

Digitalicé la Fig. 1 del paper (`fig1.png`, 8× desde el PDF, ticks calibrados):
- La línea roja está en 200.6 kpc.
- Ahí los puntos medidos están en μ ≈ 26.15 y la línea del ajuste en ≈ 26.55; a 220 kpc el ajuste está en 26.70.
- Los datos cruzan 26.5 recién hacia ~230 kpc.

Conclusión:
- R_H en el paper = cruce del ajuste logarítmico GLOBAL. Se mantiene `RH_METHOD = 'ajuste'` por defecto.
- Ese R_H queda ~13–15 % por dentro del cruce de los datos: el mismo efecto que dio la prueba Sérsic (−15 a −17 %).
- El ajuste local (26.5 ± 1 mag) reproduce el cruce verdadero en ±2 % (prueba Sérsic) y queda como opción.
- 01 guarda `R_H_local_kpc` e imprime la razón ajuste/local; en la Fig. 1 del paper esa razón es ≈ 0.87.

## SUPUESTOS ordenados por impacto en f_ICL (ver el informe)

1. Binning del perfil (lineal 5 kpc) y rango del ajuste: definen el peso de cada radio en el ajuste global, y por lo tanto R_H.
2. Límite externo de R200 en 3D (no proyectado).
3. Definición de M*200 (estrellas del FoF dentro de R200 3D).
4. Eje de proyección (línea de visión = eje menor).
5. Ajuste sin pesos.
6. Corte en M☉/h: cambia qué halos entran, no f_ICL de cada halo.

Bajo impacto: tensor de inercia reducido, exclusión de viento, sin polvo/K/(1+z)⁴.

No afectan f_ICL: estado dinámico, NFW, merger, catálogo SA, cuña, progenitores.

## Orden en el Lab

00a → 00 → 01 → 01e → 01f (opcional) → 02 → 03 (primero con `N_PRUEBA = 2`).

Validar en 01:
- el cúmulo demo (3.º más masivo) debería dar R_H ≈ 200 kpc si es el de la Fig. 1;
- mediana de f_ICL ≈ 0.10.

## 2026-10-19 · Error `header['MassTable']` en el Lab

**Causa:**
- 00 (celda de concentración NFW) leía `Header['MassTable']` del header que devuelve `il.groupcat.loadHeader`.
- En illustris_python 1.0.0 esa función hace `dict(f['Header'].attrs)` sobre `groups_NNN/fof_subhalo_tab_NNN.0.hdf5`, el Header del catálogo de grupos.
- En TNG ese Header no tiene `MassTable`, que vive en el Header del snapshot. Tampoco hay `il.snapshot.loadHeader` pública.
- El mock de la primera prueba ponía `MassTable` en el header del groupcat, y por eso no lo detectó.

**Arreglo:**
- `icl_funcs`:
  - `snapshot_header()` lee el Header de `snapdir_NNN/snap_NNN.0.hdf5`;
  - `n_stars_snapshot()`;
  - `dm_particle_mass_msun()` usa `MassTable[1]` × UM, con respaldo `P.M_DM_MSUN_FALLBACK = 7.5e6` y aviso si se usa.
- Cambios por notebook:
  - 00 usa `F.dm_particle_mass_msun`;
  - 02 usa `F.snapshot_header`;
  - 00a valida las claves de AMBOS headers (falla con la lista de claves disponibles) y comprueba m_DM entre 5e6 y 1e7 Msun.
- El mock ahora es fiel: el header del groupcat no trae `MassTable` y el Header del snapshot sí. Con el código viejo, el mock reproduce el `KeyError`; con el arreglo, toda la batería pasa.
- `common.write_nb` ahora reutiliza los ids de celda si el número de celdas no cambia, para que los diffs queden limpios.

Las demás lecturas de `il.groupcat.loadHeader` usan sólo BoxSize, Redshift y Time, que sí están en el Header del groupcat.

## 2026-10-19 (b) · Verificación contra la documentación oficial de TNG

Lo que dice la documentación:
- Doc `data/docs/specifications`, Header del snapshot: BoxSize, `MassTable` ("Masses of particle types which have a constant mass (only DM)"), NumPart_ThisFile, NumPart_Total, NumPart_Total_HighWord, Omega0, OmegaLambda, Redshift, Time, NumFilesPerSnapshot.
- Doc, Header de fof_subhalo_tab: Ngroups_*, Nsubgroups_*, NumFiles, Time, Redshift. Sin MassTable.
- Doc, PartType1: no tiene `Masses`.
- Página `data/downloads/TNG100-1/`: "dark matter particle mass [10^10 msun/h] 0.000505574296436975" (≈ 7.46e6 Msun); box 75000 ckpc/h; h = 0.6774.

Contraste con los archivos reales:
- La tabla de la doc para el Header del groupcat NO lista BoxSize, pero los outputs viejos de Mary en el Lab muestran que `il.groupcat.loadHeader(...)['BoxSize']` = 75000. La tabla de la doc está incompleta en ese punto.
- HubbleParam no aparece en ninguna de las dos tablas: 00a ya no lo exige.

Conclusión: la corrección anterior era correcta en el diagnóstico (MassTable no está en el Header del groupcat; está en el del snapshot).

Cambios:
- `dm_particle_mass_msun` ahora devuelve `(m, fuente)`: usa MassTable[1] del snapshot si es > 0; si no, el valor DOCUMENTADO `P.M_DM_INTERNAL_DOC` (ya no se trata como respaldo genérico, que antes valía 7.5e6), declarando la fuente; avisa si ambos difieren.
- 00a imprime las claves de ambos Headers, el MassTable del snapshot y la fuente usada para m_DM.
- El mock es fiel a la doc. La variante `MOCK_SIN_MASSTABLE=1` prueba el camino del valor documentado.

## 2026-10-19 (c) · Primera corrida de 00 en el Lab (respaldo con outputs: /Users/mvs/Documents/respaldo_00_catalogo_corrida_lab_20261019.ipynb, fuera del repo)

Resultados de esa corrida:
- N = 127 halos, 10 sobre 1e14 Msun/h, log M200c entre 13.00 y 14.41: igual que el paper.
- R200c entre 517 y 1529 kpc; log M*200 entre 11.34 y 12.91.
- Estado dinámico: 108 relajados, 19 intermedios, **0 perturbados**; mediana de D/R200 = 0.052.
  - El paper sí tiene cúmulos perturbados: en §5.1, 9 de 27 en la selección alta y 4 de 27 en la baja, entre otros.
  - Sospecha: el SUPUESTO de CM del gas dentro de R200 diluye los offsets. Probar el gas de todo el FoF.
- Falló en la celda de m_DM con `AttributeError: module 'icl_funcs' has no attribute 'dm_particle_mass_msun'`, aunque HEAD (60e82b1) sí tiene esa función. Probablemente el kernel tenía importado un `icl_funcs` viejo: hay que reiniciar el kernel tras cada pull.

Corrección mínima de 00:
- Sólo cambió la celda de m_DM a la v2 (diff contra HEAD: 2 líneas).
- Se limpiaron outputs y execution_count; la metadata quedó igual a HEAD; se conservaron los ids.
- Mary no había cambiado ningún código.
