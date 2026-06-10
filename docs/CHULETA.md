# DistributionTool — Chuleta de comandos

---

## COMPILAR

### distool (análisis CLI)
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
# binario: build/distool
```

### distool_viewer (visualizador OpenGL)
```bash
cmake -S viewer -B build_viewer -DCMAKE_BUILD_TYPE=Release
cmake --build build_viewer -j$(nproc)
# binario: build_viewer/distool_viewer
```

> Si falla por dependencias X11:
> `sudo apt-get install -y libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev`

---

## EJECUTAR ANÁLISIS

### HEA — dump-finalCool (bulk ~1M átomos, FCC)
```bash
./build/distool \
    data/dumps/dump-finalCool.0 data/dumps/dump-finalCool.160000 \
    --n-max 9 --l-max 9 --r-cut 5.0 \
    --threshold 0.15 --grid-spacing 1.0 \
    --pca 3 --hist \
    --output output_hea.csv \
    2>&1 | tee run_hea.log
```

### Nanopartícula BCC (~31k átomos, r_cut corto, vacías desactivadas)
```bash
# la referencia pristina vive en tests/data/; el frame dañado en data/dumps/
./build/distool \
    tests/data/dump.nanoparticula.0 data/dumps/dump.nanoparticula.8000 \
    --n-max 9 --l-max 9 --r-cut 3.5 \
    --threshold 0.15 --grid-spacing 200.0 \
    --pca 2 --hist \
    --output output_nano.csv \
    2>&1 | tee run_nano.log
```
> `--grid-spacing 200` desactiva la detección de vacancias por grilla
> (no apta para nanopartículas con superficie libre)

### Análisis rápido / prueba
```bash
./build/distool \
    test_ref.dump test_dmg.dump \
    --n-max 5 --l-max 4 --r-cut 4.5 \
    --threshold 0.15 \
    --output test_output.csv
```

---

## WORKFLOW: GUARDAR DV DE UN ÁTOMO DEFECTUOSO
```bash
# 1. Guardar el DV del átomo id=42 del frame dañado
./build/distool ref.dump damaged.dump \
    --save-dv 42 sia_ref.dat

# 2. Reclasificar con esa referencia
./build/distool ref.dump damaged.dump \
    --ref-sia  sia_ref.dat \
    --output output_con_refs.csv
```

---

## VISUALIZADOR

### Abrir un CSV o dump analizado directamente
```bash
./build_viewer/distool_viewer output_hea.csv
./build_viewer/distool_viewer analyzed_dump-finalCool.160000
```

### Abrir sin archivo (usar pestaña Analizar desde la UI)
```bash
./build_viewer/distool_viewer
```

**Controles:**
| Acción         | Control              |
|----------------|----------------------|
| Rotar          | Drag botón izquierdo |
| Pan            | Drag botón medio     |
| Zoom           | Scroll               |
| Seleccionar átomo | Click izquierdo   |

---

## SCRIPTS PYTHON

```bash
# Visualizar distribución 3D + histograma de defect_prob
python3 scripts/graficar_output.py          # usa output.csv

# PCA scatter + boxplot por tipo
python3 scripts/graficar_pca.py             # usa pca_output.csv

# Análisis completo nanopartícula vs Wigner-Seitz
python3 scripts/analisis_nanoparticula.py   # usa output_nano.csv + dump WS

# Figuras de publicación (HEA)
python3 scripts/publication_charts.py
```

---

## ARCHIVOS DE SALIDA

| Archivo                  | Contenido                                      |
|--------------------------|------------------------------------------------|
| `output.csv`             | id, tipo, xyz, dist_to_ref, defect_prob, label |
| `vacancies_output.csv`   | Vacancias: xyz, d_near_max, n_pts              |
| `hist_output.csv`        | Histograma: bin_centre, count                  |
| `pca_output.csv`         | PCA: id, tipo, pc1, pc2…, dist_to_ref, label   |
| `analyzed_<dump>`        | Dump LAMMPS con columnas extra (OVITO-ready)   |

---

## PARÁMETROS CLAVE

| Material      | Estructura | r_cut sugerido | threshold |
|---------------|------------|---------------|-----------|
| Fe, HEA FeCrNiCoCu | FCC   | 5.0 Å         | 0.15      |
| BCC (Fe-like) | BCC        | 3.5–4.5 Å     | 0.15      |
| W             | BCC        | 4.5–5.5 Å     | 0.15      |
| Cu, Ni        | FCC        | 4.0–5.0 Å     | 0.15      |

**DV size** = `n_max*(n_max+1)/2 * (l_max+1)`
→ con n=9, l=9: **450 componentes**

---

## NOTAS RÁPIDAS

- El **threshold 0.15** separa red (Lattice) de defectos; subilo a 0.20-0.30 si hay muchos FP
- **dist_to_ref alto** (~0.8–1.0) → átomo muy distorsionado (defecto confirmado)
- **Sin `--ref-sia/antv`** → todos los defectos salen como `Unknown`
- Átomos de **superficie** siempre aparecen como Unknown/defecto (CN bajo → DV distorsionado)
- Para nanopartículas: siempre usar `--grid-spacing 200` para desactivar grilla de vacancias
