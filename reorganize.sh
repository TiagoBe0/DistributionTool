#!/usr/bin/env bash
# Reorganiza el proyecto DistributionTool.
# Ejecutar desde la raíz del proyecto: bash reorganize.sh
set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

echo "=== Reorganizando DistributionTool ==="

# ── 1. data/dumps/ ─────────────────────────────────────────────────────────────
mkdir -p data/dumps data/analyzed

for f in test_dmg.dump test_ref.dump \
          training_file.0.dump training_file.7.dump training_file.15.dump \
          training_file.27.dump training_file.35.dump; do
    [ -f "$f" ] && mv "$f" data/dumps/ && echo "  mv $f → data/dumps/"
done

for f in analyzed_test_dmg.dump "analyzed_dump-finalCool.160000"; do
    [ -f "$f" ] && mv "$f" data/analyzed/ && echo "  mv $f → data/analyzed/"
done

# ── 2. scripts/ ────────────────────────────────────────────────────────────────
mkdir -p scripts

for f in graficar_output.py graficar_pca.py publication_charts.py \
          visualizer.py analisis_nanoparticula.py; do
    [ -f "$f" ] && mv "$f" scripts/ && echo "  mv $f → scripts/"
done

# ── 3. figures/ ────────────────────────────────────────────────────────────────
mkdir -p figures

for f in dist_by_type.pdf dist_by_type.png \
          analisis_nanoparticula.pdf analisis_nanoparticula.png; do
    [ -f "$f" ] && mv "$f" figures/ && echo "  mv $f → figures/"
done

# ── 4. results/ ─────────────────────────────────────────────────────────────────
mkdir -p results

for f in output.csv output_hea.csv output_nano.csv outputmoiercoles.csv \
          test_output.csv \
          hist_output.csv hist_output_hea.csv hist_output_nano.csv \
          hist_test_output.csv hist_outputmoiercoles.csv \
          pca_output.csv pca_output_hea.csv pca_output_nano.csv \
          pca_test_output.csv pca_outputmoiercoles.csv \
          vacancies_output.csv vacancies_output_hea.csv \
          vacancies_test_output.csv vacancies_outputmoiercoles.csv; do
    [ -f "$f" ] && mv "$f" results/ && echo "  mv $f → results/"
done

echo ""
echo "=== Listo. Estructura resultante:"
ls -1
