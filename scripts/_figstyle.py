#!/usr/bin/env python3
"""
Estilo compartido para las figuras del paper de la cascada PKA FeCrNi.

Fuente única de verdad para paleta + tipografía + dpi, de modo que las seis
figuras del manuscrito (survival_overview, cascade_anatomy, predictability_
horizon, ml_survival_weights, damage_evolution_*, energy_comparison) se vean
homogéneas. Importar y llamar setup() antes de crear figuras; guardar con
save(fig, "nombre_sin_extension") para emitir PNG (150 dpi) + PDF vectorial.

Uso:
    from _figstyle import setup, save, C
    setup()
    fig, ax = plt.subplots(...)
    ax.plot(..., color=C["ws"])
    save(fig, "mi_figura")
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "figures")

# Paleta semántica canónica. Dos ejes distintos, no mezclar:
#  - comparación de métodos de conteo: ws / hybrid / grid
#  - clasificación de supervivencia:    survive / transient / traj
C = {
    "ws":        "#c0392b",  # Wigner-Seitz / pares de Frenkel
    "hybrid":    "#1f6fb2",  # detector híbrido (preset robust/survival)
    "grid":      "#27ae60",  # grid / voids por volumen
    "survive":   "#2a9d8f",  # vacancias que sobreviven al relax
    "transient": "#e76f51",  # vacancias transitorias del pico
    "traj":      "#1d7d6b",  # trayectorias de sobrevivientes
    "muted":     "0.5",      # líneas auxiliares (azar, medianas, flechas)
}


def setup():
    """Aplica tipografía y dpi uniformes a todas las figuras."""
    plt.rcParams.update({
        "font.size":        11,
        "axes.titlesize":   12,
        "axes.labelsize":   11,
        "legend.fontsize":  9,
        "xtick.labelsize":  9,
        "ytick.labelsize":  9,
        "figure.dpi":       150,
        "savefig.dpi":      150,
        "legend.frameon":   False,
    })


def save(fig, name):
    """Guarda figures/<name>.png (150 dpi) y figures/<name>.pdf (vectorial)."""
    os.makedirs(FIGDIR, exist_ok=True)
    base = os.path.join(FIGDIR, name)
    fig.savefig(base + ".png", dpi=150, bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    return base + ".png"
