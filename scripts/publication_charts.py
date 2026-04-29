#!/usr/bin/env python3
"""
publication_charts.py — Gráficos de calidad de publicación para DistributionTool
==================================================================================
Genera figuras listas para publicación (PDF + PNG 300 DPI) a partir de los
archivos CSV exportados por distool.

Formatos CSV soportados:
  output*.csv           → 6 paneles: distribuciones, mapas espaciales, clasificación
  pca_output*.csv       → espacio PCA coloreado por tipo + distribución dist_to_ref
  hist_output*.csv      → histograma de distancias a referencia (lineal + semilog)
  vacancies_output*.csv → proyecciones espaciales de vacancias detectadas

Uso:
  python publication_charts.py output.csv
  python publication_charts.py pca_output.csv --outdir figures/
  python publication_charts.py output_nano.csv --threshold 0.42 --label "Nanopartícula"
  python publication_charts.py --all            (procesa todos los CSV del directorio)
  python publication_charts.py --all --outdir pub_figures/ --threshold 0.15

Dependencias:
  pip install numpy pandas matplotlib scipy
  (seaborn opcional pero recomendado)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

try:
    from scipy.stats import gaussian_kde
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


# ─── Estilo de publicación ────────────────────────────────────────────────────

def set_publication_style():
    """Aplica rcParams para gráficos de calidad de publicación (estilo Nature/ACS)."""
    mpl.rcParams.update({
        # Figura
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'figure.facecolor': 'white',
        'figure.edgecolor': 'white',
        # Fuentes
        'font.family': 'serif',
        'font.serif': ['DejaVu Serif', 'Times New Roman', 'Times', 'serif'],
        'font.size': 9,
        'axes.titlesize': 10,
        'axes.labelsize': 9,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'legend.title_fontsize': 8,
        # Ejes
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': True,
        'grid.alpha': 0.22,
        'grid.linewidth': 0.5,
        'grid.color': '#888888',
        # Ticks
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'xtick.minor.size': 2,
        'ytick.minor.size': 2,
        # Líneas y marcadores
        'lines.linewidth': 1.2,
        'lines.markersize': 4,
        # Leyenda
        'legend.framealpha': 0.88,
        'legend.edgecolor': '#cccccc',
        'legend.frameon': True,
        # Texto matemático
        'mathtext.fontset': 'cm',
        # Guardar
        'pdf.fonttype': 42,    # fuentes embebidas en PDF (requerido por muchas revistas)
        'ps.fonttype': 42,
    })


# Paleta de colores para tipos de defectos (consistente entre figuras)
DEFECT_COLORS = {
    'Lattice':      '#2166ac',   # azul
    'Unknown':      '#d73027',   # rojo
    'Interstitial': '#f46d43',   # naranja
    'VacancyAdj':   '#fdae61',   # amarillo-naranja
    'TypeA':        '#a6d96a',   # verde claro
    'Vacancy':      '#762a83',   # morado
    'Surface':      '#80cdc1',   # teal
}


# ─── Carga de datos ───────────────────────────────────────────────────────────

def detect_csv_type(path: Path) -> str:
    """Detecta el tipo de CSV por nombre de archivo y, si es necesario, por cabecera."""
    name = path.stem.lower()
    if name.startswith('vacancies'):
        return 'vacancies'
    if name.startswith('hist'):
        return 'hist'
    if name.startswith('pca'):
        return 'pca'
    if name.startswith('output'):
        return 'output'
    # Fallback: leer cabecera
    with open(path) as f:
        header = f.readline().strip().lstrip('#').split()
    if 'pc1' in header:
        return 'pca'
    if 'bin_centre' in header:
        return 'hist'
    if 'defect_prob' in header:
        return 'output'
    if set(header) <= {'x', 'y', 'z'}:
        return 'vacancies'
    return 'unknown'


def load_output_csv(path: Path) -> pd.DataFrame:
    cols = ['id', 'type', 'x', 'y', 'z', 'dist_to_ref', 'defect_prob', 'defect_type']
    df = pd.read_csv(path, sep=r'\s+', comment='#', header=None, names=cols)
    df['defect_type'] = df['defect_type'].astype(str)
    return df


def load_pca_csv(path: Path) -> pd.DataFrame:
    cols = ['id', 'type', 'pc1', 'pc2', 'dist_to_ref', 'defect_type']
    df = pd.read_csv(path, sep=r'\s+', comment='#', header=None, names=cols)
    df['defect_type'] = df['defect_type'].astype(str)
    return df


def load_hist_csv(path: Path) -> pd.DataFrame:
    cols = ['bin_centre', 'count']
    return pd.read_csv(path, sep=r'\s+', comment='#', header=None, names=cols)


def load_vacancies_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r'\s+', comment='#', header=None)
    col_names = ['x', 'y', 'z', 'd_near_max', 'n_grid_pts']
    df.columns = col_names[:len(df.columns)]   # safe for old 3-col and new 5-col files
    return df


# ─── Utilidades gráficas ──────────────────────────────────────────────────────

def get_color_map(defect_types):
    """
    Devuelve un dict {tipo: color} para los tipos dados.
    Si algún tipo no está en la paleta fija, asigna colores de tab10.
    """
    unknown_types = [t for t in defect_types if t not in DEFECT_COLORS]
    extra_colors = plt.cm.tab10(np.linspace(0, 0.9, max(len(unknown_types), 1)))
    cmap = dict(DEFECT_COLORS)
    for t, c in zip(unknown_types, extra_colors):
        cmap[t] = c
    return {t: cmap[t] for t in defect_types}


def add_colorbar(fig, ax, scatter, label='', fraction=0.046, pad=0.04):
    cbar = fig.colorbar(scatter, ax=ax, fraction=fraction, pad=pad)
    cbar.set_label(label, fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    return cbar


def kde_curve(values, x_grid, bw_method='scott'):
    """Retorna la curva KDE evaluada en x_grid, o None si falla."""
    if not HAS_SCIPY or len(values) < 5:
        return None
    try:
        kde = gaussian_kde(values, bw_method=bw_method)
        return kde(x_grid)
    except Exception:
        return None


def save_figure(fig, stem: str, outdir: Path, formats=('pdf', 'png')):
    outdir.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        p = outdir / f"{stem}.{fmt}"
        dpi = 300 if fmt == 'png' else None
        fig.savefig(p, bbox_inches='tight', dpi=dpi)
        print(f"    ✓ {p}")


def _sample(df: pd.DataFrame, n: int, seed: int = 42) -> pd.DataFrame:
    return df if len(df) <= n else df.sample(n, random_state=seed)


# ─── Figura 1: output*.csv ────────────────────────────────────────────────────

def plot_output(df: pd.DataFrame, label: str, threshold: float,
                outdir: Path, stem: str):
    """
    6-panel overview para output*.csv:
      [0,0] Histograma + KDE de dist_to_ref por tipo de defecto
      [0,1] Proyección XY coloreada por dist_to_ref
      [1,0] Violin / boxplot de dist_to_ref por tipo
      [1,1] Proyección XZ coloreada por tipo de defecto
      [2,0] Barras: conteo y porcentaje por tipo
      [2,1] Distribución de defect_prob
    """
    set_publication_style()

    defect_types = sorted(df['defect_type'].unique())
    cmap = get_color_map(defect_types)
    sample50k = _sample(df, 50_000)

    fig = plt.figure(figsize=(10, 11))
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.58, wspace=0.40)

    # ── [0,0] Histograma + KDE dist_to_ref ────────────────────────────────────
    ax0 = fig.add_subplot(gs[0, 0])
    xmax = df['dist_to_ref'].quantile(0.999)
    bins = np.linspace(0, xmax, 80)
    x_grid = np.linspace(0, xmax, 400)

    for dt in defect_types:
        vals = df[df['defect_type'] == dt]['dist_to_ref'].values
        ax0.hist(vals, bins=bins, alpha=0.55, color=cmap[dt],
                 density=True, histtype='stepfilled', edgecolor='none')
        kde_y = kde_curve(vals, x_grid)
        if kde_y is not None:
            ax0.plot(x_grid, kde_y, color=cmap[dt], lw=1.4)

    ax0.axvline(threshold, color='#333333', ls='--', lw=1.2,
                label=f'Umbral = {threshold}')
    handles = [Patch(facecolor=cmap[dt], alpha=0.7,
                     label=f'{dt}  (n={len(df[df["defect_type"]==dt]):,})')
               for dt in defect_types]
    handles.append(Line2D([0], [0], color='#333333', ls='--', lw=1.2,
                          label=f'Umbral = {threshold}'))
    ax0.legend(handles=handles, loc='upper right', handlelength=1.4)
    ax0.set_xlabel(r'$d_\mathrm{ref}$ (distancia SOAP a la referencia)')
    ax0.set_ylabel('Densidad')
    ax0.set_title('Distribución de distancia SOAP')
    ax0.set_xlim(0, xmax)

    # ── [0,1] Proyección XY coloreada por dist_to_ref ─────────────────────────
    ax1 = fig.add_subplot(gs[0, 1])
    vmax = sample50k['dist_to_ref'].quantile(0.99)
    sc = ax1.scatter(sample50k['x'], sample50k['y'],
                     c=sample50k['dist_to_ref'], cmap='plasma',
                     s=0.8, alpha=0.6, vmin=0, vmax=vmax, rasterized=True)
    add_colorbar(fig, ax1, sc, label=r'$d_\mathrm{ref}$')
    ax1.set_xlabel('x (Å)')
    ax1.set_ylabel('y (Å)')
    ax1.set_title('Proyección XY — distancia SOAP')
    ax1.set_aspect('equal', adjustable='box')

    # ── [1,0] Violin por tipo de defecto ──────────────────────────────────────
    ax2 = fig.add_subplot(gs[1, 0])
    groups = [df[df['defect_type'] == dt]['dist_to_ref'].values
              for dt in defect_types]
    # Usar boxplot cuando hay muy pocos puntos en algún grupo
    use_violin = all(len(g) >= 4 for g in groups)

    if use_violin:
        parts = ax2.violinplot(groups, positions=range(len(defect_types)),
                               showmedians=True, showextrema=False,
                               widths=0.75)
        for pc, dt in zip(parts['bodies'], defect_types):
            pc.set_facecolor(cmap[dt])
            pc.set_alpha(0.70)
        parts['cmedians'].set_color('black')
        parts['cmedians'].set_linewidth(1.5)
    else:
        bp = ax2.boxplot(groups, positions=range(len(defect_types)),
                         patch_artist=True, widths=0.55,
                         medianprops=dict(color='black', lw=1.5),
                         whiskerprops=dict(lw=0.8),
                         capprops=dict(lw=0.8),
                         flierprops=dict(marker='.', ms=2, alpha=0.4))
        for patch, dt in zip(bp['boxes'], defect_types):
            patch.set_facecolor(cmap[dt])
            patch.set_alpha(0.70)

    ax2.axhline(threshold, color='#333333', ls='--', lw=1.2,
                label=f'Umbral = {threshold}')
    ax2.set_xticks(range(len(defect_types)))
    ax2.set_xticklabels(defect_types, rotation=18, ha='right')
    ax2.set_ylabel(r'$d_\mathrm{ref}$')
    ax2.set_title('Distribución por tipo de defecto')
    ax2.legend(loc='upper left', handlelength=1.2)

    # ── [1,1] Proyección XZ coloreada por tipo de defecto ─────────────────────
    ax3 = fig.add_subplot(gs[1, 1])
    for dt in defect_types:
        sub = sample50k[sample50k['defect_type'] == dt]
        is_lattice = dt == 'Lattice'
        ax3.scatter(sub['x'], sub['z'],
                    c=cmap[dt], s=0.8 if is_lattice else 5,
                    alpha=0.35 if is_lattice else 0.85,
                    zorder=1 if is_lattice else 3,
                    rasterized=True)
    handles3 = [Line2D([0], [0], marker='o', color='w',
                       markerfacecolor=cmap[dt], markersize=5, label=dt)
                for dt in defect_types]
    ax3.legend(handles=handles3, loc='best', handlelength=1.0)
    ax3.set_xlabel('x (Å)')
    ax3.set_ylabel('z (Å)')
    ax3.set_title('Proyección XZ — tipo de defecto')
    ax3.set_aspect('equal', adjustable='box')

    # ── [2,0] Barras: conteo por tipo ─────────────────────────────────────────
    ax4 = fig.add_subplot(gs[2, 0])
    counts = df['defect_type'].value_counts().reindex(defect_types).fillna(0).astype(int)
    x_pos = np.arange(len(defect_types))
    bars = ax4.bar(x_pos, counts.values,
                   color=[cmap[dt] for dt in defect_types],
                   edgecolor='white', linewidth=0.5, zorder=3)
    for bar, val in zip(bars, counts.values):
        ax4.text(bar.get_x() + bar.get_width() / 2,
                 bar.get_height() * 1.015,
                 f'{val:,}', ha='center', va='bottom', fontsize=7.5)
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels(defect_types, rotation=18, ha='right')
    ax4.set_ylabel('Número de átomos')
    ax4.set_title('Recuento por clasificación')
    pct = '  |  '.join(
        f'{dt}: {100 * v / len(df):.1f}%'
        for dt, v in zip(defect_types, counts.values)
    )
    ax4.text(0.5, -0.30, pct, transform=ax4.transAxes,
             ha='center', va='top', fontsize=7, style='italic', color='#555555')

    # ── [2,1] Distribución de defect_prob ─────────────────────────────────────
    ax5 = fig.add_subplot(gs[2, 1])
    # Separar los que tienen prob < 1 (interesante) del resto
    non_one = df[df['defect_prob'] < 0.9999]['defect_prob']
    if len(non_one) >= 5:
        ax5.hist(non_one, bins=60, color='#4393c3', alpha=0.75,
                 density=True, histtype='stepfilled', edgecolor='none')
        x_prob = np.linspace(non_one.min(), non_one.max(), 300)
        kde_p = kde_curve(non_one.values, x_prob)
        if kde_p is not None:
            ax5.plot(x_prob, kde_p, color='#08519c', lw=1.4)
        subtitle = f'(excluyendo prob = 1.0,  n = {len(non_one):,})'
    else:
        ax5.hist(df['defect_prob'], bins=60, color='#4393c3', alpha=0.75,
                 density=True, histtype='stepfilled', edgecolor='none')
        subtitle = f'(n = {len(df):,})'
    ax5.set_xlabel('Probabilidad de defecto')
    ax5.set_ylabel('Densidad')
    ax5.set_title(f'Distribución de defect_prob\n{subtitle}')

    # ── Título global ──────────────────────────────────────────────────────────
    ntot = len(df)
    ndef = len(df[df['defect_type'] != 'Lattice'])
    fig.suptitle(
        f'{label} — Análisis SOAP\n'
        f'N = {ntot:,} átomos  |  '
        f'Defectos: {ndef:,} ({100 * ndef / ntot:.2f} %)  |  '
        f'Umbral = {threshold}',
        fontsize=11, y=1.01
    )

    save_figure(fig, f'{stem}_overview', outdir)
    plt.close(fig)


# ─── Figura 2: pca_output*.csv ────────────────────────────────────────────────

def plot_pca(df: pd.DataFrame, label: str, threshold: float,
             outdir: Path, stem: str):
    """
    3-panel figure para pca_output*.csv:
      [0] PC1 vs PC2 coloreado por tipo de defecto
      [1] PC1 vs PC2 coloreado por dist_to_ref
      [2] Histograma + KDE de dist_to_ref
    """
    set_publication_style()

    defect_types = sorted(df['defect_type'].unique())
    cmap = get_color_map(defect_types)
    sample80k = _sample(df, 80_000)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6))
    fig.subplots_adjust(wspace=0.40)

    # ── Panel 0: PCA coloreado por tipo ───────────────────────────────────────
    ax0 = axes[0]
    for dt in defect_types:
        sub = sample80k[sample80k['defect_type'] == dt]
        is_lattice = dt == 'Lattice'
        ax0.scatter(sub['pc1'], sub['pc2'],
                    c=cmap[dt],
                    s=1.5 if is_lattice else 8,
                    alpha=0.35 if is_lattice else 0.85,
                    zorder=1 if is_lattice else 3,
                    rasterized=True)
    handles = [Line2D([0], [0], marker='o', color='w',
                      markerfacecolor=cmap[dt], markersize=5,
                      label=f'{dt}  (n={len(df[df["defect_type"]==dt]):,})')
               for dt in defect_types]
    ax0.legend(handles=handles, loc='best', handlelength=1.0)
    ax0.set_xlabel(r'PC$_1$')
    ax0.set_ylabel(r'PC$_2$')
    ax0.set_title('Espacio PCA — tipo de defecto')

    # ── Panel 1: PCA coloreado por dist_to_ref ────────────────────────────────
    ax1 = axes[1]
    vmax = sample80k['dist_to_ref'].quantile(0.99)
    sc = ax1.scatter(sample80k['pc1'], sample80k['pc2'],
                     c=sample80k['dist_to_ref'],
                     cmap='plasma', s=1.0, alpha=0.6,
                     vmin=0, vmax=vmax, rasterized=True)
    add_colorbar(fig, ax1, sc, label=r'$d_\mathrm{ref}$')
    ax1.set_xlabel(r'PC$_1$')
    ax1.set_ylabel(r'PC$_2$')
    ax1.set_title('Espacio PCA — distancia SOAP')

    # ── Panel 2: Histograma + KDE dist_to_ref ─────────────────────────────────
    ax2 = axes[2]
    xmax = df['dist_to_ref'].quantile(0.999)
    bins = np.linspace(0, xmax, 80)
    x_grid = np.linspace(0, xmax, 400)
    for dt in defect_types:
        vals = df[df['defect_type'] == dt]['dist_to_ref'].values
        ax2.hist(vals, bins=bins, alpha=0.55, color=cmap[dt],
                 density=True, histtype='stepfilled', edgecolor='none')
        kde_y = kde_curve(vals, x_grid)
        if kde_y is not None:
            ax2.plot(x_grid, kde_y, color=cmap[dt], lw=1.4)
    ax2.axvline(threshold, color='#333333', ls='--', lw=1.2,
                label=f'Umbral = {threshold}')
    handles2 = [Patch(facecolor=cmap[dt], alpha=0.7, label=dt)
                for dt in defect_types]
    handles2.append(Line2D([0], [0], color='#333333', ls='--', lw=1.2,
                            label=f'Umbral = {threshold}'))
    ax2.legend(handles=handles2, loc='upper right', handlelength=1.4)
    ax2.set_xlabel(r'$d_\mathrm{ref}$')
    ax2.set_ylabel('Densidad')
    ax2.set_title('Distribución de distancia SOAP')
    ax2.set_xlim(0, xmax)

    fig.suptitle(
        f'{label} — Análisis PCA\n'
        f'N = {len(df):,} átomos  |  Umbral = {threshold}',
        fontsize=11, y=1.04
    )

    save_figure(fig, f'{stem}_pca', outdir)
    plt.close(fig)


# ─── Figura 3: hist_output*.csv ───────────────────────────────────────────────

def plot_hist(df: pd.DataFrame, label: str, threshold: float,
              outdir: Path, stem: str):
    """
    2-panel figure para hist_output*.csv:
      [0] Escala lineal con área rellena + KDE suavizada
      [1] Escala semilogarítmica
    """
    set_publication_style()

    x = df['bin_centre'].values
    y = df['count'].values.astype(float)
    total = int(y.sum())
    mean_val = float((x * y).sum() / total) if total > 0 else 0.0
    n_above = float(y[x > threshold].sum())

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9.5, 4.2))
    fig.subplots_adjust(wspace=0.38)

    # ── Panel 0: escala lineal ─────────────────────────────────────────────────
    ax0.fill_between(x, 0, y, alpha=0.50, color='#2166ac', step='mid')
    ax0.step(x, y, color='#2166ac', linewidth=1.1, where='mid')

    # Curva KDE suavizada (usando los datos del histograma como weights)
    if HAS_SCIPY:
        # Re-muestrear puntos con pesos proporcionales a los conteos
        rng = np.random.default_rng(42)
        n_samples = min(total, 50_000)
        idx = rng.choice(len(x), size=n_samples, p=y / y.sum())
        samples = x[idx] + rng.uniform(-0.5 * (x[1] - x[0]),
                                        0.5 * (x[1] - x[0]), n_samples)
        x_kde = np.linspace(x.min(), x.max(), 500)
        kde_y = kde_curve(samples, x_kde)
        if kde_y is not None:
            kde_y_scaled = kde_y * total * (x[1] - x[0])
            ax0.plot(x_kde, kde_y_scaled, color='#d73027', lw=1.4,
                     label='KDE suavizada')

    ax0.axvline(threshold, color='#333333', ls='--', lw=1.3,
                label=f'Umbral = {threshold}')
    ax0.axvline(mean_val, color='#1a9641', ls=':', lw=1.2,
                label=f'Media = {mean_val:.4f}')

    stats_text = (f'N total = {total:,}\n'
                  f'Media = {mean_val:.4f}\n'
                  f'N > {threshold} = {int(n_above):,} ({100 * n_above / total:.1f} %)')
    ax0.text(0.97, 0.97, stats_text, transform=ax0.transAxes,
             ha='right', va='top', fontsize=7.5,
             bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                       edgecolor='#cccccc', alpha=0.92))
    ax0.set_xlabel(r'$d_\mathrm{ref}$ (distancia SOAP a la referencia)')
    ax0.set_ylabel('Número de átomos')
    ax0.set_title('Distribución de distancias (escala lineal)')
    ax0.set_xlim(0, x.max())
    ax0.legend(handlelength=1.4, loc='upper right')

    # ── Panel 1: escala semilogarítmica ───────────────────────────────────────
    y_log = np.where(y > 0, y, np.nan)
    ax1.fill_between(x, 1, y_log, alpha=0.50, color='#4d9221', step='mid')
    ax1.step(x, y_log, color='#4d9221', linewidth=1.1, where='mid')
    ax1.axvline(threshold, color='#333333', ls='--', lw=1.3,
                label=f'Umbral = {threshold}')
    ax1.axvline(mean_val, color='#1a9641', ls=':', lw=1.2,
                label=f'Media = {mean_val:.4f}')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$d_\mathrm{ref}$ (distancia SOAP a la referencia)')
    ax1.set_ylabel('Número de átomos (log)')
    ax1.set_title('Distribución de distancias (escala semilog)')
    ax1.set_xlim(0, x.max())
    ax1.legend(handlelength=1.4, loc='upper right')

    fig.suptitle(f'{label} — Histograma de distancias SOAP',
                 fontsize=11, y=1.04)

    save_figure(fig, f'{stem}_hist', outdir)
    plt.close(fig)


# ─── Figura 4: vacancies_output*.csv ─────────────────────────────────────────

def plot_vacancies(df: pd.DataFrame, label: str, outdir: Path, stem: str):
    """
    4-panel figure para vacancies_output*.csv (post-clustering):
      [0,0] Proyección XY  — color = d_near_max, tamaño ∝ n_grid_pts
      [0,1] Proyección XZ  — ídem
      [1,0] Proyección YZ  — ídem
      [1,1] Distribución de vacancies: histograma de n_grid_pts por cluster
    """
    set_publication_style()

    x, y, z = df['x'].values, df['y'].values, df['z'].values
    Lx = x.max() - x.min()
    Ly = y.max() - y.min()
    Lz = z.max() - z.min()

    # Columnas opcionales (archivos viejos sin clustering no las tienen)
    has_cluster_cols = ('d_near_max' in df.columns and 'n_grid_pts' in df.columns
                        and df['d_near_max'].notna().any())

    if has_cluster_cols:
        d_near = df['d_near_max'].values
        n_pts  = df['n_grid_pts'].values
        # Tamaño de marcador: escalar entre 8 y 80 según raíz de n_grid_pts
        s_raw  = np.sqrt(n_pts.astype(float))
        s_min, s_max = s_raw.min(), s_raw.max()
        s_norm = (s_raw - s_min) / (s_max - s_min + 1e-14)
        marker_size = 8 + 72 * s_norm
        color_vals  = d_near
        cbar_label  = 'd_near_max (Å)'
        cmap        = 'hot_r'
    else:
        # Fallback para archivos sin columnas de clustering
        def norm01(v):
            vmin, vmax = v.min(), v.max()
            return (v - vmin) / (vmax - vmin + 1e-14)
        color_vals  = norm01(z)
        marker_size = 4
        cbar_label  = 'z (normalizado)'
        cmap        = 'viridis'
        n_pts       = None

    fig = plt.figure(figsize=(11, 9))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.44, wspace=0.38)

    scatter_kw = dict(alpha=0.75, rasterized=True)

    # ── XY ────────────────────────────────────────────────────────────────────
    ax0 = fig.add_subplot(gs[0, 0])
    sc0 = ax0.scatter(x, y, c=color_vals, s=marker_size, cmap=cmap, **scatter_kw)
    add_colorbar(fig, ax0, sc0, label=cbar_label)
    ax0.set_xlabel('x (Å)')
    ax0.set_ylabel('y (Å)')
    ax0.set_title('Proyección XY')
    ax0.set_aspect('equal', adjustable='box')

    # ── XZ ────────────────────────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, 1])
    sc1 = ax1.scatter(x, z, c=color_vals, s=marker_size, cmap=cmap, **scatter_kw)
    add_colorbar(fig, ax1, sc1, label=cbar_label)
    ax1.set_xlabel('x (Å)')
    ax1.set_ylabel('z (Å)')
    ax1.set_title('Proyección XZ')
    ax1.set_aspect('equal', adjustable='box')

    # ── YZ ────────────────────────────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[1, 0])
    sc2 = ax2.scatter(y, z, c=color_vals, s=marker_size, cmap=cmap, **scatter_kw)
    add_colorbar(fig, ax2, sc2, label=cbar_label)
    ax2.set_xlabel('y (Å)')
    ax2.set_ylabel('z (Å)')
    ax2.set_title('Proyección YZ')
    ax2.set_aspect('equal', adjustable='box')

    # ── Distribución de tamaño de clusters / vacancies ────────────────────────
    ax3 = fig.add_subplot(gs[1, 1])
    if has_cluster_cols and n_pts is not None:
        # Bins 1–100 individuales; todo lo que supere 100 se acumula en el
        # último bin (posición 101 en el eje, etiquetado "100+").
        OVF = 101   # posición del bin de overflow en el eje x
        n_clipped = np.where(n_pts > 100, OVF, n_pts)
        bins = np.arange(0.5, OVF + 1.5, 1)   # [0.5, 1.5, …, 101.5] → 101 bins

        ax3.hist(n_clipped, bins=bins, color='#2166ac', edgecolor='white',
                 linewidth=0.3, rasterized=True)

        # Etiquetas del eje x: cada 10 unidades + "100+" al final
        tick_pos    = list(range(10, 101, 10)) + [OVF]
        tick_labels = [str(t) for t in range(10, 101, 10)] + ['100+']
        ax3.set_xticks(tick_pos)
        ax3.set_xticklabels(tick_labels, fontsize=7)
        ax3.set_xlim(0.5, OVF + 0.5)

        # Separador visual entre la escala regular y el bin de overflow
        ax3.axvline(OVF - 0.5, color='#888888', lw=0.8, ls=':')

        # Anotar cuántas vacancias cayeron en overflow si las hay
        overflow_n = int((n_pts > 100).sum())
        if overflow_n > 0:
            ymax = ax3.get_ylim()[1]
            ax3.text(OVF, overflow_n * 1.02, str(overflow_n),
                     ha='center', va='bottom', fontsize=7, color='#d73027')

        # Línea de mediana (solo si está en el rango regular)
        med = float(np.median(n_pts))
        if med <= 100:
            ax3.axvline(med, color='#d73027', lw=1.2, ls='--',
                        label=f'mediana = {med:.0f}')
            ax3.legend(fontsize=7.5)

        ax3.set_xlabel('Puntos de grilla por vacancia (n_grid_pts)')
        ax3.set_ylabel('Número de vacancias')
        ax3.set_title('Distribución de vacancias por tamaño')
    else:
        # Fallback: distribuciones 1D por coordenada
        nbins = 60
        palette = {'x': '#1a9641', 'y': '#2166ac', 'z': '#d73027'}
        for coord_name, vals in [('x', x), ('y', y), ('z', z)]:
            counts, edges = np.histogram(vals, bins=nbins)
            centers = 0.5 * (edges[:-1] + edges[1:])
            norm_counts = counts / counts.max()
            ax3.fill_between(centers, 0, norm_counts,
                             alpha=0.25, color=palette[coord_name], step='mid')
            ax3.step(centers, norm_counts, color=palette[coord_name],
                     lw=1.2, where='mid',
                     label=f'{coord_name}  ({vals.min():.1f}–{vals.max():.1f} Å)')
        ax3.set_xlabel('Coordenada (Å)')
        ax3.set_ylabel('Densidad normalizada')
        ax3.set_title('Distribución espacial de vacancias')
        ax3.legend(fontsize=7.5, handlelength=1.4)

    vol = Lx * Ly * Lz
    density = len(df) / vol if vol > 0 else 0.0
    fig.suptitle(
        f'{label} — Vacancias detectadas  (N = {len(df):,})\n'
        f'Caja: {Lx:.1f} × {Ly:.1f} × {Lz:.1f} Å'
        r'$\quad$'
        f'Densidad: {density:.4f} Å$^{{-3}}$',
        fontsize=11, y=1.02
    )

    save_figure(fig, f'{stem}_vacancies', outdir)
    plt.close(fig)


# ─── Dispatcher y CLI ─────────────────────────────────────────────────────────

def process_file(csv_path: Path, threshold: float, label: str, outdir: Path):
    csv_type = detect_csv_type(csv_path)
    lbl = label or csv_path.stem.replace('_', ' ').title()
    print(f"  [{csv_path.name}]  tipo='{csv_type}'  label='{lbl}'")

    if csv_type == 'output':
        df = load_output_csv(csv_path)
        print(f"    {len(df):,} átomos  |  tipos: "
              f"{dict(df['defect_type'].value_counts())}")
        plot_output(df, lbl, threshold, outdir, csv_path.stem)

    elif csv_type == 'pca':
        df = load_pca_csv(csv_path)
        print(f"    {len(df):,} átomos  |  tipos: "
              f"{dict(df['defect_type'].value_counts())}")
        plot_pca(df, lbl, threshold, outdir, csv_path.stem)

    elif csv_type == 'hist':
        df = load_hist_csv(csv_path)
        print(f"    {len(df):,} bins  |  rango: "
              f"{df['bin_centre'].min():.4f} – {df['bin_centre'].max():.4f}")
        plot_hist(df, lbl, threshold, outdir, csv_path.stem)

    elif csv_type == 'vacancies':
        df = load_vacancies_csv(csv_path)
        print(f"    {len(df):,} vacancias")
        plot_vacancies(df, lbl, outdir, csv_path.stem)

    else:
        print(f"    AVISO: tipo no reconocido, se omite.")


def main():
    parser = argparse.ArgumentParser(
        description='Gráficos de calidad de publicación para DistributionTool.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('files', nargs='*', metavar='CSV',
                        help='Archivo(s) CSV a graficar')
    parser.add_argument('--all', action='store_true',
                        help='Procesa todos los *.csv del directorio actual')
    parser.add_argument('--outdir', default='figures',
                        help='Directorio de salida  (default: figures/)')
    parser.add_argument('--threshold', type=float, default=0.15,
                        help='Umbral dist_to_ref para marcar en las figuras  (default: 0.15)')
    parser.add_argument('--label', default='',
                        help='Etiqueta descriptiva para los títulos (opcional)')

    args = parser.parse_args()

    if args.all:
        csv_files = sorted(Path('.').glob('*.csv'))
    elif args.files:
        csv_files = [Path(f) for f in args.files]
    else:
        parser.print_help()
        sys.exit(0)

    if not csv_files:
        print('No se encontraron archivos CSV.')
        sys.exit(1)

    outdir = Path(args.outdir)
    print(f'\nProcesando {len(csv_files)} archivo(s) → "{outdir}/"\n')

    for p in csv_files:
        if not p.exists():
            print(f'  AVISO: "{p}" no existe, se omite.')
            continue
        try:
            process_file(p, args.threshold, args.label, outdir)
        except Exception as exc:
            print(f'  ERROR en "{p}": {exc}')
            raise

    print(f'\nListo. Figuras guardadas en "{outdir}/"')


if __name__ == '__main__':
    main()
