#!/usr/bin/env python3
"""
publication_hybrid.py — Figuras de calidad de publicación para el método híbrido.

Toma como entrada los CSV producidos por distool en modo --hybrid --ws y genera
figuras dedicadas al análisis del super-método híbrido:

  Fig.1  WS vs Hybrid: barras + reducción porcentual + breakdown de filtros
  Fig.2  Histograma 2D de consensus_score por estado (accepted/transit/frenkel)
  Fig.3  Stacked breakdown de las 6 señales positivas por candidato
  Fig.4  Mapa espacial 3D de candidatos coloreados por estado
  Fig.5  Heatmap de correlación entre señales
  Fig.6  Histograma de las penalties (transit, frenkel) y su correlación
  Fig.7  Cuasi-Frenkel: distribución de distancias vac↔int más cercana

Uso:
  python publication_hybrid.py <hybrid_vacancies_csv> <ws_sites_csv> [--outdir DIR] [--label NAME]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)


# ─── Estilo publicación ──────────────────────────────────────────────────────
def set_publication_style():
    mpl.rcParams.update({
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'figure.facecolor': 'white',
        'font.family': 'serif',
        'font.serif': ['DejaVu Serif', 'Times New Roman', 'Times', 'serif'],
        'font.size': 9,
        'axes.titlesize': 10,
        'axes.labelsize': 9,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'legend.title_fontsize': 8,
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': True,
        'grid.alpha': 0.22,
        'grid.linewidth': 0.5,
        'grid.color': '#888888',
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.major.size': 3.5,
        'ytick.major.size': 3.5,
        'lines.linewidth': 1.2,
    })


# Paleta consistente
COLOR_WS_ONLY  = '#d62728'   # rojo
COLOR_ACCEPT   = '#2ca02c'   # verde
COLOR_HYB_ONLY = '#1f77b4'   # azul
COLOR_TRANSIT  = '#ff7f0e'   # naranja
COLOR_FRENKEL  = '#9467bd'   # violeta
COLOR_NEUTRAL  = '#7f7f7f'   # gris

SIGNAL_LABELS = {
    'ws':       'WS',
    'soft_ws':  'Soft-WS',
    'vor':      'Voronoi',
    'dens':     'Densidad',
    'soap':     'SOAP',
    'topo':     'Topología',
}
SIGNAL_COLORS = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02']


# ─── Loaders ─────────────────────────────────────────────────────────────────
def load_hybrid_csv(path: Path) -> pd.DataFrame:
    cols = [
        'x', 'y', 'z', 'consensus_score',
        'ws', 'soft_ws', 'vor', 'dens', 'soap', 'topo',
        'transit_pen', 'frenkel_pen',
        'transit_filtered', 'in_recomb', 'ref_site_idx', 'accepted',
    ]
    df = pd.read_csv(path, sep=r'\s+', comment='#', header=None, names=cols)
    df['status'] = np.where(
        df['accepted'] == 1, 'accepted',
        np.where(df['transit_filtered'] == 1, 'transit',
                 np.where(df['in_recomb'] == 1, 'frenkel', 'low_score'))
    )
    return df


def load_ws_sites_csv(path: Path) -> pd.DataFrame:
    cols = ['ref_id', 'x', 'y', 'z', 'occupancy', 'min_dist', 'site_type']
    df = pd.read_csv(path, sep=r'\s+', comment='#', header=None, names=cols)
    return df


# ─── Save helper ─────────────────────────────────────────────────────────────
def save(fig, name: str, outdir: Path, formats=('pdf', 'png')):
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(outdir / f'{name}.{ext}')
    plt.close(fig)
    print(f'  → {outdir / name}.{{{",".join(formats)}}}')


# ─── Fig.1: WS vs Hybrid summary ─────────────────────────────────────────────
def fig_ws_vs_hybrid(df: pd.DataFrame, label: str, outdir: Path):
    n_total       = len(df)
    n_accepted    = int((df['accepted'] == 1).sum())
    n_transit     = int(df['transit_filtered'].sum())
    n_frenkel     = int(((df['in_recomb'] == 1) & (df['accepted'] == 0) &
                          (df['transit_filtered'] == 0)).sum())
    n_low_score   = n_total - n_accepted - n_transit - n_frenkel
    n_ws_vac      = int((df['ref_site_idx'] >= 0).sum())
    n_hyb_only    = int(((df['accepted'] == 1) & (df['ref_site_idx'] < 0)).sum())
    reduction_pct = 100.0 * (1.0 - n_accepted / max(n_ws_vac, 1))

    fig = plt.figure(figsize=(7.2, 3.0))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.2], wspace=0.35)

    # Panel A: barras WS vs Hybrid
    ax1 = fig.add_subplot(gs[0])
    bars = ax1.bar(
        ['Wigner–Seitz', 'Hybrid\n(robust)'],
        [n_ws_vac, n_accepted],
        color=[COLOR_WS_ONLY, COLOR_ACCEPT],
        width=0.55,
        edgecolor='black', linewidth=0.6,
    )
    for b, v in zip(bars, [n_ws_vac, n_accepted]):
        ax1.text(b.get_x() + b.get_width() / 2, v + 0.02 * n_ws_vac,
                 f'{v:,}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax1.set_ylabel('Vacancias detectadas')
    ax1.set_title(f'(a) Conteo total\n$\\Delta = -{reduction_pct:.1f}\\%$',
                  fontsize=9.5)
    ax1.set_ylim(0, n_ws_vac * 1.15)
    ax1.grid(axis='y', alpha=0.3)

    # Panel B: stacked breakdown del filtrado
    ax2 = fig.add_subplot(gs[1])
    categories = ['Aceptadas', 'Filtradas\n(tránsito)',
                  'Filtradas\n(Frenkel)', 'Filtradas\n(score bajo)']
    counts     = [n_accepted, n_transit, n_frenkel, n_low_score]
    colors     = [COLOR_ACCEPT, COLOR_TRANSIT, COLOR_FRENKEL, COLOR_NEUTRAL]
    bars2 = ax2.bar(categories, counts, color=colors, edgecolor='black', linewidth=0.6)
    for b, v in zip(bars2, counts):
        pct = 100.0 * v / max(n_total, 1)
        ax2.text(b.get_x() + b.get_width() / 2, v + 0.015 * max(counts),
                 f'{v:,}\n({pct:.1f}%)', ha='center', va='bottom', fontsize=8)
    ax2.set_ylabel('Candidatos')
    ax2.set_title(f'(b) Destino de los {n_total:,} candidatos',
                  fontsize=9.5)
    ax2.set_ylim(0, max(counts) * 1.25)
    ax2.grid(axis='y', alpha=0.3)

    fig.suptitle(f'Hybrid vs Wigner–Seitz — {label}', fontsize=10.5, y=1.04)
    save(fig, 'fig1_ws_vs_hybrid', outdir)


# ─── Fig.2: consensus score distribution by status ──────────────────────────
def fig_score_distribution(df: pd.DataFrame, label: str, outdir: Path):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))

    # Panel A: histograma apilado por estado
    ax = axes[0]
    bins = np.linspace(df['consensus_score'].min() - 0.05,
                       df['consensus_score'].max() + 0.05, 50)

    for status, color, name in [
        ('accepted', COLOR_ACCEPT,   'Aceptado'),
        ('transit',  COLOR_TRANSIT,  'Filtrado tránsito'),
        ('frenkel',  COLOR_FRENKEL,  'Filtrado Frenkel'),
        ('low_score',COLOR_NEUTRAL,  'Score bajo'),
    ]:
        sub = df[df['status'] == status]['consensus_score']
        if len(sub):
            ax.hist(sub, bins=bins, alpha=0.75, label=f'{name} ({len(sub):,})',
                    color=color, edgecolor='black', linewidth=0.3)
    ax.axvline(0.5, color='black', ls='--', lw=0.8, alpha=0.6,
               label='Umbral $\\tau = 0.5$')
    ax.set_xlabel('Score de consenso $S(\\mathbf{x})$')
    ax.set_ylabel('Candidatos')
    ax.set_title('(a) Distribución de scores por estado')
    ax.legend(loc='upper left', framealpha=0.9, fontsize=7)

    # Panel B: scatter score vs penalty
    ax = axes[1]
    statuses = [('accepted', COLOR_ACCEPT, 'o'),
                ('transit',  COLOR_TRANSIT, 's'),
                ('frenkel',  COLOR_FRENKEL, '^'),
                ('low_score',COLOR_NEUTRAL, 'd')]
    max_pen = df[['transit_pen', 'frenkel_pen']].max(axis=1)
    for status, color, marker in statuses:
        sub = df[df['status'] == status]
        if len(sub):
            ax.scatter(sub['consensus_score'], max_pen[sub.index],
                       s=8, c=color, marker=marker, alpha=0.6,
                       edgecolor='black', linewidth=0.15)
    ax.axhline(0.5, color='black', ls=':', lw=0.6, alpha=0.6)
    ax.axvline(0.5, color='black', ls='--', lw=0.6, alpha=0.6)
    ax.set_xlabel('Score de consenso $S(\\mathbf{x})$')
    ax.set_ylabel('$\\max(p_\\text{transit},\\ p_\\text{Frenkel})$')
    ax.set_title('(b) Score vs. penalty máxima')

    fig.suptitle(f'Régimen del consenso — {label}', fontsize=10.5, y=1.04)
    fig.tight_layout()
    save(fig, 'fig2_score_distribution', outdir)


# ─── Fig.3: stacked signal contributions ────────────────────────────────────
def fig_signal_breakdown(df: pd.DataFrame, label: str, outdir: Path):
    # Pesos del preset robust (deben coincidir con HybridParams::applyPreset)
    weights = {'ws': 1.0, 'soft_ws': 0.8, 'vor': 0.6,
               'dens': 0.6, 'soap': 0.5, 'topo': 0.7}
    w_sum   = sum(weights.values())
    sig_cols = list(weights.keys())

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))

    # Panel A: distribución de cada señal (boxplot)
    ax = axes[0]
    data = [df[col].values for col in sig_cols]
    bp = ax.boxplot(data, patch_artist=True, widths=0.6,
                    medianprops=dict(color='black', lw=1.0),
                    flierprops=dict(marker='.', markersize=1.5, alpha=0.4))
    for patch, color in zip(bp['boxes'], SIGNAL_COLORS):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
        patch.set_edgecolor('black')
        patch.set_linewidth(0.6)
    ax.set_xticks(range(1, len(sig_cols) + 1))
    ax.set_xticklabels([SIGNAL_LABELS[c] for c in sig_cols], rotation=20, ha='right')
    ax.set_ylabel('Valor de la señal $s_\\alpha$')
    ax.set_title('(a) Distribución de cada señal')
    ax.set_ylim(-0.05, 1.05)

    # Panel B: contribución promedio al consenso (w·s / Σw)
    ax = axes[1]
    n_accepted = (df['accepted'] == 1).sum()
    sub_acc = df[df['accepted'] == 1]
    sub_fil = df[df['accepted'] == 0]
    means_acc = [(weights[c] * sub_acc[c].mean() / w_sum) if len(sub_acc) else 0
                 for c in sig_cols]
    means_fil = [(weights[c] * sub_fil[c].mean() / w_sum) if len(sub_fil) else 0
                 for c in sig_cols]
    x = np.arange(len(sig_cols))
    width = 0.4
    ax.bar(x - width / 2, means_acc, width, label=f'Aceptadas ({n_accepted:,})',
           color=COLOR_ACCEPT, edgecolor='black', linewidth=0.5)
    ax.bar(x + width / 2, means_fil, width,
           label=f'Filtradas ({len(sub_fil):,})',
           color=COLOR_WS_ONLY, edgecolor='black', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([SIGNAL_LABELS[c] for c in sig_cols], rotation=20, ha='right')
    ax.set_ylabel('Contribución $w_\\alpha s_\\alpha / \\Sigma w$')
    ax.set_title('(b) Contribución promedio por señal')
    ax.legend(loc='upper right', fontsize=7)

    fig.suptitle(f'Desglose de las 6 señales positivas — {label}', fontsize=10.5, y=1.04)
    fig.tight_layout()
    save(fig, 'fig3_signal_breakdown', outdir)


# ─── Fig.4: 3D spatial map ──────────────────────────────────────────────────
def fig_spatial_3d(df: pd.DataFrame, label: str, outdir: Path):
    fig = plt.figure(figsize=(7.2, 6.4))
    ax = fig.add_subplot(111, projection='3d')

    for status, color, name, alpha, size in [
        ('low_score', COLOR_NEUTRAL,  'Score bajo',    0.20, 4),
        ('frenkel',   COLOR_FRENKEL,  'Frenkel coh.',  0.55, 8),
        ('transit',   COLOR_TRANSIT,  'Tránsito',      0.55, 8),
        ('accepted',  COLOR_ACCEPT,   'Aceptado',      0.80, 14),
    ]:
        sub = df[df['status'] == status]
        if len(sub):
            ax.scatter(sub['x'], sub['y'], sub['z'],
                       s=size, c=color, alpha=alpha,
                       edgecolor='black' if status == 'accepted' else None,
                       linewidth=0.15,
                       label=f'{name} ({len(sub):,})')

    ax.set_xlabel('x [Å]')
    ax.set_ylabel('y [Å]')
    ax.set_zlabel('z [Å]')
    ax.set_title(f'Distribución espacial de candidatos — {label}', pad=10)
    ax.legend(loc='upper left', bbox_to_anchor=(0.02, 0.98), fontsize=8)
    ax.view_init(elev=22, azim=-65)
    save(fig, 'fig4_spatial_3d', outdir)


# ─── Fig.5: signal correlation heatmap ──────────────────────────────────────
def fig_signal_correlation(df: pd.DataFrame, label: str, outdir: Path):
    sig_cols = ['ws', 'soft_ws', 'vor', 'dens', 'soap', 'topo',
                'transit_pen', 'frenkel_pen']
    pretty = ['WS', 'Soft-WS', 'Voronoi', 'Densidad', 'SOAP', 'Topología',
              '$p_\\text{transit}$', '$p_\\text{Frenkel}$']
    corr = df[sig_cols].corr()

    fig, ax = plt.subplots(figsize=(5.5, 5.0))
    im = ax.imshow(corr.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='equal')
    ax.set_xticks(range(len(sig_cols)))
    ax.set_yticks(range(len(sig_cols)))
    ax.set_xticklabels(pretty, rotation=35, ha='right')
    ax.set_yticklabels(pretty)
    for i in range(len(sig_cols)):
        for j in range(len(sig_cols)):
            v = corr.values[i, j]
            ax.text(j, i, f'{v:.2f}', ha='center', va='center',
                    fontsize=7, color='white' if abs(v) > 0.55 else 'black')
    ax.grid(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Correlación de Pearson')
    ax.set_title(f'Estructura de correlación entre señales — {label}',
                 fontsize=10.5, pad=10)
    save(fig, 'fig5_signal_correlation', outdir)


# ─── Fig.6: penalties analysis ──────────────────────────────────────────────
def fig_penalties(df: pd.DataFrame, label: str, outdir: Path):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))

    # Panel A: distribuciones de las penalties
    ax = axes[0]
    bins = np.linspace(0, 1, 50)
    ax.hist(df['transit_pen'], bins=bins, alpha=0.65, color=COLOR_TRANSIT,
            edgecolor='black', linewidth=0.3,
            label=f'Tránsito $p_T$ (>0.5: {int((df["transit_pen"] > 0.5).sum()):,})')
    ax.hist(df['frenkel_pen'], bins=bins, alpha=0.65, color=COLOR_FRENKEL,
            edgecolor='black', linewidth=0.3,
            label=f'Frenkel $p_F$ (>0.5: {int((df["frenkel_pen"] > 0.5).sum()):,})')
    ax.axvline(0.5, color='black', ls='--', lw=0.8, alpha=0.6)
    ax.set_xlabel('Valor de la penalty')
    ax.set_ylabel('Candidatos')
    ax.set_title('(a) Distribuciones de penalties')
    ax.set_yscale('log')
    ax.legend(loc='upper right', fontsize=7)

    # Panel B: scatter transit vs frenkel
    ax = axes[1]
    ax.scatter(df['transit_pen'], df['frenkel_pen'],
               s=8, c=df['consensus_score'], cmap='viridis',
               alpha=0.7, edgecolor='black', linewidth=0.15)
    ax.axhline(0.5, color='black', ls=':', lw=0.5, alpha=0.6)
    ax.axvline(0.5, color='black', ls='--', lw=0.5, alpha=0.6)
    ax.set_xlabel('Penalty de tránsito $p_T$')
    ax.set_ylabel('Penalty de Frenkel $p_F$')
    ax.set_title('(b) Correlación tránsito–Frenkel')
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    cbar = plt.colorbar(ax.collections[0], ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Score $S$', fontsize=8)

    fig.suptitle(f'Análisis de las penalties — {label}', fontsize=10.5, y=1.04)
    fig.tight_layout()
    save(fig, 'fig6_penalties', outdir)


# ─── Fig.7: Frenkel pair distance distribution ──────────────────────────────
def fig_frenkel_distances(df_hyb: pd.DataFrame, df_ws_sites: pd.DataFrame,
                           label: str, outdir: Path,
                           recomb_radius: float = 3.3):
    vacs = df_ws_sites[df_ws_sites['occupancy'] == 0][['x', 'y', 'z']].values
    ints = df_ws_sites[df_ws_sites['occupancy'] >= 2][['x', 'y', 'z']].values
    if len(vacs) == 0 or len(ints) == 0:
        print('  [skip] no WS vacancies/interstitials')
        return

    # k-d tree style: small enough (a few thousand) so brute force is fine
    from scipy.spatial import cKDTree
    tree = cKDTree(ints)
    nearest, _ = tree.query(vacs, k=1)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))

    # Panel A: histograma de distancias mínimas vac↔int
    ax = axes[0]
    bins = np.linspace(0, np.percentile(nearest, 99) * 1.05, 60)
    ax.hist(nearest, bins=bins, color=COLOR_WS_ONLY, alpha=0.75,
            edgecolor='black', linewidth=0.3)
    ax.axvline(recomb_radius, color='red', ls='--', lw=1.0,
               label=f'$r_\\text{{rec}} = {recomb_radius:.1f}$ Å')
    ax.set_xlabel('Distancia vacancia–intersticial más cercana [Å]')
    ax.set_ylabel('Pares')
    ax.set_title('(a) Distribución de distancias Frenkel')
    ax.legend(loc='upper right', fontsize=8)

    # Panel B: fracción cumulativa de Frenkel cercanos
    ax = axes[1]
    sorted_d = np.sort(nearest)
    cumfrac = np.arange(1, len(sorted_d) + 1) / len(sorted_d)
    ax.plot(sorted_d, cumfrac * 100, lw=1.5, color=COLOR_WS_ONLY)
    ax.axvline(recomb_radius, color='red', ls='--', lw=1.0,
               label=f'$r_\\text{{rec}} = {recomb_radius:.1f}$ Å')
    frac_at_rec = 100 * (nearest < recomb_radius).mean()
    ax.axhline(frac_at_rec, color='gray', ls=':', lw=0.6, alpha=0.8)
    ax.text(0.95, frac_at_rec - 4, f'{frac_at_rec:.1f}% bajo $r_\\text{{rec}}$',
            transform=ax.get_yaxis_transform(), ha='right',
            fontsize=8, color='gray')
    ax.set_xlabel('Distancia [Å]')
    ax.set_ylabel('Fracción cumulativa [%]')
    ax.set_title('(b) Cumulativa de pares Frenkel')
    ax.set_ylim(0, 105)
    ax.legend(loc='lower right', fontsize=8)

    fig.suptitle(f'Distancias vacancia↔intersticial WS — {label} '
                 f'({len(vacs)} pares)', fontsize=10.5, y=1.04)
    fig.tight_layout()
    save(fig, 'fig7_frenkel_distances', outdir)


# ─── Main ────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser()
    p.add_argument('hybrid_csv', type=Path)
    p.add_argument('ws_sites_csv', type=Path, nargs='?', default=None)
    p.add_argument('--outdir', type=Path, default=Path('figures_hybrid'))
    p.add_argument('--label', type=str, default='HEA finalCool')
    p.add_argument('--recomb-radius', type=float, default=3.3)
    args = p.parse_args()

    set_publication_style()

    print(f'[publication_hybrid] reading {args.hybrid_csv}')
    df = load_hybrid_csv(args.hybrid_csv)
    print(f'  {len(df):,} candidatos | accepted={(df["accepted"]==1).sum():,}')

    fig_ws_vs_hybrid(df, args.label, args.outdir)
    fig_score_distribution(df, args.label, args.outdir)
    fig_signal_breakdown(df, args.label, args.outdir)
    fig_spatial_3d(df, args.label, args.outdir)
    fig_signal_correlation(df, args.label, args.outdir)
    fig_penalties(df, args.label, args.outdir)

    if args.ws_sites_csv and args.ws_sites_csv.exists():
        df_ws = load_ws_sites_csv(args.ws_sites_csv)
        print(f'  WS sites: {len(df_ws):,} '
              f'({(df_ws["occupancy"]==0).sum()} vacancies, '
              f'{(df_ws["occupancy"]>=2).sum()} crowded)')
        fig_frenkel_distances(df, df_ws, args.label, args.outdir,
                                recomb_radius=args.recomb_radius)

    print(f'\n✓ Figuras guardadas en {args.outdir}')


if __name__ == '__main__':
    main()
