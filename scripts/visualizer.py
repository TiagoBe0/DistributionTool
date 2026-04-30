"""
DistributionTool Visualizer
===========================
Carga cualquier CSV generado por distool y permite inspeccionar los datos.

Formatos soportados:
  output*.csv      → visualización 3D + filtros
  vacancies*.csv   → mapa 3D de vacancias + histograma de clusters (pestaña propia)

Uso:
    python visualizer.py
    Abrir http://127.0.0.1:8050 en el browser.
"""

import base64

import dash
from dash import Input, Output, State, dcc, html
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# ---------------------------------------------------------------------------
# Estado server-side
# ---------------------------------------------------------------------------
_df: pd.DataFrame | None = None
_csv_type: str = "unknown"

MAX_RENDER = 300_000

THEMES = {
    "dark": {
        "--bg": "#0e1117",
        "--panel": "#161b22",
        "--border": "#30363d",
        "--text": "#e0e0e0",
        "--muted": "#8b949e",
        "--accent": "#58a6ff",
        "--error": "#f85149",
        "plot_bg": "#161b22",
        "canvas_bg": "#0e1117",
        "grid": "#30363d",
    },
    "light": {
        "--bg": "#f5f7fb",
        "--panel": "#ffffff",
        "--border": "#d8dee8",
        "--text": "#1f2937",
        "--muted": "#5f6b7a",
        "--accent": "#2563eb",
        "--error": "#c2410c",
        "plot_bg": "#ffffff",
        "canvas_bg": "#f5f7fb",
        "grid": "#d8dee8",
    },
}

TAB_STYLE = {
    "color": "var(--muted)",
    "background": "var(--panel)",
    "border": "none",
    "padding": "10px 24px",
    "fontSize": "14px",
}
TAB_SELECTED = {
    **TAB_STYLE,
    "color": "var(--accent)",
    "borderBottom": "2px solid var(--accent)",
    "fontWeight": "bold",
}

SIDEBAR_WIDTH = "280px"


def app_style(theme_name: str = "dark") -> dict:
    theme = THEMES.get(theme_name, THEMES["dark"])
    return {
        **{k: v for k, v in theme.items() if k.startswith("--")},
        "fontFamily": "monospace",
        "background": "var(--bg)",
        "color": "var(--text)",
        "minHeight": "100vh",
    }


# ---------------------------------------------------------------------------
# Parseo y detección de CSV
# ---------------------------------------------------------------------------
def detect_type(columns: list) -> str:
    cols = set(columns)
    if {"n_grid_pts", "d_near_max"}.issubset(cols):
        return "vacancies"
    if {"x", "y", "z"}.issubset(cols) and "dist_to_ref" not in cols:
        return "vacancies_simple"
    if "dist_to_ref" in cols:
        return "output"
    return "unknown"


def parse_upload(contents: str, filename: str):
    _, content_string = contents.split(",", 1)
    raw = base64.b64decode(content_string).decode("utf-8")
    lines = [l for l in raw.splitlines() if l.strip()]
    if not lines:
        raise ValueError("Archivo vacío.")
    header = lines[0].lstrip("# ").split()
    rows = [l.split() for l in lines[1:]]
    df = pd.DataFrame(rows, columns=header)
    for col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col])
        except (ValueError, TypeError):
            pass
    csv_type = detect_type(list(df.columns))
    return df, f"{filename}  —  {len(df):,} filas  [{csv_type}]", csv_type


# ---------------------------------------------------------------------------
# Layout (estático — todos los componentes existen desde el inicio)
# ---------------------------------------------------------------------------
app = dash.Dash(__name__, title="DistributionTool Visualizer")

app.layout = html.Div(
    id="app-root",
    style=app_style("dark"),
    children=[

        # ── Header ──────────────────────────────────────────────────────────
        html.Div(
            children=[
                html.Div("DistributionTool Visualizer"),
                dcc.RadioItems(
                    id="theme",
                    options=[{"label": "Oscuro", "value": "dark"},
                             {"label": "Claro",  "value": "light"}],
                    value="dark", inline=True,
                    style={"fontSize": "13px", "fontWeight": "normal", "color": "var(--muted)"},
                    inputStyle={"marginLeft": "14px", "marginRight": "5px"},
                ),
            ],
            style={
                "display": "flex", "alignItems": "center",
                "justifyContent": "space-between", "gap": "16px",
                "padding": "14px 24px", "fontSize": "18px", "fontWeight": "bold",
                "background": "var(--panel)", "borderBottom": "1px solid var(--border)",
                "letterSpacing": "1px", "flexWrap": "wrap",
            },
        ),

        # ── Upload ──────────────────────────────────────────────────────────
        html.Div(
            style={"padding": "16px 24px", "borderBottom": "1px solid var(--border)"},
            children=[
                dcc.Upload(
                    id="upload",
                    children=html.Div([
                        "Arrastrá o ",
                        html.A("seleccioná", style={"color": "var(--accent)", "cursor": "pointer"}),
                        " el CSV de distool  (output*.csv  o  vacancies*.csv)",
                    ]),
                    style={
                        "width": "100%", "height": "52px", "lineHeight": "52px",
                        "borderWidth": "1px", "borderStyle": "dashed", "borderRadius": "6px",
                        "borderColor": "var(--border)", "textAlign": "center",
                        "color": "var(--muted)", "background": "var(--panel)", "cursor": "pointer",
                    },
                ),
                html.Div(id="file-info", style={"marginTop": "8px", "color": "var(--muted)", "fontSize": "13px"}),
                html.Div(id="error-msg", style={"marginTop": "6px",  "color": "var(--error)", "fontSize": "13px"}),
            ],
        ),

        # ── Pestañas ────────────────────────────────────────────────────────
        dcc.Tabs(
            id="tabs", value="tab-3d",
            style={"background": "var(--panel)", "borderBottom": "1px solid var(--border)"},
            colors={"border": "var(--border)", "primary": "var(--accent)", "background": "var(--panel)"},
            children=[
                dcc.Tab(label="Vista 3D",               value="tab-3d",
                        style=TAB_STYLE, selected_style=TAB_SELECTED),
                dcc.Tab(label="Histograma de clusters", value="tab-hist",
                        style=TAB_STYLE, selected_style=TAB_SELECTED),
            ],
        ),

        # ── Vista 3D ────────────────────────────────────────────────────────
        html.Div(
            id="view-3d",
            style={"display": "flex", "flex": "1"},
            children=[
                # Sidebar
                html.Div(
                    style={
                        "width": SIDEBAR_WIDTH, "minWidth": SIDEBAR_WIDTH,
                        "padding": "20px 16px",
                        "borderRight": "1px solid var(--border)",
                        "background": "var(--panel)",
                        "overflowY": "auto", "fontSize": "13px",
                    },
                    children=[
                        html.Div([
                            html.Label("Distorsión (dist_to_ref)", style={"color": "var(--muted)"}),
                            dcc.RangeSlider(id="dist-slider", min=0, max=1, step=0.0001,
                                            value=[0, 1], marks=None,
                                            tooltip={"placement": "bottom", "always_visible": True}),
                        ], style={"marginBottom": "20px"}),
                        html.Div([
                            html.Label("Tipo de defecto", style={"color": "var(--muted)"}),
                            dcc.Checklist(id="defect-filter", options=[], value=[],
                                          labelStyle={"display": "block", "margin": "4px 0"},
                                          inputStyle={"marginRight": "6px"}),
                        ], style={"marginBottom": "20px"}),
                        html.Div([
                            html.Label("Tipo de átomo", style={"color": "var(--muted)"}),
                            dcc.Checklist(id="type-filter", options=[], value=[],
                                          labelStyle={"display": "block", "margin": "4px 0"},
                                          inputStyle={"marginRight": "6px"}),
                        ], style={"marginBottom": "20px"}),
                        html.Div([
                            html.Label("Colorear por", style={"color": "var(--muted)"}),
                            dcc.Dropdown(
                                id="color-col",
                                options=[
                                    {"label": "dist_to_ref",   "value": "dist_to_ref"},
                                    {"label": "defect_prob",   "value": "defect_prob"},
                                    {"label": "tipo de átomo", "value": "type"},
                                ],
                                value="dist_to_ref", clearable=False,
                                style={"background": "var(--bg)", "color": "var(--text)"},
                            ),
                        ], style={"marginBottom": "20px"}),
                        html.Div([
                            html.Label("Tamaño de punto", style={"color": "var(--muted)"}),
                            dcc.Slider(id="point-size", min=1, max=10, step=0.5, value=3,
                                       marks={i: str(i) for i in [1, 3, 5, 8, 10]},
                                       tooltip={"placement": "bottom"}),
                        ], style={"marginBottom": "20px"}),
                        html.Div([
                            html.Label("Opacidad", style={"color": "var(--muted)"}),
                            dcc.Slider(id="opacity", min=0.1, max=1.0, step=0.05, value=0.85,
                                       marks={v: str(v) for v in [0.1, 0.5, 1.0]},
                                       tooltip={"placement": "bottom"}),
                        ], style={"marginBottom": "20px"}),
                        html.Div([
                            html.Label("Paleta de colores", style={"color": "var(--muted)"}),
                            dcc.Dropdown(
                                id="colorscale",
                                options=[{"label": s, "value": s}
                                         for s in ["Plasma", "Viridis", "Inferno", "Turbo", "RdBu", "Hot"]],
                                value="Plasma", clearable=False,
                                style={"background": "var(--bg)", "color": "var(--text)"},
                            ),
                        ], style={"marginBottom": "24px"}),
                        html.Div(id="atom-count",
                                 style={"color": "var(--muted)", "fontSize": "12px", "lineHeight": "1.6"}),
                    ],
                ),
                # Gráfico 3D
                html.Div(
                    style={"flex": "1", "padding": "8px"},
                    children=[dcc.Graph(
                        id="scatter-3d",
                        style={"height": "calc(100vh - 200px)"},
                        config={"displayModeBar": True, "scrollZoom": True},
                    )],
                ),
            ],
        ),

        # ── Vista histograma de clusters ─────────────────────────────────────
        html.Div(
            id="view-hist",
            style={"display": "none", "flex": "1", "padding": "24px 40px"},
            children=[
                dcc.Graph(
                    id="cluster-hist",
                    style={"height": "calc(100vh - 210px)", "minHeight": "520px"},
                    config={"displayModeBar": True},
                ),
            ],
        ),

        # Store tipo CSV
        dcc.Store(id="csv-type-store", data="unknown"),
    ],
)


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------

@app.callback(
    Output("app-root", "style"),
    Input("theme", "value"),
)
def update_theme(theme_name):
    return app_style(theme_name)


@app.callback(
    Output("view-3d",   "style"),
    Output("view-hist", "style"),
    Input("tabs", "value"),
)
def switch_tab(tab):
    show    = {"display": "flex", "flex": "1", "padding": "24px 40px"}
    hide    = {"display": "none"}
    show_3d = {"display": "flex", "flex": "1"}
    if tab == "tab-3d":
        return show_3d, hide
    return hide, show


@app.callback(
    Output("file-info",       "children"),
    Output("error-msg",       "children"),
    Output("dist-slider",     "min"),
    Output("dist-slider",     "max"),
    Output("dist-slider",     "value"),
    Output("defect-filter",   "options"),
    Output("defect-filter",   "value"),
    Output("type-filter",     "options"),
    Output("type-filter",     "value"),
    Output("csv-type-store",  "data"),
    Input("upload", "contents"),
    State("upload", "filename"),
    prevent_initial_call=True,
)
def load_file(contents, filename):
    global _df, _csv_type
    if not contents:
        return ("", "", 0, 1, [0, 1], [], [], [], [], "unknown")
    try:
        df, msg, csv_type = parse_upload(contents, filename)
    except Exception as exc:
        return ("", f"Error: {exc}", 0, 1, [0, 1], [], [], [], [], "unknown")

    _df = df
    _csv_type = csv_type

    dmin, dmax = 0.0, 1.0
    defect_opts, defect_vals = [], []
    type_opts, type_vals = [], []

    if csv_type == "output":
        dmin = float(df["dist_to_ref"].min())
        dmax = float(df["dist_to_ref"].max())
        if "defect_type" in df.columns:
            defect_vals = sorted(df["defect_type"].unique().tolist())
            defect_opts = [{"label": v, "value": v} for v in defect_vals]
        if "type" in df.columns:
            type_vals = sorted(df["type"].unique().tolist())
            type_opts = [{"label": f"tipo {v}", "value": v} for v in type_vals]

    return (msg, "", dmin, dmax, [dmin, dmax],
            defect_opts, defect_vals, type_opts, type_vals, csv_type)


@app.callback(
    Output("scatter-3d", "figure"),
    Output("atom-count", "children"),
    Input("dist-slider",   "value"),
    Input("defect-filter", "value"),
    Input("type-filter",   "value"),
    Input("color-col",     "value"),
    Input("point-size",    "value"),
    Input("opacity",       "value"),
    Input("colorscale",    "value"),
    Input("theme",         "value"),
    prevent_initial_call=True,
)
def update_3d(dist_range, defect_types, atom_types,
              color_col, point_size, opacity, colorscale, theme_name):
    theme = THEMES.get(theme_name, THEMES["dark"])
    empty = go.Figure().update_layout(
        paper_bgcolor=theme["canvas_bg"], plot_bgcolor=theme["canvas_bg"],
        scene=dict(bgcolor=theme["canvas_bg"]),
    )
    if _df is None:
        return empty, "Sin datos cargados."

    df = _df

    if _csv_type in ("vacancies", "vacancies_simple"):
        sample = df if len(df) <= MAX_RENDER else df.sample(MAX_RENDER, random_state=42)
        color_vals = sample["d_near_max"].values if "d_near_max" in sample.columns else sample["z"].values
        if "n_grid_pts" in sample.columns:
            s_raw  = np.sqrt(sample["n_grid_pts"].values.astype(float))
            s_norm = (s_raw - s_raw.min()) / (s_raw.max() - s_raw.min() + 1e-14)
            msize  = (4 + 10 * s_norm).tolist()
        else:
            msize = 4
        hover = (
            "x=" + sample["x"].map("{:.2f}".format) + "<br>" +
            "y=" + sample["y"].map("{:.2f}".format) + "<br>" +
            "z=" + sample["z"].map("{:.2f}".format) +
            ("<br>d_near_max=" + sample["d_near_max"].map("{:.3f}".format) if "d_near_max" in sample.columns else "") +
            ("<br>n_grid_pts=" + sample["n_grid_pts"].astype(str) if "n_grid_pts" in sample.columns else "")
        )
        fig = go.Figure(data=[go.Scatter3d(
            x=sample["x"], y=sample["y"], z=sample["z"],
            mode="markers",
            marker=dict(size=msize, color=color_vals, colorscale="Hot_r", opacity=0.80,
                        colorbar=dict(title="d_near_max",
                                      tickfont=dict(color=theme["--text"]),
                                      titlefont=dict(color=theme["--text"]))),
            text=hover, hovertemplate="%{text}<extra></extra>",
        )])
        info = html.Div(f"Vacancias: {len(df):,}")

    else:
        mask = df["dist_to_ref"].between(dist_range[0], dist_range[1])
        if defect_types and "defect_type" in df.columns:
            mask &= df["defect_type"].isin(defect_types)
        if atom_types and "type" in df.columns:
            mask &= df["type"].isin(atom_types)
        filtered = df[mask]
        downsampled = len(filtered) > MAX_RENDER
        if downsampled:
            filtered = filtered.sample(MAX_RENDER, random_state=42)

        if color_col == "type" and "type" in filtered.columns:
            marker = dict(size=point_size, opacity=opacity,
                          color=filtered["type"].astype(str))
        else:
            cv = filtered.get(color_col, filtered["dist_to_ref"])
            marker = dict(size=point_size, opacity=opacity,
                          color=cv, colorscale=colorscale, showscale=True,
                          colorbar=dict(title=color_col, thickness=14, len=0.6,
                                        tickfont=dict(color=theme["--text"]),
                                        titlefont=dict(color=theme["--text"])))

        parts = []
        if "id" in filtered.columns:
            parts.append("<b>id " + filtered["id"].astype(str) + "</b>")
        if "type" in filtered.columns:
            parts.append("tipo: " + filtered["type"].astype(str))
        parts.append("dist: " + filtered["dist_to_ref"].map("{:.5f}".format))
        if "defect_prob" in filtered.columns:
            parts.append("prob: " + filtered["defect_prob"].map("{:.4f}".format))
        if "defect_type" in filtered.columns:
            parts.append("defecto: " + filtered["defect_type"].astype(str))
        hover = pd.Series(["<br>".join(row) for row in zip(*parts)], index=filtered.index)

        fig = go.Figure(data=[go.Scatter3d(
            x=filtered["x"], y=filtered["y"], z=filtered["z"],
            mode="markers", marker=marker,
            text=hover, hovertemplate="%{text}<extra></extra>",
        )])
        lines = [f"Visibles: {len(df[mask]):,} / {len(_df):,} átomos"]
        if downsampled:
            lines.append(f"(mostrando {MAX_RENDER:,} puntos aleatorios)")
        info = html.Div([html.Div(l) for l in lines])

    fig.update_layout(
        paper_bgcolor=theme["canvas_bg"],
        font=dict(color=theme["--text"]),
        scene=dict(
            bgcolor=theme["plot_bg"],
            xaxis=dict(title="X (Å)", color=theme["--muted"], gridcolor=theme["grid"]),
            yaxis=dict(title="Y (Å)", color=theme["--muted"], gridcolor=theme["grid"]),
            zaxis=dict(title="Z (Å)", color=theme["--muted"], gridcolor=theme["grid"]),
            aspectmode="data",
        ),
        margin=dict(l=0, r=0, b=0, t=0),
        uirevision="camera",
    )
    return fig, info


@app.callback(
    Output("cluster-hist", "figure"),
    Input("csv-type-store", "data"),
    Input("theme",          "value"),
    prevent_initial_call=True,
)
def update_cluster_hist(csv_type, theme_name):
    theme = THEMES.get(theme_name, THEMES["dark"])
    empty = go.Figure().update_layout(
        paper_bgcolor=theme["canvas_bg"], plot_bgcolor=theme["plot_bg"],
    )
    if _df is None or "n_grid_pts" not in _df.columns:
        return empty

    n_pts    = _df["n_grid_pts"].values.astype(int)
    n_total  = len(n_pts)
    median   = float(np.median(n_pts))
    mean     = float(np.mean(n_pts))
    OVF      = 101

    n_clipped    = np.where(n_pts > 100, OVF, n_pts)
    overflow_n   = int((n_pts > 100).sum())
    pct_overflow = 100.0 * overflow_n / n_total if n_total > 0 else 0.0

    counts, _ = np.histogram(n_clipped, bins=np.arange(0.5, OVF + 1.5, 1))
    x_vals    = list(range(1, OVF + 1))
    colors    = ["#d73027" if x == OVF else "#2166ac" for x in x_vals]

    text_labels = [""] * len(x_vals)
    if overflow_n > 0:
        text_labels[-1] = f"{overflow_n}<br>({pct_overflow:.1f}%)"

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=x_vals, y=counts.tolist(),
        marker_color=colors, marker_line_width=0,
        text=text_labels, textposition="outside",
        textfont=dict(size=14, color="#d73027"),
        hovertemplate="n_grid_pts: %{x}<br>clusters: %{y}<extra></extra>",
    ))

    if median <= 100:
        fig.add_vline(x=median, line_dash="dash", line_color="#ef5350", line_width=2.5,
                      annotation_text=f"  mediana = {median:.0f}",
                      annotation_font_size=14, annotation_font_color="#ef5350",
                      annotation_position="top right")

    fig.add_vline(x=min(mean, OVF), line_dash="dot", line_color="#66bb6a", line_width=2.5,
                  annotation_text=f"  media = {mean:.1f}",
                  annotation_font_size=14, annotation_font_color="#66bb6a",
                  annotation_position="bottom right")

    fig.add_vline(x=OVF - 0.5, line_color="#555555", line_width=1.5, line_dash="dot")

    stats = (
        f"<b>N total:</b> {n_total:,}<br>"
        f"<b>Mediana:</b> {median:.0f} pts<br>"
        f"<b>Media:</b> {mean:.1f} pts<br>"
        f"<b>Máx:</b> {int(n_pts.max())} pts<br>"
        f"<b>> 100 pts:</b> {overflow_n:,} ({pct_overflow:.1f}%)"
    )
    fig.add_annotation(
        xref="paper", yref="paper", x=0.99, y=0.97,
        text=stats, align="left", showarrow=False,
        font=dict(size=14, color=theme["--text"]),
        bgcolor=theme["--panel"], bordercolor=theme["--border"],
        borderwidth=1, borderpad=12,
    )

    tick_vals   = list(range(10, 101, 10)) + [OVF]
    tick_labels = [str(t) for t in range(10, 101, 10)] + ["100+"]

    fig.update_layout(
        paper_bgcolor=theme["canvas_bg"],
        plot_bgcolor=theme["plot_bg"],
        font=dict(color=theme["--text"], size=14),
        title=dict(
            text=f"Distribución de tamaño de clusters de vacancias  (N = {n_total:,})",
            font=dict(size=18), x=0.5,
        ),
        xaxis=dict(
            tickvals=tick_vals, ticktext=tick_labels, tickfont=dict(size=13),
            title=dict(text="Puntos de grilla por cluster de vacancia (n_grid_pts)", font=dict(size=15)),
            gridcolor=theme["grid"], gridwidth=0.5,
            range=[0.5, OVF + 1],
        ),
        yaxis=dict(
            title=dict(text="Número de clusters", font=dict(size=15)),
            gridcolor=theme["grid"], gridwidth=0.5,
            tickfont=dict(size=13),
        ),
        bargap=0.05,
        margin=dict(l=80, r=50, t=80, b=80),
        showlegend=False,
    )
    return fig


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Abrí http://127.0.0.1:8050 en el browser")
    app.run(debug=False)
