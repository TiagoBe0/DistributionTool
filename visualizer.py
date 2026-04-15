"""
DistributionTool Visualizer
===========================
Carga cualquier CSV generado por distool (formato: # id type x y z dist_to_ref defect_prob defect_type)
y permite inspeccionar la muestra en 3D filtrando por distorsión, tipo de defecto y tipo de átomo.

Uso:
    python visualizer.py
    Abrir http://127.0.0.1:8050 en el browser.
"""

import base64
import io

import dash
from dash import Input, Output, State, dcc, html
import pandas as pd
import plotly.graph_objects as go

# ---------------------------------------------------------------------------
# Estado server-side (app local, usuario único)
# ---------------------------------------------------------------------------
_df: pd.DataFrame | None = None

MAX_RENDER = 300_000   # Límite de puntos renderizados; muestra aleatoria si supera


# ---------------------------------------------------------------------------
# Parseo del CSV de distool
# ---------------------------------------------------------------------------
def parse_upload(contents: str, filename: str) -> tuple[pd.DataFrame, str]:
    _, content_string = contents.split(",", 1)
    raw = base64.b64decode(content_string).decode("utf-8")

    lines = [l for l in raw.splitlines() if l.strip()]
    if not lines:
        raise ValueError("Archivo vacío.")

    header = lines[0].lstrip("# ").split()
    rows = [l.split() for l in lines[1:]]
    df = pd.DataFrame(rows, columns=header)

    required = {"x", "y", "z", "dist_to_ref"}
    if not required.issubset(df.columns):
        raise ValueError(f"Columnas requeridas: {required}. Encontradas: {set(df.columns)}")

    for col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col])
        except (ValueError, TypeError):
            pass

    msg = f"{filename}  —  {len(df):,} átomos"
    return df, msg


# ---------------------------------------------------------------------------
# Hover text vectorizado (evita iterrows en 1M filas)
# ---------------------------------------------------------------------------
def build_hover(df: pd.DataFrame) -> pd.Series:
    parts = ["<b>id " + df["id"].astype(str) + "</b>"] if "id" in df.columns else []
    if "type" in df.columns:
        parts.append("tipo: " + df["type"].astype(str))
    parts.append("dist: " + df["dist_to_ref"].map("{:.5f}".format))
    if "defect_prob" in df.columns:
        parts.append("prob: " + df["defect_prob"].map("{:.4f}".format))
    if "defect_type" in df.columns:
        parts.append("defecto: " + df["defect_type"].astype(str))
    return pd.Series(["<br>".join(row) for row in zip(*parts)], index=df.index)


# ---------------------------------------------------------------------------
# App layout
# ---------------------------------------------------------------------------
SIDEBAR_WIDTH = "280px"

app = dash.Dash(__name__, title="DistributionTool Visualizer")
app.layout = html.Div(
    style={"fontFamily": "monospace", "background": "#0e1117", "color": "#e0e0e0", "minHeight": "100vh"},
    children=[
        # ── Header ──────────────────────────────────────────────────────────
        html.Div(
            "DistributionTool Visualizer",
            style={
                "padding": "14px 24px", "fontSize": "18px", "fontWeight": "bold",
                "background": "#161b22", "borderBottom": "1px solid #30363d",
                "letterSpacing": "1px",
            },
        ),
        # ── Upload ──────────────────────────────────────────────────────────
        html.Div(
            style={"padding": "16px 24px", "borderBottom": "1px solid #30363d"},
            children=[
                dcc.Upload(
                    id="upload",
                    children=html.Div(
                        ["Arrastrá o ", html.A("seleccioná", style={"color": "#58a6ff", "cursor": "pointer"}),
                         " el CSV de output de distool"],
                    ),
                    style={
                        "width": "100%", "height": "52px", "lineHeight": "52px",
                        "borderWidth": "1px", "borderStyle": "dashed", "borderRadius": "6px",
                        "borderColor": "#30363d", "textAlign": "center", "color": "#8b949e",
                        "background": "#161b22", "cursor": "pointer",
                    },
                ),
                html.Div(id="file-info", style={"marginTop": "8px", "color": "#8b949e", "fontSize": "13px"}),
                html.Div(id="error-msg", style={"marginTop": "6px", "color": "#f85149", "fontSize": "13px"}),
            ],
        ),
        # ── Body (sidebar + plot) ────────────────────────────────────────────
        html.Div(
            style={"display": "flex", "flex": "1"},
            children=[
                # Sidebar
                html.Div(
                    style={
                        "width": SIDEBAR_WIDTH, "minWidth": SIDEBAR_WIDTH, "padding": "20px 16px",
                        "borderRight": "1px solid #30363d", "background": "#161b22",
                        "overflowY": "auto", "fontSize": "13px",
                    },
                    children=[
                        # Distorsión
                        html.Div([
                            html.Label("Distorsión (dist_to_ref)", style={"color": "#8b949e"}),
                            dcc.RangeSlider(
                                id="dist-slider", min=0, max=1, step=0.0001, value=[0, 1],
                                marks=None,
                                tooltip={"placement": "bottom", "always_visible": True},
                            ),
                        ], style={"marginBottom": "20px"}),

                        # Tipo de defecto
                        html.Div([
                            html.Label("Tipo de defecto", style={"color": "#8b949e"}),
                            dcc.Checklist(
                                id="defect-filter", options=[], value=[],
                                labelStyle={"display": "block", "margin": "4px 0"},
                                inputStyle={"marginRight": "6px"},
                            ),
                        ], style={"marginBottom": "20px"}),

                        # Tipo de átomo
                        html.Div([
                            html.Label("Tipo de átomo", style={"color": "#8b949e"}),
                            dcc.Checklist(
                                id="type-filter", options=[], value=[],
                                labelStyle={"display": "block", "margin": "4px 0"},
                                inputStyle={"marginRight": "6px"},
                            ),
                        ], style={"marginBottom": "20px"}),

                        # Color
                        html.Div([
                            html.Label("Colorear por", style={"color": "#8b949e"}),
                            dcc.Dropdown(
                                id="color-col",
                                options=[
                                    {"label": "dist_to_ref", "value": "dist_to_ref"},
                                    {"label": "defect_prob", "value": "defect_prob"},
                                    {"label": "tipo de átomo", "value": "type"},
                                ],
                                value="dist_to_ref", clearable=False,
                                style={"background": "#0e1117", "color": "#e0e0e0"},
                            ),
                        ], style={"marginBottom": "20px"}),

                        # Tamaño de punto
                        html.Div([
                            html.Label("Tamaño de punto", style={"color": "#8b949e"}),
                            dcc.Slider(
                                id="point-size", min=1, max=10, step=0.5, value=3,
                                marks={i: str(i) for i in [1, 3, 5, 8, 10]},
                                tooltip={"placement": "bottom"},
                            ),
                        ], style={"marginBottom": "20px"}),

                        # Opacidad
                        html.Div([
                            html.Label("Opacidad", style={"color": "#8b949e"}),
                            dcc.Slider(
                                id="opacity", min=0.1, max=1.0, step=0.05, value=0.85,
                                marks={v: str(v) for v in [0.1, 0.5, 1.0]},
                                tooltip={"placement": "bottom"},
                            ),
                        ], style={"marginBottom": "20px"}),

                        # Colorscale
                        html.Div([
                            html.Label("Paleta de colores", style={"color": "#8b949e"}),
                            dcc.Dropdown(
                                id="colorscale",
                                options=[{"label": s, "value": s} for s in [
                                    "Plasma", "Viridis", "Inferno", "Turbo", "RdBu", "Hot"
                                ]],
                                value="Plasma", clearable=False,
                                style={"background": "#0e1117"},
                            ),
                        ], style={"marginBottom": "24px"}),

                        # Info
                        html.Div(id="atom-count", style={"color": "#8b949e", "fontSize": "12px", "lineHeight": "1.6"}),
                    ],
                ),

                # Plot
                html.Div(
                    style={"flex": "1", "padding": "8px"},
                    children=[
                        dcc.Graph(
                            id="scatter-3d",
                            style={"height": "calc(100vh - 160px)"},
                            config={"displayModeBar": True, "scrollZoom": True},
                        )
                    ],
                ),
            ],
        ),
    ],
)


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------

@app.callback(
    Output("file-info", "children"),
    Output("error-msg", "children"),
    Output("dist-slider", "min"),
    Output("dist-slider", "max"),
    Output("dist-slider", "value"),
    Output("defect-filter", "options"),
    Output("defect-filter", "value"),
    Output("type-filter", "options"),
    Output("type-filter", "value"),
    Input("upload", "contents"),
    State("upload", "filename"),
    prevent_initial_call=True,
)
def load_file(contents, filename):
    global _df
    if not contents:
        return "", "", 0, 1, [0, 1], [], [], [], []
    try:
        df, msg = parse_upload(contents, filename)
    except Exception as exc:
        return "", f"Error: {exc}", 0, 1, [0, 1], [], [], [], []

    _df = df

    dmin = float(df["dist_to_ref"].min())
    dmax = float(df["dist_to_ref"].max())

    defect_opts, defect_vals = [], []
    if "defect_type" in df.columns:
        defect_vals = sorted(df["defect_type"].unique().tolist())
        defect_opts = [{"label": v, "value": v} for v in defect_vals]

    type_opts, type_vals = [], []
    if "type" in df.columns:
        type_vals = sorted(df["type"].unique().tolist())
        type_opts = [{"label": f"tipo {v}", "value": v} for v in type_vals]

    return msg, "", dmin, dmax, [dmin, dmax], defect_opts, defect_vals, type_opts, type_vals


@app.callback(
    Output("scatter-3d", "figure"),
    Output("atom-count", "children"),
    Input("dist-slider", "value"),
    Input("defect-filter", "value"),
    Input("type-filter", "value"),
    Input("color-col", "value"),
    Input("point-size", "value"),
    Input("opacity", "value"),
    Input("colorscale", "value"),
    prevent_initial_call=True,
)
def update_plot(dist_range, defect_types, atom_types, color_col, point_size, opacity, colorscale):
    empty = go.Figure().update_layout(
        paper_bgcolor="#0e1117", plot_bgcolor="#0e1117",
        scene=dict(bgcolor="#0e1117"),
    )

    if _df is None:
        return empty, "Sin datos cargados."

    df = _df

    # Filtros
    mask = df["dist_to_ref"].between(dist_range[0], dist_range[1])
    if defect_types and "defect_type" in df.columns:
        mask &= df["defect_type"].isin(defect_types)
    if atom_types and "type" in df.columns:
        mask &= df["type"].isin(atom_types)

    filtered = df[mask]
    total_filtered = len(filtered)

    # Downsampling para el render si supera el límite
    downsampled = False
    if len(filtered) > MAX_RENDER:
        filtered = filtered.sample(MAX_RENDER, random_state=42)
        downsampled = True

    # Color
    if color_col == "type" and "type" in filtered.columns:
        color_values = filtered["type"].astype(str)
        colorbar = None
        use_categorical = True
    else:
        color_values = filtered.get(color_col, filtered["dist_to_ref"])
        colorbar = dict(title=color_col, thickness=14, len=0.6, tickfont=dict(color="#e0e0e0"))
        use_categorical = False

    hover = build_hover(filtered)

    if use_categorical:
        marker = dict(size=point_size, opacity=opacity, color=color_values)
    else:
        marker = dict(
            size=point_size, opacity=opacity,
            color=color_values, colorscale=colorscale,
            colorbar=colorbar, showscale=True,
        )

    fig = go.Figure(data=[go.Scatter3d(
        x=filtered["x"], y=filtered["y"], z=filtered["z"],
        mode="markers",
        marker=marker,
        text=hover,
        hovertemplate="%{text}<extra></extra>",
    )])

    fig.update_layout(
        paper_bgcolor="#0e1117",
        scene=dict(
            bgcolor="#161b22",
            xaxis=dict(title="X (Å)", color="#8b949e", gridcolor="#30363d"),
            yaxis=dict(title="Y (Å)", color="#8b949e", gridcolor="#30363d"),
            zaxis=dict(title="Z (Å)", color="#8b949e", gridcolor="#30363d"),
            aspectmode="data",
        ),
        margin=dict(l=0, r=0, b=0, t=0),
        uirevision="camera",   # preserva la cámara al cambiar filtros
    )

    # Info texto
    info_lines = [f"Visibles: {total_filtered:,} / {len(_df):,} átomos"]
    if downsampled:
        info_lines.append(f"(mostrando {MAX_RENDER:,} puntos aleatorios)")
    count_text = html.Div([html.Div(l) for l in info_lines])

    return fig, count_text


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Abrí http://127.0.0.1:8050 en el browser")
    app.run(debug=False)
