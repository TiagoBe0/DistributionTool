"""
Script to generate DistributionTool PowerPoint presentation.
"""

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt
import copy

# ── Color palette ──────────────────────────────────────────────────────────────
C_DARK_BG   = RGBColor(0x1A, 0x1A, 0x2E)   # dark navy
C_MID_BG    = RGBColor(0x16, 0x21, 0x3E)   # slightly lighter navy
C_ACCENT    = RGBColor(0x0F, 0x3C, 0x78)   # deep blue accent
C_HIGHLIGHT = RGBColor(0x00, 0xB4, 0xD8)   # cyan highlight
C_ORANGE    = RGBColor(0xFF, 0x6B, 0x35)   # orange accent
C_GREEN     = RGBColor(0x06, 0xD6, 0xA0)   # green accent
C_WHITE     = RGBColor(0xFF, 0xFF, 0xFF)
C_LIGHT_GY  = RGBColor(0xC8, 0xD6, 0xE5)
C_YELLOW    = RGBColor(0xFF, 0xD1, 0x66)

prs = Presentation()
prs.slide_width  = Inches(13.33)
prs.slide_height = Inches(7.5)

BLANK = prs.slide_layouts[6]   # completely blank layout


# ── Helper functions ───────────────────────────────────────────────────────────

def add_rect(slide, l, t, w, h, fill_rgb, alpha=None):
    shape = slide.shapes.add_shape(1, Inches(l), Inches(t), Inches(w), Inches(h))
    shape.line.fill.background()
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill_rgb
    return shape


def add_text(slide, text, l, t, w, h,
             font_size=18, bold=False, color=C_WHITE,
             align=PP_ALIGN.LEFT, wrap=True, italic=False):
    txBox = slide.shapes.add_textbox(Inches(l), Inches(t), Inches(w), Inches(h))
    txBox.word_wrap = wrap
    tf = txBox.text_frame
    tf.word_wrap = wrap
    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size = Pt(font_size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    return txBox


def add_multiline(slide, lines, l, t, w, h,
                  font_size=16, color=C_WHITE, bold_first=False,
                  line_spacing=None, align=PP_ALIGN.LEFT):
    """lines: list of (text, bold, color_override_or_None)"""
    txBox = slide.shapes.add_textbox(Inches(l), Inches(t), Inches(w), Inches(h))
    txBox.word_wrap = True
    tf = txBox.text_frame
    tf.word_wrap = True
    for i, item in enumerate(lines):
        if isinstance(item, str):
            text, bold, col = item, False, color
        else:
            text = item[0]
            bold = item[1] if len(item) > 1 else False
            col  = item[2] if len(item) > 2 and item[2] else color
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        run = p.add_run()
        run.text = text
        run.font.size = Pt(font_size)
        run.font.bold = bold
        run.font.color.rgb = col
    return txBox


def slide_bg(slide, color=C_DARK_BG):
    add_rect(slide, 0, 0, 13.33, 7.5, color)


def accent_bar(slide, color=C_HIGHLIGHT):
    add_rect(slide, 0, 0, 13.33, 0.08, color)
    add_rect(slide, 0, 7.42, 13.33, 0.08, color)


def section_header(slide, text, sub=""):
    slide_bg(slide)
    add_rect(slide, 0, 0, 13.33, 0.08, C_HIGHLIGHT)
    add_rect(slide, 0, 7.42, 13.33, 0.08, C_HIGHLIGHT)
    add_rect(slide, 0, 2.5, 13.33, 2.5, C_ACCENT)
    add_text(slide, text, 0.5, 2.7, 12.33, 1.2,
             font_size=40, bold=True, color=C_HIGHLIGHT, align=PP_ALIGN.CENTER)
    if sub:
        add_text(slide, sub, 0.5, 3.9, 12.33, 0.8,
                 font_size=20, color=C_LIGHT_GY, align=PP_ALIGN.CENTER)


def slide_title(slide, title, subtitle=""):
    add_rect(slide, 0, 0, 13.33, 0.06, C_HIGHLIGHT)
    add_text(slide, title, 0.5, 0.15, 12.0, 0.7,
             font_size=26, bold=True, color=C_HIGHLIGHT)
    if subtitle:
        add_text(slide, subtitle, 0.5, 0.75, 12.0, 0.4,
                 font_size=14, color=C_LIGHT_GY)
    add_rect(slide, 0, 7.4, 13.33, 0.1, C_ACCENT)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 1 — TITLE
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
add_rect(s, 0, 0, 13.33, 0.1, C_HIGHLIGHT)
add_rect(s, 0, 7.4, 13.33, 0.1, C_HIGHLIGHT)

# Central decorative block
add_rect(s, 1.0, 1.5, 11.33, 4.2, C_ACCENT)
add_rect(s, 1.0, 1.5, 0.08, 4.2, C_HIGHLIGHT)

add_text(s, "DistributionTool", 1.5, 1.7, 10.5, 1.3,
         font_size=54, bold=True, color=C_HIGHLIGHT, align=PP_ALIGN.CENTER)

add_text(s, "Detección y Clasificación Automática de Defectos\nen Materiales Cristalinos por Simulación Molecular",
         1.5, 3.0, 10.5, 1.2,
         font_size=22, color=C_WHITE, align=PP_ALIGN.CENTER)

add_text(s, "Pipeline completo · SOAP Descriptors · Clasificación estadística · Detección de vacancias",
         1.5, 4.2, 10.5, 0.6,
         font_size=14, color=C_LIGHT_GY, align=PP_ALIGN.CENTER, italic=True)

add_text(s, "C++17  ·  Eigen3  ·  OpenMP  ·  LAMMPS",
         1.5, 5.0, 10.5, 0.5,
         font_size=14, color=C_YELLOW, align=PP_ALIGN.CENTER)

add_text(s, "Investigación en Fusión y Fisión Nuclear · Ciencia de Materiales Computacional",
         1.0, 6.6, 11.33, 0.5,
         font_size=12, color=C_LIGHT_GY, align=PP_ALIGN.CENTER, italic=True)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 2 — AGENDA / ÍNDICE
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
slide_title(s, "Contenido de la Presentación")

topics = [
    ("1", "Motivación y Problema", C_ORANGE),
    ("2", "Fundamentos Científicos", C_HIGHLIGHT),
    ("3", "Pipeline Completo (4 etapas)", C_GREEN),
    ("4", "Descriptores SOAP — Núcleo del Método", C_YELLOW),
    ("5", "Clasificación Estadística de Defectos", C_ORANGE),
    ("6", "Detección de Vacancias", C_HIGHLIGHT),
    ("7", "Arquitectura del Software", C_GREEN),
    ("8", "Configuración y CLI", C_YELLOW),
    ("9", "Rendimiento y Escalabilidad", C_ORANGE),
    ("10", "Flujo de Trabajo en Práctica", C_HIGHLIGHT),
]

cols = [topics[:5], topics[5:]]
for ci, col in enumerate(cols):
    x = 0.8 + ci * 6.3
    for ri, (num, label, color) in enumerate(col):
        y = 1.4 + ri * 1.05
        add_rect(s, x, y, 0.55, 0.65, color)
        add_text(s, num, x, y, 0.55, 0.65,
                 font_size=22, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)
        add_text(s, label, x + 0.65, y + 0.1, 5.2, 0.55,
                 font_size=17, color=C_WHITE)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 3 — MOTIVACIÓN
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
slide_title(s, "Motivación y Problema", "¿Por qué necesitamos esta herramienta?")

# Left problem box
add_rect(s, 0.4, 1.2, 5.8, 5.7, C_ACCENT)
add_rect(s, 0.4, 1.2, 5.8, 0.5, C_ORANGE)
add_text(s, "⚠  El Problema", 0.4, 1.2, 5.8, 0.5,
         font_size=17, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)

prob_lines = [
    ("Simulaciones de dinámica molecular de impactos de neutrones/iones en materiales cristalinos generan cascadas de colisiones.", False, C_LIGHT_GY),
    ("", False, None),
    ("Los métodos clásicos fallan a temperaturas elevadas:", True, C_WHITE),
    ("  • Celdas de Wigner-Seitz: requieren posiciones de referencia exactas", False, C_LIGHT_GY),
    ("  • Teselación de Voronoi: sensible al ruido térmico", False, C_LIGHT_GY),
    ("  • A T > 300 K las vibraciones térmicas enmascaran los defectos reales", False, C_LIGHT_GY),
    ("", False, None),
    ("Resultado: detección errónea → estadísticas de defectos incorrectas → predicciones de daño por radiación inexactas", False, C_ORANGE),
]
add_multiline(s, prob_lines, 0.6, 1.8, 5.4, 5.0, font_size=14, color=C_WHITE)

# Right solution box
add_rect(s, 6.8, 1.2, 6.1, 5.7, C_MID_BG)
add_rect(s, 6.8, 1.2, 6.1, 0.5, C_GREEN)
add_text(s, "✓  La Solución", 6.8, 1.2, 6.1, 0.5,
         font_size=17, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)

sol_lines = [
    ("Descriptores SOAP (Smooth Overlap of Atomic Positions)", True, C_GREEN),
    ("", False, None),
    ("Invariantes rotacionales del entorno local de cada átomo → robustez frente a ruido térmico", False, C_LIGHT_GY),
    ("", False, None),
    ("No requiere posiciones de red de referencia", True, C_YELLOW),
    ("", False, None),
    ("Modelo probabilístico basado en distribución chi → umbral estadístico fundamentado", False, C_LIGHT_GY),
    ("", False, None),
    ("Detección geométrica de vacancias mediante muestreo de grilla uniforme", False, C_LIGHT_GY),
    ("", False, None),
    ("Validación independiente con método Wigner-Seitz opcional", False, C_LIGHT_GY),
]
add_multiline(s, sol_lines, 7.0, 1.8, 5.7, 5.0, font_size=14, color=C_WHITE)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 4 — FUNDAMENTOS CIENTÍFICOS
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
slide_title(s, "Fundamentos Científicos", "Base teórica del método")

refs = [
    ("SOAP Descriptors", "Bartók et al. (2013)", "Power spectrum del solapamiento de densidades atómicas gaussianas. Invariante bajo rotaciones, reflexiones y permutaciones.", C_HIGHLIGHT),
    ("Clasificación χ²", "Domínguez-Gutiérrez & von Toussaint (2019)", "Distancias en el espacio de descriptores siguen una distribución chi. Parámetros ajustados por método de momentos.", C_ORANGE),
    ("Flujo FaVaD", "von Toussaint et al. (2020)", "Vacancies and Defects: detección de vacancias sin posiciones de referencia + clustering greedy de puntos vacíos.", C_GREEN),
    ("Herramientas de Validación", "LAMMPS / OVITO", "Archivos dump de entrada; salida analizada compatible con OVITO para visualización 3D interactiva.", C_YELLOW),
]

for i, (title, ref, desc, color) in enumerate(refs):
    x = 0.4 + (i % 2) * 6.45
    y = 1.3 + (i // 2) * 2.8
    add_rect(s, x, y, 6.1, 2.5, C_ACCENT)
    add_rect(s, x, y, 6.1, 0.08, color)
    add_text(s, title, x + 0.15, y + 0.12, 5.8, 0.5,
             font_size=17, bold=True, color=color)
    add_text(s, ref, x + 0.15, y + 0.6, 5.8, 0.35,
             font_size=12, italic=True, color=C_LIGHT_GY)
    add_text(s, desc, x + 0.15, y + 1.0, 5.8, 1.3,
             font_size=13, color=C_WHITE)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 5 — PIPELINE OVERVIEW
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
slide_title(s, "Pipeline Completo — Vista General")

stages = [
    ("ETAPA 1", "Lectura\nde Entrada", C_HIGHLIGHT, "LAMMPS dump\n(2 archivos)"),
    ("ETAPA 2", "Descriptores\nSOAP", C_ORANGE, "Vectores de\nfirma SOAP"),
    ("ETAPA 3", "Clasificación\nde Defectos", C_GREEN, "Tipo de\ndefecto por\nátomo"),
    ("ETAPA 4", "Detección de\nVacancias", C_YELLOW, "Clusters de\nvacancias"),
]

box_w = 2.4
box_h = 2.8
gap   = 0.55
start_x = 0.55
y_box = 1.8

for i, (stage, name, color, out) in enumerate(stages):
    bx = start_x + i * (box_w + gap)
    # Main box
    add_rect(s, bx, y_box, box_w, box_h, C_ACCENT)
    add_rect(s, bx, y_box, box_w, 0.45, color)
    add_text(s, stage, bx, y_box, box_w, 0.45,
             font_size=12, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)
    add_text(s, name, bx, y_box + 0.5, box_w, 1.2,
             font_size=16, bold=True, color=C_WHITE, align=PP_ALIGN.CENTER)
    # Output label
    add_rect(s, bx + 0.2, y_box + 1.85, box_w - 0.4, 0.8, C_MID_BG)
    add_text(s, out, bx + 0.2, y_box + 1.85, box_w - 0.4, 0.8,
             font_size=11, color=color, align=PP_ALIGN.CENTER)

    # Arrow between boxes
    if i < 3:
        ax = bx + box_w + 0.05
        add_rect(s, ax, y_box + 1.1, gap - 0.1, 0.08, C_HIGHLIGHT)
        add_text(s, "▶", ax + gap * 0.1, y_box + 0.9, 0.4, 0.4,
                 font_size=18, color=C_HIGHLIGHT, align=PP_ALIGN.CENTER)

# Output files section
add_rect(s, 0.4, 5.1, 12.5, 1.9, C_MID_BG)
add_rect(s, 0.4, 5.1, 12.5, 0.08, C_GREEN)
add_text(s, "Archivos de Salida", 0.6, 5.15, 12.0, 0.4,
         font_size=13, bold=True, color=C_GREEN)

outputs = [
    "output.csv\n(átomos+defectos)",
    "vacancies_output.csv\n(clusters vacancias)",
    "hist_output.csv\n(histograma dist.)",
    "pca_output.csv\n(proyección PCA)",
    "analyzed.dump\n(OVITO-compatible)",
]
for i, out in enumerate(outputs):
    ox = 0.6 + i * 2.48
    add_rect(s, ox, 5.6, 2.2, 1.2, C_ACCENT)
    add_text(s, out, ox + 0.05, 5.65, 2.1, 1.1,
             font_size=11, color=C_YELLOW, align=PP_ALIGN.CENTER)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 6 — ETAPA 1: LECTURA DE ENTRADA
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
slide_title(s, "Etapa 1: Lectura de Entrada", "LammpsDumpReader — Parsing de formato LAMMPS dump")

add_rect(s, 0.4, 1.2, 5.6, 5.8, C_MID_BG)
add_rect(s, 0.4, 1.2, 5.6, 0.45, C_HIGHLIGHT)
add_text(s, "Formato LAMMPS Dump", 0.4, 1.2, 5.6, 0.45,
         font_size=15, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)

dump_content = [
    ("ITEM: TIMESTEP", True, C_YELLOW),
    ("1000", False, C_LIGHT_GY),
    ("ITEM: NUMBER OF ATOMS", True, C_YELLOW),
    ("10000", False, C_LIGHT_GY),
    ("ITEM: BOX BOUNDS pp pp pp", True, C_YELLOW),
    ("0.0  28.7  0.0  28.7  0.0  28.7", False, C_LIGHT_GY),
    ("ITEM: ATOMS id type x y z", True, C_YELLOW),
    ("1  1  2.145  0.000  0.000", False, C_LIGHT_GY),
    ("2  1  0.000  2.145  0.000", False, C_LIGHT_GY),
    ("...", False, C_LIGHT_GY),
]
add_multiline(s, dump_content, 0.6, 1.75, 5.2, 5.0, font_size=13, color=C_WHITE)

add_rect(s, 6.4, 1.2, 6.5, 5.8, C_ACCENT)
add_rect(s, 6.4, 1.2, 6.5, 0.45, C_ORANGE)
add_text(s, "Módulos Involucrados", 6.4, 1.2, 6.5, 0.45,
         font_size=15, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)

mod_lines = [
    ("LammpsDumpReader.h", True, C_ORANGE),
    ("  readNext() → un frame a la vez", False, C_LIGHT_GY),
    ("  readAll()  → todos los frames", False, C_LIGHT_GY),
    ("", False, None),
    ("AtomData.h — Estructuras de datos", True, C_GREEN),
    ("  Atom: id, type, x, y, z", False, C_LIGHT_GY),
    ("        descriptor_vector (SOAP)", False, C_LIGHT_GY),
    ("        dist_to_ref, defect_type", False, C_LIGHT_GY),
    ("", False, None),
    ("  SimBox: límites + condiciones", False, C_LIGHT_GY),
    ("          de contorno periódicas (PBC)", False, C_LIGHT_GY),
    ("", False, None),
    ("  Frame: timestep + lista de átomos", False, C_LIGHT_GY),
    ("", False, None),
    ("Entradas requeridas:", True, C_YELLOW),
    ("  1. Frame de referencia (prístino/", False, C_LIGHT_GY),
    ("     termalizado al mismo T)", False, C_LIGHT_GY),
    ("  2. Frame irradiado (dañado)", False, C_LIGHT_GY),
]
add_multiline(s, mod_lines, 6.6, 1.75, 6.1, 5.0, font_size=13, color=C_WHITE)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 7 — ETAPA 2: DESCRIPTORES SOAP (parte 1)
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
slide_title(s, "Etapa 2: Descriptores SOAP", "Smooth Overlap of Atomic Positions — Núcleo del método")

# Left: Mathematical explanation
add_rect(s, 0.4, 1.2, 6.0, 5.8, C_ACCENT)
add_rect(s, 0.4, 1.2, 6.0, 0.45, C_HIGHLIGHT)
add_text(s, "Formulación Matemática", 0.4, 1.2, 6.0, 0.45,
         font_size=15, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)

math_lines = [
    ("1. Densidad atómica gaussiana:", True, C_YELLOW),
    ("   ρ(r) = Σⱼ exp(-|r - rⱼ|² / 2σ²) · fcut(rⱼ)", False, C_LIGHT_GY),
    ("", False, None),
    ("2. Expansión en bases:", True, C_YELLOW),
    ("   ρ(r) = Σ_{nlm} c_{nlm} · φₙ(r) · Yₗᵐ(r̂)", False, C_LIGHT_GY),
    ("   φₙ: bases radiales gaussianas", False, C_LIGHT_GY),
    ("   Yₗᵐ: armónicos esféricos reales", False, C_LIGHT_GY),
    ("", False, None),
    ("3. Power Spectrum (invariante rotacional):", True, C_YELLOW),
    ("   p_{nn'l} = π√(8/2l+1) Σₘ c_{nlm}·c_{n'lm}", False, C_LIGHT_GY),
    ("", False, None),
    ("4. Normalización:", True, C_YELLOW),
    ("   q̃ = p / ||p||₂  (esfera unitaria)", False, C_LIGHT_GY),
    ("", False, None),
    ("Dimensión del vector:", True, C_GREEN),
    ("   D = (n_max×(n_max+1)/2) × (l_max+1)", False, C_LIGHT_GY),
    ("   Con n_max=9, l_max=9: D = 450", False, C_GREEN),
]
add_multiline(s, math_lines, 0.6, 1.75, 5.6, 5.2, font_size=13, color=C_WHITE)

# Right: Implementation
add_rect(s, 6.8, 1.2, 6.1, 5.8, C_MID_BG)
add_rect(s, 6.8, 1.2, 6.1, 0.45, C_ORANGE)
add_text(s, "Implementación C++", 6.8, 1.2, 6.1, 0.45,
         font_size=15, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)

impl_lines = [
    ("CellList.h — Índice espacial O(N)", True, C_ORANGE),
    ("  Lista enlazada de celdas", False, C_LIGHT_GY),
    ("  Stencil 27 celdas para PBC", False, C_LIGHT_GY),
    ("  Consulta vecinos: O(27ρ r_cut³)", False, C_LIGHT_GY),
    ("", False, None),
    ("RadialBasis.h — φₙ(r)", True, C_ORANGE),
    ("  Centros equidistantes en [0, r_cut]", False, C_LIGHT_GY),
    ("  Pre-calcula 1/(2σ²) → hot path", False, C_LIGHT_GY),
    ("  fcut(r) suave para continuidad C¹", False, C_LIGHT_GY),
    ("", False, None),
    ("SphericalHarmonics.h — Yₗᵐ(θ,φ)", True, C_ORANGE),
    ("  Recurrencia de Bonnet: O(l_max²)", False, C_LIGHT_GY),
    ("  Recurrencia trig: evita atan2", False, C_LIGHT_GY),
    ("", False, None),
    ("SOAPDescriptor.h — Orquestador", True, C_ORANGE),
    ("  compute() — un átomo", False, C_LIGHT_GY),
    ("  computeAll() — frame completo", False, C_LIGHT_GY),
    ("  #pragma omp parallel for (OpenMP)", False, C_GREEN),
]
add_multiline(s, impl_lines, 7.0, 1.75, 5.7, 5.2, font_size=13, color=C_WHITE)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 8 — ETAPA 3: CLASIFICACIÓN
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
slide_title(s, "Etapa 3: Clasificación Estadística de Defectos")

# Flow diagram
flow_items = [
    ("Construir Referencia q̄(T)", C_HIGHLIGHT, "Media de DVs del\nframe prístino"),
    ("Ajustar Dist. Chi\n(k, σ)", C_ORANGE, "Método de\nmoments → DOF efectivos"),
    ("Calcular dⁱ = ||q̃ⁱ - q̄(T)||", C_GREEN, "Distancia euclídea\nen esfera unitaria"),
    ("Clasificar Átomo", C_YELLOW, "dⁱ < umbral →\nRed cristalina\ndⁱ ≥ umbral →\nDefecto"),
]

for i, (label, color, sub) in enumerate(flow_items):
    bx = 0.4 + i * 3.15
    add_rect(s, bx, 1.3, 2.9, 2.0, C_ACCENT)
    add_rect(s, bx, 1.3, 2.9, 0.08, color)
    add_text(s, label, bx + 0.1, 1.4, 2.7, 0.9,
             font_size=14, bold=True, color=color, align=PP_ALIGN.CENTER)
    add_text(s, sub, bx + 0.1, 2.3, 2.7, 0.9,
             font_size=11, color=C_LIGHT_GY, align=PP_ALIGN.CENTER)
    if i < 3:
        add_text(s, "▶", bx + 2.9, 1.9, 0.25, 0.4,
                 font_size=18, color=C_HIGHLIGHT, align=PP_ALIGN.CENTER)

# Defect types
add_rect(s, 0.4, 3.7, 12.5, 3.4, C_MID_BG)
add_rect(s, 0.4, 3.7, 12.5, 0.4, C_GREEN)
add_text(s, "Tipos de Defectos", 0.4, 3.7, 12.5, 0.4,
         font_size=15, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)

defect_types = [
    ("Lattice", C_GREEN, "Red perfecta\ndⁱ < umbral"),
    ("SIA\n(Self-Interstitial)", C_ORANGE, "Átomo intersticial\nmás cerca de DV-SIA"),
    ("VacAdj\n(Vacancy-Adjacent)", C_HIGHLIGHT, "Átomo vecino a\nvacancia"),
    ("Type-A", C_YELLOW, "Defecto descubierto\npor PCA"),
    ("Unknown", C_LIGHT_GY, "No clasificado\n(sin referencia)"),
]
for i, (name, color, desc) in enumerate(defect_types):
    bx = 0.6 + i * 2.45
    add_rect(s, bx, 4.2, 2.2, 2.5, C_ACCENT)
    add_rect(s, bx, 4.2, 2.2, 0.4, color)
    add_text(s, name, bx + 0.05, 4.2, 2.1, 0.4,
             font_size=12, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)
    add_text(s, desc, bx + 0.05, 4.7, 2.1, 1.8,
             font_size=12, color=C_WHITE, align=PP_ALIGN.CENTER)

# Probability formula
add_text(s, "Probabilidad de defecto: P_defecto = 1 − P_chi(dⁱ; k, σ)",
         0.5, 7.0, 12.0, 0.4,
         font_size=13, italic=True, color=C_YELLOW, align=PP_ALIGN.CENTER)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 9 — ETAPA 4: VACANCIAS
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
slide_title(s, "Etapa 4: Detección de Vacancias", "Grid sampling + Greedy clustering")

# Algorithm steps
steps = [
    ("1. Muestreo\nde Grilla", C_HIGHLIGHT,
     "Puntos uniformes Δ-espaciados\nen toda la caja de simulación"),
    ("2. Campo de\nDistancia", C_ORANGE,
     "Para cada punto: distancia\nal átomo más cercano via CellList"),
    ("3. Puntos\nVacantes", C_GREEN,
     "Marcar puntos con\nd_near > d_vac (umbral de vacancia)"),
    ("4. Clustering\nGreedy", C_YELLOW,
     "Ordenar por d_near descendente.\nMerge r < r_cluster → un cluster = una vacancia"),
]

for i, (title, color, desc) in enumerate(steps):
    bx = 0.4 + (i % 2) * 6.45
    by = 1.3 + (i // 2) * 2.3
    add_rect(s, bx, by, 6.0, 2.0, C_ACCENT)
    add_rect(s, bx, by, 1.5, 2.0, color)
    add_text(s, title, bx, by + 0.5, 1.5, 1.2,
             font_size=14, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)
    add_text(s, desc, bx + 1.6, by + 0.3, 4.2, 1.4,
             font_size=14, color=C_WHITE)

# Output info
add_rect(s, 0.4, 6.0, 12.5, 1.2, C_MID_BG)
output_lines = [
    ("Salida por vacancia:", True, C_GREEN),
    ("  posición (x, y, z)  ·  d_near_max  ·  n_pts (puntos de grilla en cluster)  →  vacancies_output.csv", False, C_LIGHT_GY),
]
add_multiline(s, output_lines, 0.6, 6.05, 12.0, 1.1, font_size=13, color=C_WHITE)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 10 — ARQUITECTURA DEL SOFTWARE
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
slide_title(s, "Arquitectura del Software", "Módulos C++17 — Diseño modular header-only")

modules = [
    # (name, category, color, description)
    ("main.cpp",         "Orquestador",    C_HIGHLIGHT, "CLI, arg parsing,\npipeline completo,\nescritura de salida"),
    ("LammpsDumpReader", "I/O",            C_ORANGE,    "Lectura frame-by-frame\no bulk de archivos dump"),
    ("AtomData",         "Datos",          C_GREEN,     "Atom, SimBox, Frame\n(estructuras base)"),
    ("CellList",         "Algoritmos",     C_YELLOW,    "Índice espacial O(N)\n27-cell stencil PBC"),
    ("RadialBasis",      "SOAP",           C_HIGHLIGHT, "Bases radiales\nφₙ(r) gaussianas"),
    ("SphericalHarmonics","SOAP",          C_ORANGE,    "Yₗᵐ(θ,φ) vía\nrecurrencia Bonnet"),
    ("SOAPDescriptor",   "SOAP",           C_GREEN,     "Power spectrum,\nnormalización, OpenMP"),
    ("Statistics",       "Estadística",    C_YELLOW,    "mean, euclidean\nchiProbability, fit"),
    ("DefectClassifier", "Clasificación",  C_HIGHLIGHT, "buildReference,\nclassify, findVacancies"),
    ("PCA",              "Análisis",       C_ORANGE,    "fit (covarianza+eigen)\ntransform (proyección)"),
    ("WignerSeitz",      "Validación",     C_GREEN,     "Clasificación alternativa\npor sitios de referencia"),
    ("Scripts Python",   "Visualización",  C_YELLOW,    "graficar_output.py\ngraficar_pca.py\nvisualizer.py"),
]

cols = 4
rows = 3
for i, (name, cat, color, desc) in enumerate(modules):
    c = i % cols
    r = i // cols
    bx = 0.35 + c * 3.22
    by = 1.3 + r * 1.95
    add_rect(s, bx, by, 3.0, 1.75, C_ACCENT)
    add_rect(s, bx, by, 3.0, 0.35, color)
    add_text(s, f"{name}  [{cat}]", bx + 0.08, by + 0.04, 2.84, 0.32,
             font_size=11, bold=True, color=C_DARK_BG)
    add_text(s, desc, bx + 0.08, by + 0.42, 2.84, 1.2,
             font_size=11, color=C_WHITE)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 11 — CONFIGURACIÓN Y CLI
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
slide_title(s, "Configuración y CLI", "Todos los parámetros son argumentos de línea de comandos")

# Example command
add_rect(s, 0.4, 1.2, 12.5, 0.75, C_DARK_BG)
add_rect(s, 0.4, 1.2, 12.5, 0.08, C_GREEN)
cmd = "./distool --n-max 9 --l-max 9 --r-cut 5.5 --threshold 0.15 --grid-spacing 0.2 --vac-cluster-radius 3.0 --pca 2 --hist  pristine_Fe_300K.dump  damaged_Fe_10keV.dump"
add_text(s, cmd, 0.6, 1.28, 12.1, 0.6,
         font_size=12, color=C_GREEN, italic=True)

# Parameter table
param_groups = [
    ("SOAP Descriptores", C_HIGHLIGHT, [
        ("--n-max N",     "Bases radiales (def: 9)"),
        ("--l-max L",     "Momento angular máx (def: 9)"),
        ("--r-cut R",     "Radio de corte [Å] (def: 5.0)"),
        ("--sigma S",     "Ancho gaussiano [Å] (def: auto)"),
    ]),
    ("Clasificación", C_ORANGE, [
        ("--threshold T", "Umbral dist. defecto (def: 0.15)"),
        ("--ref-sia F",   "DV referencia SIA"),
        ("--ref-antv F",  "DV referencia VacAdj"),
        ("--save-dv ID F","Guardar DV de átomo ID"),
    ]),
    ("Vacancias", C_GREEN, [
        ("--vac-dist D",      "Dist. punto vacante [Å]"),
        ("--grid-spacing G",  "Espaciado grilla [Å] (def: 0.5)"),
        ("--vac-cluster-r R", "Radio merge vacancia [Å]"),
    ]),
    ("Análisis / Salida", C_YELLOW, [
        ("--pca [N]",     "PCA con N componentes (def: 2)"),
        ("--hist [B]",    "Histograma distancias (def: 50)"),
        ("--ws",          "Wigner-Seitz paralelo"),
        ("--output F",    "Nombre CSV de salida"),
    ]),
]

for gi, (group, color, params) in enumerate(param_groups):
    bx = 0.4 + gi * 3.22
    add_rect(s, bx, 2.1, 3.0, 5.1, C_ACCENT)
    add_rect(s, bx, 2.1, 3.0, 0.38, color)
    add_text(s, group, bx + 0.05, 2.12, 2.9, 0.36,
             font_size=13, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)
    for pi, (flag, desc) in enumerate(params):
        py = 2.6 + pi * 1.1
        add_rect(s, bx + 0.1, py, 2.8, 0.38, C_MID_BG)
        add_text(s, flag, bx + 0.15, py + 0.02, 2.7, 0.36,
                 font_size=11, bold=True, color=color)
        add_text(s, desc, bx + 0.15, py + 0.38, 2.7, 0.65,
                 font_size=10, color=C_LIGHT_GY)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 12 — RENDIMIENTO Y ESCALABILIDAD
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
slide_title(s, "Rendimiento y Escalabilidad")

# Performance metrics
add_rect(s, 0.4, 1.2, 5.8, 2.7, C_ACCENT)
add_rect(s, 0.4, 1.2, 5.8, 0.42, C_HIGHLIGHT)
add_text(s, "Métricas de Rendimiento", 0.4, 1.2, 5.8, 0.42,
         font_size=15, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)
perf_lines = [
    ("~500 ms / 10k átomos", True, C_YELLOW),
    ("  con n_max=9, l_max=9, r_cut=5 Å + OpenMP", False, C_LIGHT_GY),
    ("", False, None),
    ("~7 GB RAM / 10k átomos", True, C_ORANGE),
    ("  para matriz de descriptores (450 componentes)", False, C_LIGHT_GY),
    ("", False, None),
    ("CellList: O(N) build, O(27ρ) query", True, C_GREEN),
    ("Algoritmo de vacancias: O(N_grid × 27ρ)", True, C_GREEN),
]
add_multiline(s, perf_lines, 0.6, 1.75, 5.4, 2.0, font_size=13, color=C_WHITE)

# Complexity
add_rect(s, 0.4, 4.1, 5.8, 3.0, C_ACCENT)
add_rect(s, 0.4, 4.1, 5.8, 0.42, C_ORANGE)
add_text(s, "Complejidades Algorítmicas", 0.4, 4.1, 5.8, 0.42,
         font_size=15, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)
comp_lines = [
    ("Construcción CellList:", True, C_ORANGE),
    ("  O(N) — lineal en número de átomos", False, C_LIGHT_GY),
    ("SOAP por átomo:", True, C_ORANGE),
    ("  O(27ρ × n_max² × l_max²) — constante", False, C_LIGHT_GY),
    ("SOAP frame completo:", True, C_ORANGE),
    ("  O(N) con paralelismo OpenMP", False, C_LIGHT_GY),
    ("Detección vacancias:", True, C_ORANGE),
    ("  O(N_grid) — depende de grid_spacing", False, C_LIGHT_GY),
]
add_multiline(s, comp_lines, 0.6, 4.6, 5.4, 2.3, font_size=13, color=C_WHITE)

# Safety / Memory checks
add_rect(s, 6.6, 1.2, 6.3, 5.9, C_MID_BG)
add_rect(s, 6.6, 1.2, 6.3, 0.42, C_GREEN)
add_text(s, "Seguridad y Robustez", 6.6, 1.2, 6.3, 0.42,
         font_size=15, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)

safety_lines = [
    ("Pre-vuelo de memoria:", True, C_GREEN),
    ("  • Lee /proc/meminfo al inicio", False, C_LIGHT_GY),
    ("  • Alerta si descriptores > 75% RAM", False, C_LIGHT_GY),
    ("  • Sugiere n_max menor para 50% RAM", False, C_LIGHT_GY),
    ("  • Previene crashes OOM", False, C_LIGHT_GY),
    ("", False, None),
    ("Validación de grilla de vacancias:", True, C_GREEN),
    ("  • Aborta si > 200M puntos de grilla", False, C_LIGHT_GY),
    ("  • Sugiere espaciado mínimo", False, C_LIGHT_GY),
    ("  • Previene agotamiento de RAM", False, C_LIGHT_GY),
    ("", False, None),
    ("Optimizaciones hot-path:", True, C_GREEN),
    ("  • 1/(2σ²) pre-calculado", False, C_LIGHT_GY),
    ("  • Recurrencia trig (sin atan2)", False, C_LIGHT_GY),
    ("  • Normalización en esfera unitaria", False, C_LIGHT_GY),
    ("  • NATIVE_MARCH para binarios nativos", False, C_LIGHT_GY),
    ("", False, None),
    ("Suite de Tests:", True, C_YELLOW),
    ("  test_statistics, test_descriptors,", False, C_LIGHT_GY),
    ("  test_classifier, test_pca,", False, C_LIGHT_GY),
    ("  test_reader, test_wigner_seitz", False, C_LIGHT_GY),
]
add_multiline(s, safety_lines, 6.8, 1.75, 5.9, 5.1, font_size=13, color=C_WHITE)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 13 — FLUJO DE TRABAJO EN PRÁCTICA
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
slide_title(s, "Flujo de Trabajo en Práctica", "3 fases de uso recomendadas")

phases = [
    (
        "FASE 1\nDescubrimiento",
        C_HIGHLIGHT,
        "--pca 2",
        [
            "Ejecutar con el flag --pca 2",
            "Visualizar DVs en 2D con graficar_pca.py",
            "Identificar clusters → tipos de defectos",
            "Seleccionar átomos representativos de cada cluster",
        ]
    ),
    (
        "FASE 2\nConstrucción\nde Referencias",
        C_ORANGE,
        "--save-dv ID archivo.dv",
        [
            "Usar --save-dv para extraer DV de átomos",
            "identificados en la fase de descubrimiento",
            "Repetir para SIA, VacAdj y otros tipos",
            "Construir librería de referencias del material",
        ]
    ),
    (
        "FASE 3\nProducción",
        C_GREEN,
        "--ref-sia sia.dv --ref-antv vac.dv",
        [
            "Clasificar con referencias calibradas",
            "Estadísticas precisas de tipos de defectos",
            "Usar --ws para validación independiente",
            "Exportar a OVITO para visualización 3D",
        ]
    ),
]

for i, (phase, color, cmd, steps) in enumerate(phases):
    bx = 0.4 + i * 4.28
    add_rect(s, bx, 1.3, 3.9, 5.7, C_ACCENT)
    add_rect(s, bx, 1.3, 3.9, 0.75, color)
    add_text(s, phase, bx + 0.1, 1.32, 3.7, 0.72,
             font_size=14, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)
    # Command
    add_rect(s, bx + 0.1, 2.15, 3.7, 0.5, C_DARK_BG)
    add_text(s, cmd, bx + 0.15, 2.18, 3.6, 0.46,
             font_size=11, italic=True, color=color)
    # Steps
    for si, step in enumerate(steps):
        add_rect(s, bx + 0.15, 2.8 + si * 0.9, 0.35, 0.35, color)
        add_text(s, str(si + 1), bx + 0.15, 2.8 + si * 0.9, 0.35, 0.35,
                 font_size=12, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)
        add_text(s, step, bx + 0.6, 2.82 + si * 0.9, 3.1, 0.7,
                 font_size=12, color=C_WHITE)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 14 — STACK TECNOLÓGICO
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
slide_title(s, "Stack Tecnológico")

techs = [
    ("C++17", "Núcleo", C_HIGHLIGHT,
     "Lenguaje principal. Requiere GCC ≥7 o Clang ≥5. Templates, constexpr, structured bindings."),
    ("Eigen3 ≥3.3", "Álgebra lineal", C_ORANGE,
     "Header-only. Operaciones matriciales para PCA. SelfAdjointEigenSolver para descomposición."),
    ("OpenMP", "Paralelismo", C_GREEN,
     "Opcional. Paraleliza el cálculo de descriptores SOAP por átomo. Flag USE_OPENMP en CMake."),
    ("CMake ≥3.14", "Build system", C_YELLOW,
     "Release/Debug configs. GNUInstallDirs. AppImage support (NATIVE_MARCH=OFF). Tests integrados."),
    ("LAMMPS dump", "Formato de datos", C_HIGHLIGHT,
     "Formato de entrada estándar en MD. Frame-by-frame o bulk. PBC periódicas soportadas."),
    ("Python + Plotly", "Visualización", C_ORANGE,
     "Scripts: graficar_output.py, graficar_pca.py, publication_charts.py, visualizer.py."),
    ("OVITO", "Visualización 3D", C_GREEN,
     "Archivo analyzed.dump compatible con OVITO para inspección visual 3D de defectos."),
    ("Linux Desktop", "Integración", C_YELLOW,
     ".desktop + ícono .png. FreeDesktop spec. Soporte AppImage para distribución portable."),
]

for i, (name, cat, color, desc) in enumerate(techs):
    c = i % 4
    r = i // 4
    bx = 0.35 + c * 3.25
    by = 1.3 + r * 2.7
    add_rect(s, bx, by, 3.05, 2.4, C_ACCENT)
    add_rect(s, bx, by, 3.05, 0.08, color)
    add_text(s, name, bx + 0.1, by + 0.12, 2.85, 0.45,
             font_size=16, bold=True, color=color)
    add_text(s, f"[{cat}]", bx + 0.1, by + 0.55, 2.85, 0.3,
             font_size=11, italic=True, color=C_LIGHT_GY)
    add_text(s, desc, bx + 0.1, by + 0.85, 2.85, 1.4,
             font_size=11, color=C_WHITE)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 15 — CONCLUSIONES
# ══════════════════════════════════════════════════════════════════════════════
s = prs.slides.add_slide(BLANK)
slide_bg(s)
slide_title(s, "Conclusiones y Contribuciones Clave")

add_rect(s, 0.4, 1.2, 5.8, 5.8, C_ACCENT)
add_rect(s, 0.4, 1.2, 5.8, 0.45, C_HIGHLIGHT)
add_text(s, "Innovación Científica", 0.4, 1.2, 5.8, 0.45,
         font_size=16, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)

sci_lines = [
    ("✓ Robustez térmica", True, C_GREEN),
    ("  Los descriptores SOAP son estables a\n  temperaturas elevadas — sin posiciones\n  de referencia de red necesarias", False, C_LIGHT_GY),
    ("", False, None),
    ("✓ Modelo probabilístico fundamentado", True, C_GREEN),
    ("  Distribución chi → umbral estadístico\n  con base matemática sólida", False, C_LIGHT_GY),
    ("", False, None),
    ("✓ Detección geométrica de vacancias", True, C_GREEN),
    ("  Grid sampling → detecta voids\n  extendidos y regiones crowdion", False, C_LIGHT_GY),
    ("", False, None),
    ("✓ Validación cruzada", True, C_GREEN),
    ("  Wigner-Seitz opcional como\n  método independiente", False, C_LIGHT_GY),
]
add_multiline(s, sci_lines, 0.6, 1.75, 5.4, 5.0, font_size=13, color=C_WHITE)

add_rect(s, 6.6, 1.2, 6.3, 5.8, C_MID_BG)
add_rect(s, 6.6, 1.2, 6.3, 0.45, C_ORANGE)
add_text(s, "Calidad del Software", 6.6, 1.2, 6.3, 0.45,
         font_size=16, bold=True, color=C_DARK_BG, align=PP_ALIGN.CENTER)

sw_lines = [
    ("✓ Arquitectura modular", True, C_ORANGE),
    ("  Headers separados por responsabilidad,\n  sin dependencias circulares", False, C_LIGHT_GY),
    ("", False, None),
    ("✓ Alto rendimiento", True, C_ORANGE),
    ("  O(N) via CellList + OpenMP\n  Pre-computos en hot path", False, C_LIGHT_GY),
    ("", False, None),
    ("✓ Manejo de memoria", True, C_ORANGE),
    ("  Pre-flight check /proc/meminfo\n  Validación tamaño de grilla", False, C_LIGHT_GY),
    ("", False, None),
    ("✓ Suite de tests completa", True, C_ORANGE),
    ("  6 archivos de tests unitarios\n  Integrados en CMake", False, C_LIGHT_GY),
    ("", False, None),
    ("✓ Reproducibilidad", True, C_ORANGE),
    ("  CLI-only: cada análisis es\n  100% reproducible y scriptable", False, C_LIGHT_GY),
]
add_multiline(s, sw_lines, 6.8, 1.75, 5.9, 5.0, font_size=13, color=C_WHITE)


# ══════════════════════════════════════════════════════════════════════════════
#  SAVE
# ══════════════════════════════════════════════════════════════════════════════
out_path = "/home/user/DistributionTool/DistributionTool_Presentacion.pptx"
prs.save(out_path)
print(f"Saved: {out_path}")
