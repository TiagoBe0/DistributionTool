"""
Genera el ícono PNG 256×256 para DistributionTool.
Diseño: red BCC vista en perspectiva isométrica con un átomo "defecto" destacado.
"""
from PIL import Image, ImageDraw
import math

SIZE   = 256
MARGIN = 18

img  = Image.new("RGBA", (SIZE, SIZE), (0, 0, 0, 0))
draw = ImageDraw.Draw(img)

# ── Fondo con gradiente radial oscuro ────────────────────────────────────────
bg = Image.new("RGBA", (SIZE, SIZE), (0, 0, 0, 0))
for r in range(SIZE // 2, 0, -1):
    t   = r / (SIZE / 2)
    val = int(14 + t * 10)
    col = (val, val, int(val * 1.4), 255)
    ImageDraw.Draw(bg).ellipse(
        [SIZE//2 - r, SIZE//2 - r, SIZE//2 + r, SIZE//2 + r],
        fill=col
    )
img = Image.alpha_composite(img, bg)
draw = ImageDraw.Draw(img)

# Círculo de borde exterior
draw.ellipse([2, 2, SIZE-3, SIZE-3],
             outline=(60, 120, 200, 180), width=3)

# ── Proyección isométrica ────────────────────────────────────────────────────
# Convierte (ix, iy, iz) de una celda BCC a píxeles 2D.
CX, CY = SIZE / 2, SIZE / 2 + 4
SCALE  = 54.0

def iso(ix, iy, iz):
    x2d = (ix - iy) * math.cos(math.radians(30)) * SCALE
    y2d = (ix + iy) * math.sin(math.radians(30)) * SCALE - iz * SCALE * 0.82
    return CX + x2d, CY + y2d

# ── Nodos de la celda: 8 esquinas + 1 central (BCC) ─────────────────────────
corners = [(i, j, k) for i in range(2) for j in range(2) for k in range(2)]
center_bcc = [(0.5, 0.5, 0.5)]

# Aristas de la celda unitaria (12 aristas del cubo)
edges = []
for a in corners:
    for b in corners:
        diff = sum(abs(a[d]-b[d]) for d in range(3))
        if diff == 1:
            edges.append((a, b))

# Dibuja aristas
for a, b in edges:
    pa = iso(*a)
    pb = iso(*b)
    draw.line([pa, pb], fill=(50, 90, 160, 120), width=1)

# Radios de los átomos (en píxeles)
R_CORNER = 14
R_CENTER = 18
R_DEFECT = 17

# ── Función para dibujar un átomo con efecto 3D ──────────────────────────────
def draw_atom(draw, cx, cy, r, base_color, highlight=True):
    # Sombra
    draw.ellipse([cx-r+2, cy-r+2, cx+r+2, cy+r+2],
                 fill=(0, 0, 0, 80))
    # Cuerpo principal
    draw.ellipse([cx-r, cy-r, cx+r, cy+r],
                 fill=base_color + (255,))
    # Brillo especular
    if highlight:
        hr = max(3, r // 3)
        hx = cx - r // 3
        hy = cy - r // 3
        draw.ellipse([hx-hr, hy-hr, hx+hr, hy+hr],
                     fill=(255, 255, 255, 140))

# Átomos de esquina — verde (Lattice)
for node in sorted(corners, key=lambda n: n[0]+n[1]-n[2]):
    px, py = iso(*node)
    draw_atom(draw, px, py, R_CORNER, (58, 160, 58))

# Átomo BCC central — amarillo (TypeA / defecto)
px, py = iso(0.5, 0.5, 0.5)
draw_atom(draw, px, py, R_CENTER, (210, 170, 20))

# Átomo intersticial desplazado — rojo (Interstitial)
px, py = iso(1.0, 0.5, 1.0)
# pequeño desplazamiento para simular intersticial
px += 7; py -= 9
draw_atom(draw, px, py, R_DEFECT, (200, 50, 40))

# Átomo vacancia-adyacente — azul (VacancyAdj)
# esquina (0,0,0) marcada diferente
px, py = iso(0, 0, 0)
draw_atom(draw, px, py, R_CORNER + 3, (40, 100, 210))

# ── Leyenda pequeña de colores en la esquina inferior derecha ────────────────
legend = [
    ((58, 160, 58),   "Lat"),
    ((40, 100, 210),  "VacAdj"),
    ((200, 50, 40),   "Int"),
    ((210, 170, 20),  "TypeA"),
]
lx, ly = SIZE - 58, SIZE - 72
for color, _ in legend:
    draw.ellipse([lx-4, ly-4, lx+4, ly+4], fill=color + (200,))
    ly += 16

# ── Texto pequeño "DisTool" en la parte inferior ─────────────────────────────
try:
    from PIL import ImageFont
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 14)
    draw.text((SIZE//2, SIZE - 14), "DisTool",
              font=font, fill=(200, 220, 255, 200), anchor="mm")
except Exception:
    pass  # si no hay fuente disponible, omitir texto

# ── Guardar ──────────────────────────────────────────────────────────────────
out = "packaging/distool.png"
img.save(out)
print(f"Ícono guardado: {out}  ({SIZE}×{SIZE} px, RGBA)")
