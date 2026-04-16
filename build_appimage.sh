#!/bin/bash
# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  build_appimage.sh — Empaqueta DistributionTool como AppImage portable      ║
# ║                                                                              ║
# ║  Uso:                                                                        ║
# ║    ./build_appimage.sh                  # build normal                       ║
# ║    ./build_appimage.sh --skip-build     # solo reempaqueta (bins ya ok)      ║
# ║    ./build_appimage.sh --clean          # borra todo y empieza de cero       ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
set -euo pipefail

# ── Colores ──────────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'; NC='\033[0m'
log()  { echo -e "${GREEN}[✓]${NC} $*"; }
step() { echo -e "\n${CYAN}━━  $*${NC}"; }
warn() { echo -e "${YELLOW}[!]${NC} $*"; }
die()  { echo -e "${RED}[✗] ERROR:${NC} $*" >&2; exit 1; }

# ── Directorios ───────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_DIR="${SCRIPT_DIR}/packaging/tools"
BUILD_CLI="${SCRIPT_DIR}/_build_appimage_cli"   # distool     (pipeline)
BUILD_GUI="${SCRIPT_DIR}/_build_appimage_gui"   # distool_viewer (GUI)
APPDIR="${SCRIPT_DIR}/_AppDir"
OUTPUT="${SCRIPT_DIR}/DistributionTool-x86_64.AppImage"

# ── Argumentos ────────────────────────────────────────────────────────────────
SKIP_BUILD=0
DO_CLEAN=0
for arg in "$@"; do
    case "$arg" in
        --skip-build) SKIP_BUILD=1 ;;
        --clean)      DO_CLEAN=1   ;;
        --help|-h)
            grep '^# ║' "$0" | sed 's/^# ║  //; s/  *║$//'
            exit 0 ;;
        *) die "Argumento desconocido: $arg" ;;
    esac
done

# ── Limpieza opcional ─────────────────────────────────────────────────────────
if [[ $DO_CLEAN -eq 1 ]]; then
    warn "Borrando builds anteriores..."
    rm -rf "${BUILD_CLI}" "${BUILD_GUI}" "${APPDIR}" "${OUTPUT}" \
           "${TOOLS_DIR}/linuxdeploy"* "${TOOLS_DIR}/appimagetool"*
    log "Limpieza completa."
fi

# ─────────────────────────────────────────────────────────────────────────────
# PASO 1 — Descargar herramientas de empaquetado si no existen
# ─────────────────────────────────────────────────────────────────────────────
step "PASO 1/5  Herramientas de empaquetado"

mkdir -p "${TOOLS_DIR}"

LINUXDEPLOY="${TOOLS_DIR}/linuxdeploy-x86_64.AppImage"
APPIMAGETOOL="${TOOLS_DIR}/appimagetool-x86_64.AppImage"

download_tool() {
    local url="$1" dest="$2" name="$3"
    if [[ -f "${dest}" ]]; then
        log "${name} ya descargado."
    else
        warn "Descargando ${name} (~10 MB)…"
        curl -L --progress-bar "${url}" -o "${dest}"
        chmod +x "${dest}"
        log "${name} listo."
    fi
}

download_tool \
    "https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage" \
    "${LINUXDEPLOY}" "linuxdeploy"

download_tool \
    "https://github.com/AppImage/AppImageKit/releases/download/continuous/appimagetool-x86_64.AppImage" \
    "${APPIMAGETOOL}" "appimagetool"

# Las herramientas son AppImages; si FUSE no está disponible necesitan esta variable.
export APPIMAGE_EXTRACT_AND_RUN=1

# ─────────────────────────────────────────────────────────────────────────────
# PASO 2 — Compilar distool (pipeline CLI)
# ─────────────────────────────────────────────────────────────────────────────
step "PASO 2/5  Compilar distool (pipeline)"

if [[ $SKIP_BUILD -eq 0 ]]; then
    # -DNATIVE_MARCH=OFF es OBLIGATORIO para que el binario corra en otras máquinas
    cmake -B "${BUILD_CLI}" \
          -DCMAKE_BUILD_TYPE=Release \
          -DNATIVE_MARCH=OFF \
          -DCMAKE_INSTALL_PREFIX="${APPDIR}/usr" \
          "${SCRIPT_DIR}" \
          2>&1 | grep -v "^--" || true

    cmake --build "${BUILD_CLI}" \
          -j"$(nproc)" \
          --target distool
    log "distool compilado."
else
    warn "--skip-build activo, omitiendo compilación de distool."
    [[ -f "${BUILD_CLI}/distool" ]] || die "No existe ${BUILD_CLI}/distool — compilá antes."
fi

# ─────────────────────────────────────────────────────────────────────────────
# PASO 3 — Compilar distool_viewer (GUI)
# ─────────────────────────────────────────────────────────────────────────────
step "PASO 3/5  Compilar distool_viewer (GUI)"

if [[ $SKIP_BUILD -eq 0 ]]; then
    cmake -B "${BUILD_GUI}" \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX="${APPDIR}/usr" \
          "${SCRIPT_DIR}/viewer" \
          2>&1 | grep -v "^--" || true

    cmake --build "${BUILD_GUI}" \
          -j"$(nproc)" \
          --target distool_viewer
    log "distool_viewer compilado."
else
    warn "--skip-build activo, omitiendo compilación del viewer."
    [[ -f "${BUILD_GUI}/distool_viewer" ]] || die "No existe ${BUILD_GUI}/distool_viewer."
fi

# ─────────────────────────────────────────────────────────────────────────────
# PASO 4 — Armar el AppDir
# ─────────────────────────────────────────────────────────────────────────────
step "PASO 4/5  Armando AppDir"

rm -rf "${APPDIR}"
mkdir -p "${APPDIR}/usr/bin" "${APPDIR}/usr/lib"

# Instalar binarios + .desktop + ícono via CMake
cmake --install "${BUILD_CLI}" --prefix "${APPDIR}/usr"
cmake --install "${BUILD_GUI}" --prefix "${APPDIR}/usr"

# Archivos requeridos en la RAÍZ del AppDir (convención AppImage)
cp "${SCRIPT_DIR}/packaging/AppRun"          "${APPDIR}/AppRun"
cp "${SCRIPT_DIR}/packaging/distool.desktop" "${APPDIR}/distool.desktop"
cp "${SCRIPT_DIR}/packaging/distool.png"     "${APPDIR}/distool.png"
chmod +x "${APPDIR}/AppRun"

# Eliminar archivos de desarrollo que cmake instala pero no son necesarios en runtime
rm -rf "${APPDIR}/usr/include"
rm -rf "${APPDIR}/usr/lib/cmake"
rm -rf "${APPDIR}/usr/lib/pkgconfig"
find  "${APPDIR}/usr/lib" -name "*.a" -delete   # static libs

log "AppDir armado:"
find "${APPDIR}" -not -path "${APPDIR}/usr/lib/*" | sort | sed 's|'"${APPDIR}"'||'

# ─────────────────────────────────────────────────────────────────────────────
# PASO 5 — Resolver dependencias y generar AppImage
# ─────────────────────────────────────────────────────────────────────────────
step "PASO 5/5  Resolviendo librerías y generando AppImage"

# linuxdeploy copia las .so necesarias en AppDir/usr/lib/.
# libGL y variantes NO se incluyen: deben venir del driver del sistema (NVIDIA/AMD/Intel).
warn "Resolviendo dependencias con linuxdeploy…"
"${LINUXDEPLOY}" \
    --appdir "${APPDIR}" \
    --executable "${APPDIR}/usr/bin/distool_viewer" \
    --executable "${APPDIR}/usr/bin/distool" \
    --exclude-library "libGL*" \
    --exclude-library "libGLX*" \
    --exclude-library "libGLdispatch*" \
    --exclude-library "libOpenGL*" \
    --exclude-library "libEGL*" \
    --exclude-library "libvulkan*"

warn "Generando AppImage con appimagetool…"
ARCH=x86_64 "${APPIMAGETOOL}" \
    --no-appstream \
    "${APPDIR}" \
    "${OUTPUT}"

# ─────────────────────────────────────────────────────────────────────────────
# Resumen
# ─────────────────────────────────────────────────────────────────────────────
SIZE_MB=$(du -sm "${OUTPUT}" | awk '{print $1}')
echo ""
echo -e "${GREEN}╔══════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║  ¡AppImage listo!                                    ║${NC}"
echo -e "${GREEN}║                                                      ║${NC}"
printf  "${GREEN}║${NC}  Archivo : %-41s${GREEN}║${NC}\n" "$(basename "${OUTPUT}")"
printf  "${GREEN}║${NC}  Tamaño  : %-41s${GREEN}║${NC}\n" "${SIZE_MB} MB"
echo -e "${GREEN}║                                                      ║${NC}"
echo -e "${GREEN}║  Para ejecutar:                                      ║${NC}"
echo -e "${GREEN}║    chmod +x DistributionTool-x86_64.AppImage         ║${NC}"
echo -e "${GREEN}║    ./DistributionTool-x86_64.AppImage                ║${NC}"
echo -e "${GREEN}║                                                      ║${NC}"
echo -e "${GREEN}║  Con un archivo:                                     ║${NC}"
echo -e "${GREEN}║    ./DistributionTool-x86_64.AppImage mi_dump.dump   ║${NC}"
echo -e "${GREEN}╚══════════════════════════════════════════════════════╝${NC}"
