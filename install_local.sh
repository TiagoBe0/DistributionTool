#!/bin/bash
# install_local.sh — Registra DistributionTool en el escritorio local del usuario.
# Instala el ícono y el .desktop para que aparezca en el lanzador de aplicaciones.
#
# Uso:
#   ./install_local.sh                          # usa DistributionTool-x86_64.AppImage
#   ./install_local.sh <ruta/al/ejecutable>     # apunta a otro binario o AppImage
set -euo pipefail

RED='\033[0;31m'; GREEN='\033[0;32m'; CYAN='\033[0;36m'; NC='\033[0m'
log()  { echo -e "${GREEN}[✓]${NC} $*"; }
step() { echo -e "\n${CYAN}━━  $*${NC}"; }
die()  { echo -e "${RED}[✗] ERROR:${NC} $*" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ICON_SRC="${SCRIPT_DIR}/packaging/distool.png"
APPIMAGE="${SCRIPT_DIR}/DistributionTool-x86_64.AppImage"

if [[ $# -ge 1 ]]; then
    APPIMAGE="$(realpath "$1")"
fi

[[ -f "${APPIMAGE}" ]] || die "Ejecutable no encontrado: ${APPIMAGE}\n  Compilá primero con ./build_appimage.sh"
[[ -f "${ICON_SRC}" ]] || die "Ícono no encontrado: ${ICON_SRC}"

ICON_DIR="${HOME}/.local/share/icons/hicolor/256x256/apps"
APP_DIR="${HOME}/.local/share/applications"
mkdir -p "${ICON_DIR}" "${APP_DIR}"

step "Instalando ícono"
cp "${ICON_SRC}" "${ICON_DIR}/distool.png"
log "Ícono instalado en ${ICON_DIR}/distool.png"

step "Creando entrada de escritorio"
cat > "${APP_DIR}/distool.desktop" <<EOF
[Desktop Entry]
Version=1.0
Type=Application
Name=DistributionTool
GenericName=Defect Analyzer
Comment=SOAP-based defect classifier for crystalline materials (MD simulations)
Exec=${APPIMAGE} %f
Icon=distool
Terminal=false
Categories=Science;Physics;Education;
MimeType=application/x-lammps-dump;text/csv;
Keywords=SOAP;defect;crystal;MD;LAMMPS;vacancy;interstitial;nanoparticle;
StartupNotify=true
StartupWMClass=distool_viewer
EOF
log "Desktop entry instalado en ${APP_DIR}/distool.desktop"

step "Actualizando base de datos de aplicaciones"
if command -v update-desktop-database &>/dev/null; then
    update-desktop-database "${APP_DIR}"
    log "Base de datos actualizada."
else
    echo "  (update-desktop-database no disponible — puede que tarde en aparecer)"
fi

if command -v gtk-update-icon-cache &>/dev/null; then
    gtk-update-icon-cache -f -t "${HOME}/.local/share/icons/hicolor" 2>/dev/null || true
fi

echo ""
echo -e "${GREEN}╔══════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║  DistributionTool instalado localmente.              ║${NC}"
echo -e "${GREEN}║  Buscá 'DistributionTool' en tu lanzador.           ║${NC}"
echo -e "${GREEN}║                                                      ║${NC}"
printf  "${GREEN}║${NC}  Ejecutable: %-40s${GREEN}║${NC}\n" "$(basename "${APPIMAGE}")"
echo -e "${GREEN}╚══════════════════════════════════════════════════════╝${NC}"
