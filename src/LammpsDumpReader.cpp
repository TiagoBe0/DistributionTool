#include "LammpsDumpReader.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace DistTool {

LammpsDumpReader::LammpsDumpReader(const std::string& path)
    : path_(path), file_(path)
{
    if (!file_.is_open())
        throw std::runtime_error("Cannot open dump file: " + path);
}

bool LammpsDumpReader::hasNext() const {
    // peek() attempts to read the next byte without consuming it.
    // Returns EOF only when the stream is truly exhausted.
    return file_.peek() != std::char_traits<char>::eof();
}

Frame LammpsDumpReader::readNext() {
    if (!hasNext())
        throw std::runtime_error("No more frames in: " + path_);
    return parseFrame();
}

std::vector<Frame> LammpsDumpReader::readAll() {
    std::vector<Frame> frames;
    while (hasNext()) {
        try { frames.push_back(parseFrame()); }
        catch (const std::runtime_error& e) {
            // EOF reached cleanly between frames — normal termination.
            std::string msg = e.what();
            if (msg.find("EOF") != std::string::npos) break;
            throw;  // Re-throw real parse errors.
        }
    }
    return frames;
}

// ── Internal helpers ──────────────────────────────────────────────────────────

std::vector<std::string> LammpsDumpReader::splitColumns(const std::string& header) {
    // header is like  "ITEM: ATOMS id type x y z"
    std::vector<std::string> cols;
    std::istringstream ss(header);
    std::string tok;
    ss >> tok >> tok;           // skip "ITEM:" and "ATOMS"
    while (ss >> tok) cols.push_back(tok);
    return cols;
}

Frame LammpsDumpReader::parseFrame() {
    Frame frame;
    std::string line;

    // ── TIMESTEP ─────────────────────────────────────────────────────────────
    while (std::getline(file_, line)) {
        if (line.find("ITEM: TIMESTEP") != std::string::npos) break;
    }
    if (file_.eof()) throw std::runtime_error("EOF before TIMESTEP");

    std::getline(file_, line);
    frame.timestep = std::stoi(line);

    // ── NUMBER OF ATOMS ───────────────────────────────────────────────────────
    std::getline(file_, line);   // "ITEM: NUMBER OF ATOMS"
    std::getline(file_, line);
    int n_atoms = std::stoi(line);

    // ── BOX BOUNDS ────────────────────────────────────────────────────────────
    std::getline(file_, line);   // "ITEM: BOX BOUNDS ..."

    // Detect box format from the header tokens:
    //   "abc origin"  → 4 numbers/row: ax ay az ox | bx by bz oy | cx cy cz oz
    //                    (orthorhombic: ay=az=bx=bz=cx=cy=0; xlo=ox, xhi=ox+ax)
    //   "xy xz yz"    → 3 numbers/row: xlo_b xhi_b tilt
    //   default       → 2 numbers/row: xlo xhi
    bool box_abc    = (line.find("abc")    != std::string::npos);
    bool box_tilt   = (line.find("xy xz") != std::string::npos ||
                       line.find("xy")    != std::string::npos);

    // Parse periodic flags from the last 3 tokens of the BOX BOUNDS line.
    // LAMMPS writes them as "pp"/"p" (periodic) or "fs"/"ss"/"f"/"s" (fixed/shrink).
    // Example: "ITEM: BOX BOUNDS pp pp pp"  or  "ITEM: BOX BOUNDS xy xz yz pp pp pp"
    {
        std::vector<std::string> toks;
        std::istringstream ss(line);
        std::string token;
        while (ss >> token) toks.push_back(token);
        if (toks.size() >= 3) {
            for (int d = 0; d < 3; ++d) {
                const auto& t = toks[toks.size() - 3 + d];
                frame.box.periodic[d] = (t == "pp" || t == "p");
            }
        }
    }

    // Strip trailing \r (Windows line endings)
    auto stripCR = [](std::string& s) {
        if (!s.empty() && s.back() == '\r') s.pop_back();
    };

    if (box_abc) {
        // "abc origin" format:
        //   Row 0: ax  ay  az  ox     →  xlo = ox, xhi = ox + ax
        //   Row 1: bx  by  bz  oy     →  ylo = oy, yhi = oy + by
        //   Row 2: cx  cy  cz  oz     →  zlo = oz, zhi = oz + cz
        // (For non-orthorhombic cells ay/az/bx/bz/cx/cy ≠ 0; we use the
        //  diagonal elements and the origin to build the bounding box.)
        auto readABC = [&](std::array<double,2>& b, int diag_col) {
            std::getline(file_, line);
            stripCR(line);
            std::istringstream ss(line);
            std::vector<double> v;
            double x;
            while (ss >> x) v.push_back(x);
            // v[0..2] = lattice vector, v[3] = origin component
            double L_diag = (v.size() > (size_t)diag_col) ? v[diag_col] : 1.0;
            double origin = (v.size() > 3)                 ? v[3]         : 0.0;
            b[0] = origin;
            b[1] = origin + L_diag;
        };
        readABC(frame.box.xb, 0);
        readABC(frame.box.yb, 1);
        readABC(frame.box.zb, 2);
    } else if (box_tilt) {
        // "xy xz yz" triclinic: "xlo_b xhi_b xy" per row
        // For neighbour-list purposes we use the bounding-box hi/lo.
        auto readTilt = [&](std::array<double,2>& b) {
            std::getline(file_, line);
            stripCR(line);
            std::istringstream ss(line);
            ss >> b[0] >> b[1];  // tilt factor ignored (PBC still OK for MIC)
        };
        readTilt(frame.box.xb);
        readTilt(frame.box.yb);
        readTilt(frame.box.zb);
    } else {
        // Standard orthogonal: "xlo xhi"
        auto readBound = [&](std::array<double,2>& b) {
            std::getline(file_, line);
            stripCR(line);
            std::istringstream ss(line);
            ss >> b[0] >> b[1];
        };
        readBound(frame.box.xb);
        readBound(frame.box.yb);
        readBound(frame.box.zb);
    }

    // ── ATOMS ─────────────────────────────────────────────────────────────────
    std::getline(file_, line);   // "ITEM: ATOMS id type x y z ..."
    auto cols = splitColumns(line);

    // Map column names to indices
    auto colIdx = [&](std::initializer_list<const char*> names) -> int {
        for (auto* name : names) {
            auto it = std::find(cols.begin(), cols.end(), std::string(name));
            if (it != cols.end()) return static_cast<int>(it - cols.begin());
        }
        return -1;
    };

    int id_c   = colIdx({"id"});
    int type_c = colIdx({"type"});
    int x_c    = colIdx({"x","xu","xs"});
    int y_c    = colIdx({"y","yu","ys"});
    int z_c    = colIdx({"z","zu","zs"});

    // Detect scaled coordinates (xs/ys/zs are in [0,1] and must be unscaled).
    bool x_scaled = (x_c >= 0 && cols[x_c] == "xs");
    bool y_scaled = (y_c >= 0 && cols[y_c] == "ys");
    bool z_scaled = (z_c >= 0 && cols[z_c] == "zs");

    if (x_scaled || y_scaled || z_scaled)
        std::cout << "      [!] Scaled coordinates detected ("
                  << (x_scaled ? "xs " : "") << (y_scaled ? "ys " : "")
                  << (z_scaled ? "zs" : "")
                  << ") — converting to Å.\n";

    const double xlo = frame.box.xb[0], lx = frame.box.lx();
    const double ylo = frame.box.yb[0], ly = frame.box.ly();
    const double zlo = frame.box.zb[0], lz = frame.box.lz();

    frame.atoms.resize(n_atoms);
    for (int i = 0; i < n_atoms; ++i) {
        std::getline(file_, line);
        std::istringstream ss(line);
        std::vector<std::string> vals;
        std::string v;
        while (ss >> v) vals.push_back(v);

        Atom& a = frame.atoms[i];
        if (id_c   >= 0 && id_c   < (int)vals.size()) a.id   = std::stoi(vals[id_c]);
        if (type_c >= 0 && type_c < (int)vals.size()) a.type = std::stoi(vals[type_c]);
        if (x_c    >= 0 && x_c    < (int)vals.size()) a.x    = std::stod(vals[x_c]);
        if (y_c    >= 0 && y_c    < (int)vals.size()) a.y    = std::stod(vals[y_c]);
        if (z_c    >= 0 && z_c    < (int)vals.size()) a.z    = std::stod(vals[z_c]);

        if (x_scaled) a.x = xlo + a.x * lx;
        if (y_scaled) a.y = ylo + a.y * ly;
        if (z_scaled) a.z = zlo + a.z * lz;
    }

    // Sort by atom id for deterministic ordering
    std::sort(frame.atoms.begin(), frame.atoms.end(),
              [](const Atom& a, const Atom& b){ return a.id < b.id; });

    return frame;
}

} // namespace DistTool
