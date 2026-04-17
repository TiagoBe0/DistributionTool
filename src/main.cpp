#include "DefectClassifier.h"
#include "LammpsDumpReader.h"
#include "PCA.h"
#include "SOAPDescriptor.h"
#include "Statistics.h"
#include "WignerSeitz.h"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace DistTool;

// ─────────────────────────────────────────────────────────────────────────────
// Reference DV I/O
// ─────────────────────────────────────────────────────────────────────────────

// Load a descriptor vector from a plain-text file.
// Format: one float per line (or space-separated); comment lines start with '#'.
// Normalises to unit length to match our SOAP convention.
static std::vector<double> loadDV(const std::string& path, int expected_size) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open reference DV file: " + path);

    std::vector<double> dv;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        double v;
        while (ss >> v) dv.push_back(v);
    }

    if (dv.empty())
        throw std::runtime_error("Empty reference DV file: " + path);
    if (expected_size > 0 && (int)dv.size() != expected_size)
        throw std::runtime_error(
            "DV size mismatch in " + path + ": got " +
            std::to_string(dv.size()) + ", expected " +
            std::to_string(expected_size));

    // Normalise to unit length (consistent with SOAP normalisation)
    double norm2 = 0.0;
    for (double v : dv) norm2 += v * v;
    if (norm2 > 1e-20) {
        const double inv = 1.0 / std::sqrt(norm2);
        for (auto& v : dv) v *= inv;
    }
    return dv;
}

// Write the DV of a single atom (found by id) to a text file.
static void saveDV(const Frame& frame, int atom_id, const std::string& path) {
    const Atom* target = nullptr;
    for (const auto& a : frame.atoms)
        if (a.id == atom_id) { target = &a; break; }
    if (!target)
        throw std::runtime_error("--save-dv: atom id " +
                                 std::to_string(atom_id) + " not found in frame");
    if (target->dv.empty())
        throw std::runtime_error("--save-dv: DV not yet computed for atom " +
                                 std::to_string(atom_id));

    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write DV file: " + path);
    f << "# DV for atom id=" << atom_id
      << "  size=" << target->dv.size() << '\n';
    f << std::fixed << std::setprecision(10);
    for (double v : target->dv) f << v << '\n';
    std::cout << "  → DV saved: " << path
              << "  (atom " << atom_id << ", " << target->dv.size() << " components)\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// CLI helpers
// ─────────────────────────────────────────────────────────────────────────────

static void printBanner() {
    std::cout <<
        "╔══════════════════════════════════════════════════════════════╗\n"
        "║  DistributionTool — Defect Detection in Crystalline Solids  ║\n"
        "║  Based on: Domínguez-Gutiérrez & von Toussaint (2019)       ║\n"
        "╚══════════════════════════════════════════════════════════════╝\n\n";
}

static void printUsage(const char* prog) {
    std::cout <<
        "Usage:\n"
        "  " << prog << " [options] <reference.dump> <damaged.dump>\n\n"
        "Arguments:\n"
        "  reference.dump   LAMMPS dump of the pristine/thermalized crystal\n"
        "  damaged.dump     LAMMPS dump after irradiation/cascade\n\n"
        "SOAP options:\n"
        "  --n-max  N       Radial basis functions       (default: 9)\n"
        "  --l-max  L       Max angular momentum         (default: 9)\n"
        "  --r-cut  R       Cutoff radius [Å]            (default: 5.0)\n"
        "  --sigma  S       Gaussian width [Å] (-1=auto) (default: -1)\n\n"
        "Classification options:\n"
        "  --threshold T          Distance threshold for defects (default: 0.15)\n"
        "  --vac-dist  D          Vacancy detection distance [Å] (default: r_cut×0.4)\n"
        "  --grid-spacing G       Vacancy grid spacing [Å]       (default: 0.5)\n"
        "  --vac-cluster-radius R Cluster radius for vacancy merging [Å]\n"
        "                         (default: vac-dist; use FaVaD greedy algorithm)\n\n"
        "Reference DV files (secondary classification):\n"
        "  --ref-sia  FILE  DV file for self-interstitial atom\n"
        "  --ref-antv FILE  DV file for atom-next-to-vacancy\n"
        "  --ref-typea FILE DV file for type-A (PCA-discovered) defect\n"
        "  --save-dv ID FILE  Save DV of atom ID (from damaged frame) to FILE\n"
        "                     Use this to build reference DV files from known sites\n\n"
        "Wigner-Seitz comparison:\n"
        "  --ws             Run WS analysis in parallel with SOAP\n"
        "  --ws-r-cut R     WS search radius [Å] (default: same as --r-cut)\n\n"
        "Output options:\n"
        "  --output   FILE  Main output CSV file         (default: output.csv)\n"
        "  --pca      [N]   Run PCA with N components    (default: 2)\n"
        "  --hist     [B]   Distance histogram, B bins   (default: 50)\n\n"
        "Other:\n"
        "  --help           Print this message\n";
}

static double timerSec(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
}

// ─────────────────────────────────────────────────────────────────────────────
// Output writers
// ─────────────────────────────────────────────────────────────────────────────

static void writeAtomCSV(const Frame& frame, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write: " + path);

    f << "# id type x y z dist_to_ref defect_prob defect_type\n";
    f << std::fixed << std::setprecision(8);

    for (const auto& a : frame.atoms) {
        f << a.id    << ' '
          << a.type  << ' '
          << a.x     << ' '
          << a.y     << ' '
          << a.z     << ' '
          << a.dist_to_ref  << ' '
          << a.defect_prob  << ' '
          << defectName(a.defect_type) << '\n';
    }
    std::cout << "  → atom CSV written: " << path << '\n';
}

static void writePCAcsv(
    const std::vector<std::vector<double>>& proj,
    const Frame& frame,
    const std::string& path)
{
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write: " + path);

    int nc = static_cast<int>(proj.empty() ? 0 : proj[0].size());
    f << "# id type";
    for (int k = 0; k < nc; ++k) f << " pc" << (k+1);
    f << " dist_to_ref defect_type\n";
    f << std::fixed << std::setprecision(8);

    for (int i = 0; i < (int)frame.atoms.size(); ++i) {
        const auto& a = frame.atoms[i];
        f << a.id << ' ' << a.type;
        for (int k = 0; k < nc; ++k) f << ' ' << proj[i][k];
        f << ' ' << a.dist_to_ref
          << ' ' << defectName(a.defect_type) << '\n';
    }
    std::cout << "  → PCA CSV written: " << path << '\n';
}

static void writeVacancyCSV(
    const std::vector<VacancyCluster>& clusters,
    const std::string& path)
{
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write: " + path);
    f << "# x y z d_near_max n_grid_pts\n" << std::fixed << std::setprecision(8);
    for (const auto& c : clusters)
        f << c.center[0] << ' '
          << c.center[1] << ' '
          << c.center[2] << ' '
          << c.d_near_max << ' '
          << c.n_pts      << '\n';
    std::cout << "  → vacancy CSV written: " << path
              << "  (" << clusters.size() << " vacancies)\n";
}

// Write analyzed LAMMPS dump — same box/timestep as input, with extra columns
// defect_label is the integer value of DefectType (readable by OVITO as property)
static void writeLAMMPSDump(const Frame& frame, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write dump: " + path);

    f << "ITEM: TIMESTEP\n" << frame.timestep << "\n";
    f << "ITEM: NUMBER OF ATOMS\n" << frame.size() << "\n";
    f << "ITEM: BOX BOUNDS pp pp pp\n"
      << std::fixed << std::setprecision(10)
      << frame.box.xb[0] << " " << frame.box.xb[1] << "\n"
      << frame.box.yb[0] << " " << frame.box.yb[1] << "\n"
      << frame.box.zb[0] << " " << frame.box.zb[1] << "\n";
    f << "ITEM: ATOMS id type x y z dist_to_ref defect_prob defect_label\n";
    f << std::setprecision(8);
    for (const auto& a : frame.atoms) {
        f << a.id   << ' '
          << a.type << ' '
          << a.x    << ' '
          << a.y    << ' '
          << a.z    << ' '
          << a.dist_to_ref << ' '
          << a.defect_prob << ' '
          << static_cast<int>(a.defect_type) << '\n';
    }
    f.flush();
    f.close();
    if (!f.good())
        throw std::runtime_error("Write failed (disk full?): " + path);
    std::cout << "  → analyzed dump written: " << path << '\n';
}

static void writeHistogram(
    const std::vector<double>& dists,
    int nbins,
    const std::string& path)
{
    if (dists.empty()) return;

    double dmax = *std::max_element(dists.begin(), dists.end());
    double dmin = *std::min_element(dists.begin(), dists.end());
    if (dmax <= dmin) dmax = dmin + 1e-6;

    std::vector<int> counts(nbins, 0);
    const double width = (dmax - dmin) / nbins;
    for (double d : dists) {
        int bin = static_cast<int>((d - dmin) / width);
        if (bin >= nbins) bin = nbins - 1;
        ++counts[bin];
    }

    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write: " + path);
    f << "# bin_centre count\n" << std::fixed << std::setprecision(8);
    for (int b = 0; b < nbins; ++b)
        f << (dmin + (b + 0.5) * width) << ' ' << counts[b] << '\n';
    std::cout << "  → histogram CSV written: " << path << '\n';
}

// Per-atom WS CSV: id type x y z ws_ref_id ws_dist ws_occ ws_type
static void writeWSAtomCSV(
    const Frame&                    frame,
    const std::vector<WSAtomResult>& ws,
    const std::vector<WSSite>&       sites,
    const std::string&              path)
{
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write: " + path);
    f << "# id type x y z ws_ref_id ws_dist ws_occ ws_type\n"
      << std::fixed << std::setprecision(8);

    for (int i = 0; i < static_cast<int>(frame.atoms.size()); ++i) {
        const auto& a  = frame.atoms[i];
        const auto& r  = ws[i];
        const bool assigned = (r.ref_site_idx >= 0);
        int    ref_id  = assigned ? sites[r.ref_site_idx].ref_atom_id : -1;
        double ws_dist = assigned ? r.dist : -1.0;
        const char* ws_type = !assigned        ? "Unknown"
                            : (r.ws_occ == 1) ? "Lattice"
                                              : "Interstitial";
        f << a.id    << ' '
          << a.type  << ' '
          << a.x     << ' '
          << a.y     << ' '
          << a.z     << ' '
          << ref_id  << ' '
          << ws_dist << ' '
          << r.ws_occ << ' '
          << ws_type  << '\n';
    }
    std::cout << "  → WS atom CSV written: " << path << '\n';
}

// Per-site WS CSV: ref_id x y z occupancy site_type
static void writeWSSitesCSV(
    const std::vector<WSSite>& sites,
    const std::string&         path)
{
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write: " + path);
    f << "# ref_id x y z occupancy min_dist site_type\n"
      << std::fixed << std::setprecision(8);

    for (const auto& s : sites) {
        const char* stype = (s.occupancy == 0) ? "Vacancy"
                          : (s.occupancy == 1) ? "Normal"
                                               : "Interstitial_site";
        f << s.ref_atom_id << ' '
          << s.x           << ' '
          << s.y           << ' '
          << s.z           << ' '
          << s.occupancy   << ' '
          << s.min_dist    << ' '
          << stype         << '\n';
    }

    // Count anomalous sites
    int n_vac = 0, n_int_sites = 0;
    for (const auto& s : sites) {
        if (s.occupancy == 0) ++n_vac;
        else if (s.occupancy >= 2) ++n_int_sites;
    }
    std::cout << "  → WS sites CSV written: " << path
              << "  (" << n_vac << " vacancies, "
              << n_int_sites << " interstitial sites)\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// Path helpers
// ─────────────────────────────────────────────────────────────────────────────

// Prepend `prefix` to the filename portion of `path`, keeping the directory.
// E.g. prefixedPath("vacancies_", "/tmp/out.csv") → "/tmp/vacancies_out.csv"
static std::string prefixedPath(const std::string& prefix, const std::string& path) {
    auto sep = path.find_last_of("/\\");
    if (sep == std::string::npos) return prefix + path;
    return path.substr(0, sep + 1) + prefix + path.substr(sep + 1);
}

// ─────────────────────────────────────────────────────────────────────────────
// Summary
// ─────────────────────────────────────────────────────────────────────────────

static void printSummary(const Frame& frame,
                         const std::vector<VacancyPoint>& vac_pts,
                         const std::vector<VacancyCluster>& clusters,
                         double grid_spacing) {
    int cnt[5] = {0,0,0,0,0};
    for (const auto& a : frame.atoms)
        ++cnt[static_cast<int>(a.defect_type)];

    const double vac_vol = vac_pts.size() * grid_spacing * grid_spacing * grid_spacing;

    std::cout << "\n┌─── Defect Summary ───────────────────────────────┐\n"
              << "│  Total atoms        : " << std::setw(6) << frame.size()     << "                     │\n"
              << "│  Lattice            : " << std::setw(6) << cnt[0]           << "                     │\n"
              << "│  Interstitials      : " << std::setw(6) << cnt[1]           << "                     │\n"
              << "│  Vacancy-adjacent   : " << std::setw(6) << cnt[2]           << "                     │\n"
              << "│  Type-A defects     : " << std::setw(6) << cnt[3]           << "                     │\n"
              << "│  Unknown (distorted): " << std::setw(6) << cnt[4]           << "                     │\n"
              << "│  Vacancies (clustered): " << std::setw(4) << clusters.size()<< "                     │\n"
              << "│  Vacant grid pts    : " << std::setw(6) << vac_pts.size()   << "                     │\n"
              << "│  Void volume (est.) : " << std::setw(6) << std::fixed
                                            << std::setprecision(1) << vac_vol
                                            << " Å³                 │\n"
              << "└──────────────────────────────────────────────────┘\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
    printBanner();

    // ── Default parameters ────────────────────────────────────────────────────
    SOAPParams soap;
    soap.n_max     = 9;
    soap.l_max     = 9;
    soap.r_cut     = 5.0;
    soap.sigma     = -1.0;
    soap.normalize = true;

    double threshold          = 0.15;
    double vac_dist           = -1.0;   // -1 = auto (r_cut × 0.4)
    double grid_spacing       = 0.5;    // Å — vacancy detection grid spacing
    double vac_cluster_radius = -1.0;   // -1 = auto (= vac_dist after resolve)
    std::string sia_file, antv_file, typea_file;
    int         save_dv_id   = -1;
    std::string save_dv_path;
    std::string out_file = "output.csv";
    bool   do_pca  = false;
    int    pca_nc  = 2;
    bool   do_hist = false;
    int    hist_bins = 50;
    bool   do_ws      = false;
    double ws_r_cut   = -1.0;  // -1 = use soap.r_cut
    std::string ref_file, dmg_file;

    // ── Parse CLI ─────────────────────────────────────────────────────────────
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto nextStr = [&]() -> std::string {
            if (i+1 >= argc) throw std::runtime_error("Missing value for " + a);
            return argv[++i];
        };
        auto nextDbl = [&]() { return std::stod(nextStr()); };
        auto nextInt = [&]() { return std::stoi(nextStr()); };

        // Returns true only if the next argv looks like a plain positive integer
        // (used for optional numeric arguments to --pca / --hist).
        auto nextIsUInt = [&]() -> bool {
            if (i+1 >= argc) return false;
            std::string s = argv[i+1];
            return !s.empty() &&
                   std::all_of(s.begin(), s.end(), [](char c){ return std::isdigit(c); });
        };

        if      (a == "--help")      { printUsage(argv[0]); return 0; }
        else if (a == "--n-max")     soap.n_max   = nextInt();
        else if (a == "--l-max")     soap.l_max   = nextInt();
        else if (a == "--r-cut")     soap.r_cut   = nextDbl();
        else if (a == "--sigma")     soap.sigma   = nextDbl();
        else if (a == "--threshold")    threshold    = nextDbl();
        else if (a == "--vac-dist")           vac_dist           = nextDbl();
        else if (a == "--grid-spacing")       grid_spacing       = nextDbl();
        else if (a == "--vac-cluster-radius") vac_cluster_radius = nextDbl();
        else if (a == "--ref-sia")      sia_file     = nextStr();
        else if (a == "--ref-antv")     antv_file    = nextStr();
        else if (a == "--ref-typea")    typea_file   = nextStr();
        else if (a == "--save-dv") {
            save_dv_id   = nextInt();
            save_dv_path = nextStr();
        }
        else if (a == "--output")    out_file     = nextStr();
        else if (a == "--pca") {
            do_pca = true;
            if (nextIsUInt()) pca_nc = nextInt();
        }
        else if (a == "--hist") {
            do_hist = true;
            if (nextIsUInt()) hist_bins = nextInt();
        }
        else if (a == "--ws")       do_ws    = true;
        else if (a == "--ws-r-cut") ws_r_cut = nextDbl();
        else if (a[0] != '-') {
            if      (ref_file.empty()) ref_file = a;
            else if (dmg_file.empty()) dmg_file = a;
        }
        else {
            std::cerr << "Unknown option: " << a << "\n";
            printUsage(argv[0]); return 1;
        }
    }

    if (ref_file.empty() || dmg_file.empty()) {
        std::cerr << "Error: provide both <reference.dump> and <damaged.dump>\n\n";
        printUsage(argv[0]); return 1;
    }
    if (soap.n_max < 1 || soap.l_max < 1) {
        std::cerr << "Error: --n-max and --l-max must be >= 1\n";
        return 1;
    }
    if (soap.r_cut <= 0.0) {
        std::cerr << "Error: --r-cut must be > 0\n";
        return 1;
    }

    // ── Print configuration ───────────────────────────────────────────────────
    std::cout << "Configuration:\n"
              << "  SOAP  n_max=" << soap.n_max
              << "  l_max="       << soap.l_max
              << "  r_cut="       << soap.r_cut << " Å"
              << "  DV size="     << soap.dvSize() << "\n"
              << "  Threshold=" << threshold
              << "  grid_spacing=" << grid_spacing << " Å\n";
    if (do_ws) {
        const double wr = (ws_r_cut > 0.0) ? ws_r_cut : soap.r_cut;
        std::cout << "  WS r_cut=" << wr << " Å\n";
    }
    std::cout << '\n';

    SOAPDescriptor desc(soap);
    DefectClassifier clf(threshold);

    try {
        // ══════════════════════════════════════════════════════════════════════
        // 1. REFERENCE FRAME
        // ══════════════════════════════════════════════════════════════════════
        std::cout << "[1/4] Reading reference frame: " << ref_file << '\n';
        LammpsDumpReader ref_reader(ref_file);
        Frame ref_frame = ref_reader.readNext();
        std::cout << "      " << ref_frame.size() << " atoms, timestep "
                  << ref_frame.timestep << '\n';
        if (ref_reader.hasNext())
            std::cout << "  [!] Reference file has multiple frames; only the first is used.\n";

        // ── Pre-flight memory check ───────────────────────────────────────────
        {
            const double dv_gb = static_cast<double>(ref_frame.size())
                                 * soap.dvSize() * sizeof(double) / 1e9;

            // Read available RAM from /proc/meminfo (Linux only)
            long avail_mb = -1;
            {
                std::ifstream mi("/proc/meminfo");
                std::string key; long val;
                while (mi >> key >> val) {
                    if (key == "MemAvailable:") { avail_mb = val / 1024; break; }
                    mi.ignore(256, '\n');
                }
            }

            std::cout << std::fixed << std::setprecision(1);
            std::cout << "  Memory per frame: ~" << dv_gb << " GB";
            if (avail_mb > 0)
                std::cout << "  |  Available: ~" << avail_mb / 1024.0 << " GB";
            std::cout << '\n';

            // Abort early with a helpful message instead of crashing mid-run
            if (avail_mb > 0 && dv_gb > avail_mb / 1024.0 * 0.75) {
                // Suggest n_max that would fit in ~50 % of available RAM
                int suggest_n = 4;
                for (int n = 9; n >= 2; --n) {
                    double gb = static_cast<double>(ref_frame.size())
                                * (n*(n+1)/2) * (n+1) * 8.0 / 1e9;
                    if (gb < avail_mb / 1024.0 * 0.5) { suggest_n = n; break; }
                }
                throw std::runtime_error(
                    "Not enough RAM for SOAP descriptors.\n"
                    "  Needed per frame : ~" + [&]{ std::ostringstream s;
                        s << std::fixed << std::setprecision(1) << dv_gb;
                        return s.str(); }() + " GB\n"
                    "  Available        : ~" + [&]{ std::ostringstream s;
                        s << std::fixed << std::setprecision(1) << avail_mb/1024.0;
                        return s.str(); }() + " GB\n"
                    "  Try: --n-max " + std::to_string(suggest_n) +
                         " --l-max " + std::to_string(suggest_n));
            }
        }

        std::cout << "[2/4] Computing SOAP descriptors for reference frame…\n";
        auto t0 = std::chrono::steady_clock::now();
        desc.computeAll(ref_frame);
        std::cout << "      Done in " << std::fixed << std::setprecision(2)
                  << timerSec(t0) << " s\n";

        clf.buildReference(ref_frame.atoms);

        // Build WS spatial index before freeing reference positions.
        WignerSeitz ws;
        if (do_ws) {
            if (ws_r_cut < 0.0) ws_r_cut = soap.r_cut;
            ws.build(ref_frame, ws_r_cut);
            std::cout << "      WS index built  (r_cut=" << std::fixed
                      << std::setprecision(2) << ws_r_cut << " Å, "
                      << ws.sites().size() << " reference sites)\n";
        }

        // Liberar DVs del frame de referencia: buildReference ya extrajo todo lo
        // necesario (media DV).  Esto libera ~N×DV×8 bytes antes de cargar el
        // frame dañado, reduciendo el pico de RAM a un solo frame en lugar de dos.
        for (auto& a : ref_frame.atoms)
            std::vector<double>().swap(a.dv);

        const auto& ref = clf.reference();
        std::cout << "      q̄(T) built.  <d>="
                  << std::fixed << std::setprecision(4) << ref.mean_dist
                  << "  σ_d=" << std::sqrt(ref.var_dist)
                  << "  χ(k=" << std::setprecision(1) << ref.k_chi
                  << ", σ=" << std::setprecision(4) << ref.sigma_chi << ")\n";

        // Load optional per-defect reference DVs for secondary classification
        {
            const int dvsz = soap.dvSize();
            std::vector<double> dv_sia, dv_antv, dv_typea;
            if (!sia_file.empty()) {
                dv_sia = loadDV(sia_file, dvsz);
                std::cout << "      SIA reference loaded:   " << sia_file << '\n';
            }
            if (!antv_file.empty()) {
                dv_antv = loadDV(antv_file, dvsz);
                std::cout << "      ANtV reference loaded:  " << antv_file << '\n';
            }
            if (!typea_file.empty()) {
                dv_typea = loadDV(typea_file, dvsz);
                std::cout << "      TypeA reference loaded: " << typea_file << '\n';
            }
            if (!dv_sia.empty() || !dv_antv.empty() || !dv_typea.empty())
                clf.setDefectReferences(dv_sia, dv_antv, dv_typea);
        }

        // ══════════════════════════════════════════════════════════════════════
        // 2. DAMAGED FRAME
        // ══════════════════════════════════════════════════════════════════════
        std::cout << "\n[3/4] Reading damaged frame: " << dmg_file << '\n';
        LammpsDumpReader dmg_reader(dmg_file);
        Frame dmg_frame = dmg_reader.readNext();
        std::cout << "      " << dmg_frame.size() << " atoms, timestep "
                  << dmg_frame.timestep << '\n';
        if (dmg_reader.hasNext())
            std::cout << "  [!] Damaged file has multiple frames; only the first is used.\n";

        std::cout << "      Computing SOAP descriptors for damaged frame…\n";
        t0 = std::chrono::steady_clock::now();
        desc.computeAll(dmg_frame);
        std::cout << "      Done in " << std::fixed << std::setprecision(2)
                  << timerSec(t0) << " s\n";

        // Optional: save the DV of a specific atom for use as a reference later
        if (save_dv_id >= 0)
            saveDV(dmg_frame, save_dv_id, save_dv_path);

        // ══════════════════════════════════════════════════════════════════════
        // 3. CLASSIFY
        // ══════════════════════════════════════════════════════════════════════
        std::cout << "\n[4/4] Classifying defects…\n";
        clf.classify(dmg_frame);

        // Liberar DVs del frame dañado: classify ya extrajo todo lo necesario.
        // Los outputs (CSV, dump) sólo usan defect_type/dist_to_ref/defect_prob.
        if (save_dv_id < 0)  // conservar si el usuario pidió --save-dv
            for (auto& a : dmg_frame.atoms)
                std::vector<double>().swap(a.dv);

        // Vacancy detection via sampling grid (FaVaD §2.3.2)
        if (vac_dist < 0.0)
            vac_dist = soap.r_cut * 0.4;   // default: ~40 % of cutoff radius
        if (vac_cluster_radius < 0.0)
            vac_cluster_radius = vac_dist; // default: same radius as detection threshold

        std::cout << "      Vacancy grid: " << std::fixed << std::setprecision(2)
                  << grid_spacing << " Å spacing,  threshold " << vac_dist << " Å"
                  << ",  cluster radius " << vac_cluster_radius << " Å\n";

        auto vac_pts  = clf.findVacanciesGrid(dmg_frame, grid_spacing, vac_dist);
        auto clusters = clf.clusterVacancyPoints(vac_pts, vac_cluster_radius, dmg_frame.box);

        // Post-processing: mark Unknown atoms adjacent to vacancy cluster centers
        // as VacancyAdjacent.  These are atoms whose local environment is distorted
        // because a first-shell neighbour is missing.  We use r_cut * 0.65 as the
        // adjacency radius — enough to capture the first coordination shell in
        // typical metals (e.g. BCC Fe/W/Ni at the default 5 Å cutoff gives ~3.25 Å,
        // which is larger than any BCC/FCC nearest-neighbour distance).
        if (!clusters.empty()) {
            const double adj_r  = soap.r_cut * 0.65;
            const double adj_r2 = adj_r * adj_r;
            const double Lx = dmg_frame.box.lx();
            const double Ly = dmg_frame.box.ly();
            const double Lz = dmg_frame.box.lz();
            int n_adj = 0;
            for (auto& atom : dmg_frame.atoms) {
                if (atom.defect_type != DefectType::Unknown) continue;
                for (const auto& c : clusters) {
                    double dx = atom.x - c.center[0];
                    double dy = atom.y - c.center[1];
                    double dz = atom.z - c.center[2];
                    if (dmg_frame.box.periodic[0]) dx -= Lx * std::round(dx / Lx);
                    if (dmg_frame.box.periodic[1]) dy -= Ly * std::round(dy / Ly);
                    if (dmg_frame.box.periodic[2]) dz -= Lz * std::round(dz / Lz);
                    if (dx*dx + dy*dy + dz*dz <= adj_r2) {
                        atom.defect_type = DefectType::VacancyAdjacent;
                        ++n_adj;
                        break;
                    }
                }
            }
            if (n_adj > 0)
                std::cout << "      " << n_adj
                          << " Unknown atoms re-classified as VacancyAdjacent"
                          << " (adj_r=" << std::fixed << std::setprecision(2)
                          << adj_r << " Å)\n";
        }

        // Wigner-Seitz classification (parallel to SOAP)
        std::vector<WSAtomResult> ws_results;
        if (do_ws) {
            std::cout << "\n[WS]  Classifying via Wigner-Seitz…\n";
            ws_results = ws.classify(dmg_frame);
            std::cout << "      Vacancies: " << ws.vacancyCount()
                      << "  Interstitials: " << ws.interstitialCount() << '\n';
        }

        printSummary(dmg_frame, vac_pts, clusters, grid_spacing);

        // ══════════════════════════════════════════════════════════════════════
        // 4. OUTPUTS
        // ══════════════════════════════════════════════════════════════════════
        std::cout << "\nWriting outputs…\n";
        writeAtomCSV(dmg_frame, out_file);
        writeLAMMPSDump(dmg_frame, prefixedPath("analyzed_", dmg_file));

        // Vacancy file: one row per cluster (physical vacancy)
        if (!clusters.empty()) {
            writeVacancyCSV(clusters, prefixedPath("vacancies_", out_file));
        }

        // Wigner-Seitz outputs
        if (do_ws && !ws_results.empty()) {
            writeWSAtomCSV(dmg_frame, ws_results, ws.sites(),
                           prefixedPath("ws_", out_file));
            writeWSSitesCSV(ws.sites(),
                            prefixedPath("ws_sites_", out_file));
        }

        // Distance histogram (Fig. 4b equivalent)
        if (do_hist) {
            std::vector<double> dists;
            dists.reserve(dmg_frame.size());
            for (const auto& a : dmg_frame.atoms)
                dists.push_back(a.dist_to_ref);
            writeHistogram(dists, hist_bins, prefixedPath("hist_", out_file));
        }

        // PCA (Fig. 7 equivalent)
        if (do_pca) {
            std::cout << "  Running PCA (" << pca_nc << " components)…\n";

            // Fit on reference atoms so PCA axes reflect the pristine crystal
            PCA pca;
            pca.fit(ref_frame.atoms);

            const auto& evr = pca.explainedVarianceRatio();
            std::cout << "  Explained variance:";
            double cumvar = 0.0;
            for (int k = 0; k < std::min(pca_nc, (int)evr.size()); ++k) {
                cumvar += evr[k];
                std::cout << "  PC" << (k+1) << "="
                          << std::fixed << std::setprecision(1)
                          << evr[k]*100 << "%";
            }
            std::cout << "  (cumulative " << cumvar*100 << "%)\n";

            auto proj = pca.transform(dmg_frame.atoms, pca_nc);
            writePCAcsv(proj, dmg_frame, prefixedPath("pca_", out_file));
        }

    } catch (const std::exception& e) {
        std::cerr << "\nFatal error: " << e.what() << '\n';
        return 1;
    }

    std::cout << "\nDone.\n";
    return 0;
}
