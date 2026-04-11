#include "DefectClassifier.h"
#include "LammpsDumpReader.h"
#include "PCA.h"
#include "SOAPDescriptor.h"
#include "Statistics.h"

#include <chrono>
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
        "  --threshold T    Distance threshold for defects (default: 0.15)\n"
        "  --vac-dist  D    Vacancy detection distance [Å] (default: r_cut×0.4)\n"
        "  --grid-spacing G Vacancy grid spacing [Å]       (default: 0.5)\n\n"
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
    const std::vector<std::array<double,3>>& vacs,
    const std::string& path)
{
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write: " + path);
    f << "# x y z\n" << std::fixed << std::setprecision(8);
    for (const auto& v : vacs)
        f << v[0] << ' ' << v[1] << ' ' << v[2] << '\n';
    std::cout << "  → vacancy CSV written: " << path
              << "  (" << vacs.size() << " sites)\n";
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

// ─────────────────────────────────────────────────────────────────────────────
// Summary
// ─────────────────────────────────────────────────────────────────────────────

static void printSummary(const Frame& frame,
                         const std::vector<std::array<double,3>>& vac_pts,
                         double grid_spacing) {
    int cnt[5] = {0,0,0,0,0};
    for (const auto& a : frame.atoms)
        ++cnt[static_cast<int>(a.defect_type)];

    const double vac_vol = vac_pts.size() * grid_spacing * grid_spacing * grid_spacing;

    std::cout << "\n┌─── Defect Summary ───────────────────────────────┐\n"
              << "│  Total atoms        : " << std::setw(6) << frame.size() << "                     │\n"
              << "│  Lattice            : " << std::setw(6) << cnt[0]        << "                     │\n"
              << "│  Interstitials      : " << std::setw(6) << cnt[1]        << "                     │\n"
              << "│  Vacancy-adjacent   : " << std::setw(6) << cnt[2]        << "                     │\n"
              << "│  Type-A defects     : " << std::setw(6) << cnt[3]        << "                     │\n"
              << "│  Unknown (distorted): " << std::setw(6) << cnt[4]        << "                     │\n"
              << "│  Vacant grid pts    : " << std::setw(6) << vac_pts.size()<< "                     │\n"
              << "│  Void volume (est.) : " << std::setw(6) << std::fixed
                                            << std::setprecision(1) << vac_vol
                                            << " Å³                 │\n"
              << "└──────────────────────────────────────────────────┘\n";

    // Rough Frenkel pair estimate: interstitials only (vacancies via grid)
    if (cnt[1] > 0)
        std::cout << "  Interstitial atoms: " << cnt[1] << '\n';
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

    double threshold     = 0.15;
    double vac_dist      = -1.0;   // -1 = auto (r_cut × 0.4)
    double grid_spacing  = 0.5;    // Å — vacancy detection grid spacing
    std::string out_file = "output.csv";
    bool   do_pca  = false;
    int    pca_nc  = 2;
    bool   do_hist = false;
    int    hist_bins = 50;
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
        else if (a == "--vac-dist")     vac_dist     = nextDbl();
        else if (a == "--grid-spacing") grid_spacing = nextDbl();
        else if (a == "--output")    out_file     = nextStr();
        else if (a == "--pca") {
            do_pca = true;
            if (nextIsUInt()) pca_nc = nextInt();
        }
        else if (a == "--hist") {
            do_hist = true;
            if (nextIsUInt()) hist_bins = nextInt();
        }
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

    // ── Print configuration ───────────────────────────────────────────────────
    std::cout << "Configuration:\n"
              << "  SOAP  n_max=" << soap.n_max
              << "  l_max="       << soap.l_max
              << "  r_cut="       << soap.r_cut << " Å"
              << "  DV size="     << soap.dvSize() << "\n"
              << "  Threshold=" << threshold
              << "  grid_spacing=" << grid_spacing << " Å\n\n";

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

        // Memory estimate warning
        {
            double dv_gb = static_cast<double>(ref_frame.size()) *
                           soap.dvSize() * sizeof(double) / 1e9;
            if (dv_gb > 1.0)
                std::cout << "  [!] Estimated DV memory: "
                          << std::fixed << std::setprecision(1) << dv_gb
                          << " GB — consider smaller --n-max/--l-max if memory is tight\n";
        }

        std::cout << "[2/4] Computing SOAP descriptors for reference frame…\n";
        auto t0 = std::chrono::steady_clock::now();
        desc.computeAll(ref_frame);
        std::cout << "      Done in " << std::fixed << std::setprecision(2)
                  << timerSec(t0) << " s\n";

        // Collect reference DVs and build statistics
        std::vector<std::vector<double>> ref_dvs;
        ref_dvs.reserve(ref_frame.size());
        for (const auto& atom : ref_frame.atoms)
            ref_dvs.push_back(atom.dv);

        clf.buildReference(ref_dvs);
        const auto& ref = clf.reference();
        std::cout << "      q̄(T) built.  <d>=" << ref.mean_dist
                  << "  var(d)=" << ref.var_dist << '\n';

        // ══════════════════════════════════════════════════════════════════════
        // 2. DAMAGED FRAME
        // ══════════════════════════════════════════════════════════════════════
        std::cout << "\n[3/4] Reading damaged frame: " << dmg_file << '\n';
        LammpsDumpReader dmg_reader(dmg_file);
        Frame dmg_frame = dmg_reader.readNext();
        std::cout << "      " << dmg_frame.size() << " atoms, timestep "
                  << dmg_frame.timestep << '\n';

        std::cout << "      Computing SOAP descriptors for damaged frame…\n";
        t0 = std::chrono::steady_clock::now();
        desc.computeAll(dmg_frame);
        std::cout << "      Done in " << std::fixed << std::setprecision(2)
                  << timerSec(t0) << " s\n";

        // ══════════════════════════════════════════════════════════════════════
        // 3. CLASSIFY
        // ══════════════════════════════════════════════════════════════════════
        std::cout << "\n[4/4] Classifying defects…\n";
        clf.classify(dmg_frame);

        // Vacancy detection via sampling grid (FaVaD §2.3.2)
        if (vac_dist < 0.0)
            vac_dist = soap.r_cut * 0.4;   // default: ~40 % of cutoff radius
        std::cout << "      Vacancy grid: " << std::fixed << std::setprecision(2)
                  << grid_spacing << " Å spacing,  threshold " << vac_dist << " Å\n";
        auto vacancies = clf.findVacanciesGrid(dmg_frame, grid_spacing, vac_dist);

        printSummary(dmg_frame, vacancies, grid_spacing);

        // ══════════════════════════════════════════════════════════════════════
        // 4. OUTPUTS
        // ══════════════════════════════════════════════════════════════════════
        std::cout << "\nWriting outputs…\n";
        writeAtomCSV(dmg_frame, out_file);

        // Vacancy file
        if (!vacancies.empty()) {
            std::string vac_file = "vacancies_" + out_file;
            writeVacancyCSV(vacancies, vac_file);
        }

        // Distance histogram (Fig. 4b equivalent)
        if (do_hist) {
            std::vector<double> dists;
            dists.reserve(dmg_frame.size());
            for (const auto& a : dmg_frame.atoms)
                dists.push_back(a.dist_to_ref);
            writeHistogram(dists, hist_bins, "hist_" + out_file);
        }

        // PCA (Fig. 7 equivalent)
        if (do_pca) {
            std::cout << "  Running PCA (" << pca_nc << " components)…\n";

            // Fit on reference DVs so PCA axes reflect the pristine crystal
            PCA pca;
            pca.fit(ref_dvs);

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

            // Project damaged frame DVs
            std::vector<std::vector<double>> dmg_dvs;
            dmg_dvs.reserve(dmg_frame.size());
            for (const auto& a : dmg_frame.atoms) dmg_dvs.push_back(a.dv);

            auto proj = pca.transform(dmg_dvs, pca_nc);
            writePCAcsv(proj, dmg_frame, "pca_" + out_file);
        }

    } catch (const std::exception& e) {
        std::cerr << "\nFatal error: " << e.what() << '\n';
        return 1;
    }

    std::cout << "\nDone.\n";
    return 0;
}
