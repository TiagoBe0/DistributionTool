#!/usr/bin/env python3
"""
Defect cluster-size distribution, cascade morphology and defect chemistry for the
surviving (relaxed-plateau) state of the FeCrNi PKA cascades.

Method: Wigner-Seitz against the shared equilibrated reference (1 keV step-0, the
identical alloy lattice all runs start from).  Empty sites -> vacancies (species =
reference type at that site); over-occupied sites -> interstitials (the extra
atoms, species = their own type).  Defects are clustered under PBC by a 2nd-NN
distance criterion (<= 1.05 a0).  Type map: 1=Fe 2=Ni 3=Cr 4=Fe(PKA marker).

Outputs (results/physics/):
  defect_clusters_chem.csv          per-energy counts, clustering, morphology
  defect_chemistry.png/pdf          species enrichment of vacancies/interstitials
  defect_clustersize.png/pdf        pooled cluster-size distribution
Usage: python3 scripts/analyze_defect_clusters_chem.py
"""
import os
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/run/user/1000/gvfs/sftp:host=192.168.220.25/home/santi/pka_acero"
OUT = os.path.join(BASE, "results", "physics")
os.makedirs(OUT, exist_ok=True)

REF = os.path.join(BASE, "1kev", "dump.ballistic_1kev.FeCrNi_1kev.0")
ELEM = {1: "Fe", 2: "Ni", 3: "Cr", 4: "Fe(PKA)"}
CL_FACTOR = 1.05   # cluster cutoff in units of a0 (<= 2nd NN)

# label -> (dir, ballistic prefix).  Plateau frame = max step available.
RUNS = [
    ("1kev",      "1kev",      "dump.ballistic_1kev.FeCrNi_1kev"),
    ("2kev",      "2kev",      "dump.ballistic_2kev.FeCrNi_2kev_s"),
    ("3kev",      "3kev",      "dump.ballistic_3kev.FeCrNi_3kev"),
    ("4kev",      "4kev",      "dump.ballistic.FeCrNi_4keV"),
    ("5kev",      "5kev",      "dump.ballistic_5kev.FeCrNi_5kev_v2"),
    ("6kev",      "6kev",      "dump.ballistic.FeCrNi_6keV"),
    ("7kev",      "7kev",      "dump.ballistic.FeCrNi_7keV"),
    ("8kev_real", "8kev_real", "dump.ballistic.FeCrNi_8keV"),
]


def parse_dump(path):
    with open(path) as f:
        lines = f.readlines()
    i = 0
    box_lo = np.zeros(3); box_hi = np.zeros(3); cols = None; n = 0; start = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("ITEM: NUMBER OF ATOMS"):
            n = int(lines[i + 1]); i += 2; continue
        if ln.startswith("ITEM: BOX BOUNDS"):
            for d in range(3):
                lo, hi = map(float, lines[i + 1 + d].split()[:2])
                box_lo[d] = lo; box_hi[d] = hi
            i += 4; continue
        if ln.startswith("ITEM: ATOMS"):
            cols = ln.split()[2:]; start = i + 1; break
        i += 1
    ci = {c: k for k, c in enumerate(cols)}
    ix, it = ci["id"], ci["type"]
    cx, cy, cz = ci["x"], ci["y"], ci["z"]
    data = np.empty((n, 5))
    for k in range(n):
        p = lines[start + k].split()
        data[k] = (float(p[ix]), float(p[it]),
                   float(p[cx]), float(p[cy]), float(p[cz]))
    order = np.argsort(data[:, 0])
    data = data[order]
    L = box_hi - box_lo
    return data[:, 1].astype(int), data[:, 2:5], box_lo, L


def wrap(pos, lo, L):
    w = np.mod(pos - lo, L)
    return np.clip(w, 0.0, L - 1e-9)


def pbc_clusters(points, lo, L, cutoff):
    """Connected components under PBC; returns list of cluster sizes."""
    nP = len(points)
    if nP == 0:
        return []
    tree = cKDTree(wrap(points, lo, L), boxsize=L)
    pairs = tree.query_pairs(cutoff, output_type="ndarray")
    parent = np.arange(nP)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    roots = np.array([find(a) for a in range(nP)])
    _, sizes = np.unique(roots, return_counts=True)
    return sorted(sizes.tolist(), reverse=True)


def min_image_unwrap(points, lo, L):
    """Unwrap a localized cloud relative to its first point (single cascade)."""
    if len(points) == 0:
        return points
    ref = points[0]
    d = points - ref
    d -= L * np.round(d / L)
    return ref + d


def main():
    rtype, rpos, lo, L = parse_dump(REF)
    a0 = float(np.mean(L)) / 69.0
    cutoff = CL_FACTOR * a0
    tree = cKDTree(wrap(rpos, lo, L), boxsize=L)
    bulk = {t: float(np.mean(rtype == t)) for t in (1, 2, 3, 4)}
    bulk_merged = {"Fe": bulk[1] + bulk[4], "Ni": bulk[2], "Cr": bulk[3]}
    print(f"a0={a0:.4f} Å  cluster cutoff={cutoff:.3f} Å  ({CL_FACTOR}·a0)")
    print("bulk composition:",
          {k: f"{v*100:.2f}%" for k, v in bulk_merged.items()})

    rows = []
    vac_sizes_all, int_sizes_all = [], []
    vac_type_all, int_type_all = [], []
    pka_fate = []
    for label, d, prefix in RUNS:
        ddir = os.path.join(BASE, d)
        steps = []
        for fn in os.listdir(ddir):
            if fn.startswith(prefix + "."):
                s = fn.rsplit(".", 1)[-1]
                if s.isdigit():
                    steps.append(int(s))
        if not steps:
            print(f"  {label}: no frames"); continue
        step = max(steps)
        dtype, dpos, _, _ = parse_dump(os.path.join(ddir, f"{prefix}.{step}"))

        dist, site = tree.query(wrap(dpos, lo, L), k=1)
        occ = np.bincount(site, minlength=len(rpos))
        vac_sites = np.where(occ == 0)[0]

        # interstitials: extra atoms at over-occupied sites (farthest = extra)
        int_mask = np.zeros(len(dpos), dtype=bool)
        over = np.where(occ > 1)[0]
        if len(over):
            by_site = {}
            for ai, s in enumerate(site):
                if occ[s] > 1:
                    by_site.setdefault(s, []).append(ai)
            for s, atoms in by_site.items():
                atoms.sort(key=lambda a: dist[a])      # closest = resident
                for a in atoms[1:]:
                    int_mask[a] = True
        int_idx = np.where(int_mask)[0]

        n_vac, n_int = len(vac_sites), len(int_idx)
        vsz = pbc_clusters(rpos[vac_sites], lo, L, cutoff)
        isz = pbc_clusters(dpos[int_idx], lo, L, cutoff)
        vac_sizes_all += vsz; int_sizes_all += isz

        vt = rtype[vac_sites]; it_ = dtype[int_idx]
        vac_type_all += vt.tolist(); int_type_all += it_.tolist()

        # PKA fate: is the type-4 atom an interstitial?
        pka_is_int = bool(np.any(dtype[int_idx] == 4))
        pka_fate.append((label, pka_is_int))

        # morphology (single cascade -> unwrap clouds)
        def_pts = np.vstack([rpos[vac_sites], dpos[int_idx]]) if (n_vac + n_int) else np.empty((0, 3))
        sep = gyr_v = np.nan
        if n_vac and n_int:
            vc = min_image_unwrap(rpos[vac_sites], lo, L).mean(0)
            ic = min_image_unwrap(dpos[int_idx], lo, L).mean(0)
            dsep = ic - vc; dsep -= L * np.round(dsep / L)
            sep = float(np.linalg.norm(dsep))
        if n_vac:
            vu = min_image_unwrap(rpos[vac_sites], lo, L)
            gyr_v = float(np.sqrt(np.mean(np.sum((vu - vu.mean(0))**2, axis=1))))

        clustered_v = sum(s for s in vsz if s >= 2)
        rows.append(dict(
            label=label, step=step, n_vac=n_vac, n_int=n_int,
            n_vac_clusters=len(vsz), largest_vac_cluster=(vsz[0] if vsz else 0),
            clustered_vac_frac=(clustered_v / n_vac if n_vac else 0.0),
            vac_int_sep_A=sep, vac_gyration_A=gyr_v,
            pka_interstitial=pka_is_int))
        print(f"  {label}@{step}: vac={n_vac} int={n_int}  "
              f"vac clusters={vsz}  int clusters={isz}  sep={sep:.1f}Å")

    import pandas as pd
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(OUT, "defect_clusters_chem.csv"), index=False)
    print("\n" + tab.to_string(index=False))

    # ---------- chemistry figure ----------
    def comp(types):
        types = np.array(types)
        n = len(types)
        if n == 0:
            return {"Fe": 0, "Ni": 0, "Cr": 0}, 0
        fe = np.mean((types == 1) | (types == 4))
        return {"Fe": fe, "Ni": np.mean(types == 2), "Cr": np.mean(types == 3)}, n
    vc, nv = comp(vac_type_all)
    icc, ni = comp(int_type_all)
    els = ["Fe", "Ni", "Cr"]
    x = np.arange(len(els)); w = 0.27
    fig, ax = plt.subplots(figsize=(8, 5.2))
    ax.bar(x - w, [bulk_merged[e] for e in els], w, color="0.6", label="bulk alloy")
    ax.bar(x,     [vc[e] for e in els], w, color="tab:blue",
           label=f"vacancies (N={nv})")
    ax.bar(x + w, [icc[e] for e in els], w, color="tab:red",
           label=f"interstitials (N={ni})")
    ax.set_xticks(x); ax.set_xticklabels(els)
    ax.set_ylabel("species fraction")
    ax.set_title("Defect chemistry vs bulk composition — FeCrNi cascades\n"
                 "(surviving defects pooled over 1–8 keV)")
    ax.grid(alpha=0.3, axis="y"); ax.legend()
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"defect_chemistry.{ext}"), dpi=140)

    # ---------- cluster-size distribution ----------
    fig, ax = plt.subplots(figsize=(8, 5.2))
    mx = max(vac_sizes_all + int_sizes_all + [1])
    bins = np.arange(0.5, mx + 1.5, 1)
    ax.hist(vac_sizes_all, bins=bins, alpha=0.6, color="tab:blue",
            label=f"vacancy clusters (n={len(vac_sizes_all)})")
    ax.hist(int_sizes_all, bins=bins, alpha=0.6, color="tab:red",
            label=f"interstitial clusters (n={len(int_sizes_all)})")
    ax.set_xlabel("cluster size (defects)")
    ax.set_ylabel("count")
    ax.set_title("Surviving-defect cluster-size distribution (pooled 1–8 keV)")
    ax.grid(alpha=0.3, axis="y"); ax.legend()
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"defect_clustersize.{ext}"), dpi=140)

    npka = sum(1 for _, b in pka_fate if b)
    print(f"\nPKA ended as an interstitial in {npka}/{len(pka_fate)} cascades")
    print(f"Wrote chemistry + cluster-size figures to {OUT}")


if __name__ == "__main__":
    main()
