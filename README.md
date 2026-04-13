# DistributionTool

A C++17 command-line tool for semi-automatic identification and classification of defects in crystalline solids from LAMMPS molecular dynamics simulations. Based on the SOAP descriptor approach of Domínguez-Gutiérrez & von Toussaint and the FaVaD workflow (von Toussaint et al., 2020).

---

## Table of Contents

1. [Physical background](#1-physical-background)
2. [Method overview](#2-method-overview)
3. [Mathematical formulation](#3-mathematical-formulation)
   - 3.1 [SOAP descriptor](#31-soap-descriptor)
   - 3.2 [Radial basis functions](#32-radial-basis-functions)
   - 3.3 [Spherical harmonics](#33-spherical-harmonics)
   - 3.4 [Distance metric and probability model](#34-distance-metric-and-probability-model)
   - 3.5 [Defect classification](#35-defect-classification)
   - 3.6 [Vacancy detection](#36-vacancy-detection)
   - 3.7 [Principal Component Analysis](#37-principal-component-analysis)
4. [Code architecture](#4-code-architecture)
5. [Build instructions](#5-build-instructions)
6. [Usage](#6-usage)
7. [Output files](#7-output-files)
8. [Typical workflow](#8-typical-workflow)
9. [Parameter guide](#9-parameter-guide)
10. [References](#10-references)

---

## 1. Physical background

Energetic particle impacts (neutrons, ions) create collision cascades in crystalline materials, leaving behind point defects: **vacancies** (missing atoms), **self-interstitial atoms** (SIAs), and more complex defect clusters. Quantifying this damage is critical for fusion and fission materials research.

Classical approaches (Wigner-Seitz cell, Voronoi tessellation) work well at 0 K but fail at elevated temperatures where thermal motion blurs the lattice. This tool uses a descriptor-vector (fingerprint) approach that is robust to thermal noise and does not require the pristine lattice positions at classification time.

---

## 2. Method overview

The analysis pipeline has four stages:

```
LAMMPS dump (pristine)  ──►  Compute SOAP DVs  ──►  Build reference q̄(T)
                                                           │
LAMMPS dump (damaged)   ──►  Compute SOAP DVs  ──►  Classify atoms
                                                           │
                                                    Detect vacancies (grid)
                                                           │
                                                    Cluster vacancy points
                                                           │
                                                    PCA (optional)
```

1. **Reference frame**: compute descriptor vectors (DVs) for all atoms in the pristine/thermalized crystal; average them to obtain the reference fingerprint $\bar{q}(T)$.
2. **Damaged frame**: compute DVs for the irradiated sample; compare each atom to $\bar{q}(T)$.
3. **Classification**: threshold on distance $d^i$; secondary classification by nearest known defect reference.
4. **Vacancy detection**: uniform sampling grid; points farther than a threshold from any atom are vacancy/void regions.
5. **Vacancy clustering**: greedy algorithm merges nearby grid points into individual physical vacancy events.
6. **PCA** (optional): project DVs to 2–3 dimensions; clusters reveal unknown defect geometries.

---

## 3. Mathematical formulation

### 3.1 SOAP descriptor

The **Smooth Overlap of Atomic Positions** (SOAP) descriptor (Bartók et al., Phys. Rev. B 87, 2013) maps the local atomic environment of atom $i$ onto a rotationally-invariant fingerprint vector.

**Atomic density field**

The local environment of atom $i$ is encoded as a sum of Gaussians centered on neighbours $j$:

$$\rho^i(\mathbf{r}') = \sum_{j \neq i} f_\text{cut}(r_{ij})\, \exp\!\left(-\frac{|\mathbf{r}' - \mathbf{r}_{ij}|^2}{2\sigma^2}\right)$$

where $\sigma$ is the Gaussian width and $f_\text{cut}$ is a smooth cutoff function.

**Expansion coefficients**

The density is expanded in a basis of radial functions $\phi_n(r)$ and real spherical harmonics $Y_l^m(\hat{r})$:

$$c_{nlm}^i = \sum_{j \neq i} \phi_n(r_{ij})\, Y_l^m(\hat{r}_{ij})$$

**Power spectrum (descriptor vector)**

The rotationally-invariant power spectrum combines coefficients at fixed $l$:

$$p_{nn'l}^i = \pi\sqrt{\frac{8}{2l+1}} \sum_{m=-l}^{l} c_{nlm}^i\, c_{n'lm}^i$$

Only the upper triangle $n \leq n'$ is stored, giving a descriptor of size:

$$D = \frac{n_\text{max}(n_\text{max}+1)}{2}\,(l_\text{max}+1)$$

For the default parameters ($n_\text{max} = 9$, $l_\text{max} = 9$): $D = 450$.

**Normalisation**

The descriptor is normalised to the unit sphere:

$$\tilde{q}^i = \frac{p^i}{\|p^i\|}$$

so that distances between DVs lie in $[0, \sqrt{2}]$ and are independent of the number of neighbours.

---

### 3.2 Radial basis functions

The $n_\text{max}$ Gaussian centres are placed equidistantly in $(0,\, r_\text{cut}]$:

$$r_n = \frac{(n+1)\, r_\text{cut}}{n_\text{max}+1}, \quad n = 0, \ldots, n_\text{max}-1$$

Each basis function is:

$$\phi_n(r) = f_\text{cut}(r)\, \exp\!\left(-\frac{(r - r_n)^2}{2\sigma^2}\right)$$

with the cosine smooth cutoff:

$$f_\text{cut}(r) = \begin{cases} \dfrac{1}{2}\!\left(1 + \cos\!\dfrac{\pi r}{r_\text{cut}}\right) & r < r_\text{cut} \\ 0 & r \geq r_\text{cut} \end{cases}$$

If $\sigma$ is not specified by the user ($\sigma = -1$), it defaults to the spacing between centres: $\sigma = r_\text{cut}/(n_\text{max}+1)$.

---

### 3.3 Spherical harmonics

The real (tesseral) spherical harmonics follow the QUIP/libatoms convention:

$$Y_l^0(\hat{r}) = N_{l,0}\, P_l^0(\cos\theta)$$

$$Y_l^m(\hat{r}) = \sqrt{2}\, N_{l,m}\, P_l^m(\cos\theta)\, \cos(m\phi), \quad m > 0$$

$$Y_l^{-m}(\hat{r}) = \sqrt{2}\, N_{l,m}\, P_l^m(\cos\theta)\, \sin(m\phi), \quad m > 0$$

with normalisation:

$$N_{l,m} = \sqrt{\frac{2l+1}{4\pi}\cdot\frac{(l-m)!}{(l+m)!}}$$

Components are stored at index $l^2 + l + m$, giving $(l_\text{max}+1)^2$ values total.

**Implementation (hot path)**

A single-pass Bonnet recurrence computes all $P_l^m$ in $O(l_\text{max}^2)$ operations:

| Step | Formula |
|---|---|
| Diagonal seed | $P_0^0 = 1$ |
| Diagonal advance | $P_{m+1}^{m+1} = -(2m+1)\sin\theta\, P_m^m$ |
| Superdiagonal | $P_{m+1}^m = x(2m+1)\, P_m^m$ |
| General recurrence | $P_l^m = \dfrac{x(2l-1)\,P_{l-1}^m - (l+m-1)\,P_{l-2}^m}{l-m}$ |

The azimuthal angle is handled via a trigonometric recurrence that avoids `atan2` entirely:

$$\cos((m+1)\phi) = \cos(m\phi)\cos\phi - \sin(m\phi)\sin\phi$$

where $\cos\phi = d_x/r_{xy}$, $\sin\phi = d_y/r_{xy}$, $r_{xy} = \sqrt{d_x^2 + d_y^2}$.

This replaces $O(l_\text{max})$ transcendental calls with $O(l_\text{max})$ multiply-adds.

---

### 3.4 Distance metric and probability model

**Distance to reference**

For each atom $i$ in the damaged frame, the distance to the reference fingerprint $\bar{q}(T)$ is:

$$d^i = \|\tilde{q}^i - \bar{q}(T)\|_2$$

where $\bar{q}(T)$ is the arithmetic mean of all DVs in the pristine frame.

**Chi-distribution probability model (FaVaD, Eq. 6)**

The distribution of $d$ values in a thermalized crystal follows a **chi distribution** with $k$ degrees of freedom (corresponding to the number of active DV components):

$$P(d \mid k, \sigma) \propto d^{k-2} \exp\!\left(-\frac{d^2}{2\sigma^2}\right)$$

The distribution has a mode at:

$$d_\text{peak} = \sigma\sqrt{k-2}, \quad k > 2$$

For iron BCC at 300 K the paper reports $k \approx 7$. The tool fits $k$ and $\sigma$ automatically from the reference frame using the method of moments:

$$k = \frac{2\langle d^2 \rangle^2}{\mathrm{Var}(d^2)}, \qquad \sigma = \sqrt{\frac{\langle d^2 \rangle}{k}}$$

The probability is normalised to 1 at the mode:

$$P_\text{norm}(d) = \left(\frac{d}{d_\text{peak}}\right)^{k-2} \exp\!\left(-\frac{d^2 - d_\text{peak}^2}{2\sigma^2}\right)$$

and the defect probability output in the CSV is $1 - P_\text{norm}(d^i)$.

---

### 3.5 Defect classification

**Primary classification**

$$\text{type}^i = \begin{cases} \text{Lattice} & d^i < \tau \\ \text{(secondary)} & d^i \geq \tau \end{cases}$$

where $\tau$ is the user-supplied threshold (default 0.15).

**Secondary classification**

If per-defect reference DVs are provided (via `--ref-sia`, `--ref-antv`, `--ref-typea`), distorted atoms are assigned to the nearest reference:

$$\text{type}^i = \arg\min_{\alpha \in \{\text{SIA, ANtV, TypeA}\}} \|\tilde{q}^i - \tilde{q}^\alpha_\text{ref}\|_2$$

If no reference DVs are provided, all distorted atoms are labelled **Unknown**.

**Defect types**

| Label | Meaning |
|---|---|
| `Lattice` | Regular lattice site; $d^i < \tau$ |
| `Interstitial` | Self-interstitial atom (SIA) |
| `VacancyAdj` | Atom adjacent to a single vacancy (ANtV) |
| `TypeA` | PCA-discovered novel defect type |
| `Unknown` | Distorted but no reference DV available |

---

### 3.6 Vacancy detection

The tool implements the sampling-grid method from FaVaD (§2.3.2).

A uniform grid of $N_x \times N_y \times N_z$ points is placed inside the simulation box, with spacing $\Delta$ (default 0.5 Å):

$$N_\alpha = \left\lceil \frac{L_\alpha}{\Delta} \right\rceil, \qquad \alpha \in \{x, y, z\}$$

For each grid point $\mathbf{g}$, the distance to the nearest atom in the damaged frame is computed using the linked-cell algorithm:

$$d_\text{near}(\mathbf{g}) = \min_{j} |\mathbf{g} - \mathbf{r}_j|$$

Grid points with $d_\text{near}(\mathbf{g}) > d_\text{vac}$ are classified as vacancy/void:

$$\mathcal{V} = \{ \mathbf{g} \mid d_\text{near}(\mathbf{g}) > d_\text{vac} \}$$

The estimated void volume is $|\mathcal{V}| \cdot \Delta^3$, which converges as $\Delta \to 0$.

**Advantages over the pristine-lattice approach:**
- Does not require the reference frame for vacancy identification.
- Detects extended voids and crowdion regions, not only single vacancies.
- Volume estimate is grid-refinement consistent.

#### Vacancy clustering

Raw grid points are merged into individual physical vacancy events by a greedy algorithm (FaVaD §2.3.2):

1. Sort all vacant grid points by $d_\text{near}$ in descending order (deepest void first).
2. Take the first unassigned point as the seed of a new cluster.
3. Absorb every remaining unassigned point within $r_\text{cluster}$ (default = $d_\text{vac}$) using minimum-image PBC.
4. Repeat from step 2 until all points are assigned.

Each resulting `VacancyCluster` records:

| Field | Description |
|---|---|
| `center` | Position of the grid point with maximum $d_\text{near}$ in the cluster |
| `d_near_max` | Peak void depth $\max(d_\text{near})$ within the cluster [Å] |
| `n_pts` | Number of grid points merged into this cluster |

The cluster count equals the number of distinct physical vacancies (or void pockets) detected in the frame. The `--vac-cluster-radius` option overrides $r_\text{cluster}$ independently of `--vac-dist`.

---

### 3.7 Principal Component Analysis

PCA is applied to the DV matrix of the damaged frame (fitted on the reference DVs so that the principal axes reflect the pristine crystal).

Given the centred data matrix $\mathbf{X} \in \mathbb{R}^{N \times D}$, the covariance matrix is:

$$\mathbf{C} = \frac{\mathbf{X}^\top \mathbf{X}}{N-1}$$

Eigendecomposition (using Eigen's `SelfAdjointEigenSolver`):

$$\mathbf{C} = \mathbf{V} \boldsymbol{\Lambda} \mathbf{V}^\top$$

The projection onto the first $k$ components is:

$$\mathbf{z}^i = \mathbf{V}_k^\top (\tilde{q}^i - \bar{q})$$

Clusters of distorted atoms in the PCA plot indicate novel defect geometries (type-A defects in FaVaD terminology).

---

## 4. Code architecture

```
DistributionTool/
├── CMakeLists.txt
├── include/
│   ├── AtomData.h          # Atom, SimBox, Frame structs; DefectType enum
│   ├── CellList.h          # Linked-cell spatial index (header-only)
│   ├── RadialBasis.h       # Gaussian radial basis functions
│   ├── SphericalHarmonics.h# Real tesseral spherical harmonics
│   ├── SOAPDescriptor.h    # SOAP power spectrum + computeAll()
│   ├── DefectClassifier.h  # ReferenceSet, DefectClassifier
│   ├── Statistics.h        # mean, euclidean, chiProbability, fitChiParams
│   ├── PCA.h               # PCA fit/transform
│   └── LammpsDumpReader.h  # LAMMPS dump parser
└── src/
    ├── RadialBasis.cpp
    ├── SphericalHarmonics.cpp
    ├── SOAPDescriptor.cpp
    ├── DefectClassifier.cpp
    ├── Statistics.cpp
    ├── PCA.cpp
    ├── LammpsDumpReader.cpp
    └── main.cpp            # CLI, output writers
```

### Key data structures

**`AtomData.h`**

```cpp
struct Atom {
    int    id, type;
    double x, y, z;
    std::vector<double> dv;          // normalised SOAP descriptor
    double dist_to_ref;              // d^i = ||q̃^i − q̄(T)||
    double defect_prob;              // 1 − P_norm(d^i)
    DefectType defect_type;
};

struct SimBox {
    std::array<double,2> xb, yb, zb; // [lo, hi] per axis
    bool periodic[3];
};

struct Frame { int timestep; SimBox box; std::vector<Atom> atoms; };
```

**`CellList.h`** (header-only)

Linked-cell algorithm for $O(N)$ neighbour queries. Cells have side $\geq r_\text{cut}$; each query visits 27 cells. Supports both `neighbours(int i)` (displacement vectors for SOAP) and `nearestDist2FromPoint(x,y,z)` (for vacancy detection).

**`DefectClassifier.h`** — key types

```cpp
struct VacancyPoint {
    std::array<double,3> pos;
    double d_near;   // distance to nearest atom [Å]
};

struct VacancyCluster {
    std::array<double,3> center;     // position of grid point with max d_near
    double               d_near_max; // peak void depth in this cluster [Å]
    int                  n_pts;      // grid points merged into this cluster
};
```

`DefectClassifier` exposes four pipeline steps:

| Method | Purpose |
|---|---|
| `buildReference(dvs)` | Compute $\bar{q}(T)$ and fit chi-distribution from pristine DVs |
| `classify(frame)` | Label every atom; populate `dist_to_ref`, `defect_prob`, `defect_type` |
| `findVacanciesGrid(frame, Δ, d_vac)` | Return all empty grid points as `VacancyPoint` objects |
| `clusterVacancyPoints(pts, r, box)` | Merge grid points into `VacancyCluster` events (greedy + PBC) |

### Performance-critical paths

| Hot path | Technique |
|---|---|
| `SphericalHarmonics::computeInto` | Bonnet recurrence $O(l_\text{max}^2)$; trig recurrence (no `atan2`) |
| `RadialBasis::computeInto` | Pre-computed $1/(2\sigma^2)$; single loop |
| `SOAPDescriptor::computeAll` | `CellList` $O(N)$ neighbour search; per-thread buffers (OpenMP) |
| `CellList::nearestDist2FromPoint` | 27-cell search with PBC minimum image |
| Vacancy grid | $O(N_\text{grid} \cdot 27\rho r_\text{vac}^3)$; no heap allocation per query |
| `clusterVacancyPoints` | $O(M^2)$ in the number of vacant grid points $M$; sort + single pass |

---

## 5. Build instructions

**Dependencies**

| Library | Version | Notes |
|---|---|---|
| CMake | ≥ 3.14 | Build system |
| Eigen3 | ≥ 3.3 | Linear algebra (header-only) |
| C++17 compiler | GCC ≥ 7, Clang ≥ 5 | |
| OpenMP | optional | Parallel SOAP computation |

**Ubuntu / Debian**

```bash
sudo apt install cmake libeigen3-dev g++
```

**Build**

```bash
git clone <repo>
cd DistributionTool
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./build/distool --help
```

**Debug build** (with `-Wall -Wextra`)

```bash
cmake -B build_debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build_debug -j$(nproc)
```

---

## 6. Usage

```
./build/distool [options] <reference.dump> <damaged.dump>
```

`reference.dump` — LAMMPS dump of the pristine or thermalized crystal  
`damaged.dump`   — LAMMPS dump after irradiation or cascade

### Full option reference

| Option | Default | Description |
|---|---|---|
| `--n-max N` | 9 | Number of radial basis functions |
| `--l-max L` | 9 | Maximum angular momentum |
| `--r-cut R` | 5.0 | Cutoff radius [Å] |
| `--sigma S` | auto | Gaussian width [Å]; -1 = spacing/$1$ |
| `--threshold T` | 0.15 | Distance threshold for defect classification |
| `--vac-dist D` | r_cut×0.4 | Vacancy detection distance [Å] |
| `--grid-spacing G` | 0.5 | Vacancy grid spacing [Å] |
| `--vac-cluster-radius R` | = vac-dist | Cluster merge radius for vacancy grouping [Å] |
| `--ref-sia FILE` | — | Reference DV file for SIA atoms |
| `--ref-antv FILE` | — | Reference DV file for ANtV atoms |
| `--ref-typea FILE` | — | Reference DV file for type-A atoms |
| `--save-dv ID FILE` | — | Save DV of atom ID (damaged frame) to file |
| `--output FILE` | output.csv | Main output CSV |
| `--pca [N]` | 2 | Run PCA with N components |
| `--hist [B]` | 50 | Write distance histogram (B bins) |
| `--help` | — | Print help |

### Example — iron cascade (parameters from FaVaD)

```bash
./build/distool \
    --n-max 9 --l-max 9 --r-cut 5.5 \
    --threshold 0.15 --grid-spacing 0.2 \
    --vac-cluster-radius 3.0 \
    --pca 2 --hist \
    pristine_Fe_300K.dump damaged_Fe_10keV.dump
```

### Example — BCC tungsten

```bash
./build/distool \
    --n-max 5 --l-max 4 --r-cut 4.5 \
    --threshold 0.15 \
    --pca 2 \
    ref_W.dump damaged_W.dump
```

---

## 7. Output files

| File | Contents |
|---|---|
| `output.csv` | Per-atom: `id type x y z dist_to_ref defect_prob defect_type` |
| `vacancies_output.csv` | One row per vacancy cluster: `x y z d_near_max n_grid_pts` |
| `hist_output.csv` | Distance histogram: `bin_centre count` |
| `pca_output.csv` | PCA projection: `id type pc1 pc2 … dist_to_ref defect_type` |

**`output.csv` column descriptions**

| Column | Description |
|---|---|
| `dist_to_ref` | $d^i = \|\tilde{q}^i - \bar{q}(T)\|$ |
| `defect_prob` | $1 - P_\text{norm}(d^i \mid k, \sigma)$ from chi-distribution model |
| `defect_type` | `Lattice`, `Interstitial`, `VacancyAdj`, `TypeA`, or `Unknown` |

**`vacancies_output.csv` column descriptions**

| Column | Description |
|---|---|
| `x y z` | Position of the grid point with the largest $d_\text{near}$ in the cluster |
| `d_near_max` | Peak distance to the nearest atom within the cluster [Å] |
| `n_grid_pts` | Number of empty grid points merged into this cluster |

Each row represents one physical vacancy event (or void pocket). The total number of rows equals the vacancy count printed in the defect summary.

---

## 8. Typical workflow

### Step 1 — initial run

```bash
./build/distool ref.dump damaged.dump --pca 2 --hist
```

Inspect `pca_output.csv`: clusters of atoms far from the lattice cluster indicate unknown defect types.  
Inspect `vacancies_output.csv`: each row is one physical vacancy event; `n_grid_pts` indicates cluster size and `d_near_max` measures how deep (empty) the void is.

### Step 2 — build reference DV for a known defect site

If you know that atom 42 in `frame_with_sia.dump` is at a SIA site:

```bash
./build/distool ref.dump frame_with_sia.dump --save-dv 42 sia_dv.dat
```

### Step 3 — final classification with references

```bash
./build/distool ref.dump damaged.dump \
    --ref-sia  sia_dv.dat  \
    --ref-antv antv_dv.dat \
    --pca 2 --hist
```

Atoms that were `Unknown` will now be assigned to `Interstitial` or `VacancyAdj` according to the nearest reference DV.

### Step 4 — visualisation

The output CSVs can be loaded into OVITO, VisIt, or plotted with the provided Python scripts:

```bash
python3 graficar_output.py    # dist_to_ref distribution by defect type
python3 graficar_pca.py       # 2D PCA scatter coloured by defect type
```

---

## 9. Parameter guide

### Cutoff radius `--r-cut`

The most important parameter. Should include at least the first-nearest-neighbour shell:

| Material | Structure | $a_0$ [Å] | 1NN dist [Å] | Suggested `r-cut` |
|---|---|---|---|---|
| Fe | BCC | 2.87 | 2.48 | 4.5–5.5 |
| W  | BCC | 3.16 | 2.74 | 4.5–5.5 |
| Cu | FCC | 3.61 | 2.55 | 4.0–5.0 |
| Ni | FCC | 3.52 | 2.49 | 4.0–5.0 |

Including the 2NN shell improves discrimination between defect types at the cost of descriptor size and computation time.

### Descriptor size `--n-max` / `--l-max`

Larger values improve resolution in descriptor space but increase memory and computation:

| $n_\text{max}$ | $l_\text{max}$ | $D$ | Approx. time (per 1k atoms) |
|---|---|---|---|
| 5 | 4 | 75 | < 0.1 s |
| 9 | 9 | 450 | ~0.5 s |
| 12 | 12 | 1014 | ~2 s |

For defect detection, $n_\text{max} = 5$–$9$ and $l_\text{max} = 4$–$9$ are usually sufficient.

### Threshold `--threshold`

Separates lattice from defect atoms in $d^i$ space. Typical values:
- **0.10–0.15**: tight; reduces false positives; may miss lightly distorted atoms.
- **0.20–0.30**: loose; catches more defects; more false positives near grain boundaries or surfaces.

Calibrate by inspecting the distance histogram (`--hist`): the lattice peak is narrow near $d \approx \bar{d}_\text{ref}$; defects appear at larger $d$.

### Vacancy distance `--vac-dist`

Should be roughly half the nearest-neighbour distance. If the grid correctly recovers single vacancies but misses void interiors, reduce `--grid-spacing`. If too many grid points are flagged, increase `--vac-dist`.

---

## 10. References

1. **Domínguez-Gutiérrez & von Toussaint** (2019). *On the detection and classification of material defects in crystalline solids after energetic particle impact simulations.* Nuclear Materials and Energy, **22**, 100724. https://doi.org/10.1016/j.nme.2019.100724

2. **von Toussaint, Domínguez-Gutiérrez, Compostella & Rampp** (2020). *FaVaD: A software workflow for characterisation and visualizing of defects in crystalline structures.* arXiv:2004.08184

3. **Bartók, Kondor & Csányi** (2013). *On representing chemical environments.* Phys. Rev. B **87**, 184115. https://doi.org/10.1103/PhysRevB.87.184115

4. **Plimpton** (1995). *Fast Parallel Algorithms for Short-Range Molecular Dynamics.* J. Comp. Phys. **117**, 1–19. (LAMMPS)
# DistributionTool
