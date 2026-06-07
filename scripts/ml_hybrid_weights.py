#!/usr/bin/env python3
"""
ML prototype: LEARN the hybrid consensus weights instead of hand-tuning them.

Motivation (README §3.8 / the ML discussion): the `robust` preset weights
(transit=1.5, frenkel=0.3, ...) and the 0.5 accept cutoff are set by hand. The
hybrid consensus score is literally a linear combination of the 8 signals passed
through a threshold — i.e. an untrained logistic regression. Here we FIT it.

Well-posed supervised task (regime contrast, WS is ground truth on relaxed frames):
  - STABLE defect   (y=1): WS-empty sites from RELAXED frames
                            (5keV final survivors + FeCrNi-HEA finalCool WS sites).
  - TRANSIENT cand. (y=0): WS-empty sites from the 5keV BALLISTIC peak (t=5000),
                            which are ~98% transients (only ~6 of 424 survive).

Honest caveats (printed at runtime):
  * The ballistic candidates are weak-labeled as "transient" en masse (~2% are
    real survivors). Proper per-candidate survival labels need temporal tracking
    across frames — exactly what the Clementina energy sweep will provide.
  * Direct position-matching of peak vacancies to the 6 final survivors does NOT
    work (at 2 A only 1/6 matches): the surviving defects localize only by
    t~15000, so survival is not predictable from a single peak frame.

No sklearn dependency: logistic regression is implemented in numpy (standardized
features, full-batch gradient descent, L2). 5-fold CV AUC.

Usage:
    python3 scripts/ml_hybrid_weights.py
Outputs:
    figures/ml_hybrid_weights.png  (+ .pdf)
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))

# CSV columns: x y z consensus ws soft_ws vor dens soap topo transit_pen frenkel_pen tf recomb ref_idx acc
SIGNALS   = ["ws", "soft_ws", "vor", "dens", "soap", "topo", "transit_pen", "frenkel_pen"]
SIG_COLS  = [4, 5, 6, 7, 8, 9, 10, 11]
REF_IDX_COL = 14
# hand-tuned robust preset weights (penalties enter the score NEGATIVELY)
ROBUST_W  = {"ws":1.0, "soft_ws":0.8, "vor":0.6, "dens":0.6, "soap":0.5,
             "topo":0.7, "transit_pen":-1.5, "frenkel_pen":-0.3}


def load(path):
    rows = [[float(x) for x in ln.split()]
            for ln in open(path) if ln.strip() and not ln.startswith("#")]
    return np.array(rows)


def ws_signals(path):
    """Return signal matrix for WS-site candidates (ref_site_idx >= 0)."""
    a = load(path)
    a = a[a[:, REF_IDX_COL] >= 0]
    return a[:, SIG_COLS]


def sigmoid(z):
    return 1.0 / (1.0 + np.exp(-z))


def fit_logreg(X, y, l2=6.0, lr=0.3, iters=6000):
    n, d = X.shape
    w = np.zeros(d); b = 0.0
    # class balancing weights
    pw = np.where(y == 1, n / (2 * max(1, y.sum())), n / (2 * max(1, (y == 0).sum())))
    for _ in range(iters):
        p = sigmoid(X @ w + b)
        g = (p - y) * pw
        w -= lr * (X.T @ g / n + l2 * w / n)
        b -= lr * g.mean()
    return w, b


def auc(y, s):
    pos = s[y == 1]; neg = s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    return (pos[:, None] > neg[None, :]).mean() + 0.5 * (pos[:, None] == neg[None, :]).mean()


def main():
    f_ball = os.path.join(ROOT, "results/5kev_20_05/hybrid_vacancies_ball_robust.csv")
    f_fin  = os.path.join(ROOT, "results/5kev_20_05/hybrid_vacancies_final_robust.csv")
    f_hea  = os.path.join(ROOT, "results/hea_finalCool/hybrid_vacancies_robust.csv")

    transient = ws_signals(f_ball)                       # y=0 (ballistic peak)
    stable_5  = ws_signals(f_fin)                        # y=1 (5keV survivors)
    stable_h  = ws_signals(f_hea)                        # y=1 (HEA relaxed WS sites)

    # subsample the large HEA-relaxed pool so it doesn't dominate
    rng = np.random.default_rng(0)
    if len(stable_h) > 400:
        stable_h = stable_h[rng.choice(len(stable_h), 400, replace=False)]
    stable = np.vstack([stable_5, stable_h])

    X = np.vstack([transient, stable])
    y = np.concatenate([np.zeros(len(transient)), np.ones(len(stable))])
    print(f"dataset: {len(transient)} transient (ballistic WS)  +  "
          f"{len(stable)} stable (relaxed WS: {len(stable_5)} 5keV + {len(stable_h)} HEA)")

    # standardize
    mu, sd = X.mean(0), X.std(0) + 1e-9
    Xs = (X - mu) / sd

    # 5-fold CV AUC
    idx = rng.permutation(len(y)); folds = np.array_split(idx, 5)
    aucs = []
    for k in range(5):
        te = folds[k]; tr = np.concatenate([folds[j] for j in range(5) if j != k])
        w, b = fit_logreg(Xs[tr], y[tr])
        aucs.append(auc(y[te], Xs[te] @ w + b))
    print(f"5-fold CV AUC = {np.nanmean(aucs):.3f} ± {np.nanstd(aucs):.3f}")
    print("  [caveat] AUC is inflated by regime leakage: ballistic vs relaxed frames\n"
          "          differ globally, so the model partly detects WHICH REGIME, not\n"
          "          per-defect stability. Real survival labels need temporal tracking\n"
          "          across frames (the Clementina energy sweep). This is a proof of concept.")

    # final fit on all data
    w, b = fit_logreg(Xs, y)
    # express learned weights in the score's sign convention (positive => stable)
    learned = dict(zip(SIGNALS, w))
    print("\nlearned standardized weights (positive => predicts STABLE defect):")
    for s in sorted(SIGNALS, key=lambda s: -abs(learned[s])):
        print(f"  {s:12s} {learned[s]:+.3f}    (robust hand: {ROBUST_W[s]:+.1f})")

    # ---- figure ----
    figdir = os.path.join(ROOT, "figures"); os.makedirs(figdir, exist_ok=True)
    fig, (axW, axR) = plt.subplots(1, 2, figsize=(11.5, 4.6),
                                   gridspec_kw={"width_ratios": [1.4, 1.0]})

    order = sorted(range(len(SIGNALS)), key=lambda i: learned[SIGNALS[i]])
    names = [SIGNALS[i] for i in order]
    lw = [learned[n] for n in names]
    hw = [ROBUST_W[n] for n in names]
    yb = np.arange(len(names))
    axW.barh(yb - 0.2, lw, 0.4, color="#1f6fb2", label="learned (logistic)")
    axW.barh(yb + 0.2, hw, 0.4, color="#c0392b", alpha=0.7, label="robust hand-tuned")
    axW.axvline(0, color="0.4", lw=0.8)
    axW.set_yticks(yb); axW.set_yticklabels(names, fontsize=9)
    axW.set_xlabel("weight  (positive ⇒ stable defect)")
    axW.set_title("Learned vs hand-tuned hybrid weights")
    axW.legend(frameon=False, fontsize=9, loc="lower right")
    axW.grid(axis="x", alpha=0.25)

    # ROC
    s_all = Xs @ w + b
    thr = np.sort(np.unique(s_all))[::-1]
    P = (y == 1).sum(); N = (y == 0).sum()
    tpr = [ ((y == 1) & (s_all >= t)).sum() / P for t in thr ]
    fpr = [ ((y == 0) & (s_all >= t)).sum() / N for t in thr ]
    axR.plot([0, 1], [0, 1], "--", color="0.6", lw=1)
    axR.plot(fpr, tpr, color="#1f6fb2", lw=2)
    axR.set_xlabel("false positive rate"); axR.set_ylabel("true positive rate")
    axR.set_title(f"Survivor/stable vs transient\n(in-sample AUC = {auc(y, s_all):.3f})")
    axR.grid(alpha=0.25)
    axR.set_aspect("equal")

    fig.tight_layout()
    out = os.path.join(figdir, "ml_hybrid_weights.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
    print("\nwrote", os.path.relpath(out, ROOT))


if __name__ == "__main__":
    main()
