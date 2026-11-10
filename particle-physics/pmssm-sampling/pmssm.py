"""Learning where the pMSSM gives a physical spectrum, to sample it faster.

Most random points of the 19-parameter pMSSM are rejected by the spectrum
generator (no electroweak symmetry breaking, tachyons, charged LSP, no
convergence). A small neural network, trained on points SOFTSUSY already
judged (points.csv, from generate.py), predicts the outcome for new points, so
a scan can skip the ones likely to fail before paying for the calculation.

The network is a plain numpy multilayer perceptron (19 -> 32 -> 32 -> 1, tanh,
sigmoid output), trained with binary cross-entropy and Adam.

    python3 pmssm.py           # train + check + figures
"""
import csv
import sys
from pathlib import Path

import numpy as np

from generate import PARAMS

rng = np.random.default_rng(1)
LABELS = {"M1": "$M_1$", "M2": "$M_2$", "M3": "$M_3$", "mA": "$m_A$", "tanb": "$\\tan\\beta$", "mu": "$\\mu$",
          "At": "$A_t$", "Ab": "$A_b$", "Atau": "$A_\\tau$", "mq1L": "$m_{\\tilde Q_1}$", "mq3L": "$m_{\\tilde Q_3}$",
          "muR": "$m_{\\tilde u_R}$", "mtR": "$m_{\\tilde t_R}$", "mdR": "$m_{\\tilde d_R}$", "mbR": "$m_{\\tilde b_R}$",
          "ml1L": "$m_{\\tilde L_1}$", "mtauL": "$m_{\\tilde L_3}$", "meR": "$m_{\\tilde e_R}$", "mtauR": "$m_{\\tilde \\tau_R}$"}


def load(path="points.csv"):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    lo = np.array([p[2][0] for p in PARAMS])
    hi = np.array([p[2][1] for p in PARAMS])
    x = (np.array([[float(r[p[0]]) for p in PARAMS] for r in rows]) - lo) / (hi - lo)
    y = np.array([int(r["valid"]) for r in rows], float)
    h0 = np.array([float(r["h0"]) if r["h0"] else np.nan for r in rows])
    return x, y, h0


class MLP:
    def __init__(self, sizes):
        self.w = [rng.normal(0, np.sqrt(1 / a), (a, b)) for a, b in zip(sizes, sizes[1:])]
        self.b = [np.zeros(b) for b in sizes[1:]]

    def forward(self, x):
        acts = [x]
        for i, (w, b) in enumerate(zip(self.w, self.b)):
            z = acts[-1] @ w + b
            acts.append(np.tanh(z) if i < len(self.w) - 1 else 1 / (1 + np.exp(-z)))
        return acts

    def predict(self, x):
        return self.forward(x)[-1][:, 0]

    def fit(self, x, y, epochs=1000, lr=3e-3, l2=1e-3):
        """Full-batch Adam on the binary cross-entropy; returns the loss per epoch."""
        params = self.w + self.b
        m = [np.zeros_like(p) for p in params]
        v = [np.zeros_like(p) for p in params]
        losses = []
        for t in range(1, epochs + 1):
            acts = self.forward(x)
            p = np.clip(acts[-1][:, 0], 1e-7, 1 - 1e-7)
            losses.append(-np.mean(y * np.log(p) + (1 - y) * np.log(1 - p)))
            delta = (acts[-1][:, 0] - y)[:, None] / len(y)         # dL/dz of the sigmoid layer
            grads_w, grads_b = [], []
            for i in reversed(range(len(self.w))):
                grads_w.insert(0, acts[i].T @ delta + l2 * self.w[i])
                grads_b.insert(0, delta.sum(0))
                if i:
                    delta = (delta @ self.w[i].T) * (1 - acts[i] ** 2)
            for k, g in enumerate(grads_w + grads_b):
                m[k] = 0.9 * m[k] + 0.1 * g
                v[k] = 0.999 * v[k] + 0.001 * g**2
                params[k] -= lr * (m[k] / (1 - 0.9**t)) / (np.sqrt(v[k] / (1 - 0.999**t)) + 1e-8)
        return losses


def auc(score, y):
    """Area under the ROC curve = probability a valid point outranks an invalid one."""
    rank = np.argsort(np.argsort(score)) + 1
    n1 = y.sum()
    return (rank[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * (len(y) - n1))


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Optional house style (figstyle.py at the repository root)
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import figstyle
    figstyle.apply()
    si_unit, si_qty = figstyle.unit, figstyle.qty

    x, y, h0 = load()
    idx = rng.permutation(len(y))
    train, test = idx[: int(0.7 * len(y))], idx[int(0.7 * len(y)):]
    net = MLP([19, 32, 32, 1])
    losses = net.fit(x[train], y[train])
    score = net.predict(x[test])

    base = y[test].mean()
    a = auc(score, y[test])
    order = np.argsort(-score)
    kept = np.arange(1, len(test) + 1) / len(test)
    purity = np.cumsum(y[test][order]) / np.arange(1, len(test) + 1)
    half = purity[len(test) // 2 - 1]
    print(f"{len(y)} points, {int(y.sum())} valid ({100 * y.mean():.0f} %)")
    print(f"Test set ({len(test)} points): AUC = {a:.3f}; keeping the best-scored half, "
          f"{100 * half:.0f} % are valid against {100 * base:.0f} % at random "
          f"(x{half / base:.2f}, keeping {100 * np.cumsum(y[test][order])[len(test) // 2 - 1] / y[test].sum():.0f} % "
          f"of the valid points)")
    # Check: the network has learned something real on points it never saw
    assert a > 0.75, "network does no better than a weak guess"
    assert half / base > 1.4, "selection does not enrich valid points"
    valid_h0 = h0[y == 1]
    # Check: the MSSM bounds the light Higgs to roughly m_Z..135 GeV at tree level plus loops;
    # SOFTSUSY's two-loop masses should fall there for nearly all valid points
    inside = np.mean((valid_h0 > 90) & (valid_h0 < 150))
    print(f"m_h0 of valid points: median {np.median(valid_h0):.1f} GeV, {100 * inside:.0f} % in 90-150 GeV")
    assert inside > 0.9, "light Higgs masses off the MSSM range"

    # Figure 1: fraction of valid points among those kept, against the fraction kept
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    axes[0].plot(kept, purity, color="C0", label=f"Network (AUC $= {figstyle.num(f'{a:.2f}')}$)")
    axes[0].axhline(base, color="0.6", ls="--", label="Random sampling")
    axes[0].set(xlabel="Fraction of points kept", ylabel="Fraction of kept points that are valid",
                xlim=(0, 1), ylim=(0, 1))
    axes[0].legend()
    axes[1].semilogy(losses, color="C0")
    axes[1].set(xlabel="Epoch", ylabel="Training loss (cross-entropy)")
    fig.savefig("selection.png", dpi=200)

    # Figure 2: where valid points live, for the four most telling parameters
    bins = np.linspace(0, 1, 11)
    spread = []
    for j in range(19):
        k = np.digitize(x[:, j], bins[1:-1])
        frac = np.array([y[k == b].mean() for b in range(10)])
        spread.append((frac.max() - frac.min(), j, frac))
    fig, axes = plt.subplots(1, 4, figsize=(10, 2.8), sharey=True)
    for ax, (_, j, frac) in zip(axes, sorted(spread, reverse=True)[:4]):
        lo, hi = PARAMS[j][2]
        scale = 1 if PARAMS[j][0] == "tanb" else 1e-3                      # masses in TeV
        centres = scale * (lo + (hi - lo) * (bins[:-1] + bins[1:]) / 2)
        ax.step(centres, frac, where="mid", color="C0")
        ax.axhline(y.mean(), color="0.6", ls="--", lw=0.8)
        unit = "" if PARAMS[j][0] == "tanb" else " " + si_unit(r"\tera\electronvolt", "TeV")
        ax.set(xlabel=LABELS[PARAMS[j][0]] + unit, ylim=(0, 1))
    axes[0].set_ylabel("Valid fraction")
    fig.savefig("valid_fraction.png", dpi=200)

    # Figure 3: lightest Higgs mass of the valid points
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    ax.hist(np.clip(valid_h0, 60, 160), bins=list(range(60, 162, 2)), color="C0", alpha=0.8)
    ax.axvline(125.1, color="C1", ls="--", label="Measured, " + si_qty(125.1, r"\giga\electronvolt", "GeV"))
    in_window = np.mean(np.abs(valid_h0 - 125.1) < 3)
    ax.set(xlabel="$m_{h^0}$ " + si_unit(r"\giga\electronvolt", "GeV"), ylabel="Valid points",
           title=si_qty(f"{100 * in_window:.0f}", r"\percent", r"\%") + " within " + si_qty(3, r"\giga\electronvolt", "GeV")
           + " of " + si_qty(125.1, r"\giga\electronvolt", "GeV"))
    ax.legend()
    fig.savefig("higgs_mass.png", dpi=200)
    print("Saved selection.png, valid_fraction.png, higgs_mass.png")
