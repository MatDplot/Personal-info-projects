"""Bound states of the hydrogen atom from the radial Schrodinger equation.

With u(r) = r R(r), each l gives a 1D problem
    -1/2 u'' + [l(l+1)/(2 r^2) - 1/r] u = E u,   u(0) = u(r_max) = 0,
solved on a grid with second-order finite differences (tridiagonal matrix).
The levels are checked against Bohr's E_n = -1/(2 n^2) and <r> of 1s against 3/2 a0.

    python3 hydrogen.py        # check + figures
Atomic units: hbar = m_e = e = 4 pi eps0 = 1, lengths in a0, energies in hartree.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

HARTREE_EV = 27.211386
N, R_MAX = 6000, 80.0
r = np.linspace(0, R_MAX, N + 2)[1:-1]   # interior points, u = 0 at both ends
DR = r[1] - r[0]


def levels(l, n_states=3):
    """Lowest energies and normalised u(r) for angular momentum l."""
    diag = 1 / DR**2 + l * (l + 1) / (2 * r**2) - 1 / r
    off = np.full(N - 1, -0.5 / DR**2)
    e, u = eigh_tridiagonal(diag, off, select="i", select_range=(0, n_states - 1))
    u /= np.sqrt((u**2).sum(axis=0) * DR)
    u *= np.sign(u[np.argmax(np.abs(u), axis=0), np.arange(n_states)])   # positive lobe largest
    return e, u


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Optional house style (figstyle.py at the repository root)
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import figstyle
    figstyle.apply()
    si_unit, si_qty = figstyle.unit, figstyle.qty

    # 1. Check: E depends on n = n_r + l + 1 only, and equals -1/(2 n^2)
    states = {l: levels(l) for l in range(3)}
    for l, (e, _) in states.items():
        for n_r, en in enumerate(e):
            n = n_r + l + 1
            print(f"n = {n}, l = {l}   E = {en * HARTREE_EV:8.4f} eV   Bohr {-HARTREE_EV / (2 * n * n):8.4f} eV")
            assert abs(en * 2 * n * n + 1) < 5e-4, f"E off for n = {n}, l = {l}"
    r_mean = (r * states[0][1][:, 0] ** 2).sum() * DR
    print(f"<r>_1s = {r_mean:.4f} a0 (exact 1.5)")
    assert abs(r_mean - 1.5) < 1e-3

    # 2. Radial densities P(r) = |u|^2 = r^2 |R|^2
    fig, ax = plt.subplots(figsize=(6, 4))
    for (n, l), c in zip([(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)],
                         ["k", "C0", "C1", "C2", "C4", "0.6"]):
        ax.plot(r, states[l][1][:, n - l - 1] ** 2, color=c, label=f"{n}{'spd'[l]}")
    ax.set(xlim=(0, 25), ylim=(0, None), xlabel="$r$ $[a_0]$", ylabel="$r^2 |R_{nl}(r)|^2$ $[a_0^{-1}]$")
    ax.legend()
    fig.savefig("radial_densities.png", dpi=200)

    # 3. Spectrum: computed levels against Bohr's formula
    fig, ax = plt.subplots(figsize=(6, 3.4))
    for n in range(1, 7):
        ax.axhline(-HARTREE_EV / (2 * n * n), color="0.6", lw=0.8,
                   label=(r"Bohr, $E_n = -\frac{" + (r"\qty{13.6}{\electronvolt}" if figstyle.TEX else r"13.6\,\mathrm{eV}")
                          + "}{n^2}$") if n == 1 else None)
    for l, (e, _) in states.items():
        ax.plot(np.full(len(e), l), e * HARTREE_EV, "_", ms=30, mew=2, color=["C0", "C1", "C2"][l],
                label="Finite differences" if l == 0 else None)
    ax.set(xticks=[0, 1, 2], xticklabels=["s", "p", "d"], xlim=(-0.6, 2.6), ylabel="$E$ " + si_unit(r"\electronvolt", "eV"))
    ax.legend(loc="center right", markerscale=0.35)
    fig.savefig("spectrum.png", dpi=200)
    print("Saved radial_densities.png, spectrum.png")
