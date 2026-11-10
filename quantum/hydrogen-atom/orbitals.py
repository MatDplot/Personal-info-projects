"""Hydrogen orbitals in the xz plane, from the exact wave functions.

psi_nlm = R_nl(r) Y_lm(theta, phi) with the real spherical harmonics (p_x, p_z,
d_xz, d_z2, ...). The radial parts come from the closed form with Laguerre
polynomials, and are checked against the finite-difference solutions of
hydrogen.py: two independent routes to the same u(r) = r R(r).

    python3 orbitals.py        # check + figure
Atomic units, lengths in a0.
"""
import sys
from math import factorial
from pathlib import Path

import numpy as np
from scipy.special import genlaguerre

from hydrogen import levels, r as R_GRID, DR


def radial(n, l, r):
    """Normalised R_nl(r)."""
    rho = 2 * r / n
    norm = np.sqrt((2 / n) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return norm * np.exp(-rho / 2) * rho**l * genlaguerre(n - l - 1, 2 * l + 1)(rho)


# Real angular parts in Cartesian form (unit vector components x/r, y/r, z/r)
ANGULAR = {
    "s":     (0, lambda x, y, z: np.full_like(x, np.sqrt(1 / (4 * np.pi)))),
    "p_x":   (1, lambda x, y, z: np.sqrt(3 / (4 * np.pi)) * x),
    "p_z":   (1, lambda x, y, z: np.sqrt(3 / (4 * np.pi)) * z),
    "d_xz":  (2, lambda x, y, z: np.sqrt(15 / (4 * np.pi)) * x * z),
    "d_z^2": (2, lambda x, y, z: np.sqrt(5 / (16 * np.pi)) * (3 * z**2 - 1)),
}


def orbital(n, name, x, z):
    l, ang = ANGULAR[name]
    r = np.hypot(x, z) + 1e-12
    return radial(n, l, r) * ang(x / r, 0 * x, z / r)


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Optional house style (figstyle.py at the repository root)
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import figstyle
    figstyle.apply()
    DIVERGING = figstyle.cmap("CMAP_DIVERGING", "RdBu_r")

    # 1. Check: Laguerre closed form = finite-difference solution, for every n <= 3
    for l in range(3):
        _, u_fd = levels(l)
        for k in range(3 - l):
            n = l + 1 + k
            u_ex = R_GRID * radial(n, l, R_GRID)
            norm = (u_ex**2).sum() * DR
            gap = np.abs(np.abs(u_fd[:, k]) - np.abs(u_ex)).max() / np.abs(u_ex).max()
            print(f"n = {n}, l = {l}   int u^2 dr = {norm:.5f}   max |u_fd - u_exact| / max u = {gap:.1e}")
            assert abs(norm - 1) < 1e-4 and gap < 1e-3

    # 2. Figure: |psi|^2 in the xz plane, phase of psi as colour
    panels = [(1, "s"), (2, "s"), (2, "p_z"), (2, "p_x"), (3, "s"), (3, "p_z"), (3, "d_z^2"), (3, "d_xz")]
    fig, axes = plt.subplots(2, 4, figsize=(9, 4.8))
    for ax, (n, name) in zip(axes.flat, panels):
        size = 4 * n * n
        g = np.linspace(-size, size, 401)
        x, z = np.meshgrid(g, g)
        psi = orbital(n, name, x, z)
        # Signed density on a log-like scale, so the outer lobes and nodes show
        shown = np.sign(psi) * np.abs(psi) ** 0.5
        m = np.abs(shown).max()
        ax.imshow(shown, extent=(-size, size, -size, size), origin="lower", cmap=DIVERGING, vmin=-m, vmax=m)
        letter, _, sub = name.partition("_")
        ax.set_title(f"${n}\\mathrm{{{letter}}}" + (f"_{{{sub}}}$" if sub else "$"))
        ax.set(xticks=[], yticks=[])
    fig.suptitle(r"Hydrogen orbitals in the $xz$ plane: sign of $\psi$ in colour, $|\psi|^{1/2}$ as intensity")
    fig.savefig("orbitals.png", dpi=200)
    print("Saved orbitals.png")
