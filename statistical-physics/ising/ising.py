"""2D Ising model: ferromagnetic-paramagnetic transition with the Metropolis algorithm.

Spins s = +-1 on an L x L square lattice with periodic boundaries, nearest-neighbour
coupling J = 1, no external field. A spin flip costs dE = 2 s (sum of 4 neighbours)
and is accepted with probability min(1, exp(-dE/T)). The lattice is updated as a
checkerboard: the spins of one colour do not touch each other, so half the lattice
flips at once with numpy. The magnetisation is checked against Onsager's exact
result M(T) = (1 - sinh(2/T)^-4)^(1/8) below Tc = 2 / ln(1 + sqrt 2) = 2.269.

    python3 ising.py           # check + figures
Units: J = k_B = 1.
"""
import sys
from pathlib import Path

import functools

import numpy as np

L = 48
TC = 2 / np.log(1 + np.sqrt(2))
rng = np.random.default_rng(1)


@functools.cache
def colours(n):
    """The two checkerboard sublattices of an n x n lattice."""
    i, j = np.indices((n, n))
    return [(i + j) % 2 == c for c in (0, 1)]


def neighbours(s):
    return np.roll(s, 1, 0) + np.roll(s, -1, 0) + np.roll(s, 1, 1) + np.roll(s, -1, 1)


def sweep(s, T):
    """One Monte Carlo sweep: every spin gets one Metropolis trial."""
    for c in colours(len(s)):
        dE = 2 * s * neighbours(s)
        flip = c & (rng.random(s.shape) < np.exp(-dE / T))
        s[flip] *= -1


def measure(T, n_therm=1000, n_meas=4000, s=None):
    """Mean |m|, energy, susceptibility and heat capacity per spin at temperature T."""
    s = np.ones((L, L), dtype=np.int8) if s is None else s
    for _ in range(n_therm):
        sweep(s, T)
    m, e = np.empty(n_meas), np.empty(n_meas)
    for k in range(n_meas):
        sweep(s, T)
        m[k] = abs(s.mean())
        e[k] = -(s * (np.roll(s, 1, 0) + np.roll(s, 1, 1))).mean()
    n = L * L
    return m.mean(), e.mean(), n * m.var() / T, n * e.var() / T**2, s


def onsager(T):
    T = np.asarray(T, float)
    return np.where(T < TC, np.abs(1 - np.sinh(2 / T) ** -4) ** 0.125, 0.0)


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    # Optional house style (figstyle.py at the repository root)
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import figstyle
    figstyle.apply()
    si_unit, si_qty = figstyle.unit, figstyle.qty

    # 1. Scan in temperature, each run starting from the previous configuration
    temps = np.unique(np.round(np.concatenate([np.arange(1.0, 2.0, 0.2), np.arange(2.0, 2.6, 0.05),
                                     np.arange(2.6, 3.61, 0.2)]), 3))
    res, s = [], None
    for T in temps:
        m, e, chi, c, s = measure(T, s=s)
        res.append((m, e, chi, c))
        print(f"T = {T:.2f}   |m| = {m:.3f}   Onsager {onsager(T):.3f}   chi = {chi:7.2f}   C = {c:.2f}")
    m, e, chi, c = np.array(res).T

    # Check: |m| follows Onsager away from Tc, and chi peaks near Tc (finite size shifts it up)
    far = temps < 2.1
    assert np.all(np.abs(m[far] - onsager(temps[far])) < 0.02), "|m| off Onsager"
    t_peak = temps[np.argmax(chi)]
    print(f"chi peaks at T = {t_peak:.2f} (Tc = {TC:.3f})")
    assert abs(t_peak - TC) < 0.15

    # 2. Order parameter, susceptibility and heat capacity
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.4))
    tt = np.linspace(1, 3.6, 500)
    axes[0].plot(tt, onsager(tt), color="0.6", label="Onsager, $L = \\infty$")
    for ax, y, name in zip(axes, [m, chi, c], ["$\\langle|m|\\rangle$", "$\\chi$", "$C$"]):
        ax.plot(temps, y, "o", color="C0", ms=4, label=f"Metropolis, $L = {figstyle.num(L)}$")
        ax.axvline(TC, color="C1", lw=0.8, ls="--")
        ax.set(xlabel=r"$T$ $\left[\frac{J}{k_B}\right]$", ylabel=name)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside upper center", ncol=2)
    fig.savefig("transition.png", dpi=200)

    # 3. Snapshots below, at and above Tc, from a random start
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.2))
    for ax, T in zip(axes, [1.8, TC, 3.0]):
        s = np.where(rng.random((128, 128)) < 0.5, 1, -1).astype(np.int8)    # larger lattice, domains show
        for _ in range(1500):
            sweep(s, T)
        ax.imshow(s, cmap=ListedColormap(["white", "C0"]), interpolation="nearest")
        ax.set(title="$T \\approx T_c$" if T == TC else f"$T = {figstyle.num(f'{T:.2f}')}$", xticks=[], yticks=[])
    fig.savefig("snapshots.png", dpi=200)
    print("Saved transition.png, snapshots.png")
