"""Magnetic field of a finite solenoid from the Biot-Savart law.

The solenoid is N circular loops of radius a stacked along z. The field of one
loop is written with complete elliptic integrals (exact, no discretisation of
the wire), and the solenoid's field is the sum over loops. On the axis it is
checked against the closed form for a finite solenoid,
    B_z(z) = (mu0 n I / 2) [ (z + L/2) / sqrt((z + L/2)^2 + a^2) - (z - L/2) / sqrt((z - L/2)^2 + a^2) ].

    python3 solenoid.py        # check + figure
Units: mu0 I = 1, lengths in units of the radius a.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import ellipe, ellipk

A, N_LOOPS, LENGTH = 1.0, 40, 6.0
Z_LOOPS = np.linspace(-LENGTH / 2, LENGTH / 2, N_LOOPS)


def loop_field(rho, z):
    """(B_rho, B_z) of one loop of radius A in the plane z = 0 (mu0 I = 1)."""
    rho = np.maximum(rho, 1e-9)
    q = (A + rho) ** 2 + z**2
    m = 4 * A * rho / q
    k, e = ellipk(m), ellipe(m)
    d = np.maximum((A - rho) ** 2 + z**2, 1e-6)          # on the wire itself the field diverges
    pre = 1 / (2 * np.pi * np.sqrt(q))
    b_z = pre * (k + (A**2 - rho**2 - z**2) / d * e)
    with np.errstate(invalid="ignore"):                   # 0 * inf on the axis, masked below
        b_rho = pre * z / rho * (-k + (A**2 + rho**2 + z**2) / d * e)
    return np.where(rho > 1e-6, b_rho, 0.0), b_z           # B_rho = 0 on the axis by symmetry


def solenoid_field(rho, z):
    b = [loop_field(rho, z - zl) for zl in Z_LOOPS]
    return sum(x[0] for x in b), sum(x[1] for x in b)


def axis_exact(z):
    """Uniform current sheet equivalent to the loops: each loop stands for one spacing."""
    span = LENGTH + LENGTH / (N_LOOPS - 1)
    n, half = N_LOOPS / span, span / 2
    return n / 2 * ((z + half) / np.hypot(z + half, A) - (z - half) / np.hypot(z - half, A))


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Optional house style (figstyle.py at the repository root)
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import figstyle
    figstyle.apply()
    si_unit, si_qty = figstyle.unit, figstyle.qty
    SEQUENTIAL = figstyle.cmap("CMAP_BLUE", "viridis")

    # 1. Check: one loop on its axis gives a^2 / (2 (a^2 + z^2)^(3/2))
    z = np.linspace(-3, 3, 61)
    assert np.allclose(loop_field(0, z)[1], A**2 / (2 * (A**2 + z**2) ** 1.5), rtol=1e-6)
    # 2. Check: the solenoid on its axis against the current-sheet formula, away from the ends
    z = np.linspace(-6, 6, 241)
    bz = solenoid_field(0, z)[1]
    gap = np.abs(bz / axis_exact(z) - 1)[np.abs(z) < LENGTH]
    print(f"Axis: B_z(0) = {bz[120]:.4f}, infinite solenoid n mu0 I = {N_LOOPS / LENGTH:.4f}; "
          f"max gap to the sheet formula = {gap.max():.1e}")
    assert gap.max() < 1e-2, "axis field off the closed form"

    # Figure: field lines and |B| in a plane containing the axis
    rr, zz = np.meshgrid(np.linspace(-3, 3, 241), np.linspace(-6, 6, 481))
    b_rho, b_z = solenoid_field(np.abs(rr), zz)
    b_x = np.sign(rr) * b_rho
    fig, axes = plt.subplots(2, 1, figsize=(6.5, 6.4), gridspec_kw={"height_ratios": [1.1, 1]})
    ax = axes[0]
    mag = np.hypot(b_x, b_z)
    im = ax.pcolormesh(zz, rr, np.log10(mag), cmap=SEQUENTIAL, shading="auto", rasterized=True)
    ax.streamplot(zz[:, 0], rr[0], b_z.T, b_x.T, color="w", density=1.2, linewidth=0.5, arrowsize=0.6)
    ax.plot(Z_LOOPS, np.full(N_LOOPS, A), ".", color="k", ms=2)
    ax.plot(Z_LOOPS, np.full(N_LOOPS, -A), ".", color="k", ms=2)
    ax.set(xlabel="$z$ $[a]$", ylabel="$x$ $[a]$", aspect="equal")
    fig.colorbar(im, ax=ax, label=r"$\log_{10} |B|$ $\left[\frac{\mu_0 I}{a}\right]$")
    ax = axes[1]
    ax.plot(z, bz, "o", color="k", ms=2, label=f"Biot–Savart, {figstyle.num(N_LOOPS)} loops")
    ax.plot(z, axis_exact(z), color="C0", label="Current sheet, closed form")
    ax.axhline(N_LOOPS / LENGTH, color="0.6", ls="--", lw=0.8, label="Infinite solenoid")
    ax.set(xlabel="$z$ $[a]$", ylabel=r"$B_z$ on the axis $\left[\frac{\mu_0 I}{a}\right]$")
    ax.legend(fontsize=8)
    fig.savefig("solenoid.png", dpi=200)
    print("Saved solenoid.png")
