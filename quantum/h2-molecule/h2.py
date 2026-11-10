"""Two textbook pictures of a diatomic molecule.

1. H2+ ion in the LCAO approximation: the bonding and antibonding orbitals
   sigma_g,u = (phi_a +- phi_b) / sqrt(2 (1 +- S)) built from two 1s orbitals,
   their energy curves E(R), and the bond length and binding energy they predict.
2. Vibrations of the nuclei in a Morse potential V = De (1 - exp(-alpha x))^2,
   solved with finite differences and checked against the exact Morse levels
   E_v = hbar w (v + 1/2) - [hbar w (v + 1/2)]^2 / (4 De), next to the harmonic
   approximation E_v = hbar w (v + 1/2).

    python3 h2.py          # checks + figures
Atomic units: hbar = m_e = e = 4 pi eps0 = 1.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

HARTREE_EV = 27.211386


# --- 1. H2+ in the LCAO approximation -------------------------------------------
def overlap(R):
    return np.exp(-R) * (1 + R + R**2 / 3)


def lcao_energies(R):
    """E - E_1s for sigma_g and sigma_u, nuclear repulsion included (hartree)."""
    S = overlap(R)
    j = (1 - (1 + R) * np.exp(-2 * R)) / R     # electron on a, attracted by nucleus b
    k = (1 + R) * np.exp(-R)                   # exchange (resonance) integral
    return 1 / R - (j + k) / (1 + S), 1 / R - (j - k) / (1 - S)


def orbitals(R, x, z):
    """sigma_g and sigma_u on a grid, nuclei at z = -R/2 and +R/2."""
    phi_a = np.exp(-np.hypot(x, z + R / 2)) / np.sqrt(np.pi)
    phi_b = np.exp(-np.hypot(x, z - R / 2)) / np.sqrt(np.pi)
    S = overlap(R)
    return (phi_a + phi_b) / np.sqrt(2 * (1 + S)), (phi_a - phi_b) / np.sqrt(2 * (1 - S))


# --- 2. Morse vibrations ----------------------------------------------------------
DE, ALPHA, MU = 0.1842, 0.99852, 912.1         # course values, close to H2
OMEGA = ALPHA * np.sqrt(2 * DE / MU)
x = np.linspace(-1.5, 12, 4000)                # x = R - Re
DX = x[1] - x[0]
V_MORSE = DE * (1 - np.exp(-ALPHA * x)) ** 2


def vib_levels(V, n_states):
    diag = 1 / (MU * DX**2) + V
    off = np.full(len(x) - 1, -0.5 / (MU * DX**2))
    e, psi = eigh_tridiagonal(diag, off, select="i", select_range=(0, n_states - 1))
    return e, psi / np.sqrt(DX)


def morse_exact(v):
    return OMEGA * (v + 0.5) - (OMEGA * (v + 0.5)) ** 2 / (4 * DE)


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Optional house style (figstyle.py at the repository root)
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import figstyle
    figstyle.apply()
    si_unit, si_qty = figstyle.unit, figstyle.qty
    DIVERGING = figstyle.cmap("CMAP_DIVERGING", "RdBu_r")

    # 1a. Check: the orbitals are normalised (cylindrical integral 2 pi rho drho dz)
    rho, z = np.meshgrid(np.linspace(0, 15, 1501), np.linspace(-15, 15, 3001))
    g, u = orbitals(2.0, rho, z)
    d = (rho[0, 1] - rho[0, 0]) * (z[1, 0] - z[0, 0])
    for name, psi in [("sigma_g", g), ("sigma_u", u)]:
        norm = (2 * np.pi * rho * psi**2).sum() * d
        print(f"<{name}|{name}> = {norm:.4f}")
        assert abs(norm - 1) < 2e-3

    # 1b. Bond length and binding energy from the bonding curve
    R = np.linspace(0.6, 8, 2000)
    e_g, e_u = lcao_energies(R)
    i = np.argmin(e_g)
    print(f"LCAO: Re = {R[i]:.2f} a0, De = {-e_g[i] * HARTREE_EV:.2f} eV "
          "(exact H2+: 2.00 a0, 2.79 eV)")
    assert 2.4 < R[i] < 2.6 and 1.6 < -e_g[i] * HARTREE_EV < 1.9

    # 2. Check: finite-difference Morse levels against the exact formula
    n_bound = int(np.sqrt(2 * MU * DE) / ALPHA - 0.5) + 1
    e_m, psi_m = vib_levels(V_MORSE, n_bound)
    e_ex = morse_exact(np.arange(n_bound))
    print(f"Morse: {n_bound} bound levels, hbar w = {OMEGA * HARTREE_EV * 1000:.0f} meV, "
          f"max |E_fd - E_exact| = {np.abs(e_m - e_ex).max() * HARTREE_EV * 1000:.3f} meV")
    assert np.allclose(e_m, e_ex, atol=1e-5)

    # Figure 1: the two orbitals in a plane containing the nuclei
    xx, zz = np.meshgrid(np.linspace(-4, 4, 400), np.linspace(-4, 4, 400))
    g, u = orbitals(2.0, xx, zz)
    cmap = DIVERGING
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.4))
    for ax, psi, title in zip(axes, [g, u], [r"Bonding $\sigma_g$", r"Antibonding $\sigma_u$"]):
        m = np.abs(psi).max()
        ax.imshow(psi, extent=(-4, 4, -4, 4), cmap=cmap, vmin=-m, vmax=m, origin="lower")
        ax.contour(xx, zz, psi, levels=np.linspace(-m, m, 11), colors="k", linewidths=0.4)
        ax.plot([0, 0], [-1, 1], "o", color="k", ms=4)
        ax.set(title=title, xlabel="$x$ $[a_0]$", ylabel="$z$ $[a_0]$")
    fig.savefig("orbitals.png", dpi=200)

    # Figure 2: LCAO energy curves
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(R, e_g * HARTREE_EV, color="C0", label=r"$\sigma_g$ (bonding)")
    ax.plot(R, e_u * HARTREE_EV, color="C1", label=r"$\sigma_u$ (antibonding)")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.plot(R[i], e_g[i] * HARTREE_EV, "o", color="k")
    ax.annotate(f"$R_e = {figstyle.num(f'{R[i]:.2f}')}\\,a_0$, $D_e =$ " + si_qty(f"{-e_g[i] * HARTREE_EV:.2f}", r"\electronvolt", "eV"),
                (R[i], e_g[i] * HARTREE_EV), (3.3, -2.6))
    ax.set(xlim=(0.6, 8), ylim=(-3, 8), xlabel="$R$ $[a_0]$", ylabel=r"$E - E_{1s}$ " + si_unit(r"\electronvolt", "eV"))
    ax.legend()
    fig.savefig("lcao_energy.png", dpi=200)

    # Figure 3: Morse levels and wave functions, harmonic approximation for comparison
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ev = HARTREE_EV
    ax.plot(x, V_MORSE * ev, color="k", label="Morse")
    ax.plot(x, 0.5 * MU * OMEGA**2 * x**2 * ev, color="0.6", ls="--", label="Harmonic")
    for v in range(0, n_bound, 2):
        ax.hlines(e_m[v] * ev, x[0], x[-1], color="0.6", lw=0.4)
        ax.plot(x, (e_m[v] + 0.008 * psi_m[:, v]) * ev, color="C0", lw=1)
    ax.plot([], [], color="C0", lw=1, label=r"$\psi_v$, even $v$")
    ax.set(xlim=(-1, 6), ylim=(0, 1.15 * DE * ev), xlabel=r"$R - R_e$ $[a_0]$", ylabel="$E$ " + si_unit(r"\electronvolt", "eV"))
    fig.legend(*ax.get_legend_handles_labels(), loc="outside upper center", ncol=3)
    fig.savefig("vibrations.png", dpi=200)
    print("Saved orbitals.png, lcao_energy.png, vibrations.png")
