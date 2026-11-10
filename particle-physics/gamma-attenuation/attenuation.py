"""Attenuation of gamma rays in aluminium and copper (L3 nuclear physics lab).

A NaI scintillator counts the photopeak of a source behind absorbers of growing
thickness x. Beer-Lambert, I(x) = I0 exp(-mu x), makes ln(I/I0) a straight line
of slope -mu. Sources: 137Cs (662 keV) and 22Na (511 keV annihilation photons and
the 1275 keV line). The fitted mu are compared with NIST XCOM values.

    python3 attenuation.py     # fits + check + figure
Counts are photopeak integrals over the same run time; errors are Poisson.
"""
import sys
from pathlib import Path

import numpy as np

# NIST XCOM mass attenuation coefficients mu/rho [cm^2/g], coherent scattering included
# (physics.nist.gov/PhysRefData/XrayMassCoef), and densities [g/cm^3]
NIST_E = [500, 600, 800, 1250, 1500]       # keV
NIST = {"Al": ([8.445e-2, 7.802e-2, 6.841e-2, 5.496e-2, 5.006e-2], 2.699),
        "Cu": ([8.362e-2, 7.625e-2, 6.605e-2, 5.261e-2, 4.803e-2], 8.96)}

# Measured photopeak counts: (absorber, source, energy [keV], thickness [cm], counts,
# counts without absorber)
DATA = [
    ("Al", "$^{137}$Cs", 662, [2, 4, 6, 8, 10], [75731, 49951, 31757, 20289, 14398], 115628),
    ("Cu", "$^{137}$Cs", 662, [1, 2, 3, 4, 5], [59350, 31341, 15072, 8591, 5346], 115628),
    ("Cu", "$^{22}$Na", 511, [1, 2, 3, 4, 5], [83697, 40664, 19660, 8930, 3815], 173335),
    ("Cu", "$^{22}$Na", 1275, [1, 2, 3, 4, 5], [21115, 13743, 8736, 5683, 2989], 33624),
]


def mu_nist(material, e_kev):
    """Linear attenuation coefficient [cm^-1], log-log interpolation of the NIST table."""
    mu_rho, density = NIST[material]
    return density * np.exp(np.interp(np.log(e_kev), np.log(NIST_E), np.log(mu_rho)))


def fit_mu(x, n, n0):
    """Weighted straight-line fit of ln(n/n0) = -mu x + b. Returns mu, its error, b, y, sigma_y."""
    x, n = np.asarray(x, float), np.asarray(n, float)
    y = np.log(n / n0)
    sy = np.sqrt(1 / n + 1 / n0)
    (slope, b), cov = np.polyfit(x, y, 1, w=1 / sy, cov="unscaled")
    return -slope, np.sqrt(cov[0, 0]), b, y, sy


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Optional house style (figstyle.py at the repository root)
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import figstyle
    figstyle.apply()
    si_unit, si_qty = figstyle.unit, figstyle.qty

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for (mat, src, e, x, n, n0), c in zip(DATA, ["C2", "C0", "C1", "C4"]):
        label = f"{mat}, {src} {e} keV"
        mu, dmu, b, y, sy = fit_mu(x, n, n0)
        mu_ref = mu_nist(mat, e)
        print(f"{label.replace('$^{', '').replace('}$', ''):22s} mu = {mu:.3f} +- {dmu:.3f} cm^-1   NIST {mu_ref:.3f}   "
              f"({100 * (mu / mu_ref - 1):+.0f} %)")
        # Check: within 15 %, the spread expected from a lab setup (geometry, scattered photons)
        assert abs(mu / mu_ref - 1) < 0.15, f"mu off for {label}"
        ax.errorbar(x, y, sy, fmt="o", color=c, ms=4, capsize=2)
        xx = np.linspace(0, max(x), 2)
        ax.plot(xx, b - mu * xx, color=c, lw=1, label=f"{mat}, {src} " + si_qty(e, r"\kilo\electronvolt", "keV") + ": $\\mu =$ " + si_qty(f"{mu:.2f}", r"\per\centi\metre", "cm^{-1}"))
    ax.set(xlabel="Absorber thickness $x$ " + si_unit(r"\centi\metre", "cm"), ylabel=r"$\ln\frac{I}{I_0}$", xlim=(0, 10.5))
    ax.legend(fontsize=9)
    fig.savefig("attenuation.png", dpi=200)
    print("Saved attenuation.png")
