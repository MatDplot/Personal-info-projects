"""Three-flavour neutrino oscillations in vacuum, and what a reactor experiment sees.

The flavour states are mixed with the mass states by the PMNS matrix U (three
angles and one CP phase). After a distance L, a neutrino of energy E born as
flavour alpha is detected as beta with probability
    P = | sum_i U_beta,i U*_alpha,i exp(-i m_i^2 L / 2E) |^2.
Checks: unitarity (each row of P sums to 1), the closed-form electron
antineutrino survival used by reactor experiments, and the two-flavour limit.
Then the spectrum a JUNO-like detector would record at 52.5 km, for the two
possible orderings of the neutrino masses.

    python3 oscillations.py    # checks + figures
Oscillation parameters: PDG 2022 (normal ordering unless said otherwise).
"""
import sys
from pathlib import Path

import numpy as np

K = 1.266932   # dm^2 L / 4E [rad] = K * dm^2 [eV^2] * L [km] / E [GeV]

PDG = dict(s12=0.307, s13=0.0220, s23=0.546, delta=1.36 * np.pi, dm21=7.53e-5, dm32=2.453e-3)
PDG_IO = dict(PDG, s23=0.539, dm32=-2.536e-3)


def pmns(s12, s13, s23, delta, **_):
    """PMNS matrix from sin^2 of the three angles and the CP phase."""
    c12, c13, c23 = (np.sqrt(1 - s) for s in (s12, s13, s23))
    s12, s13, s23 = (np.sqrt(s) for s in (s12, s13, s23))
    e = np.exp(-1j * delta)
    r23 = np.array([[1, 0, 0], [0, c23, s23], [0, -s23, c23]])
    r13 = np.array([[c13, 0, s13 * e], [0, 1, 0], [-s13 / e, 0, c13]])
    r12 = np.array([[c12, s12, 0], [-s12, c12, 0], [0, 0, 1]])
    return r23 @ r13 @ r12


def probabilities(L_km, E_GeV, p, anti=False):
    """P[alpha, beta, ...] for flavours e, mu, tau, broadcast over L / E."""
    u = pmns(**p).conj() if anti else pmns(**p)
    m2 = np.array([0, p["dm21"], p["dm21"] + p["dm32"]])
    phase = np.exp(-1j * 2 * K * np.multiply.outer(m2, np.asarray(L_km) / np.asarray(E_GeV)))
    amp = np.einsum("bi,ai,i...->ab...", u, u.conj(), phase)
    return np.abs(amp) ** 2


def p_ee_closed_form(L_km, E_GeV, p):
    """Reactor electron antineutrino survival, the formula used by Daya Bay and JUNO."""
    x = K * np.asarray(L_km) / np.asarray(E_GeV)
    d21, d31, d32 = x * p["dm21"], x * (p["dm21"] + p["dm32"]), x * p["dm32"]
    s12, s13 = p["s12"], p["s13"]
    sin2_2t12, sin2_2t13 = 4 * s12 * (1 - s12), 4 * s13 * (1 - s13)
    return (1 - (1 - s13) ** 2 * sin2_2t12 * np.sin(d21) ** 2
            - sin2_2t13 * ((1 - s12) * np.sin(d31) ** 2 + s12 * np.sin(d32) ** 2))


def reactor_spectrum(E_MeV):
    """Emitted antineutrinos x inverse-beta-decay cross section, arbitrary units.

    Flux: exponential-polynomial fits of Vogel & Engel (1989) for the four
    fissile isotopes, with typical fission fractions. Cross section: zeroth
    order in 1/M, sigma ~ E_e p_e with E_e = E_nu - 1.293 MeV.
    """
    fits = {"U235": (0.58, 0.870, -0.160, -0.0910), "Pu239": (0.29, 0.896, -0.239, -0.0981),
            "U238": (0.07, 0.976, -0.162, -0.0790), "Pu241": (0.06, 0.793, -0.080, -0.1085)}
    flux = sum(f * np.exp(a0 + a1 * E_MeV + a2 * E_MeV**2) for f, a0, a1, a2 in fits.values())
    e_e = np.clip(E_MeV - 1.293, 0.511, None)
    return flux * e_e * np.sqrt(e_e**2 - 0.511**2)


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Optional house style (figstyle.py at the repository root)
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import figstyle
    figstyle.apply()
    si_unit, si_qty = figstyle.unit, figstyle.qty

    # 1. Checks
    L = np.geomspace(0.1, 1e4, 2000)
    for p in (PDG, PDG_IO):
        for anti in (False, True):
            P = probabilities(L, 0.004, p, anti)
            assert np.allclose(P.sum(axis=1), 1, atol=1e-12), "rows of P do not sum to 1"
        assert np.allclose(probabilities(L, 0.004, p, anti=True)[0, 0], p_ee_closed_form(L, 0.004, p), atol=1e-12)
    two = dict(PDG, s13=0.0)
    x = K * PDG["dm21"] * L / 0.004
    assert np.allclose(probabilities(L, 0.004, two)[0, 0], 1 - 4 * two["s12"] * (1 - two["s12"]) * np.sin(x) ** 2,
                       atol=1e-12), "two-flavour limit"
    print("Unitarity, closed-form P_ee and two-flavour limit: OK")

    # Figure 1: reactor antineutrino survival against distance at 4 MeV
    L = np.geomspace(0.1, 300, 3000)
    fig, ax = plt.subplots(figsize=(6, 3.6))
    ax.plot(L, probabilities(L, 0.004, PDG, anti=True)[0, 0], color="C0", label="Three flavours")
    ax.plot(L, probabilities(L, 0.004, two, anti=True)[0, 0], color="C1", ls="--", label=f"Two flavours ($\\theta_{{13}} = {figstyle.num(0)}$)")
    for name, d in [("Daya Bay", 1.6), ("JUNO", 52.5), ("KamLAND", 180)]:
        ax.axvline(d, color="0.6", lw=0.8)
        ax.text(d * 1.08, 0.04, figstyle.halo(name), rotation=90, fontsize=8, color="0.4")
    ax.set(xscale="log", xlabel="Distance $L$ " + si_unit(r"\kilo\metre", "km"), ylabel="$P(\\bar\\nu_e \\to \\bar\\nu_e)$ at " + si_qty(4, r"\mega\electronvolt", "MeV"), ylim=(0, 1.05), xlim=(0.1, 300))
    ax.legend(loc=(0.02, 0.4))
    fig.savefig("survival.png", dpi=200)

    # Figure 2: spectrum at 52.5 km, normal against inverted ordering
    E = np.linspace(1.8, 9, 3000)
    s = reactor_spectrum(E)
    fig, ax = plt.subplots(figsize=(6, 3.8))
    ax.plot(E, s / s.max(), color="0.6", label="No oscillation")
    for p, c, name in [(PDG, "C0", "Normal ordering"), (PDG_IO, "C1", "Inverted ordering")]:
        ax.plot(E, s * probabilities(52.5, E / 1000, p, anti=True)[0, 0] / s.max(), color=c, lw=1, label=name)
    ax.set(xlabel="$E_\\nu$ " + si_unit(r"\mega\electronvolt", "MeV"), ylabel=r"Events $[\mathrm{arb.\ units}]$", xlim=(1.8, 9), ylim=(0, 1.1))
    ax.legend()
    fig.savefig("reactor_spectrum.png", dpi=200)

    # Figure 3: accelerator appearance, nu_mu -> nu_tau over CERN - Gran Sasso (730 km)
    E = np.geomspace(0.5, 50, 1000)
    P = probabilities(730, E, PDG)
    fig, ax = plt.subplots(figsize=(6, 3.6))
    for (a, b), c, name in [((1, 1), "C0", "\\nu_\\mu \\to \\nu_\\mu"), ((1, 2), "C1", "\\nu_\\mu \\to \\nu_\\tau"),
                            ((1, 0), "C2", "\\nu_\\mu \\to \\nu_e")]:
        ax.plot(E, P[a, b], color=c, label=f"${name}$")
    ax.set(xscale="log", xlabel="$E_\\nu$ " + si_unit(r"\giga\electronvolt", "GeV"), ylabel="Probability at $L =$ " + si_qty(730, r"\kilo\metre", "km"), ylim=(0, 1.05))
    ax.legend()
    fig.savefig("accelerator.png", dpi=200)
    print("Saved survival.png, reactor_spectrum.png, accelerator.png")
