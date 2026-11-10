"""The simple pendulum beyond small angles.

theta'' = -w0^2 sin(theta) - gamma theta', integrated with scipy's adaptive RK45.
Without damping the period depends on the amplitude, T = 4 K(sin^2(theta0/2)) / w0
with K the complete elliptic integral of the first kind; small angles give back
T0 = 2 pi / w0 and Borda's correction T0 (1 + theta0^2 / 16).

    python3 pendulum.py        # check + figures
Units: w0 = 1, so T0 = 2 pi.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipk


def rhs(t, y, gamma=0.0):
    theta, omega = y
    return [omega, -np.sin(theta) - gamma * omega]


def period(theta0):
    """Measured period: time between two downward crossings of theta = 0."""
    up = lambda t, y: y[0]
    up.direction = -1
    sol = solve_ivp(rhs, (0, 200), [theta0, 0], events=up, rtol=1e-11, atol=1e-12)
    t = sol.t_events[0]
    return t[1] - t[0]


def period_exact(theta0):
    return 4 * ellipk(np.sin(theta0 / 2) ** 2)


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Optional house style (figstyle.py at the repository root)
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import figstyle
    figstyle.apply()
    si_unit, si_qty = figstyle.unit, figstyle.qty

    # 1. Check: measured period against the elliptic-integral formula
    amps = np.radians([5, 30, 60, 90, 120, 150, 170, 179])
    t_num = np.array([period(a) for a in amps])
    for a, t in zip(amps, t_num):
        print(f"theta0 = {np.degrees(a):5.0f} deg   T/T0 = {t / (2 * np.pi):.6f}   exact {period_exact(a) / (2 * np.pi):.6f}")
    assert np.allclose(t_num, period_exact(amps), rtol=1e-6), "period off the exact formula"

    # Figure 1: period against amplitude
    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    a = np.radians(np.linspace(0.5, 179.5, 400))
    ax.plot(np.degrees(a), period_exact(a) / (2 * np.pi), color="C0", label=r"Exact, $4K\left(\sin^2\frac{\theta_0}{2}\right)$")
    ax.plot(np.degrees(a), 1 + a**2 / 16, color="C1", ls="--", label=r"Borda, $1 + \frac{\theta_0^2}{16}$")
    ax.axhline(1, color="0.6", lw=0.8, label="Small angles")
    ax.plot(np.degrees(amps), t_num / (2 * np.pi), "o", color="k", label="Simulation")
    ax.set(xlabel="Amplitude $\\theta_0$ " + si_unit(r"\degree", "^\\circ"), ylabel=r"$\frac{T}{T_0}$", xlim=(0, 180), ylim=(0.9, 4.2))
    ax.legend(loc="upper left")
    fig.savefig("period.png", dpi=200)

    # Figure 2: phase portrait, with the separatrix and one damped trajectory
    fig, ax = plt.subplots(figsize=(6, 3.8))
    th, om = np.meshgrid(np.linspace(-2 * np.pi, 2 * np.pi, 400), np.linspace(-3, 3, 300))
    energy = 0.5 * om**2 - np.cos(th)
    ax.contour(th, om, energy, levels=np.linspace(-0.9, 2.5, 14), colors="0.6", linewidths=0.6, negative_linestyles="solid")
    ax.contour(th, om, energy, levels=[1], colors="C1", linewidths=1.2)
    damped = solve_ivp(rhs, (0, 40), [-2 * np.pi + 0.01, 2.5], args=(0.15,), max_step=0.02)
    ax.plot(damped.y[0], damped.y[1], color="C0", lw=1)
    ax.plot([], [], color="C1", label="Separatrix: rotation above, oscillation below")
    ax.plot([], [], color="C0", label=f"Damped, $\\gamma = {figstyle.num(0.15)}\\,\\omega_0$")
    ax.set(xlabel="$\\theta$ " + si_unit(r"\radian", "rad"), ylabel=r"$\frac{\dot\theta}{\omega_0}$", xlim=(-2 * np.pi, 2 * np.pi))
    ax.set_xticks(np.pi * np.arange(-2, 3), ["$-2\\pi$", "$-\\pi$", "0", "$\\pi$", "$2\\pi$"])
    fig.legend(*ax.get_legend_handles_labels(), loc="outside upper center", ncol=2)
    fig.savefig("phase_portrait.png", dpi=200)
    print("Saved period.png, phase_portrait.png")
