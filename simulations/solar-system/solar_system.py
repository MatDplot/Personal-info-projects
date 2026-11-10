"""The Sun and the eight planets as a full N-body problem.

Every body attracts every other one (Newton). Positions and velocities are pushed
with velocity Verlet, which is symplectic: the energy oscillates but does not
drift, unlike the explicit Euler scheme shown for comparison. Each planet starts
at perihelion with the vis-viva speed of its real orbit (a, e), in one plane.
Checks: periods against Kepler's third law T^2 = a^3, and energy conservation.

    python3 solar_system.py    # checks + figures
Units: AU, year, solar mass, so G = 4 pi^2.
"""
import sys
from pathlib import Path

import numpy as np

G = 4 * np.pi**2
NAMES = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
MASS = np.array([1.660e-7, 2.448e-6, 3.003e-6, 3.227e-7, 9.546e-4, 2.858e-4, 4.366e-5, 5.150e-5])
A = np.array([0.387, 0.723, 1.000, 1.524, 5.203, 9.537, 19.19, 30.07])      # semi-major axis [AU]
ECC = np.array([0.206, 0.007, 0.017, 0.093, 0.048, 0.054, 0.047, 0.009])
PERI_ANGLE = np.radians([77, 131, 103, 336, 14, 93, 173, 48])               # longitude of perihelion


def initial_state():
    """Sun + planets at perihelion, total momentum zero.

    Each planet is placed relative to the centre of mass of the Sun and the planets
    inside its orbit (Jacobi coordinates). Starting it from the Sun instead would
    hand Neptune the Sun's wobble around Jupiter and bias its period by ~1 %.
    """
    pos, vel, m = [np.zeros(2)], [np.zeros(2)], [1.0]
    for a, e, w, mp in zip(A, ECC, PERI_ANGLE, MASS):
        mi = np.array(m)
        com_pos = (mi[:, None] * pos).sum(0) / mi.sum()
        com_vel = (mi[:, None] * vel).sum(0) / mi.sum()
        v_p = np.sqrt(G * (mi.sum() + mp) * (1 + e) / (a * (1 - e)))     # vis-viva at perihelion
        pos.append(com_pos + a * (1 - e) * np.array([np.cos(w), np.sin(w)]))
        vel.append(com_vel + v_p * np.array([-np.sin(w), np.cos(w)]))   # anticlockwise
        m.append(mp)
    pos, vel, m = np.array(pos), np.array(vel), np.array(m)
    vel -= (m[:, None] * vel).sum(0) / m.sum()
    return pos, vel, m


def acceleration(pos, m):
    d = pos[None, :, :] - pos[:, None, :]                                  # d[i, j] = r_j - r_i
    r3 = (d**2).sum(-1) ** 1.5
    np.fill_diagonal(r3, np.inf)
    return G * (m[None, :, None] * d / r3[:, :, None]).sum(1)


def energy(pos, vel, m):
    d = np.linalg.norm(pos[None] - pos[:, None], axis=-1)
    i, j = np.triu_indices(len(m), 1)
    return 0.5 * (m * (vel**2).sum(1)).sum() - G * (m[i] * m[j] / d[i, j]).sum()


def run(t_end, dt, method="verlet", every=10):
    pos, vel, m = initial_state()
    acc = acceleration(pos, m)
    traj, e = [], []
    for k in range(int(t_end / dt)):
        if k % every == 0:
            traj.append(pos.copy())
            e.append(energy(pos, vel, m))
        if method == "verlet":
            vel += 0.5 * dt * acc
            pos += dt * vel
            acc = acceleration(pos, m)
            vel += 0.5 * dt * acc
        else:                                                              # explicit Euler
            pos, vel = pos + dt * vel, vel + dt * acc
            acc = acceleration(pos, m)
    return np.array(traj), np.array(e), dt * every


def periods(traj, dt_out):
    """Time for each planet to sweep 2 pi around the Sun (first full turn, interpolated)."""
    rel = traj[:, 1:] - traj[:, :1]
    turn = np.unwrap(np.arctan2(rel[..., 1], rel[..., 0]), axis=0)
    turn -= turn[0]
    t = np.arange(len(traj)) * dt_out
    return np.array([np.interp(2 * np.pi, turn[:, p], t) for p in range(rel.shape[1])])


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Optional house style (figstyle.py at the repository root)
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import figstyle
    figstyle.apply()
    si_unit, si_qty = figstyle.unit, figstyle.qty

    # 1. A bit more than one Neptune year with Verlet, dt = 1/2 day
    traj, e, dt_out = run(170, 1 / 730)
    T = periods(traj, dt_out)
    for name, a, t in zip(NAMES, A, T):
        print(f"{name:8s} T = {t:8.3f} yr   Kepler a^1.5 = {a**1.5:8.3f} yr")
    drift = np.abs(e / e[0] - 1).max()
    print(f"Verlet: max |dE/E| over 170 yr = {drift:.1e}")
    # Planets feel each other, so T differs from a^1.5 by about the planet/Sun mass ratio
    assert np.allclose(T, A**1.5, rtol=5e-3), "periods off Kepler"
    assert drift < 1e-5, "energy not conserved"

    # 2. Euler against Verlet over 10 years, same step
    _, e_eu, dt_eu = run(10, 1 / 730, "euler", every=1)
    _, e_ve, _ = run(10, 1 / 730, "verlet", every=1)
    print(f"10 yr: Euler |dE/E| = {abs(e_eu[-1] / e_eu[0] - 1):.1e}, Verlet {abs(e_ve[-1] / e_ve[0] - 1):.1e}")

    # Figure 1: orbits, inner planets over their first 2 years only: with 176 steps per
    # orbit, Verlet makes Mercury's ellipse precess by ~1e-3 rad per turn (numerical, not physics)
    colours = ["0.6", "C4", "C0", "C1", "C5", "C3", "C2", "k"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
    for ax, planets, lim, n in zip(axes, [range(4), range(8)], [1.8, 42], [int(2 / dt_out), None]):
        ax.plot(0, 0, "o", color="#ffb800", ms=6)                             # Sun, ATLAS yellow
        for p in planets:
            ax.plot(traj[:n, p + 1, 0], traj[:n, p + 1, 1], color=colours[p], lw=0.8, label=NAMES[p], rasterized=True)
        ax.set(xlim=(-lim, lim), ylim=(-lim, lim), aspect="equal", xlabel="$x$ " + si_unit(r"\astronomicalunit", "AU"), ylabel="$y$ " + si_unit(r"\astronomicalunit", "AU"))
    fig.legend(*axes[1].get_legend_handles_labels(), loc="outside right upper")
    fig.savefig("orbits.png", dpi=200)

    # Figure 2: energy error, Euler against Verlet
    fig, ax = plt.subplots(figsize=(6, 3.6))
    t = np.arange(len(e_eu)) * dt_eu
    ax.semilogy(t[1:], np.abs(e_eu[1:] / e_eu[0] - 1), color="C1", label="Euler")
    ax.semilogy(t[1:], np.abs(e_ve[1:] / e_ve[0] - 1), color="C0", label="Velocity Verlet")
    ax.set(xlabel="$t$ " + si_unit(r"\year", "yr"), ylabel=r"$\frac{|E(t) - E_0|}{|E_0|}$",
           title=r"Same step, $\Delta t =$ " + si_qty(0.5, r"\day", "d"))
    ax.legend()
    fig.savefig("energy.png", dpi=200)
    print("Saved orbits.png, energy.png")
