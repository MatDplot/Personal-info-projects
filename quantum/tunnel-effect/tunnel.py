"""Quantum tunnelling of a Gaussian wave packet through a square barrier.

Solves the 1D time-dependent Schrodinger equation with the split-step Fourier
method (Strang splitting), measures the transmitted probability T and checks it
against the exact plane-wave result averaged over the packet's momenta.
Method as in J. VanderPlas, "Quantum Python" (2012); code written from scratch.

    python3 tunnel.py          # check + figures + animation
Units: hbar = 1. Lengths in units of L = hbar / sqrt(2 m V0).
"""
import sys
from pathlib import Path

import numpy as np

HBAR, M = 1.0, 1.9
V0 = 1.5
L = HBAR / np.sqrt(2 * M * V0)
A = 3 * L                     # barrier width
X0, D = -60.0, 10.0           # packet start and width
N, DX = 2**13, 0.05           # grid: x in [-205, 205], periodic
DT = 0.01
x = DX * (np.arange(N) - N // 2)
k = 2 * np.pi * np.fft.fftfreq(N, DX)   # FFT order, not sorted
V = np.where((x > 0) & (x < A), V0, 0.0)


def packet(k0):
    return (D * np.sqrt(np.pi)) ** -0.5 * np.exp(-0.5 * ((x - X0) / D) ** 2 + 1j * k0 * x)


def simulate(e_ratio, n_frames=0):
    """Run until the packet has left the barrier; return T, R and frames of |psi|^2."""
    k0 = np.sqrt(2 * M * e_ratio * V0) / HBAR
    t_end = 2 * abs(X0) / (HBAR * k0 / M)
    n_steps = int(t_end / DT)
    half = np.exp(-0.5j * V * DT / HBAR)
    kin = np.exp(-0.5j * HBAR * k**2 * DT / M)
    psi = packet(k0)
    every = max(1, n_steps // n_frames) if n_frames else 0
    frames = []
    for i in range(n_steps):
        if every and i % every == 0:
            frames.append((i * DT, np.abs(psi) ** 2))
        psi = half * np.fft.ifft(kin * np.fft.fft(half * psi))
    rho = np.abs(psi) ** 2
    return rho[x > A].sum() * DX, rho[x < 0].sum() * DX, frames


def t_plane(e):
    """Exact transmission of a plane wave of energy e (any e > 0)."""
    kappa = np.sqrt(2 * M * (V0 - e + 0j)) / HBAR
    return (1 / (1 + V0**2 * np.sinh(kappa * A) ** 2 / (4 * e * (V0 - e)))).real


def t_packet(e_ratio):
    """Plane-wave T averaged over the packet's momentum distribution."""
    k0 = np.sqrt(2 * M * e_ratio * V0) / HBAR
    kk = np.linspace(k0 - 6 / D, k0 + 6 / D, 4001)
    w = np.exp(-((kk - k0) * D) ** 2)
    e = (HBAR * kk) ** 2 / (2 * M)
    tk = np.where(kk > 0, t_plane(np.where(kk > 0, e, 1.0)), 0.0)   # k < 0 goes back
    return (tk * w).sum() / w.sum()


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    # Optional house style (figstyle.py at the repository root)
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import figstyle
    figstyle.apply()
    si_unit, si_qty = figstyle.unit, figstyle.qty

    # 1. Check against theory over a range of energies, below and above the barrier
    ratios = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5]
    t_num = []
    for r in ratios:
        t, refl, _ = simulate(r)
        assert abs(t + refl - 1) < 1e-4, "T + R != 1"   # slow tail still in the barrier
        t_num.append(t)
        print(f"E/V0 = {r:.1f}   T sim = {t:.4f}   T theory = {t_packet(r):.4f}")
    for r, t in zip(ratios, t_num):
        assert abs(t - t_packet(r)) < 0.05 * t_packet(r) + 2e-3, f"T off at E/V0 = {r}"

    # 2. T(E/V0): simulation against the exact curve
    fig, ax = plt.subplots(figsize=(6, 4))
    er = np.linspace(0.05, 1.6, 300)
    ax.plot(er, [t_packet(r) for r in er], color="C0", label="Exact, averaged over the packet")
    ax.plot(er, t_plane(er * V0), color="0.6", ls="--", label="Exact, plane wave")
    ax.plot(ratios, t_num, "o", color="k", label="Simulation")
    ax.axvline(1, color="0.6", lw=0.8)
    ax.set_xlabel(r"$\frac{E}{V_0}$")
    ax.set_ylabel("Transmission $T$")
    ax.legend()
    fig.savefig("transmission.png", dpi=200)
    fig.savefig("transmission.pdf")

    # 3. |psi|^2 before, during and after the barrier, at E = 0.2 V0
    t, refl, frames = simulate(0.2, n_frames=120)
    scale = 0.8 * V0 / frames[0][1].max()
    fig, ax = plt.subplots(figsize=(7, 3.2))
    ax.fill_between(x, V, color="0.6", alpha=0.6, lw=0, label="Barrier $V(x)$")
    for (tt, rho), c in zip([frames[0], frames[len(frames) // 2], frames[-1]], ["0.6", "C1", "C0"]):
        ax.plot(x, scale * rho, color=c, label=f"$t = {figstyle.num(f'{tt:.0f}')}$")
    right = x > A
    ax.plot(x[right], 20 * scale * frames[-1][1][right], color="C0", ls="--", lw=1,
            label=f"Transmitted part $\\times {figstyle.num(20)}$")
    ax.set_xlim(-130, 100)
    ax.set_xlabel("$x$")
    ax.set_ylabel(r"$|\psi(x,t)|^2$ $[\mathrm{arb.\ units}]$")
    ax.set_yticks([])
    ax.set_title(f"$E = 0.2\\,V_0$: $T =$ " + si_qty(f"{100 * t:.2f}", r"\percent", r"\%") + ", $R =$ " + si_qty(f"{100 * refl:.2f}", r"\percent", r"\%"), fontsize=10)
    ax.legend(fontsize=9)
    fig.savefig("wave_packet.png", dpi=200)
    fig.savefig("wave_packet.pdf")

    # 4. Animation, transmitted part magnified so it shows. GIF frames need a
    # raster canvas, so this figure is drawn by Agg (mathtext) even under LaTeX.
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    fig, ax = plt.subplots(figsize=(7, 3.2))
    FigureCanvasAgg(fig)
    ax.fill_between(x, V, color="0.6", alpha=0.6, lw=0)
    line, = ax.plot([], [], color="C0")
    zoom, = ax.plot([], [], color="C1", lw=1, label=r"Transmitted part $\times 20$")
    label = ax.text(0.02, 0.9, "", transform=ax.transAxes)
    ax.set(xlim=(-130, 100), ylim=(0, 1.1 * V0), yticks=[], xlabel="$x$")
    ax.legend(loc="upper right", fontsize=9)

    def draw(i):
        tt, rho = frames[i]
        line.set_data(x, scale * rho)
        right = x > A
        zoom.set_data(x[right], 20 * scale * rho[right])
        label.set_text(f"t = {tt:.0f}")
        return line, zoom, label

    FuncAnimation(fig, draw, frames=len(frames)).save("tunnel.gif", writer=PillowWriter(fps=20), dpi=80)
    print("Saved transmission.{png,pdf}, wave_packet.{png,pdf}, tunnel.gif")
