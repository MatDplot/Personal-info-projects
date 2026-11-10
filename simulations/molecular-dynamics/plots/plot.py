"""Figures from the output of ./md (run from src/, reads src/Resultats).

    python3 ../plots/plot.py      # from src/, after one run per integration method

Energy conservation for every method found in Resultats/Energies, the speed
distribution against Maxwell-Boltzmann, and the radial distribution function g(r).
Replaces the gnuplot scripts of the original project.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AU_VELOCITY = 2.18769126e6       # m/s
OUT = Path(__file__).parent
RES = Path("Resultats")

# Optional house style (figstyle.py at the repository root)
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import figstyle
figstyle.apply()
si_unit, si_qty = figstyle.unit, figstyle.qty

# 1. Total energy drift, one curve per method: t [ps], Ep, Ec, Etot [eV], (Etot - E0)/E0 [%]
fig, ax = plt.subplots(figsize=(6, 3.8))
for method, c, ls in [("Euler", "C1", "-"), ("RK4", "C4", "-"), ("Verlet", "C0", "-"), ("Leapfrog", "C2", ":")]:
    f = RES / "Energies" / f"energies_{method}.txt"
    if f.exists():
        t, drift = np.loadtxt(f, usecols=(0, 4), unpack=True)
        ax.plot(t[1:], np.abs(drift[1:]), color=c, ls=ls, label=method)
ax.set(yscale="log", xlabel="$t$ " + si_unit(r"\pico\second", "ps"),
       ylabel=r"$\frac{|E_\mathrm{tot}(t) - E_0|}{|E_0|}$ " + si_unit(r"\percent", r"\%"))
ax.legend(ncol=4, loc="upper center")
ax.set_ylim(top=1)
fig.savefig(OUT / "energy_drift.png", dpi=200)

# 2. Speed distribution at the end of the run against Maxwell-Boltzmann
v, h = np.loadtxt(RES / "Histogrammes/Vitesses/histo_modules_vitesses.txt", unpack=True)
_, h_th = np.loadtxt(RES / "Histogrammes/Vitesses/histo_modules_vitesses_th.txt", unpack=True)
v *= AU_VELOCITY
fig, ax = plt.subplots(figsize=(5.5, 3.8))
ax.bar(v, 100 * h, width=v[1] - v[0], color="C0", alpha=0.7, label="Simulation")
ax.plot(v, 100 * h_th, color="k", label="Maxwell-Boltzmann")
ax.set(xlabel="$|v|$ " + si_unit(r"\metre\per\second", r"m\,s^{-1}"), ylabel="Atoms per bin " + si_unit(r"\percent", r"\%"), xlim=(0, v[h > 0].max()))
ax.legend()
fig.savefig(OUT / "speeds.png", dpi=200)

# 3. Radial distribution function
r, g = np.loadtxt(RES / "Histogrammes/RDF/histo_RDF.txt", unpack=True)
fig, ax = plt.subplots(figsize=(5.5, 3.8))
ax.plot(r[1:], g[1:], color="C0")
ax.axhline(1, color="0.6", lw=0.8)
ax.set(xlabel="$r$ " + si_unit(r"\text{\AA}", r"\AA"), ylabel="$g(r)$", xlim=(0, r.max()))
fig.savefig(OUT / "rdf.png", dpi=200)
print("Saved energy_drift.png, speeds.png, rdf.png in", OUT)

# 4. Simulated annealing (NV mode), if it was run: heating then cooling, and g(r) of
# the final, cold structure from its positions (the NV mode writes no g(r) itself)
nv = RES / "NV"
if (nv / "Temperature_NV.txt").exists():
    t_nv, temp = np.loadtxt(nv / "Temperature_NV.txt", unpack=True)
    lines = (nv / "geom_NV_Verlet.xyz").read_text().split("\n")
    rows = [l.split() for l in lines[2:2 + int(lines[0])]]
    pos = np.array([[float(x) for x in row[1:]] for row in rows if row[0] == "Ar"])
    box = max(float(row[1]) for row in rows if row[0] == "Xe")       # box corners are written as Xe
    d = pos[:, None] - pos[None]
    d -= box * np.round(d / box)                                       # minimum image
    dist = np.sqrt((d**2).sum(-1))[np.triu_indices(len(pos), 1)]
    edges = np.arange(0.05, box / 2, 0.1)
    h_nv, _ = np.histogram(dist, edges)
    r_nv = (edges[1:] + edges[:-1]) / 2
    g_nv = h_nv / (len(pos) / 2 * len(pos) / box**3 * 4 * np.pi * r_nv**2 * 0.1)

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    axes[0].plot(t_nv, temp, color="C0", lw=0.8)
    axes[0].set(xlabel="$t$ " + si_unit(r"\pico\second", "ps"), ylabel="$T$ " + si_unit(r"\kelvin", "K"))
    axes[1].plot(r[1:], g[1:], color="0.6", label="Liquid at " + si_qty(90, r"\kelvin", "K"))
    axes[1].plot(r_nv, g_nv, color="C1", label=f"After annealing, {len(pos)} atoms at " + si_qty(10, r"\kelvin", "K"))
    axes[1].axhline(1, color="0.6", lw=0.5)
    axes[1].set(xlabel="$r$ " + si_unit(r"\text{\AA}", r"\AA"), ylabel="$g(r)$", xlim=(0, box / 2))
    axes[1].legend(fontsize=8, loc="upper right")
    fig.savefig(OUT / "annealing.png", dpi=200)
    print("Saved annealing.png")
