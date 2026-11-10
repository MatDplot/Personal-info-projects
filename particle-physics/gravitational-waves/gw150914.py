"""Finding GW150914 in the open LIGO data with a matched filter.

32 s of strain from LIGO Hanford (H1) and Livingston (L1) around the first
detected binary black-hole merger, from the Gravitational Wave Open Science
Center. The noise PSD is estimated from the data (Welch), and the data are
matched-filtered against a bank of inspiral templates: the leading-order
(Newtonian) stationary-phase chirp, which depends only on the chirp mass
M_c = (m1 m2)^(3/5) / (m1 + m2)^(1/5). The chirp is continued up to three times
the innermost-stable-orbit frequency of an equal-mass binary, a crude stand-in
for the merger that a Newtonian template does not describe.

    python3 gw150914.py        # downloads 2 x 1.3 MB once, then check + figures
Data: GWOSC, https://gwosc.org (CC BY 4.0). Abbott et al., PRL 116 (2016) 061102.
"""
import sys
import urllib.request
from pathlib import Path

import numpy as np
from scipy.signal import spectrogram, welch
from scipy.signal.windows import tukey

G, C, M_SUN = 6.674e-11, 2.998e8, 1.989e30
FS, T0 = 4096, 1126259447                  # sampling rate [Hz], GPS start of the files
GPS_EVENT = 1126259462.4                   # GWTC-1
URL = "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/{d}-{d}1_GWOSC_4KHZ_R1-1126259447-32.txt.gz"
F_LOW = 20.0
F_CUT = 3.0                                # cut-off in units of f_ISCO
DATA = Path(__file__).parent / "data"


def load(det):
    """Strain of detector 'H' or 'L', downloaded from GWOSC on first use."""
    f = DATA / f"{det}1.txt.gz"
    if not f.exists():
        DATA.mkdir(exist_ok=True)
        urllib.request.urlretrieve(URL.format(d=det), f)
    return np.loadtxt(f)


def psd(strain):
    """One-sided noise PSD on the FFT grid of the full series (Welch, 4 s segments)."""
    f, p = welch(strain, fs=FS, nperseg=4 * FS)
    return np.interp(np.fft.rfftfreq(len(strain), 1 / FS), f, p)


def template(freqs, chirp_mass):
    """Newtonian stationary-phase chirp h(f) ~ f^(-7/6) exp(i Psi), merger at t = 0."""
    mc = chirp_mass * M_SUN * G / C**3                                   # in seconds
    f_isco = 1 / (6**1.5 * np.pi * 2 ** (6 / 5) * mc)                    # equal masses
    band = (freqs >= F_LOW) & (freqs <= F_CUT * f_isco)
    f = np.where(band, freqs, 1.0)
    h = f ** (-7 / 6) * np.exp(1j * (-np.pi / 4 + 3 / 128 * (np.pi * mc * f) ** (-5 / 3)))
    return np.where(band, h, 0)


def matched_filter(data_f, h, s, n):
    """SNR time series |<d, h(t)>| / sigma for every merger time t (both phases)."""
    df = FS / n
    sigma = np.sqrt(4 * df * (np.abs(h) ** 2 / s).sum())
    z = np.zeros(n, complex)
    z[: len(h)] = data_f * h.conj() / s
    return 4 * df * np.abs(np.fft.ifft(z) * n) / sigma


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

    det = {}
    for d in "HL":
        x = load(d)
        n = len(x)
        s = psd(x)
        det[d] = dict(x=x, s=s, xf=np.fft.rfft(x * tukey(n, 0.1)) / FS)
    n = len(det["H"]["x"])
    freqs = np.fft.rfftfreq(n, 1 / FS)
    t = np.arange(n) / FS
    search = (t > 4) & (t < 28)                                          # away from the edges

    # 1. Template bank in chirp mass, network SNR = quadrature sum of the two peaks
    bank = np.linspace(15, 45, 121)
    peaks = {d: [] for d in "HL"}
    for mc in bank:
        h = template(freqs, mc)
        for d in "HL":
            peaks[d].append(matched_filter(det[d]["xf"], h, det[d]["s"], n)[search].max())
    network = np.hypot(peaks["H"], peaks["L"])
    best = bank[np.argmax(network)]
    h = template(freqs, best)
    snr = {d: matched_filter(det[d]["xf"], h, det[d]["s"], n) for d in "HL"}
    t_peak = {d: t[search][np.argmax(snr[d][search])] for d in "HL"}
    for d in "HL":
        print(f"{d}1: SNR = {snr[d][search].max():.1f} at GPS {T0 + t_peak[d]:.3f}")
    print(f"Best chirp mass (detector frame) = {best:.1f} M_sun, network SNR = {network.max():.1f}; "
          f"L1 - H1 = {1000 * (t_peak['L'] - t_peak['H']):+.1f} ms")

    # Check: both detectors find the event at the catalogue time, within the light
    # travel time between the sites (10 ms), and the chirp mass near GWTC-1's
    # 28.6 M_sun x (1 + z = 1.09) = 31 M_sun in the detector frame
    assert all(abs(T0 + t_peak[d] - GPS_EVENT) < 0.05 for d in "HL"), "event not at the catalogue time"
    assert abs(t_peak["L"] - t_peak["H"]) < 0.011, "time delay larger than the light travel time"
    assert 25 < best < 37 and network.max() > 8, "chirp mass or SNR off"

    # Figure 1: SNR time series around the event, and SNR across the bank
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), gridspec_kw={"width_ratios": [1.4, 1]})
    near = np.abs(t - (GPS_EVENT - T0)) < 1
    for d, c in zip("HL", ["C1", "C0"]):
        axes[0].plot(t[near] - (GPS_EVENT - T0), snr[d][near], color=c, label=f"{d}1")
        axes[1].plot(bank, peaks[d], color=c, label=f"{d}1")
    axes[1].plot(bank, network, color="k", label="Network")
    axes[1].axvline(31.2, color="0.6", ls="--", lw=0.8)
    axes[0].set(xlabel="$t - t_{\\mathrm{GWTC}}$ " + si_unit(r"\second", "s"), ylabel="Matched-filter SNR")
    axes[1].set(xlabel="Template chirp mass $[M_\\odot]$", ylabel="Peak SNR")
    fig.legend(*axes[1].get_legend_handles_labels(), loc="outside upper center", ncol=3)
    fig.savefig("snr.png", dpi=200)

    # Figure 2: whitened, band-passed strain and time-frequency map in H1
    def whiten(d):
        w = np.where((freqs > 30) & (freqs < 350), det[d]["xf"] / np.sqrt(det[d]["s"]), 0)
        return np.fft.irfft(w, n)

    wh, wl = whiten("H"), whiten("L")
    shift = int(round((t_peak["L"] - t_peak["H"]) * FS))
    fig, axes = plt.subplots(2, 1, figsize=(7, 5.2), sharex=True)
    tt = t - (GPS_EVENT - T0)
    zoom = (tt > -0.35) & (tt < 0.08)
    axes[0].plot(tt[zoom], wh[zoom] / wh[zoom].std(), color="C1", lw=0.8, label="H1")
    axes[0].plot(tt[zoom], -np.roll(wl, -shift)[zoom] / wl[zoom].std(), color="C0", lw=0.8,
                 label="L1, shifted and sign-flipped")
    axes[0].set(ylabel="Whitened strain $[\\sigma]$")
    axes[0].legend(fontsize=8, loc="upper left")
    f_s, t_s, sxx = spectrogram(wh, fs=FS, nperseg=FS // 16, noverlap=FS // 16 - FS // 256)
    t_s -= GPS_EVENT - T0
    keep = (f_s > 20) & (f_s < 400)
    axes[1].pcolormesh(t_s, f_s[keep], sxx[keep], cmap=SEQUENTIAL, shading="gouraud", rasterized=True)
    # Frequency track of the best template, f = (5 / tau)^(3/8) / (8 pi Mc^(5/8)) with tau the time
    # to merger at the H1 peak, up to the template's cut-off
    mc = best * M_SUN * G / C**3
    tau = np.geomspace(0.35, 1e-4, 400)
    f_track = (5 / tau) ** 0.375 / (8 * np.pi * mc**0.625)
    f_cut = F_CUT / (6**1.5 * np.pi * 2 ** (6 / 5) * mc)
    on = f_track < f_cut
    axes[1].plot(t_peak["H"] - (GPS_EVENT - T0) - tau[on], f_track[on], color="w", ls="--", lw=1,
                 label=r"Best template, $\mathcal{M}_c = " + figstyle.num(f"{best:.0f}") + r"\,M_\odot$")
    axes[1].legend(fontsize=8, loc="upper left", labelcolor="w")
    axes[1].set(xlabel="$t - t_{\\mathrm{GWTC}}$ " + si_unit(r"\second", "s"), ylabel="Frequency " + si_unit(r"\hertz", "Hz"), xlim=(-0.35, 0.08), yscale="log",
                yticks=[32, 64, 128, 256], yticklabels=["32", "64", "128", "256"], ylim=(32, 380))
    fig.savefig("strain.png", dpi=200)
    print("Saved snr.png, strain.png")
