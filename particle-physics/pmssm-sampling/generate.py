"""Draw random pMSSM points and let SOFTSUSY decide which ones give a valid spectrum.

Each point is the 19 parameters of the phenomenological MSSM (gaugino masses,
higgsino mass mu, tan beta, trilinear couplings, sfermion masses) drawn
uniformly in the ranges below. SOFTSUSY runs them from the input scale with
the renormalisation group equations; a point fails when it has no electroweak
symmetry breaking, gives a tachyonic sfermion, a charged LSP, etc.
The result goes to points.csv, which pmssm.py reads.

    python3 generate.py 4000 /path/to/softsusy-4.1.9/.libs/softpoint.x
(the compiled binary, with its libsoft*.so next to it, not libtool's wrapper script)
SOFTSUSY: B. C. Allanach, Comput. Phys. Commun. 143 (2002) 305, softsusy.hepforge.org
"""
import csv
import os
import subprocess
import sys

import numpy as np

# name, SLHA EXTPAR code(s) (or "tanb"), range [GeV, except tan beta]
PARAMS = [
    ("M1", [1], (50, 5000)), ("M2", [2], (50, 5000)), ("M3", [3], (50, 5000)),
    ("mA", [26], (50, 5000)), ("tanb", "tanb", (2, 60)), ("mu", [23], (-1e4, 1e4)),
    ("At", [11], (-1e4, 1e4)), ("Ab", [12], (-1e4, 1e4)), ("Atau", [13], (-1e4, 1e4)),
    ("mq1L", [41, 42], (50, 5000)), ("mq3L", [43], (50, 5000)),
    ("muR", [44, 45], (50, 5000)), ("mtR", [46], (50, 5000)),
    ("mdR", [47, 48], (50, 5000)), ("mbR", [49], (50, 5000)),
    ("ml1L", [31, 32], (50, 5000)), ("mtauL", [33], (50, 5000)),
    ("meR", [34, 35], (50, 5000)), ("mtauR", [36], (50, 5000)),
]
# PDG codes of the masses kept from the output spectrum
MASSES = {"h0": 25, "gluino": 1000021, "chi10": 1000022, "chi1pm": 1000024, "stop1": 1000006}

HEADER = """BLOCK MODSEL
 1 0                 # general MSSM
BLOCK SMINPUTS
 1 1.279340000e+02   # alpha_em^-1(MZ) MSbar
 2 1.663700000e-05   # G_Fermi
 3 0.1179            # alpha_s(MZ) MSbar
 4 91.1876           # MZ pole
 5 4.18              # mb(mb) MSbar
 6 172.9             # mtop pole
 7 1.77686           # mtau pole
"""


def slha_input(values):
    v = dict(zip([p[0] for p in PARAMS], values))
    lines = [HEADER, "BLOCK MINPAR\n", f" 3 {v['tanb']}   # tan beta\n",
             "BLOCK EXTPAR\n", " 0 -1.0   # input scale = SUSY scale\n"]
    for name, codes, _ in PARAMS:
        if codes != "tanb":
            lines += [f" {c} {v[name]}   # {name}\n" for c in codes]
    return "".join(lines)


def run_point(softpoint, values):
    """Valid flag and selected pole masses (GeV) for one point."""
    env = dict(os.environ, LD_LIBRARY_PATH=os.path.dirname(os.path.abspath(softpoint)))
    out = subprocess.run([softpoint, "leshouches"], input=slha_input(values), env=env,
                         capture_output=True, text=True, timeout=60).stdout
    if "Block MASS" not in out or "problem" in out.lower():
        return False, {k: "" for k in MASSES}
    masses, in_mass = {}, False
    for line in out.splitlines():
        if line.lower().startswith("block"):
            in_mass = line.split()[1].upper() == "MASS"
            continue
        if in_mass and line.split() and line.split()[0].lstrip("-").isdigit():
            code, value = int(line.split()[0]), float(line.split()[1])
            for k, pdg in MASSES.items():
                if code == pdg:
                    masses[k] = abs(value)
    return True, masses


if __name__ == "__main__":
    n, softpoint = int(sys.argv[1]), sys.argv[2]
    rng = np.random.default_rng(2022)
    lo = np.array([p[2][0] for p in PARAMS])
    hi = np.array([p[2][1] for p in PARAMS])
    with open("points.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([p[0] for p in PARAMS] + ["valid"] + list(MASSES))
        n_valid = 0
        for i in range(n):
            values = lo + (hi - lo) * rng.random(len(PARAMS))
            valid, masses = run_point(softpoint, values)
            n_valid += valid
            w.writerow([f"{x:.6g}" for x in values] + [int(valid)] + [masses.get(k, "") for k in MASSES])
            if (i + 1) % 200 == 0:
                print(f"{i + 1} / {n} points, {n_valid} valid")
    print(f"Saved points.csv: {n_valid} valid out of {n}")
