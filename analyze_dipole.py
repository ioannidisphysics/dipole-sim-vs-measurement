"""
Analysis for the lambda/2 dipole project.

Reads, for each setup A-D:
  data/<S>_hfss.csv   HFSS export   (header carries the unit: "Freq [MHz]" or "[GHz]")
  data/<S>_nec.csv    NEC2 export   (written by sim_dipole_nec.py)
  data/<S>_meas_[1-3].s1p           repeat measurements, Touchstone
  data/<S>_meas_<tag>.s1p           cable experiments: hand, ferrite, parallel

Writes plots/<S>_s11.png per setup and prints a markdown results table.
Runs fine before any measurement exists: the measurement columns stay empty.

    pip install numpy matplotlib scikit-rf
    python analyze_dipole.py
"""
import glob
import re

import matplotlib.pyplot as plt
import numpy as np

try:
    import skrf as rf
except ImportError:                                   # measurements not needed to plot the sims
    rf = None

UNITS = {"Hz": 1, "kHz": 1e3, "MHz": 1e6, "GHz": 1e9}
SETUPS = ["A", "B", "C", "D"]
VSWR2_DB = -9.54                                      # |S11| for VSWR = 2


def load_touchstone(path):
    n = rf.Network(path)
    return n.f, 20 * np.log10(np.abs(n.s[:, 0, 0]))


def load_csv(path):
    """HFSS/NEC csv. The frequency unit is read from the header, not assumed."""
    header = open(path).readline()
    scale = UNITS[re.search(r"\[(\w+)\]", header).group(1)]
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    return d[:, 0] * scale, d[:, 1]


def f_min(f, s_db):
    """Frequency of minimum |S11|, refined by parabolic interpolation."""
    i = int(np.argmin(s_db))
    if 0 < i < len(f) - 1:
        y0, y1, y2 = s_db[i - 1:i + 2]
        denom = y0 - 2 * y1 + y2
        if denom != 0:
            return f[i] + 0.5 * (y0 - y2) / denom * (f[i + 1] - f[i]), s_db[i]
    return f[i], s_db[i]


def bandwidth(f, s_db, thr=VSWR2_DB):
    """Contiguous band around the minimum where |S11| < thr. None if never matched."""
    below = np.where(s_db < thr)[0]
    if below.size == 0:
        return None
    lo_i, hi_i = below[0], below[-1]
    lo = np.interp(thr, s_db[:lo_i + 1][::-1], f[:lo_i + 1][::-1])
    hi = np.interp(thr, s_db[hi_i:], f[hi_i:])
    return lo, hi


def r_from_depth(s_db_min):
    """Input resistance implied by the depth, assuming X = 0 at the minimum."""
    g = 10 ** (s_db_min / 20)
    return 50 * (1 + g) / (1 - g)


rows, notes = [], []
for s in SETUPS:
    hfss = glob.glob(f"data/{s}_hfss.csv")
    nec = glob.glob(f"data/{s}_nec.csv")
    if not hfss:
        continue
    fh, sh = load_csv(hfss[0])
    f_sim, _ = f_min(fh, sh)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(fh / 1e6, sh, "C0", lw=2, label="HFSS (FEM)")
    f_nec = None
    if nec:
        fn, sn = load_csv(nec[0])
        f_nec, _ = f_min(fn, sn)
        ax.plot(fn / 1e6, sn, "C3--", lw=1.6, label="NEC2 (MoM)")

    # repeat measurements
    meas = sorted(glob.glob(f"data/{s}_meas_[0-9].s1p"))
    f_meas = None
    if meas and rf is not None:
        fres, depth = [], []
        for k, p in enumerate(meas, 1):
            fm, sm = load_touchstone(p)
            fr, dp = f_min(fm, sm)
            fres.append(fr)
            depth.append(dp)
            ax.plot(fm / 1e6, sm, lw=1, alpha=0.85, label=f"measurement {k}")
        fres = np.array(fres)
        f_meas = (fres.mean(), fres.std(ddof=1) if len(fres) > 1 else 0.0, float(np.mean(depth)))

    # cable experiments
    if rf is not None:
        for p in sorted(glob.glob(f"data/{s}_meas_[a-z]*.s1p")):
            fm, sm = load_touchstone(p)
            tag = p.split("_meas_")[1].removesuffix(".s1p")
            fr, _ = f_min(fm, sm)
            notes.append(f"{s}, {tag}: f = {fr/1e6:.1f} MHz")
            ax.plot(fm / 1e6, sm, ":", lw=1.2, label=tag)

    ax.axhline(VSWR2_DB, color="grey", ls=":", lw=1)
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("|S11| (dB)")
    ax.set_title(f"Dipole, setup {s}")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"plots/{s}_s11.png", dpi=200)
    plt.close(fig)

    bw = bandwidth(fh, sh)
    rows.append((s, f_sim, sh.min(), r_from_depth(sh.min()), f_nec, f_meas, bw))

print("| Setup | f HFSS (MHz) | f NEC2 (MHz) | f measured (MHz) | |S11| min (dB) | "
      "R from depth (ohm) | BW VSWR<2 (MHz) | HFSS vs NEC (%) | HFSS vs meas (%) |")
print("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
for s, f_sim, s_min, r_in, f_nec, f_meas, bw in rows:
    nec_c = f"{f_nec/1e6:.2f}" if f_nec else "-"
    d_nec = f"{100*(f_sim-f_nec)/f_nec:+.2f}" if f_nec else "-"
    if f_meas:
        mean, sd, dp = f_meas
        meas_c = f"{mean/1e6:.2f} ± {sd/1e6:.2f}"
        d_meas = f"{100*(f_sim-mean)/mean:+.2f}"
        depth_c = f"{dp:.2f}"
    else:
        meas_c, d_meas, depth_c = "-", "-", f"{s_min:.2f}"
    bw_c = f"{(bw[1]-bw[0])/1e6:.1f}" if bw else "-"
    print(f"| {s} | {f_sim/1e6:.2f} | {nec_c} | {meas_c} | {depth_c} | {r_in:.1f} | "
          f"{bw_c} | {d_nec} | {d_meas} |")

if notes:
    print("\nCable experiment:")
    for n in notes:
        print(" ", n)
