"""
Analysis for the half-wave dipole project.

Reads, for each setup A-D:
  data/<S>_hfss.csv          HFSS export      (unit taken from the header)
  data/<S>_nec.csv           NEC2 export      (written by sim_dipole_nec.py)
  data/<S>_meas_[0-9].s1p    measurement repeats, Touchstone
  data/<S>_meas_<tag>.s1p    variations: indoor, ferrite, hand, parallel

Writes plots/<S>_s11.png, plots/all_setups.png, data/results_summary.csv and
prints the markdown tables that go into the README.

    pip install numpy matplotlib scikit-rf
    python analyze_dipole.py

Two metrics are reported for every trace, because they disagree when the trace
is not clean:

  f_min     frequency of the deepest point, parabolically interpolated
  f_centre  centre of the contiguous band where VSWR < 2

f_min is the sharper number but it jumps between dips when two are nearly
equal. f_centre is blunt but stable. Where the two agree the trace is a single
clean resonance; where they do not, something else is in the measurement.
"""
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np

try:
    import skrf as rf
except ImportError:                       # the simulations plot fine without it
    rf = None

C = 299792458.0
UNITS = {"Hz": 1, "kHz": 1e3, "MHz": 1e6, "GHz": 1e9}
SETUPS = ["A", "B", "C", "D"]
VSWR2_DB = -9.54                          # |S11| for VSWR = 2, i.e. |gamma| = 1/3


def load_touchstone(path):
    if rf is not None:
        n = rf.Network(path)
        return n.f, n.s[:, 0, 0]
    f, g = [], []                         # minimal reader, real/imag Touchstone
    for line in open(path):
        line = line.strip()
        if not line or line[0] in "!#":
            continue
        a = line.split()
        f.append(float(a[0]))
        g.append(float(a[1]) + 1j * float(a[2]))
    return np.array(f), np.array(g)


def load_csv(path):
    """HFSS/NEC2 csv. The frequency unit is read from the header, not assumed."""
    header = open(path).readline()
    scale = UNITS[re.search(r"\[(\w+)\]", header).group(1)]
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    return d[:, 0] * scale, d[:, 1]


def f_min(f, s_db):
    """Frequency of minimum |S11|, refined by parabolic interpolation."""
    i = int(np.argmin(s_db))
    if 0 < i < len(f) - 1:
        y0, y1, y2 = s_db[i - 1:i + 2]
        den = y0 - 2 * y1 + y2
        if den != 0:
            return f[i] + 0.5 * (y0 - y2) / den * (f[i + 1] - f[i]), y1
    return f[i], s_db[i]


def match_band(f, s_db, thr=VSWR2_DB):
    """
    Widest contiguous run below thr, with linearly interpolated edges.

    Returns (lo, hi, truncated). truncated is True when the run touches an end
    of the sweep, which means the real band is wider than what is reported.
    """
    below = s_db < thr
    if not below.any():
        return None

    runs, start = [], None
    for i, v in enumerate(below):
        if v and start is None:
            start = i
        elif not v and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(below) - 1))

    a, b = max(runs, key=lambda r: f[r[1]] - f[r[0]])

    def cross(i, j):
        return f[i] + (thr - s_db[i]) * (f[j] - f[i]) / (s_db[j] - s_db[i])

    lo = f[a] if a == 0 else cross(a - 1, a)
    hi = f[b] if b == len(f) - 1 else cross(b, b + 1)
    return lo, hi, (a == 0 or b == len(f) - 1)


def r_from_depth(s_db_min):
    """Input resistance implied by the depth, assuming X = 0 at the minimum."""
    g = 10 ** (s_db_min / 20)
    return 50 * (1 + g) / (1 - g)


def line_length(f, gamma):
    """
    Electrical length of transmission line left in the reference plane.

    A lossless line does not change |S11|, only its phase, at a rate set by the
    round trip. Taking the median step keeps the estimate away from the
    resonance, where the antenna's own reactance dominates the rotation.
    """
    step = np.angle(gamma[1:] / gamma[:-1])
    return float(-np.median(step) / (2 * np.pi) * C / np.median(np.diff(f)) / 2)


def describe(f, s_db, gamma=None):
    fm, depth = f_min(f, s_db)
    band = match_band(f, s_db)
    d = dict(f_min=fm, depth=depth, f_centre=np.nan, bw=np.nan, q=np.nan,
             truncated=False, cable=np.nan)
    if band:
        lo, hi, trunc = band
        d.update(f_centre=0.5 * (lo + hi), bw=hi - lo,
                 q=0.707 / ((hi - lo) / (0.5 * (lo + hi))), truncated=trunc)
    if gamma is not None:
        d["cable"] = line_length(f, gamma)
    return d


# --------------------------------------------------------------- per setup --
os.makedirs("plots", exist_ok=True)
rows, variations = [], []
fig_all, axes_all = plt.subplots(2, 2, figsize=(11, 7))

for s, ax_all in zip(SETUPS, axes_all.flat):
    hfss = glob.glob(f"data/{s}_hfss.csv")
    nec = glob.glob(f"data/{s}_nec.csv")

    fig, ax = plt.subplots(figsize=(7, 4.2))
    row = dict(setup=s)

    for target in (ax, ax_all):
        if hfss:
            fh, sh = load_csv(hfss[0])
            target.plot(fh / 1e6, sh, "C7", lw=1.4, label="HFSS, nominal geometry")
        if nec:
            fn, sn = load_csv(nec[0])
            target.plot(fn / 1e6, sn, "C0--", lw=1.6, label="NEC2, as built")

    if hfss:
        row["hfss"] = describe(*load_csv(hfss[0]))
    if nec:
        row["nec"] = describe(*load_csv(nec[0]))

    repeats = sorted(glob.glob(f"data/{s}_meas_[0-9].s1p"))
    stats = []
    for k, path in enumerate(repeats, 1):
        fm_, g = load_touchstone(path)
        sm_ = 20 * np.log10(np.abs(g))
        stats.append(describe(fm_, sm_, g))
        for target in (ax, ax_all):
            target.plot(fm_ / 1e6, sm_, "C3", lw=1.8,
                        label="measured" if k == 1 else None)
    if stats:
        row["meas"] = stats

    labelled = set()
    for path in sorted(glob.glob(f"data/{s}_meas_[a-z]*.s1p")):
        fm_, g = load_touchstone(path)
        sm_ = 20 * np.log10(np.abs(g))
        tag = path.split("_meas_")[1].removesuffix(".s1p")
        group = tag.rstrip("0123456789")            # indoor1..indoor5 share a label
        variations.append((s, tag, describe(fm_, sm_, g)))
        ax.plot(fm_ / 1e6, sm_, "C8", lw=0.9, alpha=0.55,
                label=None if group in labelled else group)
        labelled.add(group)

    for target in (ax, ax_all):
        target.axhline(VSWR2_DB, color="grey", ls=":", lw=1)
        target.set_xlabel("Frequency (MHz)")
        target.set_ylabel("|S11| (dB)")
        target.grid(alpha=0.3)
        target.set_title(f"Setup {s}")
        target.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(f"plots/{s}_s11.png", dpi=200)
    plt.close(fig)
    rows.append(row)

fig_all.tight_layout()
fig_all.savefig("plots/all_setups.png", dpi=200)
plt.close(fig_all)

# ------------------------------------------------------------------ tables --
print("| Setup | HFSS nominal (MHz) | NEC2 as built (MHz) | measured f_min (MHz) | "
      "measured f_centre (MHz) | depth (dB) | BW VSWR<2 (MHz) | FBW (%) | Q | "
      "meas vs NEC2 (%) |")
print("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
for r in rows:
    h = r.get("hfss")
    n = r.get("nec")
    ms = r.get("meas")
    if ms:
        m = ms[0] if len(ms) == 1 else None
        fmin_c = (f"{m['f_min']/1e6:.2f}" if m else
                  f"{np.mean([x['f_min'] for x in ms])/1e6:.2f} ± "
                  f"{np.std([x['f_min'] for x in ms], ddof=1)/1e6:.2f}")
        cen = np.mean([x["f_centre"] for x in ms])
        cen_c = f"{cen/1e6:.2f}" + ("*" if any(x["truncated"] for x in ms) else "")
        dep = f"{np.mean([x['depth'] for x in ms]):.2f}"
        bw = f"{np.mean([x['bw'] for x in ms])/1e6:.1f}"
        fbw = f"{100*np.mean([x['bw'] for x in ms])/cen:.1f}"
        q = f"{np.mean([x['q'] for x in ms]):.1f}"
        dev = f"{100*(np.mean([x['f_min'] for x in ms]) - n['f_min'])/n['f_min']:+.2f}" if n else "-"
    else:
        fmin_c = cen_c = dep = bw = fbw = q = dev = "-"
    print(f"| {r['setup']} | {h['f_min']/1e6:.2f} | {n['f_min']/1e6:.2f} | {fmin_c} | "
          f"{cen_c} | {dep} | {bw} | {fbw} | {q} | {dev} |")

print("\nSimulated bandwidth for comparison:")
for r in rows:
    n = r.get("nec")
    if n and np.isfinite(n["bw"]):
        print(f"  {r['setup']}: NEC2 BW {n['bw']/1e6:.1f} MHz "
              f"({100*n['bw']/n['f_centre']:.1f}%), Q {n['q']:.1f}, depth {n['depth']:.2f} dB")

print("\nTransmission line left in the reference plane, from the phase slope:")
for r in rows:
    for k, m in enumerate(r.get("meas", []), 1):
        print(f"  {r['setup']} run {k}: {m['cable']:.2f} m one way")

if variations:
    print("\nVariations:")
    for s, tag, d in variations:
        print(f"  {s} / {tag}: f_min {d['f_min']/1e6:.2f} MHz, "
              f"f_centre {d['f_centre']/1e6:.2f} MHz, depth {d['depth']:.2f} dB")

# --------------------------------------------------------------- csv summary --
with open("data/results_summary.csv", "w") as fh:
    fh.write("setup,f_hfss_MHz,f_nec_MHz,f_meas_min_MHz,f_meas_centre_MHz,"
             "depth_meas_dB,bw_meas_MHz,fbw_meas_pct,q_meas,"
             "bw_nec_MHz,fbw_nec_pct,meas_vs_nec_pct,cable_m\n")
    for r in rows:
        h, n, ms = r.get("hfss"), r.get("nec"), r.get("meas")
        if not ms:
            continue
        fm = np.mean([x["f_min"] for x in ms])
        cen = np.mean([x["f_centre"] for x in ms])
        fh.write(f"{r['setup']},{h['f_min']/1e6:.3f},{n['f_min']/1e6:.3f},{fm/1e6:.3f},"
                 f"{cen/1e6:.3f},{np.mean([x['depth'] for x in ms]):.2f},"
                 f"{np.mean([x['bw'] for x in ms])/1e6:.2f},"
                 f"{100*np.mean([x['bw'] for x in ms])/cen:.2f},"
                 f"{np.mean([x['q'] for x in ms]):.2f},"
                 f"{n['bw']/1e6:.2f},{100*n['bw']/n['f_centre']:.2f},"
                 f"{100*(fm-n['f_min'])/n['f_min']:+.2f},"
                 f"{np.mean([x['cable'] for x in ms]):.2f}\n")
print("\nwrote data/results_summary.csv")
