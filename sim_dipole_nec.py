"""
Half-wave dipole modelled with NEC2 (method of moments, thin-wire kernel).

An independent second opinion on the HFSS (finite element) model, and the
simulation that is compared against the NanoVNA measurement.

    pip install PyNEC numpy matplotlib
    python sim_dipole_nec.py

Writes data/<S>_nec.csv and the plots under plots/.

Geometry
--------
The antenna is a kit of two telescopic whips screwed into a plastic block.
Two things about it are not textbook:

  the arms are tapered, not uniform cylinders, and only partly extended
  the feed gap is 21.8 mm, which is 14% of the total length in setup D

Both are handled explicitly below rather than swept under an assumed radius.
"""
import os

import numpy as np
import matplotlib.pyplot as plt
from PyNEC import nec_context

C = 299792458.0
Z0 = 50.0
N_FREQ = 401

# ---------------------------------------------------------------- geometry --
# Telescopic elements, measured with a caliper. L_min is the collapsed length
# (the base tube), L_max the fully extended length.
LARGE = dict(d_base=6.00, d_tip=2.80, L_min=230.0, L_max=1000.0)
SMALL = dict(d_base=4.85, d_tip=2.13, L_min=50.0, L_max=130.0)

FEED_GAP_MM = 21.83          # nearest metal to nearest metal, across the block

# As built on the roof, 2026-09-21. ell is the total tip-to-tip span.
SETUPS = {
    "A": dict(ell=1000.0, element=LARGE),
    "B": dict(ell=600.0, element=LARGE),
    "C": dict(ell=240.0, element=SMALL),
    "D": dict(ell=160.0, element=SMALL),   # built at 160, not the planned 140
}


def radius_profile(arm_mm, element, n=4001):
    """
    Radius along one arm, from the feed outwards, for a partly extended whip.

    The collapsed length is the base tube, so the first L_min millimetres keep
    the base diameter whatever the extension. Beyond that the thinner sections
    are out, and the diameter falls linearly to whatever the tip section is at
    this extension.
    """
    r_base, r_tip = element["d_base"] / 2, element["d_tip"] / 2
    x = np.linspace(0.0, arm_mm, n)
    if arm_mm <= element["L_min"]:
        return x, np.full(n, r_base)

    out = (arm_mm - element["L_min"]) / (element["L_max"] - element["L_min"])
    r_end = r_base + (r_tip - r_base) * out
    r = np.where(x <= element["L_min"], r_base,
                 r_base + (r_end - r_base) * (x - element["L_min"]) / (arm_mm - element["L_min"]))
    return x, r


def equivalent_radius(arm_mm, element):
    """
    One uniform radius standing in for the tapered arm.

    The thin-wire kernel enters the integral equation as ln(segment/radius), so
    the radius that matters is the one whose logarithm is the average, not the
    radius whose value is the average. In practice the two differ by less than
    0.1% here; the log form is used because it is the one that is justified.
    """
    x, r = radius_profile(arm_mm, element)
    return float(np.exp(np.trapezoid(np.log(r), x) / arm_mm))


def n_segments(length_m, radius_m, f_max):
    """Odd segment count with delta < lambda/20 at f_max and delta > 2.5r."""
    lam_min = C / f_max
    n_min = int(np.ceil(length_m / (lam_min / 20)))
    n_max = int(np.floor(length_m / (2.5 * radius_m)))
    n = max(n_min, min(n_max, 101))
    return n if n % 2 else n + 1


def simulate(span_mm, r_mm, f_start, f_stop, n_freq=N_FREQ, gain_at=None):
    """
    S11 of a centre-fed cylinder of total length span_mm and radius r_mm.

    The feed gap is not modelled as a gap: NEC drives one segment of a
    continuous wire with an ideal delta-gap source. span_mm is therefore the
    electrical span assumed for the run, and the effect of the real 21.8 mm
    gap is bracketed by running both 2L+g and 2L (see __main__).
    """
    length = span_mm * 1e-3
    a = r_mm * 1e-3
    nseg = n_segments(length, a, f_stop)

    ctx = nec_context()
    ctx.get_geometry().wire(1, nseg, 0, 0, -length / 2, 0, 0, length / 2, a, 1.0, 1.0)
    ctx.geometry_complete(0)
    ctx.set_extended_thin_wire_kernel(True)
    ctx.gn_card(-1, 0, 0, 0, 0, 0, 0, 0)                     # free space
    ctx.ex_card(0, 1, nseg // 2 + 1, 0, 1.0, 0, 0, 0, 0, 0)  # source on the centre segment

    f = np.linspace(f_start, f_stop, n_freq)
    step = (f[1] - f[0]) / 1e6 if n_freq > 1 else 0.0
    ctx.fr_card(0, n_freq, f_start / 1e6, step)
    ctx.xq_card(0)
    z = np.array([ctx.get_input_parameters(i).get_impedance()[0] for i in range(n_freq)])

    gain = None
    if gain_at is not None:
        c2 = nec_context()
        c2.get_geometry().wire(1, nseg, 0, 0, -length / 2, 0, 0, length / 2, a, 1.0, 1.0)
        c2.geometry_complete(0)
        c2.set_extended_thin_wire_kernel(True)
        c2.gn_card(-1, 0, 0, 0, 0, 0, 0, 0)
        c2.ex_card(0, 1, nseg // 2 + 1, 0, 1.0, 0, 0, 0, 0, 0)
        c2.fr_card(0, 1, gain_at / 1e6, 0)
        c2.rp_card(0, 181, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)  # theta 0..180 deg, phi = 0
        gain = c2.get_gain_max(0)
    return f, z, nseg, gain


def s11_db(z):
    return 20 * np.log10(np.abs((z - Z0) / (z + Z0)))


def f_min(f, s_db):
    """Minimum of |S11|, refined by parabolic interpolation through three points."""
    i = int(np.argmin(s_db))
    if 0 < i < len(f) - 1:
        y0, y1, y2 = s_db[i - 1:i + 2]
        den = y0 - 2 * y1 + y2
        if den != 0:
            return f[i] + 0.5 * (y0 - y2) / den * (f[i + 1] - f[i]), y1
    return f[i], s_db[i]


def f_x0(f, z):
    """First frequency where the input reactance crosses zero from below."""
    x = z.imag
    idx = np.where((x[:-1] < 0) & (x[1:] >= 0))[0]
    if idx.size == 0:
        return np.nan, np.nan
    i = idx[0]
    fr = f[i] - x[i] * (f[i + 1] - f[i]) / (x[i + 1] - x[i])
    return fr, np.interp(fr, f, z.real)


if __name__ == "__main__":
    os.makedirs("data", exist_ok=True)
    os.makedirs("plots", exist_ok=True)

    rows = []
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.5))
    for ax, (name, p) in zip(axes.flat, SETUPS.items()):
        ell = p["ell"]
        arm = (ell - FEED_GAP_MM) / 2
        r = equivalent_radius(arm, p["element"])
        f0 = C / (2 * ell * 1e-3)

        f, z, nseg, _ = simulate(ell, r, 0.60 * f0, 1.35 * f0)
        s = s11_db(z)
        fm, sm = f_min(f, s)
        fx, rx = f_x0(f, z)

        # Bracket: the same wire with the feed gap not counted as conductor.
        f2, z2, _, _ = simulate(2 * arm, r, 0.60 * f0, 1.45 * f0)
        fm_nogap, _ = f_min(f2, s11_db(z2))

        rows.append(dict(name=name, ell=ell, arm=arm, r=r, nseg=nseg, f0=f0,
                         fm=fm, sm=sm, fx=fx, rx=rx, fm_nogap=fm_nogap))

        np.savetxt(f"data/{name}_nec.csv",
                   np.column_stack([f / 1e6, s, z.real, z.imag]),
                   delimiter=",", header="Freq [MHz],dB(S(1,1)),Re(Z) [ohm],Im(Z) [ohm]",
                   comments="", fmt="%.6f")

        ax.plot(f / 1e6, s, "C0")
        ax.axvline(f0 / 1e6, color="grey", ls=":", label=f"ideal c/2l = {f0/1e6:.0f} MHz")
        ax.axvline(fm / 1e6, color="C3", ls="--", lw=1, label=f"NEC2 = {fm/1e6:.1f} MHz")
        ax.axvspan(fm / 1e6, fm_nogap / 1e6, color="C3", alpha=0.10, label="feed-gap bracket")
        ax.set_title(f"Setup {name}: span {ell:.0f} mm, r_eq {r:.2f} mm")
        ax.set_xlabel("Frequency (MHz)")
        ax.set_ylabel("|S11| (dB)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("plots/nec_s11_all.png", dpi=200)
    plt.close(fig)

    print("| Setup | span (mm) | arm (mm) | r_eq (mm) | span/2r | segs | ideal f0 (MHz) | "
          "f(X=0) (MHz) | R at X=0 (ohm) | NEC2 f (MHz) | NEC2 f, gap not counted (MHz) | shortening (%) |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for d in rows:
        print(f"| {d['name']} | {d['ell']:.0f} | {d['arm']:.1f} | {d['r']:.2f} | "
              f"{d['ell']/(2*d['r']):.0f} | {d['nseg']} | {d['f0']/1e6:.1f} | {d['fx']/1e6:.1f} | "
              f"{d['rx']:.1f} | {d['fm']/1e6:.1f} | {d['fm_nogap']/1e6:.1f} | "
              f"{100*(d['fm']-d['f0'])/d['f0']:+.1f} |")

    # ---- How much does the assumed radius matter? Bracket it tip to base ----
    print("\nSensitivity to the equivalent radius (tip radius vs base radius):")
    for d in rows:
        el = SETUPS[d["name"]]["element"]
        edges = {}
        for tag, r_try in (("fat", el["d_base"] / 2), ("thin", el["d_tip"] / 2)):
            f_, z_, _, _ = simulate(d["ell"], r_try, 0.60 * d["f0"], 1.35 * d["f0"])
            edges[tag], _ = f_min(f_, s11_db(z_))
        print(f"  {d['name']}: {edges['fat']/1e6:.1f} MHz (all base radius) .. "
              f"{edges['thin']/1e6:.1f} MHz (all tip radius), model says {d['fm']/1e6:.1f} MHz")

    # ---- Gain at resonance, as a sanity check against the textbook 2.15 dBi ----
    b = rows[1]
    _, _, _, gmax = simulate(b["ell"], b["r"], b["fx"], b["fx"], n_freq=1, gain_at=b["fx"])
    print(f"\nSetup B peak gain at {b['fx']/1e6:.1f} MHz: {gmax:.2f} dBi")
