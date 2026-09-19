"""
Προσομοίωση διπόλου λ/2 με NEC2 (μέθοδος ροπών για λεπτά σύρματα).
Ανεξάρτητη μέθοδος από το HFSS (FEM), για σύγκριση θεωρία / NEC / HFSS / μέτρηση.

Εγκατάσταση:  pip install PyNEC numpy matplotlib
Εκτέλεση:     python sim_dipole_nec.py

Όταν μετρήσεις το κιτ με το παχύμετρο, άλλαξε L, g, r στο SETUPS και ξανατρέξε.
"""
import numpy as np
import matplotlib.pyplot as plt
from PyNEC import nec_context

C = 299792458.0
Z0 = 50.0

# ---- Διαστάσεις (mm). ΥΠΟΘΕΣΕΙΣ μέχρι να μετρηθεί το κιτ ----
# r = μέση ακτίνα σκέλους (μέσος όρος βάσης και άκρης)
SETUPS = {
    "A": dict(L=500.0, g=2.0, r=2.5),   # μεγάλα τηλεσκοπικά
    "B": dict(L=300.0, g=2.0, r=2.5),
    "C": dict(L=120.0, g=2.0, r=2.0),   # μικρά τηλεσκοπικά
    "D": dict(L=70.0,  g=2.0, r=2.0),
}
N_FREQ = 301


def n_segments(length_m, radius_m, f_max):
    """Μονός αριθμός τμημάτων: Δ < λ/20 στη μέγιστη συχνότητα και Δ > 2.5r (extended kernel)."""
    lam_min = C / f_max
    n_min = int(np.ceil(length_m / (lam_min / 20)))
    n_max = int(np.floor(length_m / (2.5 * radius_m)))
    n = max(n_min, min(n_max, 101))
    return n if n % 2 else n + 1


def simulate(L_mm, g_mm, r_mm, f_start, f_stop, n_freq=N_FREQ, gain_at=None):
    """Επιστρέφει f (Hz), Zin (Ω) και προαιρετικά μέγιστο κέρδος (dBi) σε μία συχνότητα."""
    length = (2 * L_mm + g_mm) * 1e-3
    a = r_mm * 1e-3
    nseg = n_segments(length, a, f_stop)

    ctx = nec_context()
    geo = ctx.get_geometry()
    # Σύρμα κατά μήκος του z, κέντρο στην αρχή των αξόνων
    geo.wire(1, nseg, 0, 0, -length / 2, 0, 0, length / 2, a, 1.0, 1.0)
    ctx.geometry_complete(0)
    ctx.set_extended_thin_wire_kernel(True)
    ctx.gn_card(-1, 0, 0, 0, 0, 0, 0, 0)                     # ελεύθερος χώρος
    ctx.ex_card(0, 1, nseg // 2 + 1, 0, 1.0, 0, 0, 0, 0, 0)  # πηγή τάσης στο κεντρικό τμήμα

    f = np.linspace(f_start, f_stop, n_freq)
    step = (f[1] - f[0]) / 1e6 if n_freq > 1 else 0.0
    ctx.fr_card(0, n_freq, f_start / 1e6, step)
    ctx.xq_card(0)
    z = np.array([ctx.get_input_parameters(i).get_impedance()[0] for i in range(n_freq)])

    gain = None
    if gain_at is not None:
        c2 = nec_context()
        g2 = c2.get_geometry()
        g2.wire(1, nseg, 0, 0, -length / 2, 0, 0, length / 2, a, 1.0, 1.0)
        c2.geometry_complete(0)
        c2.set_extended_thin_wire_kernel(True)
        c2.gn_card(-1, 0, 0, 0, 0, 0, 0, 0)
        c2.ex_card(0, 1, nseg // 2 + 1, 0, 1.0, 0, 0, 0, 0, 0)
        c2.fr_card(0, 1, gain_at / 1e6, 0)
        c2.rp_card(0, 181, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)  # θ 0..180° ανά 1°, φ = 0
        gain = c2.get_gain_max(0)
    return f, z, nseg, gain


def s11_db(z):
    return 20 * np.log10(np.abs((z - Z0) / (z + Z0)))


def f_min(f, s_db):
    """Ελάχιστο |S11| με παραβολική παρεμβολή (ίδια μέθοδος με το analyze_dipole.py)."""
    i = int(np.argmin(s_db))
    if 0 < i < len(f) - 1:
        y0, y1, y2 = s_db[i - 1:i + 2]
        d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
        return f[i] + d * (f[i + 1] - f[i]), s_db[i]
    return f[i], s_db[i]


def f_x0(f, z):
    """Πρώτη συχνότητα όπου η φανταστική αντίσταση περνά από 0 (-→+)."""
    x = z.imag
    idx = np.where((x[:-1] < 0) & (x[1:] >= 0))[0]
    if idx.size == 0:
        return np.nan, np.nan
    i = idx[0]
    fr = f[i] - x[i] * (f[i + 1] - f[i]) / (x[i + 1] - x[i])
    return fr, np.interp(fr, f, z.real)


if __name__ == "__main__":
    import os
    os.makedirs("data", exist_ok=True)
    os.makedirs("plots", exist_ok=True)

    rows = []
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.5))
    for ax, (name, p) in zip(axes.flat, SETUPS.items()):
        ell = (2 * p["L"] + p["g"]) * 1e-3
        f0 = C / (2 * ell)
        f, z, nseg, _ = simulate(p["L"], p["g"], p["r"], 0.7 * f0, 1.3 * f0)
        s = s11_db(z)
        fmin, smin = f_min(f, s)
        fx, rx = f_x0(f, z)
        rows.append((name, p, f0, fx, rx, fmin, smin, nseg))

        np.savetxt(f"data/{name}_nec.csv", np.column_stack([f / 1e6, s, z.real, z.imag]),
                   delimiter=",", header="Freq [MHz],dB(S(1,1)),Re(Z) [ohm],Im(Z) [ohm]",
                   comments="", fmt="%.6f")

        ax.plot(f / 1e6, s, "C0")
        ax.axvline(f0 / 1e6, color="grey", ls=":", label=f"Ιδανικό c/2ℓ = {f0/1e6:.0f} MHz")
        ax.axvline(fmin / 1e6, color="C3", ls="--", lw=1, label=f"NEC ελάχιστο = {fmin/1e6:.1f} MHz")
        ax.set_title(f"Ρύθμιση {name}: L = {p['L']:.0f} mm, r = {p['r']} mm")
        ax.set_xlabel("Συχνότητα (MHz)"); ax.set_ylabel("|S11| (dB)")
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig("plots/nec_s11_all.png", dpi=200); plt.close(fig)

    print("| Ρύθμιση | L (mm) | r (mm) | ℓ/2r | f₀ ιδανικό (MHz) | f (X=0) NEC (MHz) | R στο X=0 (Ω) "
          "| f ελαχ. S11 NEC (MHz) | S11 ελάχ. (dB) | Μετατόπιση από f₀ (%) |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for name, p, f0, fx, rx, fmin, smin, nseg in rows:
        ell = 2 * p["L"] + p["g"]
        print(f"| {name} | {p['L']:.0f} | {p['r']} | {ell/(2*p['r']):.0f} | {f0/1e6:.0f} | {fx/1e6:.1f} "
              f"| {rx:.1f} | {fmin/1e6:.1f} | {smin:.1f} | {100*(fmin-f0)/f0:+.1f} |")

    # ---- Ευαισθησία: πόσο μετακινείται η ρύθμιση B με r και g ----
    print("\nΕυαισθησία ρύθμισης B (L = 300 mm):")
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    for r in [1.0, 2.0, 3.0, 4.0]:
        f, z, _, _ = simulate(300, 0, r, 175e6, 300e6)
        s = s11_db(z); fm, sm = f_min(f, s)
        ax[0].plot(f / 1e6, s, label=f"r = {r} mm ({fm/1e6:.1f} MHz)")
        print(f"  r = {r} mm, g = 0:  f = {fm/1e6:.1f} MHz, S11 = {sm:.1f} dB")
    for g in [0.0, 5.0, 10.0]:
        f, z, _, _ = simulate(300, g, 2.5, 175e6, 300e6)
        s = s11_db(z); fm, sm = f_min(f, s)
        ax[1].plot(f / 1e6, s, label=f"g = {g:.0f} mm ({fm/1e6:.1f} MHz)")
        print(f"  r = 2.5 mm, g = {g:.0f} mm:  f = {fm/1e6:.1f} MHz, S11 = {sm:.1f} dB")
    ax[0].set_title("B: επίδραση ακτίνας σκέλους"); ax[1].set_title("B: επίδραση κενού τροφοδοσίας")
    for a_ in ax:
        a_.set_xlabel("Συχνότητα (MHz)"); a_.set_ylabel("|S11| (dB)"); a_.grid(alpha=0.3); a_.legend(fontsize=8)
    fig.tight_layout(); fig.savefig("plots/nec_sensitivity_B.png", dpi=200); plt.close(fig)

    # ---- Κέρδος στον συντονισμό (έλεγχος ~2.15 dBi) ----
    name, p, f0, fx, rx, fmin, smin, nseg = rows[1]
    _, _, _, gmax = simulate(p["L"], p["g"], p["r"], fx, fx, n_freq=1, gain_at=fx)
    print(f"\nΜέγιστο κέρδος ρύθμισης B στα {fx/1e6:.1f} MHz: {gmax:.2f} dBi")
