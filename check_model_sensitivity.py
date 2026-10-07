"""
The two sensitivity checks quoted in the README, both reproducible from here.

    python check_model_sensitivity.py

1. How much the assumed arm radius moves the resonance.

   The HFSS project carries the as-built equivalent radius for setups A and B
   and the nominal one for C and D (README §3). NEC2 models the radius too, so
   running it at both values puts a number on what that defect costs, and on
   how much of setup A's "control" was really a radius change rather than a
   feed-gap change.

2. What the 1.2 m of line in front of the reference plane is made of.

   analyze_dipole.py reports a line length from the median phase slope of
   Gamma. That estimator cannot tell a cable's rotation from the antenna's
   own, and the two add exactly, because a cable's phase slope is constant
   with frequency. Running the identical estimator on the NEC2 model, which
   has no cable at all, separates them.

Needs the measurements in data/ and PyNEC installed.
"""

import numpy as np
import skrf as rf

from sim_dipole_nec import (
    C,
    FEED_GAP_MM,
    SETUPS,
    equivalent_radius,
    f_min,
    s11_db,
    simulate,
)

Z0 = 50.0

# Arm radius each HFSS design actually used, read out of hfss/dipole.aedt.
# A and B carry the as-built value; C and D were left at the nominal one.
R_ARM_IN_AEDT_MM = {"A": 2.85, "B": 2.99, "C": 2.50, "D": 2.00}

# The nominal radius each setup was first solved at, for the §8 control.
R_NOMINAL_MM = {"A": 2.00, "B": 2.50, "C": 2.50, "D": 2.00}

# The kit's permanently attached feed cable.
PIGTAIL_M = 0.60
VF_RG174 = 0.66


def line_length(f, gamma):
    """Verbatim from analyze_dipole.py, so the two agree by construction."""

    step = np.angle(gamma[1:] / gamma[:-1])
    return float(-np.median(step) / (2 * np.pi) * C / np.median(np.diff(f)) / 2)


def as_built(name):
    """Span, arm length and equivalent radius of one setup."""

    p = SETUPS[name]
    arm = (p["ell"] - FEED_GAP_MM) / 2
    return p["ell"], arm, equivalent_radius(arm, p["element"])


def resonance(span_mm, r_mm, f0):
    f, z, _, _ = simulate(span_mm, r_mm, 0.60 * f0, 1.35 * f0)
    return f_min(f, s11_db(z))[0]


def radius_sensitivity():
    print("1. Resonance against the assumed arm radius, NEC2\n")
    print(f"{'setup':>5} {'r as built':>11} {'f':>9} | {'r in .aedt':>11} {'f':>9} "
          f"{'shift':>8} {'shift %':>8} | {'r nominal':>10} {'shift %':>8}")
    print("-" * 92)

    for name in SETUPS:
        ell, _, r_eq = as_built(name)
        f0 = C / (2 * ell * 1e-3)

        f_ref = resonance(ell, r_eq, f0)
        f_aedt = resonance(ell, R_ARM_IN_AEDT_MM[name], f0)
        f_nom = resonance(ell, R_NOMINAL_MM[name], f0)

        print(f"{name:>5} {r_eq:11.2f} {f_ref/1e6:9.2f} | {R_ARM_IN_AEDT_MM[name]:11.2f} "
              f"{f_aedt/1e6:9.2f} {(f_aedt-f_ref)/1e6:+8.2f} "
              f"{100*(f_aedt-f_ref)/f_ref:+8.2f} | {R_NOMINAL_MM[name]:10.2f} "
              f"{100*(f_nom-f_ref)/f_ref:+8.2f}")

    print("\n   C and D are the two whose HFSS exports carry the .aedt radius, so the"
          "\n   'shift %' column is what re-solving them should move the published"
          "\n   number by, with the sign reversed.")
    print("\n   For setup A the 'r nominal' column is the part of the nominal-to-as-built"
          "\n   step that was a radius change and not a feed-gap change (README §8).")


def feed_line_decomposition():
    print("\n\n2. What the line in front of the reference plane is made of\n")
    print(f"   {PIGTAIL_M:.2f} m of RG174 at VF {VF_RG174} is "
          f"{PIGTAIL_M/VF_RG174:.2f} m of free-space equivalent length\n")
    print(f"{'setup':>5} {'sweep (MHz)':>17} {'measured':>9} {'antenna':>9} "
          f"{'cable':>8} {'unexplained':>12}")
    print("-" * 68)

    for name in SETUPS:
        try:
            net = rf.Network(f"data/{name}_meas_1.s1p")
        except Exception as exc:                      # noqa: BLE001
            print(f"{name:>5}  no measurement: {exc}")
            continue

        f_meas, g_meas = net.f, net.s[:, 0, 0]
        measured = line_length(f_meas, g_meas)

        ell, _, r_eq = as_built(name)
        f_sim, z, _, _ = simulate(ell, r_eq, f_meas[0], f_meas[-1], n_freq=len(f_meas))
        antenna = line_length(f_sim, (z - Z0) / (z + Z0))

        cable = measured - antenna
        print(f"{name:>5} {f_meas[0]/1e6:8.1f}-{f_meas[-1]/1e6:8.1f} {measured:9.2f} "
              f"{antenna:9.2f} {cable:8.2f} {cable - PIGTAIL_M/VF_RG174:+12.2f}")

    print("\n   'antenna' is the same estimator run on the NEC2 model, which has no cable,"
          "\n   so it is the part of the reported length that was never line at all.")


if __name__ == "__main__":
    radius_sensitivity()
    feed_line_decomposition()
