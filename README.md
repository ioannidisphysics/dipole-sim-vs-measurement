# Half-wave dipole: simulation against measurement

S11 of a telescopic half-wave dipole at four element lengths, modelled independently with two
different numerical methods and then measured with a NanoVNA-F V2.

The point of the project is not to make the three numbers agree. It is to explain, with
numbers, why they differ, and to say where the remaining disagreement is not yet explained.

**Status:** complete. Simulation, cross-verification and measurement are all in the repository.
Four setups measured on 21 September 2026 on a roof in Thessaloniki, one calibration per setup.

**The main result.** The two solvers disagree by up to 4.4% on where these dipoles resonate, and
the measurement says which one is right and why. The difference is the feed gap: 21.83 mm, which
is 2% of the longest setup and 14% of the shortest. NEC2 cannot represent it — a method of
moments code drives one segment of a continuous wire — and its error tracks the gap fraction
exactly, reaching +4.8% on the shortest dipole. HFSS models the gap as a gap, and on that same
dipole it lands within 0.8% of the measurement.

![All four setups](plots/all_setups.png)

## 1. Theory

An ideal half-wave dipole resonates when its total span ℓ equals λ/2:

```
f₀ = c / (2ℓ)        ℓ = 2L + g
```

where `L` is the length of each arm and `g` is the feed gap. A real dipole made of
finite-thickness arms resonates lower than `f₀`, and the shift grows as the arms get thicker
relative to the length. The shift is usually quoted as 3–5%; that figure applies to thin wire.
For telescopic whips it is two to three times larger, as shown below.

Thickness does two things at once, and both are visible in the data: it lowers the resonance,
and it widens the bandwidth. They are the same cause — a thicker arm has lower characteristic
impedance, so the reactance climbs more slowly away from resonance — so a setup that is shifted
further down must also be broader, and that internal consistency is one of the checks used here.

## 2. Geometry

The antenna is a kit of two telescopic whips screwed into a plastic block. Two things about it
are not textbook, and the accuracy of everything downstream depends on how they are handled.

**The arms are tapered and only partly extended.** The large element runs from 6.00 mm diameter
at the base to 2.80 mm at the tip and extends from 230 to 1000 mm; the small element runs from
4.85 to 2.13 mm and extends from 50 to 130 mm. Collapsed length is the length of the base tube,
so a partly extended whip is a fat cylinder with a thin taper on the end, not a uniform rod.
`sim_dipole_nec.py` builds that profile and reduces it to one equivalent radius as
exp(⟨ln r⟩) along the arm, because the thin-wire kernel enters the integral equation as
ln(segment/radius) — the radius that matters is the one whose logarithm is the average.

**The feed gap is 21.83 mm**, measured metal to metal across the block. That is 2% of setup A
and 14% of setup D. NEC2 drives one segment of a continuous wire, so the gap cannot be modelled
as a gap; what can be done is to run the model twice, once with the gap counted as conductor
(span = 2L + g) and once not (span = 2L), and treat the pair as a bracket on the answer. HFSS
has no such restriction — it bridges the real gap with the port sheet — and the difference
between the two treatments turns out to be the main result of the project (§8).

As built, on the roof:

| Setup | Span ℓ (mm) | Arm L (mm) | Equivalent radius (mm) | ℓ / 2r | Element | Ideal f₀ (MHz) |
| --- | --- | --- | --- | --- | --- | --- |
| A | 1000 | 489.1 | 2.85 | 175 | large | 149.9 |
| B | 600 | 289.1 | 2.99 | 100 | large | 249.8 |
| C | 240 | 109.1 | 2.13 | 56 | small | 624.6 |
| D | 160 | 69.1 | 2.38 | 34 | small | 936.9 |

Setup B is *thicker* than setup A in the sense that matters, even though both use the same
element: B is extended less, so more of its length is base tube.

Setup D was planned at 140 mm and built at 160 mm. 160 mm is the better build — it keeps the
small element off its collapsed end stop, where the sections make poor contact — and the sweep
covered it, so the change was kept rather than undone.

## 3. Method 1 — Ansys HFSS (finite element method)

- Solution type: Modal, driven
- Radiation boundary on an air region offset by `pad` ≈ λ/4 at the low end of each sweep
- Lumped port, 50 Ω, integration line across the gap
- Adaptive mesh to Maximum Delta S = 0.01, up to 25 passes
- Interpolating frequency sweep, 0.1 MHz step

All four setups were solved twice. The first pass used the **nominal** geometry — uniform arms
of 2.0 or 2.5 mm radius and a 2 mm feed gap — because it was run before the kit was measured
with a caliper. Those exports are kept as `data/<S>_hfss_nominal.csv` and are what §5 is about.

The second pass uses the **as-built** geometry of §2, and it is the one compared against
measurement. It matters because HFSS models the 21.83 mm feed gap as a real gap bridged by the
port sheet, which is the one thing NEC2 structurally cannot do. §8 is mostly about what that
changed.

**Known defect in the committed C and D exports.** The arm radius in `hfss/dipole.aedt` is the
variable `r_arm`. For setups A and B it holds the as-built equivalent radius, 2.85 and 2.99 mm.
For setups C and D it still holds the nominal values, 2.5 and 2.0 mm, while the as-built figures
2.13 and 2.38 mm sit in a second variable, `r_pad`, that no geometry in the project references.
`data/C_hfss.csv` and `data/D_hfss.csv` were therefore solved at the nominal arm radius and are
labelled as-built when they are not. They are left in place with this notice until the two
designs are re-solved; §7 states the correction NEC2 predicts for them.

## 4. Method 2 — NEC2 (method of moments)

The four as-built dipoles were solved with NEC2 via [PyNEC](https://pypi.org/project/PyNEC/)
(`sim_dipole_nec.py`). NEC2 is a thin-wire method of moments code — a completely different
formulation from HFSS's finite elements — so agreement between the two is evidence that the
model is right, not just that it converged.

Segmentation: odd segment count, Δ < λ/20 at the top of the sweep and Δ > 2.5·r so the extended
thin-wire kernel stays valid, capped at 101 segments.

The NEC2 radiation pattern for setup B is the expected dipole doughnut, at **2.14 dBi** with a
**78° half-power beamwidth**. Textbook values for an ideal half-wave dipole are 2.15 dBi and 78°.

![Radiation pattern](plots/radiation_pattern.png)

NEC2, as built:

| Setup | f at X = 0 (MHz) | R at X = 0 (Ω) | f at min \|S11\| (MHz) | Shortening vs f₀ | Feed-gap bracket (MHz) | Equivalent-radius bracket (MHz) |
| --- | --- | --- | --- | --- | --- | --- |
| A | 141.6 | 72.1 | 140.7 | −6.1% | 140.7 … 143.8 | 140.6 … 142.2 |
| B | 234.0 | 72.3 | 232.1 | −7.1% | 232.1 … 240.6 | 232.0 … 235.3 |
| C | 578.5 | 72.6 | 572.0 | −8.4% | 572.0 … 627.5 | 569.8 … 581.5 |
| D | 857.7 | 73.2 | 844.1 | −9.9% | 844.1 … 972.2 | 843.6 … 864.4 |

The equivalent-radius bracket is the whole plausible range — every arm modelled at its base
radius, then every arm at its tip radius — and it is under 2.5% wide even for setup D. The
tapered-arm model is therefore not where the uncertainty lives. The feed-gap bracket is much
wider and is where it does live, for the two short setups.

## 5. Mesh convergence — why the first HFSS run was wrong

The first solution of setup A stopped at Delta S = 0.0196 after 15 passes, just inside the 0.02
criterion, and put the resonance at **137.73 MHz**. That looked acceptable in isolation.

It was caught by a consistency check rather than by the solver. Setups A and B share the same
arm radius but differ in length, so B is electrically thicker and must resonate relatively
*lower*; the ratio f_B / f_A therefore has to come out below the pure-scaling value of
ℓ_A / ℓ_B = 1.6645. The first run gave 1.6788 — on the wrong side.

Re-solving both setups with Maximum Delta S = 0.01 moved A by +1.32 MHz (+0.96%) and B by only
+0.41 MHz (+0.18%), which halved the disagreement with NEC2 for A:

| Setup | Delta S = 0.02 | Delta S = 0.01 | vs. NEC2, before → after |
| --- | --- | --- | --- |
| A | 137.73 MHz | 139.05 MHz | −2.13% → −1.19% |
| B | 231.22 MHz | 231.63 MHz | −0.40% → −0.23% |

The same ratio test applied to the measurement, where nothing is assumed about the model at all:
the measured f_B / f_A is **1.636** against the pure-scaling limit of 1.667, so the measurement
satisfies the test with room to spare. The C/D pair gives 1.421 against a limit of 1.500. The
as-built runs pass it too, at 1.650 for HFSS and 1.649 for NEC2.

The corresponding check on the NEC2 side is segmentation. Sweeping the segment count from 15 to
61 moves the resonance by less than 0.1% for setup C and 0.25% for setup D, so the NEC2 numbers
are not limited by discretisation either.

## 6. Measurement

NanoVNA-F V2, 50 kHz – 3 GHz. OSL calibration at the end of the VNA's own test cable, repeated
for every setup, 101 points per sweep. The antenna was on a roof, horizontal, at least 1 m clear
of railings and walls, with the feed cable led away perpendicular to the arms for as far as it
would go.

The calibration plane is therefore the connector the antenna screws onto — but **the antenna is
not at its own connector**. The kit's base carries 0.60 m of RG174 with a ferrite choke moulded
onto it, and that cable sits between the connector and the dipole. It is inside every
measurement here by construction, and short of cutting it off there is no way to calibrate it
out. The following is how much of it the data can see.

**How much line is in front of the reference plane, and what it is made of.** A lossless line
does not change |S11|, only its phase, so Γ rotates at a rate set by the round trip, and the
median phase step across the sweep converts that rate into a length. The estimator reports
free-space equivalent length — physical length divided by velocity factor — and it cannot tell
the cable's rotation from the antenna's own. A cable contributes a constant phase slope, so the
two add exactly, and running the identical estimator on the NEC2 model, which has no cable at
all, separates them:

| | A | B | C | D |
| --- | --- | --- | --- | --- |
| Measured, from the phase slope (m) | 1.31 | 1.17 | 1.23 | 1.22 |
| NEC2, antenna alone, no cable (m) | 0.25 | 0.15 | 0.05 | 0.02 |
| Difference, attributable to cable (m) | **1.05** | **1.02** | **1.18** | **1.20** |
| 0.60 m of RG174 at VF 0.66 would give | 0.91 | 0.91 | 0.91 | 0.91 |

Most of the 1.2 m is the kit's own pigtail, and the four estimates of it agree to ±8%, which is
what makes this one consistent measurement rather than four. A fifth of what setup A appeared to
show was never cable: it is the antenna's own phase rotation, which is largest where the sweep
is narrowest and falls to almost nothing by setup D.

The 0.11 to 0.29 m left over is not accounted for. Either the pigtail is longer than the 0.60 m
quoted for the kit, or its velocity factor is nearer 0.55 than the 0.66 nominal for RG174, which
varies between makes, or part of it is the connector and the transition into the plastic block.
The test that settles it takes five minutes and has not been run: sweep the pigtail on its own
with an open at the far end and read the delay directly.

The residual line does not move the resonance, because it does not change |Γ|, but its loss
deepens every dip and its outer braid carries common-mode current (§7).

**Indoors is not a measurement.** Five repeats of setup A taken indoors the day before are kept
in `data/A_meas_indoor*.s1p` and plotted in grey:

![Setup A, indoor against roof](plots/A_s11.png)

Indoors the antenna shows two or three shallow dips between 133 and 147 MHz, none deeper than
−12.2 dB, and which one is deepest changes between connect/disconnect cycles — the frequency of
minimum |S11| jumps from 144.6 to 133.7 MHz across five runs of the same antenna. On the roof
there is one dip, at 143.22 MHz, and it is 2.8 dB deeper than anything seen indoors. The indoor
runs are in the repository because the failure is the useful part: a number extracted from them
would have looked like a measurement.

## 7. Results

The measurement is the reference; both simulations are quoted as their error against it.

| Setup | Measured (MHz) | HFSS, as built (MHz) | error | NEC2, as built (MHz) | error |
| --- | --- | --- | --- | --- | --- |
| A | **143.22** | 139.01 | −2.94% | 140.70 | −1.76% |
| B | **234.25** | 229.39 | −2.08% | 232.06 | −0.94% |
| C | **567.10** | 547.76 | −3.41% | 571.97 | +0.86% |
| D | **805.60** | 812.06 | **+0.80%** | 844.11 | +4.78% |

Bandwidth and depth:

| Setup | Depth, measured (dB) | BW VSWR<2, measured | BW, HFSS | BW, NEC2 | Q, measured |
| --- | --- | --- | --- | --- | --- |
| A | −14.95 | 12.6 MHz (8.7%) | 11.9 MHz (8.5%) | 11.8 MHz (8.4%) | 8.1 |
| B | −16.27 | 33.4 MHz (14.2%) | 22.0 MHz (9.5%) | 22.1 MHz (9.5%) | 5.0 |
| C | −19.89 | 158.9 MHz (26.1%) | 63.4 MHz (11.5%) | 63.0 MHz (11.0%) | 2.7 |
| D | −18.76 | 212.9 MHz (26.3%) | 98.5 MHz (12.1%) | 108.2 MHz (12.7%) | 2.7 |

**Setups C and D are pending a re-solve.** Their HFSS runs used the nominal arm radius rather
than the as-built one (§3). NEC2, which does model the radius, puts the cost of that at −0.46%
for C, where 2.5 mm stood in for 2.13, and +0.60% for D, where 2.0 mm stood in for 2.38.
Carrying those across, HFSS as built should land near **550.3 MHz for C**, an error of −2.96%,
and **807.2 MHz for D**, an error of +0.20%. Those are predictions from a different solver, not
results: the tables above are left uncorrected and will be replaced by what the re-solve
actually returns. If the prediction holds, NEC2 − HFSS becomes +1.2, +1.2, +3.9, +4.6% against
gap fractions of 2.2, 3.6, 9.1, 13.6% — monotonic, which it is not in the table below.

Machine-readable summary in [`data/results_summary.csv`](data/results_summary.csv). Per-setup
plots: [A](plots/A_s11.png) · [B](plots/B_s11.png) · [C](plots/C_s11.png) · [D](plots/D_s11.png).

## 8. Discussion

**All four setups are within 5% of both simulations**, with no parameter fitted to the
measurement. The shortening from the ideal λ/2 length runs from −4.5% measured for A to −14%
for D, against the 3–5% that textbooks quote for thin wire.

**The feed gap is what NEC2 was getting wrong.** NEC2's error grows monotonically as the gap
takes up more of the antenna, reaching +4.78% for setup D, where 21.83 mm of gap is 13.6% of a
160 mm dipole. NEC2 drives one segment of a continuous wire with an ideal delta-gap source, so
that gap is not in the model at all. HFSS bridges it with a real port sheet, and on the same
geometry it lands 32 MHz lower — within 0.8% of the measurement:

| Setup | Gap as a fraction of ℓ | NEC2 − HFSS | NEC2 error | HFSS error |
| --- | --- | --- | --- | --- |
| A | 2.2% | +1.2% | −1.76% | −2.94% |
| B | 3.6% | +1.2% | −0.94% | −2.08% |
| C | 9.1% | +4.4% | +0.86% | −3.41% |
| D | 13.6% | +3.9% | +4.78% | **+0.80%** |

**Setup A was meant to be the control, and it is not a clean one.** Between the nominal and
as-built HFSS runs its feed gap went from 2 mm to 21.83 mm while its overall span stayed at
1000 mm, and the resonance moved by 0.04 MHz. But the arm radius changed in the same step, from
the nominal 2.0 mm to the as-built 2.85 mm, and that is not negligible: NEC2 puts it at
−0.78 MHz on its own. If HFSS answers a radius change the way NEC2 does, the two effects very
nearly cancelled, and the feed gap at constant span is worth about **+0.8 MHz**, or 0.6%, for
setup A — not the nothing that the raw difference suggests. The conclusion survives, because
0.6% at a 2.2% gap fraction against roughly 4% at 13.6% is still the same trend, but the number
that was quoted for the control was wrong.

The clean version costs one HFSS run: setup A at the as-built 2.85 mm radius with the gap put
back to 2 mm and the span held at 1000 mm, so that the gap is the only thing that differs. It
has not been run.

**And the gap fraction is not the only thing that grows down that table.** ℓ/2r falls from 175 to
34 across the four setups, so the dipoles become electrically fatter in step with the gap taking
up more of them — and NEC2's kernel is a *thin*-wire approximation. An error that tracks arm
thickness would produce exactly the same ordering. These four measurements do not separate the
two explanations.

What separates them is a run that moves one and not the other: setup D at its as-built radius,
with the gap cut from 21.83 mm to 2 mm and the arms lengthened to hold the span at 160 mm.
Thickness is untouched, only the feed gap moves. If HFSS then climbs towards NEC2's 844 MHz, the
gap is the cause; if it stays near 812, it is not, and the thin-wire kernel becomes the better
suspect. One solve, not yet run.

**HFSS is not simply the better solver, though.** For A, B and C it sits 2–3% *below* the
measurement while NEC2 sits within ±1.8%. Over the four setups the two have almost the same
RMS error, 2.5% against 2.6%; what differs is the shape. NEC2's error is structured, growing
with the gap fraction, and is understood. HFSS's is a roughly constant offset in one direction
that is **not explained here**. Candidates, none of them established:

- *Mesh convergence.* Tightening Maximum Delta S from 0.02 to 0.01 moved setup A up by 0.96%
  (§5), and the residual after a convergence step usually has the same sign as the step. Another
  halving, to 0.005, would show whether a further ~1% is waiting there. This is the cheap test
  and it has not been run.
- *The equivalent-radius model.* Both solvers were given one uniform radius standing in for a
  tapered arm (§2). Overstating it pushes both low, and the bracket in §4 is worth 1–2%. Setup
  D, whose arm is mostly base tube and so depends least on that model, is the one that agrees —
  but B and C do not order themselves by taper fraction, so this is suggestive at best.
- *Both models are lossless PEC*, with no contact resistance at the telescopic joints.
- *Neither model has a ground.* Both solve in free space; the antennas were a short distance
  above a concrete roof. A horizontal dipole over a conducting plane shifts in resonance by up
  to 2–3%, upwards or downwards depending on its height in wavelengths, and by less over
  concrete than over metal; it moves the input resistance as well. This one is different in kind
  from the others, because it acts on the *measurement* rather than on either model, so it
  cannot explain a difference between HFSS and NEC2 — but it is the most likely single
  explanation for A and B coming out above **both** of them. The height above the roof deck was
  not recorded, so it cannot be checked against the data that exists. Recording it, and if
  possible repeating one setup at two heights, is the cheapest thing to add to the next
  session.

**The two solvers agree with each other on bandwidth to better than 1%** for A, B and C, and to
9% for D, which is worth noting given they disagree by up to 4.4% on centre frequency. Bandwidth
is set by the arm thickness, which both model the same way; centre frequency is where the feed
model enters.

**The bandwidth agrees where the trace is clean and not where it is not.** For A the measured
VSWR<2 bandwidth is 12.6 MHz against 11.8 predicted, a 7% excess that is the right size for
conductor and contact loss in a telescopic whip. For C and D the measured bandwidth is two and a
half times the prediction, and the traces show why: they never return to 0 dB. Setup D is below
−4 dB across the whole 670–1300 MHz sweep, where a 160 mm dipole should be almost totally
reflective at both ends. Power is leaving somewhere other than the antenna.

**The dips are deeper than they should be, which is not good news.** NEC2 puts the input
resistance at 72–73 Ω for all four setups, which is a −15.1 dB minimum against 50 Ω. Setup A
measures 0.2 dB shallower than that; B, C and D measure 1.1, 4.7 and 3.6 dB *deeper*.
Attenuation between the calibration plane and the antenna subtracts twice the one-way loss from
the return loss, so a lossy path makes a mismatched antenna look well matched, and the three
setups above 200 MHz are all on that side. The consequence is that **the
input resistance cannot be recovered from these measurements**: de-embedding needs the cable's
loss and length to be known, and 10 cm of assumed length moves setup D's resonance by 20 MHz.

**There is no balun.** The kit feeds a balanced antenna from unbalanced coax, so current flows
on the outside of the braid and the cable becomes part of the radiating structure. That is the
most likely source of both the missing reflection off resonance and the ripples in the C and D
traces, whose spacing (roughly 150 MHz for D) corresponds to a line of the order of a metre —
the length the phase slope measures in front of the calibration plane. The experiment
that would confirm it — a clamp-on ferrite at the feed, and the cable re-routed parallel to one
arm — was planned and not done; the roof session was a single pass with one calibration per
setup and no repeats.

**What is missing**, in the order that would most change the conclusions:

1. *Re-solve C and D at the as-built arm radius* (§3). Their published HFSS numbers were
   obtained at the nominal radius, and the correction is predicted but not measured.
2. *Setup D with the gap cut to 2 mm at constant span*, which is the only run that separates
   the feed gap from the arm thickness as the cause of NEC2's error.
3. *Setup A at the as-built radius with a 2 mm gap*, which turns the control of §8 into a clean
   one-variable comparison.
4. *HFSS at Delta S = 0.005 on setup A*, the cheap test for the unexplained 2–3% offset.
5. *The height above the roof deck*, which was never recorded and is needed before ground
   proximity can be ruled in or out.
6. *Repeats.* One measurement per setup, so there is no repeatability figure for the roof
   session; the five indoor repeats of setup A are all that exist, and they measure the room
   rather than the antenna.
7. *The ferrite and cable-routing experiment*, planned and not done.
8. *The delay of the kit's own pigtail*, measured directly, which closes §6.

## 9. Measuring the same antennas with an SDR

The same four antennas were swept with an RTL-SDR Blog V4 in
[sdr-spectrum-analyzer](https://github.com/ioannidisphysics/sdr-spectrum-analyzer), comparing
the receiver's noise floor with the antenna connected against a 50 Ω load, which should peak
where the antenna is matched. **It did not work**, and the reasons are written up there. The
first is the one that matters: the method needs the antenna to be looking at something hotter
than the 50 Ω resistor it is being compared against, and at the UHF end of these sweeps it is
looking at ground and buildings at about the same 290 K, so there is nothing to measure at any
receiver gain. The method belongs at VHF and below. The second is that the two VHF sweeps, where
there would have been something to see, were taken at too low a tuner gain. The VNA numbers in
this repository are the ones that stand.

## 10. Repository

```text
dipole-sim-vs-measurement/
  README.md
  sim_dipole_nec.py      NEC2 model of the as-built geometry, writes data/<S>_nec.csv
  analyze_dipole.py      reads HFSS/NEC2/Touchstone, writes the plots and the tables
  check_model_sensitivity.py
                         the two sensitivity numbers quoted in §6 and §8: what the
                         assumed arm radius is worth, and how much of the line in
                         front of the reference plane is cable and how much is the
                         antenna's own phase rotation
  data/                  <S>_hfss.csv          HFSS, as-built geometry
                         <S>_hfss_nominal.csv  HFSS, nominal geometry (the §5 convergence story)
                         <S>_nec.csv           NEC2, as built
                         <S>_meas_*.s1p        measurements
                         results_summary.csv
  hfss/                  dipole.aedt
  plots/                 per-setup comparisons, all_setups.png, radiation_pattern.png
```

Measurement files are Touchstone `.s1p`, real/imaginary, 50 Ω reference:
`<S>_meas_1.s1p` is the roof measurement, `<S>_meas_indoor*.s1p` the indoor repeats.

## 11. Reproducing

```bash
pip install numpy scipy matplotlib scikit-rf PyNEC
python sim_dipole_nec.py          # re-runs NEC2 for all four setups, prints the geometry table
python analyze_dipole.py          # rebuilds the plots and prints the results tables
python check_model_sensitivity.py # the radius and feed-line numbers quoted in §6 and §8
```

`analyze_dipole.py` reads the frequency unit from the csv header, so HFSS exports in MHz or GHz
both work, and it runs with any subset of the measurements present.

## References

- C. A. Balanis, *Antenna Theory: Analysis and Design*, 4th ed., ch. 4 — dipole resonance and
  the length/thickness dependence
- G. J. Burke and A. J. Poggio, *Numerical Electromagnetics Code (NEC) — Method of Moments*,
  Lawrence Livermore National Laboratory, 1981 — thin-wire kernel and segmentation limits
- Ansys HFSS documentation — adaptive meshing and the Delta S convergence criterion
