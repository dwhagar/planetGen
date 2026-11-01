# Water delivery: the mass law of row b4 (GEN.510), replaced by the formation-distance law (GEN.802)

Status: GEN.510 built 2026-10-24 (lane "Source and tune atmosphere draws"); replaced 2026-10-25 by GEN.802, which wires one joint law into rows b1 to b4 (section 0). Data now `src/planetgen/physics/data/ngpps_formation_water.json` (the earlier `ngpps_water_fraction.json` is gone), tests `src/tests/test_water_delivery.py`. Decision D1 of `tunable-gaps.md` is followed: the formation model is a switch with the Bern population as the MODEL default and a flagged bracket.

## 0. GEN.802: one law in the formation distance and the mass

The mass logistic below (centre 6.1 Earth masses, slope 6.28 per dex) and the snow-line weights of b1 to b3 were two views of one process, and the gate `x >= 1` of b4 contradicted the population (`volatile-delivery-formation-distance.md`). The rows now draw their class from three nested events fitted jointly to the 10,506 planets of 0.5 to 30 Earth masses of the Bern condensed-water population (starting distance from DACE ng76): not dry (water mass fraction f of at least 0.00025), high or rich (at least 0.0016), rich (at least 0.1), each a logit `a + b ln x + c log10 M + h max(0, ln x - ln 0.6)` with ln x clipped to x in 0.03 to 4 and the mass to 0.5 to 30 Earth masses (coefficients in the data file; the separate straight fits of the research note cross each other and a shared-slope ordered logit is biased in the middle range, so the hinge at 0.6 was added: log likelihood -2557 against -2787 and -2767). The nested probabilities are held in order. b1 is the dry share, b2 the low share, b3 the high share and b4 the water-rich share times the plateau `c_b4` (the water-rich planets the plateau removes are counted as high). The mass moves the water-rich share by 0.13 per dex in the logit, so the old mass law is not applied on top. Not tuned to any outcome (tier 6 for the centre, 7 for the width and class split, tweak L0; registry row AT.b4); the thresholds that map f to Ortenzi's cases remain a design assumption. The recipes are not read by the generator yet: no generated value changes, no epoch.

Fit script: `reports/gen802-files/fit_ngpps_formation_water.py` in the project files (downloads the two inputs; the raw files never enter the repository).

Below is the record of GEN.510, whose mass law and data file are gone.

## 1. What changed (GEN.510, superseded)

Row b4 (a water-rich planet, migrated from beyond the snow line) weighed `c_b4 s^3` of the snow-line ratio alone, `c_b4` log-uniform from 0.05 to 0.4. It had no mass dependence, but the migration that brings the water grows steeply with mass. The Bern population synthesis (NGPPS, NG76 nominal, Zenodo 10527170; the group 2 worker of the tunable-gaps investigation, COMPUTED) gives the water-rich share among rocky and water planets of 0.5 to 30 Earth masses with periods 0.5 to 100 d as a logistic in log10 M: centre 0.788 (6.1 Earth masses), slope 6.28 per dex. Bin fractions: 1 to 2 Earth masses 0.02, 2 to 3 0.09, 3 to 5 0.23, 5 to 8 0.49, 8 to 12 0.75, 12 to 20 0.91.

The probability of b4 is now `q = c_b4 L(M)` with `L` that logistic. Since the weights of b1 to b3 sum to 1, the weight `q/(1 - q)` makes the probability of b4 exactly `q`. The logistic reproduces the six bins within 0.03 (a root mean square of 0.020, the tolerance a test uses at 2 sigma). The plateau `c_b4` is the flagged knob: 1 is the logistic as fitted, and it is bracketed 0 to 1, uniform, because published water-rich shares run from at most 3 percent [R3801] to 59 percent [R0715] and 87 percent by formation model, and the synthesis's own plateau is 0.92. The widest range applies where two sources are in question.

## 2. What it is not

- A MODEL label: the numbers come from a physics simulation with disc priors, not from a measurement and not from a Solar System body. The row is partly verified, with that limit written in it.
- It is not tuned to the density modes of [R2774] (0.94, 0.47, 0.24 with their selection): that comparison needs radii from composition and the mass-radius selection model (GEN.505), and mass and radius alone cannot break the degeneracy with hydrogen-helium envelopes (about 6 informative JWST spectra). It stays open under GEN.510.
- Rows b1 to b3 (dry, low, high delivery by the snow-line ratio) stay labelled stand-ins: their payloads are Earth delivery scenarios and formation distance is not observable. b5 (the continuous water fraction) and b6 (rogue water, not a separate draw) also stay stand-ins; T9 and T11 have no ensemble beyond Earth's budget.
- The stellar-mass scaling of the minimum ice mass (3 (M*/Msun)^0.48 from two figure reads) is not built.
- The recipes are not read by the generator yet: no generated value changes, no epoch.
