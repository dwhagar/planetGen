# Volatile delivery against formation distance (GEN.481, AR2)

Status: research note, written on the evening of 2026-10-24 (as the lane reported) by the Citations lane. Rows b1 to b4 of `scripts/atmosphere_recipes_data.py` were rewritten from it on 2026-10-25 (GEN.802, `water-delivery-tuning.md` section 0: one joint law in x and mass, refitted with a hinge at x = 0.6, replaces the separate fits of section 3; rows b5 and b6 are untouched). Section 7 is the reasoning that led there. The numbers are COMPUTED here from a published population-synthesis dataset; they are a model result, not a measurement.

## 1. The answer in short

- In the Bern NGPPS population (NG76, 1 Msun stars), the share of planets that hold at least 10 percent water by mass (the water-rich row b4) rises from 2 percent at a starting distance of 0.15 of the snow line to 26 percent at 0.55, 51 percent at 0.75, 85 percent at 0.95 and 99 percent at 1.75. As a logistic in ln x (x = a_form / a_snow, a_snow = 2.69 AU for 1 Lsun, part 09 row b1): centre x = 0.66 (1.78 AU at 1 Lsun), width 0.25 in ln x, all planet masses together (the same centre and width come out of the 30,581 planets of 10 Earth masses or less and of the 10,506 planets of 0.5 to 30 Earth masses).
- The mass dependence at fixed formation distance is weak (0.23 in logit per dex of mass for the water-rich share). The strong mass law of GEN.510 (a 1 Earth-mass planet 2 percent water-rich, a 12 Earth-mass planet 91 percent) is a statement about planets that sit close to the star at 5 Gyr: massive planets migrate in from beyond the snow line. The two laws are two views of one process and must not be multiplied (section 7).
- The design weights of rows b1 to b3 (centre x = 1, width 0.4, Earth delivery scenarios) sit too far out: at x = 0.75 the design gives 0.5 percent water-rich and 45 percent dry; the population gives 60 percent water-rich and 7 percent dry. The current gate of b4 (x at least 1) gives zero water-rich weight where the population gives 26 to 60 percent (x = 0.55 to 0.75), although only 7 percent of all water-rich planets start inside x = 1 (2.5 percent inside 0.75), because most embryos start far out.
- The "Low" class of Ortenzi (0.05 wt percent water) has almost no support: at most 9 percent of the planets (usually 4 to 6) land between 0.025 and 0.16 percent water; planets are either (almost) dry or hold at least 0.16 percent.
- The snow-line law over time is not recoverable from this dataset (section 5). The effective centre 1.78 AU against the instantaneous water ice line of the disc model (about 3 AU, 3.4 AU and 4.5 AU in the systems the model paper shows) says the ice line sweeps inward during planet growth, which carries ice to embryos that start inside it.

## 2. Sources, and what was read

| Source | Read | Used for |
|---|---|---|
| Burn, Mordasini, Mishra, Haldemann, Venturini, Nature Astronomy 8, 463, "A radius valley between migrated steam worlds and evaporated rocky cores" [R0600, R0601] | the arXiv text, sections on the formation model and the water treatment | what the model is (Bern model, NGPPS, condensed-water run) |
| Zenodo record 10527170, `code_upload_v2.zip` (137 MB): `condensed/ref_evo5e9.dat` and the paper's own plotting code `plot_biased_compo.py` | opened and analysed (the column names come from `rename()` in the plot code) | per-planet mass and water fraction at 5 Gyr, 31,909 planets in 1000 systems |
| DACE population ng76 (dace_query, public mode), snapshot at 100,000 yr | queried | the starting semi-major axis of each (system, planet) |
| [R1268] | text on the disc, the condensation model and the ice lines | the water ice line at about 3 AU (3.4 AU and 4.5 AU in two example systems) from an equilibrium condensation model [R4336]; planetesimal redistribution is not modelled |

Not used and why: [R4335] is about refractories (the tunable-gaps note already says so); [R4336] is cited by Emsenhuber for the condensation model and was not opened.

Method. The starting distance a_start is the semi-major axis at 100,000 yr. A 200-planet check against the first row of the planet's own track gives a ratio of median 1.000 (interquartile 0.999 to 1.000; the 5th and 95th percentile within 0.5 percent), so a(100 kyr) is a_start. The water fraction is column 51 of the condensed-water run, which the paper's plot code names `fice` and uses to flag water-rich planets (`fice > 0` with no envelope). The paper text I read does not define the column; its plateau (0.505 to 0.55 beyond 3 AU) is what a 50 percent ice share of solids beyond the ice line gives. That share is the model's input and is not the project's BC1 value (0.14 to 0.39 from the oxygen budget): use the data for the odds of each delivery class, and BC1 for the amount of water in solids.

The script that reproduces every number below is `reports/gen481-files/analyze_ngpps_formation_water.py` in the project files (it downloads the two inputs; the raw files never enter the repository).

## 3. Results

Classes by water mass fraction f (the mapping to the design classes is section 6): dry f < 0.00025, low 0.00025 to 0.0016, high 0.0016 to 0.1, water-rich f at least 0.1. x = a_start / 2.69 AU. Planets of 0.5 to 30 Earth masses at 5 Gyr, any final orbit.

| x | n | dry | low | high | water-rich | mean f |
|---|---|---|---|---|---|---|
| 0.0 to 0.1 | 352 | 0.963 | 0.011 | 0.023 | 0.003 | 0.001 |
| 0.1 to 0.2 | 486 | 0.794 | 0.058 | 0.126 | 0.023 | 0.007 |
| 0.2 to 0.3 | 392 | 0.628 | 0.087 | 0.237 | 0.048 | 0.016 |
| 0.3 to 0.4 | 293 | 0.529 | 0.055 | 0.307 | 0.109 | 0.032 |
| 0.4 to 0.5 | 243 | 0.469 | 0.049 | 0.309 | 0.173 | 0.046 |
| 0.5 to 0.6 | 179 | 0.352 | 0.039 | 0.352 | 0.257 | 0.060 |
| 0.6 to 0.7 | 126 | 0.262 | 0.079 | 0.397 | 0.262 | 0.077 |
| 0.7 to 0.8 | 136 | 0.132 | 0.059 | 0.301 | 0.507 | 0.149 |
| 0.8 to 0.9 | 160 | 0.038 | 0.006 | 0.206 | 0.750 | 0.241 |
| 0.9 to 1.0 | 235 | 0.000 | 0.004 | 0.145 | 0.851 | 0.289 |
| 1.0 to 1.2 | 589 | 0.000 | 0.000 | 0.090 | 0.910 | 0.292 |
| 1.2 to 1.5 | 831 | 0.002 | 0.000 | 0.045 | 0.953 | 0.315 |
| 1.5 to 2.0 | 1015 | 0.000 | 0.000 | 0.007 | 0.993 | 0.419 |
| 2.0 to 3.0 | 1720 | 0.000 | 0.000 | 0.002 | 0.998 | 0.480 |
| 3.0 to 15 | 3749 | 0.000 | 0.000 | 0.001 | 0.999 | 0.507 |

Split by planet mass (the full tables are in the script's output): for 0.5 to 2 Earth masses the dry share at x = 0.2 to 0.3 is 0.50 and the water-rich share 0.05; for 2 to 30 Earth masses they are 0.71 and 0.05. The water-rich share at x = 0.7 to 0.8 is 0.58 and 0.44. The sub-rich classes depend on mass (small planets are less often dry at a given x), the water-rich share almost not.

Logistic fits (maximum likelihood, 10,506 planets of 0.5 to 30 Earth masses; M in Earth masses):

| Event | logit P = a + b ln x + c log10 M | Centre at M = 1 | Brier score (fit against a constant) |
|---|---|---|---|
| f at least 0.1 (water-rich) | 1.649 + 3.997 ln x + 0.231 log10 M | x = 0.66 (1.78 AU at 1 Lsun) | 0.034 against 0.162 |
| f at least 0.0016 (high or water-rich) | 3.500 + 3.048 ln x - 0.894 log10 M | x = 0.32 | 0.038 against 0.121 |
| f at least 0.00025 (not dry) | 3.629 + 2.828 ln x - 0.888 log10 M | x = 0.28 (0.75 AU at 1 Lsun) | 0.040 against 0.113 |

The three fits cross each other at extreme x and mass (861 of 2400 grid points violate the ordering), so they are not a joint model. Use the binned table for the class probabilities, or one fit with one threshold at a time; the water-rich logistic is the headline and the only one that is safe to use alone.

## 4. Against the design weights of part 09 (b1 to b4)

Design: s = 1/(1 + exp(-ln x / 0.4)); dry (1-s)^2, low 2s(1-s), high s^2, water-rich c s^3 with c between 0.05 and 0.4 and a gate at x of at least 1. NGPPS: the table above (0.5 to 30 Earth masses).

| x | design dry | low | high | water-rich | NGPPS dry | low | high | water-rich |
|---|---|---|---|---|---|---|---|---|
| 0.3 | 0.91 | 0.09 | 0.00 | 0.00 | 0.54 | 0.07 | 0.32 | 0.07 |
| 0.5 | 0.72 | 0.26 | 0.02 | 0.00 | 0.37 | 0.04 | 0.36 | 0.23 |
| 0.75 | 0.45 | 0.44 | 0.11 | 0.01 | 0.07 | 0.02 | 0.30 | 0.60 |
| 1.0 | 0.25 | 0.50 | 0.25 | 0.02 | 0.00 | 0.00 | 0.15 | 0.85 |
| 1.5 | 0.07 | 0.39 | 0.54 | 0.06 | 0.00 | 0.00 | 0.02 | 0.98 |

(The NGPPS columns are windows of plus and minus 15 percent around each x, for planets of up to 10 Earth masses, 313 to 1347 planets each.) The design is therefore both shifted outward by a factor of about 1.5 and too gradual, and it puts the water in the "low" and "high" classes where the population puts it in the water-rich class.

## 5. The snow line over time

The dataset has no time column for the ice line. Two facts bound it. First, the disc model places the water ice line at about 3 AU at the moment of the figures in Emsenhuber et al. (3.4 AU and 4.5 AU in the two example systems; the jump at the water ice line in the solid surface density is "the most consequential one at about 3 au"). The design a_snow of 2.69 AU for the Sun is consistent with this. Second, the effective centre of the water-rich share is 1.78 AU (x = 0.66) and the dry share vanishes at x of about 0.9 (2.4 AU): inside the instantaneous ice line the planets still hold water. The model has planetesimals accreted in the feeding zone, embryo-embryo collisions, and a cooling disc whose ice line moves inward during the few Myr of growth; Emsenhuber et al. state that planetesimal redistribution (drift) is not included, so the water inside x = 1 comes from the moving ice line, from collisions with material that formed beyond it and from migration, not from drift. How much each contributes is not separable with these columns. Decision for GEN.458 (the snow-line law): keep T_snow = 170 K (bracket 161 to 213 K) for the instantaneous line and carry the inward sweep as the offset between the logistic centre (0.66) and 1; do not build a time-dependent ice line from this note. A time law needs the disc model's own ice-line track (not public in the files I reached) or an analytic viscous-disc line [R3193], neither opened here.

## 6. Mapping to Ortenzi's Dry, Low and High

[R3378] (Table 1, read) gives mantle water contents of 0, 0.05 and 0.5 weight percent for the three delivery scenarios. The population gives the bulk water mass fraction of the planet, which is not the mantle's water after core formation and outgassing. I map by nearest class in log water mass fraction: dry below the midpoint of 0 and 0.0005 (taken as 0.00025), low from 0.00025 to 0.0016 (the geometric midpoint of 0.0005 and 0.005), high from 0.0016 to 0.1, and everything of at least 0.1 as the water-rich row b4 ([R2804] hycean limit, in the part 09 text). Moving both inner thresholds by a factor of 3 up or down changes, at x = 0.3, 0.5 and 0.75, the dry share by at most 0.04, the low share by at most 0.05 and the high share by at most 0.09; the water-rich share does not change. The label is DERIVED from a model, with the mapping as an assumption: no source gives the mantle water of a planet formed at a given distance.

## 7. What this means for the rows, and why they are not rebuilt here

1. Rows b1 to b3 should take their weights from section 3 (binned in x, or a logistic in ln x for "not dry") and are then DERIVED from the model with the assumption of section 6; their payload (the Earth delivery scenarios) stays labelled.
2. Row b4 needs both laws reconciled. P(rich | x) of this note is a law of the formation distance; the GEN.510 logistic in mass is a law of close-in planets at 5 Gyr, where mass stands in for the migration distance. The generator, which does not migrate planets (planet mass is the isolation mass, GEN.650), should use the formation-distance law for a planet that formed at x and treat migration as the separate event that moves a planet from x to its final orbit. Using both multiplies one effect twice. The b4 gate at x of at least 1 is contradicted by the data (60 percent of the planets at x = 0.75 are water-rich).
3. Rows b5 (continuous water fraction) and b6 (rogue water) are not served by this dataset: b5 needs the mantle water after core formation, b6 needs an ejection model.
4. The tests of `test_atmosphere_recipes.py` pin the design numbers (for example b1 = 0.7782 at x = 0.447; b4 = 0.3638 at 5 Earth masses), and `test_water_delivery.py` pins the GEN.510 logistic. Wiring the new law changes both, so it is a code item with its own gate, not part of a research note. Proposed follow-up item: GEN.481b "Wire the formation-distance water law into rows b1 to b4 and reconcile it with GEN.510", prerequisite GEN.481, owner Lane Draws, 8.1, one registry row change per row (tier T6, tweak L0).

## 8. Limits

- One model, one star (1 Msun; the DACE populations for 0.1, 0.3, 0.5 and 0.7 Msun exist but their composition columns are not in the Zenodo file); scaling a_snow with the square root of the luminosity is the design's assumption, not tested here.
- The water fraction is the condensed-water run; in the nominal run the water is mixed into the envelope and `fice` is zero, so the nominal file cannot give a water-by-distance table.
- The population is the Bern model's disc-mass and metallicity distribution; the width 0.25 in ln x folds in the spread of the ice line across discs. A generator that draws a disc-dependent ice line would double count it.
- All numbers are a model output (source tier T6 for the centre and T7 for the width and the class split, tweak L0), no part is tuned.

## 9. Licence and data

Zenodo 10527170 carries no licence field; the DACE population ng76 is public. Nothing from either enters the repository: the tables above are derived counts. The terms are recorded in DATA-LICENSE.md.
