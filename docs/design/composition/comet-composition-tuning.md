# The comet composition rows k1 to k5 and k7 (GEN.520)

Status: built 2026-10-24 by lane M (builder B12). Rows k1 to k5 and k7 of the atmosphere recipes (`scripts/atmosphere_recipes_data.py`); 23 targets in `physics/data/tuning_targets.json`; tests in `src/tests/test_tuned_draws.py` and `src/tests/test_tuning_harness.py`. The numbers are those of the small-bodies source check (`tunable-gaps-sources/group5-small-bodies-moons.md`), computed from the Robinson et al. 2024 database (CC BY 4.0, per-comet median of the detections at heliocentric distance up to 2.5 AU) and from Harrington Pinto et al. 2022 (HP22).

## 1. k1 to k3: the coma class

The weights 0.5, 1/3 and 0.1 were rounded words of an abstract and summed to 0.93. They are now the printed counts of HP22 Table 2: 13 CO2-dominated, 9 CO-dominated and 3 comparable of 25 comets. Three proportion targets check the shares; the harness gained a `category` field for a proportion of the alternatives a group of options draws, and `recipe_sampler` returns the chosen option for an option row.

Limits kept on the rows (the state is `partly verified`): the class depends on the heliocentric distance (every coma beyond about 3.5 AU is CO-dominated), and 9 of the 13 CO2-dominated entries rest on an upper limit on CO, so the CO2 share is inflated by an amount that cannot be corrected without the list of excluded comets. The source check's better design (draw log10 CO/CO2 at the comet's own distance and derive the class from the 2/3 and 4/3 thresholds) needs a distance-dependent model and is left open (section 4).

## 2. k4 and k5: CO2/H2O and CO/H2O

| Row | Before | Now | Target |
|---|---|---|---|
| k4 CO2/H2O | log-normal, median 12 percent, sigma 0.5 dex ("by eye") | median 10.9 percent, sigma 0.51 dex | mean -0.963 +- 0.11 dex, sd 0.51 +- 0.08 (21 comets) |
| k5 CO/H2O | median 3 percent, sigma 0.6 dex, truncated to 0.5 to 23 percent | median 2.5 percent, sigma 0.65 dex, no truncation | mean -1.602 +- 0.16 dex, sd 0.65 +- 0.12 (16 comets) |

The truncation of k5 is dropped: 0.3 to 26 percent is the extremes of 16 comets, not a bound, and a cut at the sample's extremes narrows the spread under the published one. The old and the new values agree within the standard errors, so the data do not move the rows much; the test that guards them fails a wrong centre (10 percent for CO) or a spread of 0.2 dex.

## 3. k7: single-comet species

CH4, C2H6, H2CO and H2S have populations of 9 to 18 comets (geometric means 0.68, 0.41, 0.25 and 0.71 percent; sigma 0.33, 0.33, 0.35 and 0.39 dex), replacing the values of 67P, which were 2 to 10 times off the population in either direction. Eight targets check them. N2, O2, OCS, SO2, Ar, Kr and Xe are measured on one comet (67P) and keep their values with the design width of 0.6 dex, named in the row's limit; O2 presence cannot be set from one comet.

## 4. What it does not do

- **The admission rule.** The statistics are of comets that were observed. The draws do not depend on the comet's distance, so the admission model of GEN.507 (which weights distance) does not change them and the targets use none. A coma class that depends on the distance does need it.
- **The system offset.** The source check proposes a per-system offset of the log ratios, N(0, sigma_sys) with an upper edge near 1.5 dex, from the two interstellar objects (N = 2). Two objects bracket it only; it is not drawn, and the rows say so. Filed to the TODO thread.
- **Correlation.** The one-factor model (mean pair correlation 0.13) is not drawn; k6 and k7 species are independent.
- The recipes are not yet read by the generator, so no generated value changes and no epoch is taken.
