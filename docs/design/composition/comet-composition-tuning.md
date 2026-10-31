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
- **The system offset.** The source check proposes a per-system offset of the log ratios, N(0, sigma_sys) with an upper edge near 1.5 dex, from the two interstellar objects (N = 2). The live comet sampler (`physics/comet_population.py`, GEN.495) already draws an interstellar comet's CO and CO2 from the wider ranges one measured object forces, and its Solar System species agree with these rows (CH4, C2H6, H2CO, H2S within 0.05 dex, CO2 at 11 against 10.9 percent; its CO2 width 0.27 dex is at r_h up to 2 AU, the rows' 0.51 dex at 2.5 AU, the widest range rule). A shared system-level offset needs a model of the system; it is not drawn here and the rows say so.
- **Correlation.** The one-factor model (mean pair correlation 0.13) is not drawn; k6 and k7 species are independent.
- The recipes are not yet read by the generator, so no generated value changes and no epoch is taken.

## 5. The ice zones of the disc (GEN.595, lock-free step)

The population shape of rows k1 to k5 is the Solar System's (decision D6: the Solar System population as the within-system shape plus a per-system offset knob bracketed by the two interstellar objects). `physics/comet_ice_budget.py` (registry CI1 to CI5, section 6.46) derives what the disc allows:

- **Lines.** The average disc of Andrews and Williams 2007 (T = 200 K at 1 AU, index 0.62, scaled by the star's luminosity to the quarter) with the evaporation temperatures of Oberg, Murray-Clay and Bergin 2011 (H2O 135 K, CO2 47 K, CO 20 K) puts the lines of a Sun-like star at 1.9, 10 and 41 AU.
- **Budget.** Table 1 of that paper (O and C per 10^4 hydrogen, from one disc's ices) gives CO2/H2O 0.33 and CO/H2O 1.67 for a body that froze everything, scaled by the star's carbon (both) and its oxygen after the rock's share (the water). A body between the lines freezes what is below its line and traps 5 to 10 percent of the rest in water ice (Fayolle et al. 2011, as Oberg et al. cite it).
- **Where bodies formed.** Equal masses in equal intervals of ln r (a surface density as 1 / r): for a formation range of 5 to 50 AU round the Sun, 32 percent of the bodies are in the water zone, 60 percent in the CO2 zone and 9 percent in the CO zone.

**What it reproduces and what it does not.** 2I/Borisov's CO/H2O of 130 to 155 percent (Bodewits et al. 2020) is the CO zone's 167 percent to within 30 percent: an interstellar comet with a CO-rich coma is a body from beyond its star's CO line, with no offset knob. The Solar System's median comet (CO2/H2O 10.9, CO/H2O 2.5 percent) is not in any zone: it has more CO2 than CO as the CO2 zone has, at 0.5 dex less CO2 and 0.7 dex less CO, so the absolute level of the Solar System population stays a stand-in (a loss of ice by heating or the disc's own loss of CO, both outside this model). 3I/ATLAS's CO2/H2O of 7.6 (Cordiner et al. 2025) is over twenty times the most this budget gives for any carbon and oxygen: it needs a water-poor body or processing of CO to CO2 (cosmic-ray processing is proposed in the literature), which is not in the model, so `comet_population`'s CO2/H2O bracket of 10 stays a registered stand-in.

Left out: transport of icy grains, the walk of the lines with the star's age, the CO2 trapping fraction (the CO value is used), the refractory carbon.

**Wired (GEN.595, train 7).** A comet of a star's system draws where it formed from a stream of its own (keyed by the comet's orbit, so the run's stream is untouched), log-uniformly between the star's snow line and 18 times the snow line (the edge `DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER` the planets are placed to; CI6), and `comet_ice_budget.formation_ratios` gives the CO2/H2O and CO/H2O of the ice there from the star's luminosity and its carbon and oxygen over magnesium against the Sun's (`bulk_draw.ice_forming_ratios`). `comet_population.draw_volatiles(ice_ratios=...)` takes the median total and the CO share from them in place of the Solar System's median and class share, keeps the population's width (0.2 dex in the total, 1.5 in the logit of the share), and moves the observed bounds with the median (they keep their distance from the population's geometric mean, so a zone past the Solar System's 30 percent CO is drawn). The period class no longer sets the CO share; it keeps the carbon-chain flag. A star without abundances (a compact remnant) and an interstellar comet keep the population's medians and the bracket.

What it changes for a Sun-like star: a comet's median CO over water is 12.5 percent (water and CO2 zones) or 167 percent (the CO zone), and CO2 over water 2.5 or 33 percent, where the Solar System's measured medians were 2.5 and 10.9 percent. The model gives the Solar System population neither level (above): the absolute level is the model's, as the physics rule has it, and the measured medians stay tests.
