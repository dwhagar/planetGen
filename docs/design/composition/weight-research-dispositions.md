# Weight research dispositions: what GEN.480 to GEN.494 and GEN.501 found for the 43 open atmosphere rows

Status: closes the research items GEN.480 to GEN.494 and GEN.501 (lane M, 2026-10-23). Source of the findings: `docs/design/composition/tunable-gaps.md` section 3 and the worker reports in `docs/design/composition/tunable-gaps-sources/`, `docs/design/boss-papers-2026-10-17.md` and `docs/design/open-gaps-closest-laws-2026-10-17.md`; nothing in this note is a new reading of a paper. What this note adds is the decision for each row, written into the machine columns (`scripts/atmosphere_recipes_data.py`) and the data file the loader reads.

## 1. The rule applied

Before this change 43 of the 104 rows of part 09 were `open`: a `DESIGN` weight the generator would draw inside a bracket. The research items were to turn each into `SOURCED` or `DERIVED`. The tunable-gaps investigation had already read what the public record holds for each row. Boss's rule (1:23 PM EDT, 2026-10-17) settles what to do where it holds nothing: **a row with no measurement ships as a labelled stand-in at the widest bracket, replaced when a measurement exists**; and the verified-wins, newest-wins, widest-range-wins rules decide between sources. So every `open` row now has one of two states:

- **partly verified** with a stated limit, where a published statistic, gate or bracket exists and only an edge, a width or a share is a design choice. The limit says which.
- **stand-in** naming the builder item (GEN.509 to GEN.525, GEN.595) that replaces it, where no measurement exists, the value rests on Solar System bodies, or the only number is a SNIPPET or RECALLED value (the target table never tunes to those).

No row is `verified` or `derived` by this change: no source was newly opened. The loader therefore lists no `open` row, and 30 stand-ins (2 before: 28 new) and 38 partly verified rows (23 before: 15 new) carry the open work, each pointing at the builder that owns it. (Later builders move rows: GEN.516 retuned f5 to a verified source, and GEN.514 turned o4 from a stand-in into a partly verified formula; the counts are then 37 partly verified and 29 stand-ins.)

The recipes are not yet read by the generator, so this change draws nothing differently and takes no generator epoch.

## 2. Row by row

| Research item | Rows | New state | Why | Replaced by |
|---|---|---|---|---|
| GEN.480 (AR1) envelope, water-rich share, envelope Z | b4 | stand-in | published water-rich shares run from at most 3 percent [R3801] to 87 percent by formation model; decision D1: keep the NGPPS logistic as the labelled MODEL default, bracket 0 to 0.87 | GEN.510 |
| | e3, e4 | partly verified | the hycean gate (Teq below about 400 K, H2 at most 0.1 wt percent, water at least 10 wt percent) is [R2804]; the share is p[b4], a stand-in | GEN.510 |
| | e5 | partly verified | bracket edges from about 8 JWST retrievals and one 15-planet fit; no population mean exists | GEN.509 |
| GEN.481 (AR2) delivery | b1 to b3 | stand-in | the Ortenzi payloads are Earth delivery scenarios; formation distance is not observable; the only function of water against distance is population synthesis (NGPPS, logistic centre 6.1 Earth masses, slope 6.3 per dex, COMPUTED) | GEN.510 |
| | b5, b6 | stand-in | no ensemble measures sub-percent water (Boss, 1:23 PM EDT, 2026-10-17); rogue water is not a separate draw | GEN.510 |
| GEN.482 (AR3) redox | a1 | stand-in | the Earth end and the mass slope rest on Solar System bodies; the low-mass end is tunable to the white dwarfs once the scale conversion (core formation against surface magma ocean) is tested | GEN.513 |
| | a2 (already a stand-in) | stand-in | law fitted to the Moon, Mars and Earth | GEN.513 |
| | a3, a4 | partly verified | 22 polluted white dwarf values [R1152, R1153]: centre -1.25, sd 0.45; the sample carries detection and settling selection; the shipped draw is the first tunable row (see section 4) | GEN.513 |
| GEN.483 (AR4) shoreline | d6 | stand-in | only upper limits and one detection (55 Cnc e); the cap is derived as min(supply, vapour, escape) | GEN.511 |
| GEN.484 (AR5) giant C/O | f7 | partly verified | stellar C/O measured (median 0.47, sd 0.07, N=849); planets 6 of 69 above 0.8 (LExACoM) with a 0.2 to 0.3 method systematic | GEN.516 |
| GEN.485 (AR6) hot Jupiters | h2, h3 | stand-in | the regime edges are pressure-dependent outputs of equilibrium chemistry, not measurements | GEN.518 |
| | h4 | partly verified | an empirical Ti threshold between 2350 and 2640 K from N=3 to 6 planets; no population measures Kzz | GEN.519 |
| GEN.486 (AR7) lava worlds | i6 | partly verified | 2 of 2 molten-dayside ultra-short-period planets keep volatiles (55 Cnc e, TOI-561 b); N=2, Beta(3,1) interval 0.29 to 0.99 | GEN.511 |
| GEN.487 (AR8) evolved and frost | j1 | partly verified | the runaway gate is [R2446]; the retained fraction is a design residual in [0, 1] | GEN.511 |
| | j5 | stand-in | the oxygen-sink fraction is unsourced (Luger and Barnes against [R1813] span from 0 to 2400 bar per ocean) | GEN.511 |
| GEN.488 (AR9) comets | k1 to k3 | stand-in | counts 13/9/3 of 25 are a SNIPPET (Harrington [R3535: Table 4] not read); the target table never tunes to a SNIPPET | GEN.520 |
| | k4, k5 | stand-in | the Solar System population is not a prior for other systems; 2I and 3I lie 1.8 to 3.6 sigma outside the fits | GEN.520 |
| | k7 (already a stand-in) | stand-in | single-comet species; O2 presence cannot be set from one comet | GEN.520 |
| GEN.489 (AR10) white dwarfs and pulsars | l1, l2, l3 | stand-in | the exosphere threshold (log10 P = -6 +/- 0.5) is a design constant with no data; the surviving-giant share has upper limits and two detections only | GEN.513 |
| GEN.490 (AR11) trace gases | m2, m5, m6 | stand-in | abiotic O3 is a photochemistry look-up; HCl and HF are speciation of a crust halogen budget that is Earth-only | GEN.524 |
| | m3, m4, m11 | stand-in | no measurement of N2O, PH3, NO2 or CH3Cl on any exoplanet; tuning to Earth or Venus breaks NO FIXED CONTENT | GEN.525 |
| GEN.491 (AR12) clouds and haze | n6 | partly verified | A_H against Teq is measured (Fu 2017, N=34, slope 0.0008 +/- 0.00016 per K); the 950 K crossover is a model output | GEN.519 |
| | n7 | stand-in | lab rates span 250 times and every population statistic is degenerate with metallicity, Kzz and cloud cover | GEN.519 |
| GEN.492 (AR13) radioisotopes | o2 | partly verified | [Th/Mg] for thin disc (N=215) and thick disc (N=32) from Hypatia; halo, young and M hosts have no data | GEN.514 |
| | o3, o4 | stand-in | planetary K is volatile-depleted and only Solar System bodies give the depletion; Th/U calibration is Solar System only | GEN.514 |
| | o6 | partly verified | the enriched share spans 0.04 to 0.66 (clusters, white dwarf pollution); no single distribution | GEN.514 |
| GEN.493 (AR14) rotation and saturation | p6 | partly verified | percentile tracks of Tu 2015 (READ); the 0.25/0.5/0.25 weights give way to a uniform percentile | GEN.523 |
| GEN.494 (AR15) brown dwarfs | g1, g4 | partly verified | type is a function of Teff with published scatter ([R2374], N=525); the age edge and the 500 to 700 K edge are design choices | GEN.517 |

Stand-ins that already named a research item now name their builder: a2 (GEN.482 to GEN.513) and k7 (GEN.488 to GEN.520).

## 3. GEN.501: the NH3 activity factor

The photolysis lifetime of NH3 (`physics/data/photolysis_lifetimes.json`) uses an activity factor of 1.075, the middle of 1.05 to 1.10, settled with the body lane because NH3 absorbs mainly at 170 to 217 nm. The [R2040] rate tables (Planet. Space Sci. 106, 11-45) are the source that would replace it. The paper was not reachable from this lane, and the `sbpy` package that codes part of the tables (read at version 0.6.0) holds H2O, OH, HCN, CH3OH, H2CO, CO, CO2 and CN but no NH3. So the factor stays the flagged bracket 1.05 to 1.10 with its reasoning, no number changes, and GEN.501 stays open as a data request (the rate table for NH3 from the paper).

## 4. What this leaves for the builders

- Rows marked `tunable` need a target in `physics/data/tuning_targets.json` (GEN.503). Only a3 and a4 are tunable today. Each builder marks its rows tunable when it adds the targets.
- The selection models (GEN.504 to GEN.508) are the dependency for the tunable rows whose published statistic comes from a biased sample (a3, a4 need GEN.508; the radius valley needs GEN.504).
- a1: the Sossi floor (-2) and the white-dwarf centre (-1.25) are on different redox scales (core formation against surface magma ocean). The mass-law floor is therefore left as shipped, and GEN.513 tests the conversion before moving it.
- A stand-in whose value is a bracket is drawn at the widest bracket with the stand-in label (registry in `docs/design/gap-filling-plan.md` section 6); nothing here is a default until a builder registers it.
