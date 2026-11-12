# The comet species model (GEN.1005, k7)

Status: module, rows and flip, written 2026-10-30 (as the lane reported) by the Aux 3 lane. Section 4 (GEN.1021) is wired: row k7 draws the seven species log-uniform between the comets measured (the 67P calibration, rows CK1, CK2, CK5 and CK6, is retired). Section 5 (GEN.1022) is a module not yet wired: its swap rides a train with one epoch. `src/planetgen/physics/comet_species.py`; tests in `src/tests/test_comet_species.py`; registry rows CK1 to CK6 (section 6.59 of the gap-filling plan); the row is `atmosphere-k:k7`, no longer a stand-in since GEN.1021.

## 1. What replaces what

Row k7 drew N2, O2, OCS, SO2, Ar, Kr and Xe of a comet as log-normals about the values of 67P [R3835], with a design width of 0.6 dex. Those seven species have a population of one comet, so no median can be taken from a population; the other four species of the row (CH4, C2H6, H2CO, H2S, 9 to 18 comets each) are unchanged.

For N2, Ar, OCS and SO2 the median is now a chain: the molecules of the species per CO in the nebula of the comet's host (its own N/C, Ar/C and S/C over the proto-solar ones, new subject quantities `n_over_c`, `ar_over_c`, `s_over_c`, 1 for a host of solar composition), times the share of that which is in the ice, times the population's median CO/H2O (10^-1.602, row k5). The carbon in CO is the 1.5 per 10^4 hydrogen of [R3322] over the proto-solar carbon (0.52). O2, Kr and Xe have no host ratio here and keep the value of 67P. The width of every species is 0.4 dex.

| Species | Nebular, per CO (solar host) | Share in the ice (CK1, CK2) | Median per water | 67P (calibration) |
|---|---|---|---|---|
| N2 | 0.225 | 0.158 | 8.9e-4 | 8.9e-4 |
| Ar | 0.016 | 0.0145 | 5.8e-6 | 5.8e-6 |
| OCS | 0.088 | 0.187 | 4.1e-4 | 4.1e-4 |
| SO2 | 0.088 | 0.592 | 1.3e-3 | 1.3e-3 |

## 2. What this is and is not

It makes the seven species depend on the star: a host with twice the nitrogen per carbon gives twice the N2. It does not remove 67P as the anchor: the four shares and the three unscaled values are tuned to the one comet measured, so the test that the chain gives 67P at a solar host is the calibration, not a check. They are registered stand-ins (tier T7, tweak L3), with GEN.1021 to replace them with a second comet or a laboratory table of trapping by formation temperature. Bar-Nun (NASA N90-16610: the trapped gas drops six decades between 25 and 100 K; CO at 3.5 percent needs 48 +/- 5 K) says the share depends on a temperature the generator does not carry per comet; the NTRS text had no per-species tables.

The width 0.4 dex is Bar-Nun's 5 K at a decade per 12 K (CK6, T6, L1), against the old design 0.6; the four species a population measures scatter 0.33 to 0.39 dex.

## 3. Limits

- The sulfur of OCS and SO2 (0.78 of the CO-scaled nebular sulfur) plus the population's H2S exceeds what the chain allows: the sulfur is not conserved across rows.
- Krypton and xenon are not in the abundance table of the generator (`tuning.STELLAR_SOLAR_LOG_EPS` stops at the metals the star draw carries), so they keep one value.
- The star's ratios are quantities the subject must carry (`required_inputs('comet')` now lists them); `mixture_subject.subject_of_body` builds planets and moons, and a comet's subject is built by its caller until the generator draws comets (the recipes are not yet read by it). A caller with no star composition passes 1.

## 4. The brackets between the comets measured (GEN.1021, `comet_species.bracket_over_h2o`)

The 67P calibration is replaced by the range between the comets measured (the widest range where two sources are in question). Report: `reports/gen1021-comet-species-second-comet-2026-10-30.md`.

| Species | Low end | High end | Width | Sources |
|---|---|---|---|---|
| N2 (per CO, times the median CO/H2O of 2.5 percent) | 0.0287 (67P) | 0.09 (C/2016 R2) | 0.5 dex | R0109, R0836 |
| Ar (per water) | 5.8e-6 (67P) | 8.7e-3 (Hale-Bopp) | 3.2 dex | R3835, R4215 |
| OCS (per water) | 5e-4 | 4e-3, never above the nebular sulfur | 0.64 dex | range of the comets, 46P paper Table 7 |
| SO2 (per water) | 3e-4 | 2.3e-3, never above the nebular sulfur | 0.86 dex | the same |
| O2 (per water) | 0.020 | 0.0465 | 0.37 dex | Halley R3836 to 67P R0386, by their errors |
| Kr, Xe | argon times 0.089, krypton times 0.269 | the same | 3.2 dex | 67P R3834 (one comet) |

The ends scale with the host: N2 with `n_over_c`, argon, krypton and xenon with `ar_over_c`, OCS and SO2 with `s_over_c`. Argon is moved off the CO chain onto the water: Hale-Bopp's argon per CO (0.043) is 2.7 times what the nebula holds per CO, so no trapped share of the CO chain can reach it. Rows to retire when this is wired: CK1, CK2, CK5, CK6 (CK3 stays for the nebular amounts, CK4 is row k5). New rows: CB1 (the brackets, T3, L0) and CB2 (the 67P noble-gas ratios, T2, L1).

## 5. One latent ice formation temperature (GEN.1022, `comet_ice_formation`)

The data say the species move together: the two CO-rich comets are argon-rich and nitrogen-rich per CO and 67P is poor in all of them, and Bar-Nun's trapped gas falls six decades from 25 to 100 K. Independent brackets (section 4) would allow a CO-poor, argon-rich comet, which the data do not show. So a comet has one latent, the temperature its ice formed at, and it is read off CO/H2O, the quantity measured for every comet, so row k5's spread (2.5 percent, 0.65 dex) is unchanged:

1. `T = 48 K + 12.5 K x log10(0.035 / (CO/H2O))`: the laboratory law through Halley's anchor (3.5 percent CO at 48 +/- 5 K) at 12.5 K per decade (CY3). The median comet is 49.8 K, Hale-Bopp (CO/H2O 0.2, CY1) 38.5 K.
2. The position `xi` is 0 at the median comet's temperature and 1 at Hale-Bopp's, held at the ends.
3. N2, Ar, Kr and Xe are interpolated in the logarithm between the two ends of their section 4 brackets at `xi`, plus a residual of their own (CY2: N2 0.05 dex, Ar 0.4 dex from 67P's own season, Kr and Xe 0.15 dex from the lab's 40 percent). N2 is per CO, so it also follows the comet's own CO/H2O. Krypton is the comet's argon times 0.089 and xenon its krypton times 0.269.
4. OCS, SO2 and O2 are not driven by the temperature and keep their section 4 brackets (log-uniform).

Tests (never inputs): the median comet gives the warm end and Hale-Bopp the cold end of every species; the species correlate with CO across a drawn population (above 0.7) and a CO-poor, argon-rich comet is under 1 percent; argon across comets spans 67P to Hale-Bopp and N2 per CO spans 67P to R2.

Limits. The temperature is a latent: Ciesla and others 2018 show that the trapping also depends on the deposition rate and the partial pressure. Hale-Bopp's CO/H2O is recalled, not opened (CY1). Comets warmer than the median all take the warm-end values (plus residuals). The sulfur of OCS and SO2 plus the population's H2S is still not conserved across rows.

Wiring (a later part, rides a train with one epoch): rows k5 and k7 draw the CO/H2O first and then each latent species from `latent_species_over_h2o` with one standard normal each (a new draw kind in `atmosphere_recipes`), and OCS, SO2 and O2 from `bracket_over_h2o` as `loguniform` rows; `TRAPPED_OVER_NEBULAR`, `UNSCALED_OVER_H2O`, `median_over_h2o`, `sigma_dex` and rows CK1, CK2, CK5 and CK6 go, with their burn-down lines. The stored comet columns move for every comet, so the epoch is named in the note.
