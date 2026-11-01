# The comet species model (GEN.1005, k7)

Status: module, rows and flip, written 2026-10-30 (as the lane reported) by the Aux 3 lane. `src/planetgen/physics/comet_species.py`; tests in `src/tests/test_comet_species.py`; registry rows CK1 to CK6 (section 6.59 of the gap-filling plan); the row is `atmosphere-k:k7`, still a flagged stand-in, replaced by GEN.1021.

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
