# The rock fraction of an icy body (GEN.522, row T19)

Status: built 2026-11-24 (builder B14 of [tunable-gaps.md](tunable-gaps.md), row T19; the moon rows T17 and T18 went to the design-flag registry, GEN.525). Code `src/planetgen/physics/small_body_composition.py` (`ROCK_FRACTION_LAWS`, `ICY_TABLE`), targets `t19-dwarf-rock-mean`, `t19-dwarf-rock-sd` and `t19-kbo-density-share` in `physics/data/tuning_targets.json`, test `src/tests/test_kbo_rock_fraction.py`, stand-ins CP5 and CP6 of [stand-in-values.md](../research/stand-in-values.md). The data are read in [tunable-gaps-sources/group5-small-bodies-moons.md](tunable-gaps-sources/group5-small-bodies-moons.md) section 2.4 and now in full from the author copy of Bierson and Nimmo ([R0390], Table A.2 of 17 diameters and densities, the porosity model and the fit).

## 1. What was asked and what the data allow

The row was a seeded bracket (0.58 to 0.82) that the draw already replaced by a truncated normal (`icy-intact`: N(0.70, 0.08) on [0.45, 0.95]; `icy-stripped`: N(0.88, 0.05) on [0.75, 0.98]), with open weights. The item asks for the knobs to be tuned to published statistics. Every statistic there is of one system's bodies: 17 binary Kuiper belt objects, and seven dwarf planets. The rule of tunable-gaps.md section 5 (NO FIXED CONTENT: the densities of the Kuiper belt objects and the dwarf planets are a Solar System-only constraint) applies, so none of the three statistics is a target to tune to. Each is a **test point** (`role` test-point with the reason in `solar_system`, `Target.tunes` false): the draw is simulated, compared within 2 sigma and reported, and a failure is a finding and never a reason to move a knob to the value.

## 2. What the comparison shows (one draw per seed, `simulate_and_compare`)

| Test point | Published | Simulated | z |
|---|---|---|---|
| `t19-dwarf-rock-mean`: rock fraction of the seven dwarf planets (a dwarf planet takes either row at a half) | 0.78, standard error 0.045 (0.82 SD 0.12 at 2700 kg/m3 rock, 0.74 SD 0.10 at 3500) | 0.790 | 0.2 |
| `t19-dwarf-rock-sd`: its spread | 0.12, standard error 0.035 | 0.109 | -0.3 |
| `t19-kbo-density-share`: the covered objects within 2 sigma of the drawn density | 14 of 17 (0.82; Wilson interval 0.60 to 0.93) | 0.880 | 0.8 |

The density is the compaction law's (`compaction.compacted_density_kg_m3`, the tapered brittle law) at the radius of each object and a rock fraction drawn from the row, against the published density with its error carried from the diameter alone (the paper's own reading: the radius dominates, the mass comes from the binary's orbit). The law ends at 1e21 kg (the structure solver decides above it), so the test covers the nine objects up to 705 km across (radius 61 to 353 km); the eight larger ones are outside what the law covers and are flagged. At the row's center (0.70) the covered objects sit within 1.5 sigma except one: the 652 km object (820 kg/m3) is 13 sigma under it. That object is reached only at the row's lower limit, 0.45 (the law gives 819.7 kg/m3 there), which matches the published fit (about 0.5) and is more than 3 sigma under the center of the draw, so the limit stays at 0.45 and the width at 0.08.

## 3. What changed in the rows

Nothing a draw reads. The rows were open (no weight, no source); they now carry the state the data support:

- **`icy-intact`** (Kuiper belt objects, class D bodies and half the dwarf planets): partly verified. The center 0.70 is [R0390] (a constant 70 percent matches 14 of 17 within 2 sigma and 10 of 17 within 1 sigma, with the porosity of the compaction law doing the work in the small objects). The limit text says it is a fit to one system's objects and that the width 0.08 and the limits are the repository's (stand-in CP5, tier T3, tweak L0). Its weight for a dwarf planet is a half, as the open share was.
- **`icy-stripped`** (a dwarf planet whose ice a collision removed): a stand-in (CP6, tier T4, tweak L0), weight a half. The draw's mean and spread match the seven dwarf planets because a half takes the mean from 0.70 to 0.79; the weight rests on no count outside one system (two of the seven, the densest, are the collisional ones, which is also inside 0.3 to 0.6).

`Table.digest` changes (the states and weights are part of it); the draws do not: a dwarf planet drew the two rows at a half before (open rows share what the sourced rows leave), and the other kinds have the one row. No generator epoch bump.

## 4. What the data do not decide

- The two readings of the 17 densities (rock fraction rising with size, against a constant 70 percent with porosity doing the work) cannot be told apart from the table; the row keeps the second, which is the paper's and the compaction law's.
- The small bodies' porosity-free rock fraction is 0.3 to 0.45 (tunable-gaps.md section 3.5); it is a different quantity from the rock fraction of a body at its density with its porosity, and nothing here uses it.
- The width 0.08 and the lower limit 0.45 are the repository's; the three objects the constant 70 percent misses (the two densest, which the law does not cover, and the 652 km object) say the lower tail matters and the upper does not. Neither is tuned to.
- No extrasolar icy-body density exists. Replace the test points with targets when one does (CP5, CP6).
