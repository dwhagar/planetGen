# The iron branch of a rocky planet (GEN.512)

Status: built 2026-10-24 (lane "Source and tune atmosphere draws"). Code `src/planetgen/physics/bulk_draw.py` (`IRON_BRANCH_TABLE`, `iron_branch_ratios`, `branch_ratios`), the tuning sample `src/planetgen/physics/iron_branch_tuning.py`, targets `t5-outlier-share` and `t5-within-share` in `physics/data/tuning_targets.json`, tests `src/tests/test_iron_branch.py`. It is builder B4 of [tunable-gaps.md](tunable-gaps.md) and row T5 of the body tables; the data are read in [tunable-gaps-sources/group2-delivery-retention.md](tunable-gaps-sources/group2-delivery-retention.md) section 2.5.

## 1. What changed

A rocky planet's core was exactly its star's iron: `dry_partition` took the star's Fe/Mg, so two planets of one star had one core mass fraction (CMF) and no planet was ever iron-rich. Now `draw_terrestrial` draws, after the graphite and before the partition, one branch from `IRON_BRANCH_TABLE`:

- **star-matched** (weight 0.99): the planet's Fe/Mg is the star's times 10 to a normal of `IRON_SCATTER_DEX` (0.05 dex), slope 1;
- **iron-enriched** (weight `IRON_ENRICHED_WEIGHT`, 0.01): the CMF of the dry part is uniform from 0.5 to 0.8 (`RANGES["enriched_cmf"]`; the top is the cap of the observed iron-rich rocky planets, [R3545]) and never below the star-matched CMF, and the Fe/Mg that gives it is solved by bisection.

A moon and the classes D, S and C (icy, small, carbon) have no branch and draw nothing for it. A body given both a mass and a radius (the inverse route) infers its CMF from them and has no branch. The draw order is part of the generator epoch. The coreless branch has no count (no planet is known to lack a core) and is not drawn: derive it from the water fraction [R1248] when that is built.

## 2. The tuning (COMPUTED through the harness)

The published statistic is about inferred CMF, not true CMF, so the simulated side measures each planet. `iron_branch_tuning.sampler` draws a host from the generator's own star model, the planet (mass 1 to 10 Earth masses, period 0.5 to 30 d, radius from the CMF, a planet above 1.8 R_E dropped), the measurement the sample would make (`mass_radius_sample.observe`: the model of GEN.505), the CMF inferred from the measured mass and radius, and its 1-sigma error from the two errors. A planet is an outlier when its inferred CMF is more than 2 sigma above the star-matched value, and within when it is inside 1 sigma. The precisely-measured cut is the `mass-radius-sample` selection model.

| Target | Published | Simulated at the shipped knobs (3000 seeds) | Verdict |
|---|---|---|---|
| `t5-outlier-share`: planets 2 sigma more iron-rich than the star | 1 of 21 [R3546] | 0.053 | within 2 sigma (z 0.08) |
| `t5-within-share`: planets within 1 sigma of the star | 15 of 20 [R0541] | 0.61 | within 2 sigma (z -1.3) |

## 3. What the tuning shows, and its limits

- **Noise alone makes about 4 percent of planets outliers** (0.043 with no enriched branch and no scatter), against the published 4.8 percent: the data need no enriched branch. This agrees with Plotnykov's mock test (mass-radius scatter alone reproduces the outliers about 15 percent of the time). The weight is therefore an upper limit; 0.01 is the value at which the simulated share equals 1 of 21, and the data (exact interval 0.001 to 0.238 for 1 of 21) allow up to about 0.2. The widest bracket in the sources is 0 (Brinkman, 0 of 20 distinct) to 0.23 (Adibekyan, 5 of 22, heterogeneous); the newer homogeneous paper wins the centre.
- **The scatter is not measured.** The within-1-sigma share is 0.59 with no scatter and 0.57 at 0.1 dex, both inside the interval of 15 of 20 (0.53 to 0.89): the measurement error dominates. 0.05 dex is a stand-in between none and the 0.1 of CMF (about 0.2 dex) the body master suggested; replace it when a homogeneous sample of more than the 21 planets exists (Plotnykov says 150 are needed).
- **The slope stays 1.** The published slope is 0.94 (+1.02/-1.07), and 1.3 +- 1.0 or 5.6 +- 1.6 by fit on the same data: consistent with 1, not measured. A slope other than 1 needs a pivot (the reference composition the slope turns about) that the generator's star model does not state; not built.
- **Power** (GEN.665): the outlier target refuses a share of 0.30 and the within target refuses 0.5 and 1.0, but neither refuses the low side of the outlier share (the proportion cannot go below 0), and both are weak against small changes.
- **Validity:** the data cover 1 to 1.8 R_E, hosts of 0.8 to 1.2 solar masses and 1 to 10 Earth masses, FGK hosts; the generator draws the branch for every rocky planet, so M-dwarf and heavier planets are flagged extrapolations.
- The coreless branch, the CMF cap of 0.8 as a planet-formation limit (Mercury is only a label here) and the hosts' own CMF are not tuned.

## 4. Stand-ins

The scatter (0.05 dex) and the enriched CMF range (0.5 to 0.8 from the exoplanet cap) are partly verified by the limits above; the weight is tuned. Each is in the rows' `limit` text.

## 5. The larger sample and the residuals of GEN.512 (GEN.820)

GEN.512 tuned on 3,000 simulated planets per case. At 60,000 per case (the same sampler, `iron_branch_tuning.sampler`; z is the simulated share against the published one, in the units of the target's own uncertainty, `tuning_harness.simulate_and_compare`; the published shares are 0.048 outliers and 0.75 within 1 sigma):

| Case (enriched weight / scatter in dex) | Outliers 2 sigma iron-rich | z | Within 1 sigma | z | Both pass |
|---|---|---|---|---|---|
| shipped: 0.01 / 0.05 | 0.047 | -0.01 | 0.602 | -1.38 | yes |
| no scatter: 0.01 / 0 | 0.044 | -0.12 | 0.608 | -1.33 | yes |
| 0.01 / 0.1 | 0.056 | 0.12 | 0.585 | -1.54 | yes |
| 0.01 / 0.2 (the core spread of the [R0541] re-read in BC8) | 0.087 | 0.56 | 0.532 | -2.04 | no: within-1-sigma share too low |
| no enriched branch: 0 / 0.05 | 0.042 | -0.20 | 0.606 | -1.34 | yes |
| 0.05 / 0.05 | 0.072 | 0.34 | 0.586 | -1.53 | yes |
| 0.23 / 0.05 (the 5 of 22 of [R0032]) | 0.176 | 1.82 | 0.509 | -2.25 | no: within-1-sigma share too low |

What it closes:

- **The larger sample (residual 4).** The shipped knobs reproduce both targets at 60,000 planets (z -0.01 and -1.38), so the 3,000-planet tuning was not a lucky draw; a test now holds this at 20,000 planets. The sample cannot get larger than the published N (21 and 20), so the width of the target, not the simulation, limits the tuning from here.
- **The scatter (residual 3).** 0 to 0.1 dex passes both targets and 0.2 dex fails the within-1-sigma share (z -2.04), so the data do not support the 0.2 dex of the [R0541] re-read as an intrinsic spread: the registered range of BC8 stands as an upper bound and 0.05 dex stays a stand-in at tier T6 and tweak L2, with the 0.1 dex edge now measured. The within-1-sigma share sits near z -1.4 at every scatter and with no enriched branch: the simulated sample holds fewer planets within 1 sigma (0.60) than the published 0.75, which the knobs cannot move; it is measurement error, not the branch.
- **The slope (residual 2).** It stays 1: the published slopes are consistent with 1 within their errors and with each other's method (row BC8, fill 3), so no pivot is needed.
- **The coreless branch (residual 1).** No count of coreless planets exists, and the simulated sample reproduces the published shares with none, so nothing is drawn: registry row BC25 (share 0, tier T4, tweak L0) holds it with what would change the answer.
- **The enriched weight** at 0.23, the share of [R0032]'s 5 of 22, is refused by the within-1-sigma target (z -2.25), which agrees with BC8 putting that sample second to the newer homogeneous one.

No generated value changes.

