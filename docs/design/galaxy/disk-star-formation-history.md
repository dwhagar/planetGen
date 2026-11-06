# The thin disk's star formation history from its gas (GEN.469)

**Status (2026-10-17):** built (`galaxy/disk_history.py`, `tests/test_disk_history.py`). Boss chose "Build it" on the
decision card at 5:34 PM EDT (as the lane reported).

The disk's three populations (`young` 0 to 0.1 Gy, `intermediate` 0.1 to 3, `old` 3 to 10) used to share its stars in
proportion to the length of their age windows, as if the disk had formed stars at one rate for ten Gy. That is a typed
assumption. Now each population's share is the part of the disk's own star formation history that fell in its window,
and the history is worked out from the gas.

## The model

One disk, forty rings out to ten scale lengths, each following gas continuity from the disk's start (10 Gy ago):

    d(gas)/dt = inflow(t) - (1 - R) * rate,      rate = A * gas^n          (per area)

| Piece | Value | Source (read) |
|---|---|---|
| Depletion | `A = 2.5e-4` solar masses/yr/kpc2 per (solar masses/pc2)^1.4, `n = 1.4` | [R2321], eq. 4 [R2321], already in `tuning` |
| Return fraction `R` | 0.384 from the generator's own mass function (`IMF_BREAKS_SOL`, `IMF_SLOPES`, the Chabrier function of `stellar-mass-function.md`; 0.371 on the Kroupa function before GEN.927); the same code gives 0.272 for a Salpeter slope of 2.35 from 0.1 to 100 solar masses, Madau and Dickinson's 0.27 | [R2798], section 2.1: instantaneous recycling above 1 solar mass, white dwarf remnants by the [R0961] relation of `tuning.WD_IFMR_PIECES` (the old relation of `0.444 + 0.084 m` before GEN.930), 1.4 solar masses to 40, black holes keep all |
| Infall | Exponential, `exp(-t / tau(R))` from the disk's start; `tau(R) = 1.033 R - 1.267` Gy (7 Gy at the Sun) between 2 and 8 Gy, `R` the Milky Way's radius of the ring's scale lengths (2.6 kpc each) | [R5095] section 3 (the line), [R5094] section 2.4 (2 Gy inside, 8 Gy at 20 kpc) |

**The inflow** is an exponential infall whose time scale grows outward, so the inner rings are fed early and the outer ones keep being fed (inside-out
formation). Only its shape in time matters: each ring's amplitude is solved (by bisection, all rings at once) so the ring ends with its share of the
disk's stars, `Sigma*(R) = M_thin / (2 pi L^2) exp(-R/L)`. The thin disk's mass is the galaxy's stellar mass (the mass chain's, drawn from the seed) times the thin
disk's share of the density model's stars. A dense ring has a short depletion time and follows the inflow closely; a thin
one lags it. The populations' shares are the rings' star formation (weighted by ring area) inside each age window.

A shape carries the galaxy's `(stellar mass, halo mass)` in a new field `GalaxyShape.masses`. It is not stored: it is
drawn again from the seed (`db.store.get_galaxy_shape`, `shape_draw`), like the component field. A shape built without
masses (tests, toy galaxies) takes the median galaxy of its disk scale length ([R4048] sizes, the mass chain's
stellar-to-halo relation).

## What it gives

Disk-wide, for a disk scale length of 2.6 kpc (`tests/test_disk_history.py`), since GEN.1230 (the halo is no longer read):

| Stellar mass | Young (0 to 0.1 Gy) | Intermediate (0.1 to 3) | Old | Present rate over the mean | Local B-star share (between arms, 2026-10-18, halo-growth inflow) |
|---|---|---|---|---|---|
| flat history | 1.00% | 29.0% | 70.0% | 1 | 0.10% |
| 1e10 | 0.70% | 23.9% | 75.4% | 0.70 | 0.26% |
| 3e10 | 0.54% | 19.5% | 80.0% | 0.54 | 0.13% |
| 5e10 | 0.49% | 17.7% | 81.8% | 0.48 | 0.050% |
| 6.1e10 | 0.47% | 17.1% | 82.4% | 0.46 | 0.033% |

(Before GEN.1230 the inflow followed the stars a galaxy's halo holds and the rows read 2.93, 1.20, 0.39 and 0.27 percent young and a present rate over the mean of
2.95, 1.18, 0.38 and 0.27; the lighter halo of GEN.1154 moved the 5e10 row to 0.55, and the Milky Way shape's star formation to 4.0 solar masses a year against the 1.65 measured. The B-star column is from then and has not been measured again.)

The last column was recomputed on 2026-10-18 with the generator's own sampler, 100,000 draws a population, the flat row with
ages spread evenly inside the windows and the others with the ages spread by the history's rate (GEN.546, below). GEN.469
recorded 0.19, 0.38, 0.26, 0.12 and 0.085% for the same rows with even ages; the same measure gives 0.10, 0.20, 0.14,
0.069 and 0.047% today, so the sampler or the mass function has moved since.

A light spiral forms its stars later than a heavy one (downsizing), but since GEN.1230 only through the law: a thinner disk spends its gas more slowly,
so the young share is 1.4 times higher at 1e10 than at 5e10 (3 times with the halo-growth inflow, which also made the Milky Way's disc form 2.7 times the stars
it does). The infall time scale is the Milky Way's in scale lengths for every galaxy (row GI1). The mass dependence is a factor on the inflow, `exp(0.4 t (1 - V / 220 km/s))` per gigayear ([R5121]: equation 21), with the disc speed V = 220 (M / 5e10)^(1/3) km/s from a disc mass that goes as V cubed (row GI3). A disc slower than the Milky Way's is fed ever more slowly into the present and a faster one the other way, which restores the downsizing the infall alone had flattened: [R5122] has a 100 km/s disc rising to the present, a 200 km/s one about constant and a 360 km/s one falling (GEN.1246).
A galaxy's young share is still a property of the galaxy, so the magnetar rate (which follows the young stars) differs from galaxy to galaxy:
`MAGNETAR_YOUNG_STAR_SHARE` is now the share of a galaxy of 5e10 solar masses, 0.18% (it was 0.16%, and 0.44% before GEN.469).

**Test targets, none an input.** The Galaxy forms 1.65 +- 0.19 solar masses a year ([R2680, R2681], read; M* 6.08 +- 1.14e10, disk 5.17 +- 1.11e10), which against a ten-Gy mean of
5.17e10 / (1 - R) / 10 Gy is a present rate 0.20 of the mean. The derived 0.48 (5e10) and 0.46 (6.1e10) are 2.3 times it (the test allows 2.5; the observed figure counts the thick disc's stars over
the thin disc's ten Gy, so it is low for the thin disc alone). The whole Milky Way shape forms 2.8 solar masses a year now against the measured 1.65 (the census range is 1 to 3):
the outer disc (beyond 10 kpc) and the inner disc (inside 4 kpc) still form more than observed, where the model has no truncation and no bar-depleted centre. The gas left at the Sun's radius is 0.2 to 0.3 of its stars (observed about 13 against 38 solar masses per
pc2), inside the test's 0.1 to 0.5.

## The law steepens where the gas thins (GEN.999)

The rings used the plain Kennicutt-Schmidt law, so a thin outer ring spent its gas as fast as a dense one and a Milky Way mass disc ended with about a third of the gas it has: a depletion time of 1.8 Gy (3e10 solar masses of thin disc stars) against the 1.0 to 4.2 Gy measured (Sources 2 on GEN.998), and 4.3 solar masses per pc2 at the Sun's radius. The ring solve now carries the factor the present-day profile already used (GEN.450, `star_formation._thin_gas_factor`): `rate = A * gas^n * (gas / (gas + gas_s))^(m - n)`, with the saturation `gas_s` of 9 solar masses per pc2 [R0392] and the index `m` of 2.1 (Kennicutt and Evans, the middle of 1.5 to 2.7), both already in `tuning`. No new constant, and no Milky Way value is fitted.

| Milky Way shape (`tests/test_disk_history.py`) | Plain law | With the factor | Measured |
|---|---|---|---|
| Depletion time, thin disc of 3e10 stars | 1.8 Gy | 2.9 Gy | 1.0 to 4.2 Gy |
| Depletion time, the full Milky Way shape | 2.2 Gy | 3.7 Gy | 4.2 Gy (Sources 2) |
| Gas at 8 kpc | 4.3 per pc2 | 7.2 per pc2 | 11.9 +- 1.6 atomic and molecular with helium ([R3002: Table 2]: 10.9 + 1.0 within 1 kpc of the Sun; 13.7 with the ionised gas) |
| Star formation now | 1.50 a year | 1.70 a year | 1.65 +- 0.19 [R2680] |
| Gas inside the bar gap | 5.0e8 | 6.2e8 | (the ring inflow follows it) |

The gas at the Sun's radius is still 60% of the measured, and raising it further (a higher index or saturation) would push the depletion time and the present rate past their measured ranges; the plain law's coefficient ([R2321], +- 28%) is a mean over spirals, and the Milky Way is known to sit below it at the Sun's radius. That residual is left as it is and filed with the stand-in row AL1.

Populations' shares move with the rate: for stellar masses 1e10, 3e10, 5e10 and 6.1e10 the young share goes 2.99, 1.31, 0.54 and 0.37% to 3.42, 1.43, 0.60 and 0.42%, the intermediate 62.6, 44.4, 25.0 and 18.1% to 68.0, 47.4, 26.8 and 19.5%, and the old falls to match. The table under "What it gives" is as of GEN.469 and is not redone here.

(Since GEN.1230 the same disc depletes in 3.4 Gy and keeps 10.3 solar masses per pc2 at 8 kpc; the table above is the halo-growth inflow's.)

## Ages inside a window (GEN.546)

A window used to draw its stars' ages evenly, though the rate varies across it (by a factor of 1.0, 2.5 and a rise from
nothing for the young, intermediate and old windows of a 5e10 galaxy). Now the age of a young, intermediate or old disk
star is drawn from the history's own rate inside its window.

- **The law.** `galaxy.age_law.AgeLaw` is the density of ages in a window as a piecewise-linear curve through the disk's
  star formation rate (rings weighted by area) at 33 evenly spaced ages (`tuning.AGE_LAW_POINTS`, rounded to six figures so
  that equal histories give one law). It gives the weight of an interval (the window's length times the chance an age
  falls in it, which for an even spread is the overlap), and an age by inverting its cumulative.
  `galaxy.disk_history.population_age_laws(shape)` builds one per disk population (`tuning.DISK_AGE_LAW_POPULATIONS`)
  for the galaxy of a shape. The bulge, halo, clusters and satellites have no entry and draw evenly: their history is not
  the disk's.
- **One law everywhere.** The star samplers' tables, shares and draws all take the same law (`age_laws` in
  `generation.star_population`, `bright_stars.age_laws_of(shape)`), so the bright stars the scatter places, the mass
  pass, both backfills and the dim systems a sector draws for itself agree on how many stars of each kind a place holds
  and how old they are. `FillContext.apply` hands a dim system its population's law as `SystemConfig.AGE_LAW`. A draw with
  no population (a system made outside a galaxy) keeps the whole disk's even spread.
- **What changes.** Within a window the ages shift toward where the rate is high: the old window of a 5e10 galaxy
  averages 6.1 Gy instead of 6.5, the intermediate one 1.8 Gy instead of 1.55, and a 1e10 galaxy's old window 5.0 Gy.
  The shares of each window are still the history's (GEN.469), but the share of bright stars in a window changes with the
  ages: a front-loaded intermediate window holds more bright main-sequence stars.
- **The local B-star share** between the arms at the solar radius, with 100,000 draws a population, falls from 0.069% to
  0.050% (0.046 to 0.055% over four seeds), against the observed 0.04%; the table above shows the effect for each
  galaxy mass.

## What it does not do

- **The observed local B-star share (0.04%, Mamajek's census) is still 1.2 to 1.4 times lower than the Milky Way-like
  galaxy's 0.046 to 0.055% (GEN.546 closed most of the gap, below).** The test now allows 1.5 times the observed share.
  What is not in the model and could close the rest: the single-star mass function may be steeper than the cluster one above 1
  solar mass ([R2498], astro-ph/0009005, on the bias of unresolved binaries), the Sun's place above the plane and
  between the arms, and gas moving between rings.
- **Gas moves between rings and flows out** in a real galaxy ([R5094] has bar-driven radial flows and none outflows; [R5095] has neither); here each ring is a closed
  box fed by its inflow.
- **The radial pattern stays GEN.134's** (`galaxy.star_formation.recent_factor`, from the gas disk's scale length and the
  bar gap). The rings alone make the outer disk younger, the opposite of the observed fall; the factor still multiplies
  the young and intermediate shares and averages 1 over the disk, so the disk-wide shares are the history's.
- **The infall time scale is a fit to one galaxy.** [R5095] fits it to the Milky Way's abundance gradients and [R5094] takes 2 to 8 Gy; neither gives it for another size or mass, so every galaxy takes
  the Milky Way's in scale lengths (GEN.725 row GI1, level L0).

## Code

- `galaxy/disk_history.py`: infall time scale, ring inflows, ring solver, shares.
- `galaxy/star_formation.disk_population_share(name, shape)`: the share `density.population_densities`,
  `population_densities_array` and `generation/object_first._disk_share_bound` read (the bound must use the same share
  as the density it bounds).
- `tuning`: the constants above with their sources.
- A change of the generator's output: epoch bumped (OPS.56); it joins the one full rescatter (GEN.301 waits).

## The thick disc's ages in the old population (GEN.1025)

The thick disc is counted in the `old` population (`density.population_units`), but the old window ended at 10 Gy, the thin disc's age, while the thick disc formed from about 13 to about 8 Gy ago ([R1863], `tuning.THICK_DISC_AGE_RANGE_GY`, stand-in TD1). `galaxy/thick_ages.py` builds the old population's `AgeLaw` as a mixture: the thin disc's law over 3 to 10 Gy with the weight the thin disc has among the old systems, and the thick disc's law, flat over 8 to 13 Gy (stand-in TD2), with the thick disc's share (`density.thick_disk_units` over the old units). Only the shape of an age law matters, so the count of old stars per sector does not change; the stars are spread over ages up to 13 Gy, and the thick disc's share of them is older than the thin disc ever was.

Wired (part 2): `disk_history.population_age_laws` returns the mixture for `old`, `tuning.STELLAR_POPULATION_AGE_RANGES_GY["old"]` ends at 13 Gy, and `tuning.THIN_DISK_OLD_WINDOW_GY` (3 to 10) keeps the thin disc's own window where the thin disc alone is meant: its share of the disc (`population_shares`), its law, its heating (`galaxy/structure.py`) and its arm contrast. The mean system mass of `old` and the remnants formed in it follow the longer window, so the density identity fixture (`tests/fixtures/density_identity.json`) was made again: every density of a shape moves by one number (0.04 to 1.2 percent, the model's units per solar mass), the total mass not at all. Changes generated output: the old stars of a sector are drawn over 3 to 13 Gy.

## The infall replaces the halo growth (GEN.1230, GEN.1242)

The inflow followed the stars a galaxy's halo holds. GEN.1154 made the median halo of a 5e10 galaxy lighter (2.0e12 to 1.37e12 solar masses), which
moved the growth of those stars later, so the Milky Way shape's disc formed 4.0 solar masses a year (3.7 before GEN.999's gas law, 2.0 before GEN.1154) against the 1.65
measured, and 1.5 against the main sequence at that mass (Speagle's fit: 2.7 times above it). Weakening the law's coefficient made it worse (the stars are fixed, so more form
late), a cold-accretion cut at the shock-heating mass needed an exponent no source gives, and the halo's own baryon accretion gave 6.3. The one-zone Milky Way models
use an exponential infall: the line of [R5095] for its time scale, limited by [R5094]'s 2 and 8 Gy. It gives 2.8 solar masses a year for the Milky Way shape, a gas surface density
of 10.3 solar masses per pc2 at 8 kpc (measured 7 to 16) and a depletion time of 3.4 Gy at the Sun. Registry rows GI1 (the scale length mapping) and GI2 (the 2 and 8 Gy stops); the core gas
inflow time went from 1.0 to 0.6 Gy (CG1) because the inner rings keep less gas.
