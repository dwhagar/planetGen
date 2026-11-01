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
| Halo growth | `46.1 Msun/yr (M / 1e12)^1.1 (1 + 1.11 z) E(z)`, integrated back from the halo's mass today | Fakhouri, Ma and [R0508] [R1313], eq. 2 |
| Stars a halo holds | Stellar-to-halo relation with its four parameters moving as `z / (z + 1)` (log M1 11.590 + 1.195, N 0.0351 - 0.0247, beta 1.376 - 0.826, gamma 0.608 + 0.329) | [R3181], eqs. 11 to 14 and Table 1 |
| Cosmology | h 0.671, Omega_m 0.3175 (the mass chain's; the Universe is 13.83 Gy old against `UNIVERSE_AGE_GY` 13.8) | `galaxy.mass_chain`; [R3544] gives 67.4 and 0.315 |

**The inflow** is the growth of the stellar mass the galaxy's halo holds along its growth track, never below zero. The
gas comes in as fast as the stars it must build. Only its shape in time matters: each ring's amplitude is solved (by
bisection, all rings at once) so the ring ends with its share of the disk's stars, `Sigma*(R) = M_thin / (2 pi L^2)
exp(-R/L)`. The thin disk's mass is the galaxy's stellar mass (the mass chain's, drawn from the seed) times the thin
disk's share of the density model's stars. A dense ring has a short depletion time and follows the inflow closely; a thin
one lags it. The populations' shares are the rings' star formation (weighted by ring area) inside each age window.

A shape carries the galaxy's `(stellar mass, halo mass)` in a new field `GalaxyShape.masses`. It is not stored: it is
drawn again from the seed (`db.store.get_galaxy_shape`, `shape_draw`), like the component field. A shape built without
masses (tests, toy galaxies) takes the median galaxy of its disk scale length ([R4048] sizes, the mass chain's
stellar-to-halo relation).

## What it gives

Disk-wide, for a disk scale length of 2.6 kpc and the median halo of each stellar mass (`tests/test_disk_history.py`):

| Stellar mass | Halo | Young (0 to 0.1 Gy) | Intermediate (0.1 to 3) | Old | Present rate over the mean | Local B-star share (between arms, 2026-10-18) |
|---|---|---|---|---|---|---|
| flat history | | 1.00% | 29.0% | 70.0% | 1 | 0.10% |
| 1e10 | 3.2e11 | 2.93% | 62.0% | 35.1% | 2.95 | 0.26% |
| 3e10 | 8.3e11 | 1.20% | 42.3% | 56.5% | 1.18 | 0.13% |
| 5e10 | 2.0e12 | 0.39% | 19.0% | 80.6% | 0.38 | 0.050% |
| 6.1e10 | 3.1e12 | 0.27% | 12.7% | 87.1% | 0.27 | 0.033% |

The last column was recomputed on 2026-10-18 with the generator's own sampler, 100,000 draws a population, the flat row with
ages spread evenly inside the windows and the others with the ages spread by the history's rate (GEN.546, below). GEN.469
recorded 0.19, 0.38, 0.26, 0.12 and 0.085% for the same rows with even ages; the same measure gives 0.10, 0.20, 0.14,
0.069 and 0.047% today, so the sampler or the mass function has moved since.

A light spiral has a rising history and forms its stars late; a heavy one formed them early. That is the known pattern
(downsizing) and it comes out of the halo's growth and the relation, with nothing typed. It also makes a galaxy's
young share a property of the galaxy, so the magnetar rate (which follows the young stars) differs from galaxy to galaxy:
`MAGNETAR_YOUNG_STAR_SHARE` is now the share of a galaxy of 5e10 solar masses, 0.19% (it was 0.44%).

**Test targets, none an input.** The Galaxy forms 1.65 +- 0.19 solar masses a year ([R2680, R2681], read; M* 6.08 +- 1.14e10, disk 5.17 +- 1.11e10), which against a ten-Gy mean of
5.17e10 / (1 - R) / 10 Gy is a present rate 0.20 of the mean. The derived 0.38 (5e10) and 0.27 (6.1e10) are within a
factor of two. The gas left at the Sun's radius is 0.2 to 0.3 of its stars (observed about 13 against 38 solar masses per
pc2), inside the test's 0.1 to 0.5.

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
- **Gas moves between rings, flows out, and halos gain stars by mergers** in a real galaxy; here each ring is a closed
  box fed by its inflow and the relation's stellar mass counts every star the halo holds.
- **The radial pattern stays GEN.134's** (`galaxy.star_formation.recent_factor`, from the gas disk's scale length and the
  bar gap). The rings alone make the outer disk younger, the opposite of the observed fall; the factor still multiplies
  the young and intermediate shares and averages 1 over the disk, so the disk-wide shares are the history's.
- **The relation's stellar mass is that of the central galaxy**, bulge, thick disk and thin disk together. The thin disk
  takes the whole of it over its ten Gy, which puts a little too much of the early growth into the thin disk.
- The halo's mass today is the median for the galaxy's stellar mass with the chain's scatter; Moster's relation puts a
  galaxy of the Milky Way's mass in a halo about twice as heavy as dynamical estimates.

## Code

- `galaxy/disk_history.py`: cosmology, halo track, inflow, ring solver, shares.
- `galaxy/star_formation.disk_population_share(name, shape)`: the share `density.population_densities`,
  `population_densities_array` and `generation/object_first._disk_share_bound` read (the bound must use the same share
  as the density it bounds).
- `tuning`: the constants above with their sources.
- A change of the generator's output: epoch bumped (OPS.56); it joins the one full rescatter (GEN.301 waits).

## The thick disc's ages in the old population (GEN.1025)

The thick disc is counted in the `old` population (`density.population_units`), but the old window ended at 10 Gy, the thin disc's age, while the thick disc formed from about 13 to about 8 Gy ago ([R1863], `tuning.THICK_DISC_AGE_RANGE_GY`, stand-in TD1). `galaxy/thick_ages.py` builds the old population's `AgeLaw` as a mixture: the thin disc's law over 3 to 10 Gy with the weight the thin disc has among the old systems, and the thick disc's law, flat over 8 to 13 Gy (stand-in TD2), with the thick disc's share (`density.thick_disk_units` over the old units). Only the shape of an age law matters, so the count of old stars per sector does not change; the stars are spread over ages up to 13 Gy, and the thick disc's share of them is older than the thin disc ever was.

Wired (part 2): `disk_history.population_age_laws` returns the mixture for `old`, `tuning.STELLAR_POPULATION_AGE_RANGES_GY["old"]` ends at 13 Gy, and `tuning.THIN_DISC_OLD_WINDOW_GY` (3 to 10) keeps the thin disc's own window where the thin disc alone is meant: its share of the disc (`population_shares`), its law, its heating (`galaxy/structure.py`) and its arm contrast. The mean system mass of `old` and the remnants formed in it follow the longer window, so the density identity fixture (`tests/fixtures/density_identity.json`) was made again: every density of a shape moves by one number (0.04 to 1.2 percent, the model's units per solar mass), the total mass not at all. Changes generated output: the old stars of a sector are drawn over 3 to 13 Gy.
