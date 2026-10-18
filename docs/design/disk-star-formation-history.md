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
| Depletion | `A = 2.5e-4` solar masses/yr/kpc2 per (solar masses/pc2)^1.4, `n = 1.4` | Kennicutt 1998, ApJ 498, 541, eq. 4 (arXiv astro-ph/9712213), already in `tuning` |
| Return fraction `R` | 0.371 from the generator's own mass function (`IMF_BREAKS_SOL`, `IMF_SLOPES`); the same code gives 0.272 for a Salpeter slope of 2.35 from 0.1 to 100 solar masses, Madau and Dickinson's 0.27 | Madau and Dickinson 2014, ARA&A 52, 415 (arXiv 1403.0007), section 2.1: instantaneous recycling above 1 solar mass, white dwarf remnants `0.444 + 0.084 m` (Weidemann 2000), 1.4 solar masses to 40, black holes keep all |
| Halo growth | `46.1 Msun/yr (M / 1e12)^1.1 (1 + 1.11 z) E(z)`, integrated back from the halo's mass today | Fakhouri, Ma and Boylan-Kolchin 2010, MNRAS 406, 2267 (arXiv 1001.2304), eq. 2 |
| Stars a halo holds | Stellar-to-halo relation with its four parameters moving as `z / (z + 1)` (log M1 11.590 + 1.195, N 0.0351 - 0.0247, beta 1.376 - 0.826, gamma 0.608 + 0.329) | Moster, Naab and White 2013, MNRAS 428, 3121 (arXiv 1205.5807), eqs. 11 to 14 and Table 1 |
| Cosmology | h 0.671, Omega_m 0.3175 (the mass chain's; the Universe is 13.83 Gy old against `UNIVERSE_AGE_GY` 13.8) | `galaxy.mass_chain`; Planck 2018 (arXiv 1807.06209) gives 67.4 and 0.315 |

**The inflow** is the growth of the stellar mass the galaxy's halo holds along its growth track, never below zero. The
gas comes in as fast as the stars it must build. Only its shape in time matters: each ring's amplitude is solved (by
bisection, all rings at once) so the ring ends with its share of the disk's stars, `Sigma*(R) = M_thin / (2 pi L^2)
exp(-R/L)`. The thin disk's mass is the galaxy's stellar mass (the mass chain's, drawn from the seed) times the thin
disk's share of the density model's stars. A dense ring has a short depletion time and follows the inflow closely; a thin
one lags it. The populations' shares are the rings' star formation (weighted by ring area) inside each age window.

A shape carries the galaxy's `(stellar mass, halo mass)` in a new field `GalaxyShape.masses`. It is not stored: it is
drawn again from the seed (`db.store.get_galaxy_shape`, `shape_draw`), like the component field. A shape built without
masses (tests, toy galaxies) takes the median galaxy of its disk scale length (Shen et al. 2003 sizes, the mass chain's
stellar-to-halo relation).

## What it gives

Disk-wide, for a disk scale length of 2.6 kpc and the median halo of each stellar mass (`tests/test_disk_history.py`):

| Stellar mass | Halo | Young (0 to 0.1 Gy) | Intermediate (0.1 to 3) | Old | Present rate over the mean | Local B-star share (between arms) |
|---|---|---|---|---|---|---|
| flat history | | 1.00% | 29.0% | 70.0% | 1 | 0.19% |
| 1e10 | 3.2e11 | 2.93% | 62.0% | 35.1% | 2.95 | 0.38% |
| 3e10 | 8.3e11 | 1.20% | 42.3% | 56.5% | 1.18 | 0.26% |
| 5e10 | 2.0e12 | 0.39% | 19.0% | 80.6% | 0.38 | 0.12% |
| 6.1e10 | 3.1e12 | 0.27% | 12.7% | 87.1% | 0.27 | 0.085% |

A light spiral has a rising history and forms its stars late; a heavy one formed them early. That is the known pattern
(downsizing) and it comes out of the halo's growth and the relation, with nothing typed. It also makes a galaxy's
young share a property of the galaxy, so the magnetar rate (which follows the young stars) differs from galaxy to galaxy:
`MAGNETAR_YOUNG_STAR_SHARE` is now the share of a galaxy of 5e10 solar masses, 0.19% (it was 0.44%).

**Test targets, none an input.** The Galaxy forms 1.65 +- 0.19 solar masses a year (Licquia and Newman 2015,
arXiv 1407.1078, read; M* 6.08 +- 1.14e10, disk 5.17 +- 1.11e10), which against a ten-Gy mean of
5.17e10 / (1 - R) / 10 Gy is a present rate 0.20 of the mean. The derived 0.38 (5e10) and 0.27 (6.1e10) are within a
factor of two. The gas left at the Sun's radius is 0.2 to 0.3 of its stars (observed about 13 against 38 solar masses per
pc2), inside the test's 0.1 to 0.5.

## What it does not do

- **The observed local B-star share (0.04%, Mamajek's census) is still 3 times lower than the Milky Way-like
  galaxy's 0.12% (0.085% at 6.1e10 solar masses).** `test_the_local_b_share_reaches_the_observed_one` is a strict expected failure that records this; it
  turns into a failure to remove the day something closes the gap. What is not in the model and could: the single-star
  mass function may be steeper than the cluster one above 1 solar mass (Kroupa 2001, astro-ph/0009005, on the bias of unresolved
  binaries: a
  slope of 2.7 instead of 2.3 would cut B stars by about 40%, but the generator's mass function is the one the whole star
  model reads), the Sun's place above the plane and between the arms, and gas moving between rings.
- **Ages stay uniform inside a window.** The star formation rate varies across a window (by a factor of 1.0, 2.5 and a
  rise from nothing for the young, intermediate and old windows of a 5e10 galaxy), so the mean ages are somewhat off.
  Drawing an age from the history inside the window is a follow-up in the sampler.
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
