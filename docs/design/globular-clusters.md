# Globular clusters in generation

Boss's research of 2026-10-09 (18:45Z, "Procedural Generation and Dynamical Modeling of Globular Clusters in Galactic-Scale Architectures"), digested against how the generator works today, with the numbers checked by calculation. It answers Boss's question of 18:32Z, "What do we need to know in order to model globular clusters in generation?", and settles the decision he was asked: store a cluster as a density component, not as hand-placed stars.

Informs: GEN.9 and VIEW.2 (other galaxies), GEN.100 and GEN.130 (the captured-giant pulsar planet, left out until clusters exist), GEN.134 (the metallicity gradient, item 4), GEN.84 and GEN.133 (stellar populations), GEN.42 and PERF.34 (sampling and cost)
Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation. No cluster model exists in the repository today.

Evidence tags: [B] a claim from Boss's research text, its cited paper not opened by this digest (the research environment cannot reach arXiv or journal sites); [C] computed here, re-runnable with `/mnt/project-files/research/scripts/globular/calc.py` (numpy only, seeded); [S] read from the repository; [R] recalled and unconfirmed.

Related notes: [galaxy-disk-density.md](galaxy-disk-density.md) (the density model a cluster adds to), [sampling-backfill-and-resume.md](sampling-backfill-and-resume.md) (how a sector is filled), [orbital-updates.md](orbital-updates.md) (the velocity and epoch design a cluster's orbit uses), [compact-remnant-regions.md](compact-remnant-regions.md) (millisecond pulsars), [exotic-environments-planets-and-compact-binaries.md](exotic-environments-planets-and-compact-binaries.md) (pulsar-planet types), [star-types-by-galactic-radius.md](star-types-by-galactic-radius.md) (GEN.134), [multiple-galaxies.md](multiple-galaxies.md) (the galaxies that get a cluster system), [performance-eta-queue-and-caching.md](performance-eta-queue-and-caching.md) (cost per system).

## Summary

**What the generator needs to know, in one list.** (1) A cluster record: position, velocity, mass, King concentration, half-mass radius, age, metallicity, a core-collapse flag. (2) A way to turn it into stars on demand. (3) A metallicity value, which the generator does not have at all. (4) A rule that culls planets by local density and metallicity. (5) A count and mass function for clusters in a generated galaxy. Boss's text covers all five; the checks below confirm the central choice and correct five points.

1. **Density component, not stored stars: confirmed by arithmetic [C].** A Milky Way-like system of 157 clusters holds about 1.8e8 systems, 0.2% of the galaxy's 1e11. Generating them all at the measured 5 to 17 ms per system is 10 to 34 days on one worker. It is also 262,000 sectors (mean 1,670 per cluster) that would otherwise be empty. Only the bright-first backfill (generate the giants and blue stragglers, fill the dim rest on demand) is affordable.
2. **"Millions of stars in a 4 pc sector" is true only for the most massive clusters [C].** A 47 Tucanae-like cluster (8e5 Msun, half-mass radius 6 pc) puts 1.1e5 to 5.5e5 systems in its central sector depending on where the sector edge falls; a typical one (2e5 Msun, 5 pc) puts 3e4 to 8e4. Only omega Centauri-class clusters (3.6e6 Msun [R]) reach millions. The 4 pc sector holds 6.3 systems at local density, so a 47 Tuc core sector is about 2e4 to 9e4 times denser.
3. **Use King models, not Plummer [C].** Plummer's inverse is exact (the half-mass radius is 1.305 a; sampled and analytic mass fractions agree to 0.0004), but it has no edge: a typical cluster fills 5,245 sectors with at least one system against 965 for a King model with the same mass and half-mass radius. A King model is one table per concentration `W0` scaled by length, so no per-cluster table is stored (section 3).
4. **The metallicity range is too narrow [R].** Boss's text gives [Fe/H] from -2.3 to -0.5. The Milky Way's bulge clusters run to about 0 (NGC 6528 and 6553 near -0.1 and -0.2 [R]), and -2.3 is 0.5% of solar [C], metal-poor, not "nearly pristine". Sample the two-peak mixture but do not cap it at -0.5.
5. **Blue stragglers are not linear in the encounter rate.** The text's "Dynamic Filtering" step says to scale blue stragglers and millisecond pulsars "linearly with Γ", but its own section 6 says blue stragglers scale as `M_core^0.4 to 0.5`. The first is right for X-ray binaries and millisecond pulsars [R, Pooley 2003; Hui 2010], the second for blue stragglers [R, Knigge 2009 gives 0.38]; they are different laws (section 5).
6. **The cluster's own orbit is slow on game time scales [C].** At 200 km/s a cluster crosses a 4 pc sector in 19,600 years, so stars stored as offsets from the cluster centre stay valid for a century and more; the sector address of a cluster star is derived from the centre, not stored per star (section 6).

**Recommendation.** Build in four slices (section 9): the metallicity value and a cluster table with the real Milky Way catalogue; the sampling and the bright-first fill; the planet cull and pulsar planets; the synthetic cluster systems for other galaxies. Slice 1 needs Boss's decision on one input, the catalogue file (section 8).

## 1. What Boss's research says

- **Storage.** A cluster is a metadata node (position, mass, radii, metallicity) that the sector engine samples; the Milky Way's roughly 160 clusters are loaded from the Harris 1996/2010 and Baumgardt and Hilker catalogues; other galaxies get a synthetic system [B].
- **Count.** `N_GC = S_N 10^(-0.4 (M_V + 15))`; `S_N` about 1 for spirals, about 5 for ellipticals, above 10 for central giants [B].
- **Mass function.** Lognormal, turnover near 2e5 Msun (`M_V` about -7.4) in nearly all galaxies [B].
- **Metallicity.** Bimodal, peaks near -1.5 (halo, extended, accreted, high dispersion) and -0.5 (bulge or thick disc, concentrated, rotating, formed in situ) [B].
- **Structure.** Core radius about 1 pc, half-mass radius 3 to 10 pc, tidal radius 30 to 100 pc [B]. Plummer: `rho = 3M/(4 pi a^3) (1 + r^2/a^2)^(-5/2)`, `r = a (X^(-2/3) - 1)^(-1/2)`. King: dimensionless central potential `W0`, concentration `c = log10(r_t/r_c)`, a precomputed table of cumulative mass against radius. About 20% of Milky Way clusters are core collapsed, with a central cusp `rho ~ r^(-0.7 to -1.3)` [B].
- **Stars.** Age 10 to 13 Gyr, turn-off near 0.85 Msun; stars on the giant branches and the horizontal branch; metal-poor stars hotter and brighter at the same mass; horizontal-branch colour set by metallicity and a second parameter [B].
- **Blue stragglers.** From collisions and binary mass transfer; `Γ ~ rho_c^2 r_c^3 / sigma`; `N_BSS ~ M_core^0.4 to 0.5` [B].
- **Planets.** Giant planet occurrence rises with metallicity and is suppressed below [Fe/H] -1; wide orbits (above about 10 AU) are soft and ionised, tight orbits (below 1 to 2 AU) are hard and survive; millisecond pulsars are 10 to 100 times over-represented; PSR B1620-26 b in M4 is the captured-circumbinary example [B].
- **Motion.** The cluster follows an eccentric halo orbit; its stars add a Maxwellian dispersion of 1 to 15 km/s [B].

## 2. What the generator has today [S]

- Density is one smooth analytic function (`galaxy/density.py`): a bar bulge, a thin disc with spiral arms, a thick disc and a halo floor. `population_densities` splits it into four populations, `young`, `intermediate`, `old` and `bulge`, each with an age range (`tuning.STELLAR_POPULATION_AGE_RANGES_GY`: bulge 8 to 12 Gyr, old 3 to 10 Gyr). A sector's expected system count is 6.3 times `relative_density` at its centre; a sector is generated only if that is at least 1.
- `sample_living_star` draws mass and age together and redraws until the star has not collapsed; with a population it uses that population's age range. The main-sequence lifetime law gives 10 Gyr at 1 Msun and 13 Gyr at 0.9 Msun [S, run]. So a cluster population aged 10 to 13 Gyr removes the heavy stars by itself, with no separate turn-off rule: at 12 Gyr the turn-off is 0.93 Msun in the code, against Boss's 0.85 [B] for metal-poor stars. The 8% gap is metallicity (metal-poor stars burn faster [R]), which the code does not model.
- There is no metallicity anywhere (GEN.134, item 4, defers a gradient "only if planet occurrence is later tied to it"). Giant planets today do not depend on it. This note is that reason.
- Bright stars are placed first and the dim ones filled later by luminosity level (GEN.44, `star_population.py`); the cost is 5 ms per system in memory and 16.6 ms per system in a filled sector [S, performance note].
- Boss decided at 18:37Z that the captured-giant pulsar-planet type (type B) stays out until a cluster model exists [S, TODO GEN.130].

## 3. Storing a cluster: the density component

**Record.** One row per cluster: galaxy, name, source (`catalogue` or `generated`), position and velocity in the galaxy frame (the velocity design of orbital-updates.md), epoch, total mass, half-mass radius, King `W0`, age, [Fe/H], core-collapse flag, and the tidal radius derived from them. Everything else is computed.

**King tables [C].** A King model is self-similar: for a fixed `W0` the density, cumulative mass and tidal radius are one curve in units of the King radius `r0`. A solver that integrates `W'' + 2W'/x = -9 rho(W)/rho(W0)` reproduces the textbook concentrations (`c` = 0.67, 1.03, 1.53, 2.12, 2.74 for `W0` = 3, 5, 7, 9, 12), so ten tables for `W0` = 3 to 12 (integrated on a fine grid, then resampled to about 2,000 points on a log grid) cover every cluster, built once and scaled by `r0 = r_h / (r_h in units of r0)`.

| `W0` | `c` | `r_t / r0` | `r_h / r0` |
|---|---|---|---|
| 3 | 0.67 | 4.7 | 1.26 |
| 5 | 1.03 | 10.7 | 2.00 |
| 7 | 1.53 | 33.7 | 3.92 |
| 9 | 2.12 | 131 | 15.4 |
| 12 | 2.74 | 548 | 86.4 |

**What clusters look like in the sector grid [C].** Four exemplars, `MEAN_M` = 0.4 Msun per system:

| Cluster | Systems | Central density (Msun/pc^3) | King radius | `r_t` | Central sector (systems) | Sectors with at least 1 system |
|---|---|---|---|---|---|---|
| typical (2e5 Msun, `r_h` 5 pc, `W0` 5) | 5.0e5 | 1.1e3 | 2.5 pc | 27 pc | 3.0e4 to 8.4e4 | 965 |
| 47 Tuc-like (8e5, 6 pc, 9) | 2.0e6 | 1.9e5 | 0.39 pc | 51 pc | 1.1e5 to 5.5e5 | 5,848 |
| compact (5e5, 3 pc, 12) | 1.2e6 | 3.2e7 | 0.03 pc | 19 pc | 1.1e5 to 5.3e5 | 443 |
| sparse (5e4, 8 pc, 3) | 1.2e5 | 38 | 6.4 pc | 30 pc | 3.4e3 to 5.1e3 | 1,099 |

("to": the range between the cluster centred on a sector corner and in a sector's middle.) Central densities of 1e4 to 1e6 Msun/pc^3 [B] belong to the dense end only: a typical cluster is about 1e3. A galaxy's real density at the Sun is 0.099 systems/pc^3, so the 47 Tuc core is 4.9e6 times denser.

**Why this must be a density component.** Each cluster's tidal sphere (`r_t`) is a bounded region, so the sector gate (`predicted_star_count`) needs one extra term: the sum of the cluster densities at the sector, with a spatial index so a sector tests only clusters whose sphere meets it, not all 157 (or 12,000 in a giant elliptical). The skeleton's upper-bound density function (`_bound_raw_density_at`) must include the cluster peak or the sector gate rejects cluster sectors outside the disc. The population mix at a cluster sector is simply 100% "cluster" population.

**Tidal radius [C].** For an isothermal halo of 220 km/s, `r_t = R (M / (2 M_gal(<R)))^(1/3)` gives, for a 2e5 Msun cluster at 1, 4, 8, 20 and 50 kpc, 21, 52, 83, 153 and 281 pc (5e4 Msun: 13 to 177 pc; 1e6 Msun: 35 to 481 pc). The text's 30 to 100 pc is right for 8 kpc. Use the pericentre distance for an eccentric orbit. The King model's own `r_t` (19 to 51 pc above) is smaller than the Jacobi radius, so clusters underfill their tidal sphere. Rule: take `r_t` from the King solution and, if it exceeds the Jacobi radius at pericentre, lower `W0` until it does not.

**Core-collapsed clusters.** About 20% of Milky Way clusters [B]. Flag them in the catalogue; for generated clusters draw the flag with probability rising with `W0` and `t_rh` below the age. The King core is replaced by a power-law cusp `rho ~ r^(-alpha)`, alpha 0.7 to 1.3 [B] inside the radius where the cusp meets the King curve. Slice 2 can use `W0` = 12 as a stand-in; the cusp is a refinement, not a prerequisite. The compact exemplar above (`W0` 12) already puts 3.2e7 Msun/pc^3 at its centre.

## 4. Stars from the component

- **Radius by inverse transform.** Draw a mass fraction `X`, look it up in the King table (one `bisect` on a monotone array: O(log n), effectively constant), scale by `r0`, add a random direction. Plummer's closed form needs no table but has the infinite tail above; use it only truncated at `r_t`.
- **Per-sector draw.** The exact block-first draw of sampling-backfill-and-resume.md applies unchanged: the expected count in a sector is the integral of the cluster density over the cube, evaluated once per sector (a cube integral of the table, or a Monte Carlo of the table at plan time).
- **Age and mass.** Add a population `cluster` to `STELLAR_POPULATION_AGE_RANGES_GY` with the cluster's own age (10 to 13 Gyr); the redraw loop drops stars above the turn-off. White dwarfs and the giant branch come out of the existing evolution code. The horizontal branch is not modelled by the code and is a recorded gap (section 7).
- **Mean mass.** A Kroupa mass function between 0.1 and 0.85 Msun has a mean of 0.30 Msun [C]; the 0.4 Msun per system used above includes remnants and binaries [R]. The real present-day mass function is also depleted of low-mass stars by evaporation [R], so a count of `M / 0.4` is an estimate to within a factor of about 1.5.
- **Velocity.** Bulk cluster velocity plus a Maxwellian draw at the local dispersion `sigma(r)`. The 1 to 15 km/s range [B] is the central line-of-sight (1-D) value, which is what the Harris catalogue lists [R]; if the text meant the 3-D value, multiply by 1.7. The virial estimates [C] are 2 to 11 km/s for the four exemplars, so the range holds. Use the King model's own `sigma(r)` falling with radius, not a constant.
- **Bright first.** At 12 Gyr the visible stars below the turn-off are the large majority, and they are dim. The bright-first level (`GEN.44`) draws the giants, horizontal branch and blue stragglers first (hundreds to thousands per cluster) and leaves the dim main sequence to the lazy fill, exactly as for the galactic bulge.

## 5. Populations inside a cluster

**Metallicity prerequisite.** Add one number per star, [Fe/H], defaulting from the position (a mild disc gradient, GEN.134 item 4, about -0.06 dex per kpc [R]; clusters carry their own). This is not a change to every system: it feeds three things, the planet-occurrence factor (section 7), a small colour and luminosity offset, and the main-sequence lifetime (metal-poor stars live shorter). Do the first slice without the latter two.

**Two-peak sampling.** Boss's peaks are -1.5 and -0.5 [B]; real peaks move with galaxy mass and the widths are not given in the text (about 0.3 dex each [R]). Draw from a two-Gaussian mixture truncated at about [-2.4, 0.0]. Positions: draw the metal-rich clusters from the bulge-plus-thick-disc density the generator already has (as a probability density by rejection), and the metal-poor from a spherical power law near `r^-3.5` with a core of a few kpc [R], which is the halo the density model lacks. In the Milky Way table both come from the catalogue.

**Blue stragglers [R].** `N_BSS ~ M_core^0.4` (Knigge 2009: 0.38 ± 0.04), nearly independent of `Γ` at fixed core mass, with typical clusters holding tens to a few hundred. Place them as stars of 1 to 1.7 Msun on the main sequence (a merger product's mass) drawn within the core radius.

**Millisecond pulsars and X-ray binaries [R].** Counts scale with `Γ` (the number of dynamical encounters), which is why dense and core-collapsed clusters host most pulsars: 340 known in 45 clusters [S, compact-remnant-regions.md]. Relative to the 47 Tuc-like exemplar, `Γ ~ rho_c^2 r_c^3 / sigma` is 0.016 for the typical cluster, 16 for the compact one and 0.001 for the sparse one [C]; a spread of four orders of magnitude, so the rule is roughly "most pulsars come from the top few clusters". The text's 10 to 100 times over-representation per unit mass [B] agrees with the existing note's 10 times [S]; the two ends have different scopes (X-ray binaries are about 100 times, observed millisecond pulsars depend on selection [R]). Default: 10 times for generated clusters, scaled by `Γ` rank.

**Collisions [C].** A star in a 1e5 pc^-3 core meets another within its lifetime with probability about 0.1 (collision time 1.3e11 yr), 0.8 at 1e6 pc^-3 and 0.014 at 1e4 pc^-3 (focusing-dominated formula [R, Binney and Tremaine]). That is enough for the observed hundreds of blue stragglers in the densest clusters and negligible in the sparse ones.

## 6. Motion and update

- **Orbit.** A cluster is one body with a galactic velocity and an epoch, stored like a star system (orbital-updates.md). A typical halo orbit is 200 km/s; the cluster covers its own 5 pc half-mass radius in 24,000 years and a 4 pc sector in 19,600 years [C]. Inside the cluster, stars at 10 km/s drift 0.001 pc per century, and the crossing time is 0.3 to 4 million years [C].
- **Consequence.** Within a few hundred years the cluster's internal structure is frozen and the whole thing moves rigidly. Store each cluster star's offset from the cluster centre (like a moon from its planet), never an absolute position; its absolute position and sector address come from the centre at the query time. This is one level of the existing relative-to-primary hierarchy and avoids re-addressing 5e5 stars on every update.
- **Open design point.** The "every object stores a sector address" rule of the physics design assumes an object moves slowly in the sector grid. A cluster moves one sector in about 20,000 years, so the stored address of a cluster star is stale after that long. Default: a cluster star's sector is derived from the centre's sector path (the Hermite spline of orbital-updates.md), not stored. Boss to confirm.
- **Evaporation and tidal shocks.** Cluster mass loss over billions of years is in the text [B]; it is not needed for a single epoch. Store the mass at the generation epoch.

## 7. Planets in clusters

**Chemical barrier [R].** Giant planet occurrence rises with host metallicity roughly as `10^(2.0 [Fe/H])` (Fischer and Valenti 2005 form): 1 at solar, 0.1 at -0.5, 0.01 at -1.0, 0.001 at -1.5 and 3e-5 at -2.3 [C]. A search of 47 Tucanae's core for hot Jupiters (Gilliland et al. 2000 [R]) found none where field rates predicted several; the metal-poor suppression (47 Tuc is -0.7) is the usual explanation, with disc photoevaporation and encounters as alternatives. Rule: multiply the giant-planet probability by `10^(2 [Fe/H])`, capped at 1. Terrestrial planets are far less sensitive and stay near their normal rate for [Fe/H] above -1, falling below it [R].

**Dynamical barrier [C].** A planet is hard if its orbital speed exceeds the cluster's dispersion, soft otherwise. The boundary is `a_h = G M_host / sigma^2`: for a 0.5 Msun host, 440 AU at 1 km/s, 18 AU at 5, 4.4 at 10 and 2.0 at 15 km/s; the text's "about 10 AU soft, below 1 to 2 AU hard" corresponds to `sigma` 6 to 15 km/s. Soft orbits are ionised when a star passes within about `a`; the expected number of such passes in 12 Gyr is `n pi a^2 v t (1 + 2 G M / (a v^2))` with `v = sqrt(6) sigma`. With `a >= a_h` the survival probability is `exp(-that)`, which is near zero at once above `a_h` for any core density of 100 pc^-3 or more:

| Local density `n` (pc^-3) | `sigma` (km/s) | `a_h` (AU) | Planets beyond `a_h` survive 12 Gyr? |
|---|---|---|---|
| 4.8e5 (47 Tuc core) | 11.5 | 3.4 | no |
| 1e4 | 8 | 6.9 | no |
| 1e3 | 6 | 12 | no |
| 1e2 | 3 | 49 | no |

So the cull is a hard cut at `a_h`: drop every planet with `a > a_h(sigma(r))`, keep those inside (they are hardened, the rate of host collisions being about 0.1 per star at the very densest core [C]). Because `a_h` depends on `sigma` only, the cut applies from the cluster centre outward with `sigma(r)` from the King model; nothing needs a density fit. A planet in the outer part of a cluster, where `sigma` is 1 to 3 km/s, may still lie out to 50 to 440 AU. Boss's "survival decays exponentially in `a`, density and age" [B] is the soft-regime formula above and reduces to this cut.

**Free-floating planets [R].** Ionised planets are not lost from the generator's point of view: they become rogue planets in the cluster, which the interstellar-object rates already model; this note adds no rate.

**Pulsar planets, type B [C, S].** The 0.7% of millisecond pulsars with planets, 25% of them the captured-giant type [S, exotic note], gives 0.175% of cluster millisecond pulsars. PSR B1620-26 b is 1 of about 340 [C, 0.29%]. The two agree within the small-number scatter (one object). Rule: once clusters exist, enable type B, restricted to millisecond pulsars whose position lies inside a cluster. Host: a pulsar plus a white-dwarf companion with a circumbinary giant; the giant is allowed despite the metallicity cull because it formed around a different star and was captured.

**Not modelled by the code and recorded as gaps.** The horizontal branch (its colour is a metallicity and "second parameter" effect [B, R]); isochrone interpolation across metallicity (the code has analytic star models, not isochrones); the second-generation chemical populations of clusters [R]; neutron star and black hole retention after natal kicks, which is low in clusters [R] but matters for how many remnants a cluster keeps.

## 8. Clusters per galaxy

**Milky Way.** Load the catalogue as persistent rows. The authoritative inputs [B] are the Harris catalogue (1996, 2010 edition), Baumgardt and Hilker (2018, with later updates) for masses, structural parameters and dispersion profiles, and the Gaia-based velocities (Vasiliev and Baumgardt 2021). **Boss's input needed:** the repository has no catalogue file, and this environment could not download one. Supply the Harris catalogue table (public; about 157 rows with position, distance, `[Fe/H]`, `c`, `r_c`, `r_h`, `M_V`, `sigma_v`, and a core-collapse flag), or approve a one-time ingestion script that downloads it where the network allows. Catalogue positions are heliocentric (l, b, d); convert with the section 2.2 transform of galaxy-coordinate-system.md, and keep the real names.

**Generated galaxies.** `N_GC = S_N 10^(-0.4 (M_V + 15))` [B]; the generated galaxies already carry `M_V` (multiple-galaxies.md). Checks [C]: the Milky Way (`M_V` -20.9) has 157 clusters, so `S_N` = 0.69 (the text's "about 1" is round); an `M_V` -22 elliptical at `S_N` 5 has 3,155; a `M_V` -22.5 central giant at `S_N` 12 has 12,000. Draw `S_N` by morphology: spirals 0.5 to 1.5, ellipticals 2 to 8, central giants 8 to 15 [R]. The cluster system's mass is 6e-4 of the Milky Way's stellar mass [C] (3.1e7 of about 5e10 Msun); for ellipticals the fraction is several times higher [R], consistent with the larger `S_N`.

**Masses.** Lognormal with peak 2e5 Msun [B]; the width is missing from the text. Default: 0.5 dex (about 1.2 magnitudes in `M_V` [R]); truncated at 1e4 and 3e6 Msun. A 157-cluster draw with this width gives 1.8e8 systems and a median central sector of 1.1e5 systems, with the most massive cluster putting 2.8e6 systems in its central sector [C].

## 9. Slices

1. **Metallicity and the cluster table.** An [Fe/H] column on stars with the default gradient; the cluster table and the King tables; the Milky Way catalogue loaded once Boss supplies it. Tests: the King solver against the five textbook concentrations; the cluster's total expected count equals `M / 0.4` to within sampling error.
2. **Sampling and fill.** The cluster term in `relative_density` and in the skeleton's bound; the `cluster` population; bright-first level for cluster sectors. Tests: a cluster's sampled radial profile matches its table; no sector above the gate is skipped; the sector counts match the table above (to a few percent).
3. **Planets and pulsars.** The metallicity factor, the `a_h` cut, type-B pulsar planets in clusters, and the blue straggler count.
4. **Synthetic systems.** `S_N`, the mass function, the metallicity mixture, positions from the bulge and halo, for generated galaxies.

Not on the list: the collapse cusp, the horizontal branch and isochrones, evaporation, and the second chemical generation; each is a refinement with no effect on counts.

## Open questions for Boss (defaults written into the handoff)

1. **Catalogue file** for the Milky Way (section 8): supply the table, or approve a one-time download script. Default: wait for the file; build the synthetic generator first.
2. **Cluster stars' sector addresses** are derived from the cluster centre, not stored (section 6). Default: derived.
3. **Planet cut**: hard cut at `a_h(sigma)` (section 7) rather than the smooth exponential. Default: hard cut.
4. **Blue stragglers and pulsars follow two different laws**: `M_core^0.4` for blue stragglers, `Γ` for millisecond pulsars. Default: as in section 5.
5. **Mass-function width** 0.5 dex and `S_N` ranges by morphology are recalled, not in the text. Default: as in section 8.

## Corrections to the pasted text

- "Millions of stars within a 4 pc sector" is true only above about 3e6 Msun; 1e5 to 6e5 for typical and 47 Tuc-like clusters.
- The metallicity range stops at -0.5; the real range reaches about 0 [R].
- -2.3 is 0.5% of solar, not nearly pristine.
- The "Dynamic Filtering" step scales blue stragglers and millisecond pulsars "linearly with Γ, accounting for a sub-linear dependence in the most massive cores", while section 6 of the same text gives blue stragglers as `M_core^0.4 to 0.5`. Two laws are mixed: linear in `Γ` is the pulsar and X-ray-binary law; the sub-linear one belongs to blue stragglers.
- The tidal radius is given as 30 to 100 pc fixed; it follows from the cluster's mass and orbit (section 3).
- The 1 to 15 km/s dispersion is the 1-D central value, not the 3-D one [R].
- The turn-off mass is 0.93 Msun in the repository's lifetime law at 12 Gyr, not 0.85; the gap is metallicity (section 2).
- Section 7.2 says planets in hard orbits are hardened "counterintuitively"; that is Heggie's law for binaries, and it holds for a planet only while its orbital speed exceeds the perturber's.
- The cited "Flynn (2026), arXiv 2605.03099" and "Hurley and Shara 2001, Planet formation in globular clusters" could not be checked; titles and years should be confirmed before they are quoted [B].

## Evidence notes

[C] numbers all come from `calc.py` (King solver, Plummer sampling, tidal radius, hard-soft radius, survival, encounter rates, timescales, cluster counts, a 157-cluster Monte Carlo). The exemplar parameters (mass, half-mass radius, concentration) are round numbers chosen to span the real range, not a catalogue fit. [R] items to check once papers are reachable: Fischer and Valenti's exponent 2.0; Knigge's 0.38; Gilliland's null result; the mass-function width; the `S_N` ranges by morphology; the halo power-law index; the Harris catalogue's [Fe/H] ceiling; the 20% core-collapse fraction; the 100 times X-ray binary over-abundance (Clark 1975); the 0.3 dex widths of the metallicity peaks.

## Sources

Boss's research text (cited, not opened): Ashman and Zepf 1992; Baes and Camps 2015; Bahramian et al. 2013; Baumgardt and Makino 2003; Baumgardt et al. 2023; Bonatto and Bica 2010; Carretta et al. 2010; de Boer et al. 2019; Dias et al. 2016; Flynn 2026 (arXiv 2605.03099); Harris and van den Bergh 1981; Heggie and Hut 2003; Hurley and Shara 2001; Hurley et al. 2005; King 1966; McLaughlin 1999; Sigurdsson 1992; Trenti et al. 2010; Wolf et al. 2010. Recalled and used here: Fischer and Valenti 2005; Gilliland et al. 2000; Knigge 2009; Pooley 2003; Hui et al. 2010; Clark 1975; Binney and Tremaine, Galactic Dynamics (collision time); Spitzer 1987 (relaxation time). Repository: `galaxy/density.py`, `generation/star_population.py`, `physics/stellar_evolution.py`, `tuning.py` (population ages, `DEFAULT_SECTOR_EDGE_PC`), `physics/constants.py` (`LOCAL_STELLAR_DENSITY_LY3`).
