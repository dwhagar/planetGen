# A realistic galactic core: the nuclear cluster, the nuclear disc and a fill that can hold them

Boss (2026-10-16, 8:43 PM EDT), to Core 3: "I want a realistic galactic core, add that to Core 3 please. Modify our density and generation model so we keep our consistent rules and adjust sector parameters to control what is filled." At 8:48 PM EDT he added: "we may have to do a tiered system for mass and luminosity, where at the core both are 10% higher caps at the core, and then 10% lower in the halo and intergalactic space" (section 5.1). This reopens GEN.155 (now GEN.279), which [fill-order-curves-and-core.md](fill-order-curves-and-core.md) section 4.2 had closed as "leave the model alone". This note is the plan. It is a recommendation: nothing here is built, and nothing is built from it until Boss approves.

Informs: GEN.155 (reopened as GEN.279), GEN.235 (the rule set and profile), GEN.237 and GEN.244 (the component table and young clusters), GEN.223 (the orbital update), GEN.24 (the core fill), ADM.49 (the Gaussian core). Read with [galaxy-disk-density.md](galaxy-disk-density.md), [intergalactic-and-halo-contents.md](intergalactic-and-halo-contents.md) (the sector profile, section 5), [globular-clusters.md](globular-clusters.md) (the dense-component precedent), [mass-backfill.md](mass-backfill.md) (the per-sector level), [spherical-sector-coordinates.md](spherical-sector-coordinates.md) and [later-phenomena-scatter-rules.md](later-phenomena-scatter-rules.md).

Evidence tags as in the catalogue: **[S]** seen in the repository, **[C]** computed with the repository's code or a short script, **[R]** recalled and unconfirmed, **[I]** inferred.

## Summary

- **The real centre is two stellar systems and a gas ring on top of the bulge** [R, with the numbers of fill-order-curves section 4.2]: the nuclear star cluster (about 2.5e7 solar masses inside a few parsecs, Sgr A* at its centre), the nuclear stellar disc (about 1e9 solar masses inside 200 pc) and the central molecular zone (3 to 5e7 solar masses of gas). The generator has none of them: its centre is a bar bulge with a flat top, 880 to 1,060 systems per cell, against 1.1e4 to 3.4e4 in the real disc and about 1.5e7 in the cell holding Sgr A* [C].
- **Keep the rules, add terms.** The three structures enter the density model as three more terms of `relative_density`, three populations in the profile's `mix` (`nuclear cluster`, `nuclear disc`, and the existing `young` for the star-forming ring), and a gas term in `gas`. The count rule (Poisson of `lambda` times volume), the population rule and the scatter rules do not change. No special object, and no special case in `ensure_sector_generated` (section 3).
- **The one new rule is a ceiling on what a cell resolves.** A cell that expects more systems than the **cell star ceiling** (default 5,000, about 100 seconds at the measured 20 ms a system) fills only the stars at or above the **fill luminosity cap** (default 100 solar luminosities, raised for a cell that those stars would still fill past the ceiling, as the nuclear cluster's cell). The rest is *unresolved*: its count and mass are computed from the address, shown on the sector page and never stored. The per-sector state that records how deep a cell was taken (`bright_level_sol`) already exists and already supports topping a cell up to a deeper floor (section 4). Ordinary space never reaches the ceiling, so it is untouched.
- **The cost that decides it** [C]: the cell holding Sgr A* has about 1.5e7 systems; at the measured 16.4 ms a system that is 68 hours and, at the default 48 KiB a system, 720 GB for one cell. At the default cap it is 4.8e4 systems, 13 minutes. The nuclear disc inside 200 pc is 2.5e9 systems in full; a run that sets its ceiling low enough makes it 8e6, about what the current 200 pc core fill already is (7.5e6).
- **Planets nearly vanish in the nuclear cluster** [C]: with the survival rule of the cluster note, a planet outlives 12 Gyr only inside 0.03 AU at 0.5 pc from Sgr A*, 0.4 AU at 4 pc and 2 AU at 10 pc; in the nuclear disc the limit is 3 to 10 AU (section 7).
- **Density first, real-world, no fudging** (section 3.0): each term is fitted to a published number and checked in `mathcheck`. **Then tiered scatter limits** (section 5.1, five settings with sliders in 5.2): the luminosity floor (9,000) and the mass limit (14) the scatter draws its line at get a core pair (9,900 and 15.4 by default, 10% higher), the baseline in the disc and bulge, and a halo luminosity floor of 500. They touch no density. They move the pre-placed young stars by 30% to 50% for the mass limit and 11% to 14% for the floor [C], change nothing for old stars (so nothing in the halo at the defaults), and are a separate lever from the resolution cap.
- **For Boss** (defaults written in, section 10): the cap value; the tier boundaries; replacing the Gaussian core of ADM.49 by the new terms; and the two mass figures the model is tuned to.

## 1. What the real centre holds

Figures as in fill-order-curves-and-core.md section 4.2 [R], with the star counts rescaled here [C].

| Structure | Size | Mass | Systems per 64 pc^3 cell (0.4 solar masses per system) |
|---|---|---|---|
| Sgr A* | point | 4.3e6 | one object, already `add_galactic_nucleus` [S] |
| Nuclear star cluster | half-light radius about 4 pc | 2.0 to 2.5e7 | 1.5e7 in the 2 pc ball, 1.5e6 at 4 pc, 4.2e5 at 6 pc, 7.5e4 at 10 pc, 6e3 at 20 pc [C, Hernquist model, 3.1] |
| Nuclear stellar disc | radius 100 to 230 pc, scale height 30 to 45 pc | about 1e9 inside 200 pc | 3.4e4 at the centre, 1.1e4 at 100 pc, 3.7e3 at 200 pc [C, 3.2] |
| Central molecular zone | about 200 pc, a ring near 100 pc | gas 3 to 5e7 | none: gas |
| Young stars | central parsec disc (about 100 O and Wolf-Rayet stars, the S-stars), Arches (about 2e4) and Quintuplet (about 1e4) clusters | young mass about 1e4 to 2e4 each | components (section 6) |

The cell volumes are those of the spherical grid: shell 0 is a ball of 2 pc radius (33.5 pc^3) and every other cell is about one edge cube (64 pc^3) [S, spherical-sector-coordinates.md].

## 2. What the model gives today

- Density is flat through the middle: relative density 154 at the centre, 14 to 17 systems per pc^3, no cusp [S, fill-order-curves section 4.2]. ADM.49 added an optional Gaussian `core_amplitude` (default 0, scale `CORE_RADIUS_FRACTION` = 0.1 of the bulge scale radius, 158 pc at the default shape) [S, `density.py`, `tuning.py`]. It has no cluster, no disc and no cusp, and at 158 pc it is the wrong size for the cluster and a rough guess for the disc.
- The black hole is in: the potential has a softened Sgr A* (4.3e6, eps 1 pc) [S, galactic-potential.md], and the one core cell gets a quiescent black hole or, rarely, a quasar [S, `add_galactic_nucleus`].
- There is no per-cell star cap (`MAX_SYSTEMS_PER_SECTOR` only limits the estimate text) [S]. A cell is Poisson around its expected count and placed by Poisson-disk growth. That is the part that cannot meet the real centre.

## 3. The density terms

**Built (GEN.280, schema v99):** the nuclear cluster and disc are `_nuclear_cluster` and `_nuclear_disc` in `galaxy/density.py`, the `nuclear_cluster` and `nuclear_disc` populations of the profile, and `galaxy_shape.nuclear_amplitude` (`planetgen plan --nuclear-amplitude`) scales both; ADM.49's Gaussian core is removed. Departures from the text below, each so the numbers match the published ones in the default model: the cluster's mass is 6.4e-4 of the bulge's (the model's bulge holds 3.9e10 solar masses at 0.4 a system, not 1.5e10, so 1.7e-3 would put 6.5e7 in the cluster); the disc's total is 61 times the cluster's (1.5e9, 1.0e9 inside 200 pc); the cluster's cusp is held flat inside the radius where its density equals the core ball's mean (about 1.2 pc), so the ball's one cell holds exactly the mass inside 2 pc; and the age windows are single ranges with the observed mixes' means (cluster 3 to 12 Gyr, disc 4 to 12). `mathcheck` compares the mass, half-light radius, density at 1 pc, disc mass inside 200 pc and disc height with the published values. The gas term (3.4) and the young share of the disc are not part of this item.

All three are closed form, are ratios to the galaxy shape's own lengths and masses (so a toy galaxy gets a proportionate centre), and are stored as shape fields next to the bulge's. `population_densities` returns them as populations, so the skeleton bound, the Galaxy Map prisms and the profile read one source.

### 3.0 What the density must match

Boss (8:50 PM EDT): match the real world, or at least what is really possible, as closely as we can, with no fudging of density. Every term is therefore fitted to a published number, and the model's own `mathcheck` (the `reference` checks of `physics/mathcheck.py`, TEST.65) compares it with that number; a term that cannot meet its number within the stated tolerance is a bug in the term, not something a scatter setting patches.

| Observable | Observed | Model | Check |
|---|---|---|---|
| Bulge and bar mass | 1.4 to 1.7e10 of 5e10 solar masses of stars (Bland-Hawthorn and Gerhard 2016) [S] | 31% of the total | exists: `galaxy_bulge_to_total`, 0.31 +- 10% |
| Thick disc against thin | 6e9 against 3.5e10 [S] | 0.17 | exists: `galaxy_thick_to_thin_disk`, +- 30% |
| Sun's distance, local density | 8.2 kpc; about 0.1 stars per pc^3 within a few parsecs [R] | 0.099 systems per pc^3 at the Sun | `sun_galactic_radius` exists; a local-density check is a gap |
| Stellar halo | 1.1 to 1.9e9 solar masses; slopes 2.5 and 4.5, break 25 kpc [S, halo note section 4] | 1.5e9 | new with the halo (GEN.223 plan) |
| Total stars | about 1e11 [S] | about 3.4 times the old model's, close to 1e11 [S, galaxy-disk-density.md] | the plan estimate |
| Nuclear star cluster mass | 2.0 to 2.5e7 solar masses [R] | 2.5e7 | **new:** `galaxy_nuclear_cluster_mass`, +- 15% |
| Nuclear cluster size | half-light radius about 4.2 pc projected [R] | 3.8 to 4.2 pc | **new:** the projected half-light radius, +- 15% |
| Nuclear cluster density | about 1e5 to 5e5 solar masses per pc^3 at 1 pc [R, low confidence] | 2.9e5 | **new:** within a factor 3 |
| Nuclear disc mass | about 1e9 solar masses inside 200 pc [R] | 1.0e9 | **new:** `galaxy_nuclear_disc_mass_200pc`, +- 30% |
| Nuclear disc height | scale height 30 to 45 pc [R] | 35 pc | **new:** within the range |
| Central molecular zone gas | 3 to 5e7 solar masses [R] | 4e7 | **new:** component mass |
| Sgr A* | 4.3e6 solar masses, 8.2 to 8.3 kpc [S] | 4.3e6 | exists in the potential |

The values marked [R] are from memory of Schodel et al. 2014 and 2018, Launhardt et al. 2002, Sormani et al. 2022 and Bland-Hawthorn and Gerhard 2016, and must be confirmed against the papers before the checks go in (no network from this environment).

### 3.1 The nuclear star cluster

A Dehnen (1993) model with inner slope 1 (a Hernquist profile): `rho(r) = M a / (2 pi r (r + a)^3)`, `M` = 2.5e7, half-mass radius 5 pc so `a` = 2.07 pc (a projected half-light radius of 3.8 to 4.2 pc, against the observed 4.2 [R]). It has a cusp (`rho` is proportional to 1/r inside `a`), as the observed cluster does (slope 1.0 to 1.4 [R]). [C]:

| r (pc) | density (solar masses per pc^3) | mass inside | circular speed with Sgr A* (km/s) |
|---|---|---|---|
| 0.1 | 8.1e6 | 5.3e4 | 433 |
| 0.5 | 9.7e5 | 9.5e5 | 212 |
| 1 | 2.9e5 | 2.7e6 | 173 |
| 2 | 6.1e4 | 6.0e6 | 149 |
| 4 | 9.2e3 | 1.1e7 | 128 |
| 10 | 470 | 1.7e7 | 96 |
| 20 | 38 | 2.1e7 | 73 |

The 2 pc ball holds 6.0e6 solar masses, 1.5e7 systems. The sizes scale with the galaxy: `M` = 1.7e-3 of the bulge mass (2.5e7 of 1.5e10 [R]), `r_h` = 3.2e-3 of the bulge scale radius (5 of 1,580 pc). **The inner cusp below 0.1 pc is not resolved by any cell** (a cell is 2 to 4 pc); the model does not need it for generation, only for the orbit of anything placed there.

### 3.2 The nuclear stellar disc

An exponential disc with a sech^2 vertical profile, the form of the thin disc: `rho(R, z) = rho0 exp(-R / R_d) sech^2(z / 2 h)`, `R_d` = 90 pc, `h` = 35 pc, `rho0` = 215 solar masses per pc^3, which puts 1e9 inside 200 pc (65% of the total) [C]. It gives 215 at the centre, 71 at 100 pc and 23 at 200 pc, inside the 50 to 150 mean the research notes quote [R]. Scaling: `R_d` and `h` in proportion to the bulge scale radius, mass as the nuclear-cluster mass times 40 (1e9 / 2.5e7 [R]).

### 3.3 What replaces what

Remove ADM.49's Gaussian core (the field, its `--core-amplitude` option and the Gaussian in `density.py`): the nuclear terms are the core, with the real shapes. No compatibility shim; the migration is Core 1's (the shape table gains the nuclear fields). The bulge term stays: the nuclear terms add to it, and the profile's `region` label becomes `nuclear` where they dominate.

**Populations** (the profile's `mix`, section 5.3 of the halo note):

| Population | Age window (Gyr) | [Fe/H] | Dispersion |
|---|---|---|---|
| nuclear cluster | 80% 5 to 12, 20% 0.1 to 5 [R] | +0.1 to +0.3, broad [R] | 100 to 250 km/s, no net rotation to a first model |
| nuclear disc | 70% 8 to 12, 30% 0 to 8, star forming still [R] | +0.2 [R] | 70 to 90 km/s, rotates |
| young (existing) | 0 to 0.1 | as the disc's | 15 km/s |

Dispersion and metallicity are recalled and a first model, not a fit. The young share of the nuclear disc is tied to the gas term (3.4), not to the density.

### 3.4 The gas

The central molecular zone is a gas term in `gas`: about 4e7 solar masses in a ring of radius 100 pc and width 50 pc, thickness 30 pc [R]. It feeds the same `kind_factors` as the disc's gas (supernova remnants, molecular clouds, emission and reflection nebulae) and the young share of the nuclear disc (about 0.07 solar masses a year of star formation, so a young share of order 1e-3 of the disc's stars [R, I]). `nebula_field` must read this term, or the centre keeps its nebula rate of the plain disc; I did not check how it reads gas today (the halo note flagged the same question).

## 4. A cell that holds too much: the resolution cap

### 4.1 The cost [C]

| Cell | Systems | Time at 16.4 ms | Storage at 48 KiB |
|---|---|---|---|
| Sgr A* ball | 1.5e7 | 68 hours | 720 GB |
| Nuclear cluster, 4 pc | 1.5e6 | 6.8 hours | 72 GB |
| Nuclear cluster, 10 pc | 7.5e4 | 20 minutes | 3.5 GB |
| Nuclear disc, centre | 3.4e4 | 9 minutes | 1.6 GB |
| Ordinary disc | 6 | 0.1 second | 0.3 MB |

The 200 pc nuclear disc in full is 2.5e9 systems (1.3 years on one worker, 120 TB); the existing core fill of a 200 pc sphere is 7.5e6 systems [S, fill-order-curves section 4.3]. Nothing is generated unless asked (halo note 5.1), so the question is only how deep one asked cell goes. The rule must give the same answer from the same address, wherever it is.

### 4.2 The rule

**Built (GEN.281, schema v98):** steps 2 to 4 below are in `planetgen/galaxy/resolution.py` and the fill; the two settings are `galaxy_shape.cell_star_ceiling` and `fill_luminosity_cap_sol` (`planetgen plan --cell-star-ceiling/--fill-luminosity-cap`). A capped cell is taken to the cap by the existing level machinery; the dim systems are not drawn and the below-cut remnants are counted as unresolved. Step 5 is as described.

1. The profile gives `lambda` (expected systems per volume) and so the cell's expected count `N`.
2. **Ceiling.** The **cell star ceiling** is a setting of the request or run, default 5,000 systems a cell (about 100 seconds, 240 MB; Boss 1:46 PM EDT, cut from 50,000 after the nuclear cluster's cell took 14 minutes), the same for every cell in it. If `N` is at or below it the cell is filled whole, exactly as now.
3. **Fill luminosity cap.** Above the ceiling the cell fills only the stars at or above the **fill luminosity cap**, a second setting (default 100 solar luminosities, the lowest the code accepts today). Those are the pre-placed ones from the passes (at or above the scatter's floor for the band) plus, between the cap and the scatter's floor, a Poisson draw per population of `N` times its `bright_band_fraction` [S, `star_population.py`]. The draw is keyed by address and by the cap, so a cell taken down in two steps ends with the stars one step gives (the property mass-backfill.md already relies on). At the defaults the Sgr A* ball resolves 4.8e4 systems, a 4 pc cell of the cluster 4.7e3 and a 10 pc cell 240 [C]. This is a step and not a slope: a cell expecting 49,000 systems fills all of them and one expecting 51,000 fills about 160. If the step bothers him, the same two settings can feed a ladder that picks the lowest floor fitting the ceiling instead (an earlier draft of this note); the cell's recorded level is the same either way.
4. **Unresolved.** The stars below the floor are not made. Their expected count and mass follow from the same fractions (`N` times one minus the resolved share; the mass from the population's mean). The sector page shows "unresolved: about 1.5e7 systems, 6e6 solar masses, below 100 solar luminosities", computed from the address, never stored.
5. **State.** The floor reached is recorded in `sector_stats.bright_level_sol` (the level GEN.44 and GEN.187 already keep) [S]. Nothing new is stored: the existing "top up from the old floor to a new one" path is how a cell goes deeper.

Resolved counts for the two populations that matter [C, `bright_star_fraction`]: bulge-age stars at or above 100 / 300 / 1,000 solar luminosities are 3.2e-3, 1.7e-3 and 8e-5 of all systems; a young population at 100 / 1,000 / 9,000 is 1.5e-2, 3.2e-3 and 3.1e-4. So the ball resolves 4.8e4 systems at 100, 2.6e4 at 300 and 1.2e3 at 1,000; a 4 pc cell resolves 4.7e3 at 100. The code refuses a floor under 100 solar luminosities today (the brightest white dwarf [S, `_check_threshold`]), so 100 is the fill cap's lowest setting; with one cap no cell needs a lower rung. A cap under 100, or the ladder, would need white dwarfs counted in the fractions.

### 4.3 Why this and not the alternatives

- **Not a hard cap on the count.** Dropping random systems down to the cap changes the density the cell reports and breaks "same rules everywhere".
- **Not smaller cells.** A 1 pc edge in the core is a different grid; the edge is decided (4 pc, Boss) and the identifiers are 80 to 96 bits for it.
- **Not a statistical object.** A "nuclear cluster object" with a population but no cells is what was dropped in 2026-10; Boss now wants the rules kept. The unresolved part *is* a statistical population, but of the cell, derived from the profile like the dark matter on the sector page.
- **Reuses the level.** The floor is the existing per-sector depth; the cap is a rule choosing the first one.

### 4.4 Asking for more

A user, an administrator or a run can ask a cell for a lower floor; the work queue answers with the ETA and the storage (the existing progress and ETA path, halo note 5.5) and runs it as a job, the way the cluster core cells are queued. The cap is then a default and not a wall. A request for a floor below the cap's is the one place a human decision enters, and the page shows the cost before it runs.

### 4.5 Reading it back

- **Region queries** (stage C, DB.34 and GEN.276) return resolved objects like any others, with one addition: a result includes the unresolved count and mass of each cell it touched, so a count "within 50 pc of Sgr A*" is never silently short. Page and query limits stay as they are; a 50 pc sphere around the centre would be about 12,000 cells, each with at most 5e4 rows, so the query takes `LIMIT` and pagination as for any dense region [I].
- **The Galaxy Map** reads density from the model, not the rows (prisms), so the centre shows at once with no generation.
- **Runs.** `generate --core` (GEN.24) takes the cap as an option; its estimate (PERF.3) already samples cells and now reports the resolved and the unresolved counts separately, so the user sees what a cap buys before the run starts.

## 5. Populations and the profile at the core

The profile's fields (`lambda`, `mix`, per-population age and [Fe/H] and dispersion, `gas`, `kind_factors`) are filled for the new terms exactly as for the halo and clusters: add a row to the population table, no rule. Two details:

- **Pre-placed shares.** An old population holds no star at or above 9,000 solar luminosities or 14 solar masses (halo note 5.3); the nuclear cluster's old 80% adds none and its young 20% adds a few. The bright and mass passes therefore place almost nothing in the core ball by themselves; the fill and the cap floor produce its stars.
- **Compact remnants** follow `kind_factors` = 1 for every population (5.2e-3 neutron stars per living star [C, halo note]): the ball holds about 7.8e4 neutron stars and, at roughly 1e-3 per star, 1.5e4 black holes; the literature figure for the inner parsec is of order 1e4 black holes [R]. They are dark and cheap (no planets, no spectrum), so I recommend counting a remnant as one tenth of a system against the cap; the factor is a guess that the build measures, not a fact.

### 5.1 Tiered scatter limits

Boss (8:50 PM EDT): the density plan matches real-world observations as closely as we can, with no fudging; then, for the scatter sweeps, the minimum mass and minimum luminosity are tweaked by region. So the density (section 3) and the fill (section 4) are never touched by the tiers. They change only **where the scatter draws its line**.

**The limits that exist today** [S]. All three are global to the galaxy, and the same in every cell:

| Limit | Stored | Default and choices | Governs |
|---|---|---|---|
| Luminosity floor `L` | `galaxy_shape.bright_star_min_luminosity_sol` | 9,000 solar luminosities; presets from 2,500 to 4 million (`luminosity_floor.PRESETS`) | every star at least this bright is placed by the scatter |
| Stellar mass limit `M` | `galaxy_shape.bright_star_mass_limit_sol` (GEN.185), set by `--phenomenon-min-mass` (`tuning.PHENOMENON_MIN_MASS_SOLAR`) | 14 solar masses; presets 8 to 20 in steps of 2 (`PHENOMENON_MIN_MASS_PRESETS`) | every star born at or above it is placed by the scatter; also the cut of the neutron-star and black-hole below-cut draws (GEN.168) |
| Compact-object limit `C` | `--compact-min-mass` | the word "star" (use `M`), or 1, 2, 4 or 6 solar masses | the lowest neutron star or black hole the scatter places |

The backfill rings (1, 2, 5, 8 solar masses) lie below `M`, and the fill draws what is left (`FillContext.star_mass_limit_sol`).

**Proposal.** The profile gains a `tier` and a factor `t`: the core tier 1.10, the baseline 1.00, the halo tier 0.90. At a band, the limits are `L_band` = `t` x `L` and `M_band` = `t` x `M`, where `L` and `M` are whatever the galaxy was planned with (9,900 and 15.4 in the core and 8,100 and 12.6 in the halo at the defaults). `C` follows `M` when it is "star" and is left alone when the user set a number, since that is an explicit choice. The tier multiplies the user's chosen presets, and the product need not be a preset or lie inside the 2,500 to 4 million range: validation stays on the user's own value. The factors are tuning constants (`TIER_CORE_FACTOR`, `TIER_OUTER_FACTOR`), so no migration.

**Where the tier is chosen.** From the profile of the band (a shell and a band, the unit the scatter walks), at the band's centre and averaged over azimuth, by the dominant term of `lambda`:

| Tier | Factor | When |
|---|---|---|
| core | 1.10 | the nuclear terms (cluster and disc, 3.1 and 3.2) are at least half of `lambda` |
| baseline | 1.00 | everything else: thin and thick disc, bulge, arms, satellites and clusters (their own mix is a component, not a place) |
| halo | 0.90 | the halo term and the intergalactic floor are at least half of `lambda` |

**Discrete, not a blend.** The scatter's object-first sampler takes one certified bound per band and then accepts by true density, and the backfill and the fill must use the same limit for a cell. One value per band, a pure function of the address, keeps all three in step with no new state; a blend would give every cell its own limit inside a band and break the bound. Boundaries are where the dominant term changes, so they follow the density and not a radius; for the Milky Way the core tier is the inner 150 to 250 pc and the halo tier begins where the halo term passes the disc and bulge together, a few kiloparsecs above the plane [I]. A cell's tier is on its sector page.

**How the tiers sit with the ceiling and the cap.** The fill luminosity cap (section 4.2) works below the band's luminosity floor and sets how deep one dense cell is filled; the tier sets where the scatter leaves off and the fill starts. Between the cap and the band's floor a dense cell draws a Poisson number of stars per population. The two never overlap and never leave a gap.

**What it does to the scatter rows** [C, `massive_star_fraction`, `bright_star_fraction`]. The share of a population's living stars the scatter places:

| Population | mass limit 12.6 / 14 / 15.4 | luminosity floor 8,100 / 9,000 / 9,900 |
|---|---|---|
| young | 2.6e-4 / 1.8e-4 / 1.2e-4 | 3.6e-4 / 3.1e-4 / 2.8e-4 |
| intermediate | 0 / 0 / 0 | 4.7e-6 / 2.3e-6 / 2.2e-7 |
| old, bulge | 0 / 0 / 0 | 0 / 0 / 0 |

- **Fewer rows in the dense core, in proportion to its young stars.** 10% moves the young pre-placed share by +49% and -30% for the mass limit and +14% and -11% for the floor. In a cell of a young cluster such as the Arches (about 5e4 systems, all young [I]) the scatter places about 16 stars at the floor and 9 above the mass limit at baseline, and about 14 and 6 in the core tier. In the 2 pc ball around Sgr A*, where the young share is of order 0.25% [I], it is about 12 and 7 rows, becoming 10 and 5.
- **Nothing changes for old stars, which are almost all of the core and all of the halo.** At the defaults an old or bulge population places no star at 8,100 or at 12.6, so the halo tier adds no row in the halo, the intergalactic space or the old satellites. The brightest old halo star is about 2,400 solar luminosities, and about one halo star in 11,000 passes 1,000 [C, halo note]. **If Boss wants rows out in the halo, 10% will not do it**: it takes a halo tier that drops the floor by a factor of ten or more (a floor near 1,000 would place about 3e5 of the halo's 3.75e9 systems [C]). I would make the tier factors settings rather than fixed 1.10 and 0.90 so he can try that, and keep 1.10 and 0.90 as the defaults he asked for.
- **Totals** [I]. The galaxy's scatter rows change by the young share of each tier: a few percent of the bright and mass-pass rows at most, since almost all young massive stars sit in the disc and arms at baseline.

**What changes where**
- **Mass function and density:** nothing. A star's mass and age are drawn as before, and the density is the real-world model of section 3. The fill's share of a cell is `placed_star_fraction(L_band, M_band, population)`, so the scatter plus the fill give the same expected total as untiered (the halo note's rule 4).
- **Scatter rows:** the star passes read the band's `L_band` and `M_band` for their fractions and bands; the phenomena pass's below-cut draws (neutron stars and black holes under the mass cut) use the band's `M_band`. The registry needs no column (a kind does not carry the cut; the pass does); I did not check whether any builder hard-codes 14 [I], and the build greps for it.
- **Remnants:** a black hole above the mass cut is a scatter row, one below is drawn in the fill; with the tier the boundary is 15.4 in the core and 12.6 in the halo, and each side keeps its own draw, so no black hole is counted twice or lost. Neutron stars never reach 12.6 (they top out at 2.2) and are fill draws everywhere unless `C` is 1.0 or "star".
- **Metallicity (GEN.158, Core 2):** the fractions are per population and will gain a metallicity argument (a metal-poor star is hotter and brighter at the same mass [R], so more of them clear a luminosity floor). The tier and [Fe/H] are independent multipliers: the halo tier lowers the floor to 8,100, and brighter metal-poor stars mean a few halo giants can pass 8,100 that never passed 9,000. Core 2 should key its tables by population and [Fe/H] bin, and the cut value must be an argument and not a module constant.
- **Existing galaxies:** a galaxy scattered with global limits has the wrong rows in the core and halo tiers; `--redo-scatters` re-runs it. No shim.

**Tests.** (1) The tier of an address is a pure function and constant over a band and a cell. (2) With every factor 1 the scatter and fill are identical to today's (regression). (3) For sampled cells in each tier, scatter rows plus fill expectation equal the untiered expected total within the Poisson spread. (4) The mass pass places no star below `M_band` in any tier. (5) Cells either side of a boundary take the tier of their own band. (6) A cell scattered, then filled, then re-scattered gives the same rows. (7) The fraction functions are monotone in the cut and in [Fe/H]. (8) The tiers do not change `relative_density`, `lambda` or `mix` anywhere.

**Stored (DB.54).** The three settings are `galaxy_shape` columns: `core_min_luminosity_sol` and `core_mass_limit_sol` (NULL until set, then the baseline's times `CORE_TIER_FACTOR` 1.10, so 9,900 and 15.4 at the defaults) and `halo_min_luminosity_sol` (default 500, never under 50, so the range is 50 to 4,000,000); `planetgen plan --core-min-luminosity`, `--core-mass-limit` and `--halo-min-luminosity` set them, and `galaxy/scatter_limits.py` gives the values a tier reads. Which place takes which is GEN.333.

### 5.2 The sliders on the Generate page

Boss: each control a slider, clearly labelled, with a text box at the end of the slider that shows the number; use the existing layout and segments. The page already has the pattern: `generate.html`'s `mass_slider` (a `generate-slider` label holding an `<input type="range">` with a `<datalist>` of ticks and an `<output>` for the value), the luminosity floor `<select>` and the compact-limit `<fieldset>`, set out together by the `scatter_limits(id)` macro in a `search-fields scatter-limits` row on both the Plan and New galaxy forms [S]. The design reuses all of it:

- **One macro, `limit_slider(name, label, ...)`**, the old `mass_slider` generalised. Its row is: the label, the range input, then at the end a number box (`<input type="number">`, class `range-value`, same `name` as the form field) that shows the slider's value and accepts a typed one. `static/generatefolds.js`'s `data-range-output` hook keeps the two in step both ways; a typed value snaps to the slider's step and clamps to its range. The field hint underneath keeps today's wording style (what the number is, what a lower or higher one costs).
- **Luminosity sliders are logarithmic** (the floors run from 100 to 4 million): the range input holds a position and the box shows the luminosity, with ticks at 100, 300, 1,000, 3,000, 10,000 and so on. The existing preset ladder (`luminosity_floor.PRESETS`, 2,500 to 4 million) gives the baseline slider its stops; the halo, core and cap sliders step 1-2-5 from 100.
- **Mass sliders are linear**: the baseline keeps its presets (8 to 20 in steps of 2); the core mass limit runs 8 to 24 in steps of 0.2.
- **The baseline luminosity `<select>` becomes a slider of the same kind**, so all the scatter floors read alike. The compact limit keeps its checkbox and dropdown (a choice among four masses and "follow the star setting", not a continuous number).
- **Segments.** The existing row stays "Galaxy-wide" and is followed by three `<fieldset class="generate-mode">` groups under the same `scatter-limits` row, each with a legend and a one-line hint, and each folding with the rest (the `fold()` macro):

| Segment | Sliders (default) |
|---|---|
| Galaxy-wide (existing) | star mass limit (14), luminosity floor (9,000), neutron star and black hole limit |
| Core | core luminosity floor (9,900), core mass limit (15.4) |
| Halo and intergalactic space | halo luminosity floor (500, from 100) |
| Dense cells | cell star ceiling (5,000), fill luminosity cap (100, raised per cell) |

- **Linked defaults.** The core sliders start at 110% of the galaxy-wide ones and follow them until the user moves a core slider (then `data-touched` stops the following). A "Reset" link per segment returns it to the defaults. The page warns beside the core floor if it is set below the galaxy-wide floor (more rows in the densest region), and beside the ceiling with the PERF.3 estimate of the Sgr A* ball's cost.
- **Server side.** `generate_page.py` parses each field like `mass_limit_argv` does today (a number, within the range, else a form error naming the field) and passes it as a `planetgen plan` option; the settings are stored with the galaxy and shown on the Stats page beside the existing floor and limit.
- **Tests.** The macro renders the slider, the datalist and the box with matching ids; a typed value reaches the form and snaps; the defaults come from `tuning`; the server refuses out-of-range values and a core floor above 4 million; the options round-trip through the plan to `galaxy_shape` and back to the page; the existing mass slider and floor tests change to the new macro.

### 5.3 Horizontal-branch and AGB stars

Boss (9:12 PM EDT): why not build horizontal-branch stars into the mix? There is no good reason beyond the stellar evolution model not having the phases.

**What the model has** [S, `physics/stellar_evolution.py` `evolve_star`]: main sequence, then a subgiant that brightens to twice its main-sequence light, then one giant phase of 0.1 of the main-sequence lifetime (log-uniform luminosity 50 to 1,000 solar luminosities, 2% of it a tip of 1,000 to 2,500, 3,500 to 5,000 K), then a white dwarf at once. That giant phase is the red giant branch. The core helium-burning phase (the horizontal branch and the red clump) and the asymptotic giant branch are missing, so a star goes from the branch tip to a white dwarf. It was a simple chain built before the model had a metallicity or a halo, and the horizontal branch's colour depends on both [R].

**What the real phases are** [R, unconfirmed]:

| Phase | Luminosity (L_sun) | Temperature (K) | Lasts | Who has it |
|---|---|---|---|---|
| Red clump, red horizontal branch | 40 to 60 | 4,600 to 5,300 | about 0.1 Gyr | metal-rich stars, [Fe/H] above -0.7: the bulge and the disc |
| RR Lyrae (instability strip) | 40 to 60 | 5,900 to 7,400 | within the 0.1 Gyr | the intermediate halo, [Fe/H] about -1.5 |
| Blue horizontal branch | 30 to 60 | 7,500 to 20,000 | within the 0.1 Gyr | the metal-poor halo, [Fe/H] below -1.6 |
| Extreme horizontal branch (hot subdwarfs) | 10 to 40 | 20,000 to 35,000 | within the 0.1 Gyr | a few percent of old stars with heavy mass loss |
| Asymptotic giant branch | 100 to 2,500 (tip near 3,000 in the halo) | 3,500 to 5,000 | about 0.02 Gyr | all low-mass stars after the horizontal branch |

How blue the branch is depends on [Fe/H] first, with age, helium and mass loss as the "second parameters"; the share of blue against red falls from nearly all-blue at [Fe/H] -2.2 to nearly all-red above -0.5 [R]. The model would draw the colour from [Fe/H] (GEN.158, Core 2) and from the population's age; until GEN.158 lands it uses the population's default [Fe/H] (halo -1.6 inside 25 kpc and -2.2 outside, bulge and disc about -0.1 [S, halo note 5.3]).

**What it does to the scatter** [C for the base, I for the additions: AGB lasting 0.002 of the main-sequence lifetime from 100 to 2,500, log-uniform]. The halo rows at the floors Boss named:

| Floor (L_sun) | Rows now | Added by AGB | Added by horizontal branch | Rows after |
|---|---|---|---|---|
| 100 | 12 million | +2.6% | 0 | 12.2 million |
| 500 | 3.8 million | +4% | 0 | 4.0 million |
| 1,000 | 310,000 | +28% | 0 | 400,000 |

**The horizontal branch adds nothing at any floor of 100 or more**, because its stars sit at 10 to 60 solar luminosities. It shows up in the scatter only when the floor is under about 60, where it adds about 10% to the giants (0.1 Gyr against about 1 Gyr on the red giant branch) and, more to the point, makes the blue and red horizontal-branch types, RR Lyrae and hot subdwarfs exist in the halo, the bulge and the clusters at all. The sector fill draws them below any floor, so they appear in every filled cell and every cluster, which is where the real bright halo stars are. The asymptotic giant branch matters at the high floors (+28% at 1,000). The base count is itself an approximation: the model spreads giants log-uniformly in luminosity, while a real red-giant luminosity function falls roughly as L^-0.9 [R], which would roughly halve the rows at 500 and above; that is a separate refinement, noted here and not planned.

**Build** (a new item, Core 3, after the metallicity input exists or with its interim defaults):
1. `evolve_star` gains the phases after the giant phase for stars below the helium-flash mass (about 2.2 solar masses): horizontal branch (0.1 Gyr, luminosity 40 to 60 and a temperature from the colour draw) and AGB (0.02 Gyr, 100 to 2,500, tip by metallicity), then the white dwarf. The lifespan to the white dwarf lengthens by about 1.2% of the main-sequence lifetime. Constants in `tuning.py`.
2. The colour draw from [Fe/H] and age (a table of blue, RR Lyrae, red and extreme fractions), seeded from the star's own stream so a star's record does not change.
3. `star_population._bright_windows`, `_bright_table` and the mass-band tables mirror the new phases (they mirror `evolve_star`; the file says so), so the pre-placed fractions and the fill's complement stay consistent.
4. Naming and storage: the star keeps Yerkes class III (a horizontal-branch star is a giant by luminosity class) and the spectral letter its temperature gives, with "horizontal branch" or "asymptotic giant" in the type text; no column. Planets around them follow the existing giant-host rules.
5. Tests: the phase durations sum to the lifetime; the colour shares at [Fe/H] -2.2, -1.5 and -0.1; no horizontal-branch star at or above 100 solar luminosities; the halo counts at the three floors against the table above; the scatter and fill still sum to the untiered expected total.

## 6. Scatter and fill at the core

The rules of [later-phenomena-scatter-rules.md](later-phenomena-scatter-rules.md) apply unchanged. What the core changes:

| Phenomenon | At the core | Phase |
|---|---|---|
| Bright and mass passes | place few or none (old population); the central young disc and the Arches and Quintuplet clusters come in as young-cluster components (GEN.244's table, ages 2 to 6 Myr, Plummer radii 0.2 to 1 pc), which add their density and 100% `young` to the cell | scatter, as now |
| Neutron stars, black holes | kind factor 1; the cell's draw is the cap-limited one; a remnant in a resolved luminosity band is not special | fill |
| Pulsars | attribute by age; the observed shortage of pulsars near Sgr A* is a detection effect (scattering by the gas) that the generator does not model [R] | fill, no change |
| Magnetars (GEN.278) | the `young stars` gate follows the young density, so a few of the roughly 100 land in the centre by themselves; the real SGR J1745-2900 within 0.1 pc is not forced | scatter row |
| Supernova remnants, molecular clouds, nebulae | the gas term of 3.4 raises their rate across the central molecular zone; Sgr A East is not placed by name | scatter / fill, via `kind_factors` |
| Wolf-Rayet stars | the mass pass places them where young density is; the 100 or so O and Wolf-Rayet stars of the central parsec are the young component's | scatter |
| Rogue planets, comets | planet occurrence at the cell's [Fe/H] and the survival limit of section 7 | fill only |
| Multi-star companions | at scatter time as now; the cluster's soft binaries are cut by the same survival rule applied to the pair | scatter |
| Hypervelocity stars | none born here by the current rule; the centre is their natural origin and a later refinement | not built |

## 7. Planets at the core [C]

The cluster note's survival rule (expected ionising passes `n pi a^2 v t (1 + 2 G M / (a v^2))` over 12 Gyr, `v` = sqrt(6) times the dispersion, `M` = 0.5 solar masses) gives the largest orbit that survives:

| Place | systems per pc^3 | dispersion (km/s) | largest surviving orbit |
|---|---|---|---|
| Nuclear cluster, 0.5 pc | 2.4e6 | 212 | 0.03 AU |
| 1 pc | 7.1e5 | 173 | 0.06 AU |
| 4 pc | 2.3e4 | 128 | 0.4 AU |
| 10 pc | 1.2e3 | 96 | 2 AU |
| Nuclear disc, centre | 540 | 90 | 3 AU |
| 100 pc | 180 | 80 | 5.6 AU |
| 200 pc | 58 | 70 | 11 AU |

So the cluster inside a few parsecs holds essentially no planets (a hot orbit inside 0.03 AU is inside most stars' Roche limit), and the nuclear disc keeps close-in and terrestrial-zone planets of the dimmer stars only. This is the cluster note's slice 3 (the metallicity factor and the cull) applied with the dispersion and density of the new terms; the **cull rule should be the survival formula, not the hard cut at `a_h`** that note defaults to: at `sigma` = 70 to 90 km/s `a_h` = GM/sigma^2 is 0.1 AU while the formula keeps 3 to 11 AU. The hard cut was justified for densities above 100 pc^-3 at dispersions of 3 to 11 km/s, and the table reproduces its numbers (188, 11, 2.6 and 2.1 expected passes at its four rows [C]); the nuclear disc is outside that range. See section 10, item 4. Super-solar metallicity raises giant-planet occurrence by `10^(2 [Fe/H])`, 1.6 to 4 times, which the cull then mostly removes.

## 8. Motion

- **The potential needs the cluster and the disc.** The galactic potential has the softened black hole (eps 1 pc) and the bulge. The density model's bulge holds about 2e5 solar masses inside 20 pc [C, 154 times 0.0986 systems per pc^3 times 0.4], against the cluster's 2.1e7, so a star placed by this model would orbit a potential about a hundred times lighter than the matter there, apart from the black hole (the potential's own bulge term was not checked) [I]. The nuclear terms go into `galactic-potential.md`'s mass components so GEN.223's nightly update moves them correctly; at 10 pc the circular speed is 96 km/s and the period 0.6 Myr, at 1 pc about 35 kyr, at 0.1 pc about 1.4 kyr [C].
- **Rate of update.** GEN.223 moves every object and updates rows only when `next_update_due` elapses (decided). In the cap-limited core only the resolved objects exist, so the update cost is the resolved count, not the 1.5e7.
- **Unresolved stars have no orbit.** If a cell is later taken deeper, its added stars take their positions and velocities from the keyed draw at that epoch, as any late-added star would; the nightly update has not been applying to them before that.

## 9. Build items (proposed; numbers from TODO)

Order within my chain: after the free #1677 items, unless a step is wanted sooner. Core 1 owns every migration; steps 1 and 4 need one (shape fields), folded into whichever revision Core 1 names.

| Step | What | Depends on |
|---|---|---|
| 1 | Density terms: nuclear cluster, nuclear disc in `density.py`, `GalaxyShape` fields (migration), `population_densities` rows, ADM.49's Gaussian core removed, skeleton bound, prisms (**built, GEN.280, schema v99**) | GEN.235 profile |
| 2 | Ceiling and fill cap (4.2): the two settings (columns with the tier settings, one Core 1 migration), the cap applied to a cell above the ceiling, `bright_level_sol` use, unresolved count and mass on the sector page and in region-query results | step 1, stage C (GEN.276), 3b for the settings |
| 3 | Gas and young terms: the central molecular zone in `gas`, nebula rate, the nuclear disc's young share; young-cluster components for the central parsec, Arches and Quintuplet through GEN.244 | GEN.237, GEN.244 |
| 3b | Tiered limits (5.1) and the sliders (5.2): `tier` in the profile, the five settings (columns on `galaxy_shape`: Core 1 migration, `planetgen plan` options, Generate page sliders), the band's limits in the passes, the backfill and the fill, `placed_star_fraction` with the cell's cuts | GEN.235 profile; Core 2 (GEN.158) for the metallicity argument |
| 4 | Planet cull at the core: metallicity factor and the survival formula | cluster note slice 3, step 1 |
| 5 | Potential and orbits: nuclear mass components in the potential | GEN.223 |
| 7 | Horizontal-branch and AGB phases (5.3): `evolve_star`, the colour draw, the bright windows, tests | [Fe/H] (GEN.158, Core 2) or its interim defaults |
| 6 | Tests: the cell counts of sections 1 and 4; floor chosen by `N`; a cell taken down in two steps equals one taken down once; no cell above the cap resolved without a request; the survival table of section 7 | each step |

Tests for step 2 that matter most: for 30 random cells across the core the resolved count never exceeds the cap, and the unresolved count plus the resolved count equals the cell's expected count to within the Poisson spread.

## 10. Decisions and defaults

Boss's call is only needed where marked; the rest stand as written.

1. **Resolution cap: 5,000 systems a cell** (about 100 seconds, 240 MB; was 50,000, Boss 1:46 PM EDT). Higher fills more of the ring around the cluster at more cost. Default: 5,000; a galaxy setting (`--cell-star-ceiling`), to become an admin option.
2. **Replace ADM.49's Gaussian core by the nuclear terms** (no shim). Default: replace.
3. **Mass figures**: nuclear cluster 2.5e7 solar masses with half-mass radius 5 pc; nuclear disc 1e9 inside 200 pc; both scaled for other galaxies by the bulge. Default: as written.
4. **Planet cull in dense components uses the survival formula, not the hard cut at `a_h`** (section 7). Default: survival formula.
5. **Named objects (Sgr A East, SGR J1745-2900, the S-stars, the Arches and Quintuplet clusters by name)**: not placed by hand; the generator makes objects of those kinds where the profile puts them. Default: none by name.
6. **Remnants count as one tenth of a system against the cap**, to be measured. Default: one tenth.
7. **Tiers** (5.1, decided by Boss 9:09 PM EDT): a core pair of limits 10% above the baseline, the baseline, and a halo luminosity floor (default 500, slider from 100), chosen per band by the dominant term of `lambda`, with satellites and clusters on the baseline; the halo has no mass limit of its own. Default: as written.
8. **Ceiling and fill cap** (4.2, Boss 9:09 PM EDT): ceiling 50,000 (now 5,000), cap 100 solar luminosities, a step and not a slope (the cap rises per cell when the stars at the cap would still exceed the ceiling). Default: as written; the ladder is the alternative if the step matters.

## Evidence notes

- [C] figures: Hernquist cluster, exponential disc, cell masses, the survival table and circular speeds from short scripts over `star_population.bright_star_fraction`, `mass_band_fraction` and the formula of globular-clusters.md section 7; the per-system time and bytes are the cluster note's 16.4 ms and `DEFAULT_BYTES_PER_SYSTEM` = 48 KiB, not measured on the target server.
- 0.4 solar masses per system includes remnants and binaries [R, globular-clusters.md]; star counts are within a factor of about 1.5.
- Not checked: how `nebula_field` reads gas, the exact mass of the potential's bulge at 20 pc, and the real Hernquist fit (Dehnen slopes between 1.0 and 1.4 would raise the inner densities by up to a factor of 3 inside 0.5 pc).
- Observational figures (masses, radii, populations, black hole counts) are from memory of Schodel et al. 2014 and 2018, Feldmeier-Oelkers et al. 2014, Launhardt et al. 2002, Sormani et al. 2022, Genzel et al. 2010 and Hailey et al. 2018, unopened here.

## Sources

[fill-order-curves-and-core.md](fill-order-curves-and-core.md), [galaxy-disk-density.md](galaxy-disk-density.md), [globular-clusters.md](globular-clusters.md), [intergalactic-and-halo-contents.md](intergalactic-and-halo-contents.md), [mass-backfill.md](mass-backfill.md), [later-phenomena-scatter-rules.md](later-phenomena-scatter-rules.md), [galactic-potential.md](galactic-potential.md); `src/planetgen/galaxy/density.py`, `generation/star_population.py`, `generation/run_sector.py` (`add_galactic_nucleus`), `tuning.py`.
