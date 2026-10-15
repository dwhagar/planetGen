# Filling the whole sphere: the halo, the satellites and intergalactic space

Boss (2026-10-15): "investigate filling the entire spherical grid of sectors but when it's off the galactic plane or the proper galaxy, we'll have to adjust the sector contents, since clearly it'd be less ... what is known to be found in intergalactic space, stars, phenomena". This note says what is known to exist out there, what the generator does with that space today, and proposes a concrete way to give every cell of the sphere defined contents without storing the empty ones.

Informs: the whole-sphere build train (section 9), GEN.230 (the band field and the birth-field layout), GEN.232 (fast bodies), GEN.9 and VIEW.2 (neighbour galaxies), GEN.78 (the density floor), GEN.100 and GEN.223 (scatter and the nightly motion of every object), GEN.134 (metallicity), GEN.42 (sampling), the four slices of [globular-clusters.md](globular-clusters.md), which this note folds into one build, DB.32 follow-ups (object ID layout)
Status: research, 2026-10-15; nothing here is built. Boss approved the reach and 96-bit IDs and set the rule (one rule set, per-sector parameters, nothing generated unless asked, clusters in the same train) the same day; sections 5 and 9 are written to that rule and are the build plan for the TODO thread. Section 8 lists what is still open.

Evidence tags: **[S]** seen: in this repository, or in a search-result summary read during this research (sources at the end; the papers themselves were not opened, the environment cannot reach arXiv or journals); **[C]** computed here with `/mnt/project-files/research/scripts/halo/halo_numbers.py` (output beside it in `halo_numbers.txt`; numpy and the repository's own grid code only); **[R]** recalled and unconfirmed; **[I]** inferred from tagged facts, with the step shown. A figure keeps the tag of its weakest input.

Related notes: [spherical-sector-coordinates.md](spherical-sector-coordinates.md) and [galaxy-coordinate-system.md](galaxy-coordinate-system.md) (the grid), [galaxy-disk-density.md](galaxy-disk-density.md) (the density model), [interstellar-object-rates.md](interstellar-object-rates.md) and [compact-remnant-regions.md](compact-remnant-regions.md) (kinds and rates), [phenomenon-scatter-mass-cut.md](phenomenon-scatter-mass-cut.md), [sampling-backfill-and-resume.md](sampling-backfill-and-resume.md) (cost of sparse sectors), [globular-clusters.md](globular-clusters.md), [multiple-galaxies.md](multiple-galaxies.md), [galactic-potential.md](galactic-potential.md), [object-ids.md](object-ids.md).

## Summary

**What is out there.** Almost nothing, by the generator's own yardstick. A sector that "exists" today expects at least one system. The stellar halo is about 1.5e9 solar masses [S], 3% of the disk and bulge, and spread so thinly that a cell at the Sun's distance from the centre expects about 0.01 systems [C]; at 20 kpc, 4e-4; at 50 kpc a cell expects 3e-5 [C]; at the 262 kpc edge of the Milky Way's neighbourhood, 2e-8 [C]. The intergalactic average is 1e-7 systems per cell or less [C, R]. Only discrete things reach one system per cell: the Magellanic Clouds, the Sagittarius remnant, dwarf spheroidals and globular clusters. Everything else is a sparse Poisson haze of old, metal-poor stars, ejected stars and ejected neutron stars, with a gas background of fractions of an Earth mass per cell.

**What the project does today** [S, C]:

1. **The grid is already a tiling of all space, but its addresses stop short.** The 80-bit object ID cannot name a cell whose band number passes 2,047, which means every shell from 1,304 (5.2 kpc) outward is only partly addressable: above latitude 31 degrees at the disk's edge (shell 3,764), no cell has an ID. Shells end at 4,095 (16.4 kpc). The real halo, the Magellanic Clouds (50 kpc) and the Milky Way's virial radius (about 200 kpc) are outside.
2. **The "halo" is a flat floor, and it is wrong in both directions.** `MIN_RELATIVE_DENSITY` is 0.1% of the local density everywhere (0.0063 systems per cell). It is half the real halo at the Sun, 200 times too high at 50 kpc, and 340,000 times too high at 262 kpc [C]. If the grid were simply filled out to 262 kpc under the floor, it would hold 7.4e12 phantom systems, 25 times the real galaxy's 3.0e11 [C].
3. **Nothing beyond the disk outline can be generated**, the hypervelocity stars are squeezed inside it (they are redrawn when they land outside), and the bright-star floor (9,000 L_sun) finds no halo star at all [C].

**Recommendation** (section 5, decisions in section 8):

- **One rule set, a profile for every cell.** Every cell in the reach has a profile computed from its address alone: expected systems, a mix of stellar populations (age, metallicity, motion), gas, dark matter, and a factor for each rare kind. The same counting, population, scatter and fill rules run everywhere; only the profile differs. Nothing is generated unless asked; an asked sector is generated and saved wherever it is, even if empty (as GEN.76 already does). The qualifying outline stays as an optimisation (map shading, random starts), not the definition of a sector.
- **Reach 65,535 shells = 262 kpc** (1.18e15 cells), with 96-bit object IDs (approved). It contains the Magellanic Clouds, the Sagittarius remnant, the globular clusters, the satellites out to Leo I and the Milky Way's virial radius [S]. Andromeda and the rest of the Local Group stay sky objects (GEN.9 stage 1).
- **Replace the flat floor with a measured stellar halo**: a broken power law, 1.5e9 solar masses, 3% of the visible mass, old and metal-poor. It lowers the cell density out there by up to five orders of magnitude and, by itself, makes no sector qualify.
- **Four kinds of content beyond the disk**, each cheap: the halo itself (density term), discrete components (globular clusters and satellites in one table, built in the same train), wanderers (hypervelocity stars, kicked neutron stars) with analytic radial laws, and computed gas and dark-matter figures on the sector page instead of rows.
- **Uncharted sectors stay implicit.** Storing every cell of the 16 kpc sphere would take 365 years of single-worker time at the measured fixed cost per sector and 29 TB for bare rows [C]; the 262 kpc sphere is out of the question. A cell gets rows only when someone asks for it.
- **The scatters fill the sphere too**, with the same algorithm: the halo, old clusters and old satellites place no bright or massive star at the default floor and cut (a 10 to 13 Gyr population has none) [C]; the Magellanic Clouds place a few thousand; about 2.4 million extra black-hole rows appear outside the disk, 4% on top of the disk's scatter [C]. The one part that cannot reach 262 kpc as built is the phenomena pass, a band-by-band walk; it moves to the object-first sampler the star passes already use.

## 1. What the grid reaches today [S, C]

The shells tile all of space (plan, section 2: "No undefined space"). What limits them is the address, not the geometry.

| Radius | Shell S | Cells inside (sphere) | Object ID layout shell/band/wedge bits | ID bits |
|---|---|---|---|---|
| 16.4 kpc | 4,095 | 2.88e11 | 27/14/15 | 96 |
| 25 kpc | 6,250 | 1.02e12 | 25/15/16 | 96 |
| 50 kpc | 12,500 | 8.18e12 | 23/16/17 | 96 |
| 100 kpc | 25,000 | 6.5e13 | 21/17/18 | 96 |
| 262 kpc | 65,535 | 1.18e15 | 19/18/19 | 96 |
| 524 kpc | 131,071 | 9.4e15 | 17/19/20 | 96 |
| 1.05 Mpc | 262,143 | 7.5e16 | 47/20/21 | 128 |

("Layout" is what `object_uid.fit_layout` returns for the full sphere to that shell: all bands, all wedges. The sphere to shell 3,764, where the disk outline ends, holds 2.23e11 cells; the outline's 2.39e10 qualifying cells are 10.7% of them and 2e-5 of the 262 kpc sphere.)

- **The 80-bit layout was sized for the outline, not the sphere.** Its band field is 12 biased bits (bands -2048 to 2047). A shell's top band is about 1.57 S, so from shell 1,304 on the polar bands (and the caps) have no address. At shell 2,500 (10 kpc) every cell above latitude 47 degrees is unaddressable; at the outline's edge, above 31 degrees [C]. This is the bug GEN.230 already files, and the layout decision it asks for is answered here: a full sphere needs the 96-bit fit, not a 41st bit. The default galaxy never meets it because its outline reaches only band 1,456; a halo, which is not confined to the plane, would meet it at once. `provisional_sector_designation` and `parse_sector_designation` also hard-code 12/12/16 bits while `fit_layout` already widens the object ID per galaxy [S, `spherical.py`, `object_uid.py`], so the designation functions need the layout.
- **96 bits cover every sphere up to 524 kpc** and cost 2 more bytes per object ID than today [C]. Beyond that, 128 bits. The layout is chosen once per galaxy by the plan step (`outline_object_id_layout`), so a galaxy that wants no halo keeps 80 bits.
- **The default galaxy's outline** (`skeleton.build_shell_extents`) stops at the first shell whose equatorial band cannot reach one expected system. It is built from the disk and bulge terms only; no halo term changes where it ends (the halo at 15 kpc is 0.26% of the 0.159 relative density the edge needs) [C].

## 2. What is known to be out there

### 2.1 The stellar halo

| Property | Value | Tag |
|---|---|---|
| Stellar mass | 1.5 +- 0.4e9 solar masses (Deason et al. 2019 in a search summary); about 3% of the 5e10 total [R for the total, BHG16] | S, R |
| Density profile | Broken power law. SDSS fit: inner slope about 2.3, outer about 4.6, break at 27 kpc, flattening q = 0.6; others find breaks from 10 to 30 kpc, a recent Subaru fit 17 kpc with slopes 3.3 and 4.8; outer slopes agree at 4 to 5, inner slopes and breaks do not | S |
| Shape | Oblate inside the break (q about 0.6), rounder outside [R] | S, R |
| Metallicity | Inner halo peaks near [Fe/H] = -1.6, outer halo near -2.2; mean about -1.2 to -1.5 | S |
| Age | Mostly 10 to 13 Gyr [R]; the dominant accreted dwarf (Gaia-Enceladus) arrived 8 to 11 Gyr ago and makes at least half the metal-poor halo within 30 kpc | S |
| Motion | No net rotation, dispersion about 120 km/s | S (`galactic-potential.md`) |
| Local density | About 0.2% to 0.5% of the local stellar density: 0.185% from the mass and profile above [C], about 0.5% recalled (Juric et al. 2008) [R] | C, R |

The mass and the local density disagree by a factor of 2.7 and I could not settle which is better; section 7 shows it does not change any conclusion. The repository's flat floor is 0.10% of the local density.

**Stars you would actually find** [C]: sampling 60,000 living stars of age 10 to 13 Gyr with the repository's own stellar model gives a brightest star of 2,400 L_sun and 0 above 9,000 L_sun, and the repository's population tables give the same zero for the bright and the mass pass (8.9e-5 of halo stars, one in 11,000, pass 1,000 L_sun). The default bright-star scatter (9,000 L_sun) therefore places no halo star, and the dim halo is entirely the sector fill's job. The repository does not model horizontal-branch stars, which dominate halo star counts at the bright end in reality [S, `globular-clusters.md`].

### 2.2 Satellites, clusters and streams

| Thing | What is known | Tag |
|---|---|---|
| Large Magellanic Cloud | 49.6 kpc; stars 2.3 to 2.7e9 solar masses (2% of the default galaxy's systems); total mass 1 to 2e11 | S |
| Small Magellanic Cloud | 62.4 kpc; stars 5.3e8 (0.4%) | S |
| Sagittarius remnant | 22 to 25 kpc; stars 5e7 to 1.5e8 including its tails; its stream wraps the galaxy | S |
| Other dwarf satellites | About 60 known [R]; Fornax at 143 kpc [S]; most are far below 1e7 solar masses [R] | S, R |
| Globular clusters | 157 in the Harris catalogue [S]; the halo ones reach past 100 kpc (AM 1 near 120 kpc) [R]; 1.8e8 systems for a full set (0.06% of the galaxy) [S, `globular-clusters.md`] | S, R |
| Stellar streams | Tens of known streams [R]; a Sagittarius-sized stream (1e8 solar masses in a tube of 1.5 kpc radius and 150 kpc length, all [R]) is 50 times denser than the halo at 30 kpc and still expects only 0.015 systems per cell [C] | R, C |

Together the satellites and clusters are about 2.7% of the default galaxy's systems, comparable to the whole stellar halo (1.25%) [C]. They are the only far content dense enough to qualify: for example the LMC's stars, about 6e9 systems, make it a galaxy of its own inside the grid.

### 2.3 Stars and remnants that leave the galaxy

- **Hypervelocity stars** (ejected by the central black hole). Ejection rate: at most 1e-5 per year at 2 sigma, from a 2024 Gaia search that found none new and predicts 5 to 45 unbound ones in the whole catalogue; 0.7e-4 per year for a top-heavy centre; about a dozen candidates known [S]. A steady stream at speed v puts `rate x r/v` stars inside radius r, so the density falls as `r^-2` (the law `tuning` already uses): at 500 to 1,000 km/s that is 160 to 2,200 inside 16.4 kpc and 2,600 to 36,000 inside 262 kpc [C]. The repository draws 1e3 to 1e4 for the whole galaxy, consistent with the total out to about 100 kpc [S, `interstellar-object-rates.md`].
- **Neutron stars kicked out.** About 40% of neutron stars and 2% of black holes reach escape speed (one population synthesis, in `compact-remnant-regions.md`) [S]. From the usual kick distributions against the repository's 557 km/s escape speed I get 22% for either of two published kick laws [C on R inputs], so 20% to 40%. On 1e9 isolated neutron stars that is 2e8 to 7e8 unbound. Flying at 300 to 500 km/s over 10 Gyr they spread over 3 to 5 Mpc; about 1e7 are inside 262 kpc [C, assuming steady production]. Per cell that is 1e-8.
- **Hostless supernovae and runaway stars** are real but add nothing a cell-level rule needs [S, `interstellar-object-rates.md`].

### 2.4 Gas

| Thing | Value | Tag |
|---|---|---|
| High-velocity clouds | HI mass about 3e8 solar masses, ionised mass at least as much; most at 2 to 15 kpc; the Magellanic Stream at about 55 kpc, 1,500 square degrees on the sky, possibly out to 100 to 150 kpc; the number of discrete clouds is a few hundred [R] | S, R |
| Hot circumgalactic gas | About 6.5e9 solar masses within the virial radius from one fit, 5% to 15% of the missing baryons; estimates vary by factors of 3 to 4 | S |
| Hot halo at about 50 kpc | n of order 1e-4 per cm3 [R]: 2.2e-4 solar masses, 74 Earth masses, per cell [C] | R, C |
| Warm-hot intergalactic medium | 1e-6 to 1e-4 per cm3; the Local Group is 20 to 30 times the cosmic mean | S |
| Cosmic mean baryons | 1.9e-7 H per cm3, 6.2e-9 solar masses per pc3: 0.13 Earth masses of gas per cell [C]; the Local Group's 20 to 30 times that is 2.6 to 4.0 Earth masses | C |

Gas at these densities is background, like the diffuse classes A and B the generator already leaves ungenerated. High-velocity clouds are the one gas item worth a row (section 5.8).

### 2.5 Intergalactic stars

- **Cluster-scale intracluster light.** Virgo: 15.8 +- 3 +- 5% of the starlight from planetary-nebula counts, about 7% in a later reanalysis of the core; simulations give 10% to 50% rising with cluster mass; rich clusters of 20% or more are reported [S]. This is a galaxy-cluster effect, far from the Milky Way's group.
- **The Local Group** has no intragroup-light measurement that I found. What sits between the Milky Way and Andromeda (about 770 kpc [S, `multiple-galaxies.md`]) is tidal debris, ejected stars and neutron stars, and gas [I].
- **The cosmic average.** About 5.6e8 solar masses of stars per Mpc3 [R] is 9e-8 systems per 64 pc3 cell, and a tenth of that in a void [C, R]. The disk outline needs 1.1e7 times more.

### 2.6 Dark matter

It matters here as gravity only. The repository's NFW halo (`galactic-potential.md`) gives 0.0091 solar masses per pc3 at the Sun (0.345 GeV per cm3), 2.0e-4 at 50 kpc, 3.2e-5 at 100 kpc and 2.1e-6 at 262 kpc, which is 1.4e-4 solar masses per cell; 1.0e12 solar masses lie inside 262 kpc [C on S]. No object is generated from it, and the orbit code already feels it.

### 2.7 What I could not find

Nothing in these searches gave a rate of free-floating planets or brown dwarfs in the halo, a count of discrete high-velocity clouds, an intragroup-light fraction for the Local Group, or the hot gas density at a given radius from a direct measurement. The proposal does not need them; section 7 says where they would matter.

## 3. What the generator does with far space today [S]

- **Density.** `relative_density` is the disk, bulge and thick disk plus the flat floor `max(model, 1e-3)` (`galaxy/density.py`). The floor is counted as "old" population.
- **Outline.** `build_shell_extents` walks outward from the core and stops at the first shell whose equatorial band cannot qualify (expected systems at least one); `galaxy_shell` stores one row per shell with its highest band. `GalaxyBounds.contains` refuses any address outside, and `random_address` draws uniformly over the outline's cells.
- **Scatter.** The bright-star and mass passes are object-first (`object_first.py`): one Poisson count per shell from a certified bound at about 1.2 microseconds per band [S]. The phenomena pass (`phenomenon_scatter._band_rows`) still walks every band with 32 longitude bins. Every band to 262 kpc is 6.75e9 bands (the sum over shells of pi times the shell number) [C]: 2.25 hours a star pass, and too many for the phenomena walk.
- **Asking.** `ensure_sector_generated` already generates and saves a requested address even when its draw is empty (GEN.76), and an address outside the outline only with `outside_ok`, at the floor's density (GEN.81).
- **Hypervelocity stars** are placed along a random ray to the outline's reach and redrawn when they land in a band the outline excludes, so they stay near the plane and inside 15 kpc (`hypervelocity_rows`).
- **Rates** are per star (`phenomenon_rate_per_star`), so they follow whatever star count a cell has; nothing changes for a cell with few stars.
- **Cost of a sparse sector**: about 40 ms fixed plus the systems (70 to 90 ms per system for a sector of one system); a per-cell Poisson draw without the database is 0.06 ms [S, `sampling-backfill-and-resume.md`].

## 4. The halo density model, with numbers

**Form** [I, parameters from section 2.1]: a stellar halo term added to `_raw_density`, a spheroid of axis ratio `q`,

```
m          = sqrt(R^2 + (z / q)^2)               # R cylindrical, z height
halo(m)    = A * (max(m, r_core) / r_break)^-n   # n = n_in for m < r_break, n_out beyond
```

| Parameter | Baseline | Range in the sources | Tag |
|---|---|---|---|
| Mass (A is solved from it) | 1.5e9 solar masses = 3.75e9 systems at 0.4 solar masses a system | 1.1 to 1.9e9 | S |
| Inner slope | 2.5 | 2.3 to 3.3 | S |
| Outer slope | 4.5 | 4.6 to 4.8 | S |
| Break radius | 25 kpc | 17 to 27 kpc | S |
| Axis ratio q | 0.7 | 0.6 to 0.8 | S, I |
| Core | flat inside 0.5 kpc (the bulge dominates there) | | I |

As with the other terms, the radii are ratios to the shape's disk scale length (the break is 9.6 scale lengths) and the mass is a fraction of the visible mass (3%), so a galaxy of any size gets a halo that scales with it. The amplitude is solved like `k_norm`, never typed in.

**Result** [C], in the plane unless marked, per 64 pc3 cell:

| Place | Relative density | Systems per cell | Flat floor now | Floor / halo | Typical star spacing |
|---|---|---|---|---|---|
| Sun | 1.8e-3 | 1.2e-2 | 6.3e-3 | 0.5 | 18 pc |
| 8.2 kpc, 5 kpc up | 9.1e-4 | 5.8e-3 | 6.3e-3 | 1.1 | |
| 10 kpc above the centre | 4.6e-4 | 2.9e-3 | 6.3e-3 | 2.2 | |
| 20 kpc above the centre | 6.2e-5 | 3.9e-4 | 6.3e-3 | 16 | |
| 30 kpc | 5.0e-5 | 3.2e-4 | 6.3e-3 | 20 | |
| 50 kpc (LMC distance) | 5.0e-6 | 3.2e-5 | 6.3e-3 | 199 | 126 pc |
| 100 kpc | 2.2e-7 | 1.4e-6 | 6.3e-3 | 4,500 | 357 pc |
| 200 kpc | 9.8e-9 | 6.2e-8 | 6.3e-3 | 1.0e5 | 1.0 kpc |
| 262 kpc | 2.9e-9 | 1.8e-8 | 6.3e-3 | 3.4e5 | 1.5 kpc |

(Spacing is `n^(-1/3)` with n the system density. The Sun's neighbours in the disk are 2.2 pc apart.)

- **No cell qualifies.** The largest halo expectation anywhere beyond 8.7 kpc from the core is 0.01 systems per cell; the outline is unchanged [C]. The three extreme parameter sets of section 2.1 move the 50 kpc value between 3.4e-6 and 4.5e-5 and the 100 kpc value between 1.2e-7 and 1.9e-6 [C]; the conclusion that a halo cell is almost always empty holds in all of them.
- **Where halo cells are as dense as a threshold**: expectation of 0.01 out to 8.7 kpc (3.0e10 cells), 1e-3 out to 22 kpc (4.8e11 cells), 1e-4 out to 39 kpc (2.7e12), 1e-5 out to 65 kpc (1.2e13), 1e-6 out to 108 kpc (5.7e13) [C]. A lower qualifying threshold therefore buys a few more shells at an enormous cell count each; it is not a way to "fill" the halo.
- **Halo systems outside radius r**: 2.7e9 beyond 5 kpc, 1.6e9 beyond 16.4 kpc, 1.0e9 beyond 25 kpc, 3.6e8 beyond 50 kpc, 1.3e8 beyond 100 kpc, 3.0e7 beyond 262 kpc [C]. Inside 262 kpc the sphere has 1.18e15 cells for 3.7e9 halo systems: **one cell in 317,000 is occupied** [C].
- **The star spacing matters for navigation.** Beyond about 30 kpc the typical gap between halo stars passes 50 pc, the hop limit Boss set for the route planner (NAV.43 to NAV.45), so a halo star can be reached only from a satellite or cluster, or by a longer hop [I].

## 5. The proposal: one rule set, a profile for every cell

Boss's rule (2026-10-15): the program has the same rules everywhere in the sphere, and each uncharted sector carries parameters that tune those rules to where the sector is and what is expected there. Nothing is generated anywhere unless asked, as now; a sector that is asked for is generated, by the same rules, wherever it is. The first neighbourhood pull stays in the galactic region.

### 5.1 Decided and carried over

- **Reach 262 kpc, 96-bit object IDs: approved (Boss, 2026-10-15).** 65,535 shells, 1.18e15 cells. The reach is a plan-time setting stored with the galaxy; a cell beyond it is refused ("outside the grid"), as cells are today. At 500 km/s a wanderer needs 510 Myr to cross it [C], so none leaves within any plausible game clock. A galaxy with the halo is re-planned, or renumbered by a migration like the shell move's (no backward compatibility); the 80-bit layout stays for a galaxy planned without a reach.
- **Asked means generated and stored.** Today `ensure_sector_generated` already saves a requested address even when its draw comes out empty (GEN.76) and, with `outside_ok`, an address outside the outline at the halo floor's density (GEN.81). That behaviour becomes the rule for the whole reach: the outline check and `outside_ok` go, the reach is the only refusal, and every asked sector is saved and marked generated, empty or not. An uncharted cell is never stored; its profile is computed from the address when someone looks.
- **The neighbourhood pull stays in the galactic region.** The random start and the default neighbourhood still draw from the qualifying outline (`random_address`); an explicit address, span or block anywhere in the reach is accepted. A neighbourhood request out in the halo reports its expected system count (the profile's `lambda` summed over its cells) and its cost before it runs, since a cubic kiloparsec of empty halo is 1.6e7 cells at 40 ms each.

### 5.2 The rules that never change by place

1. **Address to profile.** A pure function of the galaxy shape, the component table and the address: `profile(address)`. No database read, closed form, the same code for an uncharted preview, a scatter and a fill.
2. **Count.** A sector's systems are a Poisson draw with mean `lambda x cell volume`, where `lambda` is the profile's total density (today's `args.density = relative_density(position)`, generalised).
3. **Population.** Each system picks one population in proportion to the profile's mix (`pick_population`), and its star is drawn from that population's age window and metallicity (`sample_living_star(population=...)`).
4. **Pre-placed objects.** The bright-star pass, the mass pass and the phenomena pass draw from the same profile, aggregated per shell by a certified bound. A sector's fill builds what was placed in it and draws only the rest (`FillContext`, `placed_star_fraction`), so its expected total does not change.
5. **Rare kinds.** Each kind's global rate per star times a factor the profile supplies (gas, young stars, metallicity; section 5.3).
6. **Components and wanderers add to the profile**; they do not have rules of their own (sections 5.5 and 5.6).
7. **Determinism.** A cell's draws come from its address and the galaxy seed, a shell's scatter from the seed and the shell, a component's from the seed and its id; none depends on what else was generated.

### 5.3 The sector profile: the parameters

Global settings are the same in every cell: the 4 pc edge, the bright-star floor (default 9,000 L_sun), the mass pass limit (default 14 solar masses), the compact-remnant limit, the seeds. What a place sets is below.

| Field | What it is | Where it comes from |
|---|---|---|
| `region` | disk, bulge, thick disk, halo, component (kind and id), ejecta-dominated, void | which term of the density dominates; a label for pages, not a rule |
| `lambda` | expected systems per cell | `E x relative_density`, with the halo term and the components |
| `mix` | density by population: young, intermediate, old, bulge, **halo**, **cluster**, plus one class per satellite | `population_densities`, extended |
| per population | age window, [Fe/H] mean and spread, velocity dispersion and net rotation | the table below |
| `gas` | cold-gas fraction against the local disk (0 to 1) and the ambient medium (disk ISM, hot halo, intergalactic) with its density in cm^-3 | the gas disk term, the hot-halo law, the cosmic mean (section 2.4) |
| `dark_matter` | density and mass in the cell | the NFW halo of `galactic-potential.md` |
| `kind_factors` | multiplier on each kind's rate per star | rules below |

**Population parameters** (the first four are today's; ages for the new ones follow section 2.1; [Fe/H] and dispersions are recalled, and a first model, not a fit):

| Population | Ages (Gyr) | [Fe/H] | Dispersion / rotation | Pre-placed at the default floor and cut [C] |
|---|---|---|---|---|
| young | 0 to 0.1 | about 0, falling 0.06 dex per kpc outward [R] | about 15 km/s, rotates | 3.1e-4 of its stars at 9,000 L_sun; 1.8e-4 born above 14 solar masses |
| intermediate | 0.1 to 3 | as young | about 25 km/s, rotates | 2.3e-6; none above 14 |
| old (thin and thick disk) | 3 to 10 | -0.1; thick disk about -0.6 [R] | 40 to 60 km/s, rotates | none |
| bulge | 8 to 12 | about -0.1, broad [R] | 100 to 140 km/s, bar | none |
| **halo** | 10 to 13 | -1.6 inside 25 kpc, -2.2 outside [S] | about 120 km/s, no rotation [S] | none; brightest star 2,400 L_sun |
| **cluster** | 10 to 13 | two peaks, -1.5 and -0.5 [globular-clusters.md] | 1 to 15 km/s inside, the cluster's orbit outside | none |
| **satellite** (one class per Cloud or dwarf) | its own mix of the first four | LMC about -0.4, SMC about -0.9 [R] | 20 to 30 km/s | as its mix |

The "pre-placed" column is the parameter that makes the scatters consistent: it is the share of a population's living stars that the passes already placed, so a cell's own draw is the rest. The repository's code computes it per population (`bright_star_fraction`, `massive_star_fraction`); adding a population is adding an age window and a row here. **A 10 to 13 Gyr population holds no star at or above 9,000 L_sun and none born above 14 solar masses** [C]: the bright-star and mass passes place nothing in the halo, a cluster or an old satellite, and one in 11,000 halo stars passes 1,000 L_sun [C].

**Kind factors** (rate per star times the factor):

| Kind | Factor | Why |
|---|---|---|
| Neutron stars, black holes (bound) | 1 for every population; the disk's kick thickening (`remnant_distribution`) applies to disk populations only | an old population holds 5.2e-3 neutron stars per living star, the disk 5.0e-3 [C]; a halo or cluster has no thin layer to thicken |
| Planetary nebulae, brown dwarfs | 1 | they follow stars |
| Supernova remnants | gas fraction | they need gas to be seen: 0 in the halo and in clusters |
| Molecular clouds, emission and reflection nebulae | gas fraction times the young share | no cold gas and no O to A stars in the halo, clusters or old satellites; the Clouds keep theirs |
| Rogue planets, comets | planet occurrence at the cell's [Fe/H] (`globular-clusters.md` slice 3); until that exists, 1 | metal-poor stars host far fewer planets; the generator has no metallicity yet (GEN.134) |
| Hypervelocity stars, kicked neutron stars and black holes | not per star: their own populations (section 5.6) | |

### 5.4 The scatters over the whole sphere

- **Bright-star and mass passes.** The algorithm does not change: one Poisson count per shell from a certified bound, bands in proportion, accept by true density over bound (`object_first.py`). It already costs about 1.2 microseconds per band [S]; to 262 kpc that is 6.75e9 bands, 2.25 hours per pass on one worker, or about 2 minutes with the bound taken over groups of 64 bands [C]. A group's bound is the bound of its nearest band, which the halo term gives exactly because it falls with the spheroidal radius; the efficiency loss is under 10% beyond 8 kpc [I]. Changes: the halo term in `band_majorants` (replacing the floor's share), the new populations in `fractions`, band groups beyond the outline. At the defaults the passes place nothing in the halo, clusters or old satellites; they place a few thousand stars in the Clouds (below).
- **Phenomena pass.** Today it walks every band of every shell (`_band_rows`), which is the part that cannot reach 262 kpc. Its classes are Poisson in density times a rate, the same structure as the star passes, so it moves to the same object-first sampler; the rates, the cut and the below-cut draw at sector time stay as they are.
- **Components** (clusters, Clouds, Sagittarius) are drawn by their own pass, not by the shell walk: their total star count from their mass, positions by inverse transform of their profile table (the King table of `globular-clusters.md` section 4), each star then kept or dropped by the same pass cut. The sector fill adds the component's density to the cell's `lambda` and the component's population to its `mix`, so a pre-placed cluster star and a fill-time one are never counted twice.
- **Wanderers**: section 5.6.

**What the passes yield over the sphere at the defaults** [C, with assumed mixes for the satellites]:

| Region | Systems | Bright (9,000 L_sun) | Mass pass (14) | Black holes above the cut | Planetary nebulae |
|---|---|---|---|---|---|
| Stellar halo | 3.75e9 | 0 | 0 | 7.5e5 (1.9e3 intermediate-mass) | 375 |
| LMC (1% young, 9% intermediate assumed) | 6.25e9 | 2.1e4 | 1.1e4 | 1.3e6 | 625 |
| SMC (same mix) | 1.3e9 | 4.4e3 | 2.3e3 | 2.6e5 | 132 |
| Sagittarius remnant | 2.5e8 | 0 | 0 | 5.0e4 | 25 |
| 157 globular clusters | 1.8e8 | 0 | 0 | 3.6e4 (an overcount: clusters keep few remnants [R]) | 18 |
| **Outside the disk, total** | | 2.6e4 | 1.3e4 | 2.4e6 | 1.2e3 |

Neutron stars never pass the 14 solar mass cut (they top out at 2.2), so every one is drawn at sector time; supernova remnants are zero in the halo, clusters and old satellites by the gas factor. About 2.4e6 extra phenomenon rows, 0.3 GB at the measured 138 bytes a row [S], against 6e7 for the disk: the sphere adds 4% to the scatter, not a new order of magnitude. The Magellanic rows depend entirely on the assumed young and intermediate share, which I did not find a source for (section 7).

### 5.5 Globular clusters and satellites, in the same train

Boss asked for the globular cluster work to ride with this. It does, as components under the profile rules, using the four slices of [globular-clusters.md](globular-clusters.md) and its decision to store a cluster as a density component, not as stars. What this adds to that note:

- **One component table, all kinds.** A row has a kind (cluster, satellite, stream later), position and velocity at an epoch, mass, scale radius, profile (King with `W0` for clusters, exponential for the Clouds and Sagittarius, Plummer for dwarfs), age, [Fe/H] and a core-collapse flag. The profile code, spatial index and density term are written once.
- **Milky Way default set.** About 157 clusters from the Harris catalogue [S], the LMC (49.6 kpc), SMC (62.4 kpc), Sagittarius (about 24 kpc), Fornax (143 kpc) and the other best-known dwarfs [S], with the embedding of `multiple-galaxies.md` section 2.2 so the sky agrees with the real one. **Boss's input is still needed for the catalogue file**: the repository has none, and this environment's network does not reach the Harris catalogue's host (the proxy refused `physwww.mcmaster.ca`). The options stand as in the cluster note: supply the table (public, about 157 rows), or allow that host for a one-time ingestion script. Until then the synthetic generator (slice 4) is built first and the real catalogue slots in when it arrives.
- **Clusters in the profile.** Inside a cluster's tidal sphere `lambda` is the cluster's density and the `mix` is 100% `cluster`. A typical core cell holds 3.0e4 to 8.4e4 systems and a 47 Tuc-like one 1.1e5 to 5.5e5 [S, `globular-clusters.md`]. At the measured 16.4 ms a system that is 14 minutes for the typical core cell and 0.5 to 2.5 hours for the 47 Tuc-like one [C], so an asked core cell is a queued job with an ETA (the existing work queue and progress bars), never a blocking request.
- **Cost of the pre-placed part is nil**, because no cluster star reaches the bright floor or the mass limit [C]: a cluster is entirely fill-time, bright-first backfill included.
- **Clusters need the metallicity value**: an [Fe/H] on stars, the planet cull, the pulsar-planet type and the blue-straggler count are slices 1 to 3 of the cluster note and appear unchanged in the build plan below. They also serve the halo, which has the same missing input.

### 5.6 Wanderers

- **Hypervelocity stars** keep their pre-placement. Direction uniform on the sphere, distance uniform in [0, reach] (density `r^-2`), speed 500 to 1,000 km/s; the outline test goes (they are wanderers, not disk members). With the code's 1e3 to 1e4 total read as "out to 100 kpc" (section 2.3), 84% land beyond the disk [C].
- **Kicked neutron stars and black holes** are an `ejecta` population: unbound fraction 20% to 40% of neutron stars and 2% of black holes [S, C], density `N / (4 pi r^2 v T)` out to `v T` (3 to 5 Mpc), with v the speed above escape (300 to 500 km/s) and T = 10 Gyr [I]. Against the bound halo neutron stars they take over near 50 kpc: ejecta per cell 1.3e-7 at 50 kpc, 3.3e-8 at 100 kpc, 4.8e-9 at 262 kpc (2e8 unbound, 300 km/s), 3.5 to 21 times the halo's at 100 kpc and 39 to 230 times at 262 kpc, depending on the assumptions [C]. They are not pre-placed (1e7 rows inside 262 kpc, 1e-8 per cell); they join the cell's `lambda` for neutron stars and black holes and are drawn at fill.
- GEN.232's cap on the integration of fast bodies comes first: a fast body in a sector makes the path integration run for minutes.

### 5.7 What an asked sector looks like in the far field

Examples at the default floor and cut [C]; expected systems and expected pre-placed rows in the cell:

| Sector | Systems | Bright rows | Mass-pass rows | What generation does |
|---|---|---|---|---|
| Disk, Sun-like | 6.3 | 6e-5 | 3e-5 | as today |
| Thick-disk fringe, 1.2 kpc up at 8 kpc | 0.76 | 0 | 0 | as today |
| Halo, 10 kpc above the centre | 2.9e-3 | 0 | 0 | Poisson draw; almost always saved as an empty sector, with the profile on its page |
| Halo, 50 kpc | 3.2e-5 | 0 | 0 | same |
| Halo, 262 kpc (edge) | 1.8e-8 | 0 | 0 | same |
| LMC centre (illustrative disk, 1.5 kpc by 0.4 kpc) | 35 | 1e-4 | 6e-5 | ordinary fill; young and intermediate populations present |
| Typical cluster core cell | 5e4 | 0 | 0 | queued; cluster population, metal-poor |

An empty halo sector is not an error: it is the answer, stored so that a map and a bookmark can name it, with the expected counts, the ambient medium and the dark-matter figure computed from the address (section 5.8). A route planner meets the same gaps: beyond about 30 kpc the typical spacing of halo stars passes the 50 pc hop limit (section 4).

### 5.8 Gas and dark matter on the sector page

Computed from the address, not stored, like the interstellar debris count already shown (`queryDb.interstellar_debris_count`): the dark-matter density and mass in the cell, the ambient medium (disk ISM, hot halo, intergalactic) with its density and the gas mass in the cell, and the profile's expected counts. Numbers are in sections 2.4 and 2.6. High-velocity clouds are the one gas item worth rows: a few hundred, placed like the molecular cloud field (`nebula_field`) at 2 to 15 kpc and the Magellanic Stream at 55 kpc and beyond [S, R]; a later slice. `nebula_field` must return no cold gas outside the gas disk; I did not check that it does.

### 5.9 What it changes in the code

- `galaxy/profile.py` (new, pure): `profile(address)`, built on `density.py`; `FillContext` takes a profile.
- `galaxy/density.py`: halo term, the `halo`, `cluster` and satellite populations, floor removed (its tests rewritten); the density mirror in `galaxyprisms.js` may ignore the halo (densities below 1e-3 are not drawn) but must tolerate the new population keys.
- `galaxy/skeleton.py`: `GalaxyBounds` gains the reach next to the outline; `contains` answers for the reach, `qualifies` for the outline; `random_address` stays on the outline.
- `galaxy/spherical.py`, `galaxy/object_uid.py`: designation functions take the galaxy's layout; the plan stores the reach and refuses what the layout cannot hold (this is GEN.230's layout decision).
- `generation/run_galaxy.py`: `ensure_sector_generated` loses the outline refusal and `outside_ok`; a request in the reach is generated from its profile and saved.
- `generation/object_first.py`, `bright_stars.py`: halo term in `band_majorants`, band groups, new populations.
- `generation/phenomenon_scatter.py`: phenomena on the object-first sampler; hypervelocity stars without the outline test; the ejecta population; the component pass.
- New: component table, King tables and spatial index (cluster slice 1 and 2); metallicity value, planet cull, pulsar planets, blue stragglers (slice 3).
- Pages and API: the sector page shows the profile; a read-only profile lookup for an uncharted address.
- Tests: tiling and ID round trip at the widest shell of the reach; determinism of a profile and of an empty sector; the pre-placed share of every population; King solver against the textbook concentrations; no double count of a cluster star.

### 5.10 Nebulae that need a star

Boss asked whether the required objects are generated when a nebula is made. Answer, from the code at this checkout [S]: for planetary nebulae and supernova remnants, yes, always; for emission and reflection nebulae, only when the host star already exists in the sector; otherwise the nebula carries the host as text; molecular clouds need none; there is no Wolf-Rayet class. The table says what each kind requires, what the code does, and what the sphere-wide design adds.

| Nebula | Required inside | Today | Gap | Far from the disk (the design) |
|---|---|---|---|---|
| Planetary, H to L | exactly one hot central star (exposed core 0.55 to 0.9 solar masses, 30,000 to 200,000 K), often with a close companion | **real, always.** `_add_planetary_nebula` builds a new system whose only star is a white dwarf drawn from `PLANETARY_NEBULA_CENTRAL_STAR_TYPES`; the scatter rows take the same path (`_build_scattered`). The sector can refuse it when full (the nebula is then dropped) | the draw has no companion; central stars are O3VII to B0VII, heavier than the real 0.59 solar masses (the correction table already notes this) | **kept everywhere**: old stars make them, so halo, cluster and old-satellite cells keep them at the old-star rate (a few halo ones are known [R]) |
| Supernova remnant, core collapse (R, S, T, U, V) | a neutron star (about 85%) or black hole (15%), offset by its kick, or none | **real.** `SupernovaRemnant.__init__` embeds one in 60% of core-collapse remnants (75% of all remnants); the scatter uses the same constructor | none for the hidden 40%: that is the intended "no detectable core" case | **zero in old populations** (no progenitor); needs the young share |
| Supernova remnant, Type Ia (W) | no compact object | **real**: none embedded (25% of remnants) | the Type Ia share is a constant (`SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE = 0.25`), right for the disk, wrong for old populations | **share 1 in old populations**, the only remnant they make |
| Emission (H II), C, D, E | one O or early-B star (C), a small cluster (D), a massive young cluster (E); G an early-B star | **real only through `add_star_hosted_nebulae`**: it grows the nebula around a hot main-sequence star the sector already holds (O at 0.25, B0 to B2 at 0.25; class age limits 1, 10, 20, 5 Myr). A nebula made any other way has the host in its `hosts` text only, with **no star object**; the cluster members of C, D, E, G are text only | the scatter does not place emission or reflection nebulae; a host star that is not in the sector (pre-placed before its sector exists) gets no nebula (open MAP.182, default do not build) | **none unless the cell is young**: halo, cluster and old-satellite populations have no O or B stars [C], so the host rule yields nothing; the Magellanic mix can |
| Reflection, F | a B or A star, too cool to ionize (B3 to B9 at 0.05, A at 0.01; age limit 50 Myr) | the same host-first rule; text-only host otherwise | same | same |
| Molecular and dark clouds, M, N, P, Q | none illuminating; P none or one protostar, Q embedded protostars | clouds come from the seeded `nebula_field`; **no protostar objects are generated** | embedded protostars are text only | needs the cell's `gas`; none in halo, clusters or old satellites; the Magellanic mix can |
| Wolf-Rayet bubble | a Wolf-Rayet star | **no class exists** (`anomalies.md`, Tier 2) | whole feature | young cells only |

The rule to build, one for every cell: **a nebula's required objects are decided by its kind, not by where it is, and they are made when the nebula is made.** Concretely:

1. Kinds with a point rate (planetary, supernova remnant) keep building their core (white dwarf, neutron star, black hole) at creation, in any cell.
2. Emission, reflection and Wolf-Rayet nebulae are *hosted*: the host star is chosen first, from the cell's young population, and the nebula is grown around it. The hot-star scatters (bright and mass passes) are the sphere-wide host supply: where they place nothing (halo, old clusters), no hosted nebula exists. The text-only `hosts` list for the cluster classes (D, E, G) becomes real member stars drawn from the cell's young `mix` when the nebula is created, or the class is not offered in that cell.
3. Which kinds a cell can have follows from the profile (section 5.3): `kind_factors` scale hosted nebulae and core-collapse remnants with the young share, planetary nebulae and Type Ia remnants with the old share, molecular clouds with `gas`; the Type Ia share becomes a profile field (1 where there is no young population).
4. An asked halo cell therefore never shows an H II region or a pulsar-bearing remnant, and may show a planetary nebula with its white dwarf or a Type Ia remnant with nothing inside. A Magellanic cell can show all of them, hosts included.

What I did not check: what `nebula_field` returns beyond the disk today; the build must make it read the cell's `gas`, not the disk position alone.

## 6. What it costs

| Item | Cost | Tag |
|---|---|---|
| Storing every cell of the 16.4 kpc sphere | 2.88e11 cells x 40 ms = 365 years one worker; about 29 TB for bare 100-byte rows | C |
| Storing every cell of the 262 kpc sphere | 1.18e15 cells: 1.5e6 years and 1.2e5 TB | C |
| Uncharted cell | nothing stored; a profile is a few microseconds [I] | I |
| An asked halo sector | about 40 ms and one row (GEN.76) | S |
| An asked cluster core cell | 14 minutes (typical) to 2.5 hours (47 Tuc-like), queued | C |
| Bright and mass passes to 262 kpc | 2.25 hours per pass one worker as built, about 2 minutes with band groups | C |
| Phenomena pass to 262 kpc | not possible as built (band walk); on the object-first sampler the same as a star pass | I |
| Extra scatter rows outside the disk | about 2.4e6 black-hole rows and a few thousand others, 0.3 GB | C |
| Wider object IDs | 80 to 96 bits: two more bytes on each ID column and index, only for a galaxy planned with a reach | S |

## 7. Honest uncertainty

- **The halo's normalisation** is the weakest number: the stellar mass (1.5 +- 0.4e9, from a search summary of Deason et al. 2019) and the local density (about 0.5% recalled) disagree by 2.7 times. Either way every halo cell more than 9 kpc from the core is below 0.012 systems, and beyond 20 kpc below 2e-3 [C].
- **Slopes and breaks** are still disputed: breaks from 10 to 30 kpc across studies (a 2026 DESI paper says "doubly broken", twisted) [S]. A single broken power law is a first model, not a map.
- **The kick numbers are recalled**, so the 20% to 40% unbound fraction is a range from one population synthesis [S] and a recalled kick law [R, C]. Uniform-in-time production overcounts if star formation peaked early.
- **The Magellanic populations**: the 10% young and intermediate share is an assumption; the bright and mass rows scale with it directly.
- **Satellite census and masses** are search summaries, not a catalogue check; positions, sizes and the dwarf list need McConnachie (2012) or the current Gaia census before any row ships.
- **Metallicity values and dispersions** in section 5.3 are recalled and are inputs the generator does not have yet; until slice 1 lands, halo and cluster stars would still receive full planet systems.
- **Cluster remnants**: the pre-placed cluster black holes (3.6e4) use the field's rate, which overcounts; clusters retain few.
- **I did not model** high-velocity cloud counts, the hot gas density profile, intragroup light, stellar streams or Andromeda's halo; each is either too small to alter a cell (streams) or outside the reach (Andromeda).
- **Not checked**: that `nebula_field` returns no gas outside the disk; that no page needs a pre-existing sector row to show an uncharted cell.

## 8. Decisions

Decided (Boss, 2026-10-15): the 262 kpc reach with 96-bit IDs; the same rules everywhere with per-sector parameters; nothing generated unless asked, an asked sector always generated; the first neighbourhood pull in the galactic region; globular clusters in this build. Still open, with defaults the build proceeds on:

1. **Catalogue file** for the Milky Way's clusters (section 5.5): supply the Harris table, or allow its host for a one-time ingestion script. Default: the synthetic generator first, real catalogue when it arrives.
2. **Satellites.** Real Milky Way satellites at their real places for the default shape, seeded fictional ones for other shapes later. Default: as stated.
3. **The floor.** Deleted; the halo and components replace it with no setting. Default: delete.
4. **Normalisation.** The 1.5e9 solar mass anchor, with the local density (about 2.7 times higher) kept as the stated alternative. Default: mass anchor.
5. **Magellanic mix.** 1% young and 9% intermediate until a better source is found. Default: as stated.

## 9. Build plan for TODO

One train, in build order. Titles only; the TODO thread numbers them and assigns lanes. Existing items are named where they come first.

1. **Grid reach and layout** (Core; resolves GEN.230's layout question): plan-time reach, 96-bit layout from `fit_layout`, designation functions with the layout, `GalaxyBounds` reach and outline, `random_address` on the outline, refusal text, tests at the widest shell.
2. **Sector profile** (Core): `galaxy/profile.py` with `lambda`, `mix`, per-population parameters, gas, dark matter and kind factors; `ensure_sector_generated` generates any asked cell in the reach from it and saves it; `outside_ok` and the outline refusal removed.
3. **Halo density term and populations** (Core): broken power law, `halo` population, floor removed, `fractions` and majorants updated, tests.
4. **Metallicity value** (Core; GEN.134 and cluster slice 1): [Fe/H] on stars with the gradient, the halo and cluster values, the planet cull, the sector page showing it.
5. **Component table** (Core): clusters, Clouds, Sagittarius, dwarfs in one table; King, exponential and Plummer profiles; King tables; spatial index; the density term in the profile and in the bounds; the Milky Way default set with the GEN.9 embedding (cluster slice 1 and 2; the catalogue file when supplied).
6. **Scatters over the sphere** (Core; after GEN.232): band groups in the object-first bound, phenomena pass on the object-first sampler, component pass, hypervelocity stars without the outline test, the ejecta population; tests for no double count of a pre-placed star.
7. **Globular cluster physics** (Core; cluster slices 3 and 4): the hard-orbit planet cut, type-B pulsar planets, blue stragglers, synthetic clusters and satellites for other shapes.
8. **Sector page and API** (Interface): the profile on the page, a read-only profile lookup for an uncharted address, gas and dark-matter lines, the heavy-cell queue warning with ETA, the neighbourhood request's expected-count line.
9. **Galaxy Map** (Interface): components and wanderers as points beyond the outline; shading unchanged.
10. **High-velocity clouds and the Magellanic Stream** as placed gas clouds (Core; later).
11. **Required contents per nebula kind** (Core; section 5.10; touches MAP.182, whose default stays "do not build" unless Boss asks): hosts made with the nebula (real member stars for classes D, E, G and embedded protostars for P and Q, or the class withheld), Type Ia share and kind rates read from the profile, `nebula_field` gated by the cell's `gas`, Wolf-Rayet class with the anomalies lane.
12. **Docs** (Quality): `galaxy-disk-density.md`, `galaxy-coordinate-system.md` (reach and outline), `globular-clusters.md` (status), `interstellar-object-rates.md` (populations and factors).
13. **Tests** (Quality): no emission, reflection or core-collapse remnant in a halo or old cluster cell; a planetary nebula always has its white dwarf; a core-collapse remnant has a compact core at the stated share; every hosted nebula holds its host star; tiling at the reach, ID round trip at the widest shell, profile and empty-sector determinism, the pre-placed share per population.

## Sources

All read as search-result summaries; none opened in full.

- Halo profile and breaks: [Broken and unbroken: the Milky Way and M31 stellar haloes](https://arxiv.org/pdf/1210.4929); [On the Shoulders of Giants](https://arxiv.org/pdf/1408.1787); [Milky Way Tomography with Subaru HSC II](https://arxiv.org/pdf/2604.10280); [DESI DR2 doubly broken halo](https://arxiv.org/pdf/2512.01350).
- Halo metallicity and age: [H3 survey halo metallicity](https://arxiv.org/pdf/1909.02007); [BHB age gradient and inner halo assembly](https://arxiv.org/pdf/1909.12910); [Review of current estimates of the Galaxy mass](https://arxiv.org/pdf/2305.18408) (halo mass in a summary that cites Deason et al. 2019).
- Hypervelocity stars: [Improved constraints on Galactic Centre ejection of hypervelocity stars](https://arxiv.org/pdf/2406.14134).
- Intracluster light: [Intracluster planetary nebulae in Virgo](https://arxiv.org/abs/astro-ph/0407274); [Intracluster stars in the Virgo core](https://ar5iv.arxiv.org/html/astro-ph/0502420); [Intracluster light in Virgo, large-scale distribution](https://arxiv.org/pdf/0908.3848).
- Globular clusters: [A New Catalog of Globular Clusters in the Milky Way](https://ar5iv.arxiv.org/html/1012.3224).
- High-velocity clouds: [HVCs as streams of ionized and neutral gas](https://arxiv.org/pdf/1203.2626); [High-velocity cloud](https://en.wikipedia.org/wiki/High-velocity_cloud); [Accurate distance to Complex C](https://ar5iv.arxiv.org/html/0712.0612).
- Hot circumgalactic gas: [The Structure of the Milky Way's Hot Gas Halo](https://arxiv.org/pdf/1305.2430); [The Circumgalactic Medium](https://arxiv.org/pdf/1709.09180); [Baryon budget of the hot CGM](https://ar5iv.arxiv.org/html/1802.09453).
- Warm-hot intergalactic medium and the Local Group: [Chandra and FUSE view of the WHIM](https://ar5iv.labs.arxiv.org/html/astro-ph/0311162).
- Satellites: [The Massive Satellite Population of Milky-Way Sized Galaxies](https://arxiv.org/pdf/1306.4328); [The mass of our Milky Way](https://arxiv.org/pdf/1912.02599); [LMC mass from S5 streams](https://arxiv.org/pdf/2107.13004).
- In this repository: `galaxy-disk-density.md`, `interstellar-object-rates.md`, `compact-remnant-regions.md`, `globular-clusters.md`, `multiple-galaxies.md`, `galactic-potential.md`, `sampling-backfill-and-resume.md`, `spherical-sector-coordinates.md`, `galaxy-coordinate-system.md`; `src/planetgen/galaxy/{density,skeleton,spherical,object_uid}.py`, `generation/phenomenon_scatter.py`, `tuning.py`.
- Computations: `/mnt/project-files/research/scripts/halo/halo_numbers.py` and `halo_numbers.txt`.
