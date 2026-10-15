# Filling the whole sphere: the halo, the satellites and intergalactic space

Boss (2026-10-15): "investigate filling the entire spherical grid of sectors but when it's off the galactic plane or the proper galaxy, we'll have to adjust the sector contents, since clearly it'd be less ... what is known to be found in intergalactic space, stars, phenomena". This note says what is known to exist out there, what the generator does with that space today, and proposes a concrete way to give every cell of the sphere defined contents without storing the empty ones.

Informs: GEN.230 (the band field and the birth-field layout), GEN.232 (fast bodies), GEN.9 and VIEW.2 (neighbour galaxies), GEN.78 (the density floor), GEN.100 and GEN.223 (scatter and the nightly motion of every object), GEN.134 (metallicity), GEN.42 (sampling), the globular cluster slices of [globular-clusters.md](globular-clusters.md), DB.32 follow-ups (object ID layout)
Status: research, 2026-10-15; nothing here is built. The proposal in section 5 waits on Boss's word on the decisions in section 8; everything else is a recommendation.

Evidence tags: **[S]** seen: in this repository, or in a search-result summary read during this research (sources at the end; the papers themselves were not opened, the environment cannot reach arXiv or journals); **[C]** computed here with `/mnt/project-files/research/scripts/halo/halo_numbers.py` (output beside it in `halo_numbers.txt`; numpy and the repository's own grid code only); **[R]** recalled and unconfirmed; **[I]** inferred from tagged facts, with the step shown. A figure keeps the tag of its weakest input.

Related notes: [spherical-sector-coordinates.md](spherical-sector-coordinates.md) and [galaxy-coordinate-system.md](galaxy-coordinate-system.md) (the grid), [galaxy-disk-density.md](galaxy-disk-density.md) (the density model), [interstellar-object-rates.md](interstellar-object-rates.md) and [compact-remnant-regions.md](compact-remnant-regions.md) (kinds and rates), [phenomenon-scatter-mass-cut.md](phenomenon-scatter-mass-cut.md), [sampling-backfill-and-resume.md](sampling-backfill-and-resume.md) (cost of sparse sectors), [globular-clusters.md](globular-clusters.md), [multiple-galaxies.md](multiple-galaxies.md), [galactic-potential.md](galactic-potential.md), [object-ids.md](object-ids.md).

## Summary

**What is out there.** Almost nothing, by the generator's own yardstick. A sector that "exists" today expects at least one system. The stellar halo is about 1.5e9 solar masses [S], 3% of the disk and bulge, and spread so thinly that a cell at the Sun's distance from the centre expects about 0.01 systems [C]; at 20 kpc, 4e-4; at 50 kpc a cell expects 3e-5 [C]; at the 262 kpc edge of the Milky Way's neighbourhood, 2e-8 [C]. The intergalactic average is 1e-7 systems per cell or less [C, R]. Only discrete things reach one system per cell: the Magellanic Clouds, the Sagittarius remnant, dwarf spheroidals and globular clusters. Everything else is a sparse Poisson haze of old, metal-poor stars, ejected stars and ejected neutron stars, with a gas background of fractions of an Earth mass per cell.

**What the project does today** [S, C]:

1. **The grid is already a tiling of all space, but its addresses stop short.** The 80-bit object ID cannot name a cell whose band number passes 2,047, which means every shell from 1,304 (5.2 kpc) outward is only partly addressable: above latitude 31 degrees at the disk's edge (shell 3,764), no cell has an ID. Shells end at 4,095 (16.4 kpc). The real halo, the Magellanic Clouds (50 kpc) and the Milky Way's virial radius (about 200 kpc) are outside.
2. **The "halo" is a flat floor, and it is wrong in both directions.** `MIN_RELATIVE_DENSITY` is 0.1% of the local density everywhere (0.0063 systems per cell). It is half the real halo at the Sun, 200 times too high at 50 kpc, and 340,000 times too high at 262 kpc [C]. If the grid were simply filled out to 262 kpc under the floor, it would hold 7.4e12 phantom systems, 25 times the real galaxy's 3.0e11 [C].
3. **Nothing beyond the disk outline can be generated**, the hypervelocity stars are squeezed inside it (they are redrawn when they land outside), and the bright-star floor (9,000 L_sun) finds no halo star at all [C].

**Recommendation** (section 5, decisions in section 8):

- **Fill by definition, not by storage.** Every cell inside a chosen reach gets an expected system count from the density model, realised on demand by a seeded Poisson draw. An empty cell is never stored. The qualifying outline stays as it is, now only an optimisation (map shading, pre-placed scatter, random starts), not the definition of a sector.
- **Reach 65,535 shells = 262 kpc** (1.18e15 cells), with 96-bit object IDs. It contains the Magellanic Clouds, the Sagittarius remnant, the globular clusters, the satellites out to Leo I and the Milky Way's virial radius [S]. Andromeda and the rest of the Local Group stay sky objects (GEN.9 stage 1).
- **Replace the flat floor with a measured stellar halo**: a broken power law, 1.5e9 solar masses, 3% of the visible mass, old and metal-poor. It lowers the cell density out there by up to five orders of magnitude and, by itself, makes no sector qualify.
- **Four kinds of content beyond the disk**, each cheap: the halo itself (density term), discrete components (satellites and clusters, a small table), wanderers (hypervelocity stars, kicked neutron stars) with analytic radial laws, and computed gas and dark-matter figures on the sector page instead of rows.
- **Most far sectors stay implicit.** Storing every cell of the 16 kpc sphere would take 365 years of single-worker time at the measured fixed cost per sector and 29 TB for bare rows [C]; the 262 kpc sphere is out of the question. Only occupied cells (about 3.7e9 halo systems plus components) ever get rows, and only when visited.

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

**Stars you would actually find** [C]: sampling 60,000 living stars of age 10 to 13 Gyr with the repository's own stellar model gives a brightest star of 2,400 L_sun and 0 above 9,000 L_sun (1 in 7,500 above 1,000). The default bright-star scatter (9,000 L_sun) therefore places no halo star, and the dim halo is entirely the sector fill's job. The repository does not model horizontal-branch stars, which dominate halo star counts at the bright end in reality [S, `globular-clusters.md`].

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

Gas at these densities is background, like the diffuse classes A and B the generator already leaves ungenerated. High-velocity clouds are the one gas item worth a row (section 5.5).

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
- **Scatter.** `phenomenon_scatter` and `bright_stars` walk shell by shell and band by band inside the outline. A walk of every band to 262 kpc would be about 6e9 band steps (the sum over shells of pi times the shell number), too many for the numpy walk that serves the disk.
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

## 5. The proposal

### 5.1 Reach and object IDs

| Option | Reach | Cells | ID | What it holds | Verdict |
|---|---|---|---|---|---|
| A | Disk outline only (as built, shell 4,095 limit) | 2.4e10 qualifying | 80 | Nothing beyond the disk | Does not meet the request |
| B | 262 kpc, 65,535 shells | 1.18e15 | 96 | Halo, LMC, SMC, Sagittarius, clusters, satellites out to Leo I, virial radius (199 kpc [S]) | **Recommended** |
| C | 1.05 Mpc, 262,143 shells | 7.5e16 | 128 | Adds the Andromeda region at 770 kpc, which would then need its own generated galaxy | Not now |

The reach is a plan-time setting (a new `galaxy_shape` value, default off for the 80-bit case), stored with the galaxy; `outline_object_id_layout` takes it into account so the plan refuses a reach the layout cannot hold. A cell beyond the reach is refused as today ("outside the grid"), which also settles what happens to a wanderer that flies past the edge: at 500 km/s it needs 510 Myr to cross 262 kpc [C], so none does within any plausible game clock. An existing 80-bit galaxy gets the halo by re-planning or by a renumbering migration like the shell move's (no backward compatibility).

### 5.2 Sector content outside the outline: the rule

Three tiers, one rule. A cell's contents are always a pure function of the galaxy seed and its address (`galaxy/seed.py`); what changes is how much the generator is willing to store.

1. **Qualifying (inside the outline, expectation at least one).** As today: generated when visited, filled by the queue, drawn on the map.
2. **Sparse (inside the reach, expectation below one).** Never queued, never pre-filled. On a visit the generator draws the Poisson count for the cell's expected systems (`lambda x CELL`, `lambda` from the density model with the halo term and the components below), builds that many systems and the cell's rare objects, and stores rows only if there are any. An empty visit stores nothing and costs one 0.06 ms draw [S] plus the request. The sector page computes "no objects; expected 3e-5 systems here" and the background figures of section 5.5 from the address alone.
3. **Components (anywhere in the reach).** Satellites and clusters add to `lambda` and so lift their own cells above one (section 5.4); inside one the cell is tier 1.

This drops "a sector exists where it expects at least one star" as the definition of a sector and keeps it as an optimisation, a cached bound on where the sparse tier starts. The consequence for code paths is small, and listed in section 5.7.

### 5.3 Which kinds occur, and at what rates

Rates are per star (`phenomenon_rate_per_star`), so the existing rule carries over; the table says where a halo cell differs.

| Kind | Halo rule | Count in the halo (3.75e9 systems) | Why |
|---|---|---|---|
| Ordinary stars | New population `halo`, age 10 to 13 Gyr, [Fe/H] -1.6 inside 25 kpc falling to -2.2 outside; velocity dispersion 120 km/s, no rotation | 3.75e9 | sections 2.1, 4 |
| Bright stars (scatter) | No change; the 9,000 L_sun floor finds none [C]. A lower floor would find halo giants only in the thousands | 0 at the default | 2,400 L_sun is the brightest star |
| Neutron stars | Same per-star rate. An old population has 5.2e-3 neutron stars per living star against 5.0e-3 for the 10 Gyr disk [C, Kroupa IMF recalled] | 1.5e7 | dead massive stars accumulate |
| Black holes | Same per-star rate (lower metallicity gives heavier ones [R]) | 1.9e6 | as above |
| Hypervelocity stars, kicked neutron stars | Not per star: radial laws, section 5.6 | 1e3 to 4e4; 1e7 | 2.3 |
| Planetary nebulae | Same rate; the count is only about 400, so the factor does not matter | ~400 | |
| Supernova remnants | Zero: they need gas to be seen [I]; the count at the disk rate would be about 40 | 0 | |
| Molecular clouds, emission and reflection nebulae | Zero: no cold gas, no O to A stars (except inside the Magellanic and dwarf components) | 0 | |
| Rogue planets, comets | Scale with the planet occurrence of the host: the halo has 3% of the solar metals ([Fe/H] -1.6), and giant-planet occurrence drops steeply below -1 (`globular-clusters.md`). Use the planet cull of the globular cluster slice 3 once it exists; until then any star system in the halo carries a full planet set, which is the one known wrong thing | up to 2.4e10 at the disk rate; far fewer in practice | the generator has no metallicity (GEN.134) |
| Brown dwarfs | Same per-star rate | 8e8 | they form like stars |
| Runaway flag | Keep it relative to the neighbours, not to the disk | | halo stars are fast against the disk, not against each other |
| Asteroid fields, GMC | As now (zero) | | |

### 5.4 Components: satellites and clusters

One small table of discrete density components, generalising the cluster record of `globular-clusters.md` section 3 (position, velocity, epoch, mass, scale length, profile, age, metallicity, kind). The Milky Way default would ship about 160 clusters (the Harris catalogue [S]), the LMC, the SMC, the Sagittarius remnant and the best known dwarfs (Fornax, Draco, Sculptor, ...), placed with the embedding of `multiple-galaxies.md` section 2.2. Profiles: King for clusters (already designed), exponential for the Clouds and Sagittarius, Plummer for dwarfs. A spatial index (the cluster note's) keeps the density call cheap. Draws and the block-first sampler stay as they are. A galaxy with a non-default shape gets seeded fictional satellites from a luminosity function (the GEN.9 option), later.

Mass check [C]: LMC 6.2e9 systems, SMC 1.3e9, Sagittarius 2.5e8, a cluster set 1.8e8; the Clouds carry most of it. Streams are left out of the first release: the densest one expects 0.015 systems per cell, below every threshold, so a stream tube would change what a visitor sees by almost nothing [C on R].

### 5.5 Gas and dark matter on the sector page

Computed, not stored, like the interstellar debris count already shown (`queryDb.interstellar_debris_count`): the page prints, from the address, the NFW dark-matter density and mass in the cell, a gas density class (disk, hot halo, intergalactic medium) with its figure, and the expected object counts. The numbers are in section 2.4 and 2.6. High-velocity clouds are the one gas item worth rows: a few hundred, pre-placed like the molecular cloud field (`nebula_field`) at 2 to 15 kpc and the Magellanic Stream at 55 kpc and beyond [S, R]. That is a later slice.

### 5.6 Wanderers

- **Hypervelocity stars.** Draw direction uniformly on the sphere and distance uniformly out to the reach, as the code does, but drop the outline test (they are wanderers by definition) and use the ejection rate and flight-time law of section 2.3: inside a radius r there are `rate x r/v` of them, so the one-stage draw is already right for a total and the placement is a pure function of the seed. At 500 to 1,000 km/s and the code's 1e3 to 1e4 stars, spreading them uniformly in distance out to 100 kpc rather than 15 kpc puts 84% of them beyond the disk [C], which is where the flight-time law says they are. The total should be read as "out to 100 kpc".
- **Kicked neutron stars.** 20% to 40% of the isolated ones, a radial law `r^-2`, speeds 300 to 500 km/s above escape, random direction. Not pre-placed (about 1e7 rows inside 262 kpc, 1e-8 per cell): they are drawn when a cell is visited, by the same Poisson rule as every other far object, with the expected count from the radial law. A black hole is the same at 2%.
- The orbital update (GEN.223) then moves them like everything else; positions are stored with the epoch and velocity, so nothing special is needed.

### 5.7 What it changes in the code

- `galaxy/density.py`: halo term, `halo` population, floor removed; `tuning` constants; the density mirror in `galaxyprisms.js` can ignore the halo (densities below 1e-3 are not drawn) but must not break on a new population key.
- `galaxy/skeleton.py`: `GalaxyBounds` gains a `reach` next to the outline. `contains` answers for the reach, a new `qualifies` answers for the outline. `random_address` stays on the outline (on the reach, 99.998% of draws would be empty halo).
- `galaxy/spherical.py`, `galaxy/object_uid.py`: designation functions take the galaxy's layout; the plan stores the reach and refuses what the layout cannot hold.
- `generation/run_galaxy.py` (`ensure_sector_generated`): the sparse tier: Poisson count, build, store only if occupied.
- `generation/phenomenon_scatter.py`: unchanged for the outline; the halo's discrete kinds (a few thousand black holes above the cut) get an analytic radial draw, not a band walk (the walk to 262 kpc is 6e9 steps); hypervelocity stars lose the outline test.
- Population-dependent code (`bright_stars`, `object_first`, the sector fill) gains the `halo` population key.
- Tests: the density-floor tests are rewritten; new tests for tiling at the reach, ID round trips at the widest shell, determinism of an empty visit.

## 6. What it costs

| Item | Cost | Tag |
|---|---|---|
| Storing every cell of the 16.4 kpc sphere | 2.88e11 cells x 40 ms = 365 years one worker; about 29 TB for bare 100-byte rows | C |
| Storing every cell of the 262 kpc sphere | 1.18e15 cells: 1.5e6 years and 1.2e5 TB | C |
| A draw-only pass over the 16.4 kpc sphere | 200 days one worker | C |
| One empty visit (proposal) | 0.06 ms draw plus the request; 0 rows | S |
| Occupied halo cells (proposal) | About 3.7e9 halo systems plus components; each row only when someone visits | C |
| Wider object IDs | 80 to 96 bits: two more bytes on each ID column and index, in every table that holds one; only for galaxies planned with a reach | S |
| Halo term in `relative_density` | One more exponentiation per call, a few percent on the density call [I] | I |
| Pre-placed halo scatter | A few thousand rows (black holes above the cut, hypervelocity stars); nothing like the 1.17e9 of the disk | C |

## 7. Honest uncertainty

- **The halo's normalisation** is the weakest number: the stellar mass (1.5 +- 0.4e9, from a search summary of Deason et al. 2019) and the local density (about 0.5% recalled) disagree by 2.7 times. Either way every halo cell more than 9 kpc from the core is below 0.012 systems, and beyond 20 kpc below 2e-3 [C].
- **Slopes and breaks** are still disputed: breaks from 10 to 30 kpc across studies (a 2026 DESI paper says "doubly broken", twisted) [S]. A single broken power law is a first model, not a map. Beyond 100 kpc the sources give a slope, not a measurement of the stars.
- **The kick numbers are recalled**, so the 20% to 40% figure for unbound neutron stars is a range from one population synthesis [S] and a recalled kick law [R, C]. Uniform-in-time production overcounts if the star formation peaked early.
- **Satellite census and masses** are search summaries, not a catalogue check. The positions, sizes and the number of dwarf galaxies need McConnachie (2012) or the current Gaia census before any row is shipped.
- **Metallicity** is a missing input: the proposal's halo is old and metal-poor on paper, but the generator would still give its stars full planet systems until GEN.134 and the cluster slice 3 land.
- **I did not model** high-velocity cloud counts, the hot gas density profile, intragroup light, stellar streams or Andromeda's halo; each is either too small to alter a cell (streams) or outside the reach (Andromeda).

## 8. Decisions for Boss

Defaults are the recommendation; the thread proceeds on them unless told otherwise.

1. **Reach.** B (262 kpc, 96-bit IDs, plan-time setting) over A (nothing) or C (1.05 Mpc, 128-bit). Needed before any code.
2. **Empty visits.** Store nothing (the page is computed) over storing a row for every visited cell. A row costs about 100 bytes per visit and gains only a "visited" mark; the one thing to check before building is whether any page needs a sector row to exist (sector names).
3. **The floor.** Delete it and let the halo and components replace it, with no setting.
4. **Satellites.** Real Milky Way satellites at their real places for the default shape (Magellanic Clouds, Sagittarius, clusters, a few dwarfs), seeded fictional ones for other shapes later.
5. **Normalisation.** The 1.5e9 solar mass anchor, with the local density (about 2.7 times higher) kept as the stated alternative.

## 9. Proposed TODO items

Titles only; the TODO thread numbers them and sets lanes. Order is the build order.

1. Grid reach setting (builds on GEN.230, whose layout decision this answers): plan-time `reach`, `GalaxyBounds.contains` for the reach and `qualifies` for the outline, `random_address` on the outline, designation functions with the layout, 96-bit layout for a full sphere (Core; migration only for existing galaxies).
2. Stellar halo density term: broken power law, `halo` population, floor removed, tuning constants, tests (Core).
3. Sparse sector visit: seeded Poisson on visit, no row for an empty cell, sector page text and computed gas and dark-matter figures (Core, Interface).
4. Halo rules for kinds: per-kind multipliers of section 5.3, zero for gas-dependent kinds, hookup to the metallicity planet cull (Core; depends on GEN.134 and the cluster slice).
5. Satellite and cluster components: one table, King/exponential/Plummer profiles, spatial index, Milky Way default set with the GEN.9 embedding (Core; shares slices with `globular-clusters.md`).
6. Hypervelocity stars to the reach with the `r^-2` law, and kicked neutron stars and black holes as a radial-law draw (Core; GEN.232's cap on the path integration comes first, since a fast body there runs for minutes).
7. High-velocity clouds and the Magellanic Stream as placed gas clouds (Core, later).
8. Galaxy Map: draw components and wanderers as points beyond the outline, no change to shading (Interface).
9. Docs: `galaxy-disk-density.md` "not modelled: the stellar halo beyond the floor" and `galaxy-coordinate-system.md` section 4 (reach versus outline) (Quality).

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
