# Mass backfill: four rings of heavier stars around the generated sectors

Boss (GitHub issue #952, defaults confirmed 2026-10-10 04:06Z): replace the GEN.30 luminosity
tiers with four rings of sectors around what has been generated, each ring taking the stars
born above a mass, so the sky of a generated sector is full of the stars you could see from it
and the cost falls off with distance. Items: GEN.187. Builds on GEN.183 and GEN.185 (the
scatter's mass and luminosity limits) and GEN.44 (a level per sector).

Status: built (schema v79, `sector_stats.bright_mass_sol`).

## The rule

* Ring 1 is every sector a face away from a generated sector (the face neighbours of
  `galaxy.geometry.neighbor_addresses`; no diagonals). Ring 2 is a face away from ring 1, and so on
  to ring 4, each ring counted from the previous ring's outer edge. A cell outside the galaxy's
  outline is left out.
* Ring 1 gets every star born with at least **1** solar mass, ring 2 at least **2**, ring 3 at least
  **5**, ring 4 at least **8** (`tuning.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL`). The nearest ring has the
  lowest cut. Past ring 4 the scatter's own mass limit (GEN.183) stands, and a ring whose mass is not
  below that limit has nothing to add.
* The generated sectors are the run's own (`backfill_after_run`: the sectors created since the run
  began) or the one sector an on-demand visit or `generate_sector_neighborhood` creates. A sector
  already filled is never touched, and one a later run reaches nearer is taken lower and given only
  the bands it lacks.
* The GEN.30 luminosity tiers (100/250/500/750 L_sun out to 10/25/50/100 ly) are gone.

## What a sector holds

The galaxy scatter places two kinds of star: every star born at or above the mass limit M (the mass
pass, any luminosity), and, of the lighter ones, every star at or above the luminosity floor L (the
luminosity pass). A mass backfill to m (m below M) adds the rest of the stars born in [m, M): the
living ones **dimmer than L**, since the brighter ones were placed already. A star the mass pass
or the backfill leaves out is the sector's own draw at fill time, which now caps initial mass at
m as well as luminosity at L (`FillContext.star_mass_limit_sol`, `store.bright_star_fill_mass_limit`).

White dwarfs born at 1 solar mass or more count: they are living stars of that mass, and without
them the sector's fill (which refuses initial masses at or above m) would lose them. A white dwarf
row stores `NULL` for its infinite lifespan and phase end, and `bright_stars.star_params` reads
`NULL` back as infinity for class VII.

At a galaxy with no scatter at all there is no L and no M: the backfill places every star born at or
above the ring's mass, at any luminosity.

## Per sector state

`sector_stats.bright_mass_sol` (new, `NULL` = never backfilled by mass) records how far down a
sector was taken. When a scatter exists the same transaction also sets `bright_level_sol` to the
scatter's L, so a staged scatter (`planetgen plan --bright-stars-down-to`) tops the sector up as it
does any sector with its own state: from L down to the new floor, among stars born lighter than the
sector's mass (`run_plan._draw_sector_bands` takes `min(mass limit, bright_mass_sol)`). Deleting a
filled sector puts its level back; the mass stays. `clear_bright_stars` forgets both.

## Reproducibility

Stars are drawn in fixed mass bands between the edges 1, 2, 5, 8 and the galaxy's mass limit
(`bright_stars.canonical_mass_bands`). Each band of each sector is drawn whole from its own stream,
keyed by the galaxy seed, the sector's address and the band's lower mass, so a sector taken down
in two steps ends with exactly the stars one draw down to the same mass gives, whatever order the
runs happened in. The band's Poisson mean is the sector's expected count times the band's share of
living stars (`star_population.mass_band_fraction`, the same denominator as `bright_star_fraction`, so
mass and luminosity shares subtract).

## How big

The share of living stars born at or above each mass, of the whole disk (computed from the
population model, `massive_star_fraction`): 1 Msun 9.1%, 2 Msun 3.3%, 5 Msun 0.55%, 8 Msun 0.001% (the
young disk alone: 9.3, 3.5, 0.70, 0.15%). Ring 1 therefore adds about nine stars in a hundred of what a
sector will hold, ring 2 three, ring 3 half a star, so the work is a few percent of filling the
sectors in place and grows with the surface of the run, not its volume. The luminosity shares that
remain after the scatter's L are nearly all of those (under L = 5000, 1 to 2 Msun: 5.8% of living
stars, 2 to 5: 2.8%, 5 to 8: 0.55%).

## Code

`bright_stars.backfill_mass_cells`, `star_population.mass_band_fraction` and `sample_mass_band_stars`;
`run_galaxy.backfill_ring_targets` (the ring walk) and `backfill_bright_stars_around`;
`run_plan._draw_sector_masses` (chunked under `lock_sector_stats`, like the luminosity bands);
`store.sector_bright_masses`, `set_sector_bright_masses`, `bright_star_fill_mass_limit`.
