# Galaxy coordinate system and sector grid

**Status:** built. The sector grid is spherical: shells, bands and wedges.
Boss decided the move on 2026-10-15 06:39Z; the plan is
[`spherical-sector-coordinates.md`](spherical-sector-coordinates.md) and the
study behind it is [`spherical-shell-sectors.md`](spherical-shell-sectors.md).
The grid went in with schema v91 (migration `0091_spherical_sectors.py`).
Code: `src/planetgen/galaxy/spherical.py` (the grid), `skeleton.py` (the
outline), `density.py` (the density model, see `galaxy-disk-density.md`) and
`planetgen plan` / `planetgen galaxy`. The Galaxy Map's block ladder, stages
and windows sit on the same grid (`galaxy/drill.py`, `galaxy/window.py`,
`spherical-block-ladder.md`; MAP.205, MAP.206).

The cylindrical ring, layer and slot grid this replaced is kept in
[`archive/galaxy-coordinate-system-cylinder.md`](archive/galaxy-coordinate-system-cylinder.md),
and the older spherical-shell design before it in
[`archive/galaxy-coordinate-system-shells.md`](archive/galaxy-coordinate-system-shells.md).
Words matter here: **shell**, **band** and **wedge** are the only names for
the three coordinates. "Ring", "layer", "slot" and "slab" name nothing in the
current grid and are not reused.

## 1. The galaxy frame

- Origin at the galactic center, `(0, 0, 0)`. The quasar, when a galaxy has
  one, sits in the core, shell 0.
- Right-handed axes: `X` and `Y` in the galactic plane, `+Z` toward galactic
  north. `+X` is the zero meridian (wedges are counted counterclockwise from
  it). This is a generation-space convention and does not match any real sky
  frame.
- Galaxy-scale distances are stored in parsecs (`_pc` columns); positions
  inside a sector in milliparsecs (`_mpc`); bodies inside a system in km
  (shown in AU). `planetgen/physics/units.py` has `mpc_to_pc`, `pc_to_mpc`,
  `pc_to_ly` and `ly_to_pc`.

## 2. The spherical sector grid

Every galaxy-placed sector is one cell of a grid of nested spherical shells
around the galactic origin. With edge `e` (4 pc, section 3), a cell is named
by three integers, `(shell, band, wedge)`, stored as `sectors.shell_index`,
`sectors.band_index` and `sectors.wedge_index`:

- **Shell `k`** (radius): `r` in `[(k - 1/2) * e, (k + 1/2) * e)`, so its
  centerline radius is `k * e`. A shell is always `e` deep; that is the one
  length the grid never adjusts. **Shell 0 is the core**, the ball of radius
  `e / 2` (2 pc), a single cell.
- **Band** (latitude): a ring of cells around a shell between two colatitudes.
  Each shell's two polar caps are cut off first, and the rest is cut into an
  odd number of bands of equal colatitude step, each about `e` of arc, so the
  galactic plane runs through the middle of the equatorial band. Bands are
  signed: band 0 is the equatorial band, positive bands lie north (`z > 0`),
  negative bands south.
- **Wedge** (longitude): a band is cut by meridians, pole to pole, into
  `wedge_count` equal wedges. Wedge `w` spans longitude `[2*pi*w/W,
  2*pi*(w+1)/W)` from `+X` counterclockwise.
- **Cap:** the circular segment of a shell within colatitude `alpha` of a
  pole, one cell at each pole. Six wedges meet its rim with a gap of exactly
  one edge, so the rim is `6 * e` around and the cap's chord radius is
  `6 * e / (2 * pi) = 3.82 pc` (4 pc would be 6.28 wedges, and 6 is the
  nearest whole number); `alpha = asin(rho / R)`. A cap is band `+/-(h + 1)`
  with wedge 0, where `h` is the shell's highest ordinary band (the shell
  holds `2 * h + 1` ordinary bands). A cap holds about 2.9 edge cubes.

**One wedge rule everywhere (Boss, 2026-10-15 07:17Z).** Every band of every
shell has `W = 6 * m` wedges, where `m` is the number whose prime factors are
only 2, 3, 5 and 7 that is nearest (in ratio) to
`2 * pi * R * sin(phi) / (6 * e)`, `phi` being the band's centre colatitude,
and `m` is at least 1 (`spherical.wedge_count`). The wedge's arc at the
band's centre therefore stays near one edge, and the six meridians that cut
every cap's rim are wedge boundaries in every band. The rule never looks at
the outline, the window or the neighbouring shell, so every lookup stays
closed form. Measured on shell 3,000 (plan section 2): arcs 0.981 to 1.019 of
4 pc, and 99.8% of all cells within 0.5% of an edge cube's volume; near the
poles the small counts let arcs stray further (about 0.76 to 1.06 at shell
10). A shell's band count is the odd integer nearest
`R * (pi - 2 * alpha) / e` (`spherical.shell_layout`).

**Shell 1 is cut by hand** (`HAND_CUT_SHELLS`). The rule would leave it caps
of 4.8 edge cubes against 0.67 for the cells between, so its caps end at 37.5
degrees and it holds one band of 12 wedges (14 cells, volumes about 0.9 to
1.4 edge cubes). Shells 2 and 3 follow the rule like every other shell.

**No undefined space (Boss, 2026-10-15 06:59Z).** The core ball, every shell,
every cap, band and wedge together tile all of space: every point has exactly
one address, with no gap and no overlap, whether or not that sector was ever
generated. The outline (section 4) only says which cells may be generated; it
never decides which cells exist. `sector_address_at` answers for any point,
inside or outside the outline, as far as the 40-bit birth-sector field
reaches (shell below 4,096, about 53,000 ly; the default galaxy's outline
ends at shell 3,764).

**Lookups are closed form.** `sector_address_at` maps a point to its cell
(shell from `|p|`, band from `acos(z / |p|)`, wedge from `atan2(y, x)`);
`sector_position_pc` gives a cell's center (a cap's is on the axis at
`+/-shell * e`, the core's is the origin); `cell_volume`,
`sector_cell_vertices_pc` (8 corners) and `sector_cell_test` (a predicate for
one cell) are computed. Only the address and center are stored
(`uq_sectors_address` makes the address unique). `cell_count`,
`cells_in_band_range` and `nth_cell_in_band_range` count and index the cells
of a shell without listing them.

**Neighbors.** `neighbor_addresses` lists the cells sharing a face of nonzero
area: the two wedges either side in the same band; the cells above and below
in the shell (one or several, because neighbouring bands hold different
wedge counts and are offset like brick courses); and, in the shell inside and
the shell outside, every cell whose band overlaps this one's colatitudes and
whose wedge overlaps its longitudes. A cap touches every wedge of the nearest
band, and the core touches the whole of shell 1. A cell has about 7 face
neighbours, against 6 for a cube, up to roughly 10 to 15 near a cap or a
change in wedge count. Flat height slabs and aligned vertical columns do not
exist any more: a region is an arc per (shell, band), not per ring.

**Local frame.** A sector's local axes (`sector_orientation`) are `+X`
radially outward from the galactic axis, `+Y` toward increasing longitude and
`+Z` galactic north, unchanged in meaning from before. System and phenomenon
offsets (`star_systems.position_x/y/z_mpc` and the phenomenon equivalents) use
this frame, and generation samples them inside the real cell (`SectorCell`:
cube root for the radius, uniform in the cosine of the colatitude, uniform in
longitude; a cap samples a disc), not a cube. The octant label
(`star_systems.quadrant`, and on every placed phenomenon) comes from this
frame.

**Designations.** `provisional_sector_designation` packs the address into the
40-bit birth-sector field of the object ID: `shell << 28 | (band + 2048) << 16
| wedge`, printed as uppercase hex; `parse_sector_designation` reverses it. It
allows shells 0 to 4,095, bands -2048 to 2047 and 65,536 wedges a band (the
default galaxy's widest band holds about 23,650). The core is `8000000`, and
for example `(12, 3, 5)` is `C8030005`. A sector's ID is its designation (see
[object-ids.md](object-ids.md)), the first group of every object ID born in
the sector. Object IDs were renumbered from the new sectors when the grid
changed (section 6).

**Seeds.** A unit of work is seeded from the galaxy seed and its address
(`galaxy/seed.py`): a sector's seed string is `sector:<shell>/<band>/<wedge>`,
for example `sector:12/3/0`, hashed with the galaxy seed. A sector's numbers
therefore depend only on the seed and its address, never on the order sectors
run in.

**Zones.** The Galaxy pages group sectors in Zones of `ZONE_SHELL_DEPTH =
round(100 ly / 13.05 ly) = 8` consecutive shells
(`planetgen/web/maps/galaxymap.py`, `spherical.sector_zone`).

## 3. Sector size: 4 pc

`tuning.DEFAULT_SECTOR_EDGE_PC = 4` (13.05 ly) is the shell depth and the
arc every band and wedge aims at, used by the generator, the skeleton, the
Galaxy pages and the Galaxy Map. `planetgen plan` has no edge option. It
replaced 11.5 ly (3.53 pc) in 7.0.0.

The edge sets where the galaxy ends, because a sector exists only where it
expects at least one star (`predicted_star_count >= 1`), and a bigger cell
holds more stars at the same density. At the default Milky Way shape a 4 pc
sector holds about 6.3 stars at local density, the galaxy reaches about
50,300 ly (shell 3,764) and the grid holds about 23.9 billion candidate cells
(plan section 2, measured on the real density model). 2 pc was rejected
because a local-density sector would average under one star, so ordinary
solar-neighborhood space would not qualify, and 11.5 ly was the older,
arbitrary default.

## 4. The outline: one band range per shell

Shell `k` can hold content in bands `-high .. high` around the plane (the two
caps too when `high` is the shell's top band). `planetgen plan` walks the
shells outward from the core and stops at the first shell whose equatorial
band can no longer reach the qualification threshold (a sector expected to
hold at least one star; the exponential density model never reaches zero, so
the threshold is the edge). Each shell's run is found by bisection on an
upper bound of the density over longitude
(`skeleton.build_shell_extents`, the bound is `bound_relative_density_at`:
the bar bulge at its long axis plus the disks times the largest arm factor),
with the cap checked first, because the bound falls moving north from the
plane until the last step toward the pole, where it can rise again.

`planetgen plan` stores:

- `galaxy_shape`: one row with the shape parameters, `edge_pc`,
  `expected_system_count_at_density_1`, `outer_shell_index`, the bright-star
  threshold and seed, and the phenomenon scatter seed.
- `galaxy_shell`: `(shell_index, high_band_index, cell_count)`, one row per
  shell that can hold content. A shell with no row has none.

The build takes a few milliseconds. The outline is a safe superset: a sector
inside its shell's range may still fall short (an inter-arm trough), and that
exact check is one density call when the sector is visited. A band is never
left partly in; the outline decides whole bands and the exact check decides
sectors.

**Validation before generation.** Every `planetgen galaxy` mode (an address, a
span of shells, bands and wedges, a neighborhood, a random start, a course)
and visit-time generation (`generate.ensure_sector_generated`) ask
`skeleton.GalaxyBounds.contains(shell, band)` first. An address outside is
refused with the reason (`describe_miss`), even with `--density` or
`--num-systems`. A neighborhood near the edge leaves out the sectors past it.
A random start draws a uniformly random sector from inside the outline, every
sector equally likely (`GalaxyBounds.random_address`). `planetgen galaxy`
refuses to run before `planetgen plan`.

**Bright stars (schema v43).** `planetgen plan` also places every star at or
above `BRIGHT_STAR_MIN_LUMINOSITY_SOL` at a fixed point in its sector, in
`bright_stars`, before any sector is filled (the default is 1,000 L_sun, with
the dimmer stars filled in near generated sectors by the mass backfill of
GEN.187, `sector_stats.bright_mass_sol`). Filling the sector later builds a
full system around each of them. The Galaxy Map draws them, so the arms show
before any sector is filled. How much of a planet's naked-eye sky these stars
and the backfill cover is in [sky-view.md](sky-view.md) section 3.

**Adding a system to a stored sector.** `POST /api/systems` with a
`sector_id` places a new system in an existing sector
(`store.add_system_to_sector`): at a given position inside it, or else clear
of the other systems' Hill spheres, with its location, containment and
nearest systems filled in like a generated one.

## 5. Things that move

- **Hill spheres follow the real radius.** `Star.calculate_system_perimeter`
  takes `galactic_center_dist_ly`, threaded from the owning sector's
  `galactic_radius_pc`. The fixed `GALACTIC_CENTER_DISTANCE_LY` (Sol's
  distance) is used only for a system generated outside the galaxy.
- **Galactic orbits.** `planetgen.cli.orbits` turns every system, phenomenon
  and stand-alone facility along its galactic orbit, moves anything that
  drifts into another generated sector over to it, then recomputes
  containment and the stored nearest systems (`nearest_systems`).
  Which sector holds a moved point is `spherical.sector_address_at`; the
  object keeps its ID wherever it drifts (the ID names the birth sector).

## 6. Migrations that cleared or converted sectors

Earlier grid changes deleted what no longer had a cell. `_migrate_v31_to_v32`
(cylinder) deleted every galaxy-placed sector with its systems and
phenomena; the v33 edge change did the same and rebuilt the outline at 4 pc;
the v35 slot rule deleted the sectors in every ring whose slot count changed.
Sectors never placed in the galaxy were kept each time.

The move to shells (v91) deleted nothing. Migration 0091 renames the address
columns and converts the data: every star system, phenomenon and facility
keeps its place in the galaxy, its velocity and its epoch; the sector rows are
rebuilt (one for every cell that holds an object and one for every cell whose
center lay inside a generated sector, each new sector taking the name of the
old sector holding its center); object IDs are renumbered for the new birth
sector; the bright stars and the phenomenon scatter get the address of the
cell holding their position; `sector_stats` is counted again; and
`galaxy_shell` is built again from the stored shape. The details are in
[`../database-schema.md`](../database-schema.md) under v91.

## 7. Why it works this way

- **Spherical shells again, with a round core (Boss, 2026-10-15).** The
  cylinder gave cells that were close to cubes, but it put the galactic center
  on an axis rather than at a point, made the core a column, and tied every
  region and block to flat height slabs. Boss chose "two meridian lines from
  pole to pole make a wedge (vert), shell (radius), and band (horizontal),
  together they identify a single 3D segment of the sphere", with the core a
  ball 4 pc across, the shell depth always 4 pc, and bands and wedges kept as
  close to 4 pc arcs as possible, especially at the poles. The study
  recommended keeping the cylinder; Boss decided otherwise, and the plan
  records the costs it accepted: more face neighbours, no flat slabs or
  aligned columns, and heavier cells at the poles near the core (a cap there
  holds about 2,800 systems against 1,114 for the cylinder's densest cell).
- **Caps with six rim wedges.** Six is the whole number nearest the 6.28
  wedges a 4 pc gap would give. Five (3.18 pc) and seven (4.46 pc) were
  measured and left cap volumes of 2.0 and 3.9 edge cubes against 2.9.
- **Wedges in multiples of 6 with primes up to 7 (Boss, 07:17Z).** Any wedge
  count per band gave cell volumes of 0.94 to 1.09 of a cube but boundaries
  that shift between bands; doubling 6, 12, 24 left only 23% of cells within
  10% of a cube. Multiples of 6 whose cofactor uses only 2, 3, 5 and 7 keep
  the six pole meridians as boundaries in every band and neighbouring bands
  mostly nested (95% on shell 3,000), with the same rule everywhere.
- **Equal-volume bands were not adopted (Boss, 07:22Z).**
- **The plane through the middle of a band.** An odd number of bands puts the
  galactic plane through the equatorial band's center, so band 0 is
  symmetric and the outline is symmetric about the plane.
- **Object IDs renumbered (Boss, 07:32Z).** Every object got a new ID from
  its new sector in the conversion; positions are unchanged.
- **Validate before generating.** Earlier random starts drew from a fixed
  region that could land outside the galaxy; checking afterwards wasted
  generation. Storing the outline lets every path refuse a bad address up
  front.
- **Parsecs for galaxy scale.** Milliparsecs give 8-digit numbers across the
  galaxy, and light-years would need an AU round trip to compare with the
  `_mpc` sector columns. Parsecs stay in the same unit family as
  milliparsecs.
- **Store the address and center only.** Position, density and outline are
  pure functions of the address, so storing them would only copy a formula's
  output. `galactic_radius_pc` is stored anyway so "sectors within R of the
  core" can use an index.
- **Lazy generation.** At about 23.9 billion candidate cells, generating
  everything up front is not possible. Sectors are generated when visited or
  asked for.

### Alternatives considered and rejected

- The cylindrical ring, layer and slot grid (built 6.0.0 to 8.0, replaced by
  v91; see the archive copy).
- The first spherical-shell design with Fibonacci slots and Voronoi prisms
  (built in 5.3.6, replaced in 6.0.0; see the older archive copy).
- A per-sector plan table of every qualifying sector, about 10.5 billion
  rows and about 1 TB (see `archive/galaxy-disk-density-rev2.md`).
- A k-d tree or octree over sector addresses: not needed, since every lookup
  is closed form.
- 2 pc, 3 pc and 5 pc edges (section 3).

## Open questions still standing

- **The Galaxy Map's block ladder on shells.** The old ladder (243, 81, 27,
  9, 3 and 1 sectors a side) cut blocks on shared master lines in longitude;
  shells need a latitude rule on top, and a first try with odd powers of 3 as
  band masters gave band arcs of 0.89 to 1.11 of 4 pc. This is MAP.205 and
  MAP.206, designed with Maps.
- **One galaxy per database.** There is no `galaxy_id`. Several galaxies in
  one database would need one on every galaxy table. The plan for other
  galaxies (GEN.9) recommends keeping it that way: a second galaxy is a
  second database, and the neighbours are rows of one new table. See
  [multiple-galaxies.md](multiple-galaxies.md).
