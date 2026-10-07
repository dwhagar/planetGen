# Galaxy coordinate system and sector grid

**Status:** built. The cylindrical grid arrived in 6.0.0 (schema v32), the
4 pc edge and aligned layers in 7.0.0 (schema v33), and the hybrid
master-wedge slot counts in 7.13.0 (schema v35). Code:
`src/planetgen/galaxy/geometry.py` (the grid), `galaxySkeleton.py` (the
outline), `galaxyDensity.py` (the density model, see
`galaxy-disk-density.md`) and `generate.py plan` / `generate.py galaxy`.

The spherical-shell design this replaced (shells, Fibonacci slots, Voronoi
prisms, the `sector_vertices` and `galaxy_shell_band` tables) is kept in
`archive/galaxy-coordinate-system-shells.md`.

## 1. The galaxy frame

- Origin at the galactic center, `(0, 0, 0)`. The quasar, when a galaxy has
  one, sits in ring 0, layer 0, slot 0.
- Right-handed axes: `X` and `Y` in the galactic plane, `+Z` toward galactic
  north. `+X` is the zero meridian (slot boundaries are counted
  counterclockwise from it). This is a generation-space convention and does
  not match any real sky frame.
- Galaxy-scale distances are stored in parsecs (`_pc` columns); positions
  inside a sector in milliparsecs (`_mpc`); bodies inside a system in km
  (shown in AU). `planetgen/physics/units.py` has `mpc_to_pc`, `pc_to_mpc`,
  `pc_to_ly` and `ly_to_pc`.

## 2. The cylindrical sector grid

Every galaxy-placed sector is one cell of a cylindrical grid around the
galactic axis. With edge `e`:

- **Ring `i`** (`sectors.ring_index`, `i >= 0`): cylindrical radius `R` in
  `[i*e, (i+1)*e)`. A sector's center sits on the ring's centerline,
  `(i + 1/2) * e`.
- **Layer `j`** (`sectors.layer_index`, any integer): height `z` in
  `[(j - 1/2)*e, (j + 1/2)*e)`, so layer 0 is centered on the plane.
- **Slot `k`** (`sectors.ring_slot_index`): one of `N_i` equal wedges,
  counterclockwise from `+X`.

`N_i` follows the hybrid master-wedge rule (7.13.0, schema v35, Boss's
choice of 2026-09-30):

- `c_i = 2*pi*(i + 1/2)` is the ring's centerline circumference in edges.
- `M_i` (`ring_master_count`) is 3 at the center and doubles (6, 12, ...,
  1,536 at the default galaxy's edge) at the first ring where
  `c_i >= 2 * M * 8`, that is where each doubled wedge would still hold 8
  slots. The doublings happen at rings 8, 15, 31, 61, 122, 244, 489, 978
  and 1,956.
- `N_i = M_i * round(c_i / M_i)`, never less than `M_i`: 3, 9, 15, 21, 27,
  36, 42, 48, 54, ...

Each slot's centerline arc stays within 0.94 to 1.06 of `e`, and the total
count is within 0.1% of plain rounding. Because every `N_i` is a multiple of
`M_i` and `M_i` only grows outward, every master line is a slot boundary in
every ring from where it starts out to the edge. The Galaxy Map cuts its
wedge lines and blocks on these lines. `static/galaxyprisms.js` mirrors
`ring_sector_count` and `ring_master_count`.

**Aligned columns.** `N_i` depends only on the ring, so every layer cuts
ring `i` the same way: cell `(i, j, k)` sits directly above `(i, j-1, k)`.

**Exact tiling.** The cells tile space with no gaps and no overlaps, and
every lookup is closed form: `sector_address_at` maps a point to its cell,
`neighbor_addresses` lists the face neighbors (two slots along the ring, the
cells above and below, and the one or two overlapping slots in each
neighboring ring), and `enumerate_sectors_within_radius` walks only the
cells near a point. A cell's center, corners and volume are computed, so
only the address and center are stored (`uq_sectors_address` makes the
address unique).

**Local frame.** A sector's local axes (`sector_orientation`) are `+X`
radially outward from the galactic axis, `+Y` toward increasing angle and
`+Z` galactic north. System and phenomenon offsets
(`star_systems.position_x/y/z_mpc` and the phenomenon equivalents) use this
frame, and generation samples them inside the real cell (`SectorCell`), not
a cube. The octant label (`star_systems.quadrant`, and since schema v41 the
same label on every placed phenomenon) comes from this frame.

**Designations.** `provisional_sector_designation` packs
`ring << 33 | (layer + 4096) << 20 | slot` and prints it as hex;
`parse_sector_designation` reverses it. Layers run from -4096 to 4095.

**Zones.** The Galaxy pages group sectors in Zones of
`ZONE_RING_WIDTH = round(100 ly / 13.05 ly) = 8` rings
(`planetgen/web/maps/galaxymap.py`).

## 3. Sector size: 4 pc

`program_constants.DEFAULT_SECTOR_EDGE_PC = 4` (13.05 ly), used by the
generator, the skeleton, the Galaxy pages and the Galaxy Map. `generate.py
plan` has no edge option. It replaced 11.5 ly (3.53 pc) in 7.0.0.

The edge sets where the galaxy ends, because a sector exists only where it
expects at least one star (`predicted_star_count >= 1`), and a bigger cell
holds more stars at the same density. At the default Milky Way shape:

| Edge | Stars per sector at local density | Radius reached | Layers | Candidate sectors |
|---|---|---|---|---|
| 11.5 ly (old) | 4.3 | 46,900 ly | 681 | 14.9 billion |
| 2 pc | 0.8 | 31,300 ly | 893 | 28.2 billion |
| 3 pc | 2.7 | 42,400 ly | 743 | 18.6 billion |
| **4 pc** | **6.3** | **50,300 ly** | **635** | **12.3 billion** |
| 5 pc | 12.3 | 56,400 ly | 557 | 8.6 billion |

## 4. The outline: layers and columns

The layers run from `+317` to `-317` at the default shape. Each layer runs
from ring 0 out to the last ring that could still expect a star. A sector
center's highest possible density over angle is
`bulge(r_3d) + disk(R) * f_z(z) * (1 + arm_amplitude)`, which falls as `R`
and `|z|` grow, so each layer's rings are one unbroken run from ring 0 and
each layer reaches no farther than the one below it (toward the plane).

`generate.py plan` stores:

- `galaxy_shape`: one row with the shape parameters, `edge_pc`,
  `expected_system_count_at_density_1`, `outer_ring_index`, and (schema v43)
  the bright-star threshold and seed.
- `galaxy_layer`: `(layer_index, outer_ring_index)`, one row per layer
  (`galaxySkeleton.build_layer_extents`).
- `galaxy_column`: `(ring_index, layer_index_min, layer_index_max)`, one row
  per ring.

The build takes a few milliseconds. A sector inside its layer's extent may
still fall short (an inter-arm trough); that exact check is one density call
when the sector is visited.

**Validation before generation.** Every `generate.py galaxy` mode (a ring, a
slot, a neighborhood, a random start) and visit-time generation
(`generate.ensure_sector_generated`) ask `galaxySkeleton.GalaxyBounds.contains`
first. An address outside is refused with the reason, even with `--density`
or `--num-systems`. A neighborhood near the edge leaves out the sectors past
it. A random start draws a uniformly random sector from inside the outline
(`GalaxyBounds.random_address`). `generate.py galaxy` refuses to run before
`generate.py plan`.

**Bright stars (schema v43, 7.38.0).** `generate.py plan` also places every
star at or above `BRIGHT_STAR_MIN_LUMINOSITY_SOL` (500 by default) at a
fixed point in its sector, in `bright_stars`, before any sector is filled.
Filling the sector later builds a full system around each of them. Since
7.40.1, `generate.py plan --bright-star-min-luminosity` accepts down to
100 (about 220 million stars and 35 GB in a Milky Way, against about 60
million and 10 GB at 500); white dwarfs are never pre-placed. The Galaxy
Map draws them from 7.42.0, so the arms show before any sector is filled.

**Adding a system to a stored sector (7.43.0).** `POST /api/systems` with
a `sector_id` places a new system in an existing sector
(`_db.add_system_to_sector`): at a given position inside it, or else
clear of the other systems' Hill spheres, with its location, containment and
nearest systems filled in like a generated one.

## 5. Things that move

- **Hill spheres follow the real radius.** `Star.calculate_system_perimeter`
  takes `galactic_center_dist_ly`, threaded from the owning sector's
  `galactic_radius_pc` (5.3.6). The fixed `GALACTIC_CENTER_DISTANCE_LY`
  (Sol's distance) is used only for a system generated outside the galaxy.
- **Galactic orbits.** Since 7.37.0 `planetgen.cli.orbits` turns every system,
  phenomenon and stand-alone facility along its galactic orbit, moves
  anything that drifts into another generated sector over to it, then
  recomputes containment and the stored nearest systems (`nearest_systems`,
  schema v41).

## 6. Migrations that cleared sectors

Shell addresses had no matching cell, so `_migrate_v31_to_v32` deleted every
galaxy-placed sector with its systems and phenomena. The v33 edge change did
the same and rebuilt the skeleton at 4 pc. The v35 slot rule deleted the
sectors in every ring whose slot count changed (all but 15 rings). Sectors
never placed in the galaxy were kept each time.

## 7. Why it works this way

- **Cylinders instead of spherical shells (6.0.0).** The first design
  (2026-09, "Track C") placed cube sectors on concentric spherical shells
  with Fibonacci-sphere slots. A cube cannot tile a sphere, so it needed a
  per-sector Voronoi prism (`sector_vertices`), a Fibonacci-number neighbor
  search and a per-shell band cache (`galaxy_shell_band`), and most shell
  slots sat far off the disk. The release note gives the reason for the
  switch: sectors now "follow the flat disk instead of a ball". The grid
  also tiles exactly, makes every lookup closed form, and let the
  `sector_vertices` and `galaxy_shell_band` tables be dropped (CHANGELOG
  6.0.0).
- **Rings as wide as layers are tall.** One edge sets ring width, layer
  height and slot arc, so every cell is close to an `edge_pc` cube
  (`galaxyGeometry` module docstring), like the cube sectors before it. No
  further reason is recorded.
- **Aligned layers (7.0.0).** The 6.0.0 grid used a multiple of 4 slots per
  ring; 7.0.0 made the count `round(2*pi*(i + 1/2))`, the same on every
  layer, so sectors line up in vertical columns and the skeleton is one
  row per layer (CHANGELOG 7.0.0, commit d9ed2b7).
- **4 pc edge (7.0.0).** A whole number of parsecs, and at the default shape
  it puts the galaxy's edge at about 50,000 ly, the real Milky Way disk's
  radius. 2 pc was rejected because a local-density sector would average
  under one star, so ordinary solar-neighborhood space would not qualify.
  11.5 ly was the older, arbitrary default.
- **Master wedges (7.13.0, MAP.33).** Plain rounding gave slot boundaries
  that did not line up from ring to ring, so the Galaxy Map's large blocks
  could not be cut on shared lines. Multiples of a doubling master count make
  every master line a slot boundary out to the edge, while arcs stay
  within about 6% of an edge and the total count is unchanged. Boss chose this hybrid on
  2026-09-30.
- **Validate before generating (7.0.0, commit c4b4e6e).** Earlier random
  starts drew from a fixed 15,000 pc cylinder 2,000 pc tall and could land
  outside the galaxy; checking afterwards wasted generation. Storing the
  outline lets every path refuse a bad address up front.
- **Parsecs for galaxy scale.** Milliparsecs give 8-digit numbers across the
  galaxy, and light-years would need an AU round trip to compare with the
  `_mpc` sector columns. Parsecs stay in the same unit family as
  milliparsecs (original design, section 2).
- **Store the address and center only.** Position, density and outline are
  pure functions of the address, so storing them would only copy a formula's
  output (original design, section 10). `galactic_radius_pc` is stored
  anyway so "sectors within R of the core" can use an index.
- **Lazy generation.** At about 12.3 billion candidate cells, generating
  everything up front is not possible. Sectors are generated when visited or
  asked for.

### Alternatives considered and rejected

- Spherical shells with Fibonacci slots and Voronoi prisms (built in 5.3.6,
  replaced in 6.0.0; see the archive copy).
- A per-sector plan table of every qualifying sector, about 10.5 billion
  rows and about 1 TB (see `archive/galaxy-disk-density-rev2.md`).
- A k-d tree or octree over sector addresses: not needed, since every lookup
  is closed form.
- 2 pc, 3 pc and 5 pc edges (table in section 3).
- A staggered "brick" layer pattern: rejected for direct vertical columns.
  The reason recorded is only that columns line up; nothing more is written
  down.

## Open questions still standing

- **One galaxy per database.** There is no `galaxy_id`. Several galaxies in
  one database would need one on every galaxy table.
