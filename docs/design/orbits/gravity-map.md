# Gravity map: the field of a sector as a heat map

**Status (2026-10-15):** decided; not built. GEN.198, GEN.199, PERF.72, TEST.131 and MAP.167 to MAP.170 are open. This is the third phase of 9.0 (Boss folded 9.1 into 9.0 on 2026-11-05).

Boss (2026-10-10 20:57Z): the orbits are worked out from point-mass
vectors, so show "a gravitational map of a sector, probably using colors
and calculating it per zone ... the gravitational gradients of the vector
field as kind of a heat map of the inside of the sector". Foundations are
phase 2; the feature roll-out is phase 3 and phase 3 of 9.0 (9.1 was folded into 9.0 on 2026-11-05; its foundations are built in 8.2, 8.3 and 8.5).

## What is drawn

The field has three useful scalars at every point, and the map offers each
as a mode (one at a time, switchable in the map Menu):

| Mode | Quantity | Unit | Reads as |
|---|---|---|---|
| Pull | the size of the acceleration vector, \|g\| | m/s^2 | where gravity is strong |
| Well depth | the potential, Phi | m^2/s^2 | how deep the wells are |
| Tidal strength | the size of the gradient of g (the tidal tensor's largest eigenvalue) | 1/s^2 | where space is stretched; Lagrange points show as saddles |

The gradient mode is what Boss called the gradients of the vector field.
Lagrange and saddle points (where g is nearly zero but the gradient is not)
are marked.

## How it is computed

- **Zones.** A sector (edge 4 pc, `DEFAULT_SECTOR_EDGE_PC`) is cut into
  equal zones; the default is 16 per edge (0.25 pc), a setting. The field
  is sampled at each zone's centre.
- **Sources.** Every point mass in the sector and its neighbours: **one per
  star system** (its stars, planets and moons together at the system's place,
  which changes the field at a zone's size by under 0.1 percent; MAP.386,
  `db/gravity_sources.py`), stand-alone black holes and neutron stars (the
  same table the orbital update reads); rogue planets, comets and asteroid
  fields are below the tolerance and left out. The System Map's orbital-plane
  version reads the bodies themselves (MAP.169). The galaxy's own smooth
  potential (GEN.115) is the background, joined by the rule in "The smooth
  galaxy and the generated stars" below (MAP.383). Plummer softening keeps a
  zone that holds a source finite.
- **Near and far.** Sources close to a zone are summed exactly; distant
  ones are summed as one point per cell (mass and centre of mass) of a
  hierarchy with 0.5 pc leaf cells, the way a tree code does
  (MAP.385, `galaxy/gravity_tree.py`; the MAP.151 region is the data window and
  holds no masses). The cost study (PERF.72) fixes the split, the opening
  angle (0.3 to 0.7, chosen to a pair budget) and checks it against the exact
  sum; the kernel is `galaxy/gravity_kernel.py` (MAP.384).
- **Cache.** The grid is derived data, cached like the Galaxy Map tiles and
  thrown away when the orbital update (GEN.105) moves the sources; it is
  never part of the galaxy itself, so it needs no galaxy-table migration.
- **Slices and volumes.** The map shows one slab (a plane slice) at a time
  and can show the whole volume as transparent layers; on the System Map it
  is the orbital plane.

## The smooth galaxy and the generated stars (MAP.383)

The smooth potential (`galaxy/potential.py`, GEN.115) already holds the mass of the stars the generator draws: its disks, bulge and halo are fitted to the
same density model. Adding the generated stars on top of it counts that mass twice. Measured on generated blocks (`reports/map383-background-rule-2026-11-05.md`),
the double count is small in the pull (0.02 percent of the smooth pull at the solar circle, 0.3 percent in the bulge) and in the potential, and
**large in the tidal strength: the tidal norm is off by about 36 percent at the median and 67 to 80 percent at the 95th percentile**, because the smooth tidal
trace is the mean stellar density (`4 pi G rho`) and a sum of point masses has none between the stars.

The rule (points 1 to 4 and 6 are MAP.383's; point 5 is a recommendation MAP.167 confirms):

1. **Field = smooth galaxy + generated stars - the smooth galaxy's own stars in the same sectors.** The last term is each generated sector's share of the
   smooth model's stellar mass, spread through the sector as a cube and softened by the zone edge, and summed by the same tree as the generated stars.
   Equivalent: the smooth field plus the *fluctuation* of the generated stars about their mean.
2. **The mass subtracted is the smooth potential's own stellar density** (the Poisson density of its disk and sphere terms, closed forms), not the density
   model's expected mass. The two differ by 15 to 19 percent where measured (solar circle: potential 0.026 against model 0.031 Msun per pc^3; bulge: 3.68
   against 3.10) because the potential's disks and spheres are stand-ins (GP1 to GP3). Using the potential's own mass changes the tidal norm by 5 to 7
   percent against subtracting the model's, and removes that mismatch from the map instead of showing it as stars that are not there.
3. **A sector that is not generated adds nothing and subtracts nothing** (its generated and subtracted mass would be equal), so the near region needs no
   expected-mass cells for holes.
4. **Dark halo, gas, black hole and every other smooth term stay in the background untouched.**
5. **The map has two readings of each quantity.** The *total* field (the rule above) and the *local* field (generated stars minus their mean). In the pull
   the local field is 0.03 to 0.09 percent of the total, so a total-pull map shows the galaxy and none of the stars; the tidal local field is as large as the
   smooth one (100 to 130 percent). The Pull and Well depth modes therefore draw the local field by default with the total as a switch; Tidal strength draws the total.
6. **A generated sector holds less than the model expects in the sample measured** (see the report: 85 percent of the systems and 78 percent of the mass in a
   125-sector solar-circle block, 98 and 95 percent in one bulge sector), so the local field shows a small real deficit. The rule leaves it visible on
   purpose; the generator side is filed as its own bug.

## Colours

A perceptual, colour-blind-safe ramp on a log scale, with the same rule as
the Galaxy Map tint: density of the effect sets transparency and no fill
goes above 50% opacity, so stars and routes stay readable. The legend states
the unit under the Customary or metric choice (UX.36 number formatting).

## The size of a gradient measure comes from the orbital thresholds (planned)

Boss (2026-10-10 22:14Z): "add note to TODO for the future plan of gravitational map, use the preset minimum movement values from the orbital system to build the size of each gradient measure."

The orbital update system already fixes, per scale, the smallest movement that counts (`THRESHOLDS_M` in `physics/position.py`, section 3 of [orbital-updates.md](orbital-updates.md), pinned by `test_update_due.py`):

| Scale | Preset minimum movement | In metres | Objects |
|---|---|---|---|
| Galactic | 0.01 mpc (about 2 AU, 3.1e8 km) | `0.01 * PARSEC_M / 1000` | stars, black holes, neutron stars, stand-alone bodies |
| System | 0.01 AU (about 1.5e6 km) | `0.01 * AU` | planets, companion stars |
| Planetary | 100,000 km | `1.0e8` | moons and other satellites |

The plan, to be settled in GEN.198, GEN.199 and PERF.72:

- **Step of a gradient measure.** The tidal tensor and the size of the gradient are taken by differencing the field over a step. That step is the preset minimum movement of the scale being drawn: galactic for a sector grid, system for a System Map grid, planetary for a planet or moon neighbourhood. Below it the orbital system itself says an object has not moved, so a finer difference measures nothing the rest of the galaxy can see.
- **Size of a zone.** A zone edge is a whole number of those steps, so the gradient of a zone is the difference of two samples that are a real number of steps apart. The default sector grid (16 zones of 0.25 pc over a 4 pc edge) is 250 mpc, which is 25,000 galactic steps. A zone is never finer than one step; a setting that asks for finer zones is refused.
- **When a zone changes.** A source that has not moved by its own threshold has not moved in the orbital system, so it does not change the grid: the cache (GEN.199) is thrown away only for the zones a source left or entered once it passes its threshold, which is the same event that sets its next update due time (GEN.106).
- **One set of numbers.** The grid reads the thresholds from the same constants as the orbital update, not a copy, so changing a threshold there changes the gradient step and the cache rule together.

## Items

GEN.198 field evaluator; TEST.131 its tests; PERF.72 cost study; GEN.199
per-zone grid and cache (phase 2). MAP.167 umbrella; API.24 endpoint;
MAP.168 Sector Map layer; MAP.169 System Map layer; UX.94 wording and
legend (phase 3, phase 3 of 9.0). MAP.170 Galaxy Map layer (phase 3+).
Course planning (NAV.6) does not use the evaluator: it plans around keep-out spheres.
