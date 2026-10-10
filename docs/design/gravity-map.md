# Gravity map: the field of a sector as a heat map

Boss (2026-10-10 20:57Z): the orbits are worked out from point-mass
vectors, so show "a gravitational map of a sector, probably using colors
and calculating it per zone ... the gravitational gradients of the vector
field as kind of a heat map of the inside of the sector". Foundations are
phase 2; the feature roll-out is phase 3 and the seminal feature of 9.0.

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
- **Sources.** Every point mass in the sector and its neighbours: stars,
  planets, moons, compact objects and stand-alone facilities with a mass
  (the same influence set as GEN.109's Hill-sphere rule), plus the galaxy's
  own smooth potential (GEN.115) as the background. Plummer softening keeps
  a zone that holds a source finite.
- **Near and far.** Sources close to a zone are summed exactly; distant
  ones are summed as one aggregate per region (mass and centre of mass,
  re-using the MAP.151 region pyramid), the way a tree code does. The
  cost study (PERF.72) fixes the split and checks it against the exact sum.
- **Cache.** The grid is derived data, cached like the Galaxy Map tiles and
  thrown away when the orbital update (GEN.105) moves the sources; it is
  never part of the galaxy itself, so it needs no galaxy-table migration.
- **Slices and volumes.** The map shows one slab (a plane slice) at a time
  and can show the whole volume as transparent layers; on the System Map it
  is the orbital plane.

## Colours

A perceptual, colour-blind-safe ramp on a log scale, with the same rule as
the Galaxy Map tint: density of the effect sets transparency and no fill
goes above 50% opacity, so stars and routes stay readable. The legend states
the unit under the Customary or metric choice (UX.36 number formatting).

## Items

GEN.198 field evaluator; TEST.131 its tests; PERF.72 cost study; GEN.199
per-zone grid and cache (phase 2). MAP.167 umbrella; API.24 endpoint;
MAP.168 Sector Map layer; MAP.169 System Map layer; UX.94 wording and
legend (phase 3, the 9.0 headline). MAP.170 Galaxy Map layer (phase 3+).
Course planning (NAV.6) reuses the evaluator instead of its own pull model.
