# Interstellar object IDs (GEN.64)

Boss, 2026-10-02: every object in sector space that isn't a generated
star system gets a unique ID built from where it sits, and that ID is its
name. This covers rogue planets, standalone black holes and neutron
stars, nebulae, supernova remnants, quasars, interstellar comets,
asteroid fields, and star systems built around a bright-sweep star. A
remnant's collapsed core is named `<remnant ID> Core`, as before. Ordinary star systems keep their generated names.

## Layout

64 bits, printed as 16 uppercase hex digits (`objectId.format_id`):

| Bits | Field | Size | Meaning |
|------|-------|------|---------|
| 63-60 | type | 4 | `objectId.KIND_CODES`: 1 rogue planet, 2 black hole, 3 neutron star, 4 nebula, 5 supernova remnant, 6 quasar, 7 interstellar comet, 8 asteroid field, 9 bright-sweep system. 0 and 10-15 unused. |
| 59-57 | unit | 3 | 0 mpc, 1 cpc, 2 pc, 3 kpc, 4 Mpc, 5 Gpc. 6 and 7 unused. |
| 56-40 | distance | 17 | 0 to 131,071 in that unit: the smallest unit the distance fits in. |
| 39-20 | bearing | 20 | 0-360 degrees in 2^20 steps. |
| 19-0 | mark | 20 | 0-360 degrees in 2^20 steps. |

Bearing and mark are the galactic-frame course from the core to the
object (`navigation.course_between`), the same "bearing mark mark" the
NAV page uses. Boss picked this 64-bit layout over a 35-bit one (1-degree
angles, 0-999 distance), which would have given millions of objects the
same ID.

## Resolution and clashes

One ID covers about 0.003 x 0.003 x 0.01 pc at 500 pc from the core,
0.05 x 0.05 x 1 pc at 8 kpc and 0.09 x 0.09 x 1 pc at 15 kpc. Two objects
can still sit closer than that. `_db._claim_object_ids` gives the first
one, in generation order, the plain ID and bumps each later one by one
mark step (`objectId.bump`) until it is free in its own table, so the
result doesn't depend on which worker saves first. Measured on 20 dense
sectors at ring 700 (about 2.8 kpc out, 17,368 rogue planets), 30 IDs
(0.17%) were bumped. Distance is the coarse axis (1 pc steps past 1.3
kpc), so more bits, if ever needed, should go to distance first.

## Where it applies

- `insert_sector` claims IDs for every phenomenon and bright-sweep
  system of a galaxy-placed sector in one pass, before the name registry
  sees the rest, so these objects never touch `system_name_registry`.
- A standalone save with a galaxy position (`generate.py phenomenon
  --sector-id`) claims its ID the same way.
- An object with no galaxy position (`generate.py sector`, a
  `phenomenon` with no sector) keeps its generated name or designation.
- A bright-sweep system's stars, planets and moons are named from the ID
  like any system: `<ID> A`, `<ID> II`, `<ID> IIa`.
- Rows saved before this keep their names; GEN.39 already calls for a
  fresh galaxy.
