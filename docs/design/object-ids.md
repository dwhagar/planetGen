# Interstellar object IDs (GEN.64)

Boss, 2026-10-02: every object in sector space that isn't a generated
star system gets a unique ID built from where it sits, and that ID, in
hex, is its name. This covers rogue planets, standalone black holes and
neutron stars, nebulae, supernova remnants and their collapsed cores,
quasars, interstellar comets, asteroid fields, and star systems built
around a bright-sweep star. Ordinary star systems keep their generated
names, and the stars, planets and moons inside any system are named from
the system as before (`<ID> A`, `<ID> II`).

## Layout

76 bits, printed as 19 uppercase hex digits (`objectId.format_id`):

| Bits | Field | Size | Meaning |
|------|-------|------|---------|
| 75-70 | type | 6 | `objectId.KIND_CODES`: 1 rogue planet, 2 black hole, 3 neutron star, 4 nebula, 5 supernova remnant, 6 quasar, 7 interstellar comet, 8 asteroid field, 9 bright-sweep system, 10 remnant core (black hole), 11 remnant core (neutron star). 0 and 12-63 unused. |
| 69-67 | unit | 3 | 0 mpc, 1 cpc, 2 pc, 3 kpc, 4 Mpc, 5 Gpc. 6 and 7 unused. |
| 66-48 | distance | 19 | 0 to 524,287 in that unit: the smallest unit the distance fits in. |
| 47-26 | bearing | 22 | 0-360 degrees in 2^22 steps. |
| 25-4 | mark | 22 | 0-360 degrees in 2^22 steps. |
| 3-0 | collision | 4 | 0-15: tells apart up to sixteen objects of one type at one position. |

Bearing and mark are the galactic-frame course from the core to the
object (`navigation.course_between`), the same "bearing mark mark" the
NAV page uses. Boss first picked a 64-bit layout over a 35-bit one
(1-degree angles, 0-999 distance), which would have given millions of
objects the same ID, then (06:44Z) added the 2-bit collision number in
place of bumping the mark and 2 more type bits so remnant cores get
their own type, and (06:48Z) 2 more bits each for distance, bearing
and mark and 2 more for the collision number.

## Resolution and collisions

One position covers about 0.0008 x 0.0008 x 0.001 pc at 500 pc from
the core, 0.012 x 0.012 x 1 pc at 8 kpc and 0.022 x 0.022 x 1 pc at
15 kpc. Two
objects of one type can still sit closer than that.
`_db._claim_object_ids` gives the first one, in generation order,
collision number 0 and each later one the next number
(`objectId.bump`), checking its own table, so the result doesn't depend
on which worker saves first. A seventeenth object at one position moves on one
mark step with collision number 0. Measured on 20 dense sectors at ring
700 (about 2.8 kpc out), the 68-bit draft gave 29 of 16,526 rogue
planets (0.18%) collision number 1; the 76-bit layout gives none. Distance is the coarse axis
(1 pc steps past 5.2 kpc), so more bits, if ever needed, should go to
distance first.

## Where it applies

- `insert_sector` claims IDs for every phenomenon, remnant core and
  bright-sweep system of a galaxy-placed sector in one pass, before the
  name registry sees the rest, so these objects never touch
  `system_name_registry`.
- A standalone save with a galaxy position (`planetgen phenomenon
  --sector-id`) claims its ID the same way.
- A name given by hand (`--name`) is kept.
- An object with no galaxy position (`planetgen sector`, a
  `phenomenon` with no sector) keeps its generated name or designation,
  and its core stays `<remnant> Core`.
- Rows saved before this keep their names; GEN.39 already calls for a
  fresh galaxy.

## Planned: an ID for everything, and names from IDs

Boss (2026-10-03, 2026-10-07) asked to drop word-salad names entirely
and to give everything in the galaxy a unique ID, including sectors that
haven't been filled. The plan:

1. **Research** the cheapest way to make those IDs: this document's packed
   position ID, an ID derived from a sector's address (so an unfilled
   sector has one before any row exists), a hash of the seed and address,
   or a database sequence; compared for cost, collision risk and stability
   when content moves or is regenerated in place.
2. **IDs for every object**: sectors, systems, stars, planets, moons,
   belts, comets and phenomena.
3. **A naming key** in the control database, drawn from the galaxy seed
   when the galaxy is created and changeable by an admin.
4. **Names from the codec**: `gatedPhonemeCodec.py` (repo root) turns an
   ID into pronounceable words and back, keyed by a domain and the naming
   key, so a name is unique because its ID is, and changing the key
   renames everything without rewriting rows. Planets keep the "<system>
   I" pattern. A wide binary gets a two-word name: star A's planets are
   "<word 1> I", "<word 1> II", star B's "<word 2> I", "<word 2> II",
   never "A I" (Boss, 2026-10-07 17:11Z; GEN.71).
5. **Removal** of the word lists, the nltk corpus, the name registries
   and the collision rules.

A bright star placed by the backfill shows its ID until its sector is
generated, and only then gets its name.
