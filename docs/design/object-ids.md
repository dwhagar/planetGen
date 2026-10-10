# Interstellar object IDs (GEN.64)

**Built (Boss, 2026-10-09 22:39Z):** the `uid` column is an 80-bit birth-location ID (GEN.170 layout, DB.20 schema v78 with the sector fill giving the IDs; GEN.172, GEN.176 and API.23 follow; design in [object-id-options.md](object-id-options.md)). GEN.69's hash scheme is gone; GEN.64's position IDs are only the names of interstellar objects, and a bright-sweep system no longer keeps its position ID as its `uid`.

Boss, 2026-10-02: every object in sector space that isn't a generated
star system gets a unique ID built from where it sits, and that ID, in
hex, is its name. This covers rogue planets, standalone black holes and
neutron stars, nebulae, supernova remnants and their collapsed cores,
quasars, interstellar comets and asteroid fields. Star systems keep
their generated names, including one built around a bright-sweep star,
which was once known by its position ID (GEN.72) but is named by the registry, and the stars, planets and moons inside any system are named from
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

- `insert_sector` claims IDs for every phenomenon and remnant core of a
  galaxy-placed sector in one pass, before the name registry sees the
  rest, so these objects never touch `system_name_registry`. A
  bright-sweep system goes through the registry like any system.
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
4. **Names from the codec**: `planetgen.names.gated_phoneme_codec`
   (`src/planetgen/names/gated_phoneme_codec.py`, moved there from the repo root by GEN.120) turns an
   ID into pronounceable words and back, keyed by a domain and the naming
   key, so a name is unique because its ID is, and changing the key
   renames the codec's objects without rewriting rows. Boss's decisions
   of 2026-10-08 (03:57Z and 04:00Z): stars and sectors keep the
   word-salad method, and planets, moons and belts keep the "<system> I"
   pattern, so the codec names only the objects with no star-derived name
   (the GEN.64 kinds above, whose hex ID becomes words, and
   constellations). A wide binary's planets are "<word 1> I",
   "<word 1> II", star B's "<word 2> I", "<word 2> II", never "A I"
   (Boss, 2026-10-07 17:11Z; GEN.71), where the two words are the wide
   pair's word-salad name.
5. **No removal** of the word-salad code, the word lists, the nltk
   corpus or the name registries: stars, sectors, planets, moons and
   belts still use them.

### The naming key (GEN.70, control schema v9)

The key is 8 hex digits in the control database's `galaxy_naming` (one
row per galaxy database), drawn as the first 32 bits of SHA-256(galaxy
seed || `"naming-key"`) when `planetgen plan` makes a new seed (planning
again over the same seed keeps the stored key) and changeable on the Stats
page or by `POST /api/admin/naming-key`. `names/naming_key.py` turns an
ID, a kind and the key into a name: the codec domain is `<key>:<kind>`
for the kinds in `naming_key.KINDS` (the GEN.64 kinds without the
bright-sweep system, plus `constellation`). The codec derives its
permutation from the domain's CRC32, so a kind has 128 permutations and
two keys can occasionally name an object alike; the ID stays the
identity. `CODEC_VERSION` is stored beside the key and the Stats page
warns when the code's differs; bump it whenever the codec's golden words
change. Stars, sectors, planets, moons and belts do not use the key.

### Names shown for IDs (GEN.71)

The `name` column of a codec-named object keeps its 19-digit hex ID (the
identity, and what `_claim_object_ids` checks for collisions). Every JSON
answer of the API passes through `api/naming.NamingJSONProvider`, which
replaces a `name` that is such an ID with its codec name under the
galaxy's key (the kind comes from the ID's own type bits, so a remnant
core and a black hole name themselves correctly; a bright-sweep system's
ID is not renamed). The pages call the API in-process, so the maps,
lists, info panels and the nebula page all show the codec name. The key
is read once per request (cached 5 s per process, dropped when an admin
changes it) and `query.galaxy_content_state` mixes it into its `base`, so
a key change refreshes the tile and page caches. With no key drawn, or in
the CLI, the ID is shown. Lists sorted by name order by the ID, which
groups the kinds. `naming_key.stored_name_for` turns a typed codec name
back into its ID for a search that wants it. The wide-pair rule (never
"A I"; test in `test_body_names.py`) was already in force from GEN.62.

### Unique IDs for every object (GEN.69 and GEN.170, schema v78)

Every object has a `uid` column. A sector's is its designation as an
integer (`BIGINT UNSIGNED`, unique; it needs no row: `uid.sector_uid(ring,
layer, slot)` is the ID of a sector nobody has generated). Every other
object's is the 80-bit `BINARY(10)` of `planetgen/galaxy/object_uid.py`,
unique on its own: birth sector 40 bits, serial 28 bits, body number 12
bits (layout and rules in [object-id-options.md](object-id-options.md)
section 0).

- The sector fill (`store._UidIssuer`) gives a sector's systems, then its
  phenomena, generated serials 0, 1, 2 ... in insertion order, and a
  system's stars, planets, moons, belts and comets body numbers 1, 2, 3
  ... in insertion order. Saving the same sector again from the same
  galaxy seed gives every object the ID it had; position is never an input.
- A row saved later (an admin-added body or system, a system saved on its
  own, a sector with no grid address) is a run-time birth: `store.assign_uids`
  takes the next serial of the sector's `id_counters` row (kind 01), or for
  a body the system's next body number. A system with no sector is born at
  ring 0, layer 0, slot 0.
- A row that already has a `uid` keeps it, so a system whose content is
  regenerated in place (`replace_system_content`) keeps its own ID.
- Counters only move up and `planetgen reset` keeps them, so a number is
  never given twice.

### What the codec guarantees (GEN.120)

Tests: `src/tests/test_gated_phoneme_codec.py`. Its golden words are
pinned: moving one renames objects, so it needs a naming key or version
key change first.

- **Encoding** is deterministic and uses no database or hash round: the
  ID's digits are shifted by a domain-keyed permutation, split into words
  and spelled with gated phonemes.
- **Two IDs of the same length never share a name** (checked exhaustively
  for words of up to five digits, and over 100,000 random 19-digit IDs).
- **Decoding needs the ID's length.** A word can read as more than one
  length ("bar" is the digits `08` or `00B`), so about half of the names
  of 19-digit IDs have two readings of different lengths. `decode(phrase,
  domain, length=19)` is exact; without `length` it answers only when one
  length fits and otherwise says to pass it. The original file decoded by
  taking the longest phoneme at each step, which gave the wrong digits for
  about one ID in five from five digits up; decoding now accepts only
  the reading the encoder would have written. The encoder is untouched, so
  names are the same as with the original file.
- The module needs no packages beyond the standard library (the original
  used `pygtrie` and `more_itertools`, which this project does not install).

A bright star placed by the backfill shows its ID until its sector is
generated, and only then gets its name.
