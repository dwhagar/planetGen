# Interstellar object IDs (GEN.64, replaced by GEN.1256)

**Built (Boss, 2026-10-09 22:39Z):** the `uid` column is an 80-bit birth-location ID (the first layout; the default is 96 bits, `BINARY(12)`, Boss, 2026-10-15, see intergalactic-and-halo-contents.md and numeric-precision-audit.md) (GEN.170 layout, DB.20 schema v78 with the sector fill giving the IDs; GEN.172, GEN.176 and API.23 are built too; design in [object-id-options.md](object-id-options.md)).

**GEN.1256 (Boss, 2026-11-04, answers 1, 1, 1, 1):** the 76-bit position ID of GEN.64 is gone. Every object that is not a generated star system (a rogue planet, standalone black hole or neutron star, nebula, supernova remnant and its collapsed core, quasar, interstellar comet or asteroid field) keeps its Object ID (`uid`, `galaxy/object_uid.py`) in its printed form (`FE81000A2B-0000005-000`) as the stored `name`, and the API shows it by its local name (GEN.1255) with the birth sector beside it. The position-ID layout, its type codes, the collision chain and `_claim_object_ids` were removed with no alias. Star systems keep their generated names, and the stars, planets and moons inside any system are named from the system as before (`<ID> A`, `<ID> II`).

## Where it applies

- `insert_sector` marks every phenomenon and remnant core of a galaxy-placed sector as named by its Object ID (`store._name_by_object_id`, `PENDING_ID_NAME`); the INSERT that gives the row its `uid` writes the printed ID as its `name` (`store._named_by_uid`), so these objects never touch `system_name_registry`. A bright-sweep system goes through the registry like any system.
- A standalone save with a galaxy position (`planetgen phenomenon --sector-id`) takes its name in `assign_uids`.
- A name given by hand (`--name`) is kept.
- An object with no galaxy position (`planetgen sector`, a `phenomenon` with no sector) keeps its generated name or designation, and its core stays `<remnant> Core`.
- Migration `object_id_names` rewrote the 19-digit names of an older galaxy from the row's `uid`; a row with no `uid` kept its name.

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

### Local names from the Object ID (GEN.1255, `CODEC_VERSION` 2)

Boss chose (5:14 PM EST on 2026-11-04, research in `reports/id-names-without-region-prefix-2026-11-04.md`) that an ID-derived name is
the local part of the Object ID only, with the sector shown beside it, and that the 96-bit ID stays the identity. The position-ID name above
spells a whole 19-digit ID as four words, so objects born near each other share their leading words; the Object ID's place is 56 of its
96 bits and cannot be spelled word by word on a boundary. `naming_key.local_name(uid, key)` spells only the serial and body number
(28 + 12 bits, 10 hex digits) as two words, in its own codec domain `<key>:local` (not a member of `KINDS`; no object of that kind exists),
so the sector's words never appear: `Zaidonr Jaisonn`.

Spelled directly, serials 0, 1, 2 share their first word. The 26 count bits of the serial (the two kind bits stay as they are) are first
scrambled by `y = (a * x + b) mod 2**26`, a bijection for odd `a`, with `a` and `b` drawn from SHA-256 of the key and the sector, so every
sector has its own scramble. `a` is the first candidate that keeps the first 4,096 serials of the sector at least 256 apart after the
scramble (`SPREAD_SERIALS`, `SPREAD_GAP`): a first word holds all but the lowest eight bits of the scrambled count, so those serials
have first words of their own. Inside a sector a local name is exact (the spelling is a bijection of the 40 bits) and reverses with
`naming_key.uid_of_local_name(name, sector, key)`; on its own it can be shared with an object of another sector (about 1 in 27 at
10 million objects), so a name that stands alone is shown with its sector, `local_name_with_sector`: `Zaidonr Jaisonn, sector
0000B0001803E8`. A key change renames every local name without rewriting a row. The bodies of one system share the first word and
differ in the second.

What this item does not do: it adds the functions, the version and the tests. Using the local name for rogue planets and the other
phenomena, and showing the name and sector on pages, the API and search are GEN.1256, GEN.1257, UX.259
and API.38.

### Names shown for IDs (GEN.71, GEN.1256)

The `name` column of a codec-named object keeps its printed Object ID (the identity). Every JSON answer of the API passes through `api/naming.NamingJSONProvider`, which replaces a `name` that is such an ID with its local name under the galaxy's key (`naming_key.local_name`) and adds `name_sector`, the birth sector's 14 hex digits, beside it. The pages call the API in-process, so the maps, lists, info panels and the nebula page all show the local name. The key is read once per request (cached 5 s per process, dropped when an admin changes it) and `query.galaxy_content_state` mixes it into its `base`, so a key change refreshes the tile and page caches, and nothing is stored in words, so a key change renames all at once. With no key drawn, or in the CLI, the ID is shown. Lists sorted by name order by the printed ID, that is by birth sector and then serial. `GET /api/lookup?text=` accepts a printed ID, or a local name followed by its sector (`Banhunm Finkap, sector 00FA0FFFE9E240`, `naming_key.uid_of_typed_name`); a local name without its sector is ambiguous and finds nothing by itself. The wide-pair rule (never "A I"; test in `test_body_names.py`) was already in force from GEN.62.

### Unique IDs for every object (GEN.69 and GEN.170, schema v78)

Every object has a `uid` column. A sector's is its designation as an
integer (`BIGINT UNSIGNED`, unique; it needs no row: `uid.sector_uid(shell,
band, wedge)` is the ID of a sector nobody has generated). Every other
object's is the 80-bit `BINARY(10)` (the first layout; 96 bits, `BINARY(12)`, is the default, Boss, 2026-10-15) of `planetgen/galaxy/object_uid.py`,
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
  the core (shell 0, band 0, wedge 0).
- A row that already has a `uid` keeps it, so a system whose content is
  regenerated in place (`replace_system_content`) keeps its own ID.
- Counters only move up and `planetgen reset` keeps them, so a number is
  never given twice.
- A field cloud (GEN.47, GEN.176) is stored by whichever sector reaches it
  first, but is born in the sector holding the centre of the space it
  fills: the mean of the interior points of its shape on the fixed
  24-cell grid (`NebulaShape.interior_centroid`), scaled by its radius,
  added to its centre and rounded to 1 mpc (`nebula_field.centroid_mpc`).
  Its serial is field-drawn (kind 10): its rank among every field cloud
  born in that sector, cells in index order and clouds in a cell in draw
  order (`nebula_field.cloud_object_id`), worked out from the galaxy seed
  alone, so the ID does not depend on save order and the birth sector
  need not exist yet. Nebulae and remnants a sector generates itself
  keep a generated serial in that sector.
- A container (DB.85: a globular cluster or a superbubble, a row of
  `galaxy_containers`) is born in the sector holding its centre and has a
  region serial (kind 11): its rank among the containers born in that
  sector, counted in a fixed order (clusters by name and seed, then
  bubbles by column and draw; `generation/containers.build`), worked out
  from the galaxy seed alone.

### Members of a belt (GEN.1065)

A belt's members take no body number of their own (the body field is 12 bits and a belt can hold far more bodies than that). A member is designated `<belt ID>.<rank hex>`, such as `FE81000A2B-0000005-003.0042`: the belt's printed ID, a point, and the member's rank by the mass drawn once at birth (1 is the heaviest), at least four hex digits and up to eight. `galaxy/object_uid.designation(belt_id, rank)` prints it and `parse_designation(text)` reads it back (the belt's ID by `parse_id`); both are pure, read no stream and no naming key, and give a member one text. A member has no name: `naming_key` lists the kinds it renames and belt members are not among them. A rank above the belt's individual count (`is_on_demand`) names a body drawn on demand and never stored. A belt member that a planet captures keeps this designation and gains a parent link to the planet (`captured_member(belt_id, rank, planet_id)`, the planet in the same system); a body from a reservoir that is no belt of the system prints `<planet ID>.<rank hex>` (`outside_reservoir_member`), the body field untouched; `parse_captured` reads either back (TEST.495, [captured-asteroid-moons.md](../systems/captured-asteroid-moons.md) section 4). Design: [belt-evolution.md](../systems/belt-evolution.md) section 9.

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

## Scatter rows (GEN.1254)

A star of `bright_stars` and a scattered phenomenon of `phenomenon_scatter` have an Object ID from the moment the scatter writes them: `uid`
(`BINARY(12)`, unique, NULL for a row made before migration 0213). It is the ID the object keeps when its sector is charted. The module is
`db/scatter_uids.py`.

- **Birth sector and body.** The row's own sector, body number 0, serial kind generated.
- **Serial.** The row's rank in its sector among the rows written in the same call, brightest first (`luminosity_w` of a star,
  `map_luminosity_sol` of a phenomenon), then the row's `seed`, then its row number. The rank continues from the sector's counter (`id_counters`
  kinds `scatstar` and `scatphen`, scope the sector's seven bytes), so a backfill that adds dimmer stars later takes the next numbers and no
  number is given twice. The same rows give the same IDs whatever order they were inserted in, and a cleared scatter (`clear_bright_stars`,
  `clear_phenomenon_scatter`) forgets its counters and starts from serial 0 again.
  **A scatter serial depends on the run order, by design (GEN.1302):** the rank is within one write, so a sector whose dimmer stars an earlier
  backfill wrote and whose brighter ones the first fill added numbers them differently from one whose first fill wrote them all at once. The IDs
  stay unique and never change once issued; `test_order_independence` puts a sector's scatter serials in brightness order before it compares.
  Run-order-independent scatter IDs (a serial block for each mass ring) would be a new item with an epoch.
- **Stars and phenomena do not share a counter.** Stars take serials from 0, phenomena from 2^25 (`PHENOMENON_SERIAL_BASE`), 33 million each.
- **When.** `store.insert_bright_stars` and `store.insert_phenomenon_scatter` call `scatter_uids.assign` when they have written their rows:
  one window query ranks the rows with no `uid` per sector, three joined updates move the counters and write the IDs. A row inserted by other means
  stays NULL (NULLs do not collide in the unique index) until a writer runs.
- **The fill.** A sector's fill starts its own serials past the stars' counter (`_UidIssuer`, `scatter_uids.next_star_serial`), so a system the
  fill makes is never given a serial a scatter star holds. A system built round a scatter star takes the star's `uid` (`entry.bright_star_uid`), a built
  scattered phenomenon takes its row's (`entry.scatter_uid`), and `POST /api/uncharted-systems/<id>/generate` does the same. The old ID of a
  deleted system comes back when its star is built again.
- **Not done here.** The row number of an uncharted star is still in its URL and its page title; naming from the ID is GEN.1255's codec, applied to
  these rows by GEN.1256 and GEN.1257.
- **Cost.** The column and its unique index, about 24 bytes a row: 0.65 GB at 26.9 million rows, 5.3 GB at 220 million. The rank is an extra pass
  over the rows a call wrote, not over the table.
