# Unique IDs for every object: options, measurements, recommendation

Research only. No code, TODO.md or PR was changed. Boss decides before anything is filed.
Measured on the build box (4 cores, Python 3.13, MariaDB 10.11.14, InnoDB buffer pool 256 MB); scripts and raw results are in `bench/` beside this file. Code read: `main` at 4113872 (after PERF.49).

## 0. Decision and revised design (Boss, 2026-10-09 22:39Z)

Boss: fix the defects, use birth location + serial **galaxy-wide**, up to 128 bits but shorter preferred, **the same length for every object**, and an ID that always identifies that one object. Nebulae are located by the geometric centre of the space they occupy. This supersedes section 6 where they differ (section 6 had short, system-scoped body IDs).

### 0.1 One layout for every object: 80 bits, 20 hex digits, `BINARY(10)`

```
| birth sector address | serial in that sector | body number |
|       40 bits        |       28 bits         |   12 bits   |
|     10 hex digits    |      7 hex digits     |  3 hex digits|
```

- **Birth sector address (40 bits):** ring 12, biased layer 12, slot 16. Fits the default Milky Way (3,763 rings, 2,041 layers, 23,040 slots in the outer ring) with 8% spare on rings. It re-packs today's 45-bit designation so the whole ID is nibble-aligned; a galaxy whose bounds do not fit picks wider fields at plan time (stored in `galaxy_shape`) and the ID grows to 96 or at most 128 bits. Length is fixed inside a galaxy.
- **Serial (28 bits):** a top-level object's number in its birth sector. Top two bits split the space so no allocator ever needs another's state: `00` generated (the object's rank in generation order, up to 67M a sector; omega-Centauri-class cores are a few million), `01` added at run time (admin, ejection, facility), `10` field-drawn (nebulae and other objects stored by a sector other than the one holding their centre, below).
- **Body number (12 bits):** `000` is the top-level object itself (system, rogue planet, nebula, facility, ...). `001` and up number the stars, planets, moons, belts and comets born in that system, one counter across all five kinds, so a captured moon, a merged planet or a comet that becomes a planet keeps its number. 4,095 births a system for its lifetime.
- Printed `FE81000A2B-0000005-000` (system) and `FE81000A2B-0000005-01A` (a body). Every object has the same 20 digits. Cost to store: 10 bytes, against 12 (system) or 8 (body) today.
- The kind is **not** in the ID, because kind changes (a planet ejected from its star becomes a rogue planet and keeps its ID). Lookup by ID probes the object tables by the unique `uid` index; a URL or reference may carry the kind as a hint (`planet:FE81...`), never as identity.
- A 64-bit alternative exists (36-bit sector ordinal + one flat 28-bit counter per sector, `BIGINT`), but it loses the system/body split, caps a sector at 67M objects and has no room for a bigger galaxy; I do not recommend it.

Capacity check: 2.95e11 systems and about 1.3e13 objects over about 3.3e10 sectors is 400 objects a sector on average. Densest known sector (a 47 Tuc-class core, 1e5 to 5e5 systems) fits 67M.

### 0.2 Rules that make it unique and stable

1. **Assigned once, at birth, never recomputed** from rank, position, sector, parent, table or content. The fill gives generated serials by generation order inside the sector (one worker owns a sector, so no coordination and no database read). This is also the fix for defect 1.
2. **Run-time births** (admin-added body, ejected object, facility, merger remnant, new system placed by hand) take the next value of a counter in a small `id_counters` table keyed by sector (bodies: by system) using the upsert-and-`LAST_INSERT_ID` pattern `_reserve_id_block` already uses. Counters only grow and survive deleting the sector or object, so a run-time serial is never given twice. This is the fix for defect 2.
3. **A deleted or merged object's ID is never given to another object.** The only reuse is regenerating a whole sector from the seed, which makes the same objects again at the same ranks. Regenerating a body or phenomenon in place keeps its ID: add `uid` to `_KEPT_PHENOMENON_COLUMNS` and keep the row's serial on a system's content swap. This is the fix for defect 3.
4. **Events:** move and reparent keep the ID; an ejected planet keeps its ID (now a rogue-planet row, body number unchanged); a merge keeps the heavier body's ID and retires the other; a split gives each fragment a new run-time serial.
5. **The row `id` stays as the foreign key and join key.** The ID becomes the public reference (pages, URLs, API, wiki links, `objectref`), since Boss wants it to identify the object. Row ids stay out of anything user-facing. No compatibility shim.
6. **Position IDs (GEN.64) stay as the names** of interstellar objects. The ID column no longer depends on them, so `_claim_object_ids` is needed only for name uniqueness.

### 0.3 Nebulae and other multi-sector objects

A nebula occupies many sectors, so its birth sector is the sector holding the **geometric centre of the space it occupies**: the centroid of the interior of its metaball field (`nebula_shape.py`), not the metaball origin, which the warp and the 0.55 ellipsoid offset move away from the centroid. Remnants use the centre of their shell.

- Centroid method: count interior points of the shape on a fixed grid (the existing 24-cell fit grid) in integer arithmetic, average them, scale by the radius, add the stored centre, round to 1 mpc as an integer, and take the sector from the rounded integers. Rounding first removes platform float noise except when a centroid lies within about 1e-12 of a half-mpc step of a sector face; I did not measure how often that happens.
- Today a field cloud is stored by "the first sector saved that the cloud reaches" (`_insert_field_nebulae`), so its home sector depends on save order and on which worker was quicker. Under the new rule the birth sector is the centroid's sector whoever stores it. Because that worker does not own the centre sector, the cloud takes a **field-drawn serial** (`10` prefix): the cloud's rank among the clouds of its field cell (`nebula_field.cell_clouds(seed, shape, (i,j,k))`) that have their centroid in that sector, which every worker computes identically. No shared counter, and the centre sector need not exist yet.
- The same rule covers any other object that is stored by a sector other than its own.

### 0.4 What this does to the earlier ranking

The recommendation is unchanged in kind (A, birth-address serial) but now global and fixed-width. A hash cannot meet the new requirement at 80 bits: a flat 80-bit hash over 1.3e13 objects collides with probability near 1. Section 3 numbers still hold; the 63k to 67k inserts/s figures were for a 12-byte key, and a 10-byte key will not be worse. Re-measure on the real fill order when this is built.

## 1. What is in the code today

Three IDs coexist, and only one of them is used to find anything.

| ID | Where | Form | Used for |
|---|---|---|---|
| Row `id` | every table | `BIGINT`, handed out in per-process blocks (`_allocate_id`, `id_blocks`) | Foreign keys, URLs (`system:12`, `objectref.py`), the API. Differs between runs and workers. |
| `uid` (GEN.69, PERF.44) | `uid` column on every object table | Sector: designation as `BIGINT`. System and phenomenon: `BINARY(12)`, first 96 bits of SHA-256(seed, `uid-<kind>`, parent, rank), top bit set. Star, planet, moon, belt, comet: `BIGINT` (64-bit hash), unique only under its system | Written on every save, printed in detail views. Nothing looks an object up by it. |
| Position ID (GEN.64) | `name` and `uid` of rogue planets, standalone BH/NS, nebulae, remnants, quasars, comets, asteroid fields, bright-sweep systems | 76 bits: type 6, unit 3, distance 19, bearing 22, mark 22, collision 4 | Naming (the codec turns it into words) |

`rank` in the uid hash is the row's place among its parent's rows in id order, counted at insert time by `_UidIssuer` or recounted by `assign_uids`.

### Cost of the current approach

- Compute: **1.3 to 1.4 µs per ID** (`derived_uid`, SHA-256). A sector of about 170 rows costs 0.24 ms, against a row INSERT of 70 to 90 µs. Compute is not the cost, and PERF.44 already removed the SELECT-and-UPDATE pass.
- Index: a random-valued unique index. In a 2-million-row model of `star_systems` (table: row id, sector id, 250-byte payload, sector index; filled one sector at a time in Hilbert order) adding the hashed `BINARY(12)` unique index took inserts from **70.4k to 45.7k rows/s (-35%)** and the table from **370 to 421 bytes a row (+14%)**.
- GEN.64 position IDs: **5.0 µs** to pack (trig), plus `_claim_object_ids`, which SELECTs the table for held IDs once per table per sector and again for every bump.

### Defects found, reproduced against MariaDB with the repo's own code

1. **A re-run of the uid pass can fail with a duplicate key.** Delete the first planet of a 9-planet system through `save_system_edits`, add a planet (it saves with a NULL uid), then run `assign_uids`: `IntegrityError 1062 Duplicate entry for key 'uq_planets_uid'`. The new planet is ranked by its place among the *current* rows, so it gets the hash of a rank a surviving sibling already holds. `add_system_to_sector` runs the sector-wide pass over every system in the sector (`store.py` 6040), so by code reading an admin adding a system to a sector can fail this way; I ran only the system-level call. `generation-determinism.md` section 6.2 already warns that rank can shift.
2. **Bodies added by an admin never get a uid.** `save_system_edits` inserts new bodies with no uid and nothing assigns one (the new planet's uid stayed NULL).
3. **Regenerating a phenomenon wipes its uid.** `replace_phenomenon_content` copies every column except `_KEPT_PHENOMENON_COLUMNS`, which leaves out `uid`. A rogue planet with uid `0000004986A0FFFE64000000` had `NULL` afterwards (name kept). The docs say regeneration keeps the ID.

Also: the id of a body depends on its parent (`planet uid` goes into a moon's hash). Orbital updates change parents (capture, ejection), change tables (a planet ejected from its star becomes a rogue planet, two planets merge into a star, a body disrupts into an asteroid field: `orbital-updates.md` sections 1 and 3, `collisions-and-mergers.md` section 5), and move systems between sectors. An ID that encodes the parent or the table, or is recomputed from where the object is, will not survive those.

## 2. What the ID has to do

R1 unique (galaxy-wide for top-level objects, under the system for bodies). R2 survives edits and in-place regeneration. R3 survives a move to another sector. R4 survives reparenting, reclassification, merge, split (a rule is needed for the last two whatever the scheme). R5 same seed gives same IDs with parallel workers. R6 computable before the INSERT, with no lock, sequence call or SELECT. R7 cheap in the index and clustered usefully for sector queries. R8 readable and ideally decodable to an address. R9 cheap to migrate (Boss: no compatibility needed).

## 3. Numbers

### 3.1 Collision odds (birthday bound, p = 1 - exp(-n²/2^(b+1)))

Object counts: 2.953e11 systems (default galaxy, `phenomenon-scatter-mass-cut.md`); about 45 bodies per system from the 449-system run in `generation-performance-study.md` (1.3 stars, 9 planets, 25 moons, 0.7 comets, rest rogue/other), so about 1.3e13 objects. The full galaxy will never be generated; these are the ceilings. The bound matched a 32-bit simulation (0.75 seen, 0.69 predicted).

| Namespace | n | 64 bit | 95 bit (today's 96 with fixed top bit) | 128 bit |
|---|---|---|---|---|
| All systems | 2.95e11 | certain | **1.1e-6** | 1.3e-16 |
| All objects in one flat space | 1.3e13 | certain (about 4.6e6 clashes) | 2.1e-3 | 2.5e-13 |
| Bodies of one system (about 60) | 60 | 9.8e-17 | | |
| One sector's systems (about 10) | 10 | 2.7e-18 | | |

Bits for p below 1e-9: 105 for systems, 116 for every object. So today's split is sound: 96-bit galaxy-wide for systems, 64-bit scoped to the system for bodies. A flat 64-bit hash for everything would not be unique. A hash collision is also a hard failure: the same input gives the same ID, so the sector can never save (the GEN.64 bump counter has no hash equivalent).

Positional (sector + offset, Boss's idea): at 1 mpc (12 bits per axis, 81 bits with the 45-bit designation) about 10 systems a sector gives 7e-10 clash chance per sector and about 21 expected clashing pairs across 3e10 sectors if systems were placed uniformly (placement keeps Hill-sphere spacing, so fewer; not measured). At 0.1 mpc, 0.02 expected. Never zero, so it still needs a collision counter and a check.

### 3.2 Compute per ID (Python, one thread; `bench/ids_bench.txt`)

| Scheme | µs/ID |
|---|---|
| Counter (sector key plus serial) | 0.09 |
| Snowflake-style | 0.21 |
| Truncated SHA-256 / BLAKE2b / SHA-1 | 0.6 to 0.65 |
| Current `derived_uid` | 1.3 to 1.4 |
| Sector designation plus 3x12-bit offset | 1.7 (mostly `sector_uid`'s validation) |
| UUIDv4 | 1.3 |
| UUIDv5 (SHA-1, 16 B) | 1.8 |
| UUIDv7-style | 1.0 |
| ULID string | 3.1 |
| GEN.64 position pack | 5.0 |
| Morton 3D (magic numbers) | 1.4 |
| Hilbert 3D (pure Python) | 17.9 |

### 3.3 Index and insert cost (2,000,000 rows, 10 per sector, 200,000 sectors, filled in Hilbert sector order; `bench/db_results.jsonl`)

One run each, so treat differences under about 10% as noise.

| Variant | Inserts/s | Bytes/row (data+indexes) | Notes |
|---|---|---|---|
| No uid (row id plus sector index) | 70.4k | 370 | floor |
| **Hash uid, BINARY(12), unique (today)** | **45.7k** | **421** | random index |
| UUIDv4 as clustered PK | 11.8k | 514 | table 3x the buffer pool: I/O bound; gets worse with size |
| UUIDv7 as clustered PK | 71.6k | 393 | in time order; not deterministic; point lookup 2x slower here (extra hop) |
| Sector designation + serial, `BINARY(12)` unique | 63.5k | 463 | index 268 MB: designation order is not the fill order, so pages half fill |
| Hilbert-ordinal sector + serial, `BINARY(12)` unique | 67.2k | 415 | prefix follows the fill order |
| Clustered PK (designation, serial), no row id | 71.6k | 719 | row pages half fill when filled out of key order |
| Clustered PK (Hilbert ordinal, serial), no row id | 87.5k | 305 | append-like; no secondary indexes, so not like for like |

Reading: what costs is key order against insert order, not key width. Random keys (hash, UUIDv4/5) scatter; keys that start with a sector and are inserted sector by sector are 39 to 47% faster than the hash in the same model. Designation order costs index space unless the fill follows it.

### 3.4 Does the sector key order matter for queries?

For a 3-D neighbourhood of sectors (400 random centres, 229k-sector toy grid, `bench/locality.txt`), pages touched at 4 sectors per page: radius 8 pc, designation 16.6, Hilbert 15.6, Morton 17.2; radius 20 pc, designation 162, Hilbert 155, Morton 167. The existing designation order (ring, layer, slot) is within 6% of Hilbert for reads and had fewer contiguous runs. **A curve key is not worth it for queries**; the only case for one is matching the fill order on insert.

## 4. The options

| # | Scheme | Unique | Stable on edit / move / reparent | Seed-deterministic, parallel-safe | Before INSERT | Index | Readable / decodable | Notes |
|---|---|---|---|---|---|---|---|---|
| A | **Birth-address serial**: sector designation + serial counted at birth (+ body serial) | By construction, no collisions | Yes if the serial is stored and never renumbered | Yes: serial is generation rank, one worker per sector | Yes, 0.1 µs | Sector-prefixed, clusters (63 to 67k/s) | Yes: `FE81000A2B.3F` decodes to sector and serial | Needs a stored counter for runtime additions |
| B | A plus keyed scramble (Feistel / the GEN.120 codec's permutation) | By construction (bijection) | As A | As A | Yes, a few µs | Scrambled, scatters like the hash unless the stored column is the unscrambled one | Opaque unless decoded with the key | Only if Boss wants non-guessable IDs; nothing here needs that |
| C | Hash of birth address (today, GEN.69) | Probabilistic (1e-6 systems) | Birth hash survives moves; **rank-derived hashes broke on delete+add** | Yes | Yes, 1.3 µs | Random: -35% inserts | No, one-way | Fixing defect 1 means storing a birth index, which is A; the hash then adds only risk and scatter |
| D | UUIDv5 (name-based, SHA-1) | Probabilistic, 128-bit | As C | Yes | Yes, 1.8 µs | 16 B random | No | C in a standard wrapper; MariaDB has a `UUID` type, MySQL 8.4 does not, so portable only as `BINARY(16)` |
| E | UUIDv4 / 128-bit random | Probabilistic | Yes once stored | **No** | Yes | Worst: -83% | No | Only for rows that are not generated |
| F | UUIDv7, ULID, KSUID, Snowflake, Sonyflake, TSID | Yes (clock plus worker plus sequence) | Yes once stored | **No**: time of creation | Yes | Append-like (71.6k/s) | Partly | Right for non-spatial runtime rows (accounts, jobs, API keys, audit rows); wrong for generated objects, since a re-run gives new IDs |
| G | Pure positional, GEN.64-style (galactic polar, or sector + offset) | Needs a collision counter and a DB check | **Only as a frozen birth position** (as Gaia does) | Within a sector, yes; across a sector face two workers can claim the same cell | Needs `_claim_object_ids` SELECTs | Clusters spatially | Decodes to a place that goes stale | Planets, moons and anything orbiting move every tick, so it cannot identify bodies |
| H | Curve key of position (Morton, Hilbert, HEALPix) as the ID | Not injective | Moves with the object | Yes | Hilbert 18 µs in Python | Clusters | No | A sort/cluster key, not an identity. Gaia's `source_id` pairs a curve cell at first detection with a running number, which is A with a different prefix |
| I | Per-sector counter / DB sequence alone (no sector prefix) | Yes | Yes | Not unless one allocator | Needs a round trip (`id_blocks` today) | Sequential | No | This is the row id; keep it as the join key |
| J | Content hash (Merkle leaf, fingerprint) | Collides for identical content | **No**: any edit changes it | Yes | Yes | Random | No | Right tool for change detection (`generation-determinism.md`), not identity |
| K | Natural keys (name) / composite keys with no surrogate | Names are unique via the registry | Renames break it | Yes | Registry round trip | Wide | Yes | Names change (GEN.62 cascades); composite keys widen every child index |
| L | Hierarchical path text (`S12.3.0/y4/p2`) | By construction | As A | Yes | Yes | Variable-length text | Best | The display form of A, not a storage form |

## 5. Ranking

1. **A, birth-address serial.** Only option that has zero collision risk, needs no hash, SELECT or lock, decodes without a lookup, clusters by sector, and keeps its ID through every move because the address in it is where the object was born, not where it is. In my model it is 39 to 47% faster to insert than today's hash.
2. **B, A with a keyed scramble.** Same properties plus opacity, at some cost in decodability and clustering. Only if non-guessable IDs are wanted.
3. **C, keep the hash** and fix the three defects (store a birth index, keep `uid` in the regenerate path, assign on admin add). Works, but the stored birth index needed for the fix makes the hash redundant.
4. **D, UUIDv5.** C with a standard wrapper.
5. **G, positional** for interstellar objects only, as now.
6. **F, time-ordered IDs** for non-generated, non-spatial rows.
7. **H, I, K, L, J, E**: not identities.

## 6. Recommendation: A, in detail

**Principle: an ID is assigned once at birth and is never recomputed from anything that can change** (position, sector, parent, table, content). Where an object is now is an attribute (`sector_id`, position); where it was born is its ID.

Layout, with no column-width change from today:

- **Sector:** its designation, unchanged (exists before any row, GEN.68).
- **System and every top-level object** (rogue planet, BH/NS, nebula, remnant, quasar, comet, asteroid field, facility on its own): `BINARY(12)` = 8-byte birth-sector designation + 4-byte serial from that sector's counter. The serial is the object's generation rank at fill; objects added later (admin, runtime, ejection) take the next value of `sectors.next_serial`, which only grows, so deleted serials are never reused. Four billion per sector is far above any sector's population.
- **Star, planet, moon, belt, comet:** the existing `BIGINT uid`, now a per-system serial from `star_systems.next_body_serial` (generation rank at fill, then the next value). One counter across all five kinds, so a moon that is captured, a planet that merges, or a comet that becomes a planet keeps its ID. Full form: `<system ID>.<serial>`.
- Display form: `FE81000A2B.3F` for a system, `FE81000A2B.3F.7` for a body. Hex and decodable, which suits the NAV page and URLs.
- Position IDs (GEN.64) stay as the **names** of interstellar objects, since Boss chose them; they stop doubling as the `uid`, so `_claim_object_ids` is no longer on the path for uniqueness of the uid (still for the name).
- Identity events (Boss to confirm): move and reparent keep the ID. A planet ejected becomes a rogue planet and **keeps its ID** (system ID plus serial, now living in the rogue table) with a `born_in_system` link; a merge keeps the heavier body's ID and retires the other's (never reused); a split gives the fragments fresh serials.
- Needs two counter columns (`sectors.next_serial`, `star_systems.next_body_serial`), updated in the same transaction as the insert. Bulk fills use generation rank and touch the counter once per sector, so there is no per-object round trip and no cross-worker dependency (one worker owns a sector).
- Seed independence: IDs no longer depend on the galaxy seed, only on generation order, which the seed already fixes. Two galaxies share IDs at the same address, so a cross-galaxy reference needs a galaxy number in front; the design notes already give each galaxy its own database.
- Migration: one Alembic revision that adds the two counters and rewrites `uid` values; widths do not change. Boss's no-compatibility rule applies, and GEN.39 already calls for a fresh galaxy.
- Compute drops from 1.3 µs to about 0.1 µs per ID and `_UidIssuer` becomes a counter. The hash-passing code in `uid.py`, `_assign_*_uids` and the rank recount can be removed.

Optional later: match the sector prefix to the fill order (a Hilbert ordinal of the sector) to get the last 5 to 10% of insert speed and the smaller index (415 vs 463 bytes a row). It costs the readable designation; I would not do it unless inserts are still the bottleneck after PERF.45 to 47.

## 7. What I could not measure

- Real fill order, 4 parallel workers, and contention on a shared counter. The model fills one sector at a time in Hilbert order on one connection.
- Anything near the real table sizes. 2M rows of 250-byte payload against a 256 MB pool is the I/O-bound regime, but the real tables are 100 to 1000 times larger. Single runs; differences under about 10% are noise.
- MySQL 8.4 and MariaDB 11.4 (the CI legs). Only MariaDB 10.11 ran.
- The end-to-end failure when an admin adds a system to a sector after a body delete and add. Defect 1 is reproduced at the `assign_uids` call; the sector path is by code reading.
- The 45 bodies per system is one dense run; the true average over the galaxy will differ.
- The sector-level positional clash rate with real Hill-sphere spacing, and the GEN.64 collision count (I used the 0 of 16,526 already recorded in `object-ids.md`).
- Whether any user data (bookmarks, wiki links, facilities) holds row ids or uids that a rewrite would break. Boss's no-compatibility rule says that is fine.

## 8. Sources

- Gaia DR3 `source_id`: HEALPix level-12 cell in bits 36 to 63, processing centre code, then a running number, kept for the source's life although stars move ([Gaia DR3 datamodel](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html), [G-VO note](https://blog.g-vo.org/healpix-maps-in-general-and-in-gaia.html)).
- UUID versions and database-key locality: [RFC 9562](https://www.rfc-editor.org/rfc/rfc9562.html).
- Random keys, page splits and half-full pages in InnoDB: [Percona, 2019](https://www.percona.com/blog/uuids-are-popular-but-bad-for-performance-lets-discuss/), [Percona, 2015](https://www.percona.com/blog/illustrating-primary-key-models-in-innodb-and-their-impact-on-disk-usage/). The search summaries gave the 50% and 94% fill figures from a single vendor blog; my own numbers in section 3.3 are the ones to rely on.
- Repo: `galaxy/uid.py`, `names/object_id.py`, `db/store.py` (`_UidIssuer`, `assign_uids`, `_claim_object_ids`, `replace_system_content`, `add_system_to_sector`), `db/edits.py`, `docs/design/object-ids.md`, `generation-determinism.md`, `orbital-updates.md`, `collisions-and-mergers.md`, `/mnt/project-files/notes/gen68-object-ids.md`.
