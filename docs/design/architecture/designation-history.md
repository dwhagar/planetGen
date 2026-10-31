# Designation history: names and IDs may change, the past is kept (GEN.1121)

Design note, written first because GEN.1121 is lock work (a migration and an API field). Decision: Boss, 5:28 PM EDT 2026-10-30, confirmed 5:49 PM EDT: an object may change its name and its ID, it stores the history of its previous designations, and its current ID is always completely unique. A moved rogue keeps its birth-position ID, so a move writes no history row. Region study context: [region-object-store.md](../systems/region-object-store.md) section 7.3. Evidence tags: [S] read in the code, [E] estimate, not measured.

## 1. Short answers

1. **One new table, `designations`, holds only what an object has stopped being called** (a retired name or a retired ID). Current designations stay where they are (`name` and `uid` columns), so an object that never changes costs no row. This differs from the region note ("backfill one current row per object"): see section 4.1, a question for Boss.
2. **An ID is never reused, by anyone, ever.** A retired ID stays reserved in `designations` after its object changes ID or is deleted; a link to it resolves to the current object, or says the object is gone.
3. **A name is unique among current holders and kept in history.** The collision bump (`_rank_chain`) hands a stored object's old name to a different object by design (it is how names stay independent of save order), so a retired name cannot be reserved forever. An old-name lookup prefers the current holder. Section 3.
4. **Members follow without rows of their own.** A belt member is `<belt ID>.<rank hex>`, a computed designation of the parent's current ID; a body's ID embeds its system's ID, so a system's ID change writes one retired row per body in the same transaction.
5. **A renumbering is an alias, not an in-place rewrite.** The old ID goes to `designations`, the new one is issued from the sector's counters (they only move up), and every reader of the old one is redirected.

## 2. What the code does today [S]

- **Two identifiers per object.** `id` is an auto-increment row number, never shown; `uid` is the 96-bit `BINARY(12)` (birth sector 40 bits, serial 28, body number 12) with a unique index per table (`uq_*_uid` in `db/models.py`) and `object_index` (primary key `(kind, uid)`) across kinds. A sector's ID is its address, never stored and never changes ([object-ids.md](object-ids.md)).
- **A row that has a `uid` keeps it** (`replace_system_content`); counters in `id_counters` only move up and `planetgen reset` keeps them. Today nothing changes an existing `uid`.
- **A body's ID embeds its system's.** Planets, moons, stars, belts and comets take the system's sector and serial and a body number. A system whose ID changes therefore changes every body's.
- **Belt members** print `<belt ID>.<rank hex>` (`object_uid.designation`, GEN.1065); a captured member prints `<planet ID>.<rank hex>`. Neither is stored.
- **Names change today, with no memory.** `store.rename_star_system`, `rename_star` and `rename_body` (reached by `PATCH /api/systems|stars|planets|moons/<id>`) overwrite `name`; planets and moons named after their system are renamed with it (`_rename_bodies_with_prefix`). `store.name_in_use` keeps sector, system, star and the uniquely named phenomena unique (`UNIQUE_NAME_TABLES`); a planet or moon renamed by hand may share another body's name (Boss, 2026-09-30).
- **Names move on their own.** `_rank_chain` (GEN.307) hands out the Greek or Roman decorations in the order of sector addresses; a stored object of a higher address moves up when a lower-address sector saves later. The old name goes to a different object. No record is kept.
- **Lookup by ID** is `api/ids.row_id(conn, kind, printed)` (the `UidConverter` checks only the shape). There is no lookup by a past ID or name.

## 3. The rules

**Kinds.** A designation is either a *name* or an *uid*. IDs and names are different spaces; nothing here asks that a name differ from an ID.

| Rule | IDs | Names |
|---|---|---|
| Unique among current holders | Yes (the existing unique indexes) | Yes, in `UNIQUE_NAME_TABLES` only (today's rule); bodies may share names |
| Reserved after retirement | **Forever, for every object** | **No**: a retired name may be given to another object by the collision rule only |
| Kept in history | Yes | Yes, for every object that has a name column |
| Old form resolves | To the current holder; to "gone" if deleted | To the current holder of that name if there is one, else to the latest past holder, with the others listed |

Invariants (each becomes a test):

1. A new ID equals no current ID and no retired ID of any object. It is issued from the counters; the insert into `designations` fails loudly on a duplicate (a UNIQUE key on `(designation_kind, designation)` for IDs).
2. An object changes designation in one transaction: the retired value is written to `designations` with the epoch and the reason, then the column changes.
3. A retired ID never comes back, not even to the object that held it (a revert issues a new one). A retired name may return to its object, which keeps its row, so the history shows every change in order.
4. A change of ID cascades: every body whose ID embeds the changed one gets a retired row and a new ID in the same transaction, children after the parent, and `object_index` is rewritten from the same data.
5. A belt member's or captured body's designation is derived from the parent's current ID and never stored. An old `<old belt ID>.<rank>` resolves by splitting at the point, resolving the parent through `designations`, and keeping the rank.
6. Reasons are a closed list: `collision bump`, `rename`, `regeneration`, `parent change`, `admin renumber`. A moved rogue is not among them (birth ID fixed).
7. A deleted object keeps its retired IDs and gets a tombstone row (`gone_epoch`); a link to it answers 410 with the history, never 404 and never another object. The current ID of a deleted object is retired at the same time, so it is never reissued.
8. Sector IDs never change (an address); a sector's name can.
9. Fingerprints (GEN.58, GEN.1129) digest current designations only; the history digest is separate, so a rename does not read as a different object.

## 4. Schema

`designations` (one row per retired designation):

| Column | Meaning |
|---|---|
| `id` | auto-increment |
| `designation_kind` | `uid` or `name` |
| `designation` | the retired text: the printed ID (`BINARY(12)` for IDs, held as `uid`, unique) or the name |
| `object_kind`, `object_row_id` | the object that held it (the table and its `id`) |
| `valid_from_epoch`, `retired_epoch` | when it began (NULL when it was the birth designation) and when it stopped |
| `reason` | one of the closed list |
| `gone_epoch` | NULL while the object exists, the epoch it was deleted in otherwise |

Keys: unique on `uid` for IDs (the reservation); index on `(designation_kind, designation)` for lookup; index on `(object_kind, object_row_id, retired_epoch)` for the page's "formerly" lines. The key is the row number, not the uid, because the uid is the thing that changes.

### 4.1 Why not a row for every current designation (a question for Boss)

The region note proposed one current row per existing object. I recommend against it: every body of the galaxy would add a row and two indexes for something its own table already says (the current `uid` is unique already), and the cost grows with the galaxy while the changes are rare. The retired-only table answers every question the rules ask: reservation (is it in `designations` or in a live `uid` column), redirect, "formerly" lines, and the history digest. **Option A (recommended): retired-only.** **Option B: a current row per object**, as the region note says, simpler to read and heavier. [E: sizes not measured; a row of this table is about 60 to 80 bytes plus two indexes. The count of bodies in a full galaxy is the number to ask the DB lane for before picking B.]

### 4.2 Migration

Testing mode: new epoch, Alembic revision, wiped database allowed. The table is created empty and nothing is backfilled under Option A; under Option B one `INSERT ... SELECT` per table in chunks. Number allocated at merge by `allocate_numbers`; TODO confirms before it is written (standing rule).

## 5. Behaviour

- **Writers.** `store.rename_star_system`, `rename_star`, `rename_body`, `_rank_chain` (reason `collision bump`) and a new `designations.renumber(conn, kind, row_id, reason)` for the edit flow of 8.6 write the history; there is no other writer. `renumber` is the only place an ID changes and nothing calls it until 8.6 (tests call it).
- **Lookup.** `api/ids.row_id` falls back to `designations` when the printed ID is not current. A page route answers 301 to the current URL; a `GET` API call answers 301 with the new `Location` and `X-Designation-Formerly`; a `PATCH` or `DELETE` on a retired ID answers 409 with the current ID and does nothing (a write never acts on a guess). A deleted object answers 410. Search by value (8.2) matches a past name and shows "formerly".
- **Object page.** A "Formerly" line lists previous designations newest first with the reason and the epoch's date. The JSON of an object gains `former_designations`: a list of `{designation, kind, reason, retired_epoch}` (additive, API version bump).
- **Cost.** Lookup is one indexed read, only on a miss of the live column. A renumber of a system touches its bodies (tens of rows) once.

## 6. Proposed split (TODO assigns owners; smallest parts)

| Part | Content | Lock |
|---|---|---|
| a | This note | none |
| b | `designations` table, `galaxy/designations.py` (rules, `retire`, `renumber`), the three rename writers and `_rank_chain` write history | migration, epoch if names move |
| c | Redirect and 410 in `api/ids.py` and the page routes | none after b |
| d | "Formerly" lines on the object page, `former_designations`, search by past name | API field |
| e | Cascade to bodies, belt members and captured bodies; tombstones on delete | none after b |
| f | Tests (section 3 invariants) with b to e, and the 8.6 edit flow's use | none |

## 7. Questions for Boss

1. **Option A (retired-only table, recommended) or B (a current row per object)?**
2. **Names are not reserved after retirement** (the collision bump gives a stored object's name to another object, which is how names stay independent of save order). OK, with IDs reserved forever? If you want names reserved too, the collision rule must pick a different name instead of moving the stored one, and every generated name could change once.
3. **A write against a retired ID answers 409, a read redirects.** OK?

## 8. Not verified

- The number of bodies in a full galaxy and so the real size of Option B.
- Whether `object_index` is rewritten in the same transaction by every existing rename path (read, not run).
- Whether any admin tool deletes rows outside `store.delete_*`, which would skip the tombstone (GEN.1123's audit lists the delete paths).

## 9. Built (GEN.1138)

The table is migration 0197 (number allocated at merge). It differs from section 4 in one way: the retired text is two nullable columns instead of one, `designation_uid` (`BINARY(12)`, unique, the reservation) for an ID and `designation_name` (`VARCHAR(255)`, plain index with `retired_epoch`) for a name, and a CHECK keeps exactly the one that matches `designation_kind`. `retired_at` is the wall-clock time of the change. The reasons are a CHECK of the closed list. `valid_from_epoch` is filled from the epoch the object's previous designation of the same kind was retired in, and is NULL for the birth designation.

Writers are in `db/designation_rows.py`, which `store.py` calls from every place that changes a `name` of a named table: `set_name` keeps the old name and then changes the row (a rename by a person is `rename`; the Greek and Roman decoration, the collision bump of a numbered name and the healing of a system name are `collision bump`; the bodies and comets that follow a system or star are kept under the same reason as the change that moved them; the placeholder names of a system built again keep nothing). `renumber` gives one object a new ID, refusing any ID that is live or retired, keeping the old one and moving the object's `object_index` row and the constellation and lock pointers at it. The pointers with no index (`object_index.parent_uid`, `orbital_events.containing_uid`, `sectors.max_mass_id`) and the bodies of a renumbered system are the cascade of GEN.1141.

Read in the code for section 8: no existing rename path rewrites `object_index`; the index follows the rename on its next refresh (`modified_at`), so the name in the index can lag a rename until then. The old name is searchable only after GEN.1140.

## 10. Built (GEN.1140)

`designations.former_json(rows)` is the one wording of an object's history: `[{designation, kind, reason, retired_epoch, retired_on}]`, newest first (an ID printed, a name as text; `retired_on` is the date of the `retired_at` timestamp). `designation_rows.former_designations(conn, table, row_id)` reads it and `api/ids.former_designations(conn, kind, row_id)` takes the kind of a reference. Four answers carry it as `former_designations` (additive, so the API version stays): `GET /api/objects/<ref>` (every kind, and so the page of a star, planet, moon, belt, comet or facility), `GET /api/systems/<id>`, `GET /api/sectors/<id>` and `GET /api/phenomena/<type>/<id>`; the 301/409/410 bodies of GEN.1139 use the same list.

The page line is `web/templates/partials/formerly.html`, fed by `web/formerly.py` (the reasons in words): "Formerly: Old Name (name, renamed, epoch 130, 2026-10-31); ...". It sits under the badges of the item, system, sector and phenomenon pages and is absent for an object that never changed. The epoch has no calendar date of its own, so the date is the one the change was written.

Search: each of the sector, system, star, planet and moon name searches also matches a row of `designations` that kept a name containing the words typed (`db/query._search_name`; the table is small and holds only names that were given up), and a result row found that way carries `formerly` (`_with_formerly`); the Search page writes "(formerly X)" after its name. Only objects that still exist match (`gone_epoch IS NULL`). A past ID is found by `GET /api/lookup?text=` (`api/lookup.py`, the rules of `designations.resolve_uid` and `resolve_name`), which the Search page asks when the box holds an Object ID or a name that no object has now: it writes "X is now <link>" or "X belonged to an object that has been deleted". `asteroid_belts` left `NAMED_TABLES`: a belt has no name column.
