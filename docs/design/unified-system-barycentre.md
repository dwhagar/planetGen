# One System rule: a central object, members that can be systems themselves, and the barycentre as the origin

**Status:** plan (decided). Written after 12:03 AM EDT on 2026-10-18 from `main` at 6df95ad; approved by Boss on the nine cards S1 to S9 of section 7 (12:10 to 12:15 AM EDT on 2026-10-18, as the coordinator relayed), each the recommended option. Nothing is built. It adds a section to the unified mass-object plan ([unified-star-like-objects.md](unified-star-like-objects.md), [mass-object-class.md](mass-object-class.md), [unified-unit-quantity.md](unified-unit-quantity.md)). The provisional labels of section 6 are filed as: SY1 = text edits to GEN.552, GEN.315 and GEN.264 (no new item), SY2 = TEST.241, SY3 = GEN.576, SY4 = DB.78, SY5 = GEN.577, SY6 = GEN.578, SY7 = MAP.253, SY8 = GEN.579, SY9 = GEN.580, SY10 = TEST.242, SY11 = DOC.71. **Decided (Boss, 12:10 to 12:15 AM EDT, 2026-10-18):** S1, S2, S4, S5, S6, S7, S8 and S9 tapped between 12:10 and 12:13 AM EDT on 2026-10-18, S3 at 12:15 AM EDT ('Most massive, pairs stay binary'), every one the recommended option: S1 a member's stored orbit is its barycentre about the central object (a planet's stored state is the barycentre of the planet and its moons), S2 the central object's offset is derived on read and never stored (the wobble columns go), S3 the most massive member is the central object and comparable-mass pairs stay binary nodes, S4 a system node takes a central child and any number of members, S5 a stand-in size cut of 1 km diameter for hosting satellites plus fit and stability checks, S6 moons may have satellites through the same gate with a tidal-survival gate row first (the odds start at zero; this supersedes 'moons are the lowest level'), S7 the point class's system frame measures from the top-level barycentre with a new local frame from the immediate parent's barycentre, S8 the stored values change in the 8.4 revision with DB.73 and DB.67 to DB.71, S9 a System is a view (rows with a common parent link, cached sums).

**Boss, 12:03 AM EDT (project chat):** "a System should be unified rules around the central object for the system (the same way we do for binary stars now) and each object in the system other than the central object has the ability to be another system where for that system the object is the new central object, down to the level of a small object which do not get the capability to have satellites. Since all objects in the system should follow the same physics, the center of the space will be the gravitational center, even though the 'central object' might move around the 3D space of the system, it doesn't have to be in the origin of the coordinate plane, the origin of the plane should be the gravitational center."

Evidence tags: **[M]** measured in the repository (count, grep, line read), **[C]** computed here by arithmetic, **[E]** my estimate, **[R]** recalled and unverified (a lane must verify before use). Nothing is copied from the Solar System or the Milky Way; the examples are generic mass ratios.

## 1. Short answer

Yes, and the binary-star code already does most of it. The multi-star tree (`physics/system_tree.py`) is exactly this rule for two children: the origin is the pair's barycentre, each child sits either side of it in proportion to the other's mass, and a child can be another node. What does not follow it is everything below the stars: planets and moons are stored relative to their parent's **centre**, the central object sits at the origin with its motion kept in separate "wobble" columns, a wide binary puts the first star at the origin, and a close binary puts the barycentre there. Three conventions, one rule wanted.

The design in one paragraph: **a System is not a new stored thing. It is every body that has a parent link to one host, taken together, with the host as central object. The only stored point is the top-level system's barycentre (as in a scatter row); every other position is stored relative to its parent, as it is now, and the barycentre and the central object's offset from it are derived on read, never stored.** The offset of the central object is one formula, `-sum(M_i r_i) / M_total`, which for one member is the binary formula already in the code. A member can be a body or a sub-system; a sub-system's stored orbit is the orbit of its barycentre about the parent's central object. Whether a body can host a system is decided by a gate (a size line, a stability check, an age check), not by a fixed depth.

Foundation for 8.1: five conventions and two small edits to items already filed (GEN.552, GEN.315), no code of the new rule itself (section 5). The stored-value change (planet states become barycentre states) rides the 8.4 revision with DB.73 and one epoch.

## 2. The rule in plain terms

- **System:** one host plus every body whose parent is that host. **Central object:** the host. **Member:** a body that orbits the central object, or a sub-system (a member that has members of its own).
- **A sub-system acts as one body for the system above it.** Its mass is the sum of its parts; the thing that orbits the central object is its barycentre.
- **Origin:** a system's coordinates have their origin at the system's barycentre, B = (M_c r_c + sum M_i r_i) / M_total. Positions inside it are measured from B.
- **Positions stored:** each member's state (position and velocity) is stored **relative to the central object** (for a sub-system, of its barycentre), as the scatter tree and the planet rows already do.
- **Positions derived:** with `rel_i` the stored relative position of member i, the central object is at `r_c = -sum(M_i rel_i) / M_total` and member i is at `r_i = r_c + rel_i`, both measured from B. For one member this is `r_c = -m/(M+m) rel`, the close-pair formula in `system_tree.py:16`.
- **Recursion:** the position of a body in the galaxy is the top system's barycentre plus the sum of the offsets down its chain of parents. No level needs to know about the levels above.
- **End of the recursion:** a body that fails the hosting gate (section 4.4) has no members and is a leaf.

Size of the effect **[C]**, generic ratios: a member of 1e-3 of the central mass at 5 AU moves the central object 0.005 AU (747,000 km), about a radius of a star of 7e5 km; a moon of 1e-2 of its planet at 400,000 km moves the planet 3,960 km. Real offsets, small against the orbit, large against the body: which is why they must be exact in the position and not a separate "wobble".

## 3. What the code does

Four conventions for the origin **[M]**:

| Case | Where the origin is | Where it is set |
|---|---|---|
| Single star with planets | the star's centre; its motion is a separate `reflex_offset_x/y/z_km` on `stars` and `planets` (schema v20), computed by `orbits.calculate_reflex_offset` (pairwise pulls summed) and mirrored in SQL by `store.advance_orbital_phases` | `generation/system.py:432-452` |
| Close binary | the pair's barycentre; stars at `binary_primary_position_*` and `binary_secondary_position_*` (km) | `physics/body_positions.py` via `web/maps/systemscene._star_around` |
| Wide binary | the first star; the second at `binary_mutual_position_*` | `systemscene._stars` |
| Scattered multi-star tree | the barycentre of each node; children either side in proportion to the other's mass; nodes nest | `physics/system_tree.py` (144 lines), scatter rows only |

Planets and moons (`parent_type` in star, planet or system; `parent_id`; `host_node_id` for a circumbinary planet): the state vector is relative to the parent's **centre**, and a moon is added to its planet's position (`body_positions.positions_at`). A planet's mass is never in any barycentre except through the reflex columns.

Counts **[M]**:

- Files: `physics/system_tree.py` 144 lines; `physics/body_positions.py` 116; its browser twin `html/static/orbitpositions.js` 314; `web/maps/systemscene.py` 194 (`_star_around`, three branches); `web/maps/systemmap.py` 1,853 (9 barycentre mentions); `generation/binary.py` 478 and `wide_binary.py` 371 (7 and 4 mentions); `physics/orbits.py` 510 (the reflex function, about 50 lines); `physics/insolation.py` 96 (the circumbinary centroid); `generation/system.py` 2,213; `db/store.py` 9,646.
- Tables: `star_systems` 71 columns, 27 orbit or position columns, 9 of them `binary_{primary,secondary,mutual}_position_{x,y,z}_km`; `stars` 53 (3 reflex, `parent_node_id` never written); `planets` 118 (18 orbit or position columns, 3 reflex, `host_node_id`); `moons` 115 (17); `comets` 36 (22); `scatter_nodes` 18 and `orbit_nodes` 21 (`orbit_nodes` never written).
- Names across `src/planetgen`: `parent_type` 18 lines in 7 files; `parent_id` 49 in 10; `binary_primary_position` 32 in 6; `binary_secondary_position` 32 in 6; `binary_mutual_position` 55 in 8; `secondary_mass_fraction` 40 in 7; `host_node_id` 5 in 3; `scatter_nodes` 23 in 8; `orbit_nodes` 13 in 4; `star_id` 246 in 18.
- The point class: `SpatialPosition3D` has three frames, galactic, sector and system, and the system frame is "offset from the center of the body's nearest star" (`physics/position.py`, module docstring), which is a fourth convention for the origin.

What the docstring of `calculate_reflex_offset` says is worth keeping: a lot of physics (insolation, Hill sphere, tidal locking) depends on the stored relative vector being the **true separation**, not a barycentre-reduced one. The new rule keeps that: members stay stored relative to the central object, so those consumers do not change.

## 4. The design

### 4.1 What is stored and what is derived

- **Stored truth:** the state (position, velocity, epoch) of each body relative to its parent link, in the unit and frame the column declares (the point-in-space rule: kept as given). For a top-level system, the barycentre's point, velocity and epoch in the sector grid (the scatter row, the `star_systems` row).
- **Derived:** the barycentre of every sub-system, every central object's offset, and every absolute position, all on first read and cached; marking stale on any change. This is the lazy pattern of `SpatialPosition3D._sync` and `_fresh` (PERF.46), and card Q1 of the unit design already decided it for units.
- **Gone when derived:** `stars.reflex_offset_*` and `planets.reflex_offset_*` (6 columns), `binary_primary/secondary/mutual_position_*` (9 columns), the correlated SQL updates for the wobble in `advance_orbital_phases`, and the `orbit_nodes` copy for nodes that are now the parent links. Their values become reads of one function. Nothing is lost: a read of a star's position gives the same number the wobble columns held, to the precision of the sum.

### 4.2 Keeping the barycentre in step

| Event | What happens |
|---|---|
| A member is added (a planet, a moon, a sub-system) | The host's cached sum is marked stale; every ancestor system's cache is marked stale up the parent chain. |
| A member is removed (engulfed, ejected, merged) | Same. Its row's parent link is cleared or set to the new host. |
| A body changes parent (a planet ejected becomes free, GEN.291; a capture) | Its stored relative state is **re-expressed** about the new parent at the change, from the two absolute positions, so no position jumps. Old and new hosts both marked stale. |
| A member's mass changes (mass loss, accretion, an evolved star) | The sum is stale; relative orbits are not rescaled by this rule, which leaves adiabatic orbit change to the orbital update (8.3). |
| The system moves in the galaxy | Only the top barycentre's point and velocity change (the nightly update, as in the code); nothing inside is written. |
| A sub-system's own member moves | Only its own parent-relative state changes. The ancestors' sums change only if the sub-system's total mass or barycentre-relative offset does; the barycentre of a sub-system does not move when its members orbit, which is why nothing above needs rewriting. |

The last row is the useful property: the barycentre of a sub-system is conserved under its own internal motion, so a deeper level can be advanced without touching the level above.

### 4.3 How it fits the hybrid MassObject and the point class

- **E1 (approved):** a body has `parent` as data. A System is the group of rows sharing a parent. It needs no class of its own and no table of its own: a `system(host)` function groups by the parent link and returns the cached sums. A body that has members gets a cached `SystemSummary` (total mass, sum of mass-weighted relative positions, member count), not a stored profile.
- **Point in space:** `SpatialPosition3D` keeps its three frames and gains the parent-relative truth as an option: the "system" frame of the code becomes "offset from the barycentre of the top-level system of this body" (was the nearest star), and a "local" frame is the offset from the immediate parent's barycentre, from which the others are derived by walking the chain. The point class already does the same walk for sector and galactic frames; this adds one more anchor level. The class stays the one place that knows frames and forms.
- **Units:** the offset sum is a sum of `Quantity` lengths in the declared native units; the sum is done in the unit the parent column declares (km or AU), so a 10 m body 10,000 AU from its parent keeps its metre (precision table in the unit note).
- **Precision for the 10 m floor:** a moon of a 10 m body has an offset of about a body radius or less, so the sums have to be formed from **relative** numbers, never from absolute galactic positions (a double at 8 kpc has a spacing of 32,768 m **[C]**). That is why the stored truth is parent-relative and the only large number is the top barycentre.

### 4.4 The line below which a body has no satellites

I first expected the geometry to set it, and it does not **[C]**. With the project's own stable-zone fraction (`MOON_PROGRADE_STABLE_HILL_FRACTION` = 0.4895 of the Hill radius) and a density of 2000 kg/m3, a body of any size down to the 10 m floor has a stable zone many times its own radius:

| Host diameter | Stable zone at 0.1 AU | at 1 AU | at 40 AU |
|---|---|---|---|
| 10 m | 41 m (4x host + satellite radius) | 410 m (41x) | 16,400 m (1,640x) |
| 100 m | 410 m (7.5x) | 4,100 m (75x) | 164,000 m (2,980x) |
| 10 km | 41 km (8x) | 410 km (82x) | 16,400 km (3,280x) |

(Star of 1 solar mass; the zone scales with the host's size, so the ratio is nearly the same at every size and only grows with distance.) So the Hill criterion alone never forbids a satellite of a 10 m body, and the **line has to come from something else**: the satellite floor (a satellite must be at least 10 m across, so a host of 10 m can only have an equal-size companion), the chance that a body of that size has one at all, and forces the code does not model (radiation torques, which spin small bodies up and shed satellites **[R]**).

So the gate has three parts, in this order:

1. **A satellite must fit:** the smallest satellite is the 10 m floor, and a satellite must not outweigh its host by more than the system's mass-ratio limit (a tunable; equal-mass pairs are allowed only through the node form of the binary tree).
2. **A stable orbit must exist:** inside 0.4895 Hill radius of the host at its pericentre, outside the host radius plus the satellite radius (and outside the Roche limit for a fluid body), and passing the existing hierarchical stability checks (`physics/stability.py`, `is_triple_stable`) when the satellite is comparable in mass.
3. **A chance that it exists at all:** a gated table (GEN.472 format) of the share of bodies of each size class that have one or more satellites. This is where the real line is, and its values are research: observed binary and multiple asteroid fractions are known for small bodies **[R]** (Margot et al. 2002, Science; Pravec et al. 2006, Icarus; the order of 15% of near-Earth asteroids over 0.3 km are binaries, recalled and unverified). Until a lane verifies and sources it, the table ships as a labelled stand-in (Boss's gap-plan decision 3, "ship as placeholders") with a stand-in size cut (card S5), replaced when GEN.566 (size-frequency research) gives the law.

### 4.5 Moons of moons

Boss's earlier rule, "moons are the lowest level", is superseded by this message for any body that passes the gate. A moon's stable zone is inside its planet's Hill sphere and inside the moon's own, so the second is the tighter and the same function decides it; there is no depth count. What changes for the first version is only the odds: the table starts with zero rows for a moon host, and a tidal-survival check (a sub-satellite of a moon is stripped over time by the planet's tides **[R]**, the timescale is a research item) is a gate row that must exist before any are drawn. The recursion ends because the gate says no, not because a level counter does. Needs a card (S6).

### 4.6 Validity of a system

A system is valid when every member's mass is at most the mass-ratio limit of the central object **or** the pair is expressed as a node of the existing binary tree (so a pair of comparable masses is a two-body node, which is exact, and everything else is test-mass-like and combines by superposition). The generator already makes this split for stars; planets and moons are always far below the limit. A check at generation, `valid_system(host)`, covers it.

## 5. What 8.1 needs so nothing is redone

Nothing of the rule's code has to ship in 8.1. What 8.1 has to do is not build anything that this rule would undo:

1. **GEN.552 (host registry, 8.1 phase 3, Core 4):** its text says "Moons are the lowest level, so no host class lists a moon among its orbiters". Change it to: a host class row exists for a moon (and for a small body: a row whose orbiter table is empty by the gate), and the registry carries a `may_host` gate per row. The cost is one row in a table that is being written anyway; the later edit would be a change to a registry already shipped. Also: the recursion in `populate(host, depth)` (section 3.3 of the unified note, "three working depths") ends by the gate, not by a depth of three.
2. **GEN.315 and GEN.264 (Core 3):** keep the interface host-neutral as the item already asks, and in addition return each body's state **relative to its parent link** plus the node, and compute offsets through **one** function (`physics/system_tree.member_offsets_au`, generalised so a node can hold a central child and several members, n = 1 being today's pair). Also keep the `orbit_nodes` copy as it is; it is the same table the later system view reads.
3. **Stored-position convention (a rule for every 8.1 item that stores a position):** never store an absolute sector or galactic position for anything inside a system; store the parent-relative state and let the barycentre be the system's row point. This is already how the scatter rows and the planets are stored, so the rule is a statement, not work. The 8.1 items it protects: GEN.317, GEN.130, GEN.327, GEN.273.
4. **The reflex columns and the three binary position triples:** leave them alone in 8.1; they go in the 8.4 revision with DB.73 and DB.67 to DB.71 (one revision, order agreed by Cores 1 and 4). No new code reads them in 8.1.
5. **No frame change in `SpatialPosition3D` in 8.1.** The "local" frame and the barycentre system frame arrive in 8.4 with GEN.565.

That list changes two existing item texts and adds one convention; it adds no 8.1 work item.

## 6. Items, for the TODO thread to number

Provisional labels SY1 to SY11. Version and lane are proposals; the TODO thread decides.

- **SY1 (8.1 phase 3, no behaviour change, no epoch): edit GEN.552 and GEN.315/264 as in section 5 points 1 and 2.** Core 4 for GEN.552, Core 3 for GEN.315. Cost: text edits; the generalised `member_offsets_au` is Core 3's file.
- **SY2 (8.1, convention): the stored-position rule of section 5 point 3**, written into `docs/design/` and the header of `physics/position.py`, with a test that lists every table that stores a body position and fails if one stores an absolute position of a body inside a system. Core 1 (it owns `point_columns.py`).
- **SY3 (8.4): `system(host)`, the one function** that groups by the parent link and returns the barycentre sums, offsets and valid-system check. Core 4, with the MassObject core (GEN.565). Depends on SY1, GEN.565.
- **SY4 (8.4): central-object offsets replace the wobble.** Remove `reflex_offset_*` (6 columns), `calculate_reflex_offset`, and the SQL mirror in `advance_orbital_phases`; `binary_*_position_*` (9 columns) read from `system(host)`. Part of the DB.73 revision; bumps the epoch. Core 1 for the storage, Core 4 for the physics.
- **SY5 (8.4): stored planet and moon states become the state of the barycentre of the body and its own members** (the rule of S1). A re-expression of existing rows from their moons; values change in the last digits, so it is one epoch with SY4. Core 4.
- **SY6 (8.4): the point class gets the local frame and the barycentre system frame** (section 4.3) and loses "nearest star" as the system anchor. Depends on SY3, the unit class. Core 4 or Core 5.
- **SY7 (8.4): map and scene.** `systemscene._star_around`, `body_positions.positions_at`, the browser twin `orbitpositions.js`, and the `systemmap.py` branches use the parent chain; the three conventions of section 3 go. Core 2 for the browser, Core 4 for Python. The pair must stay in step (`tests/test_js_unit.py` checks them against each other).
- **SY8 (8.4): the hosting gate.** `can_host(body)` in three parts (section 4.4) and the gated table of the satellite share by size class, as labelled stand-in rows. Core 4. Depends on GEN.552 and GEN.472 (built).
- **SY9 (research, can start at no code risk): the size-frequency and binary fraction of small bodies, and the survival of sub-satellites of moons.** Joins GEN.566 (U15a). Verify the [R] items of sections 4.4 and 4.5 first.
- **SY10 (tests): invariants for any system at any time.** `sum(M_i r_i) = 0` about the barycentre; adding a member and removing it returns the other positions exactly; a parent change leaves every absolute position unchanged; a sub-system's barycentre does not move when its members orbit; the single-member case equals the existing close-pair formula bit for bit. With SY3.
- **SY11 (docs): `database-schema.md`, the system view in `docs/design/`, and the header of `system_tree.py`.**

Estimated size **[E]**: SY3 and SY6 about 300 to 500 new lines with tests; SY4, SY5 and SY7 remove 300 to 600 lines (the wobble code, the three scene branches and their JS twin, 15 columns). Net roughly zero to a few hundred lines saved; the value is one rule instead of four, the central object's true position, and sub-satellites without a special case.

## 7. Cards

**S1. What does a member's stored orbit attach to?**
1. The member's barycentre (for a sub-system) about the parent's central object, two-body, as the scatter tree already stores it; a planet's stored state is then the barycentre of the planet and its moons. **Recommended.** One rule at every level; the moon-mass shift is tiny (section 2).
2. The member's own central body about the central object, as planets are stored in the code; the sub-system's barycentre is derived. Smallest change to stored rows; the tree and the planets stay two rules.
3. Everything relative to the top barycentre (N-body style). Truest to "origin is the gravitational centre", but every member's row then changes whenever any sibling changes, and the orbital update (8.3) loses its analytic step.

**S2. Where does the central object's own offset live?**
1. Derived on read from the members, cached, never stored (the point-in-space way). **Recommended.** The wobble columns and the SQL mirror go.
2. Stored in columns and rewritten when members change (the current `reflex_offset_*`). Fast reads; a cache that can fall out of step and needs the SQL mirror.
3. Stored for stars only, derived for planets. Keeps the star positions the map reads; two rules.

**S3. Who is the central object?**
1. The most massive member of the system (the first drawn on a tie); a pair of comparable masses is a two-body node of the existing binary tree. **Recommended.**
2. The generator's placed object (the head of the scatter row), whatever its mass. Simple; wrong when the placed object is the lighter star.
3. Always none: the barycentre only, every member offset (like the existing `group` node), for every system. Uniform; loses the meaning of "orbits the central object" for planets.

The `group` node (a young non-hierarchical system) stays as it is under any option.

**S4. How many members does a system node take?**
1. A central child and any number of members (planets, moons, test-mass-like members), the binary node being the case of one member; comparable masses nest as nodes. **Recommended.**
2. Strictly two children everywhere (a chain of nodes, one per planet). Exactly hierarchical; a system of ten planets becomes ten nested nodes.
3. N members for planets and moons, strictly two for stars, as in the code. No change; keeps two rules.

**S5. Where is the line for bodies that cannot have satellites?** (class D in the code is 50 to 500 km radius, not 500 m; the Hill zone never sets the line, section 4.4.)
1. A stand-in size cut (a tunable, labelled placeholder, set at a diameter of 1 km until the verified law replaces it), plus the fit and stability checks. **Recommended.** Ships now, is labelled, and the research item replaces the cut with a sourced law.
2. A class cut: only class D and larger can host (radius from 50 km). Simplest; wrong for the small binaries known to exist **[R]**.
3. No cut: any body of at least 20 m can host, with the share by size from a table that is zero until the research lands. Needs no number invented, but nothing below class D can have satellites before the research is done, so it behaves as option 2 first.

**S6. Do moons get satellites?** (Boss's new rule says any non-small member can host.)
1. Yes, through the same gate, with a tidal-survival gate row that must exist before any are drawn, so the odds start at zero for a moon host. **Recommended.**
2. Yes, same gate only (the Hill fraction), no survival check. Draws sub-moons the tides would strip in real systems **[R]**.
3. No: moons are leaves, Boss's earlier rule. The recursion stops at three levels again.

**S7. What does the point class's "system" frame measure from?**
1. The barycentre of the top-level system, with a new "local" frame from the immediate parent's barycentre; both derived by walking the parent chain. **Recommended.**
2. The barycentre of the top-level system only; callers derive the local frame themselves. Smaller class; every caller repeats the walk.
3. The nearest star, as in the code. Keeps the fourth convention; contradicts the new rule.

**S8. When do the stored values change?**
1. In the 8.4 revision with DB.73 and DB.67 to DB.71, after the 8.1 conventions of section 5 are in place. **Recommended.** One epoch, one rescatter group.
2. In 8.1, with the host registry. Earlier, but puts a value-changing migration on the 8.1 critical path.
3. Not at all: keep the three conventions and add the new rule for sub-systems only. No migration; the three conventions and the wobble columns stay.

**S9. Is a System an object or a view?**
1. A view: the rows with a common parent link, with cached sums; the only stored system point is the top-level row. **Recommended.**
2. An object with its own row at every level. Easy to query; a table of mostly empty rows and a second place for each parent link.
3. An object only in memory during generation, a view in the database. A hybrid with the same cache but a second code path.

## 8. Limits and risks

- **The dynamics stay Keplerian.** The rule fixes where things are; each orbit still evolves as a two-body orbit about its central object (as in the code), so mutual perturbations between members are not modelled. The superposition of the central object's offset is exact for the positions given and approximate as a dynamical model for comparable masses, which is why comparable pairs are nodes (section 4.6).
- **Mass change moves the barycentre** (section 4.2). Mass loss and evolved stars are the orbital update's problem (8.3), not this rule's.
- **Bit-identical output is not kept** for SY4 and SY5 (the stored states of planets and the star positions change in low digits); they ride one epoch.
- **Precision of the sum:** forming sums from relative numbers, never from absolute positions, is a rule of the implementation; a test with a 10 m body 10,000 AU from its parent covers it.
- **Three [R] items** (binary asteroid fraction, the order of 15% of near-Earth asteroids over 0.3 km; radiation torques shedding satellites from small bodies; tidal stripping of sub-satellites of moons) are recalled and unverified; none is used by any number in the design except as a pointer for SY9 to verify.
- **The Hill-zone table** uses the project's constant 0.4895 (`tuning.MOON_PROGRADE_STABLE_HILL_FRACTION`) and a density of 2000 kg/m3 for every size; it is a geometry check, not a stability proof.

## 9. Evidence

[M] figures were counted from the repository at 6df95ad with grep, line counts and a read of `system_tree.py`, `body_positions.py`, `systemscene.py`, `orbits.py:405-470`, `position.py`, the schema columns of `star_systems`, `stars`, `planets`, `moons`, `comets`, `scatter_nodes` and `orbit_nodes`, and the filed text of GEN.550, GEN.552, GEN.315 and GEN.264; none was run against a database. [C] figures are arithmetic from the stated constants (the Hill table, the two offset examples, the double-precision spacing). [E] figures are the weakest. The design assumes the decisions already taken on the unified plan (hybrid MassObject, host registry, GEN.552 in 8.1 phase 3, the 8.4 revision carrying DB.73 and DB.67 to DB.71) and the unit design (lazy derived forms, parent-relative storage).
