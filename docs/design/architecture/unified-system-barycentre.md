# One System rule: a central object, members that can be systems themselves, and the barycentre as the origin

**Status:** plan (decided). Written after 12:03 AM EDT on 2026-10-18 from `main` at 6df95ad and extended after 12:23, 12:24 and 12:25 AM EDT (sections 10 to 13); approved by Boss on the cards S1 to S9 of section 7 (12:10 to 12:15 AM EDT) and S10 to S18 of sections 12.3 and 13.6 (12:30 to 12:32 AM EDT on 2026-10-18), as the coordinator relayed, each the recommended option. Nothing is built. It adds a section to the unified mass-object plan ([unified-star-like-objects.md](unified-star-like-objects.md), [mass-object-class.md](mass-object-class.md), [unified-unit-quantity.md](unified-unit-quantity.md)). The provisional labels are filed as: SY1 = text edits to GEN.552, GEN.315 and GEN.264 (no new item), SY2 = TEST.241, SY3 = GEN.576, SY4 = DB.78, SY5 = GEN.577, SY6 = GEN.578 (expanded in section 12.1), SY7 = MAP.253, SY8 = GEN.579, SY9 = GEN.580, SY10 = TEST.242, SY11 = DOC.71, SY12 = GEN.581, SY13 = GEN.582, SY14 = GEN.583, SY15 = TEST.244, SY17 = GEN.584; SY16 is a note on the 8.3 item GEN.197 (no work item yet).

**Boss, 12:03 AM EDT (project chat):** "a System should be unified rules around the central object for the system (the same way we do for binary stars now) and each object in the system other than the central object has the ability to be another system where for that system the object is the new central object, down to the level of a small object which do not get the capability to have satellites. Since all objects in the system should follow the same physics, the center of the space will be the gravitational center, even though the 'central object' might move around the 3D space of the system, it doesn't have to be in the origin of the coordinate plane, the origin of the plane should be the gravitational center."

Evidence tags: **[M]** measured in the repository (count, grep, line read), **[C]** computed here by arithmetic, **[E]** my estimate, **[R]** recalled and unverified (a lane must verify before use). Nothing is copied from the Solar System or the Milky Way; the examples are generic mass ratios.

## Decided (Boss's card taps, as the coordinator relayed)

S1 to S9 (12:10 to 12:15 AM EDT on 2026-10-18), all on the recommended option: a member's stored orbit attaches to its (sub-system's) barycentre about the parent's central object; the central object's offset is derived on read; the central object is the most massive member and comparable-mass pairs stay binary nodes; a node takes a central child and any number of members; the no-satellite line is a labelled stand-in size cut of 1 km plus the fit and stability checks; moons can have satellites through the same gate plus a tidal-survival gate row; the point class's system frame is the top-level barycentre plus a local frame; stored values change in the 8.4 revision; a System is a view.

S10 to S18 (12:30 to 12:32 AM EDT on 2026-10-18), all on the recommended option: S10 the original values are the origin values (birth mass, age, [Fe/H], composition, formation distance) plus the orbit-geometry draws, and the present state is recomputed about the new primary; S11 the heaviest is found by the total mass of each branch (a sub-system counts as the sum of its parts), ties under 1e-9 relative keep the draw order; S12 when the heaviest changes in simulated time the step runs once at the event (merger, capture, ejection), state-preserving, so no body jumps; S13 trajectories are an analytic two-body orbit per frame with mu including the member's own mass, and the outer system's tide is a stability gate (Hill fraction, triple and Kozai screens); S14 a frame request names any of the named frames or any other point, found through the lowest common ancestor; S15 a point has a live parent link with a versioned cache and the copied anchors and the three `carry_*` calls go; S16 the point holds a kinematic record (position, velocity, spin axis and rate, mass with mu, epoch and next-due) with identity, radius and age beside it in the core; S17 spin is a unit axis plus an angular rate in rad/s; S18 speed, axial tilt, galactic orbital period and phase, rotation period and a* become derived and their columns are dropped in the 8.4 revision with DB.73.

## 1. Short answer

Yes, and the binary-star code already does most of it. The multi-star tree (`physics/system_tree.py`) is exactly this rule for two children: the origin is the pair's barycentre, each child sits either side of it in proportion to the other's mass, and a child can be another node. What does not follow it is everything below the stars: planets and moons are stored relative to their parent's **centre**, the central object sits at the origin with its motion kept in separate "wobble" columns, a wide binary puts the first star at the origin, and a close binary puts the barycentre there. Three conventions, one rule wanted.

The design in one paragraph: **a System is not a new stored thing. It is every body that has a parent link to one host, taken together, with the host as central object. The only stored point is the top-level system's barycentre (as in a scatter row); every other position is stored relative to its parent, as the rows hold it, and the barycentre and the central object's offset from it are derived on read, never stored.** The offset of the central object is one formula, `-sum(M_i r_i) / M_total`, which for one member is the binary formula already in the code. A member can be a body or a sub-system; a sub-system's stored orbit is the orbit of its barycentre about the parent's central object. Whether a body can host a system is decided by a gate (a size line, a stability check, an age check), not by a fixed depth.

Boss's refinements of 12:23 and 12:24 AM EDT (heaviest object becomes primary before validation, everything recomputed from original values; the point object converts between frames on demand) are in sections 10 to 12, and Boss's third (12:25 AM EDT, the point also stores motion and spin) in section 13, which also gathers the three in one table.

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
| Close binary | the pair's barycentre; stars at `binary_primary_position_*` and `binary_secondary_position_*` (km) | `physics/body_positions.py` (built; `web/maps/systemscene._star_around` does not call it yet, MAP.253 wires it) |
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
- **Gone when derived:** `stars.reflex_offset_*` and `planets.reflex_offset_*` (6 columns), `binary_primary/secondary/mutual_position_*` (9 columns), the correlated SQL updates for the wobble in `advance_orbital_phases`, and the `orbit_nodes` copy for nodes that become parent links. Their values become reads of one function. Nothing is lost: a read of a star's position gives the same number the wobble columns held, to the precision of the sum.

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
2. **GEN.315 and GEN.264 (Core 3):** keep the interface host-neutral as the item already asks, and in addition return each body's state **relative to its parent link** plus the node, and compute offsets through **one** function (`physics/system_tree.member_offsets_au`, generalised so a node can hold a central child and several members, n = 1 being the pair of the code). Also keep the `orbit_nodes` copy as it is; it is the same table the later system view reads.
3. **Stored-position convention (a rule for every 8.1 item that stores a position):** never store an absolute sector or galactic position for anything inside a system; store the parent-relative state and let the barycentre be the system's row point. This is already how the scatter rows and the planets are stored, so the rule is a statement, not work. The 8.1 items it protects: GEN.317, GEN.130, GEN.327, GEN.273.
4. **The reflex columns and the three binary position triples:** leave them alone in 8.1; they go in the 8.4 revision with DB.73 and DB.67 to DB.71 (one revision, order agreed by Cores 1 and 4). No new code reads them in 8.1.
5. **No frame change in `SpatialPosition3D` in 8.1.** The "local" frame and the barycentre system frame arrive in 8.4 with GEN.565.

That list changes two existing item texts and adds one convention; it adds no 8.1 work item.

## 6. Items, for the TODO thread to number

Provisional labels SY1 to SY11. Version and lane are proposals; the TODO thread decides.

- **SY1 (8.1 phase 3, no behaviour change, no epoch): edit GEN.552 and GEN.315/264 as in section 5 points 1 and 2.** Core 4 for GEN.552, Core 3 for GEN.315. Cost: text edits; the generalised `member_offsets_au` is Core 3's file.
- **SY2 (8.1, convention): the stored-position rule of section 5 point 3**, written into `docs/design/` and the header of `physics/position.py`, with a test that lists every table that stores a body position and fails if one stores an absolute position of a body inside a system. Core 1 (it owns `point_columns.py`). **Built (TEST.241):** the rule is in the header of `physics/position.py`, `db/point_columns.py` names the tables with a parent link (`BODY_PARENT_COLUMNS`) and the four that keep a free galactic point with their reason and way out (`ABSOLUTE_WITH_A_PARENT`: bright_stars, black_holes, neutron_stars, facilities), and `tests/test_parent_relative_positions.py` fails if a table with a parent link stores its point in another frame.
- **SY3 (8.4): `system(host)`, the one function** that groups by the parent link and returns the barycentre sums, offsets and valid-system check. Core 4, with the MassObject core (GEN.565). Depends on SY1, GEN.565.
- **SY4 (8.4): central-object offsets replace the wobble.** Remove `reflex_offset_*` (6 columns), `calculate_reflex_offset`, and the SQL mirror in `advance_orbital_phases`; `binary_*_position_*` (9 columns) read from `system(host)`. Part of the DB.73 revision; bumps the epoch. Core 1 for the storage, Core 4 for the physics.
- **SY5 (8.4): stored planet and moon states become the state of the barycentre of the body and its own members** (the rule of S1). A re-expression of existing rows from their moons; values change in the last digits, so it is one epoch with SY4. Core 4.
- **SY6 (8.4): the point class gets the local frame and the barycentre system frame** (section 4.3) and loses "nearest star" as the system anchor. Depends on SY3, the unit class. Core 4 or Core 5.
- **SY7 (8.4): map and scene.** `systemscene._star_around`, `body_positions.positions_at`, the browser twin `orbitpositions.js`, and the `systemmap.py` branches use the parent chain; the three conventions of section 3 go. Core 2 for the browser, Core 4 for Python. The pair must stay in step (`tests/test_js_unit.py` checks them against each other).
- **SY8 (8.4): the hosting gate.** `can_host(body)` in three parts (section 4.4) and the gated table of the satellite share by size class, as labelled stand-in rows. Core 4. Depends on GEN.552 and GEN.472 (built).

**As built (GEN.761).** The moon draw reads the gate. `planets.generate_moons` asks `can_host(hosting_host(planet))` first and returns for a host it refuses (the planets of the generator are on circular orbits, a gas or ice giant is fluid). `Planet._draw_moons` takes the chance that a planet has moons from `hosting_gate.satellite_share` where the table has rows for its population and size, else the even chance it had. Only a class D body (50 to 500 km in radius, the generator's small bodies) has a population: `planets.hosting_population` calls it `inner-small-body` in the hot or habitable zone, `belt-small-body` in the cold zone short of `HOSTING_TRANS_NEPTUNIAN_MIN_SNOW_LINES` (10) snow lines of the star, and `trans-neptunian-cold` or `trans-neptunian-hot` beyond, split at `HOSTING_COLD_MAX_INCLINATION_DEG` (5 degrees, the cut behind the table's 22% and 5.5%). A class D body is 100 to 1000 km across, so it takes the large-body rows: 5% inner or belt, 22% cold, 5.5% hot, against the even chance before. A directed moon count (`moon_count`, the `MOONS` directive) skips the chance and still goes through the gate. Moons never host: the draw builds no moon of a moon, and the table has no row for one. The two new knobs are stand-ins (kind P), bracketed in `tuning.py`.

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

**S5. Where is the line for bodies that cannot have satellites?** (class D in the code is a small airless body with at least 1 percent ice, intended scale 50 to 500 km radius, not 500 m, GEN.1091; the Hill zone never sets the line, section 4.4.)
1. A stand-in size cut (a tunable, labelled placeholder, set at a diameter of 1 km until the verified law replaces it), plus the fit and stability checks. **Recommended.** Ships with the first build, is labelled, and the research item replaces the cut with a sourced law.
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

## 10. The primary-first rule (Boss, 12:23 AM EDT)

**Boss:** "Any 2 objects that orbit each other will follow the same rule. First off, structurally, before a system (star system, moon system, etc...) before that system is validated, the heaviest object in the system needs to become the primary object. All other objects in the system are recalculated based on this new reality based on their original values. The center for any system is the barycentre is always the origin, wherever it is. The central object is always the heaviest and so will be closest to the center, but all positions will be updated and trajectories updated based on the gravitational gradient in their frame. So objects orbiting a planet don't need to worry about following the star, they follow the planet and the planet follows the star because we're using relative coordinates."

### 10.1 What the code does already for stars **[M]**

- The scatter tree builder ranks every pair by mass: `multistar._Builder.place` (line 620) swaps `left` and `right` when `self.mass(left) < self.mass(right)`, so **slot 0 is always the heavier branch**, and the branch mass is the sum of its leaves (a sub-system counts as one body).
- The same swap is written three more times: `generation/system.py:299-308` (the fill's two-star path, log line "Binary primary swapped"), `generation/wide_binary.py:119-121` and `generation/binary.py` (the proxy's `__init__` and its `from_dict`, "Same heavier star is the proxy's primary convention"). Four copies of one rule; GEN.316 retires the two-star classes and should leave only the tree's copy.
- Nothing does it for planets or moons. `_validate_mass` checks a planet's mass against its class range, never against its host. A planet heavier than its star cannot occur in the code only because planet classes top out far below the lightest star-like object; it can occur with brown dwarf hosts, double planets (GEN.332), ejected pairs and a free planet's moons (the unified chain), which is exactly where Boss's rule is needed.
- Validation is `StarSystem.validate_system` (`system.py:1837`, calls `validation.space_orbits`): spacing by mutual Hill radii, with `reconcile_zone_and_class` **regenerating** a moved planet's class and everything derived. That recalculation is lossy: it throws the old class away and redraws.

### 10.2 The step, before validation

1. **Collect** the members with their original values (10.3).
2. **Rank** the members by total branch mass (a sub-system counts as the sum of its parts). The heaviest is the central object. A tie (relative difference under 1e-9) keeps the draw order, so the step is stable and never flips back and forth.
3. **Nest** by distance: members that sat inside the stable zone of the demoted body stay with it as a sub-system (they then orbit it, and it orbits the new central object); members outside the pair's orbit become members of the new central object's system, around the pair's barycentre; members in between are in the unstable gap, and validation moves or drops them (the outermost-first rule that already exists).
4. **Recompute** everything that depends on the primary from the original values (10.3), never from the first pass's results.
5. **Validate** (stability, spacing, Hill and Kozai screens) as in the code.

**Two-body case.** Re-primarying a pair changes nothing physical **[C]**: the relative orbit is the same ellipse seen from the other body (relative position is negated, the phase moves by 180 degrees; semi-major axis, eccentricity, inclination, period and the barycentre are unchanged), and each body stays where it was. So the swap for a pair is a relabel plus a sign, and a test can assert that a pair's absolute positions are identical before and after.

### 10.3 What "original values" must be kept so the recalculation is lossless

The principle: **separate what a body was born with and where it was drawn from what is true of it at the present epoch about a particular primary.** The first kind is kept; the second is a function of it and of the structure, computed after the structure is final.

| Kept (origin and geometry) | Recomputed from it (present state about the primary) |
|---|---|
| Seed and birth values: initial mass, age, [Fe/H] | Present mass and radius (evolution) |
| Composition, drawn first (decided in the composition model) | Class, zone, temperature, tidal locking, spin, atmosphere retention |
| Orbit geometry: semi-major axis, eccentricity, inclination, node, argument of periapsis, phase at the epoch | Period, speed, state vector (position and velocity), Hill radius, minimum orbit distance, stability limits |
| The first host's identity and luminosity at formation (for composition and snow-line origin) | Insolation and habitable-zone placement about the new primary |

Against the code: `Planet` keeps distance, inclination and node but **drops the phase** ("a phase drawn here is where the body starts around it, and is not kept: the body is its state vector", `physics/planets.py:807`) and orbits are circular (no eccentricity column on `planets` or `moons`), so for planets and moons the geometry is recoverable from the state vector and the old primary's mass (`orbit_from_vector`, built). So no new column is needed **if the step runs in the same build, with the draw record held in memory until `validate_system` finishes**; the stored result stays the state vector (GEN.289). The scatter tree already keeps full elements per node. Because every body draws from its own seed stream, a body can also be rebuilt from its seed against a new host; that is the mechanism for the class and composition of a body whose zone changed (card S10).

### 10.4 When the heaviest changes

| Case | What happens |
|---|---|
| Equal masses (a twin pair, GEN.332's double planets) | Tie rule: draw order stays; the barycentre is midway and both sit at equal distance; nothing special. |
| At generation: a companion heavier than the first-assigned primary | The step above. For stars this is the existing swap; for planets, moons and wanderers it is new. |
| A white dwarf primary lighter than its main-sequence companion | Handled in the code by relabelling (`system.py:299`); the rule evaluates masses at the system's age, so it stays correct. |
| Simulated time (8.3): mass changes between events | None: masses are constant between events, so the heaviest cannot change; no per-step check is needed. |
| Simulated time: a merger, a capture, an ejection | The step runs once at the event, **state-preserving**: no body moves, relative states are re-expressed about the new central object, elements recomputed from them. This is the `reparent` operation (11.3) also needed by GEN.291 (a bound planet becoming free). |
| The central object is ejected | The remaining members re-rank; the heaviest remaining becomes central, again without any body moving. |

**A correction to a literal reading.** "The central object is always the heaviest and so will be closest to the centre" holds for two bodies (the heavier is nearer the barycentre by the mass ratio) and for each pair of branches at a node. For one central object and several members it does not always hold **[C]**: the central object is displaced from the barycentre by up to sum(m_i a_i)/M, so a light member on a tighter orbit than that is nearer the barycentre than the star. Example: central mass 1, a member of 1e-3 at 5 AU displaces the central object by 0.005 AU, and a second member at 0.001 AU is nearer the barycentre than the central object. I would not write a test that asserts "closest"; the tested property is that the central object is the heaviest and the sum of mass-weighted offsets is zero.

### 10.5 "Trajectories updated by the gravitational gradient in their frame"

What the code can do, and what I read Boss to mean:

- **Storage is already as he describes.** A moon's row has no reference to the star; it is stored relative to its planet, and the planet relative to the star. Nothing a moon does needs the star's motion in the galaxy; the star's motion is the top barycentre's, one stored point.
- **Each orbit is a two-body Kepler orbit in its parent's frame,** with `mu = G (M_central + M_member)` for the member's total mass (a sub-system counts as one body). The universal-variable step (`kepler.universal_step`) advances a stored state from its epoch to any time analytically, so the Kepler step of a stored orbit never integrates (the orbital update's vector change goes through the one core function for every object: Boss, 2026-10-18, 5:30 to 5:33 PM EDT, orbital-update-check-2026-10-18.md Decided 2, which amends card S13): a body's state at time t is read as `advance(relative state, mu, t - epoch)`, once per level of the chain, lazily, nothing written but the top barycentre's point.
- **The "gradient" of the outer system on an inner one** (the star's tide on a moon's orbit, the outer bodies' pull on an inner pair) is not a force in the Kepler step of this model (the vector change of the orbital update goes through the one core function, amending card S13, Boss, 2026-10-18, 5:30 to 5:33 PM EDT). It is a **validity condition**: the inner orbit must lie inside the stable zone (0.4895 of the Hill radius, `is_triple_stable`, the Kozai screen). Within it the two-body orbit is a good model; outside it the pair is not a system and validation moves it. First-order effects (precession, secular exchange) are the research item of card S13.

### 10.6 Where the step lives

In the host chain, not in `system_tree`: `system_tree` stays the pure reader of a tree (offsets from nodes). The step is the first action of validation in the pair builder generalised to any masses (GEN.555, Core 3) and in `populate(host)` (GEN.556, Core 4), with one function `primary_first(members)` that the scatter tree's `place`, the fill, and the generation of moons and free pairs all call. `place` already ranks by branch mass, so it becomes a caller of the function instead of holding its own `if left < right`.

## 11. Frames through the point-in-space object (Boss, 12:24 AM EDT)

**Boss:** "using the point in space object we then update the relative coordinates, the point in space is aware of the center's position in wider space and can then provide converted coordinates to different frames on demand."

### 11.1 What the point class does in the code **[M]**

`SpatialPosition3D` (714 lines) holds three frames (galactic, sector, system) and knows the centre's position through **copies**: each body keeps `_sector_center`, `_star_center` and `_star_velocity` as absolute galactic metres. When a star or a sector moves, every body has to be told: `carry_anchors`, `carry_star_center`, `carry_sector_center` (three call sites, `galaxy/system_position.py:118-119` and `galaxy/sector.py:869`). `.spatial` is used 37 times in 5 files. A body's copy of its star's centre is a copy that can go stale, and it can name only one level (the nearest star), so a moon's frame about its planet does not exist.

### 11.2 The design

- **A point has a parent link to its parent's point, not a copy of anchors.** A point holds its truth relative to the parent (x, y, z, unit, velocity, epoch), exactly as the rows are stored. The point of a top-level system holds the barycentre in the sector and, by the sector's centre, in the galaxy. Moving a parent's point moves every descendant by construction; the three `carry_*` calls go.
- **Each system point caches the central object's offset** `c_off = -sum(M_i rel_i) / M_total`, versioned: any change to a member (add, remove, mass, state) bumps the system point's version, and everything derived from it is recomputed on first read (the lazy `_sync` and `_fresh` pattern of the class, extended up the chain).
- **A hop** is the point's position relative to its parent's barycentre: `c_off(parent) + rel(point)`, with `rel` of the central object equal to zero. A body's coordinates in any ancestor's frame are the sum of the hops from the body up to that ancestor, summed **leaf first**, so a moon's coordinates in its planet's frame use only small numbers and never touch the large galactic ones (the precision rule of the unit note).
- **Requests:**
  - `point.coordinates(frame, form)` with `frame` one of the named frames (`galactic`, `sector`, `top` = the top-level system's barycentre, `local` = the parent's barycentre) or **any other point** (an ancestor, or any body: found through the lowest common ancestor; up to the common ancestor and down again, differences of sums); `form` as in the class (Cartesian, cylindrical, spherical).
  - `point.velocity(frame)` the same way (Galilean sum; all frames are non-rotating and share the galactic axes, as the code's sector frame already does).
  - `point.set(frame, form, values, unit)` as in the class for any frame and form: the value is converted to parent-relative truth by subtracting the chain, so setting a position from a galactic or planet-frame value stores the same parent-relative number.
  - `point.reparent(new_parent)` keeps the absolute position and velocity and re-expresses the relative state about the new parent: the operation of 10.4, of GEN.291's ejection, of a capture and of a merger.
  - `point.sector_address` is derived at the top only and read through the chain (not recomputed per body).
- **Cost:** a chain is at most a handful of hops (a moon of a planet of a star of a pair is four); each hop is a cached vector add; a read costs a few adds after the first, and a change costs one version bump.

### 11.3 How it ties to the unit design

Positions and velocities are `Quantity`-declared fields: the truth at each level has its own declared native unit (km or AU below the top, mpc at the top), the registry converts at each hop, and a request names the unit it wants (`coordinates("galactic", "cartesian", unit="pc")`). The point class's single `length_unit_m` per position is replaced by the per-field unit; `Field` declarations (UN3) cover the new fields. The sums stay in floats in native units (card Q2: plain floats in kernels).

### 11.4 Frames at a glance

| Frame | Origin | Derived how |
|---|---|---|
| galactic | the galactic centre | top barycentre in the galaxy + the chain, summed leaf first, added last |
| sector | the sector's centre | the same minus the sector centre |
| top | the top-level system's barycentre | sum of hops up to the top |
| local | the immediate parent's barycentre | one hop |
| any point P | P's position (centre of P's own frame) | up to the lowest common ancestor of body and P, and down to P |

## 12. What this adds

### 12.1 To the items (SY labels continue from section 6; Boss's S1 to S9 answers stand)

- **SY12 (8.4): `primary_first(members)`** as in 10.2 to 10.6: one function in the host chain; the four copies of the swap (tree `place`, `system.py:299`, `wide_binary.py`, `binary.py`) retire to it (GEN.316 already retires two of them). Core 3 (pair builder, GEN.555) with Core 4 (`populate`, GEN.556). Depends on SY1 and GEN.555.
- **SY13 (8.4): origin and present-state split with the draw record held until validation** (10.3): the `Planet`, `Moon` and `Star` fields marked origin or state in the `Field` table of GEN.565, the draw record held in memory until `validate_system` finishes, `reconcile_zone_and_class` re-run from the seed instead of regenerating. Core 4. Depends on GEN.565, UN3.
- **SY14 (8.4): `reparent`**, the state-preserving operation shared by re-primary at an event, GEN.291, capture and merger. Core 4. Depends on SY6 (the point class).
- **SY6 expanded (8.4): the point class gets the parent link, versioned caches, `coordinates(frame, form)` with any point as a frame, `velocity(frame)`, `set(frame, ...)`, and loses the copied anchors and the three `carry_*` entry points** (3 call sites, 37 `.spatial` uses in 5 files). Core 4 or Core 5; the largest single change in this note, about +200 and -150 lines **[E]**.
- **SY15 (tests):** a pair re-primaried gives identical absolute positions; the step from the originals is idempotent and independent of the order members are listed in; `reparent` preserves absolute position and velocity; frame round trips (body to planet frame to galactic and back) agree to the spacing of a double at the top; the mass-weighted offsets sum to zero at every level; the 10 m body 10,000 AU from its parent keeps its metre. With SY3.
- **SY16 (8.3 note, no work item yet):** the orbital update reads a body's state as `advance(relative state, mu, t - epoch)` along the chain, lazily, writes only the top barycentre's point, and runs `primary_first` and `reparent` only at events (merger, capture, ejection). Goes into the 8.3 plan text when that plan is written.

### 12.2 To the 8.1 foundation (nothing here is 8.1 work; these are the conventions that keep 8.1 from being redone)

1. **GEN.315, GEN.264, GEN.316:** the fill must build stars from the stored members **already ordered heaviest first at every node** (slot 0 = heavier branch, as the tree builder writes it) and must not add a fifth copy of the swap; GEN.316's retirement of the two-star classes removes the three old copies. State this in the items.
2. **GEN.552 (host registry):** rows express an orbiter's mass as a **ratio to the host's mass** where the source gives a ratio, so the draw cannot make an orbiter heavier than its host except through the pair builder, which is where the primary-first step lives. This is a wording edit on a table being written anyway.
3. **GEN.555 (pair builder, 8.4 phase 1) and GEN.556 (`populate`):** the primary-first step is the first action of their validation; their Done lines add "a pair re-primaried gives identical positions".
4. **Stored-position convention (SY2, 8.1):** a stored state is always parent-relative **about the parent's central object, in a non-rotating frame whose axes are the galactic axes**, with its epoch. This is how the rows are stored; stating it in SY2 keeps 8.1 items from storing a state about a barycentre or in a rotating frame.
5. **Nothing in 8.1 touches `SpatialPosition3D`'s anchors.** The parent-link change (SY6) is 8.4, with GEN.565; until then the three `carry_*` calls stay.

### 12.3 Cards

**S10. What counts as the "original values" that survive a re-primary?**
1. The origin values (birth mass, age, [Fe/H], composition, formation distance) and the orbit geometry draws are kept; everything about the present state is recomputed about the new primary. **Recommended.** Matches "recalculated from their original values", and composition stays the one drawn first.
2. All draws are re-run from each body's seed against the new primary, so a body's class and composition can change. The most consistent with the new host; the bodies change character after a swap.
3. Only masses and geometry survive; class and everything derived is regenerated, as `reconcile_zone_and_class` does in the code. Smallest change; the lossy behaviour stays.

**S11. What is compared to find the heaviest?**
1. The total mass of each branch (a sub-system counts as the sum of its parts), ties (under 1e-9 relative) keep the draw order. **Recommended.** The same comparison the tree builder makes.
2. The member's own mass only. A star with a heavy planet system and a slightly lighter star beside it compare by the stars alone.
3. The central body of each branch only. Simplest to state; ignores the mass a sub-system carries.

**S12. What happens when the heaviest changes during simulated time (8.3)?**
1. The step runs once at the event (merger, capture, ejection), state-preserving, so no body jumps. **Recommended.**
2. Never after generation: the generated primary stays primary whatever happens. Simplest; a merger product can then be central while lighter than a member.
3. Re-checked at every update step. Needs no events; masses do not change between events, so the check finds nothing and costs a pass over every system.

**S13. How far does "trajectories follow the gradient in their frame" go?**
1. Analytic two-body orbit per frame with mu including the member's own mass, and the outer system's tide as a stability condition (the Hill fraction, the triple and Kozai screens). **Recommended.** What the code does; the orbital update stays analytic.
2. Add first-order secular effects (precession and exchange between inner and outer orbits) as a research item for later. Still analytic; needs a sourced model and a verification lane.
3. Integrate each system numerically on demand. Closest to the physics; expensive, and the update stops being analytic.

**S14. What can a frame request name?**
1. The named frames (galactic, sector, top, local) or **any other point**, found through the lowest common ancestor. **Recommended.** "A body in its planet's frame" and "the planet in the star's frame" are the same call.
2. The named frames and ancestors only. Simpler to write; a sibling in another sibling's frame needs two calls.
3. The four named frames only. Smallest change; no planet frame for a moon.

**S15. How does a point know its centre's position?**
1. A live parent link to the parent's point; anchors are not copied. **Recommended.** Cannot go stale; the three `carry_*` calls go.
2. Keep the copied anchors and push updates with `carry_*` when a parent moves, extended to every level. Smaller class change; every level of a moved system must be touched, and a copy can go stale.
3. Links below the top, copied anchors at the top (sector and galactic). A halfway form with two mechanisms.

### 12.4 Cost and risk **[E]**

SY6 and SY12 to SY15 together: about +500 to 700 new lines with tests, -300 to 450 removed (four swap copies, three `carry_*`, the lossy reconcile path, scene branches already counted in SY7); net about +100 to +400. The risk is the width of the point class change (37 uses, 5 files) and the values: re-expressing planet states about barycentres changes low digits (already one epoch with SY4 and SY5). The 8.1 conventions of 12.2 cost nothing and prevent a second pass.

## 13. Motion and spin in the point (Boss, 12:25 AM EDT) and the three refinements together

**Boss:** "and point in space also stores the vector of motion as well as spin and other attributes as part of the object with mass classification."

### 13.1 The three refinements in one place

| Boss, AM EDT on 2026-10-18 | Rule | Where in this note |
|---|---|---|
| 12:23 | Any two objects that orbit each other follow one rule; before a system is validated, the heaviest becomes primary and the rest are recomputed from their original values; the barycentre is the origin; coordinates are relative, so a moon follows its planet and the planet follows the star. | Section 10 |
| 12:24 | The point-in-space object updates the relative coordinates and, aware of the centre's position in wider space, gives converted coordinates for any frame on demand. | Section 11 |
| 12:25 | The same point object stores the vector of motion, the spin and other attributes of the object-with-mass classification. | This section |

One sentence for all three: **a mass object is a core of identity and a kinematic record (position, velocity, spin, mass, epoch) held by a point that links to its parent's point; the heaviest member is the central object; every number is stored relative to the parent and converted by the chain on demand.**

### 13.2 What the code stores **[M]**

Motion (`velocity_{x,y,z}_kms`): 15 tables (stars, planets, moons, comets, belts, rogue planets, black holes, neutron stars, interstellar comets, asteroid fields, bright stars, scatter rows, star systems, remnants, nebulae), 144 lines in 27 Python files. The frame differs: the layout of `point_columns.py` says galactic axes for scatter and free objects, and "relative to the star for a `system` frame table" for planets and moons. Beside the vector, `orbital_speed_kms` is stored as well (4 tables, 124 lines in 31 files) although it is the length of the vector.

Spin: a unit axis (`spin_axis_x/y/z`, 8 tables: stars, planets, moons, comets, rogue planets, black holes, neutron stars, interstellar comets) and a rate stored **three ways**: `rotation_period_hours` (6 tables), `spin_period_ms` (neutron stars), and a dimensionless `spin` a* (black holes); plus `axial_tilt_deg` (8 tables), which is the angle between the axis and the orbit normal and so a function of two things already stored. Generated by `physics/spin.py` (281 lines, GEN.104: "a unit axis ... and a rate, so the spin vector is `axis * 2 pi / period`", right-handed, in "the same frame as its orbit", the galaxy's pole as the reference for a star). Spin is **set on the body** (`set_spin`, `spin_values`), not on `SpatialPosition3D`, which holds position, velocity, mass, mu and epoch only (`rotation_period_hours` 52 lines in 11 files; `spin_axis` 39 lines in 5; `axial_tilt` 23 in 5).

Galactic motion: `galactic_orbital_period_gy` and `galactic_orbital_phase_deg` on 8 tables each (16 columns), 178 lines in 18 files, written by the orbital update; computed from position and velocity by `galaxy/galactic_orbit.py` (124 lines).

### 13.3 What belongs in the core's point

The kinematic record of every mass object, all relative to its parent link and in the shared non-rotating axes (the galactic axes, as the code's frames already are):

| Field | Form | Frame and unit | Note |
|---|---|---|---|
| Position | 3 numbers | relative to the parent's central object, km or AU below the top, mpc at the top | section 4.1 |
| Velocity | 3 numbers | the same frame, km/s | stored once; the speed is its length |
| Spin axis | unit vector | the shared axes (a spin is the same in every non-rotating frame) | as stored |
| Spin rate | 1 number, rad/s, zero when not rotating | intrinsic, frame independent | replaces the three period forms |
| Mass | 1 number | kg or solar, declared per column | `mu` derived, as in the class |
| Epoch and next-due | 2 numbers | Unix seconds | as stored |

Not in the point, but in the core next to it (hybrid E1): identity, name, parent, age, radius. Not stored at all, derived on read from the point: **speed**, **axial tilt** (the angle between the spin axis and the orbit normal `r x v` of the parent orbit, `spin.obliquity_between`, built), **rotation period** (2 pi / rate, shown through the time ladder in hours, milliseconds or days as fits), **galactic orbital period and phase** (`galactic_orbit.py` from the top point), a black hole's **a\*** (from the angular momentum and `G M^2 / c`), and a **tidally locked** body's rate (equal to its orbital mean motion, so it follows the primary if the primary changes). A neutron star's `magnetic_field_gauss`, a planet's magnetic moment and the like are **profile** data, not core.

Not in the code, and flagged rather than invented: a moment of inertia and oblateness (needed only for spin precession and tidal evolution, which no code models); the core leaves room for them in the planetary and stellar profiles.

### 13.4 How motion and spin chain through the frames

- **Velocity.** The velocity of a point relative to its parent's barycentre is `v_hop = v_c(parent) + v_rel(point)`, where `v_c(parent) = -sum(M_i v_rel_i) / M_total` is the central object's velocity relative to the barycentre and `v_rel` of the central object is zero. A body's velocity in any ancestor's frame is the sum of the hops up to that ancestor, leaf first, and the galactic velocity adds the top point's velocity last. **A moon's velocity in the star's frame** is `v_c(S) + v_rel(planet system) + v_c(P) + v_rel(moon)`: through the planet, with the frame's own motion included at each hop. Every hop first advances its relative state from its own epoch to the requested time (`universal_step`), so the epochs need not match.
- **Spin.** A spin is intrinsic: the same axis and rate in every non-rotating frame, so a frame change is the identity; no rotation matrices are needed because every frame shares the galactic axes. What does depend on the frame is the **obliquity** and **locking**: both are derived against the parent's orbit, and so follow a re-primary automatically (the orbit normal `r x v` is unchanged when a pair is re-primaried, 10.2).
- **Angular momentum.** The spin angular momentum of a body is intrinsic and contributes nothing to a system's barycentre; the orbital angular momentum of a system about its barycentre (sum of `r x m v`) and its total linear momentum are conserved by `reparent` and by the primary-first step, which gives two exact tests (SY15).
- **Setting.** `point.set(frame, form, values, unit)` accepts a velocity in any frame like a position, and subtracts the chain's velocity to store the relative one; spin is set in the shared axes.

### 13.5 What this adds

**To U4a (GEN.565, the core).** The core is the identity fields, the hybrid's parent link and the point, and the point is the kinematic record of 13.3. Concretely: the 15 velocity triples, 8 spin triples, 6 rotation periods, `spin_period_ms`, `spin` and the 16 galactic orbit columns collapse to one set of core columns (D1 option 2, the core table plus profile tables), and `orbital_speed_kms` (4 tables), `axial_tilt_deg` (8) and the galactic orbit columns (16) are dropped as derived. About 50 table-column declarations go **[E]**, and the same names in `store.py`, `models.py` and `query.py` with them.

**To the unit design.** Velocity (km/s) and angular rate (rad/s) are `Quantity` fields; the registry gets the kinds **angular rate** (unit rad/s; period forms shown through the time ladder) and **angular momentum** (kg m^2/s), with dimensions; and the relations module (section 11 of the unit note) gets `rotation_period <-> angular_rate`, `obliquity(axis, orbit_normal)`, `breakup_rate(mass, radius)` and `galactic_orbit(position, velocity)`.

**To the SY items.**

- **SY6 (the point class)** also holds spin axis and rate and gives `velocity(frame)` through the same chain; spin is set and read in the shared axes. The same class change as section 11, about +80 lines.
- **SY17 (8.4): the core's kinematic record** as in 13.3: GEN.565 declares the fields; the derived columns are dropped in the DB.73 revision with the wobble columns (SY4); `spin.set_spin` and `spin_values` move onto the point. Core 4 with Core 1 for the storage. Depends on SY6, GEN.565.
- **SY15 (tests)** gains: total linear momentum and orbital angular momentum about the barycentre are unchanged by `primary_first` and by `reparent`; a moon's velocity in the star's frame equals the sum of its hops; a locked body's rate follows a primary change.
- **SY16 (8.3 note)** gains: the orbital update advances the relative state of every body along the chain and writes nothing but the top point; spin rate changes (tidal evolution, magnetic braking) are not modelled in the update and stay a later research item.

**To the 8.1 foundation (conventions only).** Add one line to SY2: a stored velocity is, like a stored position, relative to the parent's central object in the shared axes; a stored spin is a unit axis in the shared axes and a rate in rad/s. Nothing in 8.1 stores a spin or velocity otherwise; the 8.1 items that draw spin (GEN.104 is built) or velocity (GEN.142 is superseded by GEN.557) stay as they are.

### 13.6 Cards

**S16. What does the point object hold besides position?**
1. The kinematic record: velocity, spin (axis and rate), mass with `mu`, epoch and next-due; identity, radius and age sit in the core beside it. **Recommended.** Matches "motion, spin and other attributes" without making the point the whole object.
2. Position and velocity only (as the class does plus the parent link); spin stays on the body profile. Smallest class change; spin stays outside the one place that handles frames.
3. The whole core in the point (including radius, age, name). One object; the point class then knows things that are not kinematics.

**S17. How is spin stored?**
1. A unit axis and an angular rate in rad/s (zero when not rotating); the period in hours or milliseconds is derived for display. **Recommended.** One form for a star, a planet, a neutron star and a comet; a black hole's a\* and a locked body's rate are derived.
2. The angular velocity vector (3 numbers in rad/s). Fewer conversions in the maths; the axis and the rate are no longer separately visible.
3. As stored: the axis plus a period in hours (neutron stars in milliseconds, black holes as a\*). No change; three forms of one thing stay.

**S18. Which stored columns become derived?**
1. Speed, axial tilt, galactic orbital period and phase, rotation period and a\* are read from the point (cached per version), and their columns are dropped in the DB.73 revision. **Recommended.** Removes about 50 declarations and the chance of a derived value falling out of step.
2. Keep them as cached columns written with the state. Cheap reads; they can fall out of step and the orbital update has to rewrite them.
3. Drop speed and tilt, keep the galactic orbit columns (they advance with time and cost a rotation-curve evaluation to derive). A halfway rule.
