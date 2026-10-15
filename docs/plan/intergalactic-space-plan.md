# Intergalactic space plan

The plan for taking planetGen from a galaxy-shaped volume to a whole sphere of space: every cell out to 262 kpc has defined contents, one rule set tuned by a profile for each place, the halo, the satellites and the globular clusters in place, and a sector fill that builds every kind of object and the systems it can have. Written 2026-10-15 at Boss's request by the TODO tasks thread; [docs/TODO.md](../TODO.md) is the master file and wins where they disagree. Designs: [intergalactic-and-halo-contents.md](../design/intergalactic-and-halo-contents.md) (the sphere, 15 steps) and [phenomena-catalogue-and-fill.md](../design/phenomena-catalogue-and-fill.md) (every kind of object, the fill stages, the system rules). Lane assignments live in the Core, Interface and Quality groups of [execution-plan.html](execution-plan.html).

## 1. Goal

A visitor can ask for a sector anywhere inside the reach (the Magellanic Clouds, the halo, a globular cluster, the gap between galaxies) and get the contents the real place would hold, made by the same rules as the neighbourhood, deterministic from the galaxy seed. Nothing is generated unless it is asked for or scattered; an asked sector is always generated and saved, even when empty. The first neighbourhood pull stays in the galactic region.

## 2. Decided (Boss, 2026-10-15)

1. Reach 262 kpc (65,535 shells) and 96-bit object IDs; the research note's defaults stand (section 6 lists them).
2. One consistent rule set everywhere; each uncharted sector gets tuning parameters (its profile) from its place.
3. Early alpha (5:27 PM EDT): the database schema is not a constraint; redesign tables as the work needs, one migration or a fresh galaxy. The scatter from the origin keeps the banded, layered approach (the bands of a shell have fixed slices; Core checks this against `object_first.py` and `spherical.py` before relying on it).
4. Order (5:27 PM EDT): the spherical sector layout finishes first, then the intergalactic scatters and sector parameters, then the spherical-region database query. If a later stage is not ready when the one before it finishes, skip it and start the next; a stage waiting on research never blocks the others (5:32 PM EDT).
5. Scatters and backscatters carry no planets or moons (5:39 PM EDT). Nebulae need stars and phenomena, not planets. Planets and moons come only when a sector is filled with everything, at real-world probabilities.
6. Stellar companions of stars and star-like objects (stars, former stars, protostars, neutron stars, black holes, including neutron star and black hole pairs) are generated at scatter time, because they matter for gravity (5:51 PM EDT).
7. "All types represented" means audited frequency, never forced: accurate probabilities and no forced representation in a sector, except a kind that landed there in a scatter before generation (5:51 PM EDT).
8. Every nebula gets the stars it needs as real systems; Wolf-Rayet stars and ring nebulae are built (5:23 PM EDT).
9. Multi-star plan approved (6:08 PM EDT): plan and build every observed and theoretical type of multi-star and star-like-object system into the generator throughout, matching observed probabilities best and theoretical ones where observation is silent. The cap of stars per system is 12, not 7. Compute the best path for three build lanes and one TODO lane (6:02 PM EDT): section 11.

## 3. The stages, in build order

```
A  spherical grid ──► B1 halo + profile ──► B2 components + clusters ──► B4 scatters over the sphere
                          │                                                   │
                          ├──► B3 registry + fill machinery ──► B5 real contents (young stars, nebulae, W stars)
                          │
                          └──► B6 Interface (profile lookup, sector page, map)      B7 tests + docs beside each step
M  multi-star systems at scatter time (approved; the catalogue of types is filed)       C  region query (live; its key goes in the A migration)
L  later group: nothing waits on it
```

- **A** is first and alone on the critical path: the shared migration (GEN.230 and GEN.223 together, the 96-bit layout and the point-in-space columns) then the reach (GEN.233). GEN.232 is a bug that must land before any wide scatter.
- **B1 to B7** are the whole sphere. B1 is a strict chain (halo, then profile, then kind factors); after the profile exists, B2, B3 and B6 can run in parallel across lanes. B4 needs the components (GEN.237). B5 needs the registry and fill stages (GEN.251) and the pre-main-sequence state. B7 is written beside each step.
- **M** is the multi-star stage (section 4), partly inside A (the member tables share the migration). **C** is the region query (section 5), live since its research landed (#1659); its key and columns ride in the stage A migration. **L** is the later group (section 8).
- Two Core lanes (Boss, 2026-10-15 evening; the split approved by the TODO lane at 7:00 PM EDT, section 11): Core 1 is the only lane that adds migrations; Core and Core 2 together own B1 to B5, M and C's key; Interface owns B6 and the readers it was given; Quality owns B7 and every TEST and DOC item.

## 4. Stage M: multi-star and star-like systems at scatter time

Approved by Boss (2026-10-15, 5:51 PM EDT, widened 6:08 PM EDT to every observed and theoretical type of multi-star or star-like-object system, built throughout the generator). Design: [multi-star-systems.md](../design/multi-star-systems.md), now with the catalogue of about sixty system types (sections 13 to 16). Twenty items are filed and approved (the defaults of the note stand: compact pairs promoted to the scatter, black hole companion 1.25%, hierarchical trees only, a row's members deleted when its sector builds it, pair rates from the middle of the evidence; the cap is now 12 stars, not 7).

- **Stored form** (GEN.260): a scatter row becomes a system key: barycentre point, total mass, members and orbit nodes in two child tables, one column per value, no JSON. The same migration adds `state`, `flags`, `interaction`, a SMALLINT `slot`, the `group` and `compact` node kinds, and the offset, velocity and decay-age columns of young groups. It is the shared migration 0093 (GEN.223, GEN.230, GEN.260, DB.34, DB.35), so there is one migration.
- **Generator** (GEN.261, GEN.266, GEN.267, GEN.268): `draw_member_tree` (GEN.129's steps a to d plus member kinds and the Peters, Kozai and wide-companion checks), the compact-binary orbit laws, the cap setting and `--stars N` directive, and `SYSTEM_TYPE_RATES`, one entry per catalogue type with a completeness test.
- **Scatter and fill** (GEN.262, GEN.263, GEN.264, GEN.265): companions in the scatter passes with a rate correction so pairs are not double counted, promotion of bound compact pairs, the fill building from stored members (GEN.251's system pass shrinks to planets, moons and belts), and the nightly motion of a system with its total mass.
- **The catalogue of types** (GEN.269 to GEN.275, TEST.164): GEN.269 the star model (Roche classifier, `interaction` on every pair, flags, rejuvenated age, the kinds `stripped-helium`, `thorne-zytkow` and `wolf-rayet`); GEN.270 close binaries (contact, Algol, Be, merger products); GEN.271 white-dwarf binaries (double white dwarfs, cataclysmic variables, AM CVn, symbiotic and barium stars); GEN.272 accreting compact binaries (X-ray binaries, microquasars, spider pulsars); GEN.273 young groups and trapezia; GEN.274 runaway origins and cluster-hosted compact binaries; GEN.275 brown-dwarf, planetary-mass and circumbinary types; TEST.164 the type audit.
- **Interface and Quality**: UX.114 lists members on a scatter row's pages; TEST.162 is the multiplicity audit and TEST.164 the type audit; DOC.50 the docs.
- **Order**: GEN.260 and GEN.265 go in with the migration (stage A); GEN.268 follows the registry (GEN.250) and GEN.261 follows it; GEN.262 and GEN.264 follow once GEN.251 and the member tables exist; GEN.269 follows GEN.261; GEN.270, GEN.271, GEN.274 and GEN.275 are built by the Interface lane (section 11) once their prerequisite lands, GEN.272 and GEN.273 by Core. B1 to B3 do not wait for the rest of M. GEN.240 (scatters over the whole sphere) takes members from GEN.262 once it exists and may start without them.

## 5. Stage C: the spherical-region database query

Live: its research is finished and merged (#1659, [spherical-region-queries.md](../design/spherical-region-queries.md); scripts and results in /mnt/project-files/research/spherical-region-queries/). Boss's sequencing rule (5:32 PM EDT) said to skip it if it was not ready before the spherical layout; it is ready, so it is built. The design: every object gets one `cell_key`, a Morton key of its current position in 0.25 pc cubes (21 bits per axis, plus or minus 262,144 pc, exactly the reach), computed by MySQL as a stored generated column, with the covering index `(cell_key, x, y, z)`. A sphere query lays the sphere over the nested grid, merges the cubes into at most 256 key runs and joins them to the table, so only those runs are read; measured at 50 pc and 10 million rows it takes 42 ms against 705 ms for a plain box. The cell key is not the object ID (the ID is the permanent birth address; the key is where the object is now).

- **Key and columns**: GEN.276 (`planetgen/galaxy/cellkey.py`, numpy only), DB.34 (the column and index on every positioned table, in the shared 0093 migration), DB.35 (the scatter tables, which keep their address columns).
- **Readers**: NAV.60 (what is nearby reads by the key, for every kind, including objects with no sector and unbuilt scatter rows, with a plan-shape test) and NAV.61 (nearest-first by a radius ladder and a dense-core guard).
- **Quality**: TEST.165 (re-measure on a real generated galaxy, warm and cold), TEST.166 (the key follows an object across a cube boundary), DOC.51.
- The 21-bit constants and the outside-the-box guard were not tested in the research; DB.34 and GEN.276 test them.
- Stage C does not block A or B: GEN.276 is a small pure function that goes first, DB.34 and DB.35 ride in the migration, and the readers and tests follow on the Interface and Quality lanes.

## 6. Defaults that stand until Boss says otherwise

Open questions of the design notes; each default is a decision unless changed.

| Question | Default | Item |
|---|---|---|
| Harris catalogue file for the Milky Way clusters | synthetic generator first, the real catalogue when it arrives | GEN.237, GEN.159 |
| Satellites | real Milky Way satellites at their real places for the default shape; seeded fictional ones for other shapes later | GEN.237 |
| The density floor | deleted; the halo and components replace it with no setting | GEN.234 |
| Halo normalisation | the 1.5e9 solar mass anchor, the local density (about 2.7 times higher) kept as the stated alternative | GEN.234 |
| Magellanic mix | 1% young and 9% intermediate | GEN.235 |
| How much of an H II cluster is real | the ionizing members only (5.5 solar masses and up) | GEN.241 |
| Protostar counts and the 5x shortfall in giant molecular clouds | rates as stated; the shortfall accepted for the first release | GEN.241, GEN.258 |
| Planetary-nebula central star companion | 20% (recalled, unsourced) | GEN.251 |
| Black hole companion fraction | 1.25% (the multi-star note; the earlier 20% rested on three Gaia systems) | GEN.261, GEN.130 |
| Intermediate-mass black hole count (about 1e5 in the Galaxy against a handful of candidates) | left alone | GEN.259 (a question for Boss) |
| Unstable wide binaries already generated | count first; the repair needs Boss | GEN.229 (outside this plan) |

## 7. Items by stage

Each row is generated from the item's own entry in TODO.md on 2026-10-15; the entry wins if it has changed since.

### Stage A: The spherical grid finishes first (the shared migration also carries the scatter member tables, the cell key and the catalogue columns)

| Item | What | Owner | Needs | Unlocks |
|---|---|---|---|---|
| GEN.223 | The nightly orbital update moves every object in the galaxy, and each object is stored as a point in space in its own columns (bug) | Core | none | DB.34, DB.35, GEN.233, GEN.260, GEN.265, PERF.103, TEST.166 |
| PERF.103 | A per-night cost ceiling on the orbital update: warn with the backlog size and move the rest the next night | Core | GEN.223 | (leaf) |
| GEN.230 | Sector designations and object-ID birth fields cannot name bands at or beyond 2048 (bug) | Core | none | GEN.233 |
| GEN.260 | Scatter member columns and tables: a scatter row is a system with a barycentre, a total mass and a stored tree of members | Core | GEN.223 | DOC.50, GEN.262, GEN.264, GEN.265, GEN.269, GEN.273, UX.114 |
| GEN.233 | Grid reach: every cell out to 262 kpc has defined contents, and a plan-time reach setting widens the object ID to 96 bits | Core | GEN.230, GEN.223 | DOC.48, GEN.234, GEN.235, GEN.238, GEN.257, TEST.159 |
| TEST.156 | Rewrite the tests that name ring, layer and slot | Quality | none | (leaf) |
| TEST.158 | The math check tests the spherical shell, band and wedge grid | Quality | none | (leaf) |
| DOC.47 | Rebuild docs/design/architecture.html from the new architecture.md | Quality | none | (leaf) |

### Stage B1: The halo and the sector profile

| Item | What | Owner | Needs | Unlocks |
|---|---|---|---|---|
| GEN.234 | A measured stellar halo replaces the flat density floor, with the halo, cluster and satellite populations | Core | GEN.233 | DOC.48, GEN.235, GEN.240, GEN.248 |
| GEN.235 | One rule set for every cell: a sector profile computed from the address, and any asked sector in the reach is generated from it and saved | Core | GEN.233, GEN.234 | API.27, DOC.48, GEN.236, GEN.237, GEN.238, GEN.240, GEN.241, GEN.242, GEN.243, GEN.249, GEN.251, GEN.254, TEST.159, UX.113 |
| GEN.158 | Add a metallicity value to stars | Core | none | GEN.159, GEN.162, GEN.163, GEN.164, GEN.236, GEN.243 |
| GEN.236 | Kind factors in the profile: the rate per star of each rare kind times a factor from the cell, and the planet cull by metallicity | Core | GEN.235, GEN.158 | GEN.241, TEST.159 |

### Stage B2: Components and globular clusters

| Item | What | Owner | Needs | Unlocks |
|---|---|---|---|---|
| GEN.159 | Globular clusters: cluster table, King tables and the Milky Way catalogue | Core | GEN.158 | GEN.160, GEN.163, GEN.164, GEN.237 |
| GEN.160 | Cluster density in the sector gate, with a "cluster" population | Core | GEN.159 | GEN.161, GEN.162, GEN.237 |
| GEN.237 | Satellite and cluster components: one table for the Magellanic Clouds, Sagittarius, dwarf galaxies and globular clusters | Core | GEN.235, GEN.159, GEN.160 | DOC.48, GEN.239, GEN.240, GEN.244, GEN.248, MAP.207 |
| GEN.161 | Bright-first fill for cluster sectors | Core | GEN.160 | (leaf) |
| GEN.162 | Planet cull and blue stragglers in clusters | Core | GEN.158, GEN.160 | (leaf) |
| GEN.164 | Synthetic globular-cluster systems for generated galaxies | Core | GEN.9, GEN.158, GEN.159 | (leaf) |

### Stage B3: The registry and the fill machinery

| Item | What | Owner | Needs | Unlocks |
|---|---|---|---|---|
| GEN.251 | Sector fill stages on their own random streams, and a system pass so every object that can have a system rolls for one | Core | GEN.250, GEN.235 | DOC.49, GEN.130, GEN.241, GEN.252, GEN.264, TEST.160, TEST.161 |
| GEN.130 | Compact-remnant systems at sector fill: companions and pulsar planets for neutron stars and black holes | Core | GEN.251, GEN.266 | DOC.49, GEN.163, GEN.247, GEN.253 |
| GEN.163 | Type-B pulsar planets in globular clusters (GEN.130 follow-on) | Core | GEN.130, GEN.158, GEN.159 | (leaf) |

### Stage M: Multi-star and star-like systems at scatter time (approved by Boss, 2026-10-15 5:51 and 6:08 PM EDT)

| Item | What | Owner | Needs | Unlocks |
|---|---|---|---|---|
| GEN.261 | The member-tree generator: stars, brown dwarfs, remnants and protostars in stable hierarchies up to twelve members | Core | GEN.250, GEN.268 | GEN.262, GEN.264, GEN.266, GEN.267, GEN.269, GEN.273 |
| GEN.265 | Nightly motion of a system: the barycentre moves with the system mass, members are computed from their stored phase on read | Core | GEN.223, GEN.260 | (leaf) |
| GEN.262 | Companions in the scatter passes: the mass, luminosity and phenomenon passes draw members, with a rate correction so pairs are not double counted | Core | GEN.260, GEN.261 | DOC.50, GEN.263, GEN.264, TEST.162 |
| GEN.266 | Compact-binary orbit rules: a black hole or neutron star with a star, tight compact pairs and tertiaries | Core | GEN.261 | DOC.50, GEN.130, GEN.272 |
| GEN.263 | Promote every bound compact pair to the scatter, whatever its mass | Core | GEN.262 | (leaf) |
| GEN.264 | The fill builds a system from its stored members; the system pass keeps only planets, moons and belts | Core | GEN.260, GEN.261, GEN.262, GEN.251 | DOC.50, TEST.162 |
| GEN.267 | A MAX_STARS setting and a --stars N directive, with a prevalence knob for five or more stars | Core | GEN.261 | (leaf) |
| GEN.269 | Interaction states in the star model: a Roche classifier, `interaction` on every pair, member state and flags, rejuvenated age and the new kinds | Core | GEN.260, GEN.261 | GEN.270, GEN.271, GEN.272, TEST.164 |
| GEN.270 | Close-binary types: contact, Algol, spectroscopic and eclipsing flags, tertiary shares, merger products and Be stars | Interface | GEN.269 | TEST.164 |
| GEN.271 | White-dwarf binary types: close double white dwarfs, post-common-envelope pairs, cataclysmic variables, AM CVn, symbiotic and barium stars | Interface | GEN.269 | TEST.164 |
| GEN.272 | Accreting compact binaries: X-ray binaries, accretion luminosity, microquasars, spider pulsars and compact triples | Core | GEN.269, GEN.266 | TEST.164 |
| GEN.273 | Young groups and trapezia: the `group` node, offset and velocity columns, a decay age and the remnant draw | Core | GEN.260, GEN.261 | TEST.164 |
| GEN.274 | Cluster-hosted compact binaries and runaway bookkeeping: `runaway_origin` and the compact companion left by a supernova | Interface | GEN.268 | TEST.164 |
| GEN.275 | Substellar and circumbinary types: brown-dwarf and planetary-mass companions, circumbinary disks and planets | Interface | GEN.268 | TEST.164 |
| TEST.164 | Catalogue type audit: type counts, region gates, the N tail, interaction consistency and group decay | Quality | GEN.268, GEN.269, GEN.270, GEN.271, GEN.272, GEN.273, GEN.274, GEN.275 | DOC.50 |
| UX.114 | Scatter-row pages list the stored members and the map shows a system mass | Interface | GEN.260 | (leaf) |
| TEST.162 | Multiplicity audit: companion states, row bookkeeping, stability, determinism, no planets in a scatter, fill reuse, motion and the census anchors | Quality | GEN.262, GEN.264 | (leaf) |
| DOC.50 | Docs: the schema, the multi-star note and the catalogue say where the members live | Quality, after the code | GEN.260, GEN.262, GEN.264, GEN.266, TEST.164 | (leaf) |

### Stage C: The spherical-region database query (live: research done, #1659)

| Item | What | Owner | Needs | Unlocks |
|---|---|---|---|---|
| GEN.276 | planetgen/galaxy/cellkey.py: the cell key function and the sphere-to-key-runs cover, numpy only | Core | none | DB.34, DOC.51, NAV.60, TEST.165 |
| DB.34 | Cell key and covering index on every table that holds a position, in the shared 0093 migration | Core | GEN.223, GEN.276 | DB.35, DOC.51, NAV.60, TEST.165, TEST.166 |
| DB.35 | Scatter rows keep their address columns and gain the point columns and the cell key | Core | GEN.223, DB.34 | NAV.60 |
| NAV.60 | What is nearby reads by the cell key: one statement for every kind, including objects with no sector and unbuilt scatter rows | Interface | DB.34, GEN.276, DB.35 | DOC.51, NAV.61 |
| NAV.61 | Nearest-first by a radius ladder, and a dense-core guard with an early-stop count | Interface | NAV.60 | (leaf) |
| TEST.165 | Re-measure the region query on a real generated galaxy, warm and cold | Quality | DB.34, GEN.276 | (leaf) |
| TEST.166 | The cell key follows an object across a cube boundary, and the point columns and the key round-trip bit for bit | Quality | DB.34, GEN.223 | (leaf) |
| DOC.51 | Docs: the cell key in the coordinate system and schema docs, the region note moves to built, the API states the cap and the ladder | Quality | DB.34, GEN.276, NAV.60 | (leaf) |

### Stage B4: The scatters fill the whole sphere

| Item | What | Owner | Needs | Unlocks |
|---|---|---|---|---|
| GEN.240 | The three scatters fill the whole sphere with one algorithm: band groups, the phenomena pass on the object-first sampler, and a component pass | Core | GEN.232, GEN.234, GEN.235, GEN.237 | (leaf) |
| GEN.238 | Wanderers: hypervelocity stars reach out to the grid, and kicked neutron stars and black holes are an ejecta population | Core | GEN.233, GEN.232, GEN.235 | MAP.207 |

### Stage B5: Real contents: young stars, nebulae, Wolf-Rayet stars

| Item | What | Owner | Needs | Unlocks |
|---|---|---|---|---|
| GEN.242 | A pre-main-sequence state for young stars, so protostars can exist as real stars | Core | GEN.235 | GEN.241, GEN.245, GEN.255, GEN.258 |
| GEN.241 | Every nebula gets the stars it needs as real systems made with the nebula, inside its bounds | Core | GEN.235, GEN.236, GEN.242, GEN.251 | DOC.48, GEN.243, GEN.244, GEN.258, TEST.159 |
| GEN.243 | Wolf-Rayet stars and ring nebulae: a W state for the last 0.3 Myr of very massive stars, and nebula class X | Core | GEN.235, GEN.158, GEN.241 | DOC.48, GEN.246, GEN.256, TEST.159 |

### Stage B6: Interface: the profile and the map beyond the disk

| Item | What | Owner | Needs | Unlocks |
|---|---|---|---|---|
| API.27 | A read-only profile lookup for any address in the reach, so an uncharted sector can be previewed | Interface asks Core for the profile function | GEN.235 | UX.113 |
| UX.113 | The sector page shows the profile of any address in the reach, and heavy cells and far neighbourhoods warn before they run | Interface | GEN.235, API.27 | (leaf) |
| MAP.207 | The Galaxy Map draws satellites, clusters and wanderers as points beyond the disk outline | Interface | GEN.237, GEN.238 | (leaf) |

### Stage B7: Tests and docs

| Item | What | Owner | Needs | Unlocks |
|---|---|---|---|---|
| TEST.159 | Whole-sphere tests: required stars per nebula kind, Wolf-Rayet counts, tiling at the reach, profile and empty-sector determinism, and the pre-placed share per population | Quality | GEN.233, GEN.235, GEN.236, GEN.241, GEN.243 | (leaf) |
| TEST.160 | Representation audit: seeded sectors in each region class show every kind at its expected frequency, with systems and required contents checked | Quality | GEN.250, GEN.251 | (leaf) |
| TEST.161 | Fill-stage and system-pass tests: stream independence, same seed same system, and the probability table | Quality | GEN.251 | (leaf) |
| DOC.48 | Docs: the reach, the halo, the profile rules and the Wolf-Rayet and nebula-contents rules in the design notes | Quality, after the code | GEN.233, GEN.234, GEN.235, GEN.237, GEN.241, GEN.243 | (leaf) |
| DOC.49 | Docs: point the rate, anomaly, multi-star and nebula notes at the phenomena catalogue | Quality, after the code | GEN.250, GEN.251, GEN.130 | (leaf) |

### Stage L: The later group (nothing waits on it)

| Item | What | Owner | Needs | Unlocks |
|---|---|---|---|---|
| GEN.239 | High-velocity clouds and the Magellanic Stream as placed gas clouds (later) | Core | GEN.237 | (leaf) |
| GEN.244 | Young open clusters and OB associations as components (later) | Core, later | GEN.237, GEN.241 | (leaf) |
| GEN.245 | Debris disks around A and F stars (later) | Core, later | GEN.242 | (leaf) |
| GEN.246 | Luminous blue variables and hypergiants (later) | Core, later | GEN.243 | (leaf) |
| GEN.247 | Accreting white dwarf binaries: novae, cataclysmic and symbiotic variables, Type Ia progenitors (later) | Core, later | GEN.129, GEN.130 | (leaf) |
| GEN.248 | Tidal streams as placed components of the halo (later) | Core, later | GEN.237, GEN.234 | (leaf) |
| GEN.249 | Be stars, carbon stars and Mira variables as star attributes (later) | Core, later | GEN.235 | (leaf) |
| GEN.252 | Brown dwarfs as stars: an L, T and Y state for companions and primaries, with the brown-dwarf desert | Core, later | GEN.251 | (leaf) |
| GEN.253 | White-dwarf pollution flag and cool white dwarf presets | Core, later | GEN.130 | (leaf) |
| GEN.254 | Superbubble field beside the nebula field (the Local Bubble cavity) | Core, later | GEN.235 | (leaf) |
| GEN.255 | Drawn jets: Herbig-Haro flows, microquasar jets and pulsar jets along their spin axes | Core then Interface, later | GEN.242 | (leaf) |
| GEN.256 | A radiation hazard score derived from distance to hot stars, Wolf-Rayet stars, remnants and magnetars | Core, later | GEN.243, GEN.153 | (leaf) |
| GEN.257 | A Fermi-bubble overlay on the Galaxy Map | Interface, later | GEN.233 | (leaf) |
| GEN.258 | Protostars embedded in giant molecular clouds | Core, later | GEN.241, GEN.242 | (leaf) |
| GEN.259 | Intermediate-mass black hole share gives about 100,000 in the Galaxy against a handful of candidates (bug) | Core | none | (leaf) |

## 8. The later group

Nothing in stages A to B7 waits on these; they are planned after the train. The multi-star research may reorder GEN.252, GEN.247 and GEN.253, which share its definitions.

- Sphere: GEN.239 high-velocity clouds and the Magellanic Stream, GEN.248 tidal streams, GEN.244 young open clusters and OB associations.
- Stars and remnants: GEN.246 luminous blue variables and hypergiants, GEN.247 accreting white dwarf binaries, GEN.249 Be, carbon and Mira star attributes, GEN.252 brown dwarfs as stars, GEN.253 white-dwarf pollution, GEN.258 protostars embedded in giant molecular clouds.
- Anomalies split from GEN.114: GEN.254 superbubble field, GEN.255 jets, GEN.256 radiation hazard score, GEN.257 Fermi-bubble overlay; GEN.245 debris disks.
- Skipped on purpose (catalogue section 4.7): subspace anomalies (default off), wormholes and cosmic strings, dark matter nebulae, quark stars, bursts, blanets. Thorne-Zytkow objects and merger products came back with the 6:08 PM EDT catalogue as flagged theoretical types (GEN.269, GEN.270).

## 9. Tests and docs

| Step | Tests (Quality) | Docs (Quality) |
|---|---|---|
| A | TEST.156, TEST.158 (the grid names, the math check at the reach) | DOC.47 |
| B1 to B6 | written beside each build step | none until the code lands |
| B3 | TEST.161 (stream independence, same seed same system, the probability table, no planet or moon on a scatter row) | DOC.49 |
| M | TEST.162 (the multiplicity audit: companion states, stability, determinism, fill reuse, motion, census anchors) | DOC.50 |
| B5, B7 | TEST.159 (required stars per nebula kind, Wolf-Rayet counts, tiling at the reach, profile and empty-sector determinism, pre-placed share per population) | DOC.48 |
| B7 | TEST.160 (the representation audit: expected against observed counts per region class, the fraction with a system, required contents) | DOC.49 |

Gate for every step: a green `scripts/test.py changed`; `full --check` when a stage finishes; the PR names the run. Docs-only PRs merge themselves.

## 10. Keeping this plan true

The TODO tasks thread owns this file. It amends it whenever an item in it is filed, retired or re-ordered, whenever Boss decides one of the section 6 questions, and when the multi-star catalogue extension or the region-query research hands over its items. Merged items are retired from TODO.md and this file with their PR.

## 11. Lane queues and the critical path (computed 2026-10-15, 6:10 PM EDT)

Boss (6:02 PM EDT) approved the plan and asked for the best path for the three build lanes (Core, Interface, Quality) and the TODO lane. The full per-lane queues, with the first five items of each, are in the UPDATE 6:10 PM EDT block of `/mnt/project-files/reports/lane-briefs-three-threads-2026-10-15.md`; this is the summary.

**Method.** All 325 open items with their prerequisites (from TODO.md) form one dependency graph. Each item has a relative effort in slots (1 is about one small PR; my estimates, 2 to 5 for the big whole-sphere items). A simulation lets each lane take the first ready item in its queue; a lane never idles while it has a ready item.

**Findings.**
- The dependency critical path is 8 items, 27 slots: GEN.223, GEN.233, GEN.234, GEN.235, GEN.251, GEN.241, GEN.243, TEST.159. Nothing else is longer; Interface and Quality items hang off it.
- **Core is the bottleneck, not the chain.** Core's own plan queue (stages A, B and M) is 29 items and 84 slots, 3.1 times the critical path, because every plan item needs the generator or a migration. The multi-star items add 23 slots to it. The best path therefore keeps Core on the plan queue from the first minute and moves everything else off it.
- **Moved off Core** (no migration, no generator knowledge): the admin, remote-generation API, accounts and species pages, NAV, VIEW and the visitor-facing bugs go to Interface; test tooling, determinism helpers, benchmarks, scripts and the research items go to Quality; the planet-class chain (GEN.27, GEN.28, GEN.33, GEN.90 to GEN.92, GEN.145, GEN.146, GEN.29) is Interface's overflow after its own queue, with Core consulted before generation code is touched. Moves are listed in the lane-briefs UPDATE block.
- With these moves the whole open list finishes in about 214 slots against about 289 for the current split (26% sooner, almost all of it in the tail after the plan). The plan itself completes at about slot 84 either way, because it is Core-bound.
- **Idle time.** Core never idles. Interface has work to about slot 142. Quality has spare capacity (about 60 slots between slot 73 and its last blocked tests); Core may hand it determinism or tooling items on request.

**Handoffs** (the receiving lane keeps working its own queue until the item lands):

| From | To | What | About slot |
|---|---|---|---|
| Core GEN.260 | Interface UX.114 | members on scatter-row pages | 9 |
| Core GEN.235 | Interface API.27, UX.113 | profile lookup and sector page | 19 |
| Core GEN.250, GEN.251 | Quality TEST.160, TEST.161, DOC.49 | registry and fill-stage tests | 21, 34 |
| Core GEN.262, GEN.264 | Quality TEST.162, DOC.50 | multiplicity audit and docs | 46, 63 |
| Core GEN.236, GEN.241, GEN.243 | Quality TEST.159, DOC.48 | whole-sphere tests and docs | 40 to 66 |
| Core GEN.237, GEN.238 | Interface MAP.207 | satellites and wanderers on the map | 56, 72 |
| Core migration (GEN.223, GEN.230, GEN.260) | all lanes | Interface and Quality ask Core for any schema change; nothing else waits | 9 |

**Update (6:25 PM EDT): the catalogue and the region query are in.** Nine catalogue items (GEN.268 to GEN.275, TEST.164) and eight stage C items (GEN.276, DB.34, DB.35, NAV.60, NAV.61, TEST.165, TEST.166, DOC.51) were added to the graph (344 open items). With nothing moved, Core's plan queue grows from 29 items and 84 slots to 40 items and 113 slots, and the plan completes at slot 115. So four catalogue items that need no migration (GEN.270, GEN.271, GEN.274, GEN.275) and the two readers (NAV.60, NAV.61) go to the Interface lane, which has idle time in its queue, and Quality audits them (TEST.164) so it does not check its own work. Result: Core plan queue 36 items and 103 slots (was 84), the plan completes at slot 103, Interface plan work 10 items (23 slots), Quality 14 items (25 slots). GEN.276 goes right after GEN.232, then the migration items GEN.223, GEN.230, GEN.260, DB.34 and DB.35 (one migration, 0093); GEN.268 sits between GEN.250 and GEN.261. **Update (6:40 PM EDT):** GEN.223 gains the halo and intergalactic movement scales, a sector-crossing guard and full counting (Boss, 6:26 PM EDT), counted as 2 more slots, and PERF.103 (a per-night cost ceiling) follows GEN.265; Core's plan queue is 37 items and 106 slots and the plan completes at slot 106. Core order now begins: GEN.232, GEN.276, GEN.223, GEN.230, GEN.260, DB.34, DB.35, GEN.233, GEN.234, GEN.235, GEN.250, GEN.268, GEN.261, GEN.265, GEN.158, GEN.251, GEN.159, GEN.236, GEN.242, GEN.262, GEN.269, GEN.160, GEN.241, GEN.237, GEN.266, GEN.264, GEN.243, GEN.130, GEN.238, GEN.240, GEN.263, GEN.272, GEN.273, GEN.161, GEN.162, GEN.163.

**Two Core lanes (approved 7:00 PM EDT).** Boss asked for a second Core lane. Core and Core 2 negotiated the split from the computed one, and the TODO lane approved it; nothing already in progress moved. On my estimates the plan then completes at about slot 57 instead of 106, because the registry and member-tree chain (GEN.250, GEN.268, GEN.261) no longer waits behind the migration chain; the longest dependency chain, not lane length, is then the limit.

- **Core 1** (migration 0093 and the grid; owns every migration revision): GEN.232, then the one migration 0093 (the cell key part of GEN.276, GEN.223, GEN.230, GEN.233, GEN.260, DB.34, DB.35), GEN.234, GEN.265, PERF.103, GEN.266, GEN.236, GEN.237, GEN.264, GEN.238, GEN.263, GEN.272, GEN.161, GEN.163, then DB.30, POP.14, POP.13, DB.31, GEN.101, GEN.129.
- **Core 2** (registry, member tree, fill, nebulae): GEN.250, GEN.268, GEN.261 (pure functions first; it reads and writes the scatter member tables only after 0093 merges), GEN.158, GEN.159, GEN.262, GEN.269, GEN.160, GEN.235 (after Core 1's GEN.234), GEN.251, GEN.242, GEN.241, GEN.243, GEN.130, GEN.240, GEN.273, GEN.162, then GEN.229, DB.29, PERF.97 to PERF.100, GEN.179 to GEN.181, GEN.231 and the rest of the Core tail; plus DB.33 and the `sphere_runs()` part of GEN.276 (GEN.276 stays one item with two parts; NAV.60 and TEST.165 wait for both).
- **Rules**: Core 2 asks Core 1 before any item that adds a table or column and takes migration 0094 or later after 0093 merges. Handoffs: GEN.234 to GEN.235 (Core 1 to Core 2); GEN.223 and GEN.260 to GEN.264, GEN.265 and GEN.251 to GEN.264 (between the two); each tells the other when an item lands. A lane that wants to move an item tells the TODO lane first, and the TODO lane approves before it starts.
- The "Core" owner in the item texts and tables means whichever Core lane holds the item in these lists; the lists win.
