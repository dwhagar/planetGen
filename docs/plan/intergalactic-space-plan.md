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

## 3. The stages, in build order

```
A  spherical grid ──► B1 halo + profile ──► B2 components + clusters ──► B4 scatters over the sphere
                          │                                                   │
                          ├──► B3 registry + fill machinery ──► B5 real contents (young stars, nebulae, W stars)
                          │
                          └──► B6 Interface (profile lookup, sector page, map)      B7 tests + docs beside each step
M  multi-star definitions (research running) plugs in beside B3 and B4        C  region query (waits on research)
L  later group: nothing waits on it
```

- **A** is first and alone on the critical path: the shared migration (GEN.230 and GEN.223 together, the 96-bit layout and the point-in-space columns) then the reach (GEN.233). GEN.232 is a bug that must land before any wide scatter.
- **B1 to B7** are the whole sphere. B1 is a strict chain (halo, then profile, then kind factors); after the profile exists, B2, B3 and B6 can run in parallel across lanes. B4 needs the components (GEN.237). B5 needs the registry and fill stages (GEN.251) and the pre-main-sequence state. B7 is written beside each step.
- **M** is a reserved slot (section 4). **C** is the region query (section 5). **L** is the later group (section 7).
- Core is the only lane that adds migrations and owns B1 to B5; Interface owns B6; Quality owns B7 and every TEST and DOC item.

## 4. Slot M: the multi-star definitions (research still running)

A research thread, "Research multi-star systems (2026-10-15)", is surveying how many stable multi-object systems exist, their theoretical limits and probabilities, and proposing expanded multi-star definitions (stars, former stars, protostars, neutron stars and black holes, in pairs and hierarchies up to the limit). Its items are not filed yet. They plug in here:

- **Definitions**: GEN.129 (multi-star systems of up to seven stars) is rewritten from the research; the new object kinds and the orbit trees live in the registry (GEN.250).
- **Scatter-time companions**: the companion roll moves out of the fill-time system pass (GEN.251) to the scatter stage, drawn from the row's own seed; each scatter row (GEN.240) carries its companions' masses and orbit, placed with the primary so a pair is one gravitational source. Planets, belts and moons stay at fill.
- **Compact pairs**: GEN.130 carries neutron star and black hole companions and the pairs between them; GEN.252 (brown dwarfs) and GEN.247 (accreting white dwarf binaries) read the same definitions.
- **Order**: B1 to B3 do not wait for M. B4 can start with single-object rows and add the companion field when M lands, as long as the registry reserves the field; B5 does not depend on M. When the research arrives, the TODO thread files its items here, amends GEN.129, GEN.130, GEN.240 and GEN.251, and places them in the stage tables below.

## 5. Stage C: the spherical-region database query

Waits on its research (docs/design/spherical-region-queries.md; an indexed Morton or Hilbert cell key of the current position, separate from the permanent object ID). It is planned only once the research note is finished, never blocks A or B, and is skipped if it is not ready when the spherical layout is. Its items are filed in the Core group when the research hands them over; they attach after B1 (they need the profile and the reach).

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
| Black hole companion fraction | 20% (rests on three Gaia systems) | GEN.251, GEN.130 |
| Intermediate-mass black hole count (about 1e5 in the Galaxy against a handful of candidates) | left alone | GEN.259 (a question for Boss) |
| Unstable wide binaries already generated | count first; the repair needs Boss | GEN.229 (outside this plan) |

## 7. Items by stage

Each row is generated from the item's own entry in TODO.md on 2026-10-15; the entry wins if it has changed since.

### Stage A: The spherical grid finishes first

| Item | What | Owner | Needs | Unlocks |
|---|---|---|---|---|
| GEN.232 | Cap the sector-path integration for fast bodies so one hypervelocity star cannot run for minutes (bug) | Core | none | GEN.238, GEN.240 |
| GEN.230 | Sector designations and object-ID birth fields cannot name bands at or beyond 2048 (bug) | Core | none | GEN.233 |
| GEN.223 | The nightly orbital update moves every object in the galaxy, and each object is stored as a point in space in its own columns (bug) | Core | none | GEN.233 |
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
| GEN.250 | A phenomena registry: one table for every kind of object, with a completeness test | Core | none | DOC.49, GEN.251, TEST.160 |
| GEN.251 | Sector fill stages on their own random streams, and a system pass so every object that can have a system rolls for one | Core | GEN.250, GEN.235 | DOC.49, GEN.130, GEN.241, GEN.252, TEST.160, TEST.161 |
| GEN.130 | Compact-remnant systems at sector fill: companions and pulsar planets for neutron stars and black holes | Core | GEN.251 | DOC.49, GEN.163, GEN.247, GEN.253 |
| GEN.163 | Type-B pulsar planets in globular clusters (GEN.130 follow-on) | Core | GEN.130, GEN.158, GEN.159 | (leaf) |

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
| API.27 | A read-only profile lookup for any address in the reach, so an uncharted sector can be previewed | Interface (Core supplies the profile function) | GEN.235 | UX.113 |
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
- Skipped on purpose (catalogue section 4.7): subspace anomalies (default off), wormholes and cosmic strings, dark matter nebulae, quark stars, Thorne-Zytkow objects, mergers, bursts, blanets.

## 9. Tests and docs

| Step | Tests (Quality) | Docs (Quality) |
|---|---|---|
| A | TEST.156, TEST.158 (the grid names, the math check at the reach) | DOC.47 |
| B1 to B6 | written beside each build step | none until the code lands |
| B3 | TEST.161 (stream independence, same seed same system, the probability table, no planet or moon on a scatter row) | DOC.49 |
| B5, B7 | TEST.159 (required stars per nebula kind, Wolf-Rayet counts, tiling at the reach, profile and empty-sector determinism, pre-placed share per population) | DOC.48 |
| B7 | TEST.160 (the representation audit: expected against observed counts per region class, the fraction with a system, required contents) | DOC.49 |

Gate for every step: a green `scripts/test.py changed`; `full --check` when a stage finishes; the PR names the run. Docs-only PRs merge themselves.

## 10. Keeping this plan true

The TODO tasks thread owns this file. It amends it whenever an item in it is filed, retired or re-ordered, whenever Boss decides one of the section 6 questions, and when the multi-star or region-query research hands over its items. Merged items are retired from TODO.md and this file with their PR.
