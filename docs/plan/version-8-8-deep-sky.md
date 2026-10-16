# Version 8.8: The deep sky

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 8.8

## Goal

Nebulae, clusters, satellites and the rarest stars drawn on the map and reachable by flight. Proposed by the TODO thread; not yet approved by Boss.

## Phase 1: Nebulae, clusters and rare stars

Nebula drawing, clusters and satellites on the map, and the later stellar phenomena.

### Nebula drawing

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.132 | Overlay markers for black holes, nebulae and habitable worlds |  | Later (Boss 2026-10-08 01:59Z colors). |
| MAP.142 | Nebulae have fuzzy, fading boundaries |  | GitHub issue [#713](https://github.com/dwhagar/planetGen/issues/713). |
| MAP.176 | Nebula regions on the Galaxy Map: grouped soft sprites in the dominant family colour | MAP.173 | Phase 2. Foundations lane 2. |
| MAP.177 | Color by "Nebula cover" on the Galaxy Map | MAP.174 | Phase 2. Foundations lane 2. |
| MAP.178 | Field nebula clouds drawn one by one from the nebula table | MAP.173 | Phase 2. Foundations lane 2. |
| MAP.179 | Scattered nebulae as markers from level 8 | GEN.202, MAP.173 | Phase 2. Foundations lane 2. |
| MAP.180 | Picking, hover and fly-to for nebula regions and cover | MAP.176 | Phase 2. Foundations lane 2. |
| MAP.181 | Dust colour for the dark nebula family on both themes |  | Phase 2. Foundations lane 2. |
| MAP.183 | Optional, default not built for the first release: near-view P and Q nebulae from the seeded field below about 250 pc |  | Phase 2. Optional, default not built. |
| MAP.194 | Nebulae and globular clusters open their own zoom window from the Galaxy Map | MAP.176 |  |

### Clusters and satellites on the map

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.209 | Globular clusters are parent objects for their stars and are named when their sector is filled |  |  |
| GEN.239 | High-velocity clouds and the Magellanic Stream as placed gas clouds (later) | GEN.237 | Whole-sphere build train, Boss 2026-10-15. Step 10. Core, later. |
| GEN.248 | Tidal streams as placed components of the halo (later) | GEN.237 | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.296 | The uncharted-star information panel shows the star's metallicity [Fe/H], which is already stored on bright_stars |  |  |
| GEN.297 | Milky Way globular cluster catalogue import: load Boss's Harris file into galaxy_components, replacing the generated stand-in clusters |  |  |
| MAP.207 | The Galaxy Map draws satellites, clusters and wanderers as points beyond the disk outline | GEN.237, GEN.238 | Whole-sphere build train, Boss 2026-10-15. Step 9. Interface. |

### Rare stars and phenomena

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.245 | Debris disks around A and F stars (later) | GEN.242 | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.246 | Luminous blue variables and hypergiants (later) | GEN.243 | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.247 | Accreting white dwarf binaries: novae, cataclysmic and symbiotic variables, Type Ia progenitors (later) | GEN.129, GEN.130 | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.249 | Be stars, carbon stars and Mira variables as star attributes (later) |  | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.253 | White-dwarf pollution flag and cool white dwarf presets | GEN.130 | Phenomena catalogue, Boss 2026-10-15. Step 6, later. Core. |
| GEN.255 | Drawn jets: Herbig-Haro flows, microquasar jets and pulsar jets along their spin axes | GEN.242 | Phenomena catalogue, Boss 2026-10-15. Step 7b, later. Core and Interface. |
| GEN.256 | A radiation hazard score derived from distance to hot stars, Wolf-Rayet stars, remnants and magnetars | GEN.243 | Phenomena catalogue, Boss 2026-10-15. Step 7c, later. Core. |
| GEN.257 | A Fermi-bubble overlay on the Galaxy Map |  | Phenomena catalogue, Boss 2026-10-15. Step 7d, later. Interface. |
