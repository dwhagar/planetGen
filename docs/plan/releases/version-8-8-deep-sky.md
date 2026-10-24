# Version 8.8: The deep sky

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 8.8

Status: paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

Nebulae, clusters, satellites and the rarest stars drawn on the map and reachable by flight. Proposed by the TODO thread; not yet approved by Boss.

## Phase 1: Nebulae, clusters and rare stars

Nebula drawing, clusters and satellites on the map, and the later stellar phenomena.

### Nebula drawing

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.176 | Nebula regions on the Galaxy Map: grouped soft sprites in the dominant family colour | MAP.173 | Phase 2. Foundations lane 2. |
| MAP.179 | Scattered nebulae as markers from level 8 | GEN.202, MAP.173 | Phase 2. Foundations lane 2. |
| MAP.142 | Nebulae have fuzzy, fading boundaries |  | GitHub issue [#713](https://github.com/dwhagar/planetGen/issues/713). |
| MAP.177 | Color by "Nebula cover" on the Galaxy Map |  | Phase 2. Foundations lane 2. |
| MAP.178 | Field nebula clouds drawn one by one from the nebula table | MAP.173 | Phase 2. Foundations lane 2. |
| MAP.180 | Picking, hover and fly-to for nebula regions and cover | MAP.176 | Phase 2. Foundations lane 2. |
| MAP.181 | Dust colour for the dark nebula family on both themes |  | Phase 2. Foundations lane 2. |
| MAP.183 | Optional, default not built for the first release: near-view P and Q nebulae from the seeded field below about 250 pc |  | Phase 2. Optional, default not built. |
| MAP.222 | Nebula zoom window: choosing a nebula on the Galaxy Map flies into a view of that nebula alone | MAP.176 | Split from MAP.194. |
| MAP.215 | Nebulae tint the stars seen through them, with soft translucent boundaries | MAP.176, MAP.142 | Split from MAP.191. |
| MAP.191 | Galaxy Map default view: no density shading, only nebulae shade the stars behind them | MAP.214, MAP.215 |  |
| MAP.233 | Names for nebulae on the Galaxy Map | MAP.229, MAP.176 | Map labels, Boss 2:30 PM EDT. |

### Clusters and satellites on the map

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.359 | Sector fill links cluster stars to their cluster and names the cluster at its first fill |  | Split from GEN.209. |
| MAP.224 | Cluster pages and the Galaxy Map use the parent link | GEN.359 | Split from GEN.209. |
| GEN.209 | Globular clusters are parent objects for their stars and are named when their sector is filled | GEN.359, MAP.224 |  |
| GEN.239 | High-velocity clouds and the Magellanic Stream as placed gas clouds (later) |  | Whole-sphere build train, Boss 2026-10-15. Step 10. Core, later. |
| GEN.248 | Tidal streams as placed components of the halo (later) |  | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.296 | The uncharted-star information panel shows the star's metallicity [Fe/H], which is already stored on bright_stars |  |  |
| MAP.207 | The Galaxy Map draws satellites, clusters and wanderers as points beyond the disk outline |  | Whole-sphere build train, Boss 2026-10-15. Step 9. Interface. |
| MAP.232 | Names for the big features: spiral arms, bar, bulge, halo, satellites and clusters | MAP.229, MAP.207 | Map labels, Boss 2:30 PM EDT. |
| GEN.367 | The Magellanic Clouds' gas is not in `kind_factors` and `nebula_field`: the Clouds hold stars but no nebulae |  | Finding from Core 3. |

### Rare stars and phenomena

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.245 | Debris disks around A and F stars (later) |  | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.246 | Luminous blue variables and hypergiants (later) |  | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.247 | Accreting white dwarf binaries: novae, cataclysmic and symbiotic variables, Type Ia progenitors (later) | GEN.272 | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.249 | Be stars, carbon stars and Mira variables as star attributes (later) |  | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.255 | Drawn jets: Herbig-Haro flows, microquasar jets and pulsar jets along their spin axes |  | Phenomena catalogue, Boss 2026-10-15. Step 7b, later. Core and Interface. |
| GEN.256 | A radiation hazard score derived from distance to hot stars, Wolf-Rayet stars, remnants and magnetars |  | Phenomena catalogue, Boss 2026-10-15. Step 7c, later. Core. |
| GEN.257 | A Fermi-bubble overlay on the Galaxy Map |  | Phenomena catalogue, Boss 2026-10-15. Step 7d, later. Interface. |

## Phase 2: Windows into nebulae and clusters

The zoom windows: a nebula alone, and the stars of a globular cluster. They come after the nebula drawing and the cluster parent link of phase 1.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### Zoom windows

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.223 | Globular cluster zoom window: choosing a cluster flies into a view of its own stars | MAP.222, MAP.224 | Split from MAP.194. |
| MAP.194 | Nebulae and globular clusters open their own zoom window from the Galaxy Map | MAP.222, MAP.223 |  |
