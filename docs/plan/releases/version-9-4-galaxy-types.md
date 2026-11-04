# Version 9.4: Galaxy types

Filed on 2026-10-16 as a version map (Boss, 8:48 PM EDT: for 9.4, different galactic models for generation based on known galaxy types, as a chain of required items counted back from 9.4).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 9.4

Status: paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

A galaxy is built from one family of a known galaxy type (spiral, lenticular, elliptical, irregular or dwarf): the seed draws the galaxy's stellar mass, then the family from the mass-dependent type fractions, then every shape number from published scaling relations, so no type, size or count is a preset and no named galaxy is copied (Boss, 4:11 and 4:14 PM EDT: no fixed content). The shape, the density model, the potential, the halo, the satellites and the scatters are reused: three steps of the chain come earlier (GEN.378 in 8.1, DB.66 and GEN.379 in 8.3) so that the families only fill in the shape. Design: [docs/design/galaxy/galaxy-types.md](../../design/galaxy/galaxy-types.md). One galaxy per database stays; a second galaxy is the later note.

## Phase 1: The draw and the spheroid

The type draw and the Sérsic spheroid component, which every family reads.

### Type draw

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.380 | The galaxy type draw: stellar mass from a Schechter function, then the family from the mass-dependent fractions, then the halo mass |  | Step 4 of the galaxy types chain. |

### Spheroid component

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.381 | A Sérsic spheroid as a density component: index, effective radius, axis ratios and orientation, with a mass |  | Step 4 of the galaxy types chain. |

### Family outlines

| ID | Item | Needs | Note |
|---|---|---|---|

## Phase 2: The families

Three independent families.

### Spiral family

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.382 | The spiral family: scale length and height, bulge-to-total, bar, arm count and pitch angle from mass and Hubble stage | GEN.380 | Step 5 of the galaxy types chain; generalizes the shape row of GEN.368. |

### Elliptical family

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.383 | The elliptical family: Sérsic index and size from mass, intrinsic axis ratios, no disk, arms or bar, an old population | GEN.380, GEN.381 | Step 5 of the galaxy types chain. |

### Irregular and dwarf families

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.385 | The irregular and dwarf families: a clumpy exponential or Sérsic body with young stars, and an old dwarf spheroid | GEN.380, GEN.381 | Step 5 of the galaxy types chain; reuses the satellite laws of GEN.369. |
| GEN.619 | A variable sector edge for small galaxies: 4, 2 or 1 pc chosen by the thin disc scale length, with the one-star outline rule as a density threshold |  | Boss, 8:08 PM EDT 2026-10-21 (later item; replaces GEN.616); version is an estimate. |

## Phase 3: Lenticular, populations, clusters and the nucleus

What reads the families.

### Lenticular family

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.384 | The lenticular family: a smooth disk and a large bulge, no arms, old stars and almost no gas | GEN.382, GEN.383 | Step 6 of the galaxy types chain. |

### Populations

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.386 | Population mix, star-formation history and gas fraction by family: the young share that the nebula and profile rules read | GEN.382, GEN.383, GEN.385, GEN.405 | Step 6 of the galaxy types chain. |

### Clusters and the nucleus

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.387 | Globular clusters, the nuclear cluster and the central black hole follow the family's spheroid and halo | GEN.382, GEN.383, GEN.297 | Step 6 of the galaxy types chain. |
| GEN.635 | The cluster system has no radial profile yet, and the dwarf end of the globular cluster count uses no occupation draw (the Carlsten fits are not built in GEN.403) |  |  |
| GEN.634 | The plan outline takes a whole shell in by its highest band, so the candidate sector count is 3.5 billion where about 1.0 billion expect a system (interval outline) |  |  |


## Phase 4: Creation

The seed draws the family and the galaxy stores it.

### Galaxy creation

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.389 | Galaxy creation draws the family from the seed and stores it, with the family in the log, the API and the version check | GEN.384, GEN.386, GEN.387, GEN.385, GEN.382, GEN.383 | Step 7 of the galaxy types chain. |

## Phase 5: Pages, tests and docs

### Generate page

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.134 | Generate page: the galaxy type shows what the seed drew and can be pinned | GEN.389 | Step 8 of the galaxy types chain. |

### Galaxy Map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.242 | The Galaxy Map and its key draw what the family has: arms and bar only where they exist, a spheroid shell, a lens or a clumpy body | GEN.389 | Step 8 of the galaxy types chain. |

### Plausibility tests

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.204 | Family plausibility tests: masses, sizes, indices, bulge shares, bar fractions and arm pitch fall inside the published ranges | GEN.389 | Step 8 of the galaxy types chain. |

### Docs

| ID | Item | Needs | Note |
|---|---|---|---|
| DOC.59 | Docs: the galaxy types note moves to built and the density, potential, schema, API and Generate page notes say the shape holds a family and components | GEN.389 | Step 8 of the galaxy types chain. |

## Phase 6: The headline closes

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### The headline closes

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.390 | Headline of 9.4: generate a galaxy of any known type, with spiral, lenticular, elliptical, irregular and dwarf families from one set of rules | UX.134, MAP.242, TEST.204, DOC.59 | Headline of 9.4. |
