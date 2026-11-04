# Version later: Future features (not scheduled)

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: later

Status: paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

Features that belong to no version yet, noted as possible future work: more than one galaxy in the database and the intergalactic frame (Boss, 2026-10-16: a possible future feature, not a headline).

## Phase 1: Multiple galaxies

Neighbour galaxies, more than one galaxy in the database and the intergalactic frame.

### Multiple galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.157 | The `neighbor_galaxies` table and a verified data file of about 25 real galaxies |  | Research: stage 1 of GEN.9. |
| GEN.164 | Synthetic globular-cluster systems for generated galaxies | GEN.387 | Globular-cluster chain. |
| MAP.145 | A sky and Galaxy Map drawing rule for neighbour galaxies | GEN.157 | Research. |
| NAV.57 | The Intergalactic Frame in navigation-frames.md and `navigation.py` | GEN.157 | Research: stage 2 of GEN.9. |

## Phase 2: Waiting on other projects

Small clean-ups that wait for an upstream release; none is part of a version's build.

### Optional unit edge arithmetic

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.588 | Optional: operator arithmetic on Quantity for the edges (mass divided by volume gives a density), unscheduled |  | Optional; Boss chose R1 option 1, so nothing waits on it. |

### Paging audit

The leftover offset-paged queries, checked one by one; nothing in 8.1 waits on it.

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.219 | The Systems list with a star type filter and the sector or octant sort pages with OFFSET only (a G filter matches about a tenth of the galaxy); move it to keyset paging, with a MySQL measurement at a full galaxy against the 50 ms target |  | PERF.137 audit follow-up |

### Upstream follow-ups

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.222 | Remove the NumPy 2.5 deprecation filter around the skimage marching-cubes call once scikit-image stops raising it |  | Waits on a scikit-image release. Core 5 found it. |

## Phase 3: Speed-ups that change generated values

Boss, 12:27 PM EDT 2026-10-23 (cards): not in 8.1; each needs a generator epoch and Boss's word before it starts.

### Value-changing speed-ups

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.162 | Member-tree drawing is 53 percent of a band unit, and speeding it changes values, so it needs a generator epoch and Boss's word |  | Boss: not in 8.1 (12:27 PM EDT 2026-10-23). |
| PERF.166 | The mean-system-mass lookup bisects per age (0.2 ms each); a faster form changes the last bits and needs an epoch and Boss's word (bug) |  | Boss: leave it (12:27 PM EDT 2026-10-23); build only if a measurement after PERF.165 says it is worth an epoch. |
