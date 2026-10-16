# Version later: Future features (not scheduled)

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: later

## Goal

Features that belong to no version yet, noted as possible future work: more than one galaxy in the database and the intergalactic frame (Boss, 2026-10-16: a possible future feature, not a headline).

## Phase 1: Multiple galaxies

Neighbour galaxies, more than one galaxy in the database and the intergalactic frame.

### Multiple galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.9 | Plan for more than one galaxy in the database |  | Plan only; floats. Judgment: if a second galaxy is likely, writing it before phase 2's schema work tells those migrations whether to add a galaxy id. |
| GEN.157 | The `neighbor_galaxies` table and a verified data file of about 25 real galaxies |  | Research: stage 1 of GEN.9. |
| GEN.164 | Synthetic globular-cluster systems for generated galaxies | GEN.9 | Globular-cluster chain. |
| MAP.145 | A sky and Galaxy Map drawing rule for neighbour galaxies | GEN.157 | Research. |
| NAV.57 | The Intergalactic Frame in navigation-frames.md and `navigation.py` | GEN.157 | Research: stage 2 of GEN.9. |
