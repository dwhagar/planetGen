# Version 9.5: Life classes and sentient life

Filed on 2026-10-28 as a version map (Boss, 7:41 PM EDT: the sentient life portion is left for a 9.5 update, where the class system is expanded with a life class system).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 9.5

Status: paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then).

## Goal

The life class system for carbon and water life moved into 8.1 (Boss, 9:10 PM EDT 2026-10-28: GEN.862 and GEN.863). This version adds sentient life drawn from those classes, the census for it, and nothing else. The 8.x work stops at animals (Boss, 7:41 PM EDT 2026-10-28: near-Earth, water and carbon based worlds); any chemistry that is not carbon based needs Boss's approval before it is built. The share of worlds with sentient life stays as it is (Boss, 7:38 PM EDT: no more intelligent worlds than now, but a mix of kinds).

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

## Phase 1: Sentient life and the census

### Sentient life

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.864 | Sentient life per life class: sentient species and civilizations draw from the life class scores of a world | GEN.863 | Share of sentient worlds unchanged. |

### Census

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.435 | Prevalence tests for the life classes: the share of worlds in each life class and stage inside researched ranges | GEN.863, TEST.430 | Extends the TEST.430 suite. |

## Phase 2: The headline closes

### The headline closes

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.865 | Headline of 9.5: sentient life draws from the life classes, so every sentient world is rated for the kind of life that built it | GEN.864, TEST.435 | Headline of 9.5. |
