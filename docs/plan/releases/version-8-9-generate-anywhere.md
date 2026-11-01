# Version 8.9: Generate anywhere

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 8.9

Status: paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

Remote generation: build a galaxy on a local machine, upload it through the API and have the server verify it, plus generate by recipe. Proposed by the TODO thread; not yet approved by Boss. Boss (2026-10-16, 11:26 PM EDT: "Everything for 9.0 before splits between 9.1 and 8.9") also put here what was planned for 9.0 and needs nothing of the atmosphere colour: the sky data and research of the view from a planet (phase 3) and the backups and repair chain (ADM.91, phases 4 and 5).

## Phase 1: Upload groundwork and recipes

The reservation, staging and download pieces and the recipe format.

### Upload groundwork

| ID | Item | Needs | Note |
|---|---|---|---|
| API.10 | Reservations: claimed sectors and id blocks per run |  | Reserved id blocks build on DB.3's id-block fix (PR #347). |
| API.29 | Upload spool and unit states: received batches on disk, units assembled at seal, finalize as a queue job | API.10 | Split from API.11. |
| API.11 | Staging tables | API.29 | Galaxy schema migration (staging); after NAV.10 in the writer queue. |
| API.5 | API version and compatibility checking |  |  |
| API.12 | The download: seed, skeleton and name state | API.5 | Downloads the naming key, not name registries. Downloads the seed (what a seed means is GEN.39) and the name state (rules from GEN.46, done in PR #370). |

### Recipes

| ID | Item | Needs | Note |
|---|---|---|---|
| API.18 | Generate by recipe: JSON for sectors, systems, planets, moons and phenomena | API.28 | Narrowed: the schema and runner are API.28, built in 8.6; this is the route, the job and the 400 and 422 errors. |
| API.19 | Galaxy-scale recipes: build a whole galaxy, piece by piece, from JSON | API.18, GEN.302 | The far end; framework and capability, not a full Milky Way. |

## Phase 2: Remote generation

Generation without a database, uploads, verification and the incomplete uploads page.

### Remote generation

| ID | Item | Needs | Note |
|---|---|---|---|
| API.13 | Generation without a database | API.12 | Local generation with all workers. |
| API.14 | Upload routes, compressed, in batches | API.11 |  |
| API.8 | Verify uploaded data before it is finalized | API.11, API.14, ADM.82 | Narrowed: the checks and the finalize are ADM.82, built in 8.6 |
| ADM.13 | Incomplete uploads page | API.10, API.8 |  |
| API.30 | generate.py --remote: the client driver for sector and galaxy runs | API.5, API.12, API.13, API.14, API.8 | Split from API.3. |
| ADM.95 | A remote run is a job: it shows in the job tree and the job page with its progress and report | API.10, API.14, API.8 | Split from API.3. |
| TEST.552 | test_plan_outline expects the outline crest at 3.32 times and the density model now gives 3.18 (bug) |  |  |
| ADM.136 | A New galaxy job stalls waiting on a generation run started by another process, and the holder is not identified (bug) |  |  |
| API.3 | Remote generate: generate on a local machine, upload through the API | API.30, ADM.95 | Parent; remote mode mirrors the CLI options GEN.51/52 settle. |
| API.17 | Remote generation reproduces what the server would make | API.12, API.13, API.8 |  |

### Resources in uploads

| ID | Item | Needs | Note |
|---|---|---|---|
| API.34 | API and remote uploads carry resource locations: the resource fields are in the API, and an uploaded galaxy is checked for them | GEN.799 | Boss, 10:29 PM EDT 2026-10-24: resource locations (valuable belts and fields as points of interest), foundations 8.2 to 9.0. |


## Phase 3: Sky data and research

The sky's data table and the research that settles the view from a planet, moved here from 9.0 (Boss, 10:27 PM EDT: "Everything for 9.0 before splits between 9.1 and 8.9"): nothing in it needs the atmosphere colour, so it is ready before the view is built in 9.1.

### Sky data

| ID | Item | Needs | Note |
|---|---|---|---|
| VIEW.4 | Constellation names in the name generator |  | Constellation names from the codec under the naming key, not sliced word lists. Floats: no blockers, can run any time. |
| VIEW.6 | Settle the handedness of the generated galaxy before mapping real sky coordinates |  | Research: needed before real sky coordinates. |

### Sky research

| ID | Item | Needs | Note |
|---|---|---|---|
| VIEW.1 | View from a planet |  | Research session with Boss. |
| VIEW.7 | Sky tiers and the backfill cost measurement | VIEW.1 | Research: the sky is about 8% to 10% complete today. |
| VIEW.8 | A dust model for the sky: the galactic extinction law and `nebulae.extinction_av` | VIEW.1 | Research. |

## Phase 4: Backup parts and restore

Target 8.9 (moved from 9.0; Boss, 2:00 PM EDT): the control database and parity parts of the backup archive, and a restore for each.

### Control database

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.88 | Control database backup part: a dump of the planetgen_control schema (admin logins, API keys, run records) in the backup archive |  | Boss 1:59 PM EDT; part of the ADM.87 archive. |
| ADM.89 | Control database restore: load a backed-up planetgen_control into a scratch schema, check it, then replace the live schema, from the command and from the settings page | ADM.88, OPS.41 | Boss 2:00 PM EDT. |

### Parity data

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.40 | Parity data backup part: the parity file of the DB repair copied into the backup archive |  | Boss 1:59 PM EDT; Core 1 owns DB.9. |
| DB.41 | Parity data restore: put a backed-up parity file back and check it against the database | DB.40, DB.8 | Boss 2:00 PM EDT; Core 1. |

## Phase 5: Repair from backups, the repair page and proof

The repair of the galaxy from an archive, the one page for check, repair and restore, and the end-to-end test, guide and umbrella (moved from 9.0).

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### Repair and proof

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.42 | Repair the galaxy from an archive: `repair-db --from-archive` uses the parity part of a backup when the live parity file is missing or damaged | DB.41, DB.17 | Boss 2:00 PM EDT. |
| ADM.90 | Database health and repair page: run the check, see damaged and unprotected sectors, repair from the live parity or from an archive, restore parts, as a queue job | DB.42, ADM.89, DB.41 | Boss 2:00 PM EDT. |
| OPS.55 | Nightly full backup: a backup step in the daily maintenance run that takes every part into one archive, with a schedule setting and a nightly retention setting | OPS.16, ADM.88, DB.40 | Boss 2:06 PM EDT. |
| TEST.180 | Backup and repair end to end: take a full archive, lose the control schema and the live parity file, restore and repair, and prove the galaxy and the logins come back | ADM.89, DB.42, ADM.90, OPS.55 | Boss 2:00 PM EDT. |
| DOC.56 | Backups, restore and repair: the deployment guide and the help page | TEST.180 | Boss 2:00 PM EDT. |
| ADM.91 | Backups and repair, target 8.9: one bz2 archive of the settings, the control database and the parity data, and the ways to use it (restore the control database, restore parity, repair the galaxy, the repair page) | ADM.88, ADM.89, DB.40, DB.41, DB.42, ADM.90, TEST.180, DOC.56, OPS.55 | Umbrella of the backups and repair chain; done when its items are. |
