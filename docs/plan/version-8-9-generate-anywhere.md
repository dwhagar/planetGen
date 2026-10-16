# Version 8.9: Generate anywhere

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 8.9

## Goal

Remote generation: build a galaxy on a local machine, upload it through the API and have the server verify it, plus generate by recipe. Proposed by the TODO thread; not yet approved by Boss.

## Phase 1: Upload groundwork and recipes

The reservation, staging and download pieces and the recipe format.

### Upload groundwork

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.62 | Upload bookkeeping tables in the control database (migration only) |  | Split from API.11. Core 1 owns migrations. |
| API.10 | Reservations: claimed sectors and id blocks per run | DB.62 | Reserved id blocks build on DB.3's id-block fix (PR #347). |
| API.29 | Upload spool and unit states: received batches on disk, units assembled at seal, finalize as a queue job | DB.62, API.10 | Split from API.11. |
| API.11 | Staging tables | DB.62, API.29 | Galaxy schema migration (staging); after NAV.10 in the writer queue. |
| API.5 | API version and compatibility checking | API.4 |  |
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
| API.3 | Remote generate: generate on a local machine, upload through the API | API.30, ADM.95 | Parent; remote mode mirrors the CLI options GEN.51/52 settle. |
| API.17 | Remote generation reproduces what the server would make | API.12, API.13, GEN.57, API.8 |  |
