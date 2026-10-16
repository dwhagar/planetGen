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
| API.5 | API version and compatibility checking | API.4 |  |
| API.10 | Reservations: claimed sectors and id blocks per run |  | Reserved id blocks build on DB.3's id-block fix (PR #347). |
| API.11 | Staging tables | API.10 | Galaxy schema migration (staging); after NAV.10 in the writer queue. |
| API.12 | The download: seed, skeleton and name state | API.5 | Downloads the naming key, not name registries. Downloads the seed (what a seed means is GEN.39) and the name state (rules from GEN.46, done in PR #370). |

### Recipes

| ID | Item | Needs | Note |
|---|---|---|---|
| API.18 | Generate by recipe: JSON for sectors, systems, planets, moons and phenomena |  | Validation through the Pydantic models; 400 for nonsense, 422 for validation failures. |
| API.19 | Galaxy-scale recipes: build a whole galaxy, piece by piece, from JSON | API.18 | The far end; framework and capability, not a full Milky Way. |

## Phase 2: Remote generation

Generation without a database, uploads, verification and the incomplete uploads page.

### Remote generation

| ID | Item | Needs | Note |
|---|---|---|---|
| API.13 | Generation without a database | API.12 | Local generation with all workers. |
| API.14 | Upload routes, compressed, in batches | API.7, API.11 |  |
| API.8 | Verify uploaded data before it is finalized | API.11, API.14 |  |
| ADM.13 | Incomplete uploads page | API.10, API.8 |  |
| API.3 | Remote generate: generate on a local machine, upload through the API | API.10, API.11, API.12, API.13, API.14, API.8, ADM.13 | Parent; remote mode mirrors the CLI options GEN.51/52 settle. |
| API.17 | Remote generation reproduces what the server would make | API.12, API.13, GEN.57 |  |
