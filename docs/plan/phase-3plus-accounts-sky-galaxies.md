# Phase 3+: Accounts, the sky and more galaxies

Rebuilt on 2026-10-02 from the dependency report (every open item, its
prerequisites and the files it shares), with Boss's decisions of that
night. [docs/TODO.md](../TODO.md) is the master file: it holds each
item's full text, and its "Plan: phases" section indexes every phase.
This file gives the phase's goal, its build threads with each item's
prerequisites in order, and its open questions. Where the two disagree,
TODO.md wins; when an item ships, it leaves TODO.md and its row here is
deleted in the same PR. Research notes, the files several items share and
the judgment calls behind the placement are in [notes.md](notes.md).

## Goal

The open-ended tail: user accounts (with API.6 keys and saved courses), the view of the sky from a planet, the plan for more galaxies, and the end state of reproducible galaxies (`generate.py reproduce`).

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### Accounts

| ID | Item | Needs | Note |
|---|---|---|---|
| API.6 | User-level API keys, owned by the account that created them, that can read but not upload | API.9, USR.2 | Boss 01:31Z: every key belongs to the account that made it, so it follows USR.2. The API call logging he asked for is filed separately. |
| USR.2 | Accounts with roles: user, admin and Owner |  | Control schema; no blockers, placed late by priority. |
| USR.3 | SMTP settings in the admin config | USR.2 |  |
| USR.4 | Invite-only sign-up by unique link | USR.2, USR.3 |  |
| USR.5 | Email loop for setting and resetting passwords | USR.3 |  |
| USR.6 | Owner transfer | USR.5 |  |
| USR.7 | A user-level interface with bookmarks | USR.2, NAV.7 | Bookmarks of any object use NAV.7 references. |
| USR.1 | User accounts | USR.2, USR.3, USR.4, USR.5, USR.6, USR.7 | Parent; closes with its subitems. |
| NAV.19 | Saved courses in the account (after USR.7) | USR.7, NAV.18 |  |

### View

| ID | Item | Needs | Note |
|---|---|---|---|
| VIEW.1 | View from a planet |  | Research session with Boss. |
| VIEW.4 | Constellation names in the name generator |  | Floats: no blockers, can run any time. |
| VIEW.2 | A starmap seen from a planet. RESEARCH WITH BOSS FIRST | VIEW.1 |  |
| VIEW.3 | Render the view as a PNG, with constellations | VIEW.2, VIEW.4 |  |

### Galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.9 | Plan for more than one galaxy in the database |  | Plan only; floats. Judgment: if a second galaxy is likely, writing it before phase 2's schema work tells those migrations whether to add a galaxy id. |

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.55 | A version number and a seed reproduce the same galaxy (end goal) | OPS.12 | Parent of the chain. |
| OPS.12 | `generate.py reproduce`: a version and a seed rebuild a galaxy and check it | DB.7, GEN.57, GEN.58, TEST.77, GEN.59, OPS.14, ADM.18, GEN.61, OPS.18 | The end state. |

## Open questions for Boss

- USR.2: Accounts with roles: user, admin and Owner, see its entry in TODO.md.
- USR.3: SMTP settings in the admin config, see its entry in TODO.md.
- USR.4: Invite-only sign-up by unique link, see its entry in TODO.md.
- USR.5: Email loop for setting and resetting passwords, see its entry in TODO.md.
- USR.6: Owner transfer, see its entry in TODO.md.
- USR.7: A user-level interface with bookmarks, see its entry in TODO.md.
- VIEW.4: Constellation names in the name generator, see its entry in TODO.md.
- VIEW.3: Render the view as a PNG, with constellations, see its entry in TODO.md.
- GEN.9: Plan for more than one galaxy in the database, see its entry in TODO.md.
