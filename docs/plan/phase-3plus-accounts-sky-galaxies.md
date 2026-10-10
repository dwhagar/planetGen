# Phase 3+: Accounts, the sky and more galaxies

Rebuilt on 2026-10-07 from Boss's lists of 2026-10-03 and 2026-10-07
and the dependency tree (every open item, its prerequisites and the
files it shares). [docs/TODO.md](../TODO.md) is the master file: it holds each
item's full text, and its "Plan: phases" section indexes every phase.
This file gives the phase's goal, its build threads with each item's
prerequisites in order, and its open questions. Where the two disagree,
TODO.md wins; when an item ships, it leaves TODO.md and its row here is
deleted in the same PR. Research notes, the files several items share and
the judgment calls behind the placement are in [notes.md](notes.md).

## Goal

The open-ended tail: user accounts (with API.6 keys, saved courses and Hill-radius emails), the view of the sky from a planet, the plan for more galaxies.

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### Accounts

| ID | Item | Needs | Note |
|---|---|---|---|
| USR.2 | Accounts with roles: user, admin and Owner |  | Control schema; no blockers, placed late by priority. |
| MAP.170 | The Galaxy Map gravity layer: the galaxy potential and region aggregates, coarse and optional | MAP.168, MAP.151 | Later; only if PERF.72 says it is cheap. |
| API.6 | User-level API keys, owned by the account that created them, that can read but not upload | USR.2 | Boss 01:31Z: every key belongs to the account that made it, so it follows USR.2. The API call logging he asked for is filed separately. |
| USR.3 | SMTP settings in the admin config | USR.2 |  |
| USR.4 | Invite-only sign-up by unique link | USR.2, USR.3 |  |
| USR.5 | Email loop for setting and resetting passwords | USR.3 |  |
| USR.6 | Owner transfer | USR.5 |  |
| USR.7 | A user-level interface with bookmarks | USR.2 | Bookmarks of any object use NAV.7 references. |
| USR.9 | `seo.privacy_note` text, "download my data" and "delete my account" on the account page |  | Research suggestion; not requested. |
| USR.8 | Every signed-in user can generate a one-off system | USR.2 | Boss 05:16Z. Today /admin/generate/system is admin-only and only admin accounts exist, so no earlier step is needed. Open question: per-user limit (default 30 an hour, admins unlimited). |
| USR.1 | User accounts | USR.2, USR.3, USR.4, USR.5, USR.6, USR.7 | Parent; closes with its subitems. |
| SEC.32 | Argon2id password hashing, a 1,024-character password limit and the `__Host-` cookie prefix |  | Research: before USR.2. |
| NAV.19 | Saved courses in the account (after USR.7) | USR.7, NAV.18 |  |
| GEN.111 | Email the admin when two objects are inside each other's Hill radius | GEN.109, USR.3 | Needs USR.3's SMTP settings. |

### View

| ID | Item | Needs | Note |
|---|---|---|---|
| VIEW.1 | View from a planet |  | Research session with Boss. |
| VIEW.4 | Constellation names in the name generator |  | Constellation names from the codec under the naming key, not sliced word lists. Floats: no blockers, can run any time. |
| VIEW.2 | A starmap seen from a planet. RESEARCH WITH BOSS FIRST | VIEW.1 |  |
| VIEW.10 | Declare `pillow` in `setup.py` if it is used, and add a HYG calibration test |  | Research. |
| VIEW.9 | Replace the BC_V table with a published relation (Flower 1996 or Torres 2010) |  | Research. |
| VIEW.8 | A dust model for the sky: the galactic extinction law and `nebulae.extinction_av` | VIEW.1 | Research. |
| VIEW.7 | Sky tiers and the backfill cost measurement | VIEW.1 | Research: the sky is about 8% to 10% complete today. |
| VIEW.6 | Settle the handedness of the generated galaxy before mapping real sky coordinates |  | Research: needed before real sky coordinates. |
| VIEW.3 | Render the view as a PNG, with constellations | VIEW.2, VIEW.4 |  |

### Galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.9 | Plan for more than one galaxy in the database |  | Plan only; floats. Judgment: if a second galaxy is likely, writing it before phase 2's schema work tells those migrations whether to add a galaxy id. |
| GEN.164 | Synthetic globular-cluster systems for generated galaxies | GEN.9, GEN.158, GEN.159 | Globular-cluster chain. |
| MAP.145 | A sky and Galaxy Map drawing rule for neighbour galaxies | GEN.157 | Research. |
| NAV.57 | The Intergalactic Frame in navigation-frames.md and `navigation.py` | GEN.157 | Research: stage 2 of GEN.9. |
| GEN.157 | The `neighbor_galaxies` table and a verified data file of about 25 real galaxies |  | Research: stage 1 of GEN.9. |

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.55 | Same seed, same data: a sector's contents depend only on the seed, the version and its address (internal) |  | Parent of the chain. |

## Open questions for Boss

- USR.2: Accounts with roles: user, admin and Owner, see its entry in TODO.md.
- USR.3: SMTP settings in the admin config, see its entry in TODO.md.
- USR.4: Invite-only sign-up by unique link, see its entry in TODO.md.
- USR.5: Email loop for setting and resetting passwords, see its entry in TODO.md.
- USR.6: Owner transfer, see its entry in TODO.md.
- USR.7: A user-level interface with bookmarks, see its entry in TODO.md.
- USR.8: Every signed-in user can generate a one-off system, see its entry in TODO.md.
- VIEW.4: Constellation names in the name generator, see its entry in TODO.md.
- VIEW.3: Render the view as a PNG, with constellations, see its entry in TODO.md.
- GEN.9: Plan for more than one galaxy in the database, see its entry in TODO.md.
