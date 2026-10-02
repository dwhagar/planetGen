# Phase 4: Accounts, the view from a planet, more galaxies

This is one phase of the plan drawn up on 2026-10-01 from Boss's list
of that evening and his research notes. [docs/TODO.md](../TODO.md) is
the master file: it holds each item's full text (what's wrong, where to
look, what "done" means), and its "Plan: phases" section indexes every
phase. This file holds what the phase needs beyond that: its goal, the
order and dependencies, how to split it into build threads, and the
research notes that apply. Where the two disagree, TODO.md wins; when an
item ships, it leaves TODO.md and its row here is deleted in the same PR.

## Goal

User accounts with roles and bookmarks (and saved courses in the
account), the view of the sky from a planet after its research session
with Boss, and the plan for more than one galaxy.

## Items

### User accounts

| ID | Item | Parent |
|---|---|---|
| USR.1 | User accounts |  |
| USR.2 | Accounts with roles: user, admin and Owner | USR.1 |
| USR.3 | SMTP settings in the admin config | USR.1 |
| USR.4 | Invite-only sign-up by unique link | USR.1 |
| USR.5 | Email loop for setting and resetting passwords | USR.1 |
| USR.6 | Owner transfer | USR.1 |
| USR.7 | A user-level interface with bookmarks | USR.1 |
| NAV.19 | Saved courses in the account (after USR.7) | NAV.4 |

USR.2 (roles) first; USR.3 (SMTP) before USR.4 (invites) and USR.5
(password emails), then USR.6 (Owner transfer). USR.7 (the user
interface and bookmarks) after USR.2, and NAV.19 (saved courses in the
account) after USR.7.

### The view from a planet

| ID | Item | Parent |
|---|---|---|
| VIEW.1 | View from a planet |  |
| VIEW.2 | A starmap seen from a planet. RESEARCH WITH BOSS FIRST | VIEW.1 |
| VIEW.3 | Render the view as a PNG, with constellations | VIEW.1 |
| VIEW.4 | Constellation names in the name generator | VIEW.1 |

VIEW.2 and VIEW.3 wait for VIEW.1's research session with Boss. VIEW.4
(constellation names) can be built any time.

### More than one galaxy

| ID | Item | Parent |
|---|---|---|
| GEN.9 | Plan for more than one galaxy in the database |  |

A planning document only.

### Reproducible galaxies: the end state

| ID | Item | Parent |
|---|---|---|
| GEN.55 (new) | A version number and a seed reproduce the same galaxy (end goal) |  |
| OPS.12 (new) | `generate.py reproduce`: a version and a seed rebuild a galaxy and check it | GEN.55 |

OPS.12 needs everything before it in GEN.55's chain (phases 1 to 3)
and closes GEN.55: `generate.py reproduce --seed X --version Y`
rebuilds a galaxy or region into a fresh database and checks its
fingerprint.

## Research notes

None of Boss's 2026-10-01 notes apply to this phase.

## Build threads

Each thread is briefed with its exact item IDs and takes no others.

1. Accounts: USR.2, USR.3, USR.4, USR.5, USR.6, USR.7, NAV.19.
2. VIEW.4 any time; VIEW.1 research session, then VIEW.2 and VIEW.3.
3. GEN.9 plan.
4. Reproducible galaxies: OPS.12 (closes GEN.55).

## Open questions for Boss

- VIEW.1: the research session with Boss.
- USR.1 and GEN.9: the open questions in their TODO.md entries.
