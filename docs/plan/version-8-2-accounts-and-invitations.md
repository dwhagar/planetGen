# Version 8.2: User accounts with email invitations

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 8.2

## Goal

User accounts, with sign-up by an email invitation that only an admin can send (USR.1): roles, SMTP settings, invitations, password emails and owner transfer first, then what a signed-in user gets: bookmarks in the account, API keys, the one-off system page and a home field.

## Phase 1: Accounts, invitations and bookmarks

Roles, the Owner tier of the settings page, SMTP, invitations, password emails, owner transfer and the bookmark manager.

### Accounts and sign-up

| ID | Item | Needs | Note |
|---|---|---|---|
| SEC.32 | Argon2id password hashing, a 1,024-character password limit and the `__Host-` cookie prefix |  | Research: before USR.2. |
| DB.58 | Control-schema revision for accounts: role, email and status, the single-Owner index, account_tokens and invite_uses |  | Split from USR.2. Core 1 owns migrations. |
| USR.26 | Roles in sessions, API keys and route guards; the Owner rules; the accounts CLI | DB.58, SEC.32 | Split from USR.2. |
| ADM.86 | Owner-tier settings: Owner-only options behind a recent-password prompt | USR.26 | Boss 1:43 PM EDT: the settings page moved to 8.1; its Owner tier needs the Owner account of USR.2. |
| USR.3 | SMTP settings in the admin config | ADM.86, USR.26 |  |
| USR.28 | Account tokens: one issue-and-redeem library for invite, first-password, reset, email-change and owner-transfer links | DB.58 | Foundation found in the review. |
| USR.5 | Email loop for setting and resetting passwords | USR.3, USR.28 |  |
| USR.4 | Invite-only sign-up by unique link | DB.58, USR.26, USR.3, USR.28, USR.5 |  |
| USR.6 | Owner transfer | USR.5, USR.28 |  |

### Bookmarks

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.129 | Saved lists in the browser and in the account: one store with a local and an account backend |  | Foundation found in the review. |
| UX.46 | Wiping the galaxy wipes the bookmarks | UX.129 |  |
| UX.47 | A bookmark manager: list, rename, sort, group and delete | UX.46, UX.129 |  |
| UX.48 | Charted regions as bookmarks that frame and outline the region | UX.47, MAP.190 |  |
| UX.45 | Bookmark management | UX.46, UX.47, UX.48 | Parent; per browser until accounts (USR.7) move them. |
| UX.92 | A Bookmark button for planets, moons and belts | UX.45 |  |

## Phase 2: What a signed-in user gets

Account bookmarks, user API keys, the one-off system page, the home field and the umbrella.

### User features

| ID | Item | Needs | Note |
|---|---|---|---|
| API.6 | User-level API keys, owned by the account that created them, that can read but not upload | DB.58, USR.26 | Boss 01:31Z: every key belongs to the account that made it, so it follows USR.2. The API call logging he asked for is filed separately. |
| USR.7 | A user-level interface with bookmarks | DB.58, USR.26, UX.129 | Bookmarks of any object use NAV.7 references. |
| USR.27 | Admin Users page: list accounts, promote, demote, disable and delete | USR.26 | Split from USR.2. |
| USR.2 | Accounts with roles: user, moderator, admin and Owner | DB.58, USR.26, USR.27 | Control schema; no blockers, placed late by priority. |
| USR.1 | User accounts | USR.2, USR.3, USR.4, USR.5, USR.6, USR.7 | Parent; closes with its subitems. |
| USR.8 | Every signed-in user can generate a one-off system | USR.26 | Boss 05:16Z. Today /admin/generate/system is admin-only and only admin accounts exist, so no earlier step is needed. Open question: per-user limit (default 30 an hour, admins unlimited). |
| USR.9 | `seo.privacy_note` text, "download my data" and "delete my account" on the account page | USR.7, USR.26 | Research suggestion; not requested. |
| USR.11 | A "home" field on user accounts: where the Galaxy Map zooms to when it opens | MAP.189, MAP.190 |  |

### Halo mass research

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.399 | Research the halo mass function and compare it with the halo mass approximation GEN.368 uses | GEN.368 | Boss 10:11 PM EDT; research and comparison; Core 1 proposed (it writes GEN.368). |
