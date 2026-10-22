# Version 10.0: Roles and permissions

Filed on 2026-10-16 as a version map (Boss, 1:43 PM EDT: three roles, User, Moderator and Admin; user edit requests; moderators approve user requests; admins approve moderator requests; blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 10.0

## Goal

Three roles, User, Moderator and Admin (the Owner stays above admin): users submit edit requests, moderators approve user requests, admins approve moderator requests, and an admin's own edit needs no approval, only the preview and confirm of the plan (USR.25, Boss 1:43 PM and 1:49 PM EDT). A request is a change plan of 8.6, so it keeps the preview, the apply job, the history and the undo; the roles of 8.2 supply the accounts and one permission check `can()` replaces every admin-only guard. Design: [docs/design/roles-and-requests.md](../design/roles-and-requests.md). 10.0 is a major version and needs Boss's word for its `.major.md` note.

## Phase 1: Roles, permissions and the request tables

The permission catalogue, the moderator role, the sweep that moves every guard to it, and the one migration.

### Permissions

| ID | Item | Needs | Note |
|---|---|---|---|
| USR.12 | Permission catalogue and one `can()` check: user, moderator, admin and Owner as data | USR.2 | Roles from 8.2 (USR.2 already holds the moderator value) |
| USR.13 | Moderator role: assign and remove it on the admin accounts page, with the audit entry | USR.12, USR.27 |  |
| USR.14 | Move every route, page and command from "is admin" checks to permission checks | USR.12 |  |

### Request tables

| ID | Item | Needs | Note |
|---|---|---|---|

## Phase 2: Submitting and deciding requests

A user or moderator submits an edit as a plan; the queue shows each person their share; the role above decides.

### Request flow

| ID | Item | Needs | Note |
|---|---|---|---|
| USR.15 | Submit an edit request: a user proposes an edit and the plan engine computes what it would change | USR.12, ADM.58, ADM.59 |  |
| USR.16 | Request queue: who sees what | USR.15 |  |
| USR.17 | Approve, reject or send back a request, by role: moderators decide user requests, admins decide moderator requests | USR.16, ADM.60 |  |

## Phase 3: Reviewing, limits, notices and audit

The approver screen, comments, the request lifecycle, limits, email notices and the audit view.

### Reviewing

| ID | Item | Needs | Note |
|---|---|---|---|
| USR.18 | Approver screen: the requester's note, the plan preview and the decision buttons | USR.17, ADM.66 |  |
| USR.19 | Comments on a request between the requester and the approver | USR.17 |  |
| USR.20 | "My requests" in the user interface: state of each request and withdraw | USR.16, USR.7 |  |

### Lifecycle

| ID | Item | Needs | Note |
|---|---|---|---|
| USR.21 | Request lifecycle: stale, withdrawn and expired requests | USR.17 | 30 days is the default |
| USR.22 | Request limits: a per-user cap on open requests, a rate limit, and blocking a user from requesting | USR.15 | Limits are defaults |

### Notices and audit

| ID | Item | Needs | Note |
|---|---|---|---|
| USR.23 | Email notices: a request was submitted, decided or sent back | USR.17, USR.3 |  |
| USR.24 | Audit view for roles and requests: who changed a role, who asked, who decided | USR.17, ADM.68 |  |

## Phase 4: Close-out: proof, help and the headline

End-to-end and matrix tests, the help page and the headline.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every point release.

### Proof and help

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.178 | Roles end to end: user, moderator and admin through the whole request chain | USR.18, USR.19, USR.20, USR.21, USR.22, USR.23, USR.24 |  |
| TEST.179 | Permission matrix test: every route and page checked for every role against the catalogue | USR.14 |  |
| DOC.55 | Help: roles, edit requests and how an approval works | TEST.178 |  |
| USR.25 | Headline of 10.0: roles and permissions, user edit requests and moderators | USR.13, USR.14, TEST.178, TEST.179, DOC.55 | Headline; a major version needs Boss's word for the .major note |
