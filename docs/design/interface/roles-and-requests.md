# Roles and requests: users ask, moderators approve, admins approve moderators

**Status (2026-10-16):** plan; not built. This is the design of the 10.0 headline (USR.25). It builds on the accounts of 8.2 (USR.1 to USR.7) and the change plans of 8.6 (ADM.57, ADM.58 to ADM.68).

Boss (2026-10-16, 1:43 PM EDT): "Admins don't need approval. All other defaults fine, but add for v10.0 user edit requests and a moderator role, so User, Moderator, Admin. Moderators would be able to approve user requests, admins approve moderator requests. So we'll have to build out the plan for v10.0 to add roles, permissions, etc." Boss (1:49 PM EDT): "basically admins can approve themselves at the time, so they don't request as much as just confirm".

## 1. The idea in one paragraph

A request is a change plan of 8.6 that someone below the approving role submitted. The plan engine computes its effect, the approver sees the same preview an admin sees, and approving it applies it through the same job, journal and undo. Which role may submit which operation, and which role must approve it, are data in one permission catalogue and one `can()` check, so the rules live in one place.

## 2. Roles and who decides

| Role | May request | May approve | Own plan |
|---|---|---|---|
| User | text aspects of one object (default: rename, description, note) | nothing | needs a moderator or an admin |
| Moderator | everything a user may, plus a move, an orbit edit or a regeneration of one object | requests from users | needs an admin |
| Admin | everything | requests from users and moderators | confirm it after the preview, no second approver |
| Owner | everything an admin may | the same, plus the Owner-only settings tier (ADM.86) | the same as an admin |

Nobody decides their own request unless they are an admin confirming their own plan. A request that touches a locked aspect, a named object or a populated world goes to an admin whoever asked (default). The lists above are defaults held as data in the catalogue (USR.12), so Boss can widen them without a code change.

## 3. The parts

| Part | What it is | Items |
|---|---|---|
| Permission catalogue | Every permission the site checks, the role-to-permission map and `can(account, permission, object)` | USR.12 |
| Moderator role | Assign and remove it on the accounts page, audited | USR.13 |
| Guard sweep | Every route, page and command uses `can()`; a test finds a raw role comparison | USR.14, TEST.179 |
| Request tables | Comments and notice records; the request itself is a plan row of 8.6 (DB.36) | DB.39 |
| Submit, queue, decide | Submit as a plan in state `requested`; each role sees its share; decide by role | USR.15, USR.16, USR.17 |
| Approver screen | The plan preview of ADM.66 with the requester's note and the decision buttons | USR.18 |
| Conversation | Comments, edit and resubmit, withdraw, "my requests" | USR.19, USR.20 |
| Lifecycle and limits | Stale, withdrawn and expired requests; per-user caps and a block | USR.21, USR.22 |
| Notices and audit | Email notices; one audit view of role changes and requests | USR.23, USR.24 |
| Proof | End-to-end test of the chain, the permission matrix test, the help page | TEST.178, TEST.179, DOC.55 |

## 4. What 8.2 and 8.6 hold so nothing is rebuilt

- **8.2, accounts (USR.2).** The role is one text column with four valid values (user, moderator, admin, Owner); moderator is valid and unused until 10.0. Every admin-only page, route and command goes through one guard (`require_role`), which 10.0 replaces with `can()`; nothing else compares a role. The role change page writes its audit row through one function that 10.0 reuses.
- **8.6, change plans (DB.36, ADM.58, ADM.60, ADM.66, ADM.68).** The plan table carries the requester and their role, the role that must approve, the decision and note, and the state list already has `requested` and `withdrawn`. Apply runs under the approver's identity and the journal names both people. The approval screen shows the requester's note when there is one. The history lists requester and approver.
- **8.1, settings.** The limits (open requests, daily submissions, expiry) are options of the settings page (ADM.84).

## 5. Decisions and their reasons

- **A request is a plan, not a second mechanism.** Preview, apply, undo, history and audit are written once in 8.6; 10.0 adds who may submit and who may decide.
- **The catalogue is data.** A new role or a widened list is an edit to data and a test, not a change to routes.
- **Admins confirm, they do not request.** Boss decided an admin's own edit has no second approver; the preview stays so nothing applies unseen.
- **No self-approval below admin.** A moderator cannot approve their own request; only an admin can.
- **Stale means resubmit.** If an object changed after submission the request goes back to the requester, as stale plans do in 8.6.

## 6. Phases and threads (version 10.0)

| Phase | Theme | Threads |
|---|---|---|
| 1 | Roles, permissions and the request tables | Permissions; Request tables |
| 2 | Submitting and deciding requests | Request flow |
| 3 | Reviewing, limits, notices and audit | Reviewing; Lifecycle; Notices and audit |
| 4 | Close-out: proof, help and the headline | Proof and help |

## 7. Open questions (each with a default; work goes ahead on the default)

- **What users may request.** Default: rename, description and note on one object; no structural change. Boss may widen it.
- **What a moderator may approve.** Default: only user requests inside the user list; a locked, named or populated target goes to an admin.
- **Limits and expiry.** Default: 5 open requests, 20 submissions a day, 30 days to expiry; settings on the settings page.
- **Notices.** Default: approvers and requesters are emailed, each account can turn the emails off.
