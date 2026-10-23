# Release roadmap and version map, v8.1 to v10.0

Boss (2026-10-16, 1:07 PM EDT) turned the tiers into a version map: one version per release, starting at 8.1, which is the Phase 1 we are building today, and a version can have several phases when it holds too much for one. A version is only a PROPOSAL where its Status line says so; the headlines Boss named are used as given.

Fifteen versions from 8.1 to 10.0, then a "later" list for features that belong to no version yet. Each has one headline (8.1 has two), its phases, the open items it holds, what it needs, and how much of it is finished. Boss named 8.1 (infinite zoom, and the new generate scatters with accurate galactic generation and scattering including the halo), 8.2 (accounts with email invitations sent by an admin only), 8.5 (courses that bend around gravity wells), 8.6 (edit everything), 9.0 (composition and atmosphere colour), 9.1 (the view from a planet), 9.4 (galaxy types: different galactic models for generation, based on known galaxy types) and 10.0 (roles and permissions: user edit requests, moderators and admins). Politics and territories come before ships and stations, and many galaxies is a possible future feature, not a headline.

## How to read this

- A version is a REVISION bump (a `minor` note) and 9.0 is a `major` note, which needs Boss's word. The stamp Action numbers each merge as MAJOR.REVISION.BUILD, so a "version" here is the headline a REVISION or MAJOR earns, and the 8.0.15xx builds are the revisions inside the current line. While `REVISION_HOLD` is on, nothing stamps 8.1 until Phase 1 (this page's 8.1) is complete and Boss says so (Boss, 2026-10-09).
- 8.0 is the old Phase 0 (finished, 241 of 241). Every later item sits in exactly one version, one phase of it and one thread of that phase.
- Phases (Boss, 1:09 PM EDT): a blocker is done in the phase before, in an earlier thread position of the same thread, or in an earlier version. A version has several phases when one is too much, each with a theme. The threads inside one phase are independent work trees, so lanes can build them side by side. The page flags any item that breaks this order; there are none today.
- Release gate (Boss, 8:01 PM EDT 2026-10-22, confirming 5:10 PM EDT 2026-10-21): the full security code audit (SEC.33) and the full bug hunt audit (TEST.261) each need a dated report, with no blocking finding open, before the version number changes at its major or minor part (8.1, 8.2 and on to 9.0 and 10.0). Build numbers (8.0.NNNN) are exempt.
- Progress per version is finished over finished plus open, from the progress ledger (the version each finished item belonged to) and TODO.md on every build. Dropped items count neither way.

## Release 8.1: Infinite zoom and the new generate scatters
Status: in progress, the only version being built (Boss, 7:44 PM EDT 2026-10-22: the five core lanes work 8.1 items only until 8.1 is complete). As of 2:02 AM EDT 2026-10-23: 701 of 936 items finished (74.9 percent), 235 open, by phase 113 (foundations and open trees), 29 (infinite zoom), 28 (stars, nebulae, region query), 17 (components, clusters, audits), 33 (core, halo, interface, docs) and 15 (the headline closes). The DB lane has finished its queue (DB.85 merged as #2582) and its thread is retired.
Headline: Infinite zoom, one camera from the whole galaxy to a star system and on into its planets and moons, and the new generate scatters: accurate galactic generation and scattering out to 262 kpc, including the halo, the satellites, the nebula stars and a realistic core
Today: 8.1 infinite zoom (the fly-through Galaxy Map), MAP.146
Change: A second headline joins the first (Boss, 1:07 PM EDT): the new generate scatters, as GEN.301 (new, carries a `Headline: 8.1` tag). The Phase 1 tail (queues, reproducible galaxies, the description layer, tech levels, directives, every open bug and flaky test) stays in 8.1 as phase 1, because the stamp waits for all of it. The admin editors that sat in the Phase 1 tail (ADM.32, ADM.35, ADM.55, ADM.56, MAP.119, MAP.120) move to 8.6 and the bookmark items (UX.46, UX.47) to 8.2, so they no longer gate the 8.1 stamp. The configuration page (ADM.43, ADM.44) moves up from 8.2 into 8.1 phase 1 as its own thread (Boss, 1:43 PM EDT: "we're starting to desperately need that"), with the per-cell system cap (GEN.303) wired to it in phase 5. The settings page also holds the one backup mechanism, a bz2 tarball per backup with `settings.backups_kept` archives kept (ADM.87), and its first part, the settings (Boss, 1:59 and 2:00 PM EDT). The research, designs and data tables of the 9.0 composition headline are building blocks that go down to 8.1 (Boss, 11:57 PM EDT): phase 1 holds the four research lanes and the designs (GEN.431 to GEN.437, with GEN.409, GEN.425 and GEN.426), phase 2 the data tables (GEN.438, GEN.410, GEN.427, GEN.439), phase 3 the star composition function (GEN.440) and phase 4 its sweep (TEST.211); none changes a stored value.
Needs: the closing order of 8.1 (Boss, 6:40 PM EDT 2026-10-21, updated 10:46 PM EDT 2026-10-22): the whole-sphere audits TEST.185, TEST.184, TEST.159 and TEST.160 (merged in #2504 and re-run on main after the PERF.146 speedups as TEST.283, all passing), the performance profile PERF.146 (merged as #2559, with PERF.163 and its follow-ups PERF.164 to PERF.166), then the security audit SEC.33 and the bug hunt audit TEST.261 on the final code in their own threads (required for every major or minor version change, never for a build number), then the headline GEN.301. MAP.146 (infinite zoom) runs in Core 2's queue after the Generate form and menu bar bugs ADM.127 and UX.146 and the map selection bugs MAP.210, TEST.281, MAP.209 and MAP.249. The plan numbers bugs (GEN.626, GEN.628, GEN.369, GEN.624, GEN.403, GEN.627) run in the Bugfix 2 lane, and GEN.625 (the Auto defaults) in Core 1 waits for GEN.626 and GEN.628. Generator epochs are numbered at merge (53 is the next free number, GEN.307 took 51 and GEN.315 slice 1 took 52; the epochs file must stay contiguous).
Lanes (10:20 AM EDT 2026-10-23, queues in reports/8-1-core-queues-2026-10-22.md and the board in docs/plan/lanes.md): Core 1: GEN.633 (before 8.1, Boss 8:29 AM EDT), GEN.55, then reproducibility, storage, superbubbles and satellites; Core 2: ADM.126, MAP.174, TEST.237, UX.124 (MAP.209, MAP.210 and MAP.249 wait on Boss's retest; MAP.146 is done and retired); ADM.125 is merged (#2646); Core 3: GEN.630, GEN.631 (the split of GEN.316), GEN.332, the GEN.315 remainder and the multi-star, nebula and galaxy-law items; Core 4: GEN.450 (card with Boss) and the binaries and clusters; Core 5: PERF.123 (slice 2 merged as #2665; D1, then C, then D2), the queue and logging items and GEN.301 last; Bugfix 2: finished (all six plan numbers bugs merged; its follow-ups are GEN.636 for Core 4, TEST.294 for Bug Hunt and GEN.637, the GEN.625 sweep re-run); Blocker 2 is finished (GEN.599, GEN.416, GEN.498, GEN.496, GEN.470 and GEN.452 part 1 merged; GEN.452 part 2 is PR #2613, held for the 8.4 epoch); Blocker 3: GEN.458 (#2655) and SEC.32 (#2667) merged, then GEN.146; Bug Hunt: TEST.291 (next), TEST.284 and TEST.290 (waiting for a recurrence), TEST.287, the stable-code hunt, then the TEST.261 audit at the close; Performance: PERF.164, PERF.165 and the findings PERF.156 to PERF.162; PERF.154 (waiting list) goes to Core 5 last and PERF.155 is in the pool. Finished lanes: Fair scheduling (retired; PERF.150 to PERF.153 merged), the DB lane (retired), Admin suite, Galaxy dials.

## Release 8.2: Accounts and email invitations
Status: Boss named it; new version Progress (2:02 AM EDT 2026-10-23): 15 of 41 items finished (36.6 percent), 26 open; paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22).
Headline: User accounts, with sign-up by an email invitation that only an admin can send
Today: (new)
Change: Boss, 12:45 PM EDT. It is already planned as USR.1 to USR.7; the configuration page left this version for 8.1 (Boss, 1:43 PM EDT) and USR.2 stores four roles (user, moderator, admin, Owner) behind one guard so the roles of 10.0 fit without a rebuild; this moves it up and puts the bookmark items (UX.46, UX.47) with it. One change to USR.4 on approval: its decided answer says the invitee's email is optional and the admin may only copy the link; the email invitation becomes the main path, with the link shown only when SMTP is not set up.
Needs: the settings page (ADM.43, ADM.84), now built in 8.1, is where the SMTP form sits; its Owner tier (ADM.86) opens this version's thread after USR.2; SEC.32 (Argon2id) before USR.2 so new accounts start on it. Mail goes out as a queue job (USR.3), so the queue steps PERF.97 to PERF.100 in 8.1 phase 1 come first.

## Release 8.3: The moving galaxy
Status: proposed; the orbital update moves here from 8.0 Progress (2:02 AM EDT 2026-10-23): 33 of 67 items finished (49.3 percent), 34 open; paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22).
Headline: The orbital update: every object moves on its orbit each day, with n-body influence, collisions and ejections
Today: 8.0 the orbital update system (GEN.105)
Change: GEN.105 needs open items that take database migrations from Core 1 (GEN.289, GEN.294), so it earns its own version after the scatters. GEN.198 (the gravity field evaluator) is placed here, not with the gravity map, because the courses of 8.5 reuse it.
Needs: GEN.294 and GEN.289 (every object as one SpatialPosition3D and kept planets and moons) first, in phase 1; GEN.115 and GEN.142 have no prerequisites. GEN.111 needs 8.2's SMTP settings (USR.3). GEN.299 and GEN.300 (the Quantity class) need GEN.294 and GEN.135. OPS.16 and OPS.17 set up the nightly run.

## Release 8.4: New worlds
Status: proposed Progress (2:02 AM EDT 2026-10-23): 30 of 140 items finished (21.4 percent), 110 open; paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22).
Headline: Planet classes R to Z and the habitability-driven refactor: a world's class, atmosphere and life follow its index
Today: (new)
Change: The planet class work (GEN.33, GEN.28, GEN.27, GEN.90 and the rest) has no headline today. Boss (11:25 PM EDT) put the generation and storage half of the atmosphere colour here, as four phases ahead of the class work: the bulk composition that drives density, mass and gravity (GEN.428 to GEN.430, DB.69), the albedo records (GEN.412), the surface and small-body composition (GEN.441, GEN.442), the one mixture and its stored values (GEN.411, GEN.422 to GEN.424, DB.67), ozone, haze, dust and cloud (GEN.413, GEN.414), the stored shares (DB.70, DB.68), the text and the sweeps (GEN.415, TEST.209, TEST.210, DOC.63, DOC.65). The research, designs and tables they read moved to 8.1 (Boss, 11:57 PM EDT).
Needs: The habitability index (built). GEN.33 first, one class per PR. GEN.95 comes before GEN.93 (nebula conditions on planets). POP.13 (the species editor) goes with it because the "find me a world" idea needs species limits.

## Release 8.5: Courses that bend around gravity wells
Status: Boss named it; headline unchanged Progress (2:02 AM EDT 2026-10-23): 14 of 43 items finished (32.6 percent), 29 open; paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22).
Headline: Courses that bend around gravity wells, with waypoints, editable and saved courses (NAV.59)
Today: 8.5 courses that bend around gravity wells, NAV.59
Change: No change to the headline. NAV.19 (saved courses in the account) rides here because accounts arrive in 8.2.
Needs: MAP.122 (Select mode) from 8.1 for waypoints (NAV.49); GEN.198 (the evaluator) from 8.3; saved courses in the account need USR.7 from 8.2, while per-browser saving (NAV.18) needs nothing. NAV.6's open question on keep-out radii still stands, with its default.

## Release 8.6: Edit everything
Status: Boss named it (1:07 PM EDT); items filed Progress (2:02 AM EDT 2026-10-23): 4 of 47 items finished (8.5 percent), 43 open; paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22).
Headline: Edit everything: the admin has complete control and the world adjusts to what the admin wants, after the admin approves a generation plan
Today: (new)
Change: New headline (ADM.57), filed as 35 items in 5 phases with the design in docs/design/edit-everything.md. Boss: the admin can move stars in a sector, change orbits in a system, change densities for a sector or region and regenerate any star to any specification; changing an A planet into an M planet makes the system work out how to rearrange the star system and show that generation plan for approval before anything changes. In every regen screen the admin picks the aspects to keep and the aspects to regenerate. Every change is a plan the admin sees and approves, applies as one job and can undo. An admin approves their own plan, so an admin confirms rather than requests (Boss, 1:49 PM EDT); the plan states already include `requested` so the roles of 10.0 need no new plan table. Two pieces are shared with 8.9 and built here once (Boss, 1:33 PM EDT): the unit finalizer (ADM.82, also API.8) and the recipe schema (API.28, also API.18).
Needs: The plan model, aspects and regions come first (phase 1); then the approval screen, locks, regen screens and density overrides (phase 2); then the editors that make plans and the rearranger (phase 3); then regeneration to a specification (phase 4); then the end-to-end tests, the help and the headline (phase 5). ADM.36 (edit orbits) comes from the orbital update foundations in 8.3 and ADM.52 (the planet editor) from earlier. The override journal has to keep regeneration, parity repair and the golden-seed tests (GEN.57, DB.9, DB.17 in 8.1 phase 1) working.

## Release 8.7: The pilot's console and built-in help
Status: proposed Progress (2:02 AM EDT 2026-10-23): 0 of 17 items finished (0.0 percent), 17 open; paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22).
Headline: A built-in help section and the pilot's starmap look: every page explains itself and the interface reads like a navigation console
Today: (new)
Change: The help pages (DOC.6 to DOC.16) and the console design (UX.43) have no headline today.
Needs: The help pages describe what 8.1 to 8.6 built, so this version goes after them; DOC.6 (the template and coverage test) first.

## Release 8.8: The deep sky
Status: proposed Progress (2:02 AM EDT 2026-10-23): 1 of 33 items finished (3.0 percent), 32 open; paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22).
Headline: Nebulae, clusters, satellites and the rarest stars drawn on the map and reachable by flight
Today: (new)
Change: The nebula drawing items (MAP.176 to MAP.183), the cluster work and the "later" stellar phenomena have no headline today.
Needs: The scatters of 8.1 (components, nebula stars, Wolf-Rayet, scattered nebulae in the index) and the nebula tiles and cover (MAP.173, MAP.174). MAP.194 (nebulae and clusters open their own zoom window) needs MAP.176.

## Release 8.9: Generate anywhere
Status: proposed Progress (2:02 AM EDT 2026-10-23): 4 of 35 items finished (11.4 percent), 31 open; paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22).
Headline: Remote generation: build a galaxy on a local machine, upload it through the API and have the server verify it
Today: (new)
Change: The remote generation chain (API.3 and its parts) and recipes (API.18, API.19) have no headline today. Boss (1:33 PM EDT): 8.6 and this version are closely related, so what both need is built once in 8.6: the unit finalizer (ADM.82, used by API.8) and the recipe schema and runner (API.28, used by API.18); API.8 and API.18 here keep only the upload and route parts, and API.19 uses the region object of 8.6 (GEN.302). It depends only on the reproducibility items of 8.1 phase 1 and those 8.6 pieces, so it can move earlier if Boss wants it sooner. Also here (Boss, 11:26 PM EDT): the sky data and research of the view from a planet (VIEW.4, VIEW.6, VIEW.9, DB.64, VIEW.1, VIEW.7, VIEW.8) and the backups and repair chain (ADM.88 to ADM.91, DB.40 to DB.42, OPS.55, TEST.180, DOC.56), moved from 9.0 because none of it needs the atmosphere colour.
Needs: GEN.57, GEN.135, OPS.28 and TEST.77 (reproducibility, 8.1 phase 1) first: remote generation must reproduce what the server makes (API.17). API.10, API.11 and API.12 are the base; API.19 needs API.18.

## Release 9.0: Composition and atmosphere colour
Status: proposed; Boss named it the 9.0 feature (11:25 PM EDT) and made composition its joint headline (11:57 PM EDT) Progress (2:02 AM EDT 2026-10-23): 1 of 59 items finished (1.7 percent), 58 open; paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22).
Headline: What every body, star and atmosphere is made of, from real science, and what follows from it: the colour of the air of any world, for solid worlds and gas giants, and the radiation in the air and on the ground, with every stored value agreeing
Today: (new)
Change: Boss (11:24 PM EDT, 11:25 PM EDT): store the albedo and greenhouse values, break the molar density into per-gas partial pressures so both match, build ozone from real research, add haze, dust and cloud, keep one table of every gas and compound, and word the sky colour as atmosphere colour, absolute and including gas giants; it replaces the planet colour on the maps. The view from a planet that was planned here moved to 9.1 and the backups and repair chain to 8.9 (11:26 PM EDT: "Everything for 9.0 before splits between 9.1 and 8.9"). Boss (11:57 PM EDT): atmosphere composition and the composition of planets, moons, stars and everything is the joint headline of 9.0: every body has figures on its composition from actual science, not made up where it can be avoided, with the atmospheric radiation and the surface radiation derived from what is in the air and the ground; two research lanes (atmosphere, body composition) write the design documents first.
Needs: The building blocks go down to earlier versions (Boss, 11:57 PM EDT): 8.1 phase 1 holds the research lanes and the designs (GEN.431 to GEN.434 the surveys of atmospheres, bodies, stars and radiation, GEN.425, then the designs GEN.409, GEN.426, GEN.435, GEN.436 and the master model GEN.437), 8.1 phase 2 the data tables (the elements GEN.438, the gases GEN.410, the materials GEN.427, the stars GEN.439), 8.1 phase 3 the star composition function (GEN.440) and phase 4 its sweep (TEST.211). 8.4 phases 1 to 4 hold the generation and storage: the bulk composition that drives mass, density and gravity (GEN.428 to GEN.430, DB.69), the albedo records (GEN.412), the ground and the small bodies (GEN.441, GEN.442), the mixture and stored values (GEN.411, GEN.422 to GEN.424, DB.67), ozone, haze, dust and cloud (GEN.413, GEN.414), the stored shares (DB.70, DB.68), the description text (GEN.415) and the sweeps and docs (TEST.209, TEST.210, DOC.63, DOC.65). 9.0 holds the colour model (GEN.416 to GEN.419), the radiation (GEN.443 to GEN.445), what shows them (GEN.420, GEN.446, MAP.244, DOC.64), the radiation sweep and the docs (TEST.212, DOC.66), with the headline item GEN.421. The migrations come after the 0102 to 0105 stack of 8.1.
Alternative (Politics and territories): Territories that grow and shift on the map, with borders drawn and polities that trade. Placed at 9.2; it could be a later headline.

## Release 9.1: The view from a planet
Status: proposed; Boss chose the headline (11:42 PM EDT) Progress (2:02 AM EDT 2026-10-23): 1 of 19 items finished (5.3 percent), 18 open; paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22).
Headline: Stand on a world and look up: the sky from any planet or moon, with real constellations, rendered as a picture
Today: (new)
Change: Boss (11:26 PM EDT, 11:42 PM EDT): the view from a planet moved here from 9.0 and is the headline of 9.1, with the gravity map (MAP.167) as the second feature. Its sky data and research are in 8.9 (VIEW.4, VIEW.6, VIEW.9, DB.64, VIEW.1, VIEW.7, VIEW.8) and it builds on the atmosphere colour of 9.0 (VIEW.11 extends GEN.416). The gravity map keeps its foundations in 8.3 (GEN.198, TEST.131).
Needs: 8.3 (so the sky changes with the date), 8.4 (the worlds) and the atmosphere colour of 9.0. VIEW.2 starts with a research session with Boss; VIEW.6 settles the handedness of the galaxy before real sky coordinates. The gravity map needs GEN.199 and PERF.72 (the per-zone grid and its cost), API.24, MAP.168, MAP.169 and UX.94; MAP.170 (the coarse Galaxy Map layer) stays later and optional.

## Release 9.2: Politics and territories
Status: Boss ordered it before ships and stations; idea only, no items
Headline: Territories that grow and shift over time on the map, with borders drawn and polities that trade
Today: (new)
Change: Boss, 1:07 PM EDT: politics and territories come before ships and stations. No items exist yet beyond the population-and-politics design note; the item breakdown is for a later change after Boss's answers.
Needs: 8.3 (time) and the population items of 8.4 (species and limits). The population-and-politics design note holds the groundwork.

## Release 9.3: Ships and stations
Status: proposed; after politics and territories Progress (2:02 AM EDT 2026-10-23): 1 of 7 items finished (14.3 percent), 6 open; paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22).
Headline: Users place their own ships and stations and follow them across the live map
Today: (new)
Change: USR.10 and POP.11 have no headline today; they follow politics and territories (Boss, 1:07 PM EDT).
Needs: 8.2 (accounts), 8.3 (movement), 8.5 (courses to fly) and 9.2 (territories, so stations belong to a polity).

## Release 9.4: Galaxy types
Status: Boss named it (8:48 PM EDT); items filed Progress (2:02 AM EDT 2026-10-23): 0 of 16 items finished (0.0 percent), 16 open; paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22).
Headline: Galaxy types: the seed draws a galaxy of a known type (spiral, lenticular, elliptical, irregular or dwarf), every shape number calculated from published scaling relations
Today: (new)
Change: New headline (Boss, 8:48 PM EDT: "different galactic models for generation based on known galaxy types", broken into a chain of required items and counted back from 9.4), filed as 18 items in 6 phases with the plan in docs/design/galaxy-types.md. A galaxy type is a parameter family, never a copy of a named galaxy (no fixed content, Boss 4:11 and 4:14 PM EDT): the seed draws the stellar mass, then the family from the mass-dependent type fractions, then each number from a scaling relation. It reuses the shape, the density model, the potential, the halo, the satellites and the scatters. The old "Hubble-sequence generator" suggestion (SL-b) is this headline. 9.4 is Boss's own number and needs no .major note.
Needs: three steps of the chain come earlier so the families only fill in the shape: GEN.378 (the checked design and its scaling relation table, after GEN.368) in 8.1, then DB.66 (the family key and component columns, a migration) and GEN.379 (the density model as a list of components, which the potential GEN.338 reads) in 8.3. Then, in 9.4: the type draw and the Sérsic spheroid (phase 1); the spiral, elliptical and irregular and dwarf families (phase 2); the lenticular family, populations, clusters with the nucleus (phase 3); creation (phase 4); the Generate page, the Galaxy Map, the plausibility tests and the docs (phase 5; each family also gets a simple vector outline, UX.135, in phase 1); the headline (phase 6). It also reads GEN.368, GEN.369, GEN.297, GEN.134 and GEN.319.
Lanes: not assigned; GEN.378 is proposed for Core 1 (it writes GEN.368's table), and the rest sits after 8.3.

## Release 10.0: Roles and permissions
Status: Boss named it (1:43 PM EDT); items filed Progress (2:02 AM EDT 2026-10-23): 1 of 18 items finished (5.6 percent), 17 open; paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22).
Headline: Roles and permissions: users submit edit requests, moderators approve user requests and admins approve moderator requests
Today: (new)
Change: New headline (USR.25), filed as 18 items in 4 phases with the design in docs/design/roles-and-requests.md. Three roles, User, Moderator and Admin, with the Owner above admin. A user proposes an edit, a moderator approves it; a moderator proposes an edit, an admin approves it; an admin approves their own plan at the time, so an admin confirms with the preview and does not request (Boss, 1:49 PM EDT). A request is a change plan of 8.6, so approval, preview, history, undo and audit are written once. 10.0 is a major version; the `.major.md` note waits for Boss's word.
Needs: 8.2 (accounts: USR.2 holds the role column with the four values and one guard, USR.3 mail, USR.7 the user interface) and 8.6 (change plans: ADM.58 engine with the `requested` state, ADM.59 store, ADM.60 apply, ADM.66 approval screen, ADM.68 history). Then the permission catalogue and the guard sweep, the request tables (phase 1); submit, queue and decide (phase 2); the approver screen, comments, lifecycle, limits, notices and the audit view (phase 3); the end-to-end and matrix tests, the help page and the headline (phase 4).
Lanes: not assigned; the work sits after 8.6, so no lane queue holds it.

## Release later: Future features
Status: possible future features, not scheduled Progress (2:02 AM EDT 2026-10-23): 1 of 8 items finished (12.5 percent), 7 open; paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22).
Headline: No headline: features that belong to no version yet, noted for later
Today: (new)
Change: Boss, 1:07 PM EDT: multiple galaxies is a possible future feature, not a headline (Boss, 1:43 PM EDT: 10.0 is roles and permissions). The items (GEN.9, GEN.157, GEN.164, MAP.145, NAV.57) sit in the "later" list.
Needs: GEN.9 (the plan for more than one galaxy in the database) first; GEN.157 (a verified table of about 25 real galaxies) before MAP.145 and NAV.57; GEN.164 needs GEN.9. The infinite zoom of 8.1 would extend past the galaxy edge.

## Questions for Boss
- Decided (Boss, 1:29 PM EDT, "defaults are acceptable"): keep the rule that 8.1 stamps only when all of its phases are done; ADM.32, ADM.35, ADM.55, ADM.56, MAP.119 and MAP.120 sit in 8.6 and UX.46 and UX.47 in 8.2; 8.3 to 8.9 as proposed; 9.0 is the view from a planet; multiple galaxies stays under "later".
- 8.6 undo: an applied plan can be undone while no newer plan touched the same objects. Default: yes.
- Decided (Boss, 1:43 PM and 1:49 PM EDT): an admin needs no approval; an admin confirms their own plan after the preview and there is no second approver. Any admin may also confirm a plan another admin proposed.
- 8.6 density plans: the planner may remove generated, unnamed, unpopulated systems inside an approved plan, listed with counts per sector. Default: yes.
- 8.6 specifications: forced values the generator would not draw are allowed and marked on the plan as forced. Default: yes.
- 10.0 what users may request. Default: a user may request a rename, a description or a note on one object, and no structural change; a moderator may request those plus a move, an orbit edit or a regeneration of one object; an admin may do everything. Boss may widen the user list.
- 10.0 what a moderator may approve. Default: only user requests inside the user's list; a request that touches a locked aspect, a named object or a populated world goes to an admin whoever asked.
- 10.0 request limits and expiry. Default: 5 open requests and 20 submissions a day per user; a pending request expires after 30 days; both are settings on the settings page.
- 10.0 notices. Default: approvers and requesters are emailed, each account can turn the emails off.
