# Release roadmap and version map, v8.1 to v10.0

Boss (2026-10-16, 1:07 PM EDT) turned the tiers into a version map: one version per release, starting at 8.1, which is the Phase 1 we are building today, and a version can have several phases when it holds too much for one. A version is only a PROPOSAL where its Status line says so; the headlines Boss named are used as given.

Fifteen versions from 8.1 to 10.0, then a "later" list for features that belong to no version yet. Each has one headline (8.1 has two), its phases, the open items it holds, what it needs, how much of it is finished and a list of ideas that are not filed yet. Boss named 8.1 (infinite zoom, and the new generate scatters with accurate galactic generation and scattering including the halo), 8.2 (accounts with email invitations sent by an admin only), 8.5 (courses that bend around gravity wells), 8.6 (edit everything) 9.4 (galaxy types: different galactic models for generation, based on known galaxy types) and 10.0 (roles and permissions: user edit requests, moderators and admins). Politics and territories come before ships and stations, and many galaxies is a possible future feature, not a headline.

## How to read this

- A version is a REVISION bump (a `minor` note) and 9.0 is a `major` note, which needs Boss's word. The stamp Action numbers each merge as MAJOR.REVISION.BUILD, so a "version" here is the headline a REVISION or MAJOR earns, and the 8.0.15xx builds are the revisions inside the current line. While `REVISION_HOLD` is on, nothing stamps 8.1 until Phase 1 (this page's 8.1) is complete and Boss says so (Boss, 2026-10-09).
- 8.0 is the old Phase 0 (finished, 241 of 241). Every later item sits in exactly one version, one phase of it and one thread of that phase.
- Phases (Boss, 1:09 PM EDT): a blocker is done in the phase before, in an earlier thread position of the same thread, or in an earlier version. A version has several phases when one is too much, each with a theme. The threads inside one phase are independent work trees, so lanes can build them side by side. The page flags any item that breaks this order; there are none today.
- Progress per version is finished over finished plus open, from the progress ledger (the version each finished item belonged to) and TODO.md on every build. Dropped items count neither way.
- "Suggested, not filed" lists are ideas, each with a code (S8.1-a and so on) so Boss can say "file S8.1-a and S8.2-c". Each one gets checked against the open items first, so nothing is filed twice.

## Release 8.1: Infinite zoom and the new generate scatters
Status: Boss named both headlines; next stamp
Headline: Infinite zoom, one camera from the whole galaxy to a star system and on into its planets and moons, and the new generate scatters: accurate galactic generation and scattering out to 262 kpc, including the halo, the satellites, the nebula stars and a realistic core
Today: 8.1 infinite zoom (the fly-through Galaxy Map), MAP.146
Change: A second headline joins the first (Boss, 1:07 PM EDT): the new generate scatters, as GEN.301 (new, carries a `Headline: 8.1` tag). The Phase 1 tail (queues, reproducible galaxies, the description layer, tech levels, directives, every open bug and flaky test) stays in 8.1 as phase 1, because the stamp waits for all of it. The admin editors that sat in the Phase 1 tail (ADM.32, ADM.35, ADM.55, ADM.56, MAP.119, MAP.120) move to 8.6 and the bookmark items (UX.46, UX.47) to 8.2, so they no longer gate the 8.1 stamp. The configuration page (ADM.43, ADM.44) moves up from 8.2 into 8.1 phase 1 as its own thread (Boss, 1:43 PM EDT: "we're starting to desperately need that"), with the per-cell system cap (GEN.303) wired to it in phase 5. The settings page also holds the one backup mechanism, a bz2 tarball per backup with `settings.backups_kept` archives kept (ADM.87), and its first part, the settings (Boss, 1:59 and 2:00 PM EDT).
Needs: MAP.152 first, then MAP.146 (Core 2's next two items, in that order, Boss 12:23 PM EDT). The scatters run in the lanes' queues (Core 1, Core 2 and Core 3) from GEN.237 to GEN.243, GEN.264 to GEN.275 and GEN.280 to GEN.287. Phase 1 is independent of the later phases, so lanes can work them side by side.
Lanes: Core 1: GEN.237, GEN.264, GEN.238 and the Phase 1 trees; Core 2: MAP.152, MAP.146, PERF.105, UX.118; Core 3: GEN.280, GEN.283, UX.116, GEN.286, GEN.287, GEN.245 and the rest of the scatter work.
Suggested:
- S8.1-a: An acceptance flight and frame budget: a scripted camera flight from the whole galaxy to a named moon (galaxy, sector, system, planet, moon) with a frame-rate floor, run in the browser tests. No open item covers the last two steps as one flight.
- S8.1-b: A scale bar and distance readout that change units with the zoom (light years, parsecs, AU, kilometres) from the shared ladders.
- S8.1-c: Names fade in with the zoom: labels for the largest nearby things only, so the view stays readable at every scale.
- S8.1-d: Shareable camera positions: a link that flies a visitor to the same place, built on object IDs and the bookmark format.
- S8.1-e: A guided tour: a scripted flight past notable places (Sol, the galactic core, a nebula, a globular cluster) a visitor can replay.
- S8.1-f: Touch and keyboard flight controls and a reduced-motion setting.
- S8.1-g: A census page for the whole sphere: counts of every star type and phenomenon by region, from the representation audit (TEST.160), public and cached.
- S8.1-h: Galaxy export and import: one file holds a galaxy (seed, settings, stored sectors) to back it up or move it to another server, built on the parity file (DB.9).
- S8.1-i: Named generation presets on the Generate page (small dwarf, globular-rich, large spiral; no Milky Way preset, Boss 4:11 PM EDT) saved as recipes, a step toward API.18.

## Release 8.2: Accounts and email invitations
Status: Boss named it; new version
Headline: User accounts, with sign-up by an email invitation that only an admin can send
Today: (new)
Change: Boss, 12:45 PM EDT. It is already planned as USR.1 to USR.7; the configuration page left this version for 8.1 (Boss, 1:43 PM EDT) and USR.2 stores four roles (user, moderator, admin, Owner) behind one guard so the roles of 10.0 fit without a rebuild; this moves it up and puts the bookmark items (UX.46, UX.47) with it. One change to USR.4 on approval: its decided answer says the invitee's email is optional and the admin may only copy the link; the email invitation becomes the main path, with the link shown only when SMTP is not set up.
Needs: the settings page (ADM.43, ADM.84), now built in 8.1, is where the SMTP form sits; its Owner tier (ADM.86) opens this version's thread after USR.2; SEC.32 (Argon2id) before USR.2 so new accounts start on it. Mail goes out as a queue job (USR.3), so the queue steps PERF.97 to PERF.100 in 8.1 phase 1 come first.
Suggested:
- S8.2-a: Two-factor sign-in (TOTP) for every account, required for admins; the Owner transfer in USR.6 already assumes it.
- S8.2-b: An invitation the admin can write: an optional personal message, resend, revoke, and a notice back to the admin when it is accepted.
- S8.2-c: A sign-in activity page: devices, last seen, and sign out everywhere.
- S8.2-d: An admin user list: search, disable, last sign-in and role history.
- S8.2-e: Bounce handling: mark an address undeliverable and tell the admin who invited it.
- S8.2-f: Opt-in emails for things a user follows, such as a Hill-radius overlap in a bookmarked system (needs 8.3's GEN.111).

## Release 8.3: The moving galaxy
Status: proposed; the orbital update moves here from 8.0
Headline: The orbital update: every object moves on its orbit each day, with n-body influence, collisions and ejections
Today: 8.0 the orbital update system (GEN.105)
Change: GEN.105 needs open items that take database migrations from Core 1 (GEN.289, GEN.294), so it earns its own version after the scatters. GEN.198 (the gravity field evaluator) is placed here, not with the gravity map, because the courses of 8.5 reuse it.
Needs: GEN.294 and GEN.289 (every object as one SpatialPosition3D and kept planets and moons) first, in phase 1; GEN.115 and GEN.142 have no prerequisites. GEN.111 needs 8.2's SMTP settings (USR.3). GEN.299 and GEN.300 (the Quantity class) need GEN.294 and GEN.135. OPS.16 and OPS.17 set up the nightly run.
Suggested:
- S8.3-a: A time control on every map: pick a date and see where objects were or will be, using the stated-span run (GEN.292).
- S8.3-b: An events feed: encounters, captures, ejections and collisions as a public timeline with a page per event, from the events table in GEN.143 and the object IDs of GEN.197.
- S8.3-c: Per-object history: the track an object has followed, drawn on its page and the map.
- S8.3-d: Stars that age in real time: a star that reaches its end changes state at the update that crosses its age (giant, white dwarf, remnant).
- S8.3-e: An ephemeris export: positions at stated dates as CSV or JSON through the API.
- S8.3-f: An update dashboard: when the last nightly update ran, how many objects moved and how many encounters it found.

## Release 8.4: New worlds
Status: proposed
Headline: Planet classes R to Z and the habitability-driven refactor: a world's class, atmosphere and life follow its index
Today: (new)
Change: The planet class work (GEN.33, GEN.28, GEN.27, GEN.90 and the rest) has no headline today.
Needs: The habitability index (built). GEN.33 first, one class per PR. GEN.95 comes before GEN.93 (nebula conditions on planets). POP.13 (the species editor) goes with it because the "find me a world" idea needs species limits.
Suggested:
- S8.4-a: A planet class gallery: one page per class with the rule that makes it, an example and a link to every example in the galaxy.
- S8.4-b: Planet surface maps: a generated map image per world from its class, gravity, water and temperature (the design note for rogue planet surfaces is a start).
- S8.4-c: Ring systems and ring shadows for giant planets; tidal heating and ocean worlds for moons.
- S8.4-d: Compare two worlds side by side on the habitability index rows.
- S8.4-e: A "find me a world" search: filter by class, index and a species' limits.
- S8.4-f: A planet sheet that publishes to the wiki like a system does today.

## Release 8.5: Courses that bend around gravity wells
Status: Boss named it; headline unchanged
Headline: Courses that bend around gravity wells, with waypoints, editable and saved courses (NAV.59)
Today: 8.5 courses that bend around gravity wells, NAV.59
Change: No change to the headline. NAV.19 (saved courses in the account) rides here because accounts arrive in 8.2.
Needs: MAP.122 (Select mode) from 8.1 for waypoints (NAV.49); GEN.198 (the evaluator) from 8.3; saved courses in the account need USR.7 from 8.2, while per-browser saving (NAV.18) needs nothing. NAV.6's open question on keep-out radii still stands, with its default.
Suggested:
- S8.5-a: Ship profiles (speed, fuel range, shielding) so a course shows its fuel and time for a chosen ship.
- S8.5-b: Hazard costs: radiation (GEN.256) and nebula density as route penalties, with a safest or shortest toggle.
- S8.5-c: Read-only course sharing by link, with no account needed.
- S8.5-d: Animated travel: play a course as a ship flying along it at a chosen speed, with the bends around the wells visible.
- S8.5-e: Suggested stops at systems with the facilities a ship needs (uses the tech level and facility data).
- S8.5-f: Compare courses: straight, bent and safest, side by side.

## Release 8.6: Edit everything
Status: Boss named it (1:07 PM EDT); items filed
Headline: Edit everything: the admin has complete control and the world adjusts to what the admin wants, after the admin approves a generation plan
Today: (new)
Change: New headline (ADM.57), filed as 35 items in 5 phases with the design in docs/design/edit-everything.md. Boss: the admin can move stars in a sector, change orbits in a system, change densities for a sector or region and regenerate any star to any specification; changing an A planet into an M planet makes the system work out how to rearrange the star system and show that generation plan for approval before anything changes. In every regen screen the admin picks the aspects to keep and the aspects to regenerate. Every change is a plan the admin sees and approves, applies as one job and can undo. An admin approves their own plan, so an admin confirms rather than requests (Boss, 1:49 PM EDT); the plan states already include `requested` so the roles of 10.0 need no new plan table. Two pieces are shared with 8.9 and built here once (Boss, 1:33 PM EDT): the unit finalizer (ADM.82, also API.8) and the recipe schema (API.28, also API.18).
Needs: The plan model, aspects and regions come first (phase 1); then the approval screen, locks, regen screens and density overrides (phase 2); then the editors that make plans and the rearranger (phase 3); then regeneration to a specification (phase 4); then the end-to-end tests, the help and the headline (phase 5). ADM.36 (edit orbits) comes from the orbital update foundations in 8.3 and ADM.52 (the planet editor) from earlier. The override journal has to keep regeneration, parity repair and the golden-seed tests (GEN.57, DB.9, DB.17 in 8.1 phase 1) working.
Suggested:
- S8.6-a: Batch edits: apply one plan to a list or a region of systems, with one approval for the whole batch.
- S8.6-b: Two candidate plans for one request, side by side, and the admin picks one.
- S8.6-c: Named presets for the regen screens (a hot Jupiter system, a Sol-like system, a binary with a circumbinary planet).
- S8.6-d: Locked marks drawn on the maps so the admin sees what a regeneration will not touch.
- S8.6-e: Export and import an edited area as one file, to move a hand-made region between servers.

## Release 8.7: The pilot's console and built-in help
Status: proposed
Headline: A built-in help section and the pilot's starmap look: every page explains itself and the interface reads like a navigation console
Today: (new)
Change: The help pages (DOC.6 to DOC.16) and the console design (UX.43) have no headline today.
Needs: The help pages describe what 8.1 to 8.6 built, so this version goes after them; DOC.6 (the template and coverage test) first.
Suggested:
- S8.7-a: Context help on every control: a "?" popover linked to the matching help page.
- S8.7-b: A command palette: type a name, an ID or a coordinate to go there.
- S8.7-c: A high-contrast "cockpit" theme next to light and dark, remembered in the account.
- S8.7-d: A keyboard shortcut sheet, opened with "?".
- S8.7-e: Site-wide search over objects, help pages and docs.
- S8.7-f: Printable system sheets.

## Release 8.8: The deep sky
Status: proposed
Headline: Nebulae, clusters, satellites and the rarest stars drawn on the map and reachable by flight
Today: (new)
Change: The nebula drawing items (MAP.176 to MAP.183), the cluster work and the "later" stellar phenomena have no headline today.
Needs: The scatters of 8.1 (components, nebula stars, Wolf-Rayet, scattered nebulae in the index) and the nebula tiles and cover (MAP.173, MAP.174). MAP.194 (nebulae and clusters open their own zoom window) needs MAP.176.
Suggested:
- S8.8-a: A curated notable-objects catalogue (famous remnants, clusters, nebulae) with names, descriptions and a flight to each.
- S8.8-b: Tidal streams and the Magellanic Stream drawn as paths of stars and gas.
- S8.8-c: Fly inside a nebula and see its dust dim the stars behind it (needs the dust model in VIEW.8).
- S8.8-d: A variable-star sheet with light curves for Mira stars, Cepheids and eclipsing binaries, from the interaction states of GEN.270.
- S8.8-e: A rare-objects search: every X-ray binary, Wolf-Rayet star or magnetar within a radius of a point.

## Release 8.9: Generate anywhere
Status: proposed
Headline: Remote generation: build a galaxy on a local machine, upload it through the API and have the server verify it
Today: (new)
Change: The remote generation chain (API.3 and its parts) and recipes (API.18, API.19) have no headline today. Boss (1:33 PM EDT): 8.6 and this version are closely related, so what both need is built once in 8.6: the unit finalizer (ADM.82, used by API.8) and the recipe schema and runner (API.28, used by API.18); API.8 and API.18 here keep only the upload and route parts, and API.19 uses the region object of 8.6 (GEN.302). It depends only on the reproducibility items of 8.1 phase 1 and those 8.6 pieces, so it can move earlier if Boss wants it sooner.
Needs: GEN.57, GEN.135, OPS.28 and TEST.77 (reproducibility, 8.1 phase 1) first: remote generation must reproduce what the server makes (API.17). API.10, API.11 and API.12 are the base; API.19 needs API.18.
Suggested:
- S8.9-a: A pip and Docker package of the generator with a one-command local run.
- S8.9-b: Client libraries for Python and JavaScript generated from the API description.
- S8.9-c: Webhooks when a job finishes or an upload is verified.
- S8.9-d: Per-key rate limits and a usage page, from the API call log (API.15).
- S8.9-e: A gallery of shared recipes and notebook examples.

## Release 9.0: The view from a planet
Status: proposed; Boss does not want the gravity map here
Headline: Stand on a world and look up: the sky from any planet or moon, with real constellations, rendered as a picture
Today: 9.0 the gravity map (MAP.167)
Change: Boss (2:00 PM EDT) targeted the backups and repair chain for 9.0 (umbrella ADM.91, 10 items in phases 3 and 4): the control database and parity parts of the backup archive (ADM.88, DB.40), restoring each (ADM.89, DB.41), repairing the galaxy from an archive (DB.42), one repair page (ADM.90), a nightly full backup from the daily maintenance run (OPS.55, Boss 2:06 PM EDT: yes), the end-to-end test and the guide (TEST.180, DOC.56). The settings page and the settings backups stay in 8.1. The gravity map moved to 9.1. A 9.0 headline still needs Boss's choice; this page recommends the view from a planet. 9.0 is a `major` note and needs his word.
Needs: 8.1 (the camera that reaches a world; the settings page and the archive of ADM.87; the parity file and repair of DB.9; DB.17), 8.3 (so the sky changes with the date) and 8.4 (the worlds). The backups chain needs only 8.1. VIEW.1 and VIEW.2 start with a research session with Boss; VIEW.6 settles the handedness of the galaxy before real sky coordinates.
Alternative (Politics and territories): Territories that grow and shift on the map, with borders drawn and polities that trade. Placed at 9.2; it could be the 9.0 headline.
Alternative (Generate anywhere): Remote generation as a platform headline (API.3). Placed at 8.9; it could be the 9.0 headline instead.
Alternative (Edit everything): Already the 8.6 headline; it could be held back to a major release if Boss prefers.
Suggested:
- S9.0-a: Sky effects from the surface: atmospheric colour and twilight, moons with phases, rings in the sky and the planets that show to the eye.
- S9.0-b: A sky-at-a-date slider (needs 8.3): the sky as it was or will be, with eclipse and conjunction predictions.
- S9.0-c: Sunrise and sunset times and a day-and-night cycle for any world.
- S9.0-d: A planetarium mode for a big screen, and a PNG or poster export of any sky.
- S9.0-e: Constellation names per species or polity, from the name generator (VIEW.4).

## Release 9.1: The gravity map
Status: proposed; moved from 9.0
Headline: A heat map of the gravitational field inside a sector and across a star system, as pull, well depth and tidal strength
Today: 9.1 the gravity map (MAP.167)
Change: Moved from 9.0 (Boss does not think it the best major headline). Its foundations (GEN.198, TEST.131) are already in 8.3; the rest stays together.
Needs: GEN.198 (8.3), then GEN.199 and PERF.72 (the per-zone grid and its cost), API.24, MAP.168, MAP.169 and UX.94. MAP.170 (the coarse Galaxy Map layer) stays later and optional.
Suggested:
- S9.1-a: A Lagrange point finder: a list and map markers for every stable point in a system.
- S9.1-b: Gravity-assist hints in the courses of 8.5: a flyby that saves fuel.
- S9.1-c: A stability view for planets in multi-star systems: where an orbit survives.
- S9.1-d: Roche limit and tidal-lock zones shown around bodies.

## Release 9.2: Politics and territories
Status: Boss ordered it before ships and stations; idea only, no items
Headline: Territories that grow and shift over time on the map, with borders drawn and polities that trade
Today: (new)
Change: Boss, 1:07 PM EDT: politics and territories come before ships and stations. No items exist yet beyond the population-and-politics design note; the item breakdown is for a later change after Boss's answers.
Needs: 8.3 (time) and the population items of 8.4 (species and limits). The population-and-politics design note holds the groundwork.
Suggested:
- S9.2-a: Border drawing on the Galaxy and Sector Maps, with a layer toggle.
- S9.2-b: Territory growth over the orbital update's days, with a change log.
- S9.2-c: Trade routes between polities as courses.
- S9.2-d: Population change per system over time.

## Release 9.3: Ships and stations
Status: proposed; after politics and territories
Headline: Users place their own ships and stations and follow them across the live map
Today: (new)
Change: USR.10 and POP.11 have no headline today; they follow politics and territories (Boss, 1:07 PM EDT).
Needs: 8.2 (accounts), 8.3 (movement), 8.5 (courses to fly) and 9.2 (territories, so stations belong to a polity).
Suggested:
- S9.3-a: Fleets: group ships, give the group a course, and move it together.
- S9.3-b: A game master role between player and admin: places objects in a polity's territory (built on the permission catalogue of 10.0 later, or its own role until then).
- S9.3-c: Fog of war: each user sees only what they have charted or their ships have visited.
- S9.3-d: Station docking and a ship log that writes itself from the course.

## Release 9.4: Galaxy types
Status: Boss named it (8:48 PM EDT); items filed
Headline: Galaxy types: the seed draws a galaxy of a known type (spiral, lenticular, elliptical, irregular or dwarf), every shape number calculated from published scaling relations
Today: (new)
Change: New headline (Boss, 8:48 PM EDT: "different galactic models for generation based on known galaxy types", broken into a chain of required items and counted back from 9.4), filed as 18 items in 6 phases with the plan in docs/design/galaxy-types.md. A galaxy type is a parameter family, never a copy of a named galaxy (no fixed content, Boss 4:11 and 4:14 PM EDT): the seed draws the stellar mass, then the family from the mass-dependent type fractions, then each number from a scaling relation. It reuses the shape, the density model, the potential, the halo, the satellites and the scatters. The old "Hubble-sequence generator" suggestion (SL-b) is this headline. 9.4 is Boss's own number and needs no .major note.
Needs: three steps of the chain come earlier so the families only fill in the shape: GEN.378 (the checked design and its scaling relation table, after GEN.368) in 8.1, then DB.66 (the family key and component columns, a migration) and GEN.379 (the density model as a list of components, which the potential GEN.338 reads) in 8.3. Then, in 9.4: the type draw and the Sérsic spheroid (phase 1); the spiral, elliptical and irregular and dwarf families (phase 2); the lenticular family, populations, clusters with the nucleus (phase 3); creation (phase 4); the Generate page, the Galaxy Map, the plausibility tests and the docs (phase 5; each family also gets a simple vector outline, UX.135, in phase 1); the headline (phase 6). It also reads GEN.368, GEN.369, GEN.297, GEN.134 and GEN.319.
Lanes: not assigned; GEN.378 is proposed for Core 1 (it writes GEN.368's table), and the rest sits after 8.3.
Suggested:
- S9.4-a: Interacting and merging galaxies: a tidal tail or a merger remnant as a family variant.
- S9.4-b: Satellite galaxies drawn from the same families as the main galaxy, so a satellite's stars follow its family's laws.
- S9.4-c: Active nuclei by family: a quasar or jet drawn from the black hole mass and the gas supply.
- S9.4-d: A family comparison page: draw many seeds and plot each family's mass-size and mass-metallicity points against the published relations.

## Release 10.0: Roles and permissions
Status: Boss named it (1:43 PM EDT); items filed
Headline: Roles and permissions: users submit edit requests, moderators approve user requests and admins approve moderator requests
Today: (new)
Change: New headline (USR.25), filed as 18 items in 4 phases with the design in docs/design/roles-and-requests.md. Three roles, User, Moderator and Admin, with the Owner above admin. A user proposes an edit, a moderator approves it; a moderator proposes an edit, an admin approves it; an admin approves their own plan at the time, so an admin confirms with the preview and does not request (Boss, 1:49 PM EDT). A request is a change plan of 8.6, so approval, preview, history, undo and audit are written once. 10.0 is a major version; the `.major.md` note waits for Boss's word.
Needs: 8.2 (accounts: USR.2 holds the role column with the four values and one guard, USR.3 mail, USR.7 the user interface) and 8.6 (change plans: ADM.58 engine with the `requested` state, ADM.59 store, ADM.60 apply, ADM.66 approval screen, ADM.68 history). Then the permission catalogue and the guard sweep, the request tables (phase 1); submit, queue and decide (phase 2); the approver screen, comments, lifecycle, limits, notices and the audit view (phase 3); the end-to-end and matrix tests, the help page and the headline (phase 4).
Lanes: not assigned; the work sits after 8.6, so no lane queue holds it.
Suggested:
- S10.0-a: Two-person approval for risky plans: a density or regeneration plan above a size needs a second admin.
- S10.0-b: Custom roles: an admin builds a role from the permission catalogue (a "cartographer" who may edit names only).
- S10.0-c: Auto-approve rules: trusted users' name and description requests apply without a moderator after N approved requests.
- S10.0-d: A request digest email: one daily mail to each approver instead of one per request.
- S10.0-e: Request templates: "rename this system", "report a wrong class", each a ready-made plan the user fills in.
- S10.0-f: Moderator tools beyond approvals: lock a name, hide an object description pending review, a moderation log page.

## Release later: Future features
Status: possible future features, not scheduled
Headline: No headline: features that belong to no version yet, noted for later
Today: (new)
Change: Boss, 1:07 PM EDT: multiple galaxies is a possible future feature, not a headline (Boss, 1:43 PM EDT: 10.0 is roles and permissions). The items (GEN.9, GEN.157, GEN.164, MAP.145, NAV.57) sit in the "later" list.
Needs: GEN.9 (the plan for more than one galaxy in the database) first; GEN.157 (a verified table of about 25 real galaxies) before MAP.145 and NAV.57; GEN.164 needs GEN.9. The infinite zoom of 8.1 would extend past the galaxy edge.
Suggested:
- SL-a: Intergalactic flight: the zoom of 8.1 continues out to a galaxy group, with the cosmic web as a faint background.
- SL-c: A galaxy picker: several generated galaxies in one database, each with its seed, and a way to compare them.
- SL-d: A long-timescale merger preview such as the Milky Way and Andromeda, run as a stated-span orbital update.
- SL-e: Lookback time and travel-time readouts for objects in other galaxies.

## Questions for Boss
- Decided (Boss, 1:29 PM EDT, "defaults are acceptable"): keep the rule that 8.1 stamps only when all of its phases are done; ADM.32, ADM.35, ADM.55, ADM.56, MAP.119 and MAP.120 sit in 8.6 and UX.46 and UX.47 in 8.2; 8.3 to 8.9 as proposed; 9.0 is the view from a planet; multiple galaxies stays under "later".
- 8.6 undo: an applied plan can be undone while no newer plan touched the same objects. Default: yes.
- Decided (Boss, 1:43 PM and 1:49 PM EDT): an admin needs no approval; an admin confirms their own plan after the preview and there is no second approver. Any admin may also confirm a plan another admin proposed. Two-person approval is only a suggestion (S10.0-a).
- 8.6 density plans: the planner may remove generated, unnamed, unpopulated systems inside an approved plan, listed with counts per sector. Default: yes.
- 8.6 specifications: forced values the generator would not draw are allowed and marked on the plan as forced. Default: yes.
- 10.0 what users may request. Default: a user may request a rename, a description or a note on one object, and no structural change; a moderator may request those plus a move, an orbit edit or a regeneration of one object; an admin may do everything. Boss may widen the user list.
- 10.0 what a moderator may approve. Default: only user requests inside the user's list; a request that touches a locked aspect, a named object or a populated world goes to an admin whoever asked.
- 10.0 request limits and expiry. Default: 5 open requests and 20 submissions a day per user; a pending request expires after 30 days; both are settings on the settings page.
- 10.0 notices. Default: approvers and requesters are emailed, each account can turn the emails off.
- Which suggestions to file as items? Name the codes, for example "file S8.1-a and S8.2-c".
