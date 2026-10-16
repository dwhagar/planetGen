# Release roadmap, v8.0 to v10.0 (proposal)

PROPOSAL, asked for by Boss on 2026-10-16 (12:45 PM EDT). Nothing in the repository changes because of this page: no `Headline:` line in `docs/TODO.md` moves, no `.major.md` note is added and none of the suggestions is filed as an item until Boss approves. Where Boss named a headline, it is used as given.

Thirteen releases from the current 8.0 line to 10.0. Each has one headline, the open items it holds, what it needs, how much of it is built and a list of ideas that are not filed yet. Boss's picks are 8.1 infinite zoom, 8.2 accounts with email invitations and 8.5 courses that bend around gravity wells; the 9.0 headline is open.

## How to read this

- A point release is a REVISION bump (a `minor` note); 9.0 and 10.0 are `major` notes, which need Boss's word. The stamp Action numbers each merge as MAJOR.REVISION.BUILD, so a "release" here is the headline a REVISION or MAJOR earns, and the 8.0.15xx builds are the revisions inside the current line.
- While `REVISION_HOLD` is on nothing stamps 8.1 until Phase 1 is complete (Boss, 2026-10-09). 8.0 therefore carries everything that merges until then: the whole Phase 1 tail and the Phase 2 plan stages the lanes are on.
- Each open item is placed once, by the feature it serves, not by its phase. An item sits in the same release as, or an earlier release than, everything that needs it; the page flags any item that breaks this, and there are none in this proposal.
- Progress share is finished over finished plus open. Finished counts the items the same ledger threads built before this roadmap (a fixed number, counted on 2026-10-16) plus every listed item that retires afterwards. Phase tiers keep their own bars on the tier plan; a release that draws on several phases shows its mix.
- "Suggested, not filed" lists are ideas, each with a code (S8.1-a and so on) so Boss can say "file S8.1-a and S8.2-c". Each one gets checked against the open items first, so nothing is filed twice.

## Release 8.0: The whole-sphere galaxy
Status: current line
Headline: The whole-sphere galaxy: spherical sectors and 96-bit IDs, a reach of 262 kpc with halo, satellites and a realistic core, and every kind of star and phenomenon represented
Today: 8.0 the orbital update system (GEN.105)
Change: Moves the orbital update to 8.3. GEN.105 needs eight open items (GEN.115, GEN.142, GEN.289, GEN.294, GEN.109, GEN.110, GEN.290, GEN.291), two of which take database migrations from Core 1, and the lanes are on the whole-sphere stages first. So 8.0 closes with what the lanes are building in this period, and the orbital update earns its own release.
Needs: nothing outside itself. This line closes when Phase 1 is complete (the hold), which also stamps 8.1.
Lanes: Core 1, Core 2 and Core 3 hold most of it: GEN.237 to GEN.243, GEN.264 to GEN.275, GEN.280 to GEN.285, PERF.97 to PERF.100, PERF.105 and UX.118 (Core 2, after MAP.146).
Items (Phase 1 tail, which gates the 8.1 stamp): ADM.32, ADM.50, ADM.55, ADM.56, API.4, API.7,
  API.15, API.26, DB.9, DB.16, DB.17, DB.29, DB.30, DB.31, DB.33, DOC.18, DOC.19, DOC.20, DOC.24,
  DOC.45, DOC.47, DOC.52, DOC.53, GEN.24, GEN.57, GEN.99, GEN.102, GEN.103, GEN.135, GEN.136,
  GEN.179, GEN.180, GEN.181, GEN.202, GEN.218, GEN.219, GEN.220, GEN.221, GEN.229, GEN.231, MAP.172,
  MAP.173, MAP.174, OPS.28, OPS.29, OPS.35, OPS.37, PERF.29, PERF.30, PERF.41, PERF.47, PERF.75,
  PERF.88, PERF.91, PERF.97, PERF.98, PERF.99, PERF.100, PERF.105, POP.14, TEST.77, TEST.111,
  TEST.116, TEST.126, TEST.127, TEST.128, TEST.129, TEST.130, TEST.134, TEST.135, TEST.137,
  TEST.138, TEST.140, TEST.141, TEST.144, TEST.148, TEST.155, TEST.156, TEST.167, TEST.170,
  TEST.171, TEST.172, TEST.173, TEST.174, UX.3, UX.46, UX.47, UX.49, UX.105, UX.107, UX.108, UX.109,
  UX.110, UX.118
Items (Whole-sphere galaxy and stars, Phase 2 plan stages): API.27, GEN.55, GEN.129, GEN.130,
  GEN.134, GEN.161, GEN.162, GEN.163, GEN.210, GEN.237, GEN.238, GEN.240, GEN.241, GEN.242, GEN.243,
  GEN.244, GEN.252, GEN.258, GEN.259, GEN.263, GEN.264, GEN.267, GEN.269, GEN.270, GEN.271, GEN.272,
  GEN.273, GEN.274, GEN.275, GEN.276, GEN.280, GEN.282, GEN.283, GEN.284, GEN.285, GEN.286, GEN.287,
  GEN.298, NAV.45, NAV.60, NAV.61, UX.102, UX.113, UX.114, UX.115, UX.116
Items (Generation speed and the queue, Phase 2): DB.18, GEN.40, GEN.42, OPS.15, PERF.18, PERF.20,
  PERF.35, PERF.36, PERF.48, PERF.86
Items (Tests and docs that follow the whole sphere, Phase 2): DOC.4, DOC.48, DOC.49, DOC.50, DOC.51,
  TEST.159, TEST.160, TEST.161, TEST.162, TEST.164, TEST.165, TEST.166, TEST.168, TEST.169
Built before: 590 (every finished Phase 0 and Phase 1 item)
Suggested:
- S8.0-a: A census page for the whole sphere: counts of every star type and phenomenon by region, from the representation audit (TEST.160), public and cached.
- S8.0-b: Galaxy export and import: one file holds a galaxy (seed, settings, stored sectors) to back it up or move it to another server, built on the parity file (DB.9).
- S8.0-c: Named generation presets on the Generate page (Milky Way, small dwarf, globular-rich) saved as recipes, a step toward API.18.

## Release 8.1: Infinite zoom
Status: Boss named it; next stamp
Headline: Infinite zoom: one camera flies from the whole galaxy to a star system and on into its planets and moons
Today: 8.1 infinite zoom (the fly-through Galaxy Map), MAP.146
Change: No change to the headline (Boss, 2026-10-16: 8.1 gets infinite zoom). MAP.152 and MAP.146 are Core 2's next two items, in that order (Boss, 12:23 PM EDT). The existing items cover the galaxy to star system flight and the 3D system scene (MAP.125, built) covers planets and moons; no open item names an end-to-end check that one camera flies from the galaxy to a moon (checked against the open items on 2026-10-16), which is suggested as S8.1-a.
Needs: The 8.0 hold: 8.1 stamps when Phase 1 is complete, so the Phase 1 tail above must finish too. MAP.152 first, then MAP.146 (MAP.189 and MAP.171 need MAP.152). MAP.122 (Select mode) feeds NAV.49 in 8.5.
Items: ADM.35, MAP.59, MAP.75, MAP.119, MAP.120, MAP.121, MAP.122, MAP.139, MAP.140, MAP.141,
  MAP.146, MAP.147, MAP.152, MAP.156, MAP.160, MAP.171, MAP.189, MAP.190, MAP.191, MAP.195, MAP.209,
  MAP.210
Suggested:
- S8.1-a: An acceptance flight and frame budget: a scripted camera flight from the whole galaxy to a named moon (galaxy, sector, system, planet, moon) with a frame-rate floor, run in the browser tests. No open item covers the last two steps as one flight.
- S8.1-b: A scale bar and distance readout that change units with the zoom (light years, parsecs, AU, kilometres) from the shared ladders.
- S8.1-c: Names fade in with the zoom: labels for the largest nearby things only, so the view stays readable at every scale.
- S8.1-d: Shareable camera positions: a link that flies a visitor to the same place, built on object IDs and the bookmark format.
- S8.1-e: A guided tour: a scripted flight past notable places (Sol, the galactic core, a nebula, a globular cluster) a visitor can replay.
- S8.1-f: Touch and keyboard flight controls and a reduced-motion setting.

## Release 8.2: Accounts and email invitations
Status: Boss named it; new release
Headline: User accounts, with sign-up by an email invitation that only an admin can send
Today: (new)
Change: New release (Boss, 2026-10-16). It is already planned as USR.1 to USR.7, in Phase 3+ today; this moves it up. Nothing needs filing: USR.1 holds the accounts, USR.3 the SMTP settings and USR.4 the invitations. One change to USR.4 on approval: its decided answer says the invitee's email is optional and the admin may only copy the link; Boss's 12:45 PM EDT message makes the email invitation the main path, with the link shown only when SMTP is not set up.
Needs: ADM.43 (the full admin configuration page) is where the SMTP form sits, so it comes first; SEC.32 (Argon2id) before USR.2 so new accounts start on it. ADM.55 and ADM.56 (admin streamline, in the 8.0 tail) make the admin invitation pages easier to build. Mail goes out as a queue job (USR.3), so the queue steps PERF.97 to PERF.100 should be in first.
Items: ADM.43, ADM.44, API.6, SEC.32, USR.1, USR.2, USR.3, USR.4, USR.5, USR.6, USR.7, USR.8, USR.9,
  USR.11, UX.45, UX.48, UX.92
Built before: 5 (the Bookmarks thread of Phase 1)
Suggested:
- S8.2-a: Two-factor sign-in (TOTP) for every account, required for admins; the Owner transfer in USR.6 already assumes it.
- S8.2-b: An invitation the admin can write: an optional personal message, resend, revoke, and a notice back to the admin when it is accepted.
- S8.2-c: A sign-in activity page: devices, last seen, and sign out everywhere.
- S8.2-d: An admin user list: search, disable, last sign-in and role history.
- S8.2-e: Bounce handling: mark an address undeliverable and tell the admin who invited it.
- S8.2-f: Opt-in emails for things a user follows, such as a Hill-radius overlap in a bookmarked system (needs 8.3's GEN.111).

## Release 8.3: The moving galaxy
Status: proposed; moved from 8.0
Headline: The orbital update: every object moves on its orbit each day, with n-body influence, collisions and ejections
Today: 8.0 the orbital update system (GEN.105)
Change: Moved here from 8.0. GEN.105 is the headline item. GEN.198 (the gravity field evaluator) is placed here, not in 9.0, because the courses of 8.5 reuse it, so NAV.59's note "the pull model of the gravity map (9.0)" becomes "of the evaluator in 8.3".
Needs: GEN.294 and GEN.289 (every object as one SpatialPosition3D and kept planets and moons) take migrations from Core 1 first; GEN.115 and GEN.142 have no prerequisites. GEN.111 needs 8.2's SMTP settings (USR.3). GEN.299 and GEN.300 (the Quantity class) need GEN.294 and GEN.135. OPS.16 and OPS.17 set up the nightly run.
Items: ADM.36, GEN.105, GEN.109, GEN.110, GEN.111, GEN.115, GEN.139, GEN.140, GEN.141, GEN.142,
  GEN.143, GEN.144, GEN.197, GEN.198, GEN.212, GEN.289, GEN.290, GEN.291, GEN.292, GEN.293, GEN.294,
  GEN.295, GEN.299, GEN.300, MAP.192, MAP.193, OPS.16, OPS.17, OPS.30, TEST.131
Built before: 35 (the Orbital updates and Daily maintenance threads)
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
Change: New release. The planet class work (GEN.33, GEN.28, GEN.27, GEN.90 and the rest) is Phase 2 today with no headline.
Needs: The habitability index (built). GEN.33 first, one class per PR. GEN.95 comes before GEN.93 (nebula conditions on planets). POP.13 (the species editor) goes with it because the "find me a world" idea needs species limits.
Items: ADM.52, GEN.27, GEN.28, GEN.29, GEN.33, GEN.90, GEN.91, GEN.92, GEN.93, GEN.95, GEN.112,
  GEN.145, GEN.146, MAP.143, POP.13, UX.32, UX.100, UX.101
Built before: 15 (the Habitability, Nebula planets and Classes threads)
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
Items: NAV.4, NAV.5, NAV.6, NAV.17, NAV.18, NAV.19, NAV.21, NAV.22, NAV.23, NAV.25, NAV.26, NAV.27,
  NAV.28, NAV.36, NAV.39, NAV.49, NAV.51, NAV.54, NAV.55, NAV.56, NAV.58, NAV.59
Built before: 16 (the Routing, Courses and Nearby search threads)
Suggested:
- S8.5-a: Ship profiles (speed, fuel range, shielding) so a course shows its fuel and time for a chosen ship.
- S8.5-b: Hazard costs: radiation (GEN.256) and nebula density as route penalties, with a safest or shortest toggle.
- S8.5-c: Read-only course sharing by link, with no account needed.
- S8.5-d: Animated travel: play a course as a ship flying along it at a chosen speed, with the bends around the wells visible.
- S8.5-e: Suggested stops at systems with the facilities a ship needs (uses the tech level and facility data).
- S8.5-f: Compare courses: straight, bent and safest, side by side.

## Release 8.6: The pilot's console and built-in help
Status: proposed
Headline: A built-in help section and the pilot's starmap look: every page explains itself and the interface reads like a navigation console
Today: (new)
Change: New release. The help pages (DOC.6 to DOC.16) and the console design (UX.43) have no headline today.
Needs: The help pages describe what 8.1 to 8.5 built, so this release goes after them; DOC.6 (the template and coverage test) first.
Items: DOC.6, DOC.7, DOC.8, DOC.9, DOC.10, DOC.11, DOC.12, DOC.13, DOC.14, DOC.15, DOC.16, UX.43,
  UX.82
Suggested:
- S8.6-a: Context help on every control: a "?" popover linked to the matching help page.
- S8.6-b: A command palette: type a name, an ID or a coordinate to go there.
- S8.6-c: A high-contrast "cockpit" theme next to light and dark, remembered in the account.
- S8.6-d: A keyboard shortcut sheet, opened with "?".
- S8.6-e: Site-wide search over objects, help pages and docs.
- S8.6-f: Printable system sheets.

## Release 8.7: The deep sky
Status: proposed
Headline: Nebulae, clusters, streams and the rarest stars drawn on the map and reachable by flight
Today: (new)
Change: New release. The nebula drawing items (MAP.176 to MAP.183), the cluster work and the "later" stellar phenomena have no headline today.
Needs: MAP.173 and MAP.174 (nebula tiles and cover) and GEN.202 (scattered nebulae in the index), all in the 8.0 tail, and GEN.237, GEN.241 and GEN.243 (components, nebula stars, Wolf-Rayet) in 8.0. MAP.194 (nebulae and clusters open their own zoom window) needs MAP.176.
Items: GEN.209, GEN.239, GEN.245, GEN.246, GEN.247, GEN.248, GEN.249, GEN.253, GEN.255, GEN.256,
  GEN.257, GEN.296, GEN.297, MAP.132, MAP.142, MAP.176, MAP.177, MAP.178, MAP.179, MAP.180, MAP.181,
  MAP.183, MAP.194, MAP.207
Suggested:
- S8.7-a: A curated notable-objects catalogue (famous remnants, clusters, nebulae) with names, descriptions and a flight to each.
- S8.7-b: Tidal streams and the Magellanic Stream drawn as paths of stars and gas.
- S8.7-c: Fly inside a nebula and see its dust dim the stars behind it (needs the dust model in VIEW.8).
- S8.7-d: A variable-star sheet with light curves for Mira stars, Cepheids and eclipsing binaries, from the interaction states of GEN.270.
- S8.7-e: A rare-objects search: every X-ray binary, Wolf-Rayet star or magnetar within a radius of a point.

## Release 8.8: Generate anywhere
Status: proposed
Headline: Remote generation: build a galaxy on a local machine, upload it through the API and have the server verify it
Today: (new)
Change: New release. The remote generation chain (API.3 and its parts) and recipes (API.18, API.19) are Phase 3 today with no headline. It could move ahead of 8.6 and 8.7 if Boss wants it sooner, since it depends only on the 8.0 determinism items.
Needs: GEN.57, GEN.135, OPS.28 and TEST.77 (reproducibility, 8.0) first: remote generation must reproduce what the server makes (API.17). API.10, API.11 and API.12 are the base; API.19 needs API.18.
Items: ADM.13, API.3, API.5, API.8, API.10, API.11, API.12, API.13, API.14, API.17, API.18, API.19
Built before: 6 (the API threads)
Suggested:
- S8.8-a: A pip and Docker package of the generator with a one-command local run.
- S8.8-b: Client libraries for Python and JavaScript generated from the API description.
- S8.8-c: Webhooks when a job finishes or an upload is verified.
- S8.8-d: Per-key rate limits and a usage page, from the API call log (API.15).
- S8.8-e: A gallery of shared recipes and notebook examples.

## Release 9.0: The view from a planet
Status: proposed; Boss does not want the gravity map here
Headline: Stand on a world and look up: the sky from any planet or moon, with real constellations, rendered as a picture
Today: 9.0 the gravity map (MAP.167)
Change: The gravity map moves to 9.1. A new 9.0 headline needs Boss's choice; the page recommends the view from a planet, listed first below. Whichever is chosen, 9.0 is a `major` note and needs his word.
Needs: 8.1 (the camera that reaches a world), 8.3 (so the sky changes with the date) and 8.4 (the worlds). VIEW.1 and VIEW.2 start with a research session with Boss; VIEW.6 settles the handedness of the galaxy before real sky coordinates.
Alternative (Ships and stations on the live map): Users place their own ships and stations and track them moving across the map (USR.10, POP.11), the strongest hook for the game. Needs 8.2 accounts, 8.3 movement and 8.5 courses, so it follows them; shown as 9.2.
Alternative (Generate anywhere): Remote generation as a platform headline (API.3). It is placed at 8.8; it could be the 9.0 headline instead, with 8.8 dropped.
Alternative (Many galaxies): Neighbour galaxies and the intergalactic frame (GEN.9, GEN.157). It is placed at 10.0 as the end of this roadmap.
Alternative (The living universe): The events feed and time control of 8.3 as a major release. It is already inside 8.3, so it is only a candidate if 8.3 is split.
Items: VIEW.1, VIEW.2, VIEW.3, VIEW.4, VIEW.6, VIEW.7, VIEW.8, VIEW.9, VIEW.10
Built before: 0
Suggested:
- S9.0-a: Sky effects from the surface: atmospheric colour and twilight, moons with phases, rings in the sky and the planets that show to the eye.
- S9.0-b: A sky-at-a-date slider (needs 8.3): the sky as it was or will be, with eclipse and conjunction predictions.
- S9.0-c: Sunrise and sunset times and a day-and-night cycle for any world.
- S9.0-d: A planetarium mode for a big screen, and a PNG or poster export of any sky.
- S9.0-e: Constellation names per species or polity, from the name generator (VIEW.4).

## Release 9.1: The gravity map
Status: proposed; moved from 9.0
Headline: A heat map of the gravitational field inside a sector and across a star system, as pull, well depth and tidal strength
Today: 9.0 the gravity map (MAP.167)
Change: Moved here from 9.0 (Boss does not think it the best major headline). Its foundations (GEN.198, TEST.131) are already in 8.3; the rest stays together.
Needs: GEN.198 (8.3), then GEN.199 and PERF.72 (the per-zone grid and its cost), API.24, MAP.168, MAP.169 and UX.94. MAP.170 (the coarse Galaxy Map layer) stays later and optional.
Items: API.24, GEN.199, MAP.167, MAP.168, MAP.169, MAP.170, PERF.72, UX.94
Built before: 0
Suggested:
- S9.1-a: A Lagrange point finder: a list and map markers for every stable point in a system.
- S9.1-b: Gravity-assist hints in the courses of 8.5: a flyby that saves fuel.
- S9.1-c: A stability view for planets in multi-star systems: where an orbit survives.
- S9.1-d: Roche limit and tidal-lock zones shown around bodies.

## Release 9.2: Ships and stations
Status: proposed
Headline: Users place their own ships and stations and follow them across the live map
Today: (new)
Change: New release. USR.10 and POP.11 are Phase 3+ with no headline.
Needs: 8.2 (accounts), 8.3 (movement) and 8.5 (courses to fly).
Items: POP.11, USR.10
Built before: 0
Suggested:
- S9.2-a: Fleets: group ships, give the group a course, and move it together.
- S9.2-b: A game master role between player and admin: places objects in a polity's territory.
- S9.2-c: Fog of war: each user sees only what they have charted or their ships have visited.
- S9.2-d: Station docking and a ship log that writes itself from the course.

## Release 9.3: Politics and territories
Status: idea only; no items
Headline: Territories that grow and shift over time on the map, with borders drawn and polities that trade
Today: (new)
Change: New idea; no items exist yet beyond the population design note.
Needs: 8.3 (time) and 9.2 (ships) if trade routes use them. The population-and-politics design note holds the groundwork.
Built before: 0
Suggested:
- S9.3-a: Border drawing on the Galaxy and Sector Maps, with a layer toggle.
- S9.3-b: Territory growth over the orbital update's days, with a change log.
- S9.3-c: Trade routes between polities as courses.
- S9.3-d: Population change per system over time.

## Release 10.0: Many galaxies
Status: proposed; the end of the roadmap
Headline: A universe, not one galaxy: neighbour galaxies with their own seeds, and flight between them
Today: (new)
Change: New release, and the next `major` after 9.0. GEN.9, GEN.157 and the rest are Phase 3+ today.
Needs: GEN.9 (the plan for more than one galaxy in the database) first; GEN.157 (a verified table of about 25 real galaxies) before MAP.145 and NAV.57; GEN.164 needs GEN.9. The infinite zoom of 8.1 is extended past the galaxy edge.
Items: GEN.9, GEN.157, GEN.164, MAP.145, NAV.57
Built before: 0
Suggested:
- S10.0-a: Intergalactic flight: the zoom of 8.1 continues out to a galaxy group, with the cosmic web as a faint background.
- S10.0-b: A Hubble-sequence generator: elliptical, spiral and irregular galaxies from one set of rules.
- S10.0-c: A galaxy picker: several generated galaxies in one database, each with its seed, and a way to compare them.
- S10.0-d: A long-timescale merger preview such as the Milky Way and Andromeda, run as a stated-span orbital update.
- S10.0-e: Lookback time and travel-time readouts for objects in other galaxies.

## Questions for Boss
- 8.0 headline: move the orbital update to 8.3 and call the current line the whole-sphere galaxy? Default if no answer: yes.
- 8.1 stamp: keep the rule that 8.1 stamps only when all of Phase 1 is complete, or stamp it when the zoom items are done and let the Phase 1 tail ride into 8.2? Default: keep the rule.
- 8.2: make the admin's email invitation the main path (the link is shown only when SMTP is not set up) and change USR.4 to say so? Default: yes.
- 9.0 headline: the view from a planet (recommended), ships and stations, generate anywhere, or another? Default: the view from a planet.
- 9.1: the gravity map as a point release after 9.0, with its foundation GEN.198 built in 8.3 for the courses? Default: yes.
- Which suggestions to file as items? Name the codes, for example "file S8.1-a and S8.2-c".
