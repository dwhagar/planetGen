# Changelog

## [8.0.911] - 2026-10-10

### Added
- Generate has one "Redo scatters" box in place of "Rebuild the bright stars" (GEN.196). Tick the scatters to redo (the massive stars, the bright stars, the phenomena), set each one's new limit, and they run as one job with every stage listed and the unticked ones shown as skipped. A scatter clears and rewrites only its own rows, and the central black hole or quasar is always rebuilt. `planetgen plan --redo-scatters mass luminosity phenomena` does the same from the command line.
- Random neighborhoods have a "Keep away from filled space" option (`planetgen galaxy --avoid-filled-space`, random-start mode, and a checkbox on the Generate page). A start whose neighborhood touches a sector that is already generated is dropped, so each neighborhood grows in empty space, and the run says how many candidate centres qualified. (GEN.186)
- **A sector that nothing was generated in opens if a scatter placed something there, marked uncharted (MAP.162).** On the Galaxy Map, clicking a sector that is not generated yet, or opening its address (`?sector=<designation>&open=1`), now flies to it and shows the stars, black holes, neutron stars, quasars, nebulae, remnants and hypervelocity stars the scatters left in it, each in the colour the map gives its kind and each saying in its info panel that its sector is uncharted. The sector carries the mark in three places: its crumb in the title reads "(uncharted)", its info panel is headed "Uncharted sector" with the count of stars and other objects waiting, and a dashed "Uncharted sector ..." tag sits in the corner of the map. Generating the sector replaces all of it with the normal sector view. An empty ungenerated sector is selected as before.
- **`GET /api/galaxy/uncharted?ring=&layer=&slot=`.** One cell's place (designation, centre, edge) and what the scatters left in it (the waiting bright stars and the unbuilt scattered phenomena), with `sector_id` when the cell has been generated; the web server's `/galaxy/uncharted/<ring>/<layer>/<slot>/scene` turns it into the map's scene.
- `planetgen check-db --deep` now also loads and validates every star system. It prints how long it expects to take (from the speed this server recorded, or a rough guess it says so about), asks for a yes unless `--yes` is given, shows a progress bar, and lists each failing system. `--estimate-only` prints just the estimate. The Generate page's "Check the database" section has a "Deep check" option that shows the time and asks you to confirm. (DB.21)
- Neutron stars and black holes get their own mass limit (GEN.195): `--compact-min-mass` (1, 2, 4 or 6 solar masses, or `star` to use the stellar mass limit, the default) and a control beside the stellar mass limit in Plan and New galaxy, with the database cost shown. The central black hole or quasar is always created.
- A molecular cloud of the galaxy's cloud field (GEN.176) is born in the sector holding the centre of the space it fills, with a field-drawn serial worked out from the galaxy seed alone, so its object ID no longer depends on which sector is saved first. Run-time births (GEN.172) have tests for every path that exists today: a system added to a sector, a phenomenon or system saved on its own, a facility, and a body added by an edit; a number is never given again after its object is deleted.
- A "Habitability levels" page under Classes (`/classes/habitability`, UX.90) explains the five equipment levels, the four PHI-4 domains with their Blue, Green, Yellow and Red thresholds, the pulsar, neutron star and black hole host rule, and the Habitable, Habitable moon and Inhabited chips. It is linked from the Classes list, from a visible line above a system's body list and from the Equipment tags on the search page.
- Object IDs in the database (DB.20, GEN.171): `uid` is now the 80-bit birth-location ID (`BINARY(10)`, unique on its own) on every object table, and `facilities` gains the column. The sector fill gives a sector's systems and phenomena their generated serials and a system's bodies their body numbers as the rows are inserted; rows saved later take run-time serials and body numbers from the new `id_counters` table, whose numbers are never given twice. The schema is v78; existing rows are numbered by row order. The hashed IDs, the bright-sweep systems' position-ID `uid` and the per-system uniqueness are gone.
- API keys have scopes (API.9): `read`, `generate`, `upload` and `admin` (admin implies all; generate and upload imply read), an optional expiry, and a visible prefix. Existing keys become admin keys. A key short of a route's scope gets `403` with `required_scope`; every key has its own rate limit bucket and skips the per-address defaults; `last_used_at` is written at most once a minute. The admin page's new-key form takes scopes and a lifetime. Control schema v13.
- **The habitability score of every planet and moon (GEN.89).** Each rocky body now stores the PHI-4 display (Pressure, Temperature, Chemistry and Radiation, each a 0 to 1 score and a Blue, Green, Yellow or Red tier, and their overall mean), the microbial-life, complex-life and human-operability numbers behind them, and the equipment a human needs: shirtsleeve, mask, mask and scrubber, sealed suit or full life support. A pulsar's or neutron star's planet always needs full life support. The system page shows the equipment as a coloured chip on every planet and moon row (the tooltip gives the score and what limits it), and the search page has Equipment a Human Needs filters for planets and moons, with the equipment in their result tables. A cosmic-ray change from a nebula or remnant re-scores the body. Schema v76; no new random draws, so seeded output keeps its earlier values and gains these.
- The object ID layout (GEN.170): `galaxy/object_uid.py` packs, unpacks, prints and parses the 80-bit ID of birth sector, serial and body number (20 hex digits, `BINARY(10)`), and picks a longer layout (96 or 128 bits) for a galaxy whose bounds do not fit. Nothing uses it yet; the schema (DB.20) and the fill (GEN.171) come next.
- The Generate page shows the galaxy's layer specs (ADM.28): how many layers, how tall each is, how far the widest reaches, how many sectors are charted, and a row per layer (an even sample past 41). `GET /api/galaxy/shape` gives the same as `layers`.
- **Surface radiation dose, UV and galactic hazards (GEN.87).** Every rocky planet and moon now stores its yearly surface dose in mSv (cosmic rays through the air and any dipole, stellar particles, crust and radon, and a cosmic-ray boost while a nebula or remnant presses on its star), a DNA-weighted UV index after the ozone layer, and an ozone-loss flag. Each star stores how many lethal supernovae per Gyr its place in the galaxy brings. Schema v75. Seeded output differs from earlier versions (the crust draws one more value per rocky body).
- The mass limit (GEN.183) is now picked from the presets 8, 10, 12, 14, 16, 18 and 20 solar masses, 20 by default: a slider in the Generate page's new galaxy and plan forms, and `--phenomenon-min-mass` on the command line, which refuses any other value. Every star, neutron star and black hole at or above it is placed across the whole galaxy; lighter ones are drawn when their sector is made.
- Every rocky planet and moon now has a hydrosphere (GEN.88): its water share of its mass, where the water is (dry, in the air, frozen, under ice, open ocean, or a hycean ocean under hydrogen), how much of the surface is ocean and land, the ocean's depth, any ice lid and high-pressure ice beneath, and the ocean's chemistry class (ice-sealed, chloride brine, acid sulfate, soda or neutral) with its pH, water activity and phosphorus supply. Schema v74.

### Changed
- The "Rebuild the bright stars" action is gone from the Generate page; Redo scatters replaces it.
- TODO list: PERF.33 records the audit gaps in the ETA estimates (scatter floors, backfill count, command line, job tree).
- The Generate page's whole-job bar now estimates a `galaxy` or `plan` step from the stored time of each stage it will run, looked up by the settings that stage runs with (mass limit, luminosity floor, workers); a step with a stage not yet recorded still uses what the step took before.
- TODO list: GEN.186 retired (merged, PR #1088).
- TODO list: PERF.62 recorded as built (PR #1087); PERF.58 and PERF.61 note the early stop.
- TODO list: MAP.162 retired (merged, PR #1085).
- TODO list: DB.21 retired (merged, PR #1082).
- TODO list: GEN.195 retired (merged, PR #1081).
- The stellar mass limit (`--phenomenon-min-mass`) now governs stars only.
- TODO list: PERF.62 (interim early stop) dropped, superseded by PERF.58.
- Filed PERF.62 (interim 100-empty-layer early stop, already built by lane 1, superseded by PERF.58) and noted on PERF.58 and PERF.61 that it goes away with the old walk.
- Retired PERF.56 (stage timings stored with their settings, PR #1077).
- Each stage of a galaxy or plan run is now stored with how long it took and the settings it ran with (mass limit, luminosity floor, workers, radius ...) and what it did (layers visited and changed, objects placed); skipped stages are stored with their reason. The admin Stats page lists the latest stage times, and the estimate helper `stage_seconds` reads the runs with matching settings (PERF.56, control schema v14).
- Retired GEN.172, GEN.176, TEST.110 and DOC.5 (object IDs for run-time births, nebula birth sector, ID tests and docs, PR #1075). Filed the leftovers as GEN.197 (ID rules on ejection, merger and split) and DOC.17 (api.md object ID, after API.23).
- PERF.58: the shunting amendment is withdrawn; objects over a sector's capacity are dropped (Boss 10:09Z).
- PERF.58: objects over a sector's capacity are shunted to a face-adjacent sector with room instead of dropped (Boss 10:07Z).
- PERF.58: Boss dropped the stack-of-layers rule; the sampler runs per layer only.
- Retired MAP.166 (honest Dimmest star shown label, PR #1070).
- PERF.58: final stack rule (no layer grouping; per-layer sampler) and several objects per sector by capacity tiers, from Research lane 3 (PR #1068); PERF.61 follows the same tiers.
- Filed PERF.58 to PERF.61 (object-first scatter sampler, shared ring inputs, large-mean Poisson helper, phenomena sampler; top priority) from the scatter study; PERF.57 is superseded by PERF.58.
- Retired UX.90 (habitability explanation page and Ideal rename, PR #1065).
- The lowest PHI-4 equipment level is now called "Ideal" instead of "Shirtsleeve", on the system page, in the search tags and column and in the docs. The stored tier numbers are unchanged.
- PERF.56: GEN.187 needs no new stage entry; the existing backfill stage is relabelled (PR #1062).
- **The backfill stage is named for what it does now (GEN.187).** On a galaxy run's numbered stage list, "Backfill the bright stars" becomes "Scatter the massive stars from the neighborhood", the four mass rings around the generated sectors. It is still shown as skipped, with its reason, when `--backfill-from none` is used or the run generated no sector, and its timing and per-layer counts are logged as before.
- TEST.124 is first in Bugfixes lane 1 (priority: main is red).
- Retired GEN.187 (bright-star back scatter, PR #1059). Filed TEST.124 (Phenomena page tests fail with KeyError scattered, a bug).
- **The bright-star backfill goes by mass, in four rings (GEN.187).** Around the sectors a run generates, the sectors a face away (no diagonals) now get every star born with at least 1 solar mass, the next ring out 2, then 5, then 8; each ring counts from the previous ring's outer edge, and past the fourth ring only the scatter's own mass limit applies. The old luminosity tiers (100 to 750 L_sun out to 10, 25, 50 and 100 ly) are gone, as is the `--backfill-from` radius wording. A sector filled later builds its own stars only below the mass and luminosity already placed, a staged scatter tops backfilled sectors up with the lighter stars they lack, and white dwarfs born at 1 solar mass or more are placed too (stored with a `NULL` lifespan). Schema v79 (`sector_stats.bright_mass_sol`). Design: `docs/design/mass-backfill.md`.
- Retired UX.89 (stage lists and numbering, PR #1057); noted the stage entry GEN.187 must add.
- Every staged job now numbers its stages across the whole job (a New galaxy shows "Stage 7 of 12", not "Step 4 of 4") and lists them; stages a run will not do are listed as skipped with the reason (UX.89). The galaxy and plan commands print "Stage N of M" lines and record the stage in the progress file.
- Retired DB.20 and GEN.171 (object IDs in the schema and the sector fill, PR #1055); noted what is already built for GEN.172.
- PERF.57: after a group places something the next group is half the size (Boss 09:19Z).
- Rewrote PERF.57 to Boss's layer-grouping rule for galactic scatters (replaces the 100-empty-layer stop).
- Retired GEN.194 (defaults 14 Msun and 9,000 Lsun, PR #1051).
- The default mass limit is now 14 solar masses (was 8) and the default luminosity floor for the brightest stars is 9,000 solar luminosities (was 5,000); both are presets on their lists.
- Retired UX.88 (the "uncharted" wording sweep, PR #1049).
- Text shown to visitors says "uncharted" instead of "unbuilt" or "not generated" (the Phenomena table's rows, the nearby search and route notes, the Sector Map and Galaxy Map labels and hovers, the command-line query note); the Generate system and admin panels keep the technical words (UX.88).
- Filed UX.91 (full PHI-4 explanation in the planet and moon description) from Boss's request.
- Added Boss's layers-modified statistic to PERF.56 and a pointer on PERF.57.
- Added Boss's clarification to PERF.57: star passes count stars, the phenomena pass counts phenomena.
- Filed PERF.57 (stop a layer-walking scatter early after 100 empty layers) from Boss's request.
- **Docs only:** GEN.195 and GEN.196 gain details from Bugfixes lane 2 (every scatter path, redo scope).
- **Docs only:** GEN.196 (one Redo scatters box on Generate) filed for Bugfixes lane 2.
- **Docs only:** MAP.166 (honest "Dimmest star shown" label) filed; GEN.195 gains the CLI option and storage figures.
- **Docs only:** GEN.195 (separate compact-object mass limit; central black hole or quasar always created) filed for Bugfixes lane 2.
- The phenomena scatter's final line lists only the kinds it created, like the star scatter's, with no zero counts.
- **Docs only:** the plan notes schema v77 (Phenomena table class totals, PR #1038) and next Alembic revision 0078.
- The Phenomena table counts and pages the scattered, unbuilt phenomena from stored per-class totals (schema v77), after the built ones, so it stays fast with hundreds of millions of scatter rows.
- **Docs only:** UX.90 now covers the Shirtsleeve to Ideal rename and an explanation page under Classes, owner Bugfixes lane 2; docs use "Ideal".
- **Docs only:** MAP.148 notes the follow-ups to MAP.163 and MAP.164 (PRs #1032, #1035).
- **Docs only:** API.9 (key scopes, expiry, prefix, per-key rate bucket) retired; control schema is v13.
- **Docs only:** UX.90 (explain the habitability chips in the web interface) filed as an unassigned Phase 1 item; docs/html-interface.md notes the chips.
- **Docs only:** MAP.163 and MAP.164 retired; MAP.165 (store a mass for scattered phenomena) filed; MAP.148 notes the 400-star tile cap.
- Galaxy Map: black holes are drawn purple and neutron stars dark blue (MAP.164), at the top of the star scale with full core brightness so they out-shine brighter stars. Scattered black holes and neutron stars that no sector has built yet are now listed in the map tiles too, sized by mass class, and the biggest classes (the nucleus and intermediate-mass black holes) show from the whole-galaxy view.
- Galaxy Map: the "Dimmest star shown" slider now runs from the dimmest to the brightest star at the current zoom (an open sector's stars, else the loaded tiles', down to 2,500 solar luminosities at galaxy scale), 0 still showing every star (MAP.163).
- **Docs only:** GEN.194 (defaults 14 Msun and 9,000 Lsun), UX.89 (stage counts bug) and PERF.56 (per-stage timing with settings) filed for Bugfixes lane 1.
- **Docs only:** GEN.89 and its umbrella GEN.83 (habitability index) retired; schema is v76.
- **Docs only:** UX.88 moves to Bugfixes lane 1.
- A Generate-page job with several steps now shows one overall bar above the step bar, with the time left across all steps: the running step's own estimate plus what earlier runs of each later step took ("at least" when a later step has no record yet).
- **Docs only:** PERF.33 (job ETA adds unstarted steps, PR #1021) and PERF.55 (overall bar on the web pages, PR #1023) noted as partly built.
- **Docs only:** UX.88 (say "uncharted" for anything not yet generated, outside Generate and admin) filed; the wording rule is in docs/html-interface.md.
- A multi-step job's time left on the Queue page now includes every step that has not started yet, using how long earlier runs of that step took; a step with no earlier run marks the estimate as partial.
- **Docs only:** GEN.193 (Phenomena table empty after the scatter) filed and retired with PR #1019.
- The Phenomena table now lists every phenomenon the scatter has placed but no sector has built yet (black holes, neutron stars, nebulae, remnants, hypervelocity stars), as greyed "Uncharted ..." rows without a page. Before, a fresh scatter left the table empty until sectors were filled.
- **Docs only:** GEN.170 retired; TEST.123 (a flaky spatial-position test) filed.
- **Docs only:** PERF.53 retired; PERF.55 notes the missing whole-job layers-per-second stat.
- The generation-speed stats no longer count a sector, bright-star layer or phenomena layer that produced nothing, so empty ones stop skewing the per-sector, per-layer and per-unit times (the progress bar still counts every layer of the job).
- **Docs only:** GEN.192 retired (phenomena scatter in Generate-page jobs and its log output, PRs #1008, #1012, #1013).
- The phenomena scatter's output now matches the star scatter's: black holes are split into stellar, intermediate and supermassive, a "landed in N of M layers" (or "none landed") line, the sectors already filled that it leaves out, and a summary that lists every class including the ones that drew none.
- The phenomena scatter logs each layer's phenomena by kind as it finishes ("Phenomena, layer 3: placed 120: 90 neutron_star, ..."), and the special ones (nucleus, hypervelocity stars) on one line, as the star scatter does; before, only the final total was listed.
- **Docs only:** GEN.192 (phenomena scatter log and Phenomena table check, remainder of GEN.190) filed for Bugfixes lane 1.
- **Docs only:** GEN.190 and GEN.191 (Generate-page phenomena scatter, New galaxy mass limit) retired.
- **Docs only:** OPS.40 and GEN.188 retired; TEST.122 (browser map tests fail on plain main in one container) filed; MAP.163 notes the 400-star tile cap; the scatter design note now says the defaults are 8 solar masses and 5,000 L_sun.
- The Generate page's New galaxy, Plan and Rebuild the bright stars jobs now scatter the phenomena too (black holes, neutron stars, nebulae and the rest) after the stars; they ran only the star scatter, so a new galaxy had no phenomena scatter at all. New galaxy also passes the form's mass limit to its scatters (it used 20 whatever was chosen).
- On the Generate page, the mass limit slider and the luminosity floor dropdown now sit side by side, in Plan the galaxy and in the New galaxy section (GEN.188). The defaults are 8 solar masses and 5,000 solar luminosities.
- **Docs only:** GEN.191 (New galaxy ignores the mass limit slider, bug) is filed for Bugfixes lane 1; GEN.190 gets the cause found.
- **Docs only:** PERF.55 (one global progress bar with an ETA across the phases of a generation job) is filed as an unassigned Phase 1 item.
- **Docs only:** GEN.190 (Phenomena table empty after a web-generated galaxy, bug) and MAP.164 (Galaxy Map phenomena colors and visibility) are filed, both ASAP.
- **Docs only:** TEST.121 (PR #1002) and TEST.120 (fixed by PR #975) are retired; PERF.54 files separate stats rows for the star scatter passes as a Phase 2 item.
- test_open_map_menus_hold_no_overlap closes each menu in the page instead of clicking its button, so a busy machine no longer fails it on the closing click.
- A map Menu whose content arrives after it opens (the Galaxy Map's kinds and star filters) is placed again when it grows, so it no longer hangs over its own button. Menus opened above their button covered the button when the content was late.
- The default luminosity floor (GEN.184) is now 5,000 solar luminosities (was 3,000) and the default mass limit for massive stars and phenomena (GEN.183) is now 8 solar masses (was 20, the bottom preset). Both are still on their preset lists and the ranges are unchanged. A galaxy already scattered keeps the levels it was scattered at.
- **Docs only:** PERF.53 (timing stats skewed by empty layers and sectors) is filed.
- **Docs only:** MAP.153 (rank birth-radius fade, PR #998) is retired; GEN.189 files gamma-ray burst and AGN ozone loss as a Phase 2 item.
- The Galaxy Map's stars now fade in with the zoom (MAP.153). Each star gets a birth radius from its place in its tile's list (most luminous first), and its opacity rises smoothly over one halving of the camera distance, cross-fading from the coarser tile's rank across the octave a tile level serves. Zooming in adds stars a few percent at a time instead of up to eight times as many in one frame at a tile level change, zooming out removes them as smoothly, a late tile changes nothing visible, and a bookmarked view always draws the same picture. A dense sector now shows its dimmest stars only near sector zoom, brightest first. No server or database change.
- **Docs only:** ADM.49 (user-set density of spiral arms, inter-arm space, core and bulge) is filed as an unassigned Phase 1 item.
- **Docs only:** MAP.163 (Galaxy Map brightness scale) is filed; GEN.188 and MAP.163 are one-offs on Bugfixes lane 2; TEST.120 is owned by Bugfixes lane 1.
- **Docs only:** TEST.119 (PR #993) and PERF.52 (PR #994) are retired; TEST.121 files a load-dependent browser test failure.
- The Stats page's generation-speed table says what each row counts: a sector fill is seconds per sector and systems per sector, a bright-star or phenomena layer is seconds per layer and stars (objects) per layer, and so on. Layers showed "1,262 systems per sector" in a galaxy of 7,663 systems; they were stars per layer all along.
- A reload keeps a highlighted kind when a hidden-by-default kind (rogue planets) is also shown; taking the hidden kinds up rewrote the address before the highlight was read, and dropped it. The browser tests for the grouped Galaxy Map Menu (UX.86) find its button again and show rogue planets before checking a highlight.
- **Docs only:** ADM.28 (simpler Generate page, closes issue #736) and ADM.45 (star mix, PR #991) are retired.
- The Generate page's binary-system and wide-pair prevalence boxes become a star mix (ADM.45): the share of systems with one star, a close binary and a wide pair, which must total exactly 100%. The page keeps a running total and says which way to move; the server refuses any other total. The command line has `--star-mix SINGLE CLOSE WIDE` with the same rule (it replaces `--prevalence` for `binary_system` and `wide_binary`).
- The Generate page keeps its common actions on the page and moves the less common settings into a Customize window with one tab each (ADM.28): prevalence and override for Generate sectors; galaxy shape, prevalence and bright stars for a new galaxy. Without JavaScript the groups stack in the form as before.
- **Docs only:** PERF.52 (admin generation-stats table is wrong) is filed.
- **Docs only:** GEN.173 and GEN.174 (PR #987) are retired; TEST.119 and TEST.120 file two test failures on main.
- GEN.173, GEN.174: a planet, moon or belt added by an admin edit now gets a uid (it was saved with none), and after a delete the new body's uid skips ones its siblings already carry (it failed with IntegrityError 1062 on `uq_planets_uid`).
- **Docs only:** GEN.87 (surface radiation dose, PR #985) is retired from the TODO list and the plans.
- **Docs only:** Foundations lane 1 queue notes (ADM.28 and ADM.45 first, then the object-ID block and API.9); ADM.45 star mix decision.
- **Docs only:** the lane queue puts GEN.170, API.9 (moved to Foundations lane 1) and MAP.153 first.
- **Docs only:** OPS.39 (Windows support removed, PR #981), TEST.115, OPS.34 (superseded) and TEST.118 (PR #975) are retired.
- **Docs only:** OPS.40 (update.sh step 8 bug) is filed.
- **Docs only:** GEN.188 is owned by Foundations lane 2.
- **Docs only:** GEN.188 (mass limit default 8; mass slider and luminosity dropdown side by side on Generate and New galaxy) is filed, and the mass-cut design note records the new default.
- **Docs only:** ADM.31 (every generate action offers the Galaxy Map, PR #976) is retired from the TODO list and the plans.
- `GET /api/jobs/<id>` now gives `made_url` for a finished neighborhood or regenerate job: the Galaxy Map fitted to the sectors it made (ADM.31). The map menus' Generate buttons already reach the Generate page's job, which offers the same link.
- **Docs only:** OPS.39 (remove Windows support) is assigned to Foundations lane 3.
- **Docs only:** OPS.39 (remove Windows support, keep only docs/WINDOWS.md) is filed; TEST.115 and OPS.34 are marked superseded by it.
- TEST.116: when test_bughunt_end_to_end finds a star type other than K2V, the message lists every star of the system with its role and type, so a rare failure names its cause.
- TEST.115: the Windows CI job keeps Redis in WSL alive (it ran in the foreground of a wsl.exe it holds open), waits until Windows can connect and passes the working URL on to the tests. A job's lock is removed with a few retries, and three Windows-only test failures now say what state they saw.
- **Docs only:** GEN.183 (mass limit presets, PR #969) is retired from the TODO list and the plans; TEST.118 files two test failures that follow GEN.184's luminosity floor.
- A scatter or a phenomena-only re-scatter run without `--phenomenon-min-mass` now keeps the mass limit already stored with the galaxy instead of going back to 20.
- **Docs only:** OPS.38 (committed Redis dumps, PR #967), TEST.117 (already fixed by TEST.113, PR #920) and GEN.175 (regenerate keeps the uid, PR #960) are retired from the TODO list and the plans.
- Removed the Redis `dump.rdb` snapshots that had been committed (root and `src/`) and ignore `*.rdb`.
- **Docs only:** GEN.184 (luminosity floor presets, PR #965) is retired from the TODO list and the plans.
- The bright-star luminosity floor (GEN.184) is now chosen from presets: 2,500 to 4,000,000 solar luminosities on an exponential ladder (steps of 100 near 2,500, about 400,000 near the top), default 3,000 (was 1,000). Nothing below 2,500 is accepted by `--bright-star-min-luminosity` or the Generate page, which offers the presets in a list. A galaxy already scattered keeps the level it was scattered at.
- **Docs only:** GEN.88 (hydrosphere and ocean chemistry, PR #963, schema v74) is retired from the TODO list and the plans; OPS.38 (committed Redis dump files) and TEST.117 (generatejobs.test.mjs failing since PERF.33) are filed.
- Rogue planet oceans now stop at the depth where high-pressure ice forms (the rest is stored as high-pressure ice), their ice lid is compared with the water in matching units and melts lower under its own weight, and an ocean under a hydrogen envelope is shown as a hycean ocean. New seeded output differs from earlier versions.
- **Docs only:** Boss confirmed the defaults on GEN.183, GEN.184, GEN.187 and UX.87; their open questions are now decisions.

### Fixed
- The Galaxy Map's "Dimmest star shown" slider says "every star" only when the view holds every star there is (an open sector). At galaxy scale, where tiles are capped and some sectors are uncharted, it names the dimmest star the view actually carries (MAP.166).
- Galaxy Map: the "Dimmest star shown" slider now runs from the dimmest to the brightest star the current view's tiles actually carry, at galaxy scale too, instead of stopping at a fixed 2,500 solar luminosities (MAP.163 follow-up). 0 still shows every star.
- Galaxy Map: the galactic nucleus now shows from the whole-galaxy view whether it is a quasar or a supermassive black hole, and stays there once its sector is filled (MAP.164 follow-up). The nucleus quasar was missing from the scattered objects the map lists.
- `./update.sh` step 8 ("Checking the cache, jobs and debug log locations") no longer fails with `FileNotFoundError` on `util/appconfig.py`; `setup-debug-log.sh` now loads the log paths from `util/logpaths.py`, so the logs and their rotation are set up again. Two new tests (the log-location check runs against the checkout; deploy scripts only name repo paths that exist) fail on the old script.
- The star scatter pass tests use a luminosity floor the toy galaxy can fill with lighter stars, since the real floor now starts at 2,500 solar luminosities (GEN.184).

### Removed
- Windows support (OPS.39). `install.ps1`, `update.ps1`, `scripts/deploy-common.ps1`, `examples/maintenance/install-maintenance-task.ps1`, `examples/windows/` and `docs/deployment/windows.md` are gone, as are the Windows CI jobs (`windows-jobs`, `windows-installers`), the Windows branches in the code (the detached-process, `taskkill` and `OpenProcess` handling of Generate page jobs and the no-Redis fallback to run a job directly, the CPU-percent load reading on the admin queue page, `SpawnWorker`, the below-normal worker priority class, the checkout-relative log and settings folders, drive-letter disk measuring) and the Windows-only tests. The admin Generate page's jobs now always run on Redis, and the queue page shows only the load average. `waitress` leaves the `server` extra and `requirements-server.lock`.
- `docs/WINDOWS.md` gives basic instructions for a typical Windows setup (WSL2 and the Linux guides); anyone who wants native Windows does that work themselves.

## [8.0.866] - 2026-10-10

### Changed
- GEN.175: regenerating a phenomenon keeps its uid; it was set to NULL.
- **Docs only:** filed GEN.187 (bright-star back scatter by mass, issue #952), MAP.162 (open sectors that hold scattered objects but were never generated, issue #928) and UX.87 (uncharted systems in the system list, issue #929).
- **Docs only:** filed TEST.115 (Windows CI leg failures: Redis in WSL unreachable) and TEST.116 (a one-off K2V failure in test_bughunt_end_to_end) from the Bugfixes lane 1 CI findings.
- GitHub CI (`ci.yml`, every test leg) now runs only by hand: Actions > CI > Run workflow. It no longer runs on a push or a pull request. `stamp-version.yml` and `release-note.yml` are unchanged.
- CI: the `linux-update` job's "needs migrating" and "failed migration" steps rebuild the database from the v61 baseline fixture (the oldest schema `update.sh` upgrades from) instead of faking v48, which the code has refused since the Alembic cleanup. The old steps failed with "database is at schema v48, older than v61".
- **Docs only:** TEST.112 (the 8 s wait for queued API edits, PR #955) is retired from the TODO list and the plans; TEST.111 stays open with a note.
- TEST.112: tests now wait up to 120 s for a queued API edit instead of 8 s, so a slow worker start in a busy parallel run no longer turns a 200 into a 202 (the cause of the intermittent `test_regenerate_phenomenon_keeps_id_name_and_place` failure). TEST.111's assertions now print the results they checked, so its next failure names the cause.
- **Docs only:** GEN.185 (the five-pass scatter, PR #953, schema v73) is retired from the TODO list and the plans; the next free Alembic revision is 0074.
- The star scatter runs in passes (GEN.185). The mass pass places every star born at or above the mass limit (the same 8 to 20 solar mass limit as the phenomenon scatter's, default 20), whatever its luminosity. Each sector holding one at least as bright as the luminosity floor is marked, and the luminosity pass then places the lighter stars at least that bright and skips the marked sectors, which already hold a star that bright. A sector's own draw is lighter than the limit and dimmer than the floor, and the backfill and the staged bands below the floor draw only lighter stars too. This changes the stars a given seed gives: **re-plan the galaxy** (`planetgen plan`) to get the new scatter.
- `galaxy_shape.bright_star_mass_limit_sol` (schema v73) records the mass limit the star scatter used; a galaxy scattered before it keeps NULL and behaves as before.
- **Docs only:** GEN.97 (random neighborhoods, PR #950) is retired from the TODO list and the plans; GEN.186 files the leftover "keep away from filled space" option.
- **Docs only:** GEN.165 (a test bug about a comet round a second star, PR #948) is retired from the TODO list and the plans.
- GEN.165: the epoch-position test compared a comet round a second star with the barycenter; it now compares against what the comet goes round, like planets, so it no longer fails in some random systems. No product code was wrong.
- **Docs only:** TEST.114 (exact COUNT(*) in the check test, PR #946) is retired from the TODO list and the plans.
- TEST.114: the "check writes nothing" test compares exact row counts instead of MySQL's background-refreshed row estimates, so it no longer fails under load.
- The System Map no longer shows a comet on an open (parabolic) orbit once its pass is over (or before it arrives): beyond the drawn path, 2,000 AU out, it disappears instead of flying on forever. Closed-orbit comets always show.
- **Docs only:** recorded Boss's answers on NAV.8 (stars only get pages), NAV.11 (stay per stop defaults to 0 minutes, user-changeable) and API.23 (the object ID replaces row ids, API break accepted).
- **Docs only:** ADM.31 is partly built (PR #942); the remaining half stays open and MAP.151 notes the 2,000-sector cap.
- **Docs only:** the scatter-preset rush job is split across the Foundations lanes (GEN.183 to lane 2, GEN.184 to lane 1, GEN.185 to lane 3) and Bugfixes lane 1 is running again.
- **Docs only:** filed GEN.183 (mass cut presets, 8 to 20 solar masses), GEN.184 (luminosity floor presets, 2500 to 4 million L_sun, default 3000) and GEN.185 (the five-pass scatter order) for Foundations lane 3, and noted that Foundations lanes 1 and 2 are paused.
- **Docs only:** GEN.182 (comet ejection is by design; System Map panel says bound or unbound, PR #938) is retired from the TODO list and the plans.
- GEN.182: comets that never come back are by design (parabolic, about 30%). The System Map's info panel now says whether a comet's orbit is bound (and its period) or unbound and not returning, with perihelion and eccentricity; a test checks every closed comet orbit stays inside the star's Hill sphere.
- **Docs only:** UX.86 (compact grouped Galaxy Map Menu, PR #936) is retired from the TODO list and the plans.
- UX.86: the Galaxy Map's Menu is more compact. The "Show on the map" toggles sit side by side at their natural width in their own collapsible group, the action buttons wrap in a row, and the map's controls use smaller buttons.
- **Docs only:** ADM.30 (radial fills, PR #934) is retired from the TODO list and the plans, and TEST.114 (a load-sensitive test flake) is filed.
- **Docs only:** UX.85 (menus open where they can be seen, PR #932) is retired from the TODO list and the plans.
- TODO: UX.84 retired (PR #930); DB.21 notes the step registry.
- TODO: ADM.29 retired (PR #926); MAP.151 notes the span module it reuses.
- TODO: DB.15 retired (PR #924).
- TODO: PERF.50 retired (PR #922).
- TODO: PERF.51 and TEST.113 retired (PRs #919, #920); PERF.50 notes the worker relay.
- Every generation and maintenance step reports its progress through one helper (`planetgen.generation.steps`): a bar is drawn by itself at the start of a step predicted to take longer than 15 seconds from the speeds this server recorded for its kind, nothing is drawn for a shorter one, and a step that runs past 15 seconds gets its bar at that moment. A step with no recorded speed yet is treated as long, so a first run still shows its bar. Finished steps record their speed back, so the next prediction has history; the Stats page lists each kind.
- The sector bars, neighbour linking, sector paths, bright-star backfill and top-up, the phenomenon scatter (now timed in the bar's own units, so its bar starts from a recorded rate too), the population pass, `reset` and the orbit update draw through it. Steps inside a worker can report to the parent through a channel (`worker_step` and `Relay`).
- TODO: GEN.182 notes that the comet ejection may be a rendering issue.
- TODO: bug GEN.182 for comets ejected into space instead of returning.
- TODO: GEN.96 retired (PR #915); GEN.179 to GEN.181 and bug TEST.113 filed.
- A rocky planet that no tide has slowed now turns in 8 to 48 hours, drawn log-uniform as giant-impact formation models give (Kokubo and Genda 2010), instead of anywhere from 10 to 1,400 hours. Slow days now come only from tidal locking. Most unlocked rocky planets with a working dynamo now get a dipole field instead of a multipolar one (GEN.86 follow-up, Boss's decision of 2026-10-10). Rogue planets draw the same days. Seeded output changes for planet days and fields.
- TODO: ADM.42 retired (PR #912); ADM.43 records the read side that exists.
- A wrong value in `config.json` (a port out of range, a negative proxy count, an unknown `log_rotation`) now stops the program with a message naming the field, instead of being read as given. `PLANETGEN_ADMIN_COOKIE_INSECURE` reads text like `PLANETGEN_DEBUG` does (`false`, `0`, `no`, `off` and empty mean off) instead of only `1`.
- TODO: PERF.33 records what PR #910 built and what remains.
- The time left on every progress bar now starts from the rate this server recorded for the same kind of work and worker count (PERF.32), so a bar has an estimate before its first unit finishes, and blends in the live rate as units finish (weight n / (n + 15), held back for the first 5 units and 20 seconds). The decay time constant follows the task length.
- The Generate and Queue pages show the time left as a range ("about 1 m 28 s to 1 m 54 s left"), "estimating" while there is none, and hold it when nothing has finished for a minute or more; the Queue page's time left blends the recorded task time into the job's own pace.
- TODO: GEN.86 retired (PR #908); GEN.177 and GEN.178 filed for its two unmodelled parts.
- TODO: TEST.112 filed for an intermittent phenomenon-regeneration test; API.9 control migration number moved past v12.
- TODO: PERF.32 retired (PR #905); PERF.33 records what PERF.32 left for it.
- Recorded rates are deleted by the first run of a new version, since they describe the release, Python and machine that measured them. Control schema v12 replaces the old `generation_stats` table.
- TODO: UX.85 and UX.86 scope confirmed as the Galaxy Map.
- TODO: UX.85 and UX.86 for menus opening out of sight and oversized map controls.
- TODO: items DOC.6 to DOC.16 for a static Help section and one help page per web interface feature area.
- TODO: the progress-bar chain (PERF.32, PERF.33, PERF.51, PERF.50, DB.15, UX.84) moves to Bugfixes lane 1.
- TODO: UX.84 records the corrected survey of which generation steps have progress bars.
- TODO: the automatic progress bar prediction (PERF.51, UX.84) and the deep check estimate (DB.21) read the recorded generation statistics (generation_stats, PERF.32), with a conservative fallback only for a step with no history.
- TODO: replaced the hand-picked progress bars with the rule Boss stated (a bar starts by itself on any sub-step predicted over 15 seconds): UX.84 (bug), PERF.51 (the mechanism), PERF.50 reworked as its first application, PERF.33 and DB.15 moved into Phase 1.
- TODO: filed DB.21, a deep pass for the database check (every star system validated) with the estimated time shown first.
- TODO: retired UX.83 (PR #895); filed PERF.50 (a bar inside one sector save, from the worker-to-parent channel) and TEST.111 (a parallel-run test flake) in Phase 1.
- Linking new sectors to their neighbours now counts three steps per sector (containment, nearest systems, the neighbours' lists) and names the one under way, so the bar moves during a long batch.
- The phenomenon scatter has a progress bar of its own, weighted by each layer's expected work, in the terminal and on the Generate and Queue pages.
- The population pass shows a bar over its steps.
- TODO: retired DB.8 (PR #893, planetgen check-db and the Generate-page job).
- TODO: moved MAP.161 (first-visit script load) into Phase 1, on Foundations lane 2 next to MAP.157 and MAP.158.
- TODO: added the wire format detail from Research Lane 3 and Research Lane 2 to MAP.157 to MAP.160 and filed the first-visit load item MAP.161.
- TODO: recorded the measured Galaxy Map wire format facts on MAP.147 and split its recommendation into MAP.157 (trim the tile JSON), MAP.158 (gentler prefetch, IndexedDB), MAP.159 (packed binary tiles, with MAP.154) and MAP.160 (deferred GPU buffer quantising); fixed prerequisite lines on MAP.146, MAP.148, MAP.151 and MAP.152.
- TODO: recorded Boss's approval of the MAP.156 default.
- TODO: filed the top-down view of one layer or a range of layers as a secondary Galaxy Map option (Phase 2).
- TODO: UX.83 lists the generation steps found without a progress bar.
- TODO: filed the missing progress bars for neighbour linking, the phenomenon scatter and other long generation steps as a bug.
- TODO: recorded that the rank fade (MAP.153) is stage 1 and the magnitude law (MAP.148) the end state, as decided by Boss.
- TODO: recorded Boss's answers on the zoom star visibility defaults (MAP.153, MAP.155).
- TODO: moved the fly-through Galaxy Map (MAP.146 and MAP.147 to MAP.155) into Phase 1 as the headline of the 8.1 release; Phase 1 is now complete only when they are done, so the 8.0 hold stays until then. Placed in Foundations lane 2's queue (step 4) in dependency order.
- TODO: retired UX.35 (PR #879); recorded the dependencies between the zoom visibility items (MAP.148 builds after MAP.153, MAP.152 needs MAP.154, MAP.151 needs MAP.147) and added section 8a to the fly-through design note.
- TODO: filed MAP.153 to MAP.155 from the zoom star visibility research note (rank birth-radius fade, nested tile lists, other objects fading), and marked them as the first stages of MAP.148's visibility rule.
- Rewrote MAP.146 as the fly-through umbrella (scroll-zoom, double-click flight, distance-based visibility, see-through near field) and split it into MAP.148 to MAP.152; added the two research reports under docs/design.
- The NAV page's route now runs left to right with the distance of each hop after its stop, wrapping onto more lines as the panel narrows, and a route of more than nine stops shows its first stop, last three stops, longest hop and every jump through unknown space, with the whole route under "All N stops" (UX.35).
- Retired OPS.14 (PR #876) and closed DB.19, which the 20 solar mass cut (PR #866) made unnecessary.
- GEN.170 records Boss's decision: the object ID is 80 bits (20 hex digits).
- Retired PERF.49 (PR #870) and corrected the measurement note on PERF.31 and in the generation performance study.
- Parallel generation workers no longer queue for the system-name registry: each save claims its names in a short transaction of its own instead of holding the registry's row locks until the whole sector commits. In a dense core, 8 sectors on 4 workers took 54 s instead of about 90 s. A claim is given back if the save fails, and a first holder that a later save wanted to rename renames itself when it commits (PERF.49). Names, bodies and registry rows are the same as before on a seeded run.
- Name reservation looks a name's candidates up by key instead of scanning the whole sector for each name, and the offensive-word check is one compiled pattern (PERF.49).
- The database connection pool allows 10 connections plus 10 overflow, since a save now uses two.
- Retired API.20, API.21, GEN.138, GEN.147, MAP.144, OPS.32, PERF.37 (PR #865) and OPS.36 (PR #868).
- Retired GEN.166, GEN.167 and GEN.168 (PR #866); DB.19 is now conditional on lowering the mass cut to 10 solar masses or less.
- The phenomenon scatter places only the neutron stars and black holes of at least 20 solar masses (GEN.166 to GEN.168), so `planetgen plan` writes about 2.7e5 phenomenon rows instead of about 1.17e9. A sector draws the lighter ones itself when it is filled, from its own stream, so the galaxy holds the same number of each. `--phenomenon-min-mass` sets the cut and `--phenomena-only` re-scatters at a new one; the cut is stored with the scatter (schema v71) and in the settings file. Stellar-mass and intermediate-mass black holes are now drawn as separate kinds. This is the second half of the one-time reseed that began with lazy names (PERF.43): the same seed now gives a different galaxy than before both changes.
- Retired PERF.44, PERF.45 and ADM.48 (PR #863); PERF.46 keeps only its finite_domain question.
- A saved sector's rows now get their unique ID in the INSERT instead of being selected back and updated afterwards, which saves a query pass per sector (PERF.44). The IDs are the same.
- A `galaxy` run now links its new sectors to their neighbours (containment in nebulae and remnants, nearest systems) once at its end instead of one sector at a time while it saves, so the workers no longer queue for it (PERF.45). The stored links are the same. Sectors made on demand still link as they are saved.
- A planet's or moon's position is worked out once when first read instead of after every move (PERF.46); the values are unchanged.
- Retired PERF.43 (PR #861) and filed ADM.48, the failing auth-sweep tests for the galaxy-settings download.
- Stars and phenomena draw their generated name only when it is first read (PERF.43). A placed phenomenon is named by its object ID and never draws one, which saves about 4 percent of a sector fill; building a rogue planet is about four times faster. This is the first half of a one-time reseed: each name now takes one draw at construction instead of many, so the same seed gives a different galaxy. The second half is the 20 solar mass phenomenon cut (GEN.166 to GEN.168).
- OPS.37 and API.22 record Boss's decision on their defaults.
- Recorded the decisions on API.17 (keep the fingerprint check), OPS.16 (keep the daily positional update), DB.17 (keep repair from the seed) and the new-galaxy form (no seed field).
- The TODO list records the seeding requirements for lazy names and the mass cut (PERF.43, GEN.167, GEN.168) and the ordering note on GEN.57.
- The plan drops the user-facing rebuild of a galaxy from a seed and a version: OPS.12, GEN.59, GEN.61, OPS.18, ADM.17, ADM.19, ADM.20, API.16 and DB.10 are removed, and GEN.55 becomes the internal same-seed umbrella. The seed stays an internal mechanism for parallel workers, fills, backfills and settle.
- The TODO list and plan record Boss's decision that lazy names (PERF.43) and the 20 solar mass cut (GEN.166 to GEN.168) share one combined reseed.
- The TODO list corrects the name-reservation timing quoted on PERF.31, PERF.43 and PERF.49 (the 30 s figure was a contaminated measurement; alone it was 1.7 s of a 17.8 s sector), and the execution plan moves PERF.49 behind PERF.44, PERF.45, PERF.47, DB.19 and OPS.14, with a re-measure as its first step.
- The TODO list files OPS.36, a bug: the size and free-space checks measure the boot drive instead of the drive that holds the database.
- The execution plan moves the generation speed items to the front of Foundations lane 1 (PERF.44, PERF.45, PERF.49, then PERF.47 and DB.19) and records that Boss approved the PERF.49 plan.
- The TODO list retires GEN.85 (atmosphere species and mantle redox, PR #848).
- The TODO list files PERF.49 (batch the system-name reservation: 30 s of an 84 s dense core sector) with the earlier naming-cost findings re-checked against main.
- A planet's atmosphere text is now written from its gases ("a mix of nitrogen, oxygen, and argon, with traces of water vapor and carbon dioxide"), and its surface pressure drops by any gas too cold to stay in the air.
- The TODO list retires ADM.47 (PR #846) and records Bugfixes lane 1's measurement that name reservation takes 30 s of a dense sector's 84 s (PERF.31, PERF.43).
- The TODO list records that Boss accepted the 20 solar mass cut for the phenomenon scatter (DB.19, GEN.166 to GEN.168), and the execution plan moves those items to the front of Foundations lane 2.
- The TODO list retires GEN.137, NAV.53, OPS.33, UX.79 and UX.80 (PR #843).
- The TODO list retires NAV.12 (unbounded routes, PR #838) and PERF.42 (warm queue worker, PR #841).
- The queue worker loads the generation code and the bright-star sampling tables once, before it forks a work horse per job, so each job no longer spends about 2.6 s importing and rebuilding them (PERF.42). A worker also restarts itself when an update changes the release.
- The TODO list records Boss's decisions on the generation performance items: the nearest-system and containment work moves to its own phase (PERF.45), and skipping the finite-domain check stays open (PERF.46).
- The TODO list files the phenomenon scatter mass cut (GEN.166 to GEN.169) and reworks DB.19 around it: at the recommended 20 solar masses the scatter table falls from 1.17 billion rows (161 GB) to about 2.7e5 rows. The design notes drop the unverified 1.6e8 rows and 21 GB figures.
- Two systems in one sector are routed by their nearest stars even when those lie in the sector next door.
- The Generate page progress-line to-do item (ADM.46) was withdrawn at Boss's word.

### Added
- Random neighborhoods (GEN.97): `planetgen galaxy --neighborhoods N` (random-start mode) generates N neighborhoods instead of one. Each start is inside the galaxy with its whole neighborhood, and at least twice the radius from every other start, so they never overlap; one estimate covers them all. `--neighborhood-gamma G` biases the starts toward dense space (a start is kept with probability min(1, density) ** G; 0, the default, is uniform by volume). The Generate page's "Around a random start" has Neighborhoods and Density bias fields.
- A finished galaxy-generating job on the Generate page offers "Show on Galaxy Map" (ADM.31): the map opens fitted to the sectors that run made, each ringed, with a line saying how many (the first 2,000 are ringed on a bigger run). `/galaxy?made=<since>,<until>` and `GET /api/galaxy/made` do the work, by when the sectors were created.
- Radial fills (ADM.30): `planetgen galaxy --center-sector ID --cylinder-sectors X [--cylinder-layers H]` (or with `--ring --layer --slot`) generates a round disc of radius X + 0.385 sector edges around the centre, through H layers either side (default X). X = 1 is the centre and its in-plane face neighbours. The Generate page's "Around a sector" mode takes a radius in sectors and layers either side as an alternative to the radius in parsecs.
- Every long step of a generation or maintenance run now draws its own bar by itself when the recorded speeds predict more than 15 seconds (or it has none yet): the phenomenon scatter's clearing, special rows, their insert and the epoch stamp, the containment and nearest-systems refreshes outside a galaxy run, the bright-star scatter's main bar, the Galaxy Map warm-up and the name-dedupe passes. `steps.STEP_KINDS` lists every timed kind (the Stats page uses it for its names), and `tests/test_step_registry.py` fails when a step is built with a kind that isn't registered (UX.84).
- Span fills (ADM.29): `planetgen galaxy --rings 3:5`, `--layers=-1:1` and `--slots 50:5` (an arc, wrapping through slot 0, inside one ring) generate every missing sector in a range of rings, layers or both, inside the galaxy's outline. Ranges are inclusive, and the sector count comes from prefix sums so the size warning is instant. The Generate page has a "A span of rings, layers or slots" mode. `galaxy/span.py` holds the shape for the Galaxy Map's region data layer to reuse.
- The migration script shows a bar over its steps and, for a revision that works in batches (`alembic_runner.report_progress`), a second bar for that step, both with the time left; when output is not a terminal it prints a line per step and, for a step in batches, at most one line every 30 seconds and at each 10 percent. Migration statements wait up to an hour for a lock (DB.15).
- A dense sector's save now has a bar of its own, under the sector bar, in the terminal and on the Generate and Queue pages: a worker reports the systems it has saved to the run, which draws the bar when the recorded save speed predicts more than 15 seconds (PERF.50).
- Generation directives (GEN.96): `--directive systems>=N`, `habitable>=N` or `type:X>=N` on `sector` and `galaxy` (with `--directive-attempts`, default 200) redraw each sector until it holds at least that much, from a repeatable per-attempt seed. If no draw meets it, the closest is kept and the log says which minimum it missed. The Generate page's "Generate sectors" form has an Override section for the same.
- One settings model describes every `config.json` option (ADM.42): type, default, help text, unit, secret, restart and editable flags, and the environment variable that overrides it. `python -m planetgen.cli.config check` validates the files (a misspelt option name is an error there) and `config docs` writes the option table in `docs/config.md`, `config.json.example` and `config.schema.json`. A web-owned `settings.json` overlay is read when present, holding only options marked editable from the web. Tests fail when the generated files drift or when code reads a `PLANETGEN_*` variable the model doesn't name.
- Stars store their activity (GEN.86): the coronal X-ray share of their light, their X-ray plus EUV output, whether the corona is still saturated, how often they flare above 1e33 erg and the XUV they have given off over their life. M dwarfs stay saturated for billions of years. Hot stars and white dwarfs add their photosphere's ionizing output. Neutron stars and black holes give off their thermal, spin-down or accretion X-rays.
- Planets and moons store a magnetic dipole moment from their mass, density, age and rotation, a dipole class (none, weak, earth-like, strong or multipolar), the magnetopause standoff against their star's wind, and the star's XUV flux, lifetime XUV exposure and flare irradiation at their distance (schema v72). These are the inputs of the radiation dose (GEN.87) and the habitability index.
- Generation rates are recorded per kind of work and worker count (PERF.32): sector fills, bright-star layers and now phenomenon layers. The Stats page shows the worker count of each rate and has a Reset stats button; the time estimate reads the rate for the run's own worker count, blending the neighbouring counts when it has none.
- `planetgen check-db` and a "Check the database" section on the Generate page (DB.8): a read-only check of the schema version and models, table health, rows whose parent is gone, ids that would clash, impossible values, sector counts and version keys. It ends with a pass or fail line per check and exits 1 on damage and 2 when a check could not run.
- A warning when the running code differs from the code the galaxy was planned with (OPS.14): release, Python, platform, version key or a changed `requirements.lock`, read from the galaxy's settings file. It appears with the mixed-sector warning in `planetgen galaxy`, now also in `planetgen fingerprint`, and on the Generate page.
- Filed the object ID items from Boss's decision of 2026-10-09 22:39Z: GEN.170 to GEN.176 (layout, fill, run-time births, nebula birth sector and three bugs), DB.20 (schema), API.23 (public reference), TEST.110 and DOC.5, all in Phase 1. Added docs/design/object-id-options.md.
- Filed MAP.147 (the Galaxy Map wire format: measure the payload, compare options) in Phase 2.
- Filed MAP.146 (zoom drill-down centred on the clicked point) in Phase 2.
- Filed OPS.37 (a Generator version number) and API.22 (an API version number), both plain sequential integers due by the end of Phase 1.
- Planets and moons store their mantle redox (reduced, intermediate or oxidized, with its offset from the iron-wustite buffer) and the partial pressures of O2, CO2, CO, N2, Ar, H2, H2O, CH4, H2S and SO2 (GEN.85, schema v70). Each class's mix shifts with the redox, and no gas exceeds its vapour pressure at the surface temperature.
- A route now reports its longest hop and flags each hop whose line crosses sectors that have not been generated as unknown space (NAV.12). `/api/nav` returns `route.hops` and `route.longest_hop_ly`, and the NAV page states the longest hop and the unknown-space jumps.
- Eight to-do items from the generation performance study (a warmed worker, lazy names, uids in Python, the phenomenon-row size decision, a later pass for links, one position per body, benchmark records and a faster INSERT).
- A to-do item (generating a neighbourhood from the Generate page shows no per-sector stats).
- A to-do item (the Generate page shows a progress line and per-layer counts instead of one line per sector).

### Fixed
- The Galaxy Map's Menu, its Steps list and the Bookmarks list now open inside the visible window: a panel that would run off the bottom opens above its button when there is more room there, or is capped to the room left and scrolls inside, and one that runs off a side slides back in (UX.85). A browser test opens each in a short and a narrow window.
- The Generate page's job-status script test used the old `remaining_text` field, so it failed on main since the time-left range text (`remaining_label`) arrived (TEST.113).
- The disk-space check and the admin Stats tile measure the drive that actually holds the database's data directory (asked of the server with `SELECT @@datadir`, symlinks and mounts resolved), not the boot drive. They show which path and drive were measured, and say "unknown" when a remote server's disk can't be reached.
- A moon's Hill sphere (and so its minimum orbit spacing and the facility orbit slider) is measured about its planet, not its star (GEN.138); the Hill sphere uses the pair's total mass.
- Class N (a Venus analog) no longer carries a life chemical or a life timeline; Class Q is a habitable class capped at microbial life (GEN.147).
- The progress ETA is a ratio of decayed sums, so the early estimate of a run on several workers is no longer up to twice too long (PERF.37).
- The phenomenon render uses `THREE.Timer` in place of the deprecated `THREE.Clock` (MAP.144).
- The macOS update daemon's plist is well-formed XML again (OPS.32).
- A 100,000-deep nested JSON body is a 400 ("nested too deeply") instead of a 500 (API.20), and `Retry-After` is only sent on a 429, no longer on every response Flask-Limiter counts (API.21).
- The admin page's creation-settings download answers an anonymous or non-admin caller with a plain 403, like the other file and data views, instead of a redirect (ADM.18).
- A neighbourhood or single-address run now reports the first sector (the one named outright) with the same stats as every other sector, each summary has a "Totals" line (star systems, stars, planets, phenomena), and a batch says how many sectors it runs with how many workers before the first report (ADM.47).
- Hypervelocity stars now move in a straight line (position plus velocity times time) when the orbit update runs, keeping their velocity, instead of being turned about the galactic axis like bound objects (GEN.137). `phenomenon_scatter` gains `epoch_unix` (schema v69), the time a scattered star's position holds at, and the built system inherits it.
- The sector search reached by a radius now covers every cell that touches the sphere, not only cells whose centers lie inside it (`cells_touching_sphere`, NAV.53); the center-based listing's misses are documented.
- Python and the browser now round half-way numbers the same way (away from zero, on the shortest decimal), so 9.995 reads "10" and 1.005 reads "1.01" in both; Python used to give "9.99" and "1".
- A negative number that rounds to zero prints "0", not "-0", in both.
- A checkout with `core.autocrlf=true` no longer changes the bytes of the lock files and the word list, so their hashes agree between Windows and Linux.

### Removed
- `planetgen.util.appconfig`; the log locations moved to `planetgen.util.logpaths` and every other option to `planetgen.util.settings`.

## [8.0.783] - 2026-10-09

### Added
- A to-do item (a comet's scene position at time 0 occasionally disagrees with its stored position).
- Seven to-do items for modelling globular clusters (metallicity, cluster table, density, bright-first fill, planet cull, pulsar planets, synthetic clusters).
- Every rotating object stores a spin axis and an axial tilt (GEN.104, schema v68): stars, planets, moons, comets, rogue planets, interstellar comets, and standalone black holes and neutron stars. Stars, comets, rogue planets and interstellar comets also store a rotation period. Cool stars spin by gyrochronology, hot stars by a log-normal speed held under breakup, small bodies never faster than the 2.2-hour spin barrier, and black hole spin follows Beta(1.4, 3.6).
- 64 to-do items from the research handoff (12 of them bugs), about 100 research notes on open items, and a Documentation section in the to-do list.

### Changed
- The globular-cluster items link their design note and record the sector-address answer.
- The globular-cluster items record Boss's answers on the catalogue file and the planet cut.
- The multi-star systems item (GEN.129) carries the exact triple-stability equation from Vynatheya et al. 2022 and its tests.
- The to-do items for nebula planets, multi-star systems and exotic star systems cite the verified orbit-expansion law and triple-stability criteria, with unit tests.
- The spin vector and axial tilt item (GEN.104) is finished.
- Planets now tidally lock to their star by the same rule moons use; a locked body turns once per orbit, upright. Seeds give different systems than before.
- The nebula planet-formation study (GEN.94) is finished and its rule table stands.
- The to-do items for nebula planets, multi-star systems and exotic star systems carry the exotic-environments research.
- The version stays on 8.0 until Phase 1 is complete: patch and minor changes keep the revision and join the 8.0 entry.

## [8.0.711] - 2026-10-09

### Added
- A to-do item (MAP.134) to warm the Galaxy Map's opening view ahead of time on every update.
- The path each star system, rogue planet and interstellar comet takes through its sector is saved (`sector_paths` and `sector_path_knots`, schema v62): worked out at each orbit update for every sector holding one of them (not at generation, where they would depend on the order workers fill neighbouring sectors). Run `sudo ./update.sh` to migrate.
- The Galaxy Map's opening view (its first tiles and the galaxy stage) is built into the tile cache in the background at the end of every `update.sh` / `update.ps1`, so the first visit after an update no longer waits for the database. `python -m planetgen.cli.warm_map` does the same by hand, e.g. after clearing the tile cache (MAP.134).
- A to-do item (GEN.126) to run an orbital update as the last step of every generation run.
- A `galaxy` run now ends by saving the sector paths of the sectors it created and of the sectors around them, once the neighbours are final, so a freshly generated galaxy has every body's velocity, sector address and sector path filled in without waiting for the first orbit update. The same seed still makes the same galaxy at any worker count. About 1.8 ms a body; `--no-settle` skips the step.
- The same final step now runs after every edit that changes the masses in a sector: `planetgen phenomenon --sector-id` (skip with `--no-settle`) and the admin API's regenerate-sector, delete-sector, change-star and regenerate- or delete-phenomenon jobs. The two deletes now run on the queue like the other edits (answered at once if they finish in 8 s, else `202` with a job id).
- A sector generated on the spot (`ensure_sector_generated`) queues a job that saves the paths of that sector and its neighbours (run inline when no Redis answers); the neighbourhood job settles every sector it made at its end.
- To-do items for the open GitHub issues: 26 new items (11 of them bugs) in phases 1 and 2, with #736 folded into ADM.28 and #677 into NAV.45.
- A TODO tree page (`docs/plan/todo-tree.html`, built by `scripts/build_todo_docs.py`) that shows every open item under the prerequisite that finishes last, and a refresh of the architecture map for the modules added since 2026-10-07.
- A stand-alone facility stores a galactic velocity (`facilities.velocity_x_kms`, `_y_kms`, `_z_kms`, schema v64): the rotation curve's tangent at its place, filled when it is added and turned with its position at each orbit update, so every object in space carries a vector. Existing stand-alone facilities get theirs in the migration; the facility API returns the vector.
- A to-do item (ADM.45) for the Generate page prevalence fields: type the override percent, the form forces the set to total 100%, and the page explains what to do.
- A to-do item (GEN.131) for the bright-star scatter to log how many stars it added to each layer, by type.
- `planetgen plan` now places the galaxy's black holes, neutron stars, planetary nebulae, supernova remnants and hypervelocity stars, and its nucleus, galaxy-wide just before the bright-star scatter (schema v64, `phenomenon_scatter`). A sector's fill builds them from their rows at the stored points, no longer rolls those kinds itself, and keeps them where they were put. Molecular clouds stay on the galaxy's seeded cloud field; rogue planets, comets, asteroid fields and brown dwarfs are still rolled per sector.
- Every 3D view can be re-centered (MAP.138): C, or Menu's Center on selection on the Galaxy Map and Sector Map, Center on the system page's 3D view, and Re-center (or C) on the nebula page, put the center the view turns and zooms about on the selection (or the middle of the view), keeping the angle and zoom. The mouse then turns and zooms about it.
- The Work Queue's Jobs list shows a "Time left" column for running jobs, from the job's measured pace or the estimate a Generate page job publishes; it stays blank when none can be told.
- A "Clear finished jobs" button on the Work Queue page (and `POST /api/admin/work/clear-finished`) deletes every finished job's record in one step; running jobs stay.
- The habitability index's reference maths (`planetgen.physics.habitability`, GEN.84): the PHI-4 display (Pressure, Temperature, Chemistry and Radiation, each a score and a Blue, Green, Yellow or Red tier), PHI_bio, PHI_cpx, Phi_tech and the equipment a visitor needs, with every constant and threshold fixed and sourced in `docs/design/habitability-index.md`. Nothing calls it yet; GEN.89 stores the scores.
- `GET /api/near`, `planetgen.db.near.objects_within` and `python -m planetgen.cli.query near` (NAV.43): everything generated within up to 50 parsecs of a place (an object reference such as `system:12`, or a galaxy-frame point), of every kind: systems with their stars, planets, moons, belts and comets, facilities and the standalone phenomena. Nearest first, paged, with a kind filter and a count of the sectors in range that are not generated yet. The search never generates anything.
- A "What's nearby" page (`/nearby`, NAV.44): pick a place (a sector then a system, or a point), type a distance in parsecs, and list what is there, with a kind filter, paging, and a count of the sectors in range not generated yet. Linked from the menu and from every system, sector and phenomenon page.
- To-do items (GEN.132, GEN.133) for per-sector regional rates for neutron stars and black holes, and an analysis of every star type's rate against its distance from the galactic core.
- The bright-star scatter, the backfill and the band top-up now log, for each layer, how many stars they added and of what kind (for example "Bright stars, layer 0: added 1,204 stars: 310 B-type, 55 K-type giants, ..."), in the same words as the sector fill's summary.
- Neutron stars and black holes are now placed by region (GEN.132): both thin out in the galactic plane and spread to greater heights than the stars, black holes gain toward the core, and the share of neutron stars that are active pulsars follows the Lorimer radial profile. These are estimates; see `docs/design/compact-remnant-regions.md`.
- `docs/design/star-types-by-galactic-radius.md` and `scripts/star_type_by_radius.py` (GEN.133): how every star type's share changes with distance from the core, with proposals.
- `planetgen fingerprint` (GEN.58): a canonical SHA-256 digest of each sector's generated content and one for the region (`--ring`, `--sector`, or the whole galaxy with its plan), skipping row ids, timestamps, update clocks and what is rebuilt from positions, so two builds of one galaxy can be compared sector by sector.
- Every update now records the version key it ran under, one row per planned galaxy in the control database (`version_key_history`, control schema v11): the galaxy's seed, the key, the release and the SHA-256 of `requirements.lock`. The last 10 rows per galaxy are kept. `planetgen versions` (or `python -m planetgen.cli.version_history --list`) lists them (OPS.13). The galaxy seed itself never changes on an update.
- Every sector now records the PlanetGen version, Python and platform that generated it (schema v67, Alembic 0067: `sectors.version_key`, `planetgen_version`, `python_version`, `platform`). `planetgen galaxy` and the Generate page warn before extending a galaxy whose sectors came from a different version, naming what differs; `GET /api/galaxy/shape` carries the same text as `version_warning` (DB.7).
- `corridor.objects_near_segment`: every generated system, star and phenomenon within a distance of a line segment, in order along it (for NAV.6's course steering).
- A to-do item (GEN.134) to tune the star populations to the observed star-formation profile by galactic radius, from the proposals of the star-types-by-galactic-radius analysis.
- Planning a galaxy writes its creation settings to a JSON file (ADM.18): every plan option, the seed, the version key with its parts, the naming key, the `requirements.lock` hash and the name generator's word lists. Planning again with other settings keeps the earlier file as a dated backup. The Admin dashboard has a new "Galaxy settings" panel that lists the files and downloads any of them.

### Changed
- System and wiki Markdown is rendered by the `markdown` package (cut down to headers, pipe tables, paragraphs and `<sup>`), replacing the hand-written `mdconvert.py`. The pages look the same; table alignment colons now give an `align` attribute instead of being ignored (UX.39).
- Kepler's equation is solved with scipy (Brent's method, and a vectorised solver the comet orbit update uses for all comets at once); the hand-written Newton and bisection solvers are gone.
- Physical constants (G, c, the Boltzmann and Stefan-Boltzmann constants, solar, Earth and Jupiter figures, the AU, parsec and light-year) come from astropy. Solar mass is now 1.98841e30 kg (was 1.989e30) and solar luminosity 3.828e26 W (was 3.82e26), so derived star figures shift by up to 0.2%.
- Values that were JSON blobs now live in columns: a generation run's arguments are rows in `generation_run_arguments`, a work job's argv is rows in `work_job_args`, and the unused task `result` column is gone (schema v63, control schema v10). Adds indexes on the star-system quadrant and the star, planet and moon radius and star type columns.
- **Phase 0 is complete: the groundwork the new architecture builds on is in.** The package layout, SQLAlchemy models with Alembic migrations, values in columns instead of JSON blocks, Pydantic request models, one point-in-space object with stored velocity, orbital elements and saved sector paths, codec names from IDs, the RQ job queue with streamed logs, Shoelace and TanStack web components, one 3D map engine from the galaxy down to a moon, scipy and astropy for the physics, and a UX sweep. Version 8.0 starts the work built on top of it.
- Every random draw in generation goes through one wrapper (`planetgen/util/draw.py`) whose helpers are built on Python's `random()` alone and draw from the running unit's seeded stream, so a seed gives the same galaxy on any Python release, hash seed or locale, and changing the random source later means editing one file. This changes the numbers a given seed generates (GEN.56).
- The bright-star backfill after a `galaxy` run (GEN.98) goes out from the run's edge: the backfill distance past the farthest generated sector in every direction, not a radius around the starting sector. `--backfill-from` is now `edge` (the default) or `none`; the `requested` and `all` choices and the Generate page's "backfill from every generated sector" box are gone.
- Sector Map and an opened sector: rogue planets start hidden and the stars on; once rogue planets are shown each is ringed by default (MAP.137). A shown rogue planet is kept in the URL as `show=roguePlanet`.
- The Jobs list is minimal: job, status, started, duration, time left, progress and its controls. The Generate page's recent-jobs list drops its "By" column, is paged like every list, and points to the Work Queue for every job, including ones started from the command line.
- The orbit update (`python -m planetgen.cli.orbits`) moves and counts only objects that have moved far enough since their own last update: 0.01 mpc on a galactic orbit, 0.01 AU in a system, 100,000 km round a planet. Every moving object stores when its position holds and an indexed next-update-due time worked out from its speed, so a run looks only at what is due, and an object it skips catches up the whole time when it next comes due. Schema v66 (GEN.106).
- Retuned the phenomenon rates to the central observed values: neutron stars 0.4% and black holes 0.05% of stars, planetary nebulae 1e-7 per star, rogue planets 5.8 per star. Intermediate-mass black holes are now 0.1% of black holes (was 2%), accretion disks 0.1% (was 15%), 2% of neutron stars pulse (was 70%) and a tenth of those are millisecond pulsars (was 30%).
- `update.sh` now reloads Apache itself when it is running (and restarts it when the update just enabled a module), instead of printing the command; if Apache isn't running or the command fails it prints the command as before (OPS.8). macOS and Windows are unchanged.
- A route across sectors now searches only the systems near the straight line between its ends (a corridor that widens when the route has long hops) with A*, instead of rebuilding the nearest-neighbour graph over every placed system, so routing scales to large galaxies (NAV.10). `scripts/bench_nav.py` times it on a synthetic galaxy of any size.
- The orbital maths guards its edge cases (GEN.108): very eccentric comets past e = 0.999 move by the universal variable instead of Kepler's equation, a near-parabolic radius keeps its digits at perihelion, and an orbit update many orbits long drops the whole orbits before turning the phase. New building blocks for the n-body update: `kepler.universal_step` (any conic), modified equinoctial elements, the Roche limit and `classify_encounter`. Docs: orbital-updates.md section 7.
- The orbit update (`python -m planetgen.cli.orbits`) now also says how many objects entered or left a nebula and a supernova remnant (GEN.107), next to the objects moved and the sector changes it already reported.

### Removed
- `GET /api/systems/<id>/near` and `queryDb.systems_within_radius` (same-sector systems only, in light-years): `GET /api/near` replaces them.

### Fixed
- A finished Generate job no longer leaves its own queue for a moment after its lock is released.
- The databases-listing test no longer depends on which schema sorts first on the server.
- The two-factor tests no longer fail when a slow run crosses a 30-second code boundary.
- A Windows job whose runner was slow to start (more than 15 seconds) was called interrupted, because the process-creation check used the short grace period; it now uses the start limit.
- Releasing the job lock on Windows failed silently (the lock file was removed while still open), so a job that failed to start kept the lock.
- Tests: the octant test finds a point inside a warped nebula even when no metaball centre is inside it; the end-to-end K2V test reads the primary star, not whichever star the database returns first; the Windows no-Redis job test no longer swaps in POSIX process options.
- Galaxy Map: picking an arc of the galaxy, or the first slab inside it, now shows where it lies (distance from the core, bearing, height) and how many sectors it holds and has generated, as a picked block does (MAP.135).
- The Classes index lists the class categories in alphabetical order (UX.77).
- Maps: a binary star system picks as one system (its companion's dot picks the system, named for the system), not as two stars (MAP.136).
- 3D system view: a close pair's two orbits are each drawn around the barycentre at the star's own share of the separation, so a large star with a small companion no longer looks on a collision course (MAP.136).
- The Generate page no longer says "Connection lost; reconnecting..." every time the job log stream rolls over. The server ends each stream every 40 seconds on purpose and the browser resumes at the right line; the note now appears only if the stream stays down for 8 seconds, and clears when it is back.
- A visited secondary link-button (outlined, such as Navigate or Show on Galaxy Map) kept the filled button's dark text on a clear background, so it was dark on dark in dark mode. Secondary links now keep their own colors visited or not.
- The outlined button's hover no longer drops its text below readable contrast in light mode; it thickens its outline instead of shading the background.
- A browser test checks every button look (filled, outlined, active, pressed, danger, map toggles) at rest, hovered and focused at 4.5:1 text contrast, and an icon from the sprite at 3:1, with the OS in dark or light and the site's own theme choice following or overriding it.
- A test now checks that generating a sector around a backfilled bright star gives that star its planets (GEN.72 had already fixed this).
- The Galaxy Map no longer times out while a bright-star fill or a big sector run is writing. Every minute the web layer's tile cache saw "so much changed" (more than 1,000 sectors, or new bright stars) and threw every cached tile away, so every map request recomputed its tiles on a database that was already busy. Now a change that big is marked busy and the cache keeps serving its tiles for up to ten minutes before one refresh, and a freshness check that times out serves the cache instead of failing the request.

## [7.379.678] - 2026-10-09

### Changed
- The Galaxy Map's star class buttons and dimmest-star slider now also thin the stars drawn from the map's tiles at galaxy scale, and the Menu offers them even with no sector open (MAP.123).

## [7.378.678] - 2026-10-09

### Changed
- The hand-written `_migrate_vN_to_vM` steps (v8 to v61), their old-schema test fixtures and the tests of each step are removed; Alembic carries every migration from the v61 baseline. A database older than v61 is refused with `SchemaTooOldError` (run `update.sh` from an earlier checkout first). `PLANETGEN_MIGRATIONS_DIR` points the migrations elsewhere for tests (DB.11 complete).

## [7.377.678] - 2026-10-09

### Changed
- Every API request body and the Generate page's form fields are checked by Pydantic models with the old limits, and a bad request lists every problem field at once (an `errors` list in the API response, one joined message on the page) instead of stopping at the first. The hand-written checks are deleted (ADM.21).

## [7.376.678] - 2026-10-09

### Added
- The path of a body with no closed orbit through a sector (`physics/sector_path.py`): a test particle is followed from where it enters, with its velocity, against the sector's heavy masses until it leaves, and kept as a few cubic Hermite spline knots (two for a straight crossing, more where a mass bends it, at most 48). The exit is the next sector's entry. Saving it with the sector comes next.

## [7.375.678] - 2026-10-09

### Added
- `planetgen/db/models.py`: every table in `schema.sql` as SQLAlchemy `Table` objects, generated by `scripts/generate_db_models.py`, with a test that they match a database built from `schema.sql` (DB.11, third step).

## [7.374.678] - 2026-10-09

### Added
- Alembic now carries database migrations from schema v61 on (`src/planetgen/db/migrations/`, run by `planetgen.cli.migrate` after the older steps). A revision's id is its schema version, `SCHEMA_VERSION` is the newest revision's number, and each revision is still recorded in `schema_migrations`. New schema changes are revisions instead of `_migrate_vN_to_vM` steps; `migrations/README.md` lists the steps (DB.11, second step).

## [7.373.678] - 2026-10-09

### Added
- A planet, moon or comet can tell the orbit it is on right now from its position and velocity (`orbit_from_vector`: periapsis, eccentricity, inclination, node, argument of periapsis, true anomaly, semi-major axis, period) and give its projected course, the closed ellipse round its primary (`projected_orbit_au`). The orbit is never stored, so it follows the vector whenever the vector changes.

## [7.372.678] - 2026-10-09

### Changed
- Rogue planets, standalone black holes and neutron stars, nebulae, supernova remnants and their cores, quasars, interstellar comets and asteroid fields show a pronounceable name made from their ID and the galaxy's naming key where they showed a 19-digit hex ID. The ID stays stored and nothing is rewritten, so changing the key renames them at once (GEN.71, folds in GEN.73).

## [7.371.678] - 2026-10-09

### Added
- A to-do item (GEN.125) for stand-alone facilities to store a velocity like stars and systems do.

## [7.370.678] - 2026-10-09

### Added
- A star system's galactic velocity is stored (`star_systems.velocity_x_kms`, `velocity_y_kms`, `velocity_z_kms`, schema v61, worked out for systems saved before): the galaxy's rotation at its place plus, for a runaway or hypervelocity star, its speed along a random direction. The galactic orbit update turns it with the position, and a loaded sector gives it to the system's entry, its stars and its bodies. Run `sudo ./update.sh` to migrate.

## [7.369.678] - 2026-10-09

### Changed
- The database connection pool is now SQLAlchemy's `QueuePool` (5 kept, up to 10 open, a dead connection replaced when checked out) instead of DBUtils; `dbutils` is no longer a dependency. Id blocks, batching, statement timeouts, UTC sessions and the never-fail-over rule behave as before (DB.11, first step).

## [7.368.678] - 2026-10-09

### Fixed
- Static scripts and styles are linked with a fingerprint of the files themselves as well as the release (`?v=<release>-<fingerprint>`). The release version only changes when the post-merge stamp lands, so an update taken just before it (or an Apache that was not restarted) served new scripts under an old `?v=` that browsers keep for a year; modules from two versions then met, an import failed, and the Galaxy and Sector Maps stayed blank with dead controls while the System Map kept working. A changed file now always gets a new URL.

## [7.367.678] - 2026-10-09

### Added
- The galaxy has a naming key (8 hex digits) in the control database, drawn from the galaxy seed when `planetgen plan` makes a new seed and changeable by an admin on the Stats page or with `POST /api/admin/naming-key` (control schema v9, GEN.70). Run `update.sh` to add the table.

## [7.366.678] - 2026-10-09

### Added
- Every planet, moon and comet has a velocity relative to what it orbits (km/s, `velocity_x_kms`, `velocity_y_kms`, `velocity_z_kms`): set when its orbit is, kept in the database (schema v60, worked out for bodies saved before) and moved with its position by each orbital update. A comet's comes from its Kepler or Barker state (`comet_orbital_state` returns it).
- Stars, systems and exotic objects move at their galactic orbital speed along the galaxy's rotation, and the bodies in a system carry their star's velocity with their own. A sector loaded from the database knows its objects' velocities, sector address (ring, layer, slot) and the time the orbits were last advanced (`store.get_orbit_epoch_unix`).

## [7.365.678] - 2026-10-09

### Fixed
- No block or sector fill on the Galaxy Map is more opaque than 50% now, so what is inside it always shows, even when every sector is filled. Before, a filled sector could reach 85% and a fully filled block drew solid (Boss's lone one-star sector was a solid tan wedge).
- A generated sector's opacity follows what it holds (systems per sector on a log scale, up to 50% at about 600) instead of ranking against the other sectors in view, so a sector with one star is faint even when it is the only one on screen.
- When only one to three cells hold stars, the brightness and Color by scales put their mean in the middle and pad both ends (at least a decade of luminosity, 3 Gy of age, and so on), instead of stretching from the lowest to the highest of so few values; the legend shows the padded ends.

## [7.364.678] - 2026-10-09

### Added
- On a system's page, a star, planet, moon, asteroid belt or comet picked on the System Map (diagram or 3D) now has Start Here and End Here buttons that open the NAV page with that body as the start or destination. While a start or destination is being picked (`?pick=...`), the one button for that end returns to the NAV page with the course (NAV.50).

## [7.363.678] - 2026-10-09

### Added
- **A TODO item for asteroid-field avoidance.** NAV.51: courses go around asteroid fields (Boss, 2026-10-09).

## [7.362.678] - 2026-10-09

### Fixed
- The Galaxy Map showed no bright stars in the bulge, and none beyond the first 120 or so layers off the plane (MAP.133). Its galaxy-wide sample and each tile's picks ranked stars by luminosity alone, so the old giants of the bulge and the thick disk (no brighter than about 2,500 Lsun) never made the cut against young blue stars a hundred times brighter. Both now take an equal share from each population (young, intermediate, old, bulge), a population with fewer than its share leaving the rest to the others.
- Schema v59 replaces the `idx_bright_stars_off_plane` index with `idx_bright_stars_population`; run `sudo ./update.sh` to build it (reads every bright star once). The scattered stars themselves are unchanged.

## [7.361.678] - 2026-10-09

### Added
- The map Menu's "Highlight" buttons draw a kind of phenomenon larger and brighter on an opened sector, kept in the address as `mark=` (MAP.123).

### Changed
- A bookmark of a Galaxy Map view keeps what the map shows (`hide=`, `mark=`, `stars=`, `lum=`, `color=`), so opening it restores them (MAP.123).

## [7.360.678] - 2026-10-09

### Added
- Light-travel positions: `planetgen.physics.light_travel.apparent_position` gives where an observer sees an object, its position at the time the light left it (the retarded time, found by iteration), and `apparent_position_of` does the same for a `SpatialPosition3D` moving at its galactic velocity. A star 1 ly away at rest is seen a year ago where it is; one receding at 0.1 c is seen at 1/1.1 ly (VIEW.5).

## [7.359.678] - 2026-10-09

### Added
- One keep-out radius for every kind of object, returned as `keep_out` by `GET /api/objects/<ref>`: a Hill radius for planets and moons, the stored system perimeter for stars and systems, the galactic Hill radius (never below the object's own size) for black holes, neutron stars, quasars and rogue planets. Nebulae, remnants and asteroid fields have none: a course passes through them, with a note (NAV.24).

## [7.358.678] - 2026-10-09

### Added
- The Galaxy Map's course view draws the straight line from start to end apart from the route (dashed and muted, against the solid route), with a legend, the course's distance, bearing and mark, warp times and the stops as links above the map. Ringed stops on the map open their page when clicked (NAV.20).

## [7.357.678] - 2026-10-09

### Changed
- The web interface pass (UX.21) found nothing left to fix: with every page at 390, 600, 820 and 1280 px in light and dark, no buttons or labels overlap or run off their panel (the maps' open menus included), and every plain button changes something when clicked. Two new browser tests keep it that way: one clicks every button on every page, the other checks the controls inside each open map menu.

## [7.356.678] - 2026-10-09

### Added
- The map Menu can hide stars by spectral class and set the dimmest star shown (a luminosity floor); both ride in the address (`stars=`, `lum=`) on the Galaxy Map's opened sector and the sector page (MAP.123, the star classes and luminosity slider).

## [7.355.678] - 2026-10-09

### Added
- **A TODO item for each object's sector address.** GEN.124: ring, layer and slot as a coordinate that recalculates whenever the position changes.

## [7.354.678] - 2026-10-09

### Added
- `SpatialPosition3D` keeps a velocity in the galactic frame and relative to the nearest star (the system frame), the star's own velocity, and an `epoch_unix`: when the position and velocity hold. A body's speed relative to its star is read separately from the star's speed round the galaxy.
- `planetgen.physics.state_vectors`: `state_from_elements` and `elements_from_state` convert between a position and velocity and orbital elements for any ellipse, parabola or hyperbola, with `mean_anomaly_from_true` and `closed_orbit_points` for drawing the ellipse; `orbits.circular_orbital_velocity_au_per_year` gives the velocity on the circular orbits planets and moons follow.
- `SpatialPosition3D.sector_address`: the `(ring, layer, slot)` of the sector cell the position is in, worked out again whenever the position changes, with `set_sector_address` to move to another cell. Sector entries and every star, planet, moon and comet of their systems know theirs once the sector is placed in the galaxy.

### Changed
- `SpatialPosition3D.get_time_to_observable_movement` measures a system-scale or planetary move by the speed relative to the star, and a galactic one by the galactic speed.

## [7.353.675] - 2026-10-08

### Fixed
- The edit menu's browser tests no longer fail on the `?view=3d` the system page adds when it opens in 3D.

## [7.352.675] - 2026-10-08

### Added
- NAV takes any object: `/nav` and `GET /api/nav` accept a star, planet, moon, asteroid belt or comet as either end, written as an object reference (`planet:40`). A course to or from a body gets legs: out of its system to the heliopause, between the systems, into the destination body (all in the System Local Frame inside a system). Two objects in one system are one in-system leg, using the same warp and fold tables with a note (NAV.16).

### Changed
- `GET /api/nav` takes `from`/`to` as object references; the `from_kind`/`to_kind`/`from_type`/`to_type` parameters are gone (`from=nebula:2` replaces `from=2&from_kind=phenomenon&from_type=nebula`).

## [7.351.675] - 2026-10-08

### Changed
- The Galaxy Map opens as close as its zoom range allows with every charted sector in view, centered on them; Reset and zooming out still reach the whole galaxy (MAP.124).

## [7.350.675] - 2026-10-08

### Added
- **A Galaxy Map bug in the TODO.** MAP.133: after the bright-star scatter the map shows no bright stars in the bulge, because its sample ranks by luminosity alone.

## [7.349.675] - 2026-10-08

### Added
- **Three physics items in the TODO.** GEN.121 (a velocity and epoch on every object), GEN.122 (orbital elements kept in step with the state vector) and GEN.123 (the projected path through a sector as a spline), approved by Boss on 2026-10-08.

## [7.348.675] - 2026-10-08

### Added
- The Galaxy Map's Menu has a Color by choice: the default (age hue, density opacity, luminosity brightness) or one statistic (density, mean age, luminosity, star count) on a single ramp, with a legend that says what the colors mean. The choice rides in the address (`color=`) and every choice of a stage is ranked on one scale (MAP.131).

## [7.347.675] - 2026-10-08

### Fixed
- The planet-offset tests no longer assume the primary star is the system frame's origin: a close pair's origin is its barycenter and a wide pair's is its primary star, so a random system used to fail them about one run in four (TEST.108).
- The Sector Map nebula click test waits until the pointer, with the real mouse on the spot, still reads the nebula before it clicks, so a pick that settles after load no longer fails it (TEST.109).

## [7.346.675] - 2026-10-08

### Added
- One reference form for every object, `<kind>:<id>` (`sector:3`, `system:12`, `star:7`, `planet:40`, `moon:41`, `belt:5`, `comet:8`, `nebula:2`, ...; a bare number is a system). `GET /api/objects/<ref>` resolves one to its name, parent chain up to the galaxy, sibling references and its position in each frame (galaxy parsecs, sector-local light-years, system-local kilometers). `planetgen.galaxy.objectref` and `static/objectref.js` parse and print the form (NAV.7).

## [7.345.675] - 2026-10-08

### Added
- Picking a body in a 3D system view (the system page, or a system opened on the Galaxy Map) draws its orbit bright and full in the frame of what it goes round, a moon's around its planet, and fades the other orbits back (MAP.126).

## [7.344.675] - 2026-10-08

### Added
- On the Galaxy Map the wheel carries the zoom from a sector into a star's system (the star picked, else the one nearest the middle of the view) and back out again, and the address names the body picked (`&object=planet:12`) so a reload flies back to it (MAP.125, part 2).

## [7.343.675] - 2026-10-08

### Changed
- The system page opens in the 3D view by default; a visitor who last chose the diagram gets the diagram (MAP.74 follow-up).

## [7.342.675] - 2026-10-08

### Changed
- The Admin menu of a sector page and of a phenomenon page is now the last button of the action bar, as on a system page, instead of a separate button above the page content.

## [7.341.675] - 2026-10-08

### Added
- A star of an opened sector on the Galaxy Map can be opened in place: its system's stars, planets, moons, belts and comets are drawn where it sits, can be picked and flown to (down to a moon, followed as it orbits), and shown at true scale. The address keeps it (`&system=<id>`) and Up, Back and Forward step through it (MAP.125, first part).

## [7.340.675] - 2026-10-08

### Added
- A new test bug is on the list: the Sector Map nebula click test fails on main (TEST.109).

## [7.339.675] - 2026-10-08

### Added
- The system page has a 3D view beside the Diagram: stars with a glow, lit planets, moons and comets on 3D orbit lines, belts as particle rings, a free camera (turn, pan, zoom, fly keys, fly to and follow a body), true and compressed scale, and a time control. The address keeps the view and the selected body (`?view=3d&object=planet:12`), a list of bodies stands in for the canvas, and without WebGL the Diagram stays (MAP.71 to MAP.74, MAP.67's system and body forms).

## [7.338.675] - 2026-10-08

### Added
- A new test bug is on the list: a point-in-space test fails about one run in four (TEST.108).

## [7.337.675] - 2026-10-08

### Changed
- Search shows a card only for the groups that have results, and one line ("No sectors, systems, stars, planets or moons match.") when none do. The name you searched for is no longer echoed as a chip, and the header search box is hidden on the Search page.
- Home is the front door: the counts of sectors, systems, phenomena and standalone systems as links to their lists, plus the Galaxy Map, Navigate and Search. The Systems page no longer has a second Standalone Systems card; standalone is a filter on its list (the "standalone" badge opens it).

## [7.336.675] - 2026-10-08

### Added
- Positions of a system's bodies at any time, in Python (`physics/body_positions`) and its JavaScript twin, with the time clock the 3D system view will use (MAP.70). The scene endpoint now carries a pair's mass fraction and a comet's primary mass.

## [7.335.675] - 2026-10-08

### Added
- `GET /api/systems/<id>/scene` returns a system's stars, planets, moons, belts and comets with real orbit elements, colours and positions at the epoch the stored phases were last advanced, for the 3D system view (MAP.69).

## [7.334.675] - 2026-10-08

### Changed
- Documented how the stored position, mass and mu columns map to each object's `SpatialPosition3D`, with a save and load test covering stars, planets and moons (GEN.74 part 3).

## [7.333.675] - 2026-10-08

### Changed
- List and map cards no longer repeat the page title as a visible heading (Galaxy Map, Sector Map, System Map, NAV Map, All Systems, All Phenomena, Sectors): the heading stays for screen readers. "(3D)" is gone from the Galaxy Map's name.
- A single star's description and its spectral and luminosity class links, which the System list no longer repeats, now sit under the Stars table.

## [7.332.675] - 2026-10-08

### Changed
- Planets, moons, comets and stars now hold a `SpatialPosition3D` (in AU) with their mass and mu, anchored under their system, and the system moves them with it; stored position columns map to it (GEN.74 part 2b).

## [7.331.675] - 2026-10-08

### Changed
- A backfilled bright star gets a normal word-salad name when its sector is generated, and keeps the position ID it was known by as its unique ID (GEN.72). Before, the generated system was named by that 19-digit ID.

## [7.330.675] - 2026-10-08

### Changed
- The sector's Contents table column "Location" is now "Nearest" and lists only the three nearest neighbours; the sector name is no longer repeated on every row. The system page's header line reads "Nearest: ..." the same way (the breadcrumb already names the sector).
- The system page lists a single star only in the Stars table (not again as the first row of the System list), and the Stars table hides its Role column when there is only one star.

## [7.329.675] - 2026-10-08

### Changed
- A table that scrolls sideways shows a shadow on each edge with more to scroll to and keeps a thin visible scrollbar. On a phone the Quadrants and a system's Stars tables stack: the name is the row's title and every other column is a labelled line under it, so no column is cut off (UX.69).

### Fixed
- The layout test hides the page's own controls while it checks a header menu, since the menu drops over the page's action bar by design (TEST.106).

## [7.328.675] - 2026-10-08

### Fixed
- The layout test no longer opens the maps' popover menus (Menu, Bookmarks, history) when it opens every folded section; they drop over the page by design, like the header's menus, which it already checks one at a time (TEST.106).

## [7.327.675] - 2026-10-08

### Added
- A new test bug is on the list: the Sector Map's buttons overlap the Contents filters (TEST.106).

## [7.326.675] - 2026-10-08

### Changed
- The NAV result page is titled "Course: A → B" with no From/To chips, shows Optimal Route only when it has stops beyond the two ends, and is ordered summary, map, route, then the travel times in one collapsed section with a Warp/Dimensional Fold switch. Reverse course, Another destination and Show on Galaxy Map are a button row (UX.70).
- The NAV landing page is one "Plan a course" card: start from a map pick or a sector (UX.71).

### Fixed
- The links in the Admin hub list are underlined, so they no longer fail the link-in-text-block contrast check in either theme (TEST.105).

## [7.325.675] - 2026-10-08

### Changed
- The System card's Wikitext and Markdown buttons are now a quiet "View source" menu in the card header (UX.67).

## [7.324.675] - 2026-10-08

### Added
- A new item is on the list: the sector and phenomenon pages' Admin menus move into the shared action bar (UX.75).

## [7.323.675] - 2026-10-08

### Changed
- The system page's Admin menu moves into the action bar and lists only the system's own actions (upload to wiki, regenerate, change star, delete, place or remove a facility). Planet, moon and asteroid belt actions are in an Admin menu on each of their rows in the System list.

## [7.322.675] - 2026-10-08

### Changed
- Every system and phenomenon placed in a sector now holds one `SpatialPosition3D` (in light-years) with its mass and mu; `SpaceSector.place_in_galaxy` carries a sector's entries to their galactic place. `SpatialPosition3D` now takes a length unit and can carry its sector or star anchor (GEN.74 part 2a).

## [7.321.675] - 2026-10-08

### Added
- A new test bug is on the list: the Admin hub's links fail the contrast check (TEST.105).

## [7.320.675] - 2026-10-08

### Changed
- The request rate limits (Flask-Limiter) and the login lockouts now count on the Redis server (`redis.url`) by default, so every worker process shares them and a restart keeps them. `ratelimit.storage_uri` is empty by default (the Redis server); `memory://` still counts per process, and a Redis outage falls back to counting in memory (SEC.30).
- The failed-login counts moved out of the control database: its `login_throttle` table is dropped (control schema v8), so the counts start over once on the first `update.sh` after this change. `python -m planetgen.cli.lockouts` lists and lifts lockouts in Redis (SEC.30).

## [7.319.675] - 2026-10-08

### Added
- `planetgen.physics.position.SpatialPosition3D`: one body's position in the galactic, sector and system frames, each in Cartesian, cylindrical and spherical form, kept in step when any coordinate in any of them changes, with velocity, the next-due time from the design thresholds, and mass with its gravitational parameter (GEN.74 part 1; Boss's prototype `spacial-position.py` is removed).

## [7.318.675] - 2026-10-08

### Changed
- System Map, NAV Map and the phenomenon diagram: the caption paragraph is gone; the same text is behind a Map help button (UX.50). The diagram's Reset view is now Re-center.
- System Map: Measure distance sits in a button row under the map, with Map help (UX.61).
- Galaxy page: the map's breadcrumb replaces the page's own, so one trail reads Home, Galaxy, the levels, with the star at its end; on a phone it shows the current level and a labelled Steps button (UX.59).

## [7.317.675] - 2026-10-08

### Changed
- Galaxy Map and Sector Map: the how-to text under the map is gone; the gestures, keys and legend are behind a new Map help item in the Menu (UX.50, Galaxy and Sector maps).
- The toolbar's Reset is now Whole galaxy and the Menu's Reset view is Re-center (UX.57).
- Current moved into the round Steps menu as Jump to newest, and the Steps menu shows at every width (UX.58).
- The Slabs rail only appears while the stage has slab buttons (UX.60).
- On the Galaxy page only the map card goes wide; the title and crumb keep the page's usual left edge (UX.62).

## [7.316.665] - 2026-10-08

### Changed
- The sector header's chips hold facts only (edge, systems, stars, phenomena). The Quadrant link is plain text in the header line, and the interstellar comet estimate moves to a Details line under the map.

## [7.315.665] - 2026-10-08

### Fixed
- The System Map marks itself ready (`data-ready`) once its click handlers are on, and its browser test waits for that instead of a fixed time, so it no longer fails now and then under load (TEST.104).

## [7.314.665] - 2026-10-08

### Changed
- The gear menu now lists only Theme, Account, Admin and Logout (UX.72). The Admin page is the hub: it links Generate, Queue, Stats and Account, and every admin page carries a small tab row (Overview, Generate, Queue, Stats). The "Signed in" card is gone from the Admin page, and Log out and Account are in the gear.
- A sector's wiki link is set or cleared from "Set wiki link" in the Sector page's Admin menu, which starts from the current link and needs no sector ID (UX.73). The "Sector wiki link" form left the Admin page.

## [7.313.665] - 2026-10-08

### Fixed
- Choosing a NAV start or destination on the Galaxy Map and the Sector Map, every control now works as when browsing: the "Charted only" button was the last one locked on while picking. Only a choice holding something generated can still be taken (NAV.32).

## [7.312.665] - 2026-10-08

### Changed
- The Sector Map no longer tints the space around the sector: the amber and blue-grey block fills are gone and only the faint block outlines remain (MAP.61; Boss: sector space is not tinted).

## [7.311.665] - 2026-10-08

### Fixed
- A failed task's worker traceback rides with its exception without `add_note`, so it also works on Python 3.9 and 3.10.

## [7.310.665] - 2026-10-08

### Fixed
- The Galaxy Map's Menu panel no longer opens past the right edge of a narrow screen: below 700 px it hangs from the controls row instead of its button. The layout test also waits for the page to settle before it counts overlaps (TEST.103).

## [7.309.665] - 2026-10-08

### Added

- TODO item ADM.38: worker exceptions fail on Python 3.9 because `_picklable_error` uses `Exception.add_note`.

## [7.308.665] - 2026-10-08

### Added

- Every sector, star system, star, planet, moon, belt, comet and phenomenon
  has a unique ID (GEN.69, schema v58): a `uid` column on each table, written
  when a sector, system or phenomenon is saved. A sector's is its designation,
  so a sector nobody has generated already has one; a system's or
  phenomenon's is 96 bits and a star's, planet's, moon's, belt's or comet's 64
  bits (unique under its system), hashed from the galaxy seed, the parent's ID
  and the object's slot in it, so saving the same sector again gives the same
  IDs. An interstellar object or bright-sweep system keeps its position ID.
  Rows saved before this have none until their sector or system is saved
  again. The design is in `docs/design/object-ids.md`.

## [7.307.665] - 2026-10-08

### Changed
- Navigation buttons, links and the Wikitext/Markdown toggles are now outlined (secondary), so the primary action of a page or panel stands out. The bookmark is one icon-only ☆ / ★ toggle with the same look on the Galaxy Map breadcrumb, the map panels and the object pages.

## [7.306.665] - 2026-10-08

### Changed
- A system, a phenomenon and a sector page now share one action bar (Navigate, Show on Galaxy Map, Bookmark, in that order). A phenomenon gains Show on Galaxy Map; the sector's Show on Galaxy Map moves from the facts row into the bar.

## [7.305.665] - 2026-10-08

### Changed
- The Sector page's chip and the Search table's column now say "Edge" instead of "Cube edge" (the sector is an arc-shaped cell), and so does the sector's wiki text (UX.66).
- The admin Generate page shows its Current job card only while a job runs; with nothing running a "No job running" chip sits in the header row (UX.74).

## [7.304.664] - 2026-10-08

### Changed

- Choosing a NAV course on the Galaxy Map and the Sector Map stays on the map
  (NAV.29, NAV.33). An object's panel offers "Start Here" and "End Here"
  instead of "Nav from here" and "Nav to here". The first one pressed keeps
  the view and zoom where they are, puts the chosen end in a banner and in the
  URL (`?pick=to&from=system:12`), and asks for the other end; the user zooms
  out, in and across to find it with the same picker, and its button opens the
  NAV page with both ends. Cancel on the banner clears a pick begun on the map.
  The same buttons end a pick the NAV page began, and the system and
  phenomenon pages' pick buttons carry the same words. The sector scene JSON
  no longer carries NAV links: the page's script builds them (`navpick.js`).

## [7.303.664] - 2026-10-08

### Added

- TODO items TEST.103 (sector-page control overlap at 600 px) and TEST.104 (a flaky System Map browser test).

## [7.302.664] - 2026-10-08

### Changed

- Two-step sign-in codes are checked with `pyotp` and the setup QR code is
  drawn by `segno` (SEC.29). Secrets, the 30-second step, the one-step
  clock window and the no-reuse rule are unchanged, so enrolled admins keep
  working. The two hand-written modules (`totp.py`, `qrcode.py`, with a
  copy of Nayuki's QR generator) are deleted.
- The setup page's `otpauth://` link now reads `planetGen:<name>` where it
  read `planetGen%3A<name>` (the same label, written the usual way).

## [7.301.664] - 2026-10-08

### Fixed
- The population test that plants a species without a civilization no longer fails now and then: it cleared nothing from the first planet, which could already be a real species' homeworld (TEST.102).

## [7.300.664] - 2026-10-08

### Added

- TODO items UX.50 to UX.74: the 25 approved findings of the UX audit (UX.37), each a removal, merge or move that lets the interface get out of the way.

## [7.299.637] - 2026-10-08

### Added

- The admin API keys list, the duplicate-names list on the stats page, the
  job queue's Jobs list and every Search result panel are data tables like
  the others (UX.41): they scroll through every row and keep the 50-row pages
  for scripts-off visitors. API keys sort by label, dates and status and
  filter by status; Revoke keeps working from the scrolled table.
- `GET /api/search` takes `panels` (comma-separated panel names) to run just
  those result panels.
- Rows come from `/table/api-keys`, `/table/duplicate-names`,
  `/table/queue-jobs` (admins only) and `/table/search-<panel>`.

### Changed

- Times in the API keys, duplicate-names and Jobs tables read in UTC (they no
  longer switch to the viewer's time zone).
- Revoking an API key returns to the API keys list, not to a numbered page.

## [7.298.637] - 2026-10-08

### Added

- A sector's Contents and a Galaxy Map Quadrant's sectors are data tables
  like the others (UX.41): sort from the headers, filter Contents by Type and
  Octant and a Quadrant by Zone, and scroll through every row. Many rogue
  planets still read as one folded row; it links to the Contents filtered to
  Rogue Planet, which lists them one by one.
- The Contents "Show on map" buttons keep working as rows scroll in and out.
- `GET /table/sector-contents?sector=<id>` and
  `GET /table/galaxy-quadrant?quadrant=<I-IV>` serve their rows to the page.

### Changed

- The Contents table's folded rogue-planet group is no longer an expandable
  `<details>` row; see above.

## [7.297.637] - 2026-10-08

### Changed
- A Galaxy Map view now has a stated cap on the stars it can fetch (70,000, about 8 MB before compression), built from the per-tile caps, with a test that fails if a budget is raised past it; a full view of the test galaxy is also timed (MAP.109).

## [7.296.637] - 2026-10-08

### Changed
- The Galaxy Map lists fewer generated stars per tile the farther out you zoom (a table of budgets per level), and each sector may give only a share of a tile's budget, its brightest, so a dense filled region stays readable two or three zoom levels out from a sector while a sparse one keeps all its stars. Black holes, neutron stars and quasars follow a budget too, so they show only where it reaches them. Comets, rogue planets and asteroid fields stay on the Sector Map only (MAP.116, MAP.115).

## [7.295.637] - 2026-10-08

### Changed
- The galaxy-wide sample of bright stars the biggest Galaxy Map tiles draw from now holds the brightest stars off the plane as well as on it, found through a new index on `bright_stars` (schema v57, a virtual `off_plane` flag), so old giants above and below the plane show at the widest zooms (GEN.117).

## [7.294.637] - 2026-10-08

### Added

- The Species list, the Polities list and a polity's Systems are data tables
  like the others (UX.41): sort from the headers, filter species by
  spacefaring and era and polities by government and era, and scroll through
  every row. The old "All / Spacefaring / Not spacefaring" links are now the
  Spacefaring menu (`?spacefaring=yes|no`).
- `GET /api/species`, `GET /api/polities` and `GET /api/polities/<id>` take
  `sort`, `order`, filters and `facets=1`.

## [7.293.637] - 2026-10-08

### Fixed
- A tall Galaxy Map tile no longer fills its bright-star list with the plane's luminous stars alone: the list now reserves half its room for stars above and below the thin disk, so old giants in the thick disk and halo show up too (GEN.117).

## [7.292.637] - 2026-10-08

### Added

- The Sectors table (home page and `/sectors`), the All Systems table
  (`/systems`) and the Standalone Systems table are data tables like the
  Phenomena list (UX.41): click a header to sort, filter from the menus
  above (sectors by Quadrant; systems by where they are, single or binary,
  and octant) and scroll through every row, 50 fetched at a time. Each table
  on a page keeps its own sort, filters and place in the address bar
  (`sectors_sort`, `systems_sort`, `standalone_sort`, ...). Sectors are
  still nearest the core first until you sort them.
- `GET /api/sectors` and `GET /api/systems` take `sort`, `order`, filters and
  `facets=1`.

### Changed

- The Sectors table shows its density and distance in plain text with
  Unicode superscripts instead of HTML ones.

## [7.291.637] - 2026-10-08

### Added

- The Phenomena list is the site's first data table on TanStack Table and
  TanStack Virtual (UX.41). Click a column header to sort by it (click again
  to reverse), filter by type and by descriptor (the nebula class, remnant
  shape, rogue planet kind and so on) from the menus above the table, and
  scroll through every row: the next 50 arrive as they come into view, so
  the page never holds more than a few dozen rows. The address bar follows
  the sort and filters, so a reload or a shared link shows the same table.
  Without scripts the table still sorts and filters with plain links and a
  form, and pages with the usual pager.
- `GET /api/phenomena` takes `sort`, `order`, `type`, `descriptor`, `placed`
  and `facets=1` (the option counts for the filter menus).

## [7.290.637] - 2026-10-08

### Added
- A system inside a nebula gets a faint wash of the nebula's color over its System Map, with a line in the map's hint saying so (MAP.103, System Map part: the Galaxy Map, Sector Map and System Map now all show nebulae). A system's `inside` cloud also carries its `descriptor` (a nebula's type, a remnant's morphology).

## [7.289.637] - 2026-10-08

### Fixed
- Nebulae are drawn over space where no sector has been generated yet, as well as over filled sectors, in close views and far ones (MAP.104; the Galaxy Map draws every placed nebula from its shape since MAP.103, and a browser test now checks it over unfilled space).

## [7.288.637] - 2026-10-08

### Changed
- A nebula's page now shows the nebula itself: its irregular shape in 3D, with the galaxy's brightest stars dimmed around it for reference (drag to turn, scroll to zoom), instead of the flat AU diagram (MAP.105). New `GET /api/nebulae/<id>/surroundings` supplies the stars.

### Fixed
- The "-" button on a supernova remnant's AU diagram works from the first view: a remnant too big to open inside the 1 ly limit may zoom out to twice its opening view (UX.38).

## [7.287.637] - 2026-10-08

### Added

- A TODO item for an intermittent failure in `test_a_pass_removes_species_stored_without_a_civilization` (TEST.102).

## [7.286.637] - 2026-10-08

### Added
- The Sector Map (and a sector opened in place on the Galaxy Map) draws each nebula from its shape mesh once it arrives, in place of the plain sphere, so a nebula looks like the irregular cloud it is (MAP.103, Sector Map part).

## [7.285.637] - 2026-10-08

### Added
- A browser test that walks the Galaxy Map down to a sector, out with Escape, Home, Back and Forward and checks the breadcrumb steps match what a fresh load of the same URL shows. The out-of-sync breadcrumb MAP.106 reported no longer happens since the shared breadcrumb and URL work (NAV.14, MAP.67, MAP.66); this keeps it that way.

## [7.284.637] - 2026-10-08

### Fixed
- On a sector's map, Escape clears what is selected, and a click on the selected object clears it, so a nebula covering the whole sector can be unselected; the panel goes back to its opening text.

## [7.283.637] - 2026-10-08

### Added
- The Galaxy Map draws a nebula from its shape (the mesh served by `GET /galaxy/nebula/<id>/shape`) once it is big enough on screen, instead of only a flat sprite that faded out as you came close, so you can fly into one (MAP.103, Galaxy Map part).

## [7.282.637] - 2026-10-08

### Fixed

- A sector page (and the Galaxy page) no longer embeds the tile cache's
  count of cached tiles in its first-frame data, so the page's text stays the
  same between two loads until something changes. It started to differ when
  the sector page took the Galaxy Map's tiles (MAP.68), which failed
  `test_sector_page_and_tiles_show_a_phenomenon_the_cli_added`.

## [7.281.637] - 2026-10-08

### Changed
- A nebula holds a star system, phenomenon or smaller nebula only when the point lies inside its shape (GEN.75), not anywhere in its bounding sphere. Supernova remnants stay spheres. Containment is recomputed when a sector or nebula is generated; to refresh existing data, regenerate the affected sectors.

## [7.280.637] - 2026-10-08

### Added
- Each nebula now keeps a seeded irregular shape (schema v56: `nebulae.shape_*` columns and a `nebula_shape_balls` table), and `GET /api/nebulae/<id>/shape` serves it as a triangle mesh at a low-poly or full level of detail. Nebulae already saved get the same shape drawn from their own properties (GEN.75, part 2 of 3). **Run `sudo ./update.sh` to migrate.**

## [7.279.637] - 2026-10-08

### Added

- The Sector Map, on the sector page and in the Galaxy Map's opened sector,
  has a button per kind of object it draws (stars, nebulae, supernova
  remnants, asteroid fields, black holes, neutron stars, quasars, rogue
  planets, interstellar comets, neighboring sectors) in the Menu under "Show
  on the map" (MAP.79). All are on to begin with; turning one off removes
  those points, rings and labels, and a hidden kind can't be hovered or
  picked. The choice is kept in the URL (`hide=`), so a reload or a shared
  link keeps it, and showing something with "Show on map" or the screen-reader
  list brings its kind back. Rogue planets were already drawn faint by default.

## [7.278.637] - 2026-10-08

### Added
- `planetgen.galaxy.nebula_shape`: a nebula's shape from seeded metaballs in an anisotropic ellipsoid, warped by gradient noise and read at an isovalue, with a containment test and a marching-cubes mesh at a low-poly and a full level of detail. Nothing uses it yet; storage, the API and containment follow (GEN.75, part 1 of 3).

## [7.277.637] - 2026-10-08

### Changed

- The sector page's Sector Map is now the Galaxy Map's own engine locked to
  that sector (MAP.68, with MAP.67's one URL scheme): the same stars,
  clouds and bodies, hover, tooltip, ring, info panel and NAV links, with
  zoom buttons, Reset view, "Mark rogue planets" and the Contents table's
  "Show on map" buttons as before. It has no steps or history of its own, and a
  click on a neighboring sector opens that sector's page. Because it is
  the Galaxy Map, its scale line now reads in sectors and parsecs, the
  neighbors are the galaxy's own, and the browser's tile cache is shared
  with the Galaxy Map.
- The old Sector Map code (`sectormap.js`, `render_map_panel` and its
  no-script link list) is gone; `starmap.py` now only builds the scene data,
  `GET /sector/<id>/scene`. A sector with no place in the galaxy has no map.

## [7.276.637] - 2026-10-08

### Changed
- The Galaxy Map draws its stars measured from where the camera is looking instead of from the galaxy's center, so zoomed in they keep their exact positions rather than snapping by a fraction of a pixel far from the center (MAP.102, first part).

## [7.275.637] - 2026-10-08

### Changed
- Every page's breadcrumb is now the Galaxy Map's one-line breadcrumb (`static/breadcrumb.js`): when the steps don't fit, the middle ones fold into a "…" menu, at any width. The Galaxy Map uses the same component for its stages (NAV.14).

## [7.274.636] - 2026-10-08

### Changed

- On the Galaxy Map, clicking a generated sector no longer loads its
  page: the sector opens in place as the drill-down's last stage (MAP.66).
  Its stars, black holes, nebulae and other bodies are drawn where the
  sector is in the galaxy, the camera flies to fit it, and hovering or
  clicking one gives the same tooltip, ring, info panel, ☆ and NAV links
  as on the sector page. The breadcrumb, Back, Forward, Up and Reset view
  work as on any other stage, and the stage has its own URL,
  `/galaxy?sector=<designation>&open=1`. Clicking a neighboring sector
  steps sideways into it. A sector's panel keeps its "View sector →"
  link to the page.
- The Sector Map's scene is built by a new module (`static/sectorscene.js`)
  that both the sector page and the Galaxy Map use, and `starmap.py` hands
  the sector's place and size in the galaxy to the script. The scene is
  also available as JSON at `/sector/<id>/scene`.

## [7.273.636] - 2026-10-08

### Added
- A browser test that, choosing a NAV end on the Galaxy Map, clicks through the arcs, slabs and blocks down to a sector. The pick-mode wedge clicking NAV.46 reported already works since the shared picking layer (MAP.65, NAV.13, NAV.15); this keeps it that way.

## [7.272.636] - 2026-10-08

### Added
- Choosing a NAV start or destination now works on the system and phenomenon pages too: the same "Choosing a destination" banner with Cancel and bookmarks, and one "Use as start/destination" button in place of the Navigate buttons (NAV.15).

## [7.271.636] - 2026-10-08

### Added

- Planned NAV.50 (pick a planet, moon or belt as a NAV endpoint), split from NAV.15 because the NAV page only takes system and phenomenon endpoints until NAV.16.

## [7.270.636] - 2026-10-08

### Added
- `static/picker.js`, one shared selection (select, step out, step in, step sideways, one change event, and the trail for a breadcrumb). The Galaxy Map's Up and Reset buttons, Escape and Backspace and its breadcrumb now move it, and the map follows (NAV.13).

## [7.269.635] - 2026-10-08

### Added
- The Galaxy Map has a Current button beside Forward that jumps straight to the newest view in the map's history; it is disabled when you are already there, and works in the phone layout (MAP.95).

## [7.268.635] - 2026-10-08

### Changed
- Every page with admin actions now has one Admin menu (hidden from visitors) listing only that page's actions: regenerate, delete, change class or star, upload to a wiki, generate a neighborhood, and place or remove a facility. The inline admin panels and per-row Remove buttons are gone (ADM.34).

## [7.267.635] - 2026-10-08

### Fixed
- A failed Generate-page job no longer reloads away from its error. The log stays open with the error in view until you click **Continue** (ADM.24). At an interactive terminal a failed `planetgen` run waits for Enter before it exits, so the output isn't lost with the console window; runs with redirected input or output exit at once as before.

## [7.266.635] - 2026-10-08

### Fixed
- An unexpected error, a database error or a file error in a running action now prints its full traceback to the console and into the job's web log (and the debug log), not only a one-line message. A worker's own traceback comes back with its error. The job panel has a **Copy log** button that puts the whole output on the clipboard (ADM.25).

## [7.265.635] - 2026-10-08

### Fixed
- The bright-star backfill shows its progress bar on the Generate page from the moment it starts (unmeasured while it finds its sectors, then counting them with an ETA), and the "add a dimmer layer" run shows a bar for the phase that tops up already-backfilled sectors, which had none (ADM.26).

## [7.264.635] - 2026-10-08

### Changed
- The sector page's Admin panel (generate neighborhood, upload to the wiki) is one Admin menu button whose items open their forms in a dialog, like the Edit menu; the neighborhood's size-and-time check comes back as a dialog (UX.26).

## [7.263.635] - 2026-10-08

### Changed
- The Galaxy Map colors a generated sector and the blocks holding it by what its stars are, not by their average color (which came out red nearly everywhere). The fill is translucent: denser sectors are more solid (never fully opaque), the mean star age sets the hue (blue young, slate at the 4.5 Gy disk average, amber old), and the summed luminosity sets the brightness. A block takes the same three from its sectors, with ages weighted by star count.
- The `sector_stats` table keeps the raw numbers instead of a baked color: `mean_age_gy` and `total_luminosity_sol` replace `fill_share` and `color_r`/`color_g`/`color_b` (schema v55; the migration works the new two out for sectors already generated). `GET /api/galaxy/stage` returns `stats` (systems, expected systems, stars, mean age, luminosity) in place of `look`.

## [7.262.635] - 2026-10-08

### Fixed
- The parallel-galaxy interrupt test starts its run with Ctrl+C at its default, since a pytest worker can leave it ignored and an ignored Ctrl+C is inherited, so the run finished with status 0 in the full suite (TEST.101).

## [7.261.635] - 2026-10-08

### Added
- The ID-to-words codec (`gatedPhonemeCodec.py`, which sat in the repo root and was used by nothing) is now `planetgen.names.gated_phoneme_codec`, with tests that pin its words. It needs only the standard library, and its decoding is exact when the identifier's length is given (GEN.120).

## [7.260.635] - 2026-10-08

### Changed

- **The Generate page's job log streams live in a terminal (ADM.22).** A running job's output and state now arrive over Server-Sent Events (`/admin/generate/jobs/<id>/stream`) instead of a poll every two seconds, and the output shows in an Xterm.js terminal (vendored under `static/vendor/xterm/`), so progress redraws and colours look as in a console. The stream ends every 40 seconds and the browser reconnects, resuming from the byte offset it reached, so a dropped connection loses and repeats nothing; the progress bar and step update from the same stream. Every job also has a **Download log** link (`/admin/generate/jobs/<id>/log`) with the full output. Without EventSource the page polls as before. The two pages that show the terminal (`/admin/generate`, `/admin/generate/jobs/<id>`) allow inline styles, which Xterm.js writes; their scripts and every other page keep the strict Content-Security-Policy. The deployment notes say why the app server must be threaded.

## [7.259.635] - 2026-10-08

### Fixed
- The lunar-spacing and class-change tests set their own random seed and retry until the generated body suits them, so their result no longer depends on which tests ran before (TEST.95, TEST.97).

## [7.258.635] - 2026-10-08

### Fixed
- The job tests treat a job as finished only once its runner has written a final status and released the job lock, so starting the next job right after one no longer meets a busy lock under load (TEST.94).

## [7.257.635] - 2026-10-08

### Fixed
- The parallel-galaxy interrupt test generates bigger sectors so a loaded machine can no longer finish the run before the test's signal arrives (TEST.101).

## [7.256.635] - 2026-10-08

### Fixed
- The Galaxy Map's browser tests wait until the map has drawn the stage its breadcrumb names (new `galaxyReady` on the canvas) instead of a clock, so the camera, scale line, slab button and hover tests no longer fail under load (TEST.96, TEST.98, TEST.99, TEST.100).

## [7.255.635] - 2026-10-08

### Added

- Planned a follow-up to the Shoelace migration: form fields as Shoelace components (UX.49), and closed UX.40 and ADM.14 in the plan.

## [7.254.635] - 2026-10-08

### Added

- **A progress bar on the database reset and the orbital update.** `python3 -m planetgen.cli.reset` shows the table being wiped and how many of them are done, and the New galaxy and Reset actions on the Generate page show the same bar in the job's status. `python3 -m planetgen.cli.orbits` shows its four steps (orbital phases, comets and facilities, galactic orbits, containment/nearest systems/locations) with the tables or sectors of the current step beneath, and writes the same progress for any job that runs it.

### Fixed

- **Resetting a big database is much faster.** The reset used to count every row of every table (`COUNT(*)`) before it wiped anything, which reads nearly the whole database from disk; it now takes the storage engine's own estimates, so its row counts read "about N rows" and the reset goes straight to the wipes.

## [7.253.635] - 2026-10-08

### Fixed

- **The Generate pages' text boxes line up (ADM.14).** On the admin Generate page (every form: New galaxy, Generate sectors and its modes, Plan, the bright-star and layer forms) and on the one-off system page, each group of fields is a grid of equal columns whose fields share three rows (label, box, hint), so the boxes have one width, sit level with each other and share left edges however long or wrapped a heading is, at phone and desktop widths in both themes.

## [7.252.634] - 2026-10-08

### Fixed
- With "Charted only" on, the Galaxy Map's blocks, slabs and wedges holding nothing generated can be hovered and picked just as with it off; only choosing a NAV start or destination still needs something generated (MAP.112).

## [7.251.634] - 2026-10-08

### Fixed

- **The system page's buttons sit on one row (UX.27).** Navigate from here, Navigate to here, Show on Galaxy Map and Bookmark no longer wrap or overlap: where the four don't fit, the two navigate buttons fold into one Navigate menu (From here, To here), and on a phone Show on Galaxy Map moves into that menu too. The layout follows the row's own width (container queries), not the device.

## [7.250.634] - 2026-10-08

### Fixed
- Picking on the Galaxy Map while the stage just picked was still loading its data no longer says "There is no layer x here" and leaves the breadcrumb ahead of the drawing: a choice is now taken on the stage that is drawn (MAP.107).
- A test now checks that the slab buttons never cover any part of the Galaxy Map at any width (MAP.108, the buttons-block-clicks part).

## [7.249.634] - 2026-10-08

### Added

- TODO item GEN.120: move the gated phoneme codec from the repo root into the naming package (Boss, 2026-10-08).

## [7.248.634] - 2026-10-08

### Changed

- **Admin edit actions are one Edit button that opens a menu (UX.26, UX.31).** On the system, sector and phenomenon pages the inline Regenerate, Delete, Change star and Change class folds are replaced by an Edit menu per row; picking an action opens its confirm step in a dialog on top of the page (Escape or Cancel closes it, the primary button submits), with the menu and dialog working from the keyboard. The posted forms and the server's `edit_action` handling are unchanged.

## [7.247.634] - 2026-10-08

### Added

- TODO items TEST.98 to TEST.101: four load-only test flakes seen in the full suite run for PR #531.

## [7.246.629] - 2026-10-08

### Fixed
- The Galaxy Map's "Generated only" toggle is now "Charted only": besides dimming blocks with nothing generated it dims the stars outside charted sectors and outlines each charted block, so the charted space stands out in every wedge (MAP.111).

## [7.245.629] - 2026-10-08

### Changed

- The Galaxy Map and the Sector Map now pick, hover and fill their info
  panels through one shared module (`static/mappick.js`, MAP.65). The
  Sector Map gains a hover tooltip and ring, and a ☆ Bookmark button in
  its info panel for a system, a phenomenon or a generated neighboring
  sector. On the Galaxy Map, black holes, neutron stars, quasars and
  nebulae gain a hover tooltip and ring, "Nav from here" and "Nav to
  here" links (or the pick button while a NAV end is being picked), and a
  ☆; a generated sector's panel gains a ☆ too.

## [7.244.629] - 2026-10-08

### Added

- **Shoelace components for buttons and menus (UX.40, first step).** The site now ships a subset of [Shoelace](https://shoelace.style) 2.20.1 as vendored ES modules under `static/vendor/shoelace/` (no CDN, no bundler, same-origin files under the existing Content-Security-Policy), registered on every page by `static/components.js` and coloured from the site's own light and dark tokens (`static/shoelace-theme.css`). `scripts/vendor_shoelace.py` rebuilds the vendored folder from the npm package.

### Changed

- **The header's Menu and settings gear are Shoelace dropdowns (UX.2).** They close on an outside click or Escape and give focus back to their button, replacing the `<details>` script in `theme.js`; each panel is as wide as its longest entry (at most 22 rem, never past the screen) instead of a fixed 14 rem minimum.

## [7.243.629] - 2026-10-08

### Added

- TODO items DB.14 and MAP.128 to MAP.132: how sectors and blocks are colored (Boss, 2026-10-08).

## [7.242.629] - 2026-10-08

### Added

- TODO item TEST.97: a load-only test flake seen in the full suite run for PR #525.

## [7.241.629] - 2026-10-08

### Fixed

- The Generate page's jobs folder defaults to `/var/lib/planetGen/jobs`, in the checkout's own folder, instead of `/var/lib/planetgen/jobs`, which differed from it only by case and left two folders; `update.sh` moves jobs from the old folder into the new one (a running job stays put) and removes the old folder once empty. A `jobs.dir` set in `config.json` is left alone (OPS.19).

## [7.240.629] - 2026-10-08

### Fixed

- Galaxy Map: the slab buttons beside the map, and their leader lines,
  are always in slab-number order (MAP.110). The order runs the same way
  as the slabs on screen (the highest slab first when the stack's top
  shows above its bottom, the lowest first when the view has turned
  under the plane), and each line ends on its slab's outline no higher
  than the line above it, so the lines still don't cross.

## [7.239.622] - 2026-10-08

### Added

- TODO items TEST.95 and TEST.96: two load-only test flakes seen in the full suite run for PR #522.

## [7.238.622] - 2026-10-08

### Fixed

- `GET /api/databases` counts each schema's sectors and systems over the one connection it lists them with, instead of opening a connection pool per schema and keeping it; under the parallel test suite that ran MariaDB out of connections (TEST.93).
- A Generate-page job in its first seconds shows as starting, not interrupted, while its runner process hasn't yet started running the job (TEST.92).

## [7.237.622] - 2026-10-08

### Changed

- Wiki uploads (`POST /api/systems/<id>/wiki`, `POST /api/sectors/<id>/wiki`) now run on the Redis queue and wait up to eight seconds, so a quick upload still answers `201` with the page; a slower one answers `202` with a job id (PERF.24, step 4e). The worker reads the wiki's login details from the environment and `config.json` itself, so they never pass through Redis.

## [7.236.620] - 2026-10-08

### Added

- TODO item TEST.94: test_old_jobs_are_pruned raises JobBusy again because the job lock outlives the finished job (phase 0, Bugfixes: ops and flakes).

## [7.235.620] - 2026-10-08

### Fixed

- The Generate page's prevalence fields start at each feature's usual share (habitable worlds 24.2% of systems, asteroid belts 59%, and so on, each naming what it is a share of) instead of a 0% change, and take the share wanted; the run gets the percentage that moves the usual share there (ADM.37).

## [7.234.620] - 2026-10-08

### Changed

- The Generate page's size and time estimate and the one-off system page's generator now run on the Redis queue and the page waits for them (PERF.24, step 4d). Without a Redis server (Windows without WSL's Redis) they run in the web process as before.

## [7.233.620] - 2026-10-08

### Added

- TODO item ADM.37: the Generate page's prevalence fields should show each feature's real default share instead of "0% change" (phase 0, Bugfixes: prevalence).

## [7.232.618] - 2026-10-08

### Changed

- Creating a system, regenerating a system, a planet, a moon, a belt or a phenomenon, changing a planet's or moon's class, and changing a system's star now run on the Redis queue (PERF.24, step 4c). The request waits up to eight seconds for the job, so a quick edit still answers in the same response; a slower one answers `202` with a job id for `GET /api/jobs/<id>`. Deletes and plain renames still happen in the request. Refusals (404, 409, 400) come back unchanged, and `503` means no Redis server answered.
- A queued job that the work refuses (an `ApiError`) reports `error_status` with its HTTP status in `GET /api/jobs/<id>`.

## [7.231.618] - 2026-10-08

### Fixed

- Log lines printed while a progress bar is up are no longer hard-wrapped at the console's width (80 columns in a web job's log), and text in square brackets is printed as written instead of being read as `rich` markup, which could drop or crash a line (ADM.23).

## [7.230.618] - 2026-10-07

### Fixed

- The Generate page's "New galaxy" and "Generate sectors" forms (every mode, "around a sector" included) have a folded Prevalence section: one percentage per feature, blank or 0 for the usual chance, passed to the run as `--prevalence FEATURE=PERCENT` (ADM.16).

## [7.229.618] - 2026-10-07

### Added

- TODO item TEST.93: timing tests fail and MariaDB drops connections under full-suite load (phase 0, Bugfixes: ops and flakes).

## [7.228.618] - 2026-10-07

### Changed

- `POST /api/sectors/<id>/regenerate` now queues the delete-and-regenerate on Redis and answers `202` with a job id, instead of holding a web worker (PERF.24, step 4b). Unknown (404) and off-grid or unplanned (409) sectors still answer before anything is deleted. The admin page's Regenerate button says the job is queued.

## [7.227.618] - 2026-10-07

### Fixed

- Sector and galaxy runs can make a feature more or less common instead of forcing it on every system: `--prevalence comets=+50` gives 1.5 times as many systems with comets, `--prevalence habitable_world=-100` none. It works for every former forcing option (habitable world, asteroid belt, comets, large star, moons, max planets, intelligent life, binary, wide binary, planets), and each system stores the setting with its config (database schema v54).
- `-large_star` on a single system now really keeps out a large star; it used to make no difference.

## [7.226.618] - 2026-10-07

### Changed

- `POST /api/sectors/<id>/generate-neighborhood` now queues the run on Redis and answers `202` with a job id, instead of holding a web worker for the whole run (PERF.24, step 4a). The math check, the unknown-sector and missing-skeleton errors (404, 409) and the disk-space refusal (507) still answer before anything is queued; `estimate_only` still answers at once.
- Added `GET /api/jobs/<id>`, the state, result and error of a queued API job.

## [7.225.618] - 2026-10-07

### Fixed
- **The galaxy had almost no bulge, so edge-on views showed a flat disk
  (GEN.118), and no thick disk (GEN.119).** The density model's bulge was
  a 200 pc sphere holding 0.6% of the stars; the Milky Way's holds about
  31%. The model is now three published Milky Way components:
  - a thin disk (scale length 2.6 kpc, scale height 300 pc) carrying the
    spiral arms. `--disk-scale-height-pc` is now the exponential scale
    height surveys quote, so the disk is twice as thick as before;
  - a thick disk (2.0 kpc, 900 pc, 4% of the plane's density at the Sun);
  - a boxy bar bulge, Dwek et al. 1995's fit to the COBE/DIRBE image,
    angled 27 degrees from the Sun-center line.

  New `planetgen plan` defaults: `--disk-scale-length-pc 2600`,
  `--disk-scale-height-pc 300`, `--bulge-scale-radius-pc 1580` (along the
  bar), `--bulge-amplitude 3.11`, and the calibration point at 8.2 kpc
  (3.15 scale lengths). The young and intermediate populations' scale
  heights are now 50 and 150 pc. The Galaxy Map's prisms draw the same
  model. `tests/test_galaxy_milky_way.py` and two new math checks compare
  it with Bland-Hawthorn & Gerhard 2016.

  At the default shape a full galaxy has about 23.9 billion candidate
  sectors (was 12.3), 21.0 billion qualifying (was 10.3) and 299 billion
  systems (was 88), reaching 4.1 kpc above the plane (was 1.3). An
  existing galaxy keeps its stored shape numbers but the new formulas read
  them differently, so plan and generate it again.
- **The Galaxy Map's slab lines could cross after the view turned.** The
  buttons were ordered by each slab's middle but each line ends where the
  slab's outline comes nearest its button; when those ends come out of
  order the buttons now follow them.

## [7.224.618] - 2026-10-07

### Changed
- **The admin Generate page's jobs run on Redis with RQ (PERF.24, step 3).** Starting a job queues it on a queue of its own and starts one burst worker for it, detached from the web server; the worker exits when the job is done. The job's files (`state.json`, `output.log`, progress, the Cancel file) and the one-job lock work as before. Without a Redis server at `redis.url` no job starts and the page says so, except on Windows, where Redis runs in WSL: there the job runs in its own process as before. `planetgen.cli.job` is now `planetgen.web.job_runner`.

## [7.223.617] - 2026-10-07

### Added

- TODO item TEST.92: the web job runner's first-failure test sometimes reports the job as interrupted under full-suite load (phase 0, Bugfixes: ops and flakes).

## [7.222.617] - 2026-10-07

### Changed
- **Parallel generation runs on Redis with RQ (PERF.24, step 2).** A run with more than one worker queues its sectors, layer scatters and backfill blocks as RQ jobs on a queue of its own and starts that many burst workers for it (`python -m planetgen.cli.worker`), which exit when it's done. Results are the same for the same seed at any worker count, as before. A task whose worker dies now runs once more on a fresh worker, and the run fails with `WorkerDied` only if that worker dies too. Without a Redis server at `redis.url`, a run says so and generates one sector at a time. The control database's lease still keeps two runs' workers off the machine at once.

## [7.221.616] - 2026-10-07

### Added
- **The work-queue audit (PERF.19).** `docs/design/work-queue-audit.md` lists every path where the API or website generates or writes something, where each runs today (a web job, the generation work queue, or inside the request), and where each goes once the queue moves to Redis and RQ (PERF.24). Long generation leaves the request, short generation is queued and waited on, and plain row writes and logins stay in the request.

## [7.220.616] - 2026-10-07

### Fixed

- The System Map browser test no longer fails now and then: it clicks each body on its dot instead of the middle of its marker, which can be empty map when the body's label sits to one side, and it measures from the star to the planet with moons, so a system with only one planet still works.
- The Galaxy Map drill-down browser test no longer clicks the block it is already in when, under load, the map's tooltip still names it.

## [7.219.616] - 2026-10-07

### Fixed
- The Windows installer no longer stops when planetGen isn't installed yet: the check for an existing editable install asks Python where the package is without importing it, so it prints no traceback.

## [7.218.616] - 2026-10-07

### Fixed

- Every rogue planet in a sector's Contents table, alone or in the rogue planet group, shows a small map icon beside its name in place of the "Show on map" button. The icon has a tooltip and a screen reader label, and it selects the planet on the Sector Map as the button did.
- The site has one icon sprite (`static/icons.svg`) with the approved set: edit, delete, show on map, filter, navigate, menu, settings, back, forward, up, reset and bookmark. Pages use it through `icon(name)`.

## [7.217.616] - 2026-10-07

### Added

- TODO item TEST.91: the System Map drill-and-measure browser test sometimes cannot click the first moon (phase 0, Bugfixes: ops and flakes).

## [7.216.615] - 2026-10-07

### Fixed

- A sector's Contents table lists its star systems first, then its other phenomena and facilities, then its rogue planets, each group nearest the center first. The rogue planet group opens as a row the table's full width, and each rogue planet in it shows its octant and location.

## [7.215.615] - 2026-10-07

### Fixed

- Every comet on a system page links its class, single-apparition (parabolic) comets included: they have a class page of their own now (UX.29).

## [7.214.615] - 2026-10-07

### Added

- TODO items GEN.118 (the galaxy bulge is about 40 times too light) and GEN.119 (no thick disk) in phase 0, Bugfixes: generation.

## [7.213.615] - 2026-10-07

### Fixed

- A Generate page job whose runner is slow to start (a busy server) is no longer shown as interrupted while the runner is still alive, which could let a new job start and then find the old one running again (TEST.90).
- A test now checks that bright stars spread far above the galactic plane, by population (GEN.117: the thin band seen on a server was the code from before 7.194.608).
- Several intermittently failing tests no longer depend on random draws or timing (TEST.82, TEST.84, TEST.86, TEST.88, TEST.89, TEST.71).

## [7.212.613] - 2026-10-07

### Added
- TODO: GEN.117 (bug): Boss's bright-star sweep still shows every bright star in a thin band on the galactic plane after the GEN.79 fix.

## [7.211.612] - 2026-10-07

### Changed
- **The pages' in-memory API cache runs on cachetools (PERF.25).** `planetgen.web.lib.pagecache` keeps its entries in a `cachetools.TTLCache` sized by body length, with least recently used out first and nothing served past `max_age_seconds`. That replaces its hand-rolled ordered dict, byte count and age checks. The settings, the stamp check, clearing on writes and what gets cached are unchanged. The Galaxy Map's tile cache (`tilecache.py`) stays as it is.

## [7.210.612] - 2026-10-07

### Added
- TODO: GEN.116 (bug) keeps watch for Boss's error when generating a neighbourhood near the galaxy edge, until his error text arrives. TEST.90 (bug) records a JobBusy flake in the old-job pruning test under parallel runs.

## [7.209.610] - 2026-10-07

### Fixed

- When the database disk has no room for a run, the Generate page and a sector's "Generate more sectors around this one" now offer "Generate anyway" instead of only refusing (ADM.33). Choosing it starts the job and records the override in the activity log (`job.generate_anyway`).

## [7.208.610] - 2026-10-07

### Changed
- **Redis on Windows means Redis in WSL2 (OPS.27).** The Windows guide, `config.md` and the warning `install.ps1`/`update.ps1` give when no Redis answers now point only at Redis in WSL2, with the setup steps (`sudo apt install redis-server`, `systemctl enable --now redis-server`, the default `redis.url` through WSL2's localhost forwarding, and keeping WSL running). Memurai is no longer suggested.

## [7.207.610] - 2026-10-07

### Changed
- **planetGen is installed as an editable package, and the package move is done (OPS.24, step 14, second half).** install and update now install the checkout with `pip install -e` (no dependencies; they still come from the lock), so `planetgen` imports from anywhere and a `git pull` still takes effect without a reinstall. The web app, the job runner and the tests no longer add `src/` to `sys.path` or `PYTHONPATH`, every command-line tool runs from any directory (`python3 -m planetgen.cli.<name>`, on Windows with the venv's Python), and the `planetgen` command is pip's own console script (linked into `/usr/local/bin` from the venv on macOS). The wiki upload client moved from `src/wikiClient` to `planetgen.wiki`, the last module outside the package. Run `update.sh` (or `update.ps1`) once after pulling this so the editable install is made.

## [7.206.610] - 2026-10-07

### Changed
- **`generate.py` is split into `planetgen.generation.run_*` and `planetgen.cli.generate` (OPS.24, step 14, first half).** Each command's work is in its own module (`run_system`, `run_sector`, `run_galaxy`, `run_plan`, `run_phenomenon`, `run_population`, with what they share in `run_common`), and the command line is `planetgen.cli.generate`. The root `generate.py` is gone: on Linux and macOS run `planetgen <command>` (the installer's `planetgen` command now runs `planetgen.cli.generate` from the checkout); on Windows run `python -m planetgen.cli.generate <command>` from the checkout's `src\`. The Generate page, the one-off system page, the job queue's retries and the install and update scripts run the generator the same way. Nothing a command does changes.

## [7.205.610] - 2026-10-07

### Fixed

- The size estimate before a bulk run is closer to what the run adds (PERF.26). Fills measured on MariaDB add 37 to 52 KB per star system, so the default before anything is measured is now 48 KB instead of 70 KB. The size measured after a run leaves out the galaxy-wide tables a plan or the bright-star scatter writes, and on MySQL 8 it reads the tables' current sizes instead of a copy cached for a day.

## [7.204.610] - 2026-10-07

### Changed
- **`stellarObjects/utils.py` is split into the `planetgen` packages, and `stellarObjects` is gone (OPS.24, step 13 of 14).** The number, distance, speed, duration, temperature, pressure and age text is `planetgen.util.format`; `finite_domain` is `planetgen.util.checks`; the power-law, bounded-bell and log-uniform draws are `planetgen.util.random`; unit conversions are `planetgen.physics.units`; the orbital helpers are `planetgen.physics.orbits`; the disk physics is `planetgen.physics.formation`; the galactic orbit is `planetgen.galaxy.galactic_orbit`; and the word-salad name generator is `planetgen.names.wordsalad`. The few helpers with one caller went to that caller's module. `planetgen.web.lib.fmt` now wraps `planetgen.util.format` instead of keeping its own copies, the Galaxy pages' Quadrant and Zone come from `planetgen.galaxy.geometry`, and the copies of the log-uniform draw are one `log_uniform`. Every draw is unchanged, so a seeded galaxy generates the same.

## [7.203.610] - 2026-10-07

### Fixed

- The console no longer tells the user no (GEN.81). A run the database disk is too small for, a ring, shell or block past 2,000 sectors without `--limit` or `--yes`, an address named outside the galaxy's outline, `--min-habitable` above a sector's drawn count, and a forced body no system had room for each print a warning, and the run goes ahead: the named address is generated at the halo density, the sector gets that many systems, the last system tried is kept. Bad arguments are still refused. `--strict` (on `system`, `sector` and `galaxy`) keeps the old stop for scripts.
- The installer no longer mixes apt's NumPy 1.x builds with pip's NumPy 2 (OPS.26). On Ubuntu 24.04 pip put NumPy 2.5.3 into /usr/local as a dependency of a newer scipy, and apt's astropy (erfa) and scikit-image, built for NumPy 1.x, stopped importing. NumPy, scipy, astropy with pyerfa and scikit-image now always come from one place: when pip must provide any of them, or would pull any of them in, it installs all of them from `requirements.lock`, and the report says so.

## [7.202.610] - 2026-10-07

### Changed
- TODO: Boss's answers to the open decisions are recorded. Wide-binary planets are named from two words (GEN.71). Orbital updates follow real time, one day per day, with an option to advance more (GEN.105). The orbital neighbour search is checked against the ring, layer and slot sectors (GEN.109). GEN.65 becomes a test across neighbourhood centres. New OPS.27 points the Windows installer at Redis in WSL.

## [7.201.610] - 2026-10-07

### Changed
- **The HTML pages and the app factory move into `planetgen.web` (OPS.24, step 12 of 14).** `src/html/web/`, with its templates, is now `src/planetgen/web/`, and `create_app` moved from `planetgen.api.app` to `planetgen.web.app`. `src/html/` now holds only `wsgi.py` and `static/`, so the Apache, gunicorn and waitress setup doesn't change. `wsgi.py` adds only `src/` to `sys.path`, and the tests no longer need `src/html` on it.

## [7.200.610] - 2026-10-07

### Added
- TODO: OPS.26 (bug) records Boss's Ubuntu 24.04 install failure. apt's NumPy-1 builds of astropy, erfa and scikit-image cannot load alongside the NumPy 2 that pip pulls in, so the requirements probe reports them as unusable.

## [7.199.608] - 2026-10-07

### Changed

- `update.sh` and `update.ps1` no longer ask about or run the population pass, so an update that just wiped the database doesn't offer to fill it; the closing message says how to run `generate.py population` by hand. `POPULATION=1` and `-Population` are gone from the update; the installers keep their prompt (OPS.7).

## [7.198.608] - 2026-10-07

### Fixed

- Signing in no longer ends on a "form expired" error when the sign-in (or two-step code) form is sent a second time after the first one already signed in; the second send goes back to the login page, which forwards the signed-in admin onward (SEC.31).

## [7.197.608] - 2026-10-07

### Fixed

- Search: a Phenomenon Class tag only narrows its own phenomenon type, so picking Black Hole and Nebula with a nebula class keeps the black holes; tags in one group combine with OR and groups with AND (UX.44).

## [7.196.608] - 2026-10-07

### Fixed

- The Galaxy Map no longer offers a planned sector outside the galaxy's stored outline (such as a layer above the top), which generation would refuse to fill (MAP.118).

## [7.195.608] - 2026-10-07

### Fixed

- Admin scripts reject a `--mysql-port` outside 1 to 65535 with a usage error before connecting, instead of a connection error later (OPS.6).
- Whole numbers stay in plain digits up to 999,999 and go scientific from 7 digits; numbers shown with decimals still go scientific from 5 whole digits (UX.36).
- The System Map side panel shows a planet's or moon's surface pressure under its surface temperature (MAP.117).
- The NAV page's course map spans the page's width, and its stop names, bearing and scale labels stay at body text size on phones and desktops; hop names that would overlap are left out (the route list below names every stop) (NAV.41).

## [7.194.608] - 2026-10-07

### Fixed

- Every sector inside the galaxy's outline now has some chance of a star and of a bright star (GEN.78). The density model has a halo floor (`tuning.MIN_RELATIVE_DENSITY`, a thousandth of the local density, made of old stars), and the bright-star scatter and backfill no longer skip cells that expect under one star.
- Bright stars follow the density model in every layer, including a large bulge (GEN.79). Old disk and bulge giants could never reach 1000 Lsun, so every bright star came from the young thin disk near the plane. A low-mass giant now spends a short bright tip (2% of its giant phase) between 1000 and 2500 Lsun, and the scatter reports how many layers actually drew stars.
- A neighborhood started from a sparse sector at the galaxy's edge fills every sector in range (GEN.77, fixed by GEN.76; now tested), and skip notes no longer mention the old one-star threshold.

## [7.193.608] - 2026-10-07

### Changed
- TODO: MAP.127 (neighbouring blocks and slabs on the Galaxy Map) is folded into MAP.121, which now covers seeing and stepping into neighbouring regions on every map, built as one map engine shared by the Galaxy, Sector and System Maps.

## [7.192.608] - 2026-10-07

### Changed
- TODO: Boss's three Galaxy Map problems are filed. MAP.109 (slow updates from the volume of star data) and MAP.116 (which stars and objects show at each zoom, now including MAP.115's rule for comets, rogue planets and asteroid fields) carry his words, and MAP.127 is new: see into and step to the neighbouring blocks and slabs without zooming out.

## [7.191.608] - 2026-10-07

### Added
- **A searchable TODO reference page**: `docs/plan/todo-reference.html`
  explains every TODO ID ever issued (full text, status, phase, thread,
  prerequisites and what each unblocks) with search and filters.
  `python scripts/build_todo_docs.py` rebuilds it and the tier plan from
  `docs/TODO.md` and the plan files, and a test fails when either page is
  out of date.

## [7.190.607] - 2026-10-07

### Changed
- **The JSON API moves into `planetgen.api` (OPS.24, step 11 of 14).** `src/html/api/` is now `src/planetgen/api/`, with the same module names (`app`, `routes`, `auth`, `config`, ...). The Apache example drops its deny rule for the old folder, since the package sits outside the DocumentRoot. Every caller moved too.

## [7.189.607] - 2026-10-07

### Fixed
- **The bright-star backfill after a `galaxy` run has its own progress
  bar with a working ETA** (PERF.28). A `--slot` run (the Sector Map's
  "Generate neighborhood" among them) used to backfill inside the
  sector's bar, which sat at 0 of 1 with no ETA for minutes; the backfill
  now runs at the end of the run, counting the sectors it visits.

## [7.188.607] - 2026-10-07

### Added
- **Two browsable development pages in the docs**: `docs/design/architecture.html`
  (modules, flows, data and deployment, as of 2026-10-07 before the
  package move) and `docs/plan/tier-plan.html` (every open TODO item by
  phase, lane and dependency), linked from the README's documentation
  table.

## [7.187.607] - 2026-10-07

### Changed
- **The web pages' shared modules move into `planetgen.web.lib` and `planetgen.web.maps` (OPS.24, step 10 of 14).** apiclient, fmt, pagination, pagecache, tilecache, classref, tabledisplay, mdconvert, privatedir and systempage are now in `planetgen.web.lib`. The map renderers (starmap, systemmap, navmap, galaxymap, galaxymap3d, phenomenonmap, phenomenonrender) are now in `planetgen.web.maps`. `src/html/lib/` is gone, along with the `sys.path` pushes that reached it, and the Apache example no longer needs a deny rule for it.

## [7.186.607] - 2026-10-07

### Changed
- **The work queue and job runner move into `planetgen.queue` and `planetgen.cli.job` (OPS.24, step 9 of 14).** `workQueue`, `progressFile`, `progressRate` and `systemLoad` are now `planetgen.queue.work`, `progress_file`, `progress_rate` and `load`. `src/jobRunner.py` is gone; the Generate page starts `python -m planetgen.cli.job` from the checkout's `src/`. Every caller moved too.

## [7.185.607] - 2026-10-07

### Changed
- **The population modules move into `planetgen.population` (OPS.24, step 6 of 14).** `population` and `facilities` are now `planetgen.population.model` and `planetgen.population.facilities`. Every caller moved with them.

## [7.184.607] - 2026-10-07

### Changed
- **The database modules and command-line scripts move into `planetgen.db` and `planetgen.cli` (OPS.24, step 8 of 14).** `_db`, `editStore`, `systemRender`, `queryDb`'s queries and `adminStats` are now `planetgen.db.store`, `edits`, `render`, `query` and `stats`, with `schema.sql` and `control_schema.sql` beside them. The `src/` scripts are gone: run `python3 -m planetgen.cli.query`, `migrate`, `reset`, `orbits`, `dedupe`, `lockouts` or `render_parity` from the checkout's `src/` folder. install, update, the Generate page's reset step and the maintenance examples do this already. Run `update.sh` (or `update.ps1`) after pulling.

## [7.183.607] - 2026-10-07

### Changed
- **The admin modules move into `planetgen.admin` (OPS.24, step 7 of 14).** `adminAuth`, `loginThrottle`, `totp`, `qrcodegen`, `activitylog` and `adminEdits` are now `planetgen.admin.auth`, `throttle`, `totp`, `qrcode`, `activity_log` and `edits`. The common-password list and its license moved with them, and every caller moved too.

## [7.182.607] - 2026-10-07

### Fixed
- **install.ps1 and update.ps1 no longer fail when no Redis answers**:
  the Redis check stays a warning instead of leaving its exit code as
  the script's.

## [7.181.607] - 2026-10-07

### Fixed
- **A black hole without a disk shows its Hawking temperature and
  luminosity** (GEN.82) instead of zero, from its mass.
- **The sector summary names star kinds plainly and logs one line per
  entry** (UX.34, OPS.9): white dwarfs, neutron stars, black holes,
  giants and supergiants instead of raw spectral codes.
- **Species are kept only for worlds with a technological civilization**
  (GEN.80); a population pass removes any stored species without one.
- **Changing a planet's or moon's class re-generates it as that class**
  (ADM.27): radius, density, composition, atmosphere, temperature,
  pressure and life, keeping its orbit, mass and name. A class its mass
  can't fit (a gas giant class for a small rocky world) is refused with
  a message and nothing changes.

## [7.180.607] - 2026-10-07

### Changed
- **The generation modules move into `planetgen.generation` (OPS.24, step 5 of 14).** Stars, planets, systems, binaries, belts, comets, life, evolution, the star population, bright stars, limits, stats, validation and the two plausibility engines now live in `planetgen.generation`. The six interstellar phenomena live in `planetgen.generation.phenomena` (`asteroid_field`, `compact_remnant`, `nebula`, `quasar`, `rogue`, `supernova_remnant`). Every caller moved with them.

## [7.179.607] - 2026-10-07

### Changed
- **The name modules move into `planetgen.names` (OPS.24, step 4 of 14).** `names`, `nameUniqueness`, `bodyNames` and `objectId` are now `planetgen.names.wordlists`, `.uniqueness`, `.bodies` and `.object_id`, with `offensive_words.txt` beside them. Every caller moved with them.

## [7.178.607] - 2026-10-07

### Changed
- **The galaxy modules move into `planetgen.galaxy` (OPS.24, step 3 of 14).** The sector grid, density model, skeleton, drill-down blocks, viewport, seed, version key, sectors, sector colors, molecular cloud field and navigation now live in `planetgen.galaxy` (`geometry`, `density`, `skeleton`, `drill`, `viewport`, `seed`, `version_key`, `sector`, `sector_look`, `nebula_field`, `navigation`, `nav_graph`). `stellarObjects/__init__.py` no longer imports its classes eagerly, which would have created import cycles with the moved modules. Every caller moved with them.

## [7.177.607] - 2026-10-07

### Changed
- **The physics modules move into `planetgen` (OPS.24, step 2 of 14).** These modules moved:
  - `program_constants` is now `planetgen.tuning`.
  - `physical_constants`, `keplerMotion`, `planetPhysics`, `rogueSurface`, `stellarEvolution` and `mathCheck` are now `planetgen.physics.constants`, `.kepler`, `.planets`, `.rogue_surface`, `.stellar_evolution` and `.mathcheck`.

  Every caller moved with them. The math check runs as `python -m planetgen.physics.mathcheck`.

## [7.176.607] - 2026-10-07

### Fixed
- **Every sector inside the galaxy now generates, however sparse
  (GEN.76).** The predicted density only sets how many systems a sector's
  draw expects; it no longer decides whether the draw happens. A sector
  below one expected star is saved and marked generated like any other,
  in single-address, ring, column, shell, neighborhood and random-start
  runs and on a map visit. Only addresses outside the galaxy's stored
  outline are refused, and that error now says so.
- **`generate.py` gives a usage error for an option whose value is a lone
  `--`** (`--workers=--`) on Python 3.9 too, instead of crashing
  (PERF.27).
- **CI's update.sh migration check rolls a database back to v48 again**
  (OPS.25): it dropped `bright_star_blocks`, which v53 removed.
- **The -0.0 round-trip test accepts either sign of zero** (DB.12):
  MariaDB 11 reads a stored -0.0 back as +0.0.
- **Changing a star keeps every planet's class** (TEST.80): re-spacing
  the orbits afterwards could reclassify a planet it nudged into another
  zone.
- **Two processes reserving ids for the same table no longer deadlock**
  (TEST.81, TEST.87): the reservation is one upsert that takes its row
  lock exclusively from the start. The two-process test reports a
  worker's error instead of timing out on an empty queue.
- **The Sector Map can always zoom out past its opening view** (MAP.114):
  a crowded sector opened at the zoom floor, so the - button and Reset
  view did nothing.

## [7.175.606] - 2026-10-07

### Fixed
- **`pip install` failed after the first package move.** `setup.py` mapped `planetgen.util` to a folder named `src/planetgen.util`; subpackages now map to their real folders.

## [7.174.606] - 2026-10-07

### Changed
- **The code starts moving into the `planetgen` package (OPS.24, step 1 of 14).** `log`, `appconfig` and `serialization` are now `planetgen.util.log`, `planetgen.util.appconfig` and `planetgen.util.serialization`. The version file is now `src/planetgen/_version.py`. Every caller moved with them, with no compatibility stubs left behind.

## [7.173.606] - 2026-10-07

### Added
- **The migration libraries are pinned, and install, update and CI provide Redis (OPS.21).** `setup.py` and the hash-pinned lock files now carry redis, rq, pyotp, segno, markdown, cachetools, SQLAlchemy, Alembic, Pydantic, numpy, scipy, astropy and scikit-image, so each library swap only changes code. `config.json` has a new `redis.url` (default `redis://127.0.0.1:6379/0`). `install.sh` and `update.sh` check that Redis answers there; on Linux with a local address they install and start `redis-server` when it doesn't. On macOS they say to run `brew install redis`. `install.ps1` and `update.ps1` check and point to Memurai or Redis in WSL. CI's test job runs a Redis service. Nothing uses Redis yet, so a missing server only warns. diskcache is left out: its last release has an unfixed advisory (PYSEC-2026-2447).

## [7.172.606] - 2026-10-07

### Changed
- **A package layout plan for the code reorganization (OPS.23).** `docs/design/library-migration.md` section 6 names the `planetgen` packages, which module goes where, the shared helpers merged into utility modules, and the order of the move (OPS.24). No compatibility layer: each move updates every caller.

## [7.171.606] - 2026-10-07

### Changed
- **The TODO list and phase plans were rebuilt from Boss's lists of
  2026-10-03 and 2026-10-07.** 130 new items were filed and duplicates
  merged. Every bug is now in phase 0, which has two lanes: bug fixes and
  groundwork. The groundwork covers the package layout, the move to
  third-party libraries with Redis, the data model, names from IDs, the
  job queue and logs, web components, the map engine and nebula shapes.
  Bugs that the groundwork fixes are folded into it. Phases 1 to 3+ were
  rebuilt to match. New design notes: `docs/design/library-migration.md`,
  `docs/design/habitability-index.md` and `docs/design/orbital-updates.md`.
  The other design documents gained the planned changes. Docs only; no
  code changed.

## [7.170.474] - 2026-10-03

### Fixed
- **Stars on the Galaxy Map no longer take the click (MAP.101).** They are still drawn, but a click or tap on one picks the block, slab or sector under it, so a sector thick with stars can be picked at any zoom. A star's details are on its sector's page. Black holes, neutron stars, quasars and small clouds still open their panel when clicked.

## [7.169.474] - 2026-10-02

### Changed
- **Galaxy Map sectors and blocks are colored by what is in them (MAP.86).** A generated sector is now translucent, only a little more solid than unfilled space, in its own color: the hue of its stars' average color, saturated by how full it is (from empty to the densest a sector can be) and lit by its stars' average luminosity against the Sun. A generated sector with no stars is the unfilled color, a shade more saturated and solid. A block or mega block averages the color and opacity of every sector in it, unfilled ones counted as unfilled. Filled blocks keep their amber edges. `GET /api/galaxy/stage` gives each child block, and each listed sector, a `look` (`share`, `color`, `colored`) from the stored sector stats.

## [7.168.474] - 2026-10-02

### Fixed
- **The Galaxy Map shows every star of a sector it is zoomed to, and its black holes, neutron stars and quasars (MAP.80).** Zoomed in to about a sector, the map also fetches the finest (16 pc) tiles around it, which list every generated star, so the sector is no longer thinned to its brightest stars. Tiles of 64 pc and finer now list the placed black holes, neutron stars (pulsars) and quasars (`points` in `GET /api/galaxy/tiles`). The map draws them among the stars in their own colors (violet, mint green and pink), and a click names one and links to its page.

## [7.167.474] - 2026-10-02

### Fixed
- **Courses between separately generated areas now find a route (NAV.34).** The NAV route graph links each system to its 6 nearest, so any separately generated area with 7 or more systems used to be an island with no way out, and a course from it to anywhere else said "No route via adjacent systems could be found". The islands are now joined, each to its nearest few islands by the closest pair of systems they have, until the whole graph is one piece, so any two placed endpoints always get a route. Same-sector routes, the hop-length limit and the unknown-space flag are unchanged (NAV.12).

### Added
- **Route edge-case tests (TEST.79).** Tests for the cases routing must handle, from the hop-length study: an isolated system, empty sectors between the endpoints, the galaxy edge and the halo, a graph in hundreds of pieces, and a same-sector course whose best route leaves the sector. The cases NAV.12 still has to build (the longest hop, the unknown-space flag, same-sector routes leaving the sector) are marked expected to fail.

## [7.166.474] - 2026-10-02

### Added
- **One stats row per sector (GEN.44, PERF.11, schema v53).** The new `sector_stats` table, keyed by grid address so an unfilled sector can have a row, holds each sector's bright-star level, its expected density from the galaxy model and, once it is generated, the systems and stars it got, the mean temperature and luminosity of those stars, and the color the Galaxy Map will draw it in (MAP.86: hue from its stars' temperature, saturation from how full it is, lightness from their luminosity). The galaxy keeps a decaying average of the systems fills got against the systems expected. The admin Stats page shows it as "Sector density", `GET /api/admin/stats` returns it as `sector_stats`, and `GET /api/sectors/<id>` returns a sector's own row as `stats`.

### Changed
- **The bright-star backfill works sector by sector (GEN.44).** Each sector keeps its own level: -1 untouched, the dimmest L_sun a backfill drew it down to, or 0 once generated. A backfill gives each sector the stars from its own distance tier's floor up to its level, and skips a sector already that deep. A band run (`plan --bright-stars-down-to`) tops up each backfilled sector from its own level, so a sector filled to 1000 L_sun that now needs 500 gets only the stars from 500 to 1000. A sector still at -1 that holds stars is left over from a failed run, so it is wiped and drawn again. A sector's stars come out the same whether it was taken down to 500 in one step or in two. Deleting a generated sector puts its old level back. This replaces the per-block levels of v49 (`bright_star_blocks`); the migration moves each block's level onto its sectors.

## [7.165.474] - 2026-10-02

### Fixed
- **Galaxy Map slab buttons: lines end at the slab's nearest edge, one-line labels, and they fit the window (MAP.98, MAP.100, MAP.99).** Each slab button's line now ends on its slab's outline as drawn on the map, at the point nearest the button, and keeps doing so as the view turns and zooms. Each button reads on one line, the slab number and how much of it is charted: "#4 Unknown" with nothing generated, "#2 < 0.01% charted", or "#6 ≈ 2.43% charted" (the full name and counts stay in the button's tooltip and screen-reader label). When one column of buttons is taller than the map, the buttons split into two columns, one each side of the map; if that still doesn't fit they shrink to the slab number alone, and when even that won't fit the buttons are left out and the box says to pick a slab on the map. On a phone, where the buttons sit below the map, they shrink and then give way the same way.

## [7.164.473] - 2026-10-02

### Fixed
- **Molecular clouds now appear across a generated region (GEN.47).** Dark nebulae (giant molecular clouds, dark clouds, Bok globules and star-forming cores) used to be rolled once per sector and placed inside it, about one per three thousand sectors, so they almost never appeared even where a generated region sat inside one. In a galaxy they now belong to the galaxy itself: they are drawn from the galaxy seed, more on the spiral arms and near the plane, and every sector a cloud reaches shows it, with the cloud stored once whichever of those sectors is generated first. On an arm at the Sun's distance about one sector in eight lies inside a dark cloud. Sectors already generated keep what they have; newly generated sectors next to them pick up the clouds that reach them.

## [7.163.473] - 2026-10-02

### Fixed
- **Rogue gas giants get a radius that follows their mass (GEN.60).** A free-floating gas giant used to be drawn around Jupiter's radius whatever its mass, so a 0.05 Jupiter-mass rogue was as big as a 10 Jupiter-mass one. Rogue gas giants now use the same giant mass-radius relation as gas giants around stars: a Neptune-mass rogue is about four Earth radii across, and from about Saturn's mass up the radius stays near Jupiter's. Rogue planets already in a database keep their stored size.

## [7.162.473] - 2026-10-02

### Added
- **Class S, the barren rocky super-Earth (GEN.38).** A new planet class for rocky worlds 1.2 to 1.8 times Earth's size and 2 to 10 Earth masses, with no atmosphere and no life, found close to a star, in the habitable zone and farther out. Class V stays the super-Earth that can carry life. Class S also appears on the class reference pages and the System Map.

### Fixed
- **Rocky rogue planets over 10,000 km are no longer classed C (GEN.38).** Class S can be a free-floating planet, so a rocky rogue too big for Class C's 10,000 km ceiling is now Class S. Rogue planets already in a database keep their stored class.

## [7.161.473] - 2026-10-02

### Fixed
- **The Galaxy Map turns freely, and each zoom step has its own camera (MAP.96, MAP.97).** Dragging, Shift and the arrow keys, or a one-finger touch drag now turn the map any way by any amount, past edge-on and round under the galactic plane, trackball style, without flipping at the poles; hovering and picking keep working at any angle. Each zoom step flies the camera to a preset for what it shows: the whole galaxy and a slab straight down (so the spiral arms show; this replaces the 35-degree opening tilt), and a block of several slabs (an arc, an entered block, the cube of sectors) at the isometric slant, going in and coming back out. A manual turn holds only within its step, and Reset view flies back to the step's preset.

## [7.160.473] - 2026-10-02

### Changed

- Every interstellar object in a galaxy-placed sector is now named by a
  76-bit position ID, shown as 19 hex digits, instead of a generated name
  or designation (GEN.64). The ID packs the object's type, distance from
  the galactic center (with its unit, mpc to Gpc), bearing, mark and a
  collision number that tells apart up to sixteen objects at the same spot.
  This covers rogue planets, standalone black holes and neutron stars,
  supernova remnants and their collapsed cores, nebulae, quasars,
  interstellar comets, asteroid fields, and star systems built around a
  bright-sweep star, whose stars, planets and moons are named from the
  ID. Ordinary star systems keep their names, and names given with
  `--name` are kept. These objects no longer go through the name
  registry, which made dense sectors about 1.8 times faster to generate
  with one worker and 2.3 times with two. Objects saved before this keep
  their names. See `docs/design/object-ids.md`.

## [7.159.473] - 2026-10-02

### Fixed
- **Every Galaxy Map stage is framed whole and fitted to the map's size (MAP.53, MAP.78).** The galaxy, an arc, a slab or a block now fills the map by its actual width and height, centred on its picture, so a bigger window shows it bigger and nothing is cut off at the edges; the view refits when the window is resized or turned, and stays whole while you turn it. Shift and the arrow keys now turn the view about the middle of what it shows, alongside dragging and a one-finger touch drag.

## [7.158.473] - 2026-10-02

### Fixed
- **While picking a slab the Galaxy Map draws only the lines between slabs (MAP.77).** Before, every block inside every slab was outlined while you were choosing a slab. Now the map draws just the boundaries between the slabs, and once you are on one slab it draws the divisions between that slab's blocks, which are the segments you pick next; this repeats at every level down to a sector. The whole galaxy still shows no lines.

## [7.157.473] - 2026-10-02

### Fixed
- **Slab buttons with leader lines instead of the slab slider (MAP.54, MAP.76).** While the Galaxy Map asks for a slab, each slab now has its own button beside the map (in one column below it on a phone), with a line from the button to its slab that is redrawn as the view turns, zooms, pans or the window resizes. The buttons are ordered by their slabs' height on screen so the lines don't cross; hovering or focusing a button lights its slab and its line, and clicking it picks the slab. A slab that falls outside the map gets a line ending in an arrow at the map's edge.

## [7.156.473] - 2026-10-02

### Fixed
- **The Galaxy Map drill-down is arc, slab, segment again (MAP.56).** Below an arc you now pick a slab, then click one block of that slab (a segment) to zoom into it, then a slab and a segment inside that block, and so on down to a sector. The 3 by 3 "region" pick between them is gone. A segment's URL token is `s<ring>.<wedge>`; an older link holding a region (`r4`) opens at the stage before it. Up and Down arrows move to the nearest block a ring further out or in.

## [7.155.469] - 2026-10-02

### Fixed
- **The System Map's side panel shows a planet's or moon's radius and mass (MAP.92).** Clicking or tapping a planet or moon (or the planet at the center of its moon view) now lists its radius in km and Earth radii and its mass in kg and Earth masses (Jupiter masses for a gas giant), alongside everything it showed before. An unknown value shows a dash.

## [7.154.469] - 2026-10-02

### Fixed
- **The whole star system fits on the System Map (MAP.88).** An outer planet, its ring or moons, a belt, a facility, a wide pair's companion star or a name no longer runs past the edge of the map. Once everything is placed, the map measures how far the drawn scene reaches and zooms out evenly around the star (or the planet, in a moon view) just enough to hold it all with a small margin. A scene that already fits looks as before.

## [7.153.469] - 2026-10-02

### Fixed
- **The System Map never draws a broken number (MAP.57).** A NaN or infinite value stored for a star, planet, moon, belt or facility no longer ends up in the map's SVG. A planet or moon whose position was lost is drawn at its orbit distance, due east of what it orbits, with a note in its info panel; one with no distance either is left out.

## [7.152.469] - 2026-10-02

### Fixed

- A sector save no longer fails with "existing_count must be >= 0, got
  -1" (TEST.85). When a system name with no decoration left was drawn
  again during name reservation, the new name was also counted as a
  holder of the registry row it matched later in the same pass, leaving
  that row one holder short. Each name's key is now taken once per pass,
  so a redrawn name only counts in the next pass. Names that never
  collide come out exactly as before.

## [7.151.469] - 2026-10-02

### Fixed
- **Hover lights the galaxy's arcs while picking a course (NAV.31).** While choosing a NAV start or destination on the Galaxy Map, hovering now lights and outlines the arc under the pointer, and at every later stage the choice under it, just as when browsing the map. Before, with "Generated only" forced on for the pick, only the few arcs holding generated sectors reacted at all. An arc or block with nothing generated still can't be taken, and its tooltip says so. The same applies whenever "Generated only" is on.

## [7.150.469] - 2026-10-02

### Fixed
- **The Galaxy Map's breadcrumb stays on one line (MAP.93, MAP.94).** A deep drill-down no longer wraps the breadcrumb onto several lines: when the steps don't fit, it shows the first step, a "…" button whose menu lists the hidden steps, and as many of the last steps as fit, then the current one, and it fits itself again when the window is resized. On a phone the breadcrumb line gives way to a round Steps button between Back and Forward, whose menu lists every step with the current one marked; Reset stays beside the arrows.

## [7.149.469] - 2026-10-02

### Fixed
- **Bookmarks work while picking a course (NAV.40).** While choosing a NAV start or destination, every bookmark keeps the pick and the end already chosen: in the Galaxy Map's Bookmarks menu (and its 1 to 9 keys) a system or phenomenon bookmark sets that end and opens the course, a sector bookmark opens that sector's page in pick mode, and a saved map view opens the Galaxy Map there, still picking. The sector page in pick mode now has a Bookmarks menu that works the same way, the NAV page's Bookmarks select also lists saved map views, and the course page offers bookmarks to change either end. The Galaxy Map also keeps the pick in its own URL as it moves, so Back or a reload no longer drops it.

## [7.148.469] - 2026-10-02

### Changed

- The `+name`/`-name` forcing options (`+habitable_world`, `-planets`,
  `+comets` and the rest) now work only for a single system: `generate.py
  system` and the one-off system page (GEN.51). `generate.py sector` and
  `galaxy` refuse them with an error that names the option and points to
  `system`, so a saved or queued command line that still has one gets a
  clear message. `--min-habitable` still works for sectors. Prevalence
  controls for sector and galaxy runs come later (GEN.52).

## [7.147.469] - 2026-10-02

### Fixed
- **Picking a slab lights the whole slab (MAP.91).** On the Galaxy Map, whenever the next pick is a slab (a layer of blocks, or a layer of sectors in the 3 by 3 by 3 view at the bottom of the drill-down), hovering any cube of it now lights and outlines that whole slab, round its full height, and fades the others; the slider beside the map marks the same slab, and a click picks it. Before, the 3 by 3 by 3 view lit and outlined a single sector, and the bigger slab picks gave the hovered slab no outline. A single cube is highlighted only when the pick really is one sector.

## [7.146.467] - 2026-10-02

### Changed

- Binary star names stay two words (GEN.62). A close pair's stars are now
  the system's A and B (`Voranthis A`, `Voranthis B`, shown as
  `Voranthis A / B` in the pair's table), and its planets are named for the
  system. A wide pair's stars share the system name's first word and add
  their own (`Voranthis Kelmoor`, `Voranthis Pikkita`); the second star's
  word is drawn from words for small, little, daughter, son and child, so
  it sounds like a diminutive. The primary's planets take the shared word
  (`Voranthis I`) and the secondary's take its own (`Pikkita I`). Planet
  names may repeat across systems. Only new systems get these names.

## [7.145.467] - 2026-10-02

### Added
- **The seed and version at the top of every run (OPS.10).** Every `generate.py` run now writes one line first, at normal level, to the console, the `--debug` file and the debug log: the galaxy seed, the PlanetGen release with its 22-hex-digit version key, and the command line (without the `--mysql-*` and `--debug` options), for example `Galaxy seed 3F2A...C901, PlanetGen 7.127.352 (0007007F000160030C0300), run: galaxy --ring 3`. A Generate page job's log starts with the same line for the job, before its first step.

## [7.144.463] - 2026-10-02

### Added
- **What made the galaxy, and every run (DB.6).** The galaxy now records the version key of the code that made it: 22 hex digits packing the PlanetGen release (MAJOR 4, REVISION 4, BUILD 6), the Python version (2+2+2), the OS and the architecture (1 each), for example `0007007F000160030C0300`, next to the full release, Python version and platform (schema v52). Every `generate.py` run that changes the galaxy adds a row to the new `generation_runs` table: its command line, its own seed, the galaxy seed, the version key, start and end times and outcome. A galaxy is built by a series of runs, not by its seed alone, so these rows are what a rebuild replays.

## [7.143.452] - 2026-10-02

### Added
- **One seed, one galaxy (GEN.39).** The galaxy now has a 128-bit seed, shown as 32 hex digits and stored with it (`galaxy_shape.galaxy_seed`, `BINARY(16)`, schema v51). `generate.py plan --seed <32 hex digits>` sets it; without it the first plan draws one and every later plan keeps it, and a different seed is refused once any sector exists. Every sector, the bright-star scatter, each band and each backfill block draws from its own seed, the SHA-256 of the galaxy seed and its address (for example `sector:12/3/0`), so the same seed fills a sector the same way at any `--workers` count and in any run that reaches it, on the same PlanetGen release. The galaxy has to be wiped and planned again: a galaxy made before has no seed.

### Changed
- Generation no longer draws from the operating system's random source: planet classes, life, moons, names and star positions all come from the sector's seed. A save retried after a deadlock replays the same draws, and the bright-star backfill picks its center sector by address, not save order.

## [7.142.452] - 2026-10-02

### Fixed

- `generate.py system` no longer saves a system that misses a forced body
  (GEN.49). A forced habitable world or asteroid belt that a star can't
  host (around an O5V star, a habitable world failed about two systems in
  five) is now tried on up to five whole systems; if none has it, the run
  exits with an error and saves nothing. Before, it saved the system anyway
  and reported success.
- `-planets` with `+asteroid_belt` is now refused, and for a single system
  `-planets` is also refused alongside the same keys, `num_orbits` or
  `slots` in a `--system-file` (GEN.50). Before, the belt was placed
  anyway.

## [7.141.445] - 2026-10-02

### Changed
- Galaxy Map: the whole galaxy is a 3D disk you can turn, tilt and zoom (to about twice as close), with no sector, block or wedge lines drawn on it, so its stars and spiral show through (MAP.85).
- Galaxy Map: the first pick is an arc instead of a quarter: about 40° of bearing (45°, snapped to wedge lines every block ring shares) by a third of the radius, the full height of the disk, 24 in all. Hovering one outlines it along its blocks' real sides, top and bottom, with its neighbors' outlines faint; the breadcrumb reads "Arc 90°–135° (middle)" and the URL `p=a1.90`. Older `q<n>` links still open (MAP.85, MAP.52).
- Galaxy Map: only Back, Forward, Up, Reset and Bookmarks stay in the controls row; Reset view, Generated only and Territories moved into a Menu (Escape or a click outside closes it). "Whole galaxy" is now Reset, the Wedges button is gone, and the info panel sits beside the slab slider when the row has room (MAP.55).
- Galaxy Map: the scale readout is one line, a bar and its length in sectors, pc and ly (MAP.60).

### Fixed
- Galaxy Map: hover outlines start and end on real wedge lines instead of crossing blocks (MAP.52), and an arc label no longer reads "360.0°" for a bearing a hair under 0.

## [7.140.445] - 2026-10-02

### Fixed

- Generation with more than one worker no longer stalls for 50 seconds at a time. A worker waiting for the neighbour lock kept the rows it had written, and when the lock's holder needed one of them (a name collision renames an existing system) neither moved until MySQL's lock wait timeout. The holder now gives up after 3 seconds and starts its sector save over. A seven-ring test run at 2 workers went from 266 s to 18 s (PERF.21).
- A generation run no longer hangs forever on Python 3.12 when a worker process dies. The pool stops the other workers with SIGTERM, but a worker caught it in the middle of a task, ended only that task and then waited for another one, so Python 3.12's pool cleanup waited on it for good. A worker sent SIGTERM now exits as soon as its task has been reported (PERF.22).
- Cancelling or interrupting a parallel run (SIGTERM, Ctrl+C) can no longer be lost. When the signal landed inside a database call it became a database error the work queue logged and carried on past, and the run finished as if nothing had happened. The queue now remembers the signal, stops at the next step and ends the job as cancelled (TEST.73).
- The bright-star progress bar no longer ends at 101%: a progress report that arrives after its layer has finished is ignored, and the terminal, `progress.json` and the Generate page never show a share above 100% (PERF.23).
- Re-running `generate.py plan --bright-stars-down-to` after a run of it stopped part way no longer puts the band in twice: the re-run first removes what the stopped run left below the star-fill level, then draws the whole band (GEN.32).
- A bright-star placement test now works on Python 3.9 and 3.10 (TEST.76).

### Added

- CI runs the generation tests again with 2 and 4 worker processes (Python 3.12 and 3.13), and tests can patch the workers as well as their own process (`tests/worker_patches.py`, described in `docs/testing.md`) (TEST.74, PERF.21).

## [7.139.445] - 2026-10-02

### Fixed

- New star system names are at most two words (GEN.46). A name collision
  is now resolved within two words: a one-word name gets a Greek letter
  ("Beta Vor") or, against a sector of the same name, one diminutive
  ("Little Vor"); where a decoration would add a third word (a two-word
  name, the Roman-numeral tier, a Greek letter plus a diminutive) a fresh
  name is drawn instead. Planet and moon names built on them stay short.
  Names already in the database are unchanged; sector names are not
  affected.

## [7.138.440] - 2026-10-02

### Fixed

- Both stars of a binary now share one age (GEN.53). A `--star-type`
  primary's companion is born at the primary's age, and the age adjustment
  for planets now settles one age for the pair: old enough for the planets
  around either star, and no older than either star's generated state
  allows. Before, a wide pair's stars were aged separately (an M2V of 13.8
  Gy next to an M8V of 0.45 Gy) and a `--star-type` companion rolled its
  own age.
- A `--star-type` primary's companion is now a real star for its mass
  (GEN.54). It comes from the population model at a fraction of the
  primary's mass, so its type, temperature and luminosity follow from that
  mass. Before, it kept the requested type's temperature and luminosity
  with an unrelated mass (a "G2V" of 0.17 Msun).

## [7.137.440] - 2026-10-02

### Fixed

- The Galaxy Map's tile-level helper (`galaxyViewport.tile_level_for_view_radius`
  and its copy in `lib/galaxymap3d.py`) no longer crashes with
  `OverflowError` on a subnormal view radius; any tiny positive radius gets
  the finest level (MAP.90).

## [7.136.440] - 2026-10-02

### Changed
- **The reproducible-galaxies design note names the "merge now" button
  (ADM.20, phase 3, low priority)** in place of its placeholder, with
  its lock and backup-slot rules. Documentation only.

## [7.135.440] - 2026-10-02

### Changed
- **Design and reference docs match the plan of 2026-10-02.** A new
  design note, `docs/design/reproducible-galaxies.md` (OPS.11), lays out
  the planned 128-bit galaxy seed, the 22-digit version key, the update
  history, the settings JSON file with its daily merge and 18 backups,
  the consistency check, parity repair and `generate.py reproduce`, each
  marked with its TODO item and phase. A new `docs/design/course-routing.md`
  describes today's routing and the planned routes with no hop limit and
  unknown-space jumps (NAV.12). The Galaxy Map drill-down design gains
  the planned arc pick (MAP.85) and an up-to-date build table; the
  database schema, API and command-line references list their planned
  changes, and the stale schema version numbers are corrected (galaxy
  v50, control v7). Documentation only.

## [7.134.440] - 2026-10-02

### Added

- `galaxyGeometry.sectors_along_segment` lists every sector a straight
  segment passes through, in the order the segment enters them, exact at
  sector faces and edges (NAV.38). Its browser twin is
  `galaxyprisms.sectorsAlongSegment`, and tests check the two agree.
  Course planning (NAV.12) uses it to flag hops through unknown space.

## [7.133.440] - 2026-10-02

### Fixed
- **Faint stars are brighter on the Sector Map and the Galaxy Map (MAP.87).** One shared curve draws the dimmest red dwarfs with four times the halo light, tapering to no change at 1000 L☉ and up, slowly enough that a brighter star is never drawn fainter than a dimmer one (a Sun gets about 2.2 times). Display only: nothing stored changes.
- **Rogue planets on the Sector Map are faint until marked (MAP.82, MAP.84).** Unmarked, a rogue planet is a dim speck with no glow, smaller than any star and only picked by a click right on it; "Mark rogue planets" makes each one bigger, fully lit, glowing and ringed, with a wide pick reach.
- **"Mark rogue planets" shows when it is on (MAP.83).** It starts off and stays highlighted while on, in both themes.
- **No link out of a course being picked (NAV.30).** While choosing a NAV start or destination, the Sector Map's info panel offers only the pick button, and the Galaxy Map's star and cloud panels drop "View system" and "View phenomenon".
- **Galaxy Map bookmark keys no longer clash with the browser's tab keys (MAP.81).** The first nine bookmarks open with plain 1 to 9 while the map has focus, instead of Ctrl+1 to Ctrl+9.

### Added
- **Browser tests for the maps that need no database (TEST.70).** `test_web_browser_fixture_maps.py` drives the Galaxy Map and the Sector Map, served by the real Flask views over fixture data, through picking, hover, keys, Back and Forward, URL state, bookmarks and the scale line.
- **One module for the maps' shared helpers (MAP.63).** `static/mapcore.js` holds what the Galaxy Map, the Sector Map and the System Map each had their own copy of: reading the scene data, theme colors, info-panel fields, the highlight ring, the scale bar's numbers, fitting the canvas and picking a point of light on screen. No visible change.
- **One camera and input controller for the Galaxy Map and the Sector Map (MAP.64).** `static/mapcontrol.js` turns, moves and zooms both maps, with a zoom policy each view sets (free, a short range, or locked), and tells a drag from a click in one place. No visible change.

## [7.132.433] - 2026-10-02

### Fixed

- A point one float step under layer 0's top face is now in layer 0, not
  layer 1 (GEN.31). The sector lookup (`galaxyGeometry.sector_address_at`
  and the Galaxy Map's `sectorAddressAt`) checks its ring, layer and slot
  against the cell's own bounds, so a point inside a cell's range always
  maps to that cell, in Python and JavaScript alike. The bright-star box
  query uses the same lookup.

## [7.131.433] - 2026-10-02

### Fixed

- Planets no longer pile up in the cold zone (GEN.37). The first planet now
  sits at 5-40% of the habitable zone's inner edge and each next one 1.3-1.8
  times farther out, stopping at the protoplanetary disk's outer edge, so a
  red dwarf's planets stay close in. About 30% of planets are now hot, 6% in
  the ecosphere and 64% cold (it was about 1%, 2% and 96%).
- Gas and ice giants have real densities, and super-Jupiters occur (GEN.34).
  A giant's mass is drawn first (dN/dlogM ~ M^-0.31 within its class) and its
  radius follows the Chen & Kipping (2017) mass-radius relation, so the median
  giant is about 1.2 g/cm^3 (it was 0.25) and about 8% of giants pass a
  Jupiter mass. A giant is at most 1% of its star's mass, so red dwarfs don't
  get super-Jupiters.
- Rocky planets get the moon classes their size and zone allow, not only
  Class D (GEN.35): a class only has to fit at its smallest size, and a moon's
  radius is capped so it stays under the planet's moon size and a tenth of its
  mass.
- A moon regenerated after its planet moves only gets a class and size a moon
  of that planet may have, never a gas giant or a blacklisted class (GEN.36),
  and is never too large for its planet (GEN.25). Moons are re-spaced after one
  is regenerated, and a moon heavier than a tenth of its planet is dropped
  when the planet is reclassified.
- A planet reclassified after a spacing push keeps the Hill sphere of its real
  orbit, not of the spot its new class first drew.
- Rogue planets follow one mass function, dN/dlogM ~ M^-0.65 from 0.1 Earth
  masses to 13 Jupiter masses (GEN.45): about 96% terrestrial and 4% gas
  giants (it was 91% and 9%), still 6.5 per star.

## [7.130.426] - 2026-10-02

### Fixed
- **Wiping the database (resetDb.py, or "New galaxy" on the Generate page)
  while the web app or a generation worker is still running no longer
  risks "Duplicate entry" errors (DB.3).** The id counters are now kept
  through a reset, so a new galaxy's ids carry on from where the old one
  stopped instead of starting again at 1, and a program still holding ids
  from before the reset can never be handed the same ones as a program
  started after it.

## [7.129.426] - 2026-10-02

### Fixed
- **Asteroid field and interstellar comet pages now show their saved
  composition (DB.2).** Each component was saved in its own row but never
  read back; the page showed the one-line summary stored beside it. The
  page now builds the Composition line from those rows, and
  `GET /api/phenomena/asteroid_field/<id>` and
  `/api/phenomena/interstellar_comet/<id>` return them as `composition`
  (a field's as `{component, concentration}`, a comet's as a list of
  components). No schema change.

## [7.128.426] - 2026-10-02

### Fixed
- **Several programs opening a new, empty database at the same moment no
  longer trip over each other (DB.5).** Each one used to create the tables
  itself and all but one could fail with "Duplicate entry ... for key
  'PRIMARY'". Now one creates them while the others wait, then carry on.
  The same goes for the control database (admin logins and the work
  queue), and for a migration running while another program connects.
- **A database whose schema version record was emptied or lost is no
  longer treated as up to date (DB.4).** Its version is now worked out
  from which tables and columns it has, so `migrateDb.py` (and update.sh)
  still runs the steps it is missing. No schema change.

## [7.127.352] - 2026-10-01

### Added
- **Tests for parallel generation, population passes, names and navigation
  (TEST.19, TEST.22, TEST.27, TEST.37, TEST.38, TEST.39).** One, two and
  three workers fill the same sectors with the same counts; every bulk mode
  (`--shell`, `--block`, `--column`, `--center-sector`, random start,
  `sector --num-sectors`) on two workers counts what it saved and never fills
  a sector twice; the ETA and progress file under a clock stepping back, NaN
  and infinite values and a tiny rate; sectors and systems with one name saved
  by several workers at once; population rescans after new fills; and the
  NAV k-d tree checked against brute force.

### Fixed
- **One worker didn't seed its tasks the way a pool does.** With `--workers 1`
  each task now draws from the same per-task seed a worker would use, and the
  run's own random stream carries on unchanged afterwards.
- **Two population passes at once could fail on a duplicate species.** A pass
  now holds a per-database lock, so a second pass waits for the first and then
  only adds what is new.
- **A NaN position in a NAV route gave NaN distances.** Systems without a
  finite position are left out of the route graph instead.
- **One bad progress value could break the ETA for the rest of a run.** NaN or
  infinite amounts and times are ignored, a clock stepping back no longer adds
  time to the estimate, units finished at a bar's very first instant are no
  longer dropped, and the web job's progress file never holds NaN or Infinity
  (which the Generate page's browser can't read).

## [7.126.352] - 2026-10-01

### Added

- Database tests (TEST.6-9, TEST.11-18): SQL portability and strict `sql_mode`, every released schema migrated and compared with a new database, crashed migrations re-run, every column written and read back, boundary values, collation collisions, CHECK constraints, failed sector saves, id-block edges, batch limits and search edge cases.

### Fixed

- Databases from v8-v37 can migrate again (steps re-run safely after a crash, and MariaDB column CHECKs are dropped correctly). Galaxy schema v50 brings migrated databases to exactly the new shape, keeps nebulae and asteroid fields when their sector is deleted, and stores the `--comets`/`--wide-binary` choices with a system's recipe.
- Saving one-off systems at the same time now names them Alpha/Beta like sector saves, and retries on a deadlock.
- Names that differ only by accent ("Vega"/"Véga") are treated as the same name everywhere, and the first holder keeps its own spelling.
- Very large batched saves are split to fit the server's packet limit.
- Search finds words longer than the full-text index's longest token.
- The Systems list works under MariaDB's strict GROUP BY mode.
- A failed sector save no longer leaves objects renamed.
- System and star names are limited to 200 characters, so their planets and moons always fit.

## [7.125.352] - 2026-10-01

### Added
- **Generation tests (TEST.4, 5, 23-26, 28-36).** New tests for resuming an
  interrupted ring, shell or block fill, bright-star scatter edge cases and
  interrupted scatters, `--force` scatter then fill, every `generate.py`
  argument error by its message with each limit at its maximum, the limits
  staying consistent, grid seams and the nucleus, sector placement running
  out of room, and direct tests of the system builder, moon, Kepler, star
  and phenomenon helpers. The old known-bug sweeps run over hundreds of
  seeds, the Tier 2 reports are now hard asserts or strict xfails, and the
  grid boundary tests also run at the real 4 pc sector edge.

### Fixed
- **The Kepler solver could return an unconverged answer.** When Newton's
  method runs out of iterations it now falls back to bisection, and it
  refuses an eccentricity outside [0, 1).
- **A float rounding slip could put a bright star in a zero-weight bin.**
- **A fill after a failed bright-star scatter could place the same bright
  stars twice.** The backfill now clears unbuilt leftovers in the cells it
  draws for when no scatter threshold was recorded.
- **A habitable world with no viable life chemistry could still found a
  civilization.** It got an evolutionary timeline anyway, so the population
  pass gave it a species that its system page never showed. Such a world now
  gets no timeline.
- **Three colony and population tests failed on some draws (TEST.69).**
  Two searched a page for a generated name with an apostrophe without
  escaping it; the third let a second capital keep random ages, which
  sometimes founded a polity with no systems.

## [7.124.339] - 2026-10-01

### Added

- **Install and update set up and check both logs on every OS (OPS.5).**
  On Linux and macOS, `install.sh` and `update.sh` (through
  `examples/apache/setup-debug-log.sh` and the new
  `examples/apache/log-locations.py`) now prepare the debug log whether or
  not `debug` is on, and the activity log's folder, wherever
  `PLANETGEN_LOG_FILE`/`log_file` and `PLANETGEN_LOG_DIR`/`log_dir` put
  them, and check that the web server's user and its group can really
  write each one. On Windows, `install.ps1` and `update.ps1` check that the
  app's account can write both logs' folders. Anything they can't set up
  (no rights, a read-only or missing drive, a folder in the way that is a
  file, a user that doesn't exist yet) only warns, with the exact
  `mkdir`/`chown`/`chmod` or `New-Item`/`icacls` commands that fix it, or
  how to point the setting somewhere writable; a log never stops an
  install or update.
- **CI runs `install.sh` and `update.sh` against a live database (TEST.62).**
  A new `linux-update` job runs both for real on Linux: a fresh install,
  an update with nothing new, a database that needs migrating, one newer
  than the code, an unreachable server and a failed migration. The shared
  database step is also tested on every database engine in pytest.
- **Command-line tests for the admin scripts (TEST.60):** `queryDb`,
  `adminStats`, `checkRenderParity`, `dedupeNames`, `resetDb`,
  `updateOrbits` and `loginLockouts`, including a bad port, an unknown
  database and an empty password against the environment for each.

### Changed

- **A database newer than the code is refused (TEST.10).** `migrateDb.py`
  (and so install and update), `migrateDb.py --status` and every
  read-write connection now stop with a clear message when the galaxy
  database or the control schema is at a higher version than this code
  knows, instead of carrying on (and before this code's older
  `schema.sql` touches it). `/api/health` says the database is newer than
  the code rather than telling you to run `migrateDb.py`.
- **`migrateDb.py --status` no longer changes anything:** it reads the
  version without laying down the schema, so it works with a read-only
  account, and reports a new empty database as current.

### Fixed

- `set-permissions.sh`, `create-cache-dir.sh` and `setup-debug-log.sh`
  stopped silently (exit 1, no message) on a Debian or Ubuntu server with
  Apache installed: reading `/etc/apache2/envvars` under `set -u` hit its
  unset `$APACHE_CONFDIR` and killed the user/group lookup. It is read
  safely now, and a lookup that still fails says so.

- `checkRenderParity.py`, `dedupeNames.py`, `resetDb.py` and
  `loginLockouts.py` print `error: ...` and exit 1 on a database they
  can't reach instead of a traceback.
- `updateOrbits.py` no longer fails every run with "elapsed_years must be
  >= 0" after the database server's clock went back; it moves nothing and
  restarts the clock from now.
- `loginLockouts.py --ip ""` is a usage error instead of quietly listing
  the lockouts.

### Removed

- `src/migrateSqliteToMysql.py` (TEST.61). It only accepted a SQLite file
  already at the current schema version, which no SQLite database ever
  reached (SQLite stopped at v5 and the MySQL migrations start at v8), so
  it could never import anything.

## [7.123.339] - 2026-10-01

### Changed
- **The tests now check that HEAD on the Generate page answers like a
  normal page load.** The bug was fixed in ADM.4; the test that pinned it
  as a known failure now checks the fix instead.

## [7.122.329] - 2026-10-01

### Fixed

- Galaxy Map: moving the view (right-drag or Shift-drag) now stops one and a half views from where the stage opened, as intended; before, nothing held it and the view could slide away for good.
- Galaxy Map: reloading the page after the map's Back button keeps Forward working.

### Added

- Tests for the page scripts under node (`src/tests/js`, run by `test_js_unit.py`): the Galaxy Map's drill-down, history, address bar, zoom, pan and tilt limits and buttons; the Sector Map's zoom, turning and picking; the phenomenon diagram's zoom; the Generate page's job panel; the facility form (TEST.57, TEST.58).
- Browser tests: every map button changes the view, the System Map's selection and measuring, the Galaxy Map drill-down by clicks with Back and Forward and the free camera (TEST.55, TEST.58, TEST.59), and no overlapping or off-screen controls on any page at 390, 600, 820 and 1280 px in both themes (TEST.56). The large-nebula diagram's dead "-" button is pinned as a known failure for UX.21.

## [7.121.284] - 2026-10-01

### Fixed

- A parallel run no longer stops when the control database drops mid-run: the task rows are only a record, so the run finishes its sectors and its lease goes stale on its own (TEST.20).
- When a worker process dies, the tasks it never got are recorded as cancelled and the ones in flight as failed, instead of being left as running (TEST.20).
- Cancelling a parallel run from the Generate page (or any SIGTERM to its process group) now ends it as cancelled: each worker rolls back its unfinished sector and the pool stays whole, instead of the workers being killed and the run recorded as failed (TEST.21).

## [7.120.278] - 2026-10-01

### Changed

- "Generate the neighborhood" on the Sector Map now starts a background job (followed on the Generate page and Admin, Queue) instead of running inside the page request, so closing the browser no longer matters (ADM.11).

### Fixed

- Two admins starting a job at the same moment could both start one: the job lock was briefly empty and the second caller cleared it as stale. The lock now appears with the job id in it, and only one caller clears a stale lock (TEST.40).
- A job id drawn twice in the same second gets a fresh one, and pruning old jobs never removes the running one (TEST.41).

## [7.119.271] - 2026-10-01

### Changed

- Bright stars now come after the sectors (GEN.30). A new galaxy generates its first sectors and then scatters the bright stars galaxy-wide, leaving those sectors out (`generate.py galaxy --then-scatter`).
- The bright-star backfill runs once, after a run has generated every sector it was asked for. By default it backfills only around the requested sector: the random start, the center sector or the slot address, or for ring, column, shell and block runs the generated sector nearest the middle. `--backfill-from all` (on the Generate page, "Backfill from every generated sector (farthest out)") backfills around every generated sector instead, and `--backfill-from none` skips the backfill.
- The scatter always leaves filled sectors out. The "Leave filled sectors out" checkbox is gone, and `--force` is accepted but no longer needed. A generated sector never gets scattered or backfilled stars.

## [7.118.271] - 2026-10-01

### Added
- **The Generate page's sections fold (ADM.4).** Click a section's
  heading (or focus it and press Enter or Space) to open or close it.
  Current job starts open; every other section opens or closes the way
  you last left it in this browser. A form shown again with an error or
  its size-and-time estimate stays open.
- **"Around a sector" can find its center.** Pick a filled sector by
  searching for its name (or a star system's name, which finds the sector
  it is in), or from a paged list of every filled sector, nearest the
  core first; or give a sector address (ring, layer and slot), or a
  galaxy-frame position in parsecs. An address or position that isn't
  generated yet is generated first, then its neighborhood.

### Fixed
- **A HEAD request to the Generate page no longer runs its form.** It
  is answered like a GET, so it can't skip the form's CSRF check.

## [7.117.271] - 2026-10-01

### Changed

- A new galaxy's galaxy-wide bright-star scatter now places every star of 1,000 solar luminosities and up (was 500). A galaxy already scattered keeps its level and its stars.
- The bright-star backfill around each generated sector is tiered by distance (GEN.30): down to 100 solar luminosities within 10 ly, 250 within 25 ly, 500 within 50 ly and 750 out to 100 ly. A block takes the tier of its nearest sector, and a block a nearer sector reaches later is topped up with only the band it lacks.
- The Generate page has a "Bright stars from (solar luminosities)" field on New galaxy, Plan and Rebuild the bright stars, passed to the scatter as `--bright-star-min-luminosity`.

## [7.116.271] - 2026-10-01

### Added

- An admin Queue page (`/admin/queue`, in the settings menu) to view and
  manage generation jobs (ADM.10): workers active and the server's load
  as "x / x / x" (CPU percent over 1, 5 and 15 minutes on Windows), the
  job trees with timing, progress and ETA on every node, and Pause,
  Resume, Cancel, Retry and Delete on any job, part of a job or failed
  sector. A paused job finishes its running tasks and stands by without
  holding up other jobs; "Pause the queue" stops every job from starting
  or taking another task until it is resumed. Each control is confirmed
  first and written to the activity log.

## [7.115.271] - 2026-10-01

### Added
- Tests for pages staying fresh after a command-line write (TEST.42) and
  the caches under real threads (TEST.54).
- Tests for who may call what across the whole API, generated from the
  route list (TEST.43), what an API key may do (TEST.44), more than one
  admin (TEST.45), trusted-device and two-factor edge cases (TEST.46),
  oversized requests (TEST.47) and security headers on every kind of
  response (TEST.48), the thinly tested API routes (TEST.49) and Galaxy
  Map URLs combined (TEST.50).

### Changed
- A request body over 2 MB is refused with a 413 (JSON under `/api`, the
  error page elsewhere) before the login checks or the database see it.
  Before, any size was read into memory.
- An API key can no longer make API keys, change credentials, set up or
  turn off two-factor sign-in, or log out: those answer 403 and need a
  browser sign-in, so a leaked key can't make itself a replacement.
- Turning two-factor sign-in off forgets every trusted device of that
  admin; the browser that turned it off gets a new one.

### Fixed
- An admin route whose control database can't be reached answers 503
  "database unavailable" instead of a 500.
- A phenomenon added from the command line (`generate.py phenomenon
  --sector-id N`) now shows up on cached sector pages and Galaxy Map
  tiles; before, it marked nothing as changed, so cached tiles never
  showed it.
- `/galaxy/locate` answers JSON when the API says "not found" (for
  example an unknown database), like `/galaxy/tiles` and `/galaxy/stage`
  already did, instead of the HTML 404 page.

## [7.114.271] - 2026-10-01

### Added
- **Bulk generation checks the math first (TEST.68).** `generate.py
  check-math` runs the math check by hand (`-v` lists every check). Every
  bulk run (`galaxy`, `plan`, `population`, and `sector` with more than one
  sector), every Generate page job (its new first step, "Check the math")
  and the Sector page's neighbourhood generation run it first and refuse
  to start if a check fails, naming the failed checks and writing nothing.
  `update.sh` and `update.ps1` run it after updating (a new step 4) and
  warn if it fails, skipping the population pass.

## [7.113.271] - 2026-10-01

### Added

- Every generation job is now a tree, with timing on every node (ADM.12):
  a Generate page job, its steps, each `generate.py` run, its phases
  (skeleton, bright stars, population), its work queues and their sectors
  or layers, each with its own state, start, end and duration, and
  totals, timings and an ETA added up from the children. One-worker runs
  are recorded too. Control schema v7: run `update.sh` (or `update.ps1`)
  after updating.

## [7.112.265] - 2026-10-01

### Added
- **A math check that runs first (TEST.63 to TEST.67).** A new module,
  `planetgen/physics/mathcheck.py`, checks the generator's math against 49
  fixed answers before anything trusts it: known values from real
  astronomy (the Sun's luminosity, lifetime and temperature, Earth's and
  Jupiter's orbits, Earth's Hill sphere, the habitable zone and snow line,
  white dwarf sizes, the Sun's Schwarzschild radius and galactic orbit,
  Holman & Wiegert's stability limits, the Kepler and Barker equations),
  identities that hold for any input (unit conversions, constants that
  agree with each other, energy conservation around an orbit, the sector
  grid), and seeded draws from every sampler (star masses, ages, sector
  counts, planet sizes and classes) checked against their intended shares.
  Each check says where its expected value comes from. It takes under a
  second. The test suite runs it before any test and stops if it fails,
  CI runs it as its own first job, and the website runs it at startup and
  shows admins a warning if it fails. `python -m planetgen.physics.mathcheck -v`
  prints the report.

### Fixed
- **The speed of light disagreed with the light-year.** `SPEED_OF_LIGHT_M_S`
  was 2.998e8 m/s while the light-year used the exact 299,792,458 m/s; it
  is now exact too, so black hole and quasar event horizons come out
  0.005% larger.

## [7.111.254] - 2026-10-01

### Changed
- Tests: a test that runs past 10 minutes now fails with a stack trace of where it was stuck (pytest-timeout, in the `test` extra), and CI's test jobs stop after an hour and list the 25 slowest tests.

## [7.110.254] - 2026-10-01

### Added
- **Bright stars around every generated sector (GEN.23).** As soon as a
  galaxy sector is generated, anywhere and by any route (the command
  line, the Generate page, the Sector page, or a Galaxy Map visit), every
  sector block (3x3x3 sectors) within 100 ly of it gets every star from
  100 L_sun up to what was already placed there. Each block remembers how
  dim it has gone (new `bright_star_blocks` table, schema v49), so a
  block is only drawn once, and sectors that are already filled are
  never touched. The new stars show on the Galaxy Map like the plan's
  bright stars, and later sectors in those blocks build their systems
  around them.

### Changed
- **The default generate-around sphere is 12 pc (about 39 ly), not
  100 ly.** "About 10 pc, rounded up" to the next whole 4 pc sector. It
  applies to `generate.py galaxy`'s random start, the Sector page's
  neighborhood button, `POST /api/sectors/<id>/generate-neighborhood`
  with no radius, and the Galaxy Map's neighborhood dialog. A radius you
  give is still used as is.

Run `update.sh` (or `update.ps1`) after updating: it migrates the
database to schema v49, adding one empty table. No regeneration is
needed; the backfill starts with the next sector generated.

## [7.109.254] - 2026-10-01

### Added

- Rogue planets have surface conditions (GEN.26). With no star, a rogue's only heat
  is its own: radioactive decay and leftover formation heat in a rocky
  rogue, slow cooling in a giant or brown dwarf. From that the generator
  works out its age, heat flow, effective temperature and surface: bare
  frozen rock, a frozen-out atmosphere, an ice shell over a liquid ocean,
  ice to the rock, a thick hydrogen envelope warming the ground (with an
  ocean beneath when it is warm enough), or, for a giant, the temperature
  at 1 bar. The rogue's page and description show it all; "Internal Heat"
  is now "Geologically Active", computed rather than rolled. The model and
  its defaults are in `docs/design/rogue-planet-surface.md`. Schema v48
  adds the columns and fills them for stored rogue planets, so run
  `update.sh` (or `migrateDb.py`).

## [7.108.254] - 2026-10-01

### Added
- **Change a planet's or moon's class, or a system's star (ADM.6,
  ADM.7).** On the system page's Edit panel an admin can pick a new class
  for any planet or moon: the menu lists the classes that fit where it
  is without moving any other planet first, and every other class under
  "Force", which is kept even where it couldn't form. A single-star
  system can have its star replaced by any spectral type: planets,
  moons and belts keep their classes and their orbits scale with the new
  star's light, and the page lists anything that had to move or that no
  longer had a stable orbit and was removed. New API routes: `GET
  /api/systems/<id>/class-options`, `POST /api/planets/<id>/class`,
  `POST /api/moons/<id>/class` and `POST /api/systems/<id>/star`.

## [7.107.254] - 2026-10-01

### Changed

- The bright-star scatter's progress bar and its time left now count each layer's expected stars (a quick sample of the density model before the scatter starts, plus a little for every ring walked) rather than layers, and move as layers in progress report their stars, so the thin layers at the edges of the disk no longer throw the estimate off (PERF.9). It shows a share done, "Bright stars (12 of 1,271 layers) 34%".
- While layers take longer than 30 seconds each, a second bar under it shows the stars of the layers being drawn, done of their estimate, with its own time left; it goes once layers finish faster than one every 20 seconds (PERF.4). The Generate page shows both, and each layer's time goes into the speed stats as a "scatter" task.

## [7.106.254] - 2026-10-01

### Added
- **Delete and Regenerate buttons on sectors, systems, planets, moons,
  asteroid belts and phenomena (ADM.8).** An admin sees an Edit panel on
  the sector, system and phenomenon pages; each button asks for
  confirmation first. Regenerating rolls the object again in the same
  place and keeps its name; after a body changes, the rest of the system
  is re-checked and re-spaced until it is stable, and the page says what
  moved and anything still unstable. Deleting a sector removes everything
  in it and leaves its place in the galaxy free to generate again. New
  API routes: `DELETE`/`POST .../regenerate` on `/api/planets`,
  `/api/moons`, `/api/belts` and `/api/phenomena/<type>`, plus
  `DELETE /api/sectors/<id>/contents` and `POST /api/sectors/<id>/regenerate`.

## [7.105.254] - 2026-10-01

### Added

- Every bulk generation now shows its size and time before it writes
  anything, and refuses one the database disk can't hold (PERF.3).
  `generate.py galaxy` (every mode) and `generate.py sector
  --num-sectors` print the expected sectors, star systems, database
  growth (+10%) and time, and on a terminal ask `Generate these N
  sectors? [y/N]` (`--yes` skips it). A run that would take more than a
  quarter of the database disk, or leave under 5 GB free, is refused
  with what it needs. `--estimate-only` prints the estimate and stops.
  The Generate page, the Galaxy Map's and Sector Map's Generate buttons,
  and a sector page's "Generate more sectors around this one" show the
  estimate and ask before starting.
- The server records how fast it generates, per star density (PERF.10):
  every filled sector adds its time to a log-scale density bucket (two
  per decade from 0.01) as a decaying average, and each galaxy's bytes
  per star system are measured after every run. The estimates use them;
  the admin Stats page and `GET /api/admin/generation-stats` show them.
  Control schema v6 (`generation_stats`, `generation_size`): run
  `update.sh`.

## [7.104.254] - 2026-10-01

### Added

- Bright stars can now be scattered in stages (PERF.5). `generate.py plan
  --bright-stars-down-to 100` keeps every bright star already placed and
  adds only those from 100 up to the galaxy's current star-fill level (500
  by default), then lowers the level. The Generate page shows the level
  and has a new "Add a dimmer layer of bright stars" panel. Sectors already
  filled are left out, since their own systems already include stars that
  bright. No database change.

## [7.103.254] - 2026-10-01

### Added

- Rogue planets have a planet class (GEN.8). Each planet class now says
  whether a rogue planet can have it (a new `"r"` zone flag in
  `PLANET_CLASSES`): dead worlds (C), icy bodies (D), ice giants (I), gas
  giants (J) and gas dwarfs (T). A rogue draws its class from those that fit
  its type, radius and mass; a brown dwarf has none. The class shows on the
  rogue's page (linked to the class page), in the sector's Contents and on
  the Sector Map, and the class reference lists "Interstellar space" as a
  zone. Schema v47 adds `rogue_planets.planet_class` and gives stored rogue
  planets their class, so run `update.sh` (or `migrateDb.py`).

## [7.102.253] - 2026-10-01

### Added
- **Test suite markers and more database engines in CI.** Every test is
  marked `db`, `slow` or `browser` as it applies, so `pytest -n auto -m
  "not db and not slow"` is a one-minute loop. CI's test job now runs on
  MySQL 8.0, MySQL 8.4 and MariaDB 11.4 (local runs cover MariaDB 10.11).
  `test_todo_tags.py` takes its category list from `bump_version.py`, so
  `TODO(TEST.N)` tags are accepted.
- **A guide to running CI on your own computers** (`docs/ci-runners.md`):
  what each runner needs, how to register it, security settings,
  troubleshooting. Pull requests from forks now always run on
  GitHub-hosted runners, never on self-hosted ones, and the browser job
  no longer needs passwordless sudo.

## [7.101.253] - 2026-10-01

### Removed
- **Ten tests (2,819 parametrized cases) that no longer checked anything.**
  Assert-free "generates
  without error" tests whose very next test runs the same generation and
  asserts on it (planets, the full star-by-class matrix, every star type,
  the example files), three comet-orbit validation tests repeated word for
  word in `test_kepler_motion.py`, three diminutive-prefix tests covered by
  the fuzz walk over every prefix, a check that a constant equals its own
  definition, a `None` Markdown check already in `test_mdconvert.py`, and a
  one-off TODO-marker check `test_todo_tags.py` now covers. The guard that
  `spaceSector.py` never imports a root script checked for `systemGen`,
  which no longer exists; it now checks `generate` too.

## [7.100.253] - 2026-10-01

### Fixed

- The Species and Polities pages answered 404 on MySQL 8 even after a
  population pass: the population status probe aliased a column as
  `generated`, a reserved word in MySQL 8 (not in MariaDB), so the query
  failed and every population page stayed hidden. The alias is now quoted.

## [7.99.253] - 2026-10-01

### Changed
- **The slowest tests run in seconds, and two tests that skipped at
  random now always run.** The radius-neighborhood generation test builds
  a sparser neighborhood (70 s to 3 s) and also checks every generated
  sector is inside the radius; the sector-enumeration fuzz test checks
  tolerances with `math.isclose` (30 s to 2 s). Tests hash passwords with
  1,000 PBKDF2 rounds instead of 600,000, except the tests of the hashing
  setting itself (marked `real_password_hashing`). The moon life-data
  round trip used to skip almost every run (a moon with its own life data
  is too rare to wait for) and the black-hole test skipped on one draw in
  ten; both are now seeded and always run. The full suite on 4 workers
  went from 5 min 13 s to 4 min 27 s.

## [7.98.253] - 2026-10-01

### Changed
- **Surface temperatures and pressures show customary units alongside
  metric.** A planet's or moon's surface temperature now reads in K with
  °C and °F ("288 K (15 °C, 59 °F)") in its description and on the
  System Map's info panel, and atmospheric, internal and core pressures
  read on a Pa, kPa, MPa, GPa ladder with atm and psi alongside
  ("101 kPa (1 atm, 14.7 psi)"), through the new
  `stellarObjects.utils.format_temperature_k` and `format_pressure_pa`.
  Star temperatures stay in K alone.

## [7.97.253] - 2026-10-01

### Changed
- **Stars and glowing phenomena are points of light on the Sector Map
  (MAP.15).** Every star is now drawn the way the Galaxy Map draws its
  bright stars: a tiny bright core in a soft aura a fixed number of
  pixels across, the core sized by the star's radius, the aura's width
  and brightness by its luminosity, and the color by its temperature,
  instead of a textured ball with a glow shell. Quasars, neutron stars
  and accreting black holes are points of light too; quiescent black
  holes, rogue planets (still ringed by "Mark rogue planets") and
  interstellar comets keep their spheres, and nebulae, supernova
  remnants and asteroid fields stay clouds. Points grow a little (up to
  1.5 times) as the camera closes in, a click within a few pixels of one
  still selects it and shows its details, and on the light theme a
  point's core is a darker shade of its color so it stays visible.
  `lib/starmap.py` now sends each star's point (`light`: color, core
  and aura size in pixels, aura strength, core brightness) in place of
  the old glow-shell numbers.

## [7.96.253] - 2026-10-01

### Changed
- **Speeds and time periods are always shown in a meaningful unit**
  (UX.13, UX.14). Speeds go through one shared ladder, km/h, km/s, Mm/s,
  then multiples of c from a tenth of light speed ("29.8 km/s",
  "4.5 Mm/s", "0.25 c"); periods and durations through another, µs, ms,
  s, minutes, hours, days, years, ky, My and Gy ("27.3 days",
  "1.88 years", "236 My"), each to three significant figures. Orbital
  periods and speeds of planets, moons, comets, binaries, wide binaries
  and facilities, galactic orbits, pulsar spin periods, NAV travel
  times, the facility form's live readout and the admin pages' uptimes
  all use them, in Python (`stellarObjects.utils.format_speed_kms`,
  `format_duration_seconds`, `format_period_years`) and in the browser
  (`static/speed.js`, `static/period.js`). The old "x years y days z
  hours" period text (`years_to_time_string`) is gone.

## [7.95.253] - 2026-10-01

### Changed

- Galaxy Map: slabs are picked with a vertical slider to the right of the
  map instead of a list under it (MAP.30). It has one step per slab, so
  any slab is one pick; dragging fades the other slabs and shows the
  slab's generated share, and letting go (or Enter, or Open) takes it.
  The map is now 4:3 and no taller than the window (1:1 on a phone).
- Galaxy Map: inside a block, whenever a slab is to be picked, the view
  opens from an isometric slant so the layers can be told apart.

## [7.94.253] - 2026-10-01

### Added

- Bookmarks (MAP.23, which also finishes the NAV page's map picks,
  MAP.22, and so the Galaxy Map drill-down, MAP.2). A ☆ on the Galaxy
  Map's breadcrumb saves the view (its stage URL) or the selected sector,
  and a ☆ Bookmark button on system, phenomenon and sector pages saves
  that page; it turns into ★ Bookmarked, and pressing it again offers to
  remove the bookmark. A Bookmarks menu in the map's controls opens,
  renames and deletes them, and Ctrl+1 to Ctrl+9 open the first nine
  while the map page has the focus (not while typing in a box). The NAV
  page offers the system and phenomenon bookmarks as a start or a
  destination, and a sector bookmark opens that sector's system picker.
  Bookmarks are kept in this browser (`localStorage`, one list per
  database, up to 100), in the new `static/bookmarks.js`; with storage
  blocked there are just none.

## [7.93.253] - 2026-10-01

### Changed
- **The test suite runs in parallel.** `pytest -n auto` (pytest-xdist, now
  in the `test` extra) runs one worker per core: the full suite went from
  17 min 21 s to 5 min 57 s on 4 cores, and CI and the weekly deep fuzz
  run use it. Each test process now gets its own throwaway control
  database instead of falling back to `planetgen_control`, so workers
  never share one and the tests never touch a real control schema.
  Plain `pytest` still runs serially.

## [7.92.253] - 2026-10-01

### Fixed

- Galaxy Map: picking a quarter or an arc now zooms into the wedge it highlights (the bearings its blocks really cover, set by the meridians that run to the core) instead of an even 90° or one-third share of the view.

### Changed

- Galaxy Map: zoomed into a wedge, only that wedge shows (wedge lines, stars and clouds outside it are hidden), and every view below the whole galaxy opens at an isometric slant so layers can be clicked on the map as well as picked from the list.

## [7.91.253] - 2026-10-01

### Added
- **One module to validate a planet, a lunar system and a star system
  (ADM.5).** `planetgen/generation/validation.py` holds checks that report what
  is wrong without changing anything, the orbit-spacing pass generation
  already ran (moved there unchanged), and a stabilize pass that re-spaces
  an edited system from its moons outward, ready for the admin overrides
  (ADM.6, ADM.7).

## [7.90.177] - 2026-10-01

### Changed
- **The roadmap lists only open work again.** `docs/TODO.md` drops the
  finished items still in it (the login security items SEC.1 and SEC.20
  to SEC.28, the database-call items PERF.6 to PERF.8 and PERF.12
  to PERF.17, and the finished Galaxy Map and Sector Map parents MAP.5,
  MAP.11, MAP.14 with their fixed bugs), and its plan now shows the
  current round of work (PERF.5 first; then units and the Galaxy Map
  layout, generation estimates and progress, and admin editing).
  `docs/design/todo-number-map.md` marks the items fixed since its last
  update as done and adds the IDs it was missing (MAP.19, MAP.51, GEN.9,
  ADM.4). No IDs changed.

## [7.89.177] - 2026-10-01

### Added
- **Bright stars are drawn several layers at a time (PERF.7).**
  `generate.py plan` hands each layer of the galaxy to a worker process,
  like sector generation; three workers placed the same stars about
  three times faster (313 s down to 105 s at 20,000 L_sun on a 4-core
  machine). `--workers` sets the count. Each layer has its own random
  stream, so the same seed gives the same stars on any number of
  workers.
- **The Generate page shows the time left** ("about 4 m 10 s left") for
  a running job, the same estimate the terminal shows.

### Changed
- **Steadier time-remaining estimates (PERF.7).** Every generation
  progress bar's remaining time now comes from a decaying average of
  sectors (or layers) finished per second, weighted toward the last
  minute, instead of rich's own short-window estimate, so it holds
  steady while several workers report at once and stays on screen until
  the bar is done. A run stopped early by `--limit` now ends its bar
  finished rather than part way.

## [7.88.177] - 2026-10-01

### Added
- **Sectors are generated several at a time (PERF.8).** `generate.py
  sector` and every `generate.py galaxy` mode, from the command line or
  the Generate page, fill sectors in parallel worker processes: by
  default 80% of the machine's cores (one fewer when MySQL runs on the
  same machine), at a lower priority than everything else on it.
  `--workers N` (or `PLANETGEN_WORKERS`) sets the count, and
  `--workers 1` works one sector at a time as before. On a 4-core
  machine with MySQL local, 60 sectors took 9.8 s instead of 16.5 s.
  Only one run's workers use the machine at a time: a second run waits
  for the first, using a lease in the control database (schema v5,
  `work_jobs`/`work_tasks`/`work_lease`).

### Changed
- **Linking a sector to its neighbors is about four times faster.**
  The nearest-systems search no longer walks a mostly empty grid for
  every object next to a newly saved sector; five neighboring sectors
  went from 13.9 s to 3.4 s.

Run `update.sh` (or `update.ps1`) after updating: it adds the work
queue's tables to the control database.

## [7.87.177] - 2026-10-01

### Changed

- Galaxy Map: the layers of a small cube of sectors now touch instead of being pulled apart, so no space shows between blocks or layers anywhere in the drill-down.

## [7.86.177] - 2026-10-01

### Changed
- **Search matches whole words (PERF.16).** Name search uses a FULLTEXT
  index (schema v46) instead of scanning every row for the text, so
  "ara" no longer finds "Kemaral", while "Mu" or "Tobar IV" still find
  "Ossiran Mu" and "Tobar IV". The Galaxy Map's locate box still finds a
  name as you type the start of its last word. Result counts stop at
  300 and show "300+".
- **Web pages run fewer queries (PERF.15).** System pages load moons in
  one query, sector pages load their stars in one query, system lists
  read star types once per page, and search's filter lists are cached
  until a sector or system changes.

### Added
- **A time limit on web database queries (PERF.17).**
  `mysql.statement_timeout_seconds` in `config.json` (default 10, 0
  turns it off) stops any one query on the web interface and API; the
  page says "Took too long" (HTTP 504) instead of hanging.

Run `update.sh` (or `update.ps1`) after updating: it migrates the
database to schema v46, adding full-text indexes to five tables. Writes
to those tables pause while each index builds, which can take a few
minutes on a large galaxy.

## [7.85.177] - 2026-10-01

### Changed
- **Sectors save about 2.5 times faster, with 30 times fewer database
  statements (PERF.12, PERF.13).** A filled sector's rows are now written
  as one multi-row INSERT per table instead of one INSERT per star,
  planet, moon, belt, comet and child row, and generation checks the
  schema once per process instead of on every connection. Measured on a
  40-system sector (median of five saves): 1.35 s and about 5,650
  statements before, 0.54 s and about 186 after. New ids come from a
  small `id_blocks` table (schema v45), so each row knows its parent's id
  before anything is written.
- **Several sector writers can run at once (PERF.14).** A sector reserves
  all its system and phenomenon names in one locked statement instead of
  a `SELECT ... FOR UPDATE` per name, which deadlocked when four
  `generate.py sector` runs started together. A sector save now runs at
  READ COMMITTED and is retried from the start on a deadlock or lock
  wait timeout.

### Fixed
- **A deadlock no longer leaves a half-saved transaction carrying on.**
  The connection pool used to re-run a statement that hit a deadlock on
  a fresh session, so the rest of the save continued in a transaction
  MySQL had already rolled back and failed later with a foreign key
  error (1452). The pool now raises the error instead, and the save
  retries cleanly.

Run `update.sh` (or `update.ps1`) after updating: it migrates the
database to schema v45.

## [7.84.177] - 2026-10-01

### Added

- Admin sign-in hardening (SEC.22 to SEC.25): a wrong current password on
  the account page now counts as a failed login and is limited like one;
  a successful login gives the browser a 90-day trusted-device cookie
  that keeps it past its username's lockout (the per-address lockout
  still applies); new passwords are checked against a bundled list of
  about 47,000 common and breached passwords and may not be just the
  username or the site's name with a few characters added; passwords are
  hashed with PBKDF2-SHA256 at 600,000 rounds and older hashes are
  upgraded at the next login. `src/loginLockouts.py --forget-devices
  <user>` revokes an admin's trusted devices.
- A fail2ban filter and jail for failed admin sign-ins
  (`examples/fail2ban/`, guide in `docs/deployment/fail2ban.md`; SEC.27).
- Optional two-factor sign-in for admins (SEC.26): turn it on from the
  account page by scanning a QR code with any authenticator app; sign-in
  then asks for the 6-digit code after the password (or one of ten
  single-use recovery codes). Wrong codes count as failed logins and are
  logged as `AUTH totp.failed` (matched by the fail2ban filter). API keys
  are unaffected. `src/loginLockouts.py --reset-two-factor <user>` turns
  it off for a lost phone. Control schema v4 (rerun
  `update.sh`).

## [7.83.177] - 2026-10-01

### Security
- **Per-address login lockout (SEC.1).** Three failed logins from one
  address lock it for 5 minutes; each further lockout doubles that, up
  to a day. Checked before the password, so a guess during a lock learns
  nothing. A successful login clears the address's count, but its
  doubling level only halves for each day without a lockout. IPv6
  addresses are counted by their /64. Loopback and the new
  `login_allowlist` in `config.json` are never locked. If a private
  address gets locked while `proxy_fix` is off, the app warns in the
  error log and on the admin Stats page, because behind a reverse proxy
  every visitor shares that address.
- **Lockouts are shared and survive restarts (SEC.21).** The
  per-username backoff (10 free failures, then 1 s doubling to 15
  minutes) and the new per-address lockout live in the control
  database's new `login_throttle` table (control schema v2) instead of
  each worker's memory. Until `update.sh` creates the table, the counts
  are kept in memory as before.
- **Lifting a lockout.** The admin Stats page lists current lockouts
  with Lift buttons (`GET /api/admin/lockouts`, `POST
  /api/admin/lockouts/lift`), and `python3 src/loginLockouts.py` lists or
  lifts them from a shell (`--ip`, `--user`, `--all`) for an admin locked
  out of the site itself.

## [7.82.177] - 2026-10-01

### Security
- **An always-on activity log (SEC.28).** planetGen now always writes a
  short record of who did what to `planetgen.log` in the platform's
  standard log folder (`/var/log/planetgen/` on Linux,
  `/Library/Logs/planetgen/` on macOS, `logs\` under the checkout on
  Windows; `"log_dir"` in `config.json` moves it): sign-ins, logouts,
  credential and API key changes, refused requests (missing, expired or
  forged sessions, bad API keys, admin-only actions, failed form tokens),
  every database write made through the web interface or API, each
  schema migration step, and each `generate.py` run or Generate page job
  as one start and one finish line with its counts. One fixed line format
  (`2026-10-01T08:00:00Z planetgen[1234]: AUTH login.failed ip=203.0.113.5 user="admin"`),
  documented in `docs/config.md`, with the address always before any
  text a visitor typed, which is quoted and escaped. It is separate from
  the debug log, which keeps working as before and gets a copy of each
  line while it is on. `install.sh`/`update.sh` create the folder and
  install logrotate (`/etc/logrotate.d/planetgen-log`) or newsyslog
  (`/etc/newsyslog.d/planetgen-log.conf`) rotation, daily or past 100 MB,
  30 copies; on Windows, or without those files, the program rotates the
  file itself.
- **Failed and locked logins are recorded with their address (SEC.20).**
  Each one is an `AUTH login.failed`/`login.locked` line in the activity
  log and a row in the control database's audit log (kept 90 days), as is
  a wrong current password on `/account`. The admin Stats page lists the
  newest ones under "Failed sign-ins" (new `GET /api/admin/login-failures`).
  A wrong password on the `/login` page now answers 401 instead of 200,
  so the web server's access log shows it too.

## [7.81.177] - 2026-10-01

### Changed

- Galaxy Map: once you've picked an arc (and inside every block), the view can be turned, moved and zoomed again: drag to turn it, right-drag or Shift-drag to move it, scroll or pinch to zoom, and a new Reset view button brings it back. The whole galaxy and its quarters stay locked top-down. Where the view can be turned, a layer can also be picked by clicking it on the map. Each step of the drill-down still opens on its own view.

## [7.80.177] - 2026-10-01

### Fixed

- Galaxy Map: stars in generated sectors now show in the drill-down. A
  generated sector's block is solid and hid the stars inside it; stars now
  draw over the blocks, and the faintest ones are a little brighter
  (MAP.51 follow-up).

## [7.79.177] - 2026-10-01

### Changed

- Every page, map panel and generated description now shows a number with
  5 or more digits before the decimal point in scientific notation, to 3
  significant figures: "384,400 km" reads "3.84 × 10⁵ km", "12,345
  systems" reads "1.23 × 10⁴ systems" (UX.20). Counts and measurements
  both follow it. Left as they were: IDs, seeds, years in dates, page
  numbers, ring and slot numbers, coordinates and galaxy positions, and
  the raw numbers in the API's JSON. One shared formatter does it
  (`stellarObjects.utils.format_number`, `num()` in the templates, and
  `static/numberformat.js` for the maps). No regenerate is needed: pages
  render descriptions from the stored numbers.

## [7.78.177] - 2026-10-01

### Changed

- Galaxy Map: once the drill-down is down to a small cube of sectors (a 3×3×3 block), it is shown at a fixed slant with its layers pulled apart, so any one of its sectors can be hovered and clicked straight on the map; the list beside the map still picks a layer.
- Galaxy Map: a block whose sectors are all in one layer (one sector tall, at the disk's top and bottom) shows its sectors right away instead of making you pick its smaller blocks first.

## [7.77.177] - 2026-10-01

### Fixed

- Galaxy Map: filled sectors now show their stars (MAP.51). Each tile also
  serves its generated systems' stars, down to a luminosity floor that drops
  as you zoom in (about 260 L☉ seen from kiloparsecs away, 1 L☉ at a few
  hundred parsecs, every star, red dwarfs included, at sector depth), so
  filled regions show up more and more as you close in. Every star, bright
  or generated, draws as a point of glowing light sized by the star's
  radius, brighter the more luminous, colored by its temperature; clicking
  a generated star shows it and links to its system.

## [7.76.177] - 2026-10-01

### Fixed

- System page: an asteroid belt's row in the object list shows only its
  density, its range ("2.1 AU to 3.3 AU") and its three largest minerals,
  no longer the belt's distance twice (UX.19). No row (planets, moons,
  belts) shows the zone any more; the zone stays in each body's own
  description.
- API: `GET /api/systems/<id>` gives each belt its `composition`
  (`{component, concentration}`, largest share first).

## [7.75.177] - 2026-10-01

### Fixed

- System page, "Place a facility" (ADM.9): the form asks for the name,
  then the placement, then the host, and the host list holds only the
  stars, planets, moons or belts that take that placement.
- An orbital facility's distance is now a logarithmic slider from just
  above the host's surface to the edge of its sphere of influence (a
  planet's or moon's Hill sphere, a star's heliosphere), with the
  distance, period and speed read out as it moves. It shows only for "in
  orbit"; a surface facility has no distance control. The API refuses an
  orbit outside the host's sphere of influence.
- A facility in an asteroid belt gets a random spot in the belt and a
  circular orbit around its star from there, and `updateOrbits.py` moves
  it along with orbital facilities. No schema change; belt facilities
  saved before this have no orbit and stay where they are.

## [7.74.177] - 2026-10-01

### Fixed

- The Galaxy Map's grid guide lines are now faint: the wedge lines (meridians) are mostly see-through and the block edges marking the ring, wedge and layer boundaries are fainter, just visible enough to follow. The bearing labels stay readable.
- Zoomed in to an arc, the wedge lines again stop just past the part in view instead of running across the whole map (a block reaching a little before the arc's first bearing made the clip wrap all the way round).

## [7.73.177] - 2026-10-01

### Fixed

- The Galaxy Map no longer has a free camera, rotation, zoom buttons, Slice or Free look (MAP.17). It is always seen from above: click a quarter of the galaxy, pick a layer of the disk from the list beside the map, click an arc of the ring band in view (about a third of it) to zoom in, and so on down to a sector (MAP.19, Boss's "Layer + arc"). Every pick is a big target.
- Hovering a quarter or arc dims everything else to 25% and outlines it (MAP.18); hovering a layer in the list dims the other layers.
- Zoomed in, the wedge lines stay within the part in view and a small margin past it (MAP.44).
- "Show on Galaxy Map" (`?sector=`) opens at the sector's own layer with the sector highlighted, and the map has its own Back, Forward, Up and Whole galaxy buttons (MAP.26). Stage links are now `/galaxy?at=…&p=…`; older `?slab=` links still open.

## [7.72.174] - 2026-10-01

### Fixed
- **Test suite green again.** The route fuzz test that tries ids like `0001` or `١` on every page now expects the "Show on Galaxy Map" links (`/sector/<id>/galaxy`, `/system/<id>/galaxy`) to answer with their redirect to the Galaxy Map instead of failing on it. The pages themselves were already correct.

## [7.71.174] - 2026-10-01

### Fixed
- **Names stay on screen on the System Map (MAP.50).** A planet's or moon's name near the edge of the map now slides back inside the frame (or takes a side that fits) instead of running off it, and a star's name does the same. The page also checks the real text width once the map is shown and pulls any name that still crosses an edge back in.

## [7.70.174] - 2026-10-01

### Fixed
- **Asteroid belts on the System Map no longer cover planet orbits (MAP.49).** A belt's ring now runs from its inner edge to its outer edge on the same scale as the orbits, instead of being centered on its inner edge with a width that ignored that scale. A very thin belt is still widened so it can be seen, but never over a neighboring orbit. A planet's orbit is drawn at its real distance from its star, so a planet on a tilted orbit just past a belt no longer looks like it sits inside it. Generated systems were checked too: no planet actually orbits inside a belt.

## [7.69.174] - 2026-10-01

### Fixed

- Sector Map: a neighboring sector's rogue planets, comets, black holes and
  neutron stars are no longer drawn in this sector, outside its wireframe
  (MAP.45). A neighbor's cloud that reaches in is still drawn, fainter, and
  says where it comes from. Stored positions were always correct; no
  regenerate is needed.
- Sector Map: rogue planets are easy to find (MAP.46): a brighter violet, a
  ring around each that stays visible zoomed out (with a "Mark rogue
  planets" button to hide the rings), and a "Show on map" button for each
  one in the sector's Contents.

## [7.68.174] - 2026-10-01

### Fixed
- **Bright stars on the Galaxy Map show at every zoom (MAP.47).** Zoomed out, each tile's bright-star query sorted every star in the tile (millions in a big tile) before picking the brightest, so a zoomed-out view ran past the 30-second API timeout and the stars vanished. The query now walks the luminosity index and stops at the tile's 400, and a tile too big (or too empty) to query on its own takes the brightest stars inside it from one galaxy-wide sample of the 100,000 most luminous, read once per request.
- **No more stars popping in after a zoom (MAP.48).** While a zoom's new tiles load, the stars already on screen stay, along with any a cached coarser tile holds there; stars new to the view fade in (not with reduced motion). The faster query above does most of the rest.
- **Wedge lines stop at the galaxy's edge (MAP.43).** They used to run 7% past it (the padded view radius plus 2%); now every wedge line, including the finer ones clipped to the view, ends at the outside of the outermost ring, and the bearing labels sit just past the ends.
- **Generated systems are easy to find on the Galaxy Map (MAP.37).** Blocks not yet filled are drawn much more transparent (10-30% opaque by density, was 50-80%), and any block holding a generated sector is at least 60% of the way to solid and painted a saturated amber (deeper on the light theme) instead of a faint warm white, in the free view and the drill-down stages alike.

## [7.67.174] - 2026-10-01

### Fixed
- **Buttons always have space between them (UX.16).** One shared rule in
  `static/style.css` gives every group of buttons on every page (`.btn`,
  `.btn-small`, `.starmap-btn`, plain `<button>`s, and side-by-side forms
  that each hold one) the same gap, across and between wrapped lines:
  `--btn-gap`, 0.5rem with a mouse and 0.75rem on touch screens
  (`pointer: coarse`), so 44-48 px touch targets never sit edge to edge.
  Groups that had their own smaller gap (the Galaxy Map's address
  matches, 0.3rem) now use it too.
- **An object's data sits beside its 3D render when there's room
  (UX.15).** On the System Map the info panel moves to the right of the
  map once the panel is at least 46rem wide (the map shrinks to fit and
  stays square); on a phenomenon page with a 3D view (neutron stars,
  black holes, quasars, rogue planets, interstellar comets) the data
  table sits beside the view from 50rem. Narrower screens keep the
  stacked layout. Both switches are container queries, so the layout is
  set before the render loads and nothing jumps.

## [7.66.174] - 2026-10-01

### Changed
- **Version numbers are now MAJOR.REVISION.BUILD** (OPS.1). A `major`
  release note bumps MAJOR, any other note bumps REVISION, and BUILD is
  the sum of the TODO category counters in the new "Next free IDs" table
  of `docs/design/todo-number-map.md`, so it tracks how many TODO items
  have ever been filed. `bump_version.py --check` fails when `docs/TODO.md`
  uses an ID that table hasn't counted. See `changes/README.md`.

## [7.65.174] - 2026-10-01

### Changed
- **README and INSTALL split.** `README.md` now says what planetGen is,
  what it does, what it needs and how to use the website and the command
  line. The new `INSTALL.md` walks from nothing to a running site on
  Linux, Windows or macOS with the provided scripts, says where every
  example config lives, and covers updates and scheduled maintenance.
  The old `pip install .` setup step is gone: the install scripts set up
  the libraries and everything runs the checkout's code. The full
  command-line reference moved to `docs/cli.md`.
- **TODO items have permanent category IDs** (`UX.1`, `MAP.16`, ...),
  a plain running count in each category like the schema version,
  replacing the running numbers. Bugs and features are listed under the
  item they belong to, `docs/TODO.md` is the one place that links an
  item to its design document, and code tags read `TODO(MAP.16)` (a test
  checks each names an open item). `docs/design/todo-number-map.md` maps
  every old number, by date, and the short-lived dotted IDs (`MAP.2.1`)
  to the new IDs, for the changelog, commits and PRs that cite them.
- **Every reference doc checked against the code.** `database-schema.md`
  now describes schema v44 and its tables; `api.md`,
  `html-interface.md`, `config.md`, `system-file-format.md`, the
  deployment guides and the rest have their errors fixed (for example,
  `update.sh` resets to the branch tip rather than refusing to run over
  local changes).

### Added
- `docs/design/architecture.md`: how the program fits together, mapping
  every script, package and module to what it holds and tracing the main
  flows with diagrams.
- `docs/design/design-decisions.md`: the big design choices, when they
  were made, why, and what was rejected. The design documents were
  brought up to date with their reasons; superseded ones moved to
  `docs/design/archive/` and `docs/analysis/archive/`.

## [7.64.174] - 2026-10-01

### Added
- **Show on Galaxy Map.** Sector pages, system pages and sector and system search results link to the sector on the Galaxy Map (`/galaxy?sector=<designation>`), which opens the drill-down stage holding it. The sector page's old Quadrant link stays as its own badge.
- **Pick on a map from the NAV page.** At each step the NAV page offers "Pick on Galaxy Map" and, once the other endpoint is known, "Pick in this sector", which opens the Sector Map's pick mode.

## [7.63.174] - 2026-10-01

### Added
- The Galaxy Map has a NAV pick mode (`/galaxy?pick=from|to`). It shows a "Choosing a start/destination · Cancel" banner and keeps "Generated only" on. Clicking a sector opens it in the Sector Map's pick mode.

## [7.62.174] - 2026-10-01

### Added
- **Species and polity pages.** `/species` lists every species (with a spacefaring filter), `/species/<id>` shows one, `/polities` lists every polity and `/polities/<id>` lists the systems it holds, nearest its capital first. A life world's planet row names its dominant species, and an owned system's page says whose territory it is in.
- **Hidden until there is population data.** The header's Species section, these pages, "Dominant species" and "Territory of ..." only appear once a population pass has made species (or polities); before that the pages are 404s.

## [7.61.174] - 2026-10-01

### Added
- The Galaxy Map gives an admin "Generate this block" on a 3-sector block and "Generate this layer" on one of its layers (drill-down stages 7 and 8). Each button starts the Generate page's block job.

## [7.60.174] - 2026-10-01

### Added
- **Generate a Galaxy Map block.** `generate.py galaxy --block M.I.S.SLAB [--block-layer J]` fills one drill-down block (or one of its layers), and the admin Generate page has a matching form.
- **Neighborhood radius in light-years.** The Generate page's single-sector neighborhood takes a radius in ly (13 ly and up), converted to parsecs.
- **JSON answers from the Generate page.** A post with `Accept: application/json` returns the started job's id and status URL, so the Galaxy Map can start jobs without leaving the map.

## [7.59.174] - 2026-10-01

### Added
- **More on the bright stars.** The Stats page counts the pre-placed bright stars exactly: how many were placed, how many are built into systems and how many are still waiting for their sectors. The Generate page says whether the bright stars have been scattered, at what brightness and with what seed. On the Sector Map, a neighboring sector that hasn't been generated yet lists the bright stars waiting in it. New `GET /api/galaxy/bright-stars` lists one sector cell's bright stars.

## [7.58.2] - 2026-10-01

### Changed

The population pass (species, civilizations, territories) is now optional and off by default. `generate.py sector` and `generate.py galaxy` run it after saving only with the new `--population` flag (this replaces `--no-population`). `install.sh` and `update.sh` (and `install.ps1`/`update.ps1`) ask whether to run it after the database step, defaulting to No after 30 seconds and skipping it with no terminal; `POPULATION=1` (`-Population` on Windows) runs it without asking. `generate.py population` still runs it by hand.

### Added

`GET /api/population` (and `population.population_status`) says whether a population pass has run and whether any species, polity or owned system exists, so the pages and the Galaxy Map's Territories button can hide themselves when there is no population data.

## [7.58.1] - 2026-10-01

### Changed
- The Galaxy Map leaves out its Territories button and legend until population has made at least one polity.

## [7.58.0] - 2026-10-01

### Added
- **Pick NAV endpoints on the Sector Map.** Clicking a system or phenomenon on the Sector Map now offers "Nav from here" and "Nav to here". Opening a sector with `?pick=to&from=system:12` (or `pick=from&to=...`) shows a "Choosing a destination" banner with Cancel, and each system or phenomenon offers "Use as destination" (or start), which goes straight to the plotted course.

## [7.57.0] - 2026-10-01

### Added
- Phenomenon pages show a view that suits the object: neutron stars spin with their radio beams in rough time with their real spin period (slowed, and the caption says by how much), black holes and quasars show 3D accretion disks (quasars add a dusty torus, and jets when radio-loud), rogue planets and interstellar comets are rendered bodies (a coma and tail for an active comet), and asteroid fields have no view. Nebulae and supernova remnants keep the AU-scale diagram. With reduced motion a still frame is drawn, and without JavaScript a simple drawing shows.

## [7.56.0] - 2026-10-01

### Added
- **Faster pages.** The web pages now keep the answers they get from the API in memory, so a repeat visit to a sector, system or list page doesn't query the database again. Any edit made through the site or the API clears it at once, sectors and systems added by generation jobs are noticed within 15 seconds, and nothing is kept longer than 5 minutes. It can be tuned or turned off with `page_cache` in `config.json` (see `docs/config.md`).

## [7.55.0] - 2026-10-01

### Changed
- **A bigger Galaxy Map.** The map now spans the full width of a wider page and as much of the window's height as fits under the header, with its buttons in a row underneath (wrapping on narrow screens), the block details beside them, and the how-to text below. On a phone it stays about square so the page still scrolls past it.

## [7.54.0] - 2026-10-01

### Added

- A Territories button on the Galaxy Map shows who holds what: each polity's reach as a soft ball of its own color around its capital, and the systems it owns as dots in that color, with a list under the map naming each polity, its government and how many systems it holds. It draws the same way at every zoom, so turning it on over the whole galaxy costs nothing.

## [7.53.0] - 2026-10-01

### Changed

- Generate neighborhood now asks how far it should reach. The buttons on an ungenerated sector (on the Galaxy Map and the Sector Map, for an admin) take a radius in light years, between 13 and 652, and say how many sectors that covers as the number changes; past about 5,000 sectors it asks before starting. It used to be a fixed 100 light years with no warning.

## [7.52.0] - 2026-10-01

### Added

- A plotted course can be seen on the Galaxy Map: the NAV result offers "Show on Galaxy Map", which opens the smallest view holding both ends and draws the course through its stops, each ringed and the two ends named. A course that stays inside one sector opens that sector instead.

## [7.51.0] - 2026-10-01

### Added
- Bright-star queries for the web pages: `adminStats.bright_star_counts` (placed, filled and unfilled pre-placed bright stars, also in `GET /api/admin/stats` as `database.bright_stars`), `queryDb.bright_stars_in_sector` (a cell's bright stars, unfilled only by default) and an `unfilled_only` option on `queryDb.galaxy_bright_stars_in_box`.
- `GET /api/galaxy/shape` now returns `bright_stars`: whether the bright-star scatter has run, its threshold and seed, and the default threshold (`queryDb.bright_star_scatter_status`).

## [7.50.0] - 2026-10-01

### Added

- The Galaxy Map takes an address: a field over its path accepts a sector designation, `ring/layer/slot` (or `ring 312 layer -3 slot 1042`), `x, y, z` in parsecs, or a sector or star system name, and flies to that sector. A name with several matches lists them to pick from, and anything that can't be a sector says why.

## [7.49.0] - 2026-10-01

### Added

Worlds with complex life now have a named dominant species, and technological civilizations have an age and an era (Industrial through Elder). Every spacefaring species founds one polity that claims the generated systems around its homeworld, out to a reach that grows with its age (up to 100 ly), so groups of systems form territories in 3D. `generate.py population` builds all of this from the stored galaxy with no regenerate; `generate.py sector`/`galaxy --population` also run it after saving. New read endpoints: `/api/species`, `/api/polities`, `/api/systems/<id>/owner`, `/api/planets/<id>/species` and `/api/territories`. Schema v44 adds the `species`, `polities`, `system_owners` and `population_state` tables; run `update.sh` (or `migrateDb.py`), then `generate.py population` when you want it. Design: docs/design/population-and-politics.md (TODO 51-54).

## [7.48.0] - 2026-10-01

### Added
- Sector Contents rows for phenomena show their octant and their three nearest star systems, linked. Systems list their nearest systems across sector boundaries, and phenomenon pages show their octant and nearest systems. `GET /api/phenomena/<type>/<id>` gains `nearest`.

## [7.47.0] - 2026-10-01

### Added
- **Place facilities from the system page.** Admins can add a starbase, station, outpost, colony or mining colony to a star system: pick the star, planet, moon or asteroid belt it goes on or around, preview the orbit's distance, period and speed (worked out from the host's mass, like every other orbit) and any placement rule it breaks, then save it. Each facility can be removed again.
- Facilities now show on the system page (in their own panel and in their host's row), as small diamonds on the System Map, and, for stand-alone ones and those on asteroid fields, in the sector page's Contents. A colony makes its world show as Inhabited.

## [7.46.0] - 2026-10-01

### Added
- **Class reference pages.** A new Classes section (`/classes`) lists every kind of class the generator gives out: star spectral and luminosity classes, planets, nebulae, supernova remnants, asteroid fields, black holes, rogue planets and comets. Each type has a page listing its classes, and each class has its own page of facts, all read from the generator's own tables when the site starts, so they always match what it generates.
- Class labels now link to these pages: a star's type and a planet's or comet's class on the system page, and the class of a nebula, supernova remnant, asteroid field, black hole or rogue planet on its phenomenon page.

## [7.45.0] - 2026-10-01

### Added
- **Bright stars on the Generate page.** Planning a galaxy (on its own or as part of a new galaxy) now shows the bright-star scatter as its own step, with its progress bar and the number of stars it placed. A "Skip the bright-star scatter" box leaves it out, and a new "Rebuild the bright stars" form scatters them again on the current plan, optionally leaving already generated sectors out.
- The admin Stats page shows about how many bright stars have been placed.

### Changed
- **Sector Map stars stay small and glow by brightness.** Every star is now a small point, and how bright it is shows in its glow: a supergiant has a big soft halo, the Sun a modest one, and a white dwarf is a tiny dot with almost none. Stars are still easy to click.

## [7.44.0] - 2026-10-01

### Added

- The Galaxy Map opens on the drill-down: the whole galaxy in blocks 243 sectors a side, where hovering a slab highlights it and clicking pulls it out to a view from above; clicking a block there flies into it and shows its contents as blocks a ninth the size, down to single sectors, where a click opens the sector (or, for a sector not generated yet, shows where it is and its Generate buttons for an admin). A breadcrumb with sibling menus, a slab list with generated counts, a hover tooltip, keys (arrows, Enter, Escape, Home), touch taps and a "Generated only" toggle come with it. Each stage has its own URL (`/galaxy?slab=`, `?at=`, `?sector=<designation>`), so Back and Forward work and a stage can be linked. The old free camera stays behind a Free look button.

## [7.43.0] - 2026-10-01

### Added

- **The API can add a system to an existing sector and regenerate a
  system in place.** `POST /api/systems` takes an optional `sector_id`
  (and `position`): the new system is placed clear of the sector's other
  systems' Hill spheres, with its location, containment and nearest
  systems filled in. `PATCH /api/systems/<id>` takes `{"regenerate":
  recipe}` to replace a system's stars, planets, moons, belts and comets
  while keeping its id, name, place and links.

## [7.42.1] - 2026-10-01

### Changed
- **Lighter Galaxy Map meshes where sectors are generated.** Blocks whose sectors are all generated (drawn solid) no longer draw the faces they share with each other, which nobody can see. A fully generated neighborhood now needs about a twentieth of the vertices it did, which matters most on phones. Translucent blocks keep every face, since those faces draw the block grid you see through the glass.

## [7.42.0] - 2026-10-01

### Added

The Galaxy Map draws the pre-placed bright stars (every star of 500 L☉ or more, placed by `generate.py plan` before any sector is filled), so the spiral arms show before anything is generated. Each star is a tiny point with a big soft glow in its own color, the same few pixels across at every zoom. Clicking one shows its type, luminosity and sector, and links to its system once that sector is filled. `/api/galaxy/tiles` lists each tile's most luminous 400 as `stars`.

## [7.41.3] - 2026-10-01

### Added
- **The drill-down's stage API.** `GET /api/galaxy/stage?at=m.ring.wedge.slab`
  (and the site's cached `/galaxy/stage`) returns how many generated
  sectors each block inside a drill-down block holds, and the sectors
  themselves at the smallest level. It feeds the Galaxy Map's coming
  drill-down navigation; nothing on the map changes yet.

## [7.41.2] - 2026-10-01

### Added
- **The drill-down's block ladder.** The Galaxy Map's coming drill-down
  navigation has its geometry: blocks 243, 27 and 3 sectors a side, each
  sitting wholly inside one block of the next size up, with the same
  rules on the page (`galaxyprisms.js`) and the server
  (`planetgen/galaxy/drill.py`). Nothing on the map changes yet.

## [7.41.1] - 2026-10-01

### Changed
- The Sector Map draws nebulae and supernova remnants as see-through volumes: densest through the middle and fading at the edge, with a remnant showing as a bright shell. A cloud far larger than the sector, or one centered in another sector, still tints the view from inside it, and stars inside a cloud stay visible and clickable.

## [7.41.0] - 2026-10-01

### Added
- System and phenomenon pages show an "Inside" badge linking the nebula or supernova remnant they sit in, and the sector Contents list says "Inside <name>" for those systems. `GET /api/systems/<id>` gains `inside`.
- Sector Contents rows show a nebula's, remnant's or asteroid field's class.

## [7.40.1] - 2026-10-01

### Changed

- **Bright stars can go down to 100 solar luminosities.** The default
  stays at 500 (about 60 million stars in a Milky Way, about 10 GB);
  `generate.py plan --bright-star-min-luminosity 100` now works too
  (about 220 million stars, about 35 GB). White dwarfs, which are never
  pre-placed, always stay in a sector's own draw.

## [7.40.0] - 2026-10-01

### Changed

- **Bright stars down to 100 solar luminosities by default.** `generate.py
  plan` now pre-places every star of 100 solar luminosities or more
  (about 220 million in a Milky Way, roughly an hour and a quarter plus
  the database load) instead of 500. `--bright-star-min-luminosity 500`
  still gives a quick test galaxy. The hottest white dwarfs, which sit
  right at 100, are never pre-placed and now always stay in a sector's
  own draw.

## [7.39.0] - 2026-10-01

### Added

- **A nebula squeezes the heliosphere of every system inside it.** The
  gas of a nebula or supernova remnant pushes on a star's wind bubble
  much harder than open space does, so a system inside a dense cloud
  now reports a much smaller heliosphere (a Sun-like star's drops from
  about 85 AU to well under 1 AU in a dense cold cloud). The system
  page text says so, and the system API returns `inside`,
  `heliopause_au` and `heliopause_open_space_au`; navigation uses the
  squeezed heliopause as the edge of a system's local frame.

## [7.38.0] - 2026-10-01

### Added

- **Every bright star is placed before its sector is filled (schema
  v43).** `generate.py plan` now ends by drawing every star of 500 solar
  luminosities or more across the whole galaxy and storing each at a
  fixed point in its sector, in a new `bright_stars` table, while the
  sectors themselves stay unfilled. Filling a sector builds a full system
  around each of its bright stars first and draws the rest from dimmer
  stars, so a sector's expected count is unchanged. New plan options:
  `--bright-star-min-luminosity`, `--no-bright-stars`,
  `--bright-stars-only` and `--force`.
- **Star ages follow where a sector sits.** Each system in a galaxy
  sector draws its star from the young, intermediate, old or bulge
  population in proportion to their density there, so O and B stars and
  supergiants gather in the spiral arms near the plane.

## [7.37.0] - 2026-10-01

### Changed

- **The correlative update moves everything.** `updateOrbits.py` now
  turns every star system, phenomenon and stand-alone facility along its
  galactic orbit, moves anything that drifts into another generated
  sector over to it (sector, position, octant and location text), and
  then recomputes containment and the stored nearest systems. Orbital
  facilities advance around their hosts like moons.

## [7.36.0] - 2026-10-01

### Added
- **Nebulae and supernova remnants on the Galaxy Map.** They are drawn as
  soft translucent clouds their real size, in the Sector Map's colors,
  and fade out as the camera gets close. Clicking one shows its type,
  class and radius, with a link to its page.

## [7.35.0] - 2026-10-01

### Added

- **Schema v42: starbases, colonies and outposts.** Facilities can stand
  on a planet or moon, orbit a star, planet or moon, sit in an asteroid
  belt or field, or park in open space, following Boss's placement rules
  (gas giants take orbital facilities only). An orbital facility's
  period and speed come from its host's mass. The API adds, lists,
  previews and removes them (`/api/facilities`).

## [7.34.0] - 2026-10-01

### Changed
- **The Galaxy Map has Generate buttons.** Clicking a single sector that
  isn't generated yet now gives a logged-in admin the same four buttons as
  the Sector Map: generate this sector, its neighborhood, its column, or
  its whole shell (after a confirm). They start a job on the Generate
  page. Visitors see the sector's address and designation only. The
  "Copy CLI command" button is gone.

## [7.33.0] - 2026-10-01

### Added

- **Schema v41: octants and nearest systems.** Every placed phenomenon
  records the sector octant it sits in, and the database stores the 3
  nearest star systems to every placed system and phenomenon, found
  across sector boundaries. The sector data now returns both, ready for
  the pages to show.

## [7.32.0] - 2026-10-01

### Changed
- **The Galaxy Map zooms smoothly.** The blocks for each view are now built
  in a Web Worker (`static/galaxyblocks.js`), so a zoom step no longer
  freezes the page while hundreds of milliseconds of block listing runs.
  Zoom steps glide over 160 ms instead of jumping (they still jump with
  "reduce motion" turned on), and a change of block size crossfades
  instead of popping.
- **Zoom steps you've already seen are instant.** Built views are kept (up
  to about 4 million vertices), and while the map is idle it prepares the
  views and fetches the tiles one zoom step in and out. The first frame is
  still built on the page, and so is everything in a browser where the
  worker can't start.

## [7.31.0] - 2026-10-01

### Changed

- **Schema v40: names that follow one standard.** Nebulae, supernova
  remnants, neutron stars, black holes, quasars and rogue planets are
  named the way star systems are, so no two objects in the galaxy share
  a name. Comets and asteroid fields get designations: `P/<star>-<n>`
  for a periodic comet, `C/<star>-<n>` for a long-period one,
  `I/<sector>-<n>` for an interstellar comet and
  `AF <class>-<sector>-<nn>` for an asteroid field. The migration
  renames existing rows.

## [7.30.0] - 2026-10-01

### Added

- **Nebulae and supernova remnants are generated with the stars they
  need.** Sectors now make molecular clouds (dark classes M-Q) at the
  research rate. Every planetary nebula comes with its own new hot white
  dwarf system at its center. Every O star, and half the B0-B2 stars,
  sits in an H II region (classes C-E or G), and a few later B and A
  stars light a reflection nebula (`NEBULA_HOST_RULES`). A core-collapse
  remnant's neutron star or black hole has drifted off-center by its
  birth kick (`SUPERNOVA_KICK_SPEED_RANGE_KMS`) times the remnant's age.
  Diffuse gas (classes A-B) is the background and isn't generated.
- The sector's nearby phenomena now carry each nebula, supernova remnant and asteroid field's letter `class`.

## [7.29.2] - 2026-10-01

### Changed

- **Sector Map star dots are sized on a log scale.** White dwarfs, red dwarfs, the Sun, giants and supergiants now draw at visibly different sizes; before, every star past about 8 solar radii hit the same 14 px cap. White dwarf and giant systems are now covered by tests on the system page and System Map.

## [7.29.1] - 2026-10-01

### Changed

- The Sector Map draws the sector's own cell again as a faint wireframe, with its inner and outer faces following the ring's curve.

## [7.29.0] - 2026-10-01

### Added

- **Search phenomena.** The search page's Browse by Tag gains Phenomenon (nebula, black hole, rogue planet and the rest) and Phenomenon Class (for example "D: Classical H II region" or "C3 asteroid field") groups, and a Phenomena results table that links to each one. `GET /api/search` takes `phenomenon` and `phenomenon_class` (written `<type>:<class>`) and returns a `phenomena` result panel.

## [7.28.0] - 2026-10-01

### Added

- **Phenomenon classes on their pages.** A nebula or supernova remnant page shows its class letter and name (for example "D: Classical H II region") and what it holds: dominant species, gas density, temperature and extinction. An asteroid field shows its class (for example "C3") and composition family. The Sector Map colors the new diffuse nebulae a pale pink.

## [7.27.0] - 2026-10-01

### Added

- **Generate a column or a shell, and generate from the Sector Map.** `generate.py galaxy` gains `--ring I --slot K --column` (one slot through every layer the galaxy reaches) and `--ring I --shell` (a whole ring through every layer; needs `--limit` or `--yes` past 2,000 sectors), and `--ring --slot` now takes `--radius-pc` to generate that sector's neighborhood too. The admin Generate page offers the column and shell modes and the neighborhood radius. On the Sector Map, a logged-in admin who clicks a not-yet-generated neighbor gets Generate this sector, Generate neighborhood, Generate column, and Generate the entire shell (not recommended) buttons that start the job; visitors see only its address and designation, with no command line.

## [7.26.0] - 2026-10-01

### Added

- **New database fields on the web pages.** A black hole's page shows its Class (stellar, intermediate or supermassive) and a rogue planet's its Mass Class. A runaway or hypervelocity star shows its speed as a badge on its system page and in its sector's Contents, and `GET /api/systems/<id>` returns `runaway_class` and `runaway_speed_kms`. The sector page adds its star count and the estimated number of interstellar comets and planetesimals drifting through it, and folds two or more rogue planets into one expandable Contents row.

## [7.25.0] - 2026-10-01

### Changed

- The Galaxy Map shows generated and not-yet-generated sectors through its blocks alone; the marker dots are gone. Unfilled space is see-through (more so where it's sparse), and a block grows more solid and warmer the more of its sectors are generated, fully solid once all are, so generated sectors can be found by zooming at every scale. Zoomed in to single sectors, a generated sector takes its real density's color and links to its page.
- Clicking a block shows its ring, layer and slot ranges, its exact sector count and how many are generated; a click prefers the nearest block holding generated sectors along the line of sight, and centers on them, so a double-click zooms toward them.
- `/api/galaxy/tiles` tiles carry a `filled` summary counting every placed sector in the tile (per sector, or per cell for large or crowded tiles), which the map sums into its blocks.

## [7.24.0] - 2026-10-01

### Added

- **Schema v39: what sits inside a nebula.** Star systems, rogue planets,
  interstellar comets, black holes, neutron stars, asteroid fields and
  nebulae now record the innermost nebula or supernova remnant they sit
  inside (`inside_nebula_id` / `inside_remnant_id`). This is set by a 3D
  distance test when a sector is generated and whenever a nebula or
  remnant is placed, so sectors generated later inside an existing cloud
  see it. The migration fills it for existing data, and the sector detail
  query reports each system's cloud.

## [7.23.0] - 2026-10-01

### Changed
- **Random stars are now drawn as physics: a mass, an age, then how the
  star has evolved.** Masses come from the Kroupa (2001) initial mass
  function and ages from a 0-10 Gy star-formation history; each star is
  main sequence, subgiant, giant, supergiant or white dwarf according to
  how far its age is into its own lifetime, and its spectral letter and
  Yerkes class follow from its temperature and luminosity. Every O, B and
  A star used to be a supergiant or subgiant and no white dwarfs or red
  giants were ever made; now about 74% of stars are M dwarfs, 14% K,
  5% white dwarfs, a quarter of a percent giants, and supergiants are
  about one in a million, as in the real galaxy. `+large_star` uses the
  same model above 1.4 Msun. A specified `--star-type` draws its
  luminosity log-uniformly within its class instead of linearly.
- **A binary's secondary is born with its primary.** It takes a mass
  ratio of 0.1-1.0 of the primary's initial mass and the primary's age,
  and evolves by the same model, so it's never brighter or hotter than
  its mass allows.
- **Planets follow their star's age and history.** A star younger than
  10 Myr keeps belts only; no habitable-class world or moon forms around
  a star younger than 0.1 Gy; a giant has engulfed everything inside
  twice its radius, and a white dwarf's progenitor everything inside
  1.5 AU. A required habitable world gets a star old enough to host one.

### Added
- **Stellar populations by position.** `galaxyDensity.population_densities`
  splits a point's star density into young (0-0.1 Gy, a thin disk that
  crowds the spiral arms), intermediate (0.1-3 Gy), old (3-10 Gy, the
  thick disk) and bulge (8-12 Gy) stars, summing to the same total as
  before. A system config's new `POPULATION` draws its star's age from
  one, and `MAX_STAR_LUMINOSITY_SOL` keeps it dimmer than a threshold.
- **Bright-star sampling for pre-placement.** `stellarPopulation` gives
  the share of a population's stars at or above a luminosity
  (`bright_star_fraction`), draws stars conditional on being that bright
  (`sample_bright_stars`, exact and about 10 microseconds a star) or
  dimmer (`sample_dim_star`), and `Star.from_params` /
  `StarSystem(primary_star_params=...)` rebuild a stored star without
  re-rolling it. `SpaceSector.add_preplaced_system` places one at its
  stored position before the rest of the sector fills around it.

## [7.22.1] - 2026-10-01

### Changed
- **Search's tag groups fold away.** Each "Browse by Tag" group is a
  collapsible section showing how many tags it has; a group with a
  selected tag starts open (and says how many are selected), the rest
  start closed. Works without script.

## [7.22.0] - 2026-10-01

### Added

- The System Map's Measure distance now draws its path on the map and routes it around every planet, moon and star in the way, keeping a wide berth from stars and never threading between the two stars of a close binary. The result shows how much longer the route is than the straight line.

## [7.21.1] - 2026-10-01

### Fixed

- System Map names no longer overlap: the browser measures each name once a view is shown and moves or hides any that would collide. A hidden name still shows when its marker is hovered or focused.

## [7.21.0] - 2026-10-01

### Changed
- **Times show in the viewer's own time zone.** Pages write every time
  as UTC in a `<time>` element (labelled "UTC", so they read correctly
  without script), and the new `static/localtime.js` rewrites each in the
  browser's zone with its abbreviation: API key created/last used/revoked
  times, the Stats page's activity times, and a Generate job's start
  time. The database connection's session zone is now pinned to UTC, so
  `TIMESTAMP` columns read back the same whatever the server's own zone
  is, and the API's key times and the stats times end in `Z`.

## [7.20.0] - 2026-10-01

### Added
- **The Systems page lists every system.** A new All Systems table pages
  through every system 50 at a time, with its sector (linked) and
  octant; standalone systems keep their own table below, each paging on
  its own. `GET /api/systems` rows now carry `sector_name` and
  `quadrant` too.

## [7.19.1] - 2026-09-30

### Fixed
- **A binary's star word could be two words** (from the base name
  "El Nath"), so the star's name no longer ended in its own word and a
  naming test failed at random. Star words with a space are now redrawn.

## [7.19.0] - 2026-09-30

### Added

- **Schema v38: letter classes for nebulae, supernova remnants and
  asteroid fields.** Nebulae are classed A-Q and supernova remnants R-W
  (`program_constants.NEBULA_CLASSES`), each with what fills it: dominant
  species, particle density, temperature and optical extinction. A
  remnant's class follows its progenitor and core (a Type Ia remnant is
  always W; a pulsar wind nebula needs a pulsar). Asteroid fields get a
  class made of a composition-and-density letter and a size digit, like
  `C3`, with composition drawn from real asteroid families
  (`ASTEROID_FIELD_COMPOSITIONS`). Nebulae gain a `diffuse` family.
  Existing rows get the class they most likely are.

## [7.18.1] - 2026-09-30

### Fixed
- **Moons could orbit outside their planet's Hill sphere, or inside the
  planet itself.** Moons now orbit between the planet's surface (plus
  room for the largest moon it can hold and its atmosphere) and the
  prograde stability limit of about half the Hill radius (Domingos,
  Winter & Yokoyama 2006). A planet whose class is regenerated after a
  move drops the moons that no longer fit (TODO items 40 and 41).

## [7.18.0] - 2026-09-30

### Changed

- **Schema v37: interstellar objects at real-world rates.** Each sector
  now draws its phenomena per star from the research densities in
  `docs/design/interstellar-object-rates.md`
  (`program_constants.PHENOMENON_DENSITY_PC3`, with a
  `PHENOMENON_RATE_SCALE` dial per type). Isolated asteroid fields are no
  longer generated (they disperse), and planetary nebulae now appear.
- Rogue planets are drawn from four mass bins (terrestrial, sub-Neptune,
  Saturn-class, Jupiter-mass), so terrestrial rogues are now the most
  common kind instead of almost never appearing.

### Added

- Free-floating brown dwarfs, stored as rogue planets with
  `rogue_planets.mass_bin = 'brown-dwarf'`; the migration fills `mass_bin`
  for existing rogues from their mass.
- Runaway and hypervelocity stars: `star_systems.runaway_class` and
  `runaway_speed_kms`, flagged on ordinary systems (hypervelocity stars
  grow rarer with distance from the galactic center).
- `queryDb.sector_detail` reports each sector's star count and its
  estimated count of interstellar comets and planetesimals.

## [7.17.0] - 2026-09-30

### Changed

- The Galaxy Map sizes its sector blocks from the screen: each block is the smallest power-of-3 cube of whole sectors that is at least 4 pixels across at the focus, and blocks line up with the sector grid's master wedges wherever that keeps them near one block long. Only the solid's visible surface is built, so views zoom in to single sectors much sooner, and a clicked block shows its exact sector count, leaving out sectors the galaxy's outline doesn't allow.
- The Galaxy Map's wedge lines now follow every master wedge (3 from the core, doubling outward), each zone shown once its lines are far enough apart on screen.

## [7.16.1] - 2026-09-30

### Fixed
- **A planet's Hill sphere could overlap the asteroid belt inside it.**
  A planet after a belt now keeps 5 Hill radii clear of the belt's outer
  edge, the same rule a belt after a planet already followed, and a
  planet moved to fix spacing now gets its Hill radius recomputed for
  its new distance (TODO item 43).

## [7.16.0] - 2026-09-30

### Added
- **Windows installers.** `install.ps1` and `update.ps1` do what
  `install.sh` and `update.sh` do, with the Windows guide's layout: a venv
  with the libraries and waitress from `requirements-server.lock`
  (checked by hash), the NLTK corpus with `NLTK_DATA` set machine-wide,
  `config.json` from the Windows example when there is none, the
  migrate-or-delete question (y/N, 30 seconds), the runtime folders, and
  `icacls` permissions for the app's account. `update.ps1` pulls and
  installs only what's missing. `examples/maintenance/install-maintenance-task.ps1`
  schedules the monthly orbit update and `update.ps1` with Task
  Scheduler.
- **The bash installers run on macOS.** `install.sh`, `update.sh`,
  `install-maintenance-timer.sh` and their helpers run under macOS's
  bash 3.2: a venv from Homebrew's `python3` with gunicorn, `_www`
  ownership, `newsyslog` for the debug log, and gunicorn, the orbit update
  and `update.sh` as launchd daemons (`examples/macos/org.planetgen.update.plist`
  is new). `install.sh --skip-database` leaves out the database step.
- **`requirements-server.lock`** and a `server` extra in `setup.py`
  (gunicorn, or waitress on Windows) for those venvs.
- CI runs `install.sh` on macOS and `install.ps1` on Windows.

## [7.15.0] - 2026-09-30

### Added

- **Schema v36: a supermassive black hole at every galaxy's center.** When
  the nucleus roll finds no quasar (90% of galaxies), a quiescent
  supermassive black hole (1e6-1e8 solar masses, like Sagittarius A*) is
  placed at the galactic center instead. It has no galactic orbit, a faint
  accretion flow far below its Eddington limit, and a sphere of influence
  of `G*M/sigma^2`.
- `black_holes.mass_class` (`stellar`, `intermediate`, `supermassive`);
  the migration fills it from each row's mass.

### Changed

- Intermediate-mass black holes now span 1e2-1e5 solar masses
  (log-uniform) instead of 100-1,000, still 2% of black holes
  (`BLACK_HOLE_INTERMEDIATE_MASS_CHANCE`), and their text says they grew
  in a dense star cluster rather than in one supernova.

## [7.14.1] - 2026-09-30

### Security
- **Failed logins are counted per username, not only per address.** After
  10 failures for one username, from any mix of addresses, each further
  failure locks that username for twice as long as the last (1 s, 2 s,
  4 s, ... up to 15 minutes). A locked login is refused with a 429 and
  `Retry-After` before the password is checked, and the login page says
  how long to wait. Unknown usernames are counted the same way, so the
  lock doesn't reveal which usernames exist; a successful login clears
  the count (`src/html/api/loginbackoff.py`).

## [7.14.0] - 2026-09-30

### Changed

- NAV courses read "bearing mark mark" (for example "045 mark 330") on Boss's nested reference frames: bearing 000 points at the sector's center for a course inside one sector and at the galactic core between sectors, and the mark is the elevation mod 360 (270 is straight down). The NAV page shows the course and its frame in place of the Azimuth and Altitude rows, and the NAV Map's compass arrow marks bearing 000.
- `/api/nav`'s `direct` carries `bearing_deg`, `mark_deg`, `elevation_deg` and `frame` instead of `azimuth_deg` and `altitude_deg`.
- The NAV page's distances and the NAV Map's scale bar use the shared distance ladder (for example "1.07 pc (3.5 ly)").

## [7.13.1] - 2026-09-30

### Fixed
- **A close binary's planets could orbit inside the binary.** A close
  pair's innermost planet or belt now starts at the Holman & Wiegert
  (1999) circumbinary stability limit, about twice the stars'
  separation. When a habitable world is required and the pair's whole
  habitable zone lies inside that limit, the pair is made wide instead,
  unless the binary type was forced (TODO item 42).

## [7.13.0] - 2026-09-30

### Changed

- **Schema v35: hybrid master-wedge sector slots.** Each ring now holds
  the multiple of its master wedge count nearest `2*pi*(i + 1/2)`
  (3, 9, 15, 21, 27, 36, ...), with 3 master wedges at the center that
  double (6, 12, ... 1,536) once each would hold 8 slots. Slot
  boundaries line up on the master lines from the center to the edge,
  sector arcs stay within about 6% of 4 pc, and the total sector count
  is unchanged. `galaxyGeometry.ring_master_count` is new, and
  `galaxyprisms.js` mirrors both functions.
- The migration deletes every sector in a ring whose slot count changed
  (all but 15 rings) with its systems and phenomena. Run `update.sh`,
  then regenerate the galaxy; the skeleton (`generate.py plan`) is kept.

## [7.12.0] - 2026-09-30

### Changed
- **The system page lists everything in orbit once, in order.** The
  Planets & Moons, Asteroid Belts and Comets tables are gone; each list
  row now shows what they did (class, zone, distance, period, gravity;
  a belt's distance and inner-to-outer range; a comet's perihelion and
  period). Comets sort in among the planets and belts by semi-major axis
  (unbound ones last), and the Stars table sits first, right under the
  System Map.
- **One type chip per world:** "Habitable", "Terrestrial" or "Gas
  Giant", never two, plus a "Habitable moon" chip when one of its moons
  is habitable and "Inhabited" when it is.
- **A planet's moons get their own collapsed "N moons of ..." group**
  right under the planet's row, instead of sitting after its whole
  description.

## [7.11.0] - 2026-09-30

### Changed
- **A lighter top bar with a settings gear.** The theme switch, search
  (when the bar has no room for it) and the account links moved under a
  gear at the upper right: Account, Admin, Generate, Stats and Logout
  for an admin, Login for a visitor (who never sees Stats). Galaxy,
  Sectors, Systems, Phenomena and Nav stay buttons while they fit and
  fold into a Menu when they don't, and the bar's search box shows only
  while its text entry is at least twice the width of its button. The
  layout follows the header's own width with container queries, so it
  needs no script; with script, an open menu closes on an outside click
  or Escape. The logged-in header no longer collapses at a wider width
  than a visitor's.

## [7.10.3] - 2026-09-30

### Changed

- Galaxy Map: the scale readout and the cell info panel show distances on the
  shared distance ladder (`static/distance.js`), so parsec values carry ly in
  parentheses and large or small ones switch to kpc, cpc or ly like every
  other page.

## [7.10.2] - 2026-09-30

### Security
- **Admin generation inputs have upper bounds.** The Generate page, the
  one-off system page, the API and `generate.py` now reject a radius over
  200 pc (about 652 ly for the API's generate-neighborhood `radius_ly`,
  down from 1,000,000 ly), a ring or highest ring over 100,000, a ring
  `--limit` over 628,322 (ring 100,000's slot count), and more than 500
  orbital slots (also in a `--system-file`). All four share the constants
  in `src/planetgen/generation/limits.py`, and the page's number
  inputs carry them as `max`.

## [7.10.1] - 2026-09-30

### Changed

- Renames check only uniquely named objects (sectors, systems and stars)
  for a clashing name. Planet and moon names come from their star's, so
  `PATCH /api/planets|moons|stars|systems/<id>` no longer searches the
  planet and moon tables; a planet or moon may share another body's name.
  Generation already skipped them (since v34). Profiling a galaxy run
  showed name reservation at about 1% of generation time; row-by-row
  inserts are the main cost.

## [7.10.0] - 2026-09-30

### Changed
- **Every distance is shown in its most meaningful unit.** One helper,
  `stellarObjects.utils.format_distance_m` (with `_km`, `_au`, `_ly` and
  `_pc` wrappers, re-exported by `html/lib/fmt.py`), and its browser copy
  `html/static/distance.js` pick the largest of km < AU < mpc < cpc < ly
  < pc < kpc < Mpc < Gpc the value is at least 1 of. Parsec values add
  ly in parentheses ("4.2 pc (13.7 ly)"), or AU below 0.01 ly ("2.4 mpc
  (495 AU)"). The system, sector, galaxy and phenomenon pages, the
  System, Sector and phenomenon map readouts, the wiki sector page and
  the generated text (planet, star, binary, belt and heliosphere
  distances) all use it.
- **Planet, moon and star radii are always km in scientific notation**
  (`utils.format_body_radius_km`), including a rogue planet's.
- **The distance constants are exact:** AU = 149,597,870,700 m,
  lightyear = 9,460,730,472,580,800 m, parsec = 3.085677581491367e16 m,
  with every conversion derived from them. Stored values shift by about
  one part in 70,000.

### Removed
- The unused `LY_THRESHOLD`, `HELIOSPHERE_DISPLAY_THRESHOLD_LY` and
  `ROUND_HABITABLE_ZONE_AU(_SMALL)` display constants.

## [7.9.2] - 2026-09-30

### Fixed
- **The admin Generate page's jobs work on native Windows.** Cancel now
  writes a `cancel` file that the job runner checks while a step runs,
  and stops the step's whole process tree (`os.killpg` on POSIX,
  `taskkill /T /F` on Windows), so a cancelled generation step no longer
  keeps writing to the database. Liveness uses `OpenProcess` and
  `GetExitCodeProcess` on Windows (with a creation-time check against
  reused pids), the runner starts detached in its own process group
  (breaking away from IIS's job object where allowed), `state.json`
  writes retry while the page has the file open, and the private
  fallback directory check no longer needs `os.geteuid`. CI runs the job
  tests on a Windows runner.

## [7.9.1] - 2026-09-30

### Fixed
- **A binary's secondary star could outweigh its primary.** When the
  secondary's own class pushes its mass above the primary's, the two
  swap roles, so the primary is always the heavier star (TODO item 44).

## [7.9.0] - 2026-09-30

### Added

- Galaxy Map: wedge lines run out from the galactic core in the plane, each
  labelled with its bearing (degrees counterclockwise from +X, ring slot 0),
  with a Wedges button to hide them.
- Galaxy Map: the scale readout has three lines, what one screen pixel
  spans, how big one block is, and a bar, each in sectors, pc and ly.

### Changed

- Galaxy Map: the density prisms are shaded by each prism's arm factor
  (its density over the ring's mean) as well as its density, so the spiral
  arms stand out at every zoom.
- Galaxy Map: the density blocks fill their whole cells, so the galaxy is one
  solid made of blocks with no gaps. A Slice button (on by default) cuts the
  solid at the focus's layer, so the view looks down on its cut face;
  turning it off shows the whole solid.

### Removed

- The server's leftover density point clouds: `/api/galaxy/tiles` and
  `/galaxy/tiles` no longer take `density=` or return `density`, and
  `galaxyViewport.density_sample_points` / `density_points_for_tile` are
  gone. The map has drawn density itself since the prisms.

## [7.8.0] - 2026-09-30

### Changed

- Warp travel times follow Boss's warp curve, which matches warp^(10/3) up to about warp 9 and then climbs toward warp 10 (warp 9.995 is about 12,200 c). The NAV warp table covers warp 1, 2, 4, 8, 9, 9.5, 9.9 and 9.995.

### Added

- Dimensional fold travel times (6F^4 / (10 - F) times c) at fold 4 to 8.5 on the NAV page and in `/api/nav`'s `fold_times`.

## [7.7.0] - 2026-09-30

### Security
- **pip installs only locked, hash-checked files.** `requirements.lock`
  pins every runtime library and its dependencies to one version (per
  Python range) with the sha256 hashes of its files, and
  `scripts/install-python-deps.sh` installs from it with
  `--require-hashes` on both the ordinary and the externally managed
  path. apt-provided libraries are still used as they are.
  `scripts/lock-requirements.sh` (uv) regenerates it; a test fails when
  it no longer meets `setup.py`'s floors, and CI audits it with
  pip-audit.

## [7.6.1] - 2026-09-30

### Fixed
- **Sector growth could place star systems inside a black hole's or
  neutron star's Hill sphere.** `SpaceSector._fine_tune_position` now
  checks every massive neighbor (systems and placed compact remnants),
  not just other systems (TODO item 45).

## [7.6.0] - 2026-09-30

### Added
- **Deployment guides for every major platform.** `docs/deployment/` has
  a comparison table and a guide for each setup: Apache2 + mod_wsgi (the
  existing guide, moved from `docs/apache-deployment.md`), nginx +
  gunicorn and Caddy + gunicorn on Linux, three native Windows setups
  (IIS + HttpPlatformHandler, Caddy, Apache Lounge, each with waitress)
  plus WSL2, Homebrew nginx + gunicorn under launchd on macOS (macOS
  Server is discontinued), and a VPS and platform-as-a-service note.
  Example configs are in `examples/nginx/`, `examples/caddy/`,
  `examples/systemd/`, `examples/windows/` and `examples/macos/`.
- **`proxy_fix` setting for running behind a reverse proxy.** Behind
  nginx, Caddy, IIS or Apache's `mod_proxy`, the app saw the proxy's
  address for every visitor (so everyone shared one rate-limit budget)
  and never sent `Strict-Transport-Security`. `config.json`'s
  `proxy_fix` (`x_for`, `x_proto`, `x_host`; environment variables
  `PLANETGEN_PROXY_FIX_X_FOR`, `_X_PROTO`, `_X_HOST`) sets how many
  proxies to trust for each `X-Forwarded-*` header, applied with
  werkzeug's `ProxyFix`. All 0 by default, so Apache + mod_wsgi is
  unchanged. A value that isn't a whole number 0 or above stops the app
  at startup.

### Changed
- `docs/server-checklist.md` works for every platform, and its schema
  and index checks match the current schema.
- `docs/api.md` no longer suggests a read-only database account for the
  web app: the app writes through its one account (admin pages, write
  endpoints), as `docs/config.md` already said.
- Docs and code comments that named the old `sectorGen.py`,
  `systemGen.py`, `phenomenonGen.py`, `galaxyGen.py` and `galaxyPlan.py`
  scripts now name the `generate.py` subcommands.
- `docs/TODO.md` item 54 tracks the admin Generate page's POSIX-only job
  handling, which doesn't work reliably on native Windows.

## [7.5.0] - 2026-09-30

### Security
- **No more published first login.** A new install's admin gets a random
  password, printed once by `migrateDb.py` (so by `install.sh` and
  `update.sh`); it must still be changed at first login. An existing
  install whose admin is still on the old `admin`/`password` login gets a
  random password the same way on the next update. See `docs/api.md`
  ("The first admin login", "Resetting the admin login").
- **Changing the admin username or password logs out every other
  browser** signed in as that admin. API keys keep working.
- **Login takes the same time whether or not the username exists.**
- **Rate limits on the web pages**, per IP: search 30 a minute, the Galaxy
  Map 60, its tiles 600, `/api/health` 60 and every other page 300,
  configurable as `ratelimit.pages` in `config.json` (an empty value
  turns one off).
- **CSRF tokens are tied to the login session**, and a new API key's
  one-time display cookie is only sent to `/admin`.
- **Installer permissions:** the code under `src/html` is now owned by
  root and read-only for Apache; only the runtime directories belong to
  Apache's user. `config.json` is set to `root:<apache group>`, mode 640.
  The debug log is mode 0660 (CLI users must be in Apache's group to
  append). The root-run cache-directory helper no longer imports code the
  web user could change.
- **The `/tmp` fallbacks** for the jobs and tile-cache directories are
  created private and refused when another user owns them.
- **A sector's wiki link must be an `http` or `https` URL**, and one
  stored earlier with another scheme is no longer shown as a link.
- **`/api/databases` no longer lists the control database**, and `?db=`
  can't select it.
- **Database errors no longer reach visitors**: `/api/health` and the
  read routes answer "database unavailable" and log the detail.
- **HSTS** is sent on HTTPS responses, and `%`/`_` in the `star_type`
  filter are matched literally.

## [7.4.0] - 2026-09-30

### Added
- **Brute-force (property-based) tests across the whole codebase.** New
  `test_fuzz_*.py` files use Hypothesis to hunt for inputs that break the
  galaxy grid, system and planet generation, sector placement, names, text
  rendering, the physics helpers, config loading, every `generate.py`
  command, and every web page and API endpoint. New `test_edge_admin_scripts.py`
  covers the job runner, `resetDb`, `migrateDb` and `updateOrbits` against
  empty, broken and unreachable databases. Every normal `pytest` run uses
  a fixed, repeatable `ci` profile; a weekly **Deep fuzz** workflow runs
  3000 random examples per test. See `docs/testing.md`.

### Fixed
- **Bugs the new tests found.** Among them: a point a sector cell sampled
  could count as outside that cell; negative ring or out-of-range slot
  designations were accepted; the galaxy density could exceed 1; zero,
  NaN or infinite map radii picked the wrong tile level; a planet pinned
  by mass alone crashed; star types with trailing junk (`G2Vjunk`) were
  accepted; NaN or infinite values in the physics helpers came back as
  NaN instead of an error; NaN sector edges and positions poisoned
  placement; a sector's grid cell was lost on save and load; `generate.py`
  accepted NaN or infinite numbers, bad star types, broken
  `--system-file`s and out-of-range ports, and printed tracebacks for
  database and file errors; `generate.py plan` accepted shapes that could
  not be built; a huge `--radius-pc` tried to enumerate far more sectors
  than the galaxy has; decoration words and offensive words hidden by an
  apostrophe or space slipped into generated names; and the debug log did
  not redact every password or token option.
- **Web and API inputs that crashed or misbehaved.** A non-ASCII CSRF
  token, a query string that is not valid UTF-8, a page number past
  2^63, NaN or infinite search bounds and radii, a Unicode digit in a
  sector address, a huge galaxy cell coordinate, a non-JSON request body,
  non-string logins or key labels, and over-long names, URLs, usernames
  and key labels now get a clear 400 instead of a 500 or a traceback.
  `/api/health?db=` with an unknown database answers 404, and
  `/api/databases` lists a schema it can't open instead of failing.
- **A flaky job runner test.** It checked for the job lock in the instant
  between the runner writing its final state and releasing the lock.

## [7.3.0] - 2026-09-30

### Changed
- **Planets and moons are named after their star.** Only the system draws
  a generated name now. Planets are numbered in orbit order with roman
  numerals (`Voranthis I`, `Voranthis II`), moons add a letter
  (`Voranthis IIa`, `Voranthis IIb`), and asteroid belts take no number.
  A binary's two stars each put their own word after the system name
  (`Voranthis Kelmoor`, `Voranthis Ostra`) instead of `Voranthis` and
  `Voranthis B`; a close pair's planets use the system name and a wide
  pair's use their own star's (`Voranthis Kelmoor I`). See
  `src/planetgen/names/bodies.py`.
- **Planet and moon names are no longer searched for duplicates.** They
  derive from the system name, which is already unique, so saving a
  system skips a registry lookup and write per planet and moon, and the
  "Kin"/"Ami" companion suffixes are gone. When a system is renamed
  (`PATCH /api/systems/<id>`, or the Alpha/Beta and "Little" decorations
  a name collision applies), every star, planet and moon still carrying
  its name follows. Schema v34 drops `body_name_registry`; run
  `update.sh`. Existing rows keep their old names until the galaxy is
  regenerated.
- **`PATCH /api/systems/<id>` refuses a name already in use** by any
  sector, system, star, planet or moon, with a `409`.

### Added
- **Rename stars, planets and moons through the API.**
  `PATCH /api/stars/<id>`, `/api/planets/<id>` and `/api/moons/<id>` take
  `{"name": str}` and need an admin login, like the other writes.
  Renaming a single star renames its system; renaming a binary's star
  carries the planets and moons named after it. See `docs/api.md`,
  "Renaming".

## [7.2.1] - 2026-09-30

### Fixed
- **An asteroid belt could overlap the next planet out.** When spacing
  pushed a planet into the habitable zone and it was reclassified there,
  some classes redrew its distance anywhere in the zone, which could put
  it back inside the belt it had just been moved past. It happened most
  often around giant stars, whose habitable zones are wide, and made
  `test_systems.py` fail at random. A reclassified planet now keeps the
  distance it was moved to.

## [7.2.0] - 2026-09-30

### Changed
- **`update.sh` no longer reinstalls anything.** It used to re-run
  `install.sh` after every pull that brought new commits, which
  force-reinstalled the Python package with pip, rebuilt the fallback
  venv and re-fetched the NLTK corpus. Now it checks instead: each Python
  library is imported with the system Python and compared with
  `setup.py`'s floor (`scripts/install-python-deps.sh --check`), and only
  one that is missing, too old or broken is installed, the same way
  `install.sh` would on that host. The NLTK corpus and mod_wsgi are
  installed, and Apache's modules enabled, only when missing. It prints
  one line per library (present, installed, upgraded, repaired or failed)
  with where it came from, and finishes by importing the web app as
  Apache's user, so a library www-data can't use fails the update instead
  of the site.
- **No more venv: libraries go into the system Python.** On an
  externally managed Python (Ubuntu 24.04+), everything apt packages at
  or above `setup.py`'s floor comes from apt. Only what apt lacks or ships
  too old is pip-installed system-wide into `/usr/local`, alongside apt's
  copy and never over it. The report says which ones and why. An existing
  `/opt/planetgen/venv` and its `planetgen-venv.pth` are removed on the
  next run, and their libraries are installed system-wide.
- The NLTK and Apache module steps now live in `scripts/deploy-common.sh`,
  shared by `install.sh` and `update.sh`, and `install.sh` installs
  mod_wsgi when apt can instead of only warning about it.
- `/usr/local/bin/planetgen` is now the checkout wrapper on every host
  (pip's console script ran the pip-installed copy, which would go stale
  once updates stopped reinstalling it).
- **Checking which Python Apache really uses.** `install.sh`/`update.sh`
  now warn when mod_wsgi is built for a different Python version than the
  one they set the libraries up for, and both take `PYTHON=` to pick
  another interpreter. The admin Stats page shows the web app's Python
  prefix and the directory it imports its libraries from. See
  `docs/apache-deployment.md`.

### Added
- **Migrate or delete the database on update.** When the database is
  behind the current schema, `update.sh` (and `install.sh`) asks whether
  to delete the galaxy data instead of migrating it: y/N, with a
  30-second timeout that defaults to N (keep the data and migrate it).
  The same happens with no terminal to ask on. Deleting wipes every
  generated sector and system (`resetDb.py`) and keeps admin logins.
- **Progress bar and ETA for database migrations.** `migrateDb.py` shows
  each migration step as it runs, with the elapsed time and an estimate
  of the time left. `migrateDb.py --status` reports the current and
  target schema versions without changing anything.

## [7.1.0] - 2026-09-30

### Removed
- **The CGI pages are gone.** Every page is served by the Flask app, so the
  `src/html/*.py` redirect shims, the CGI page shell (`lib/page.py`),
  `fmt.post_link`/`data_nav_params`, `static/navform.js`, the POST mode of
  the pager and the helpers' `LEGACY_PAGES` are removed, along with the
  CSS for the old side rail and link buttons.
- **The example Apache vhost has no CGI rules.** `ScriptAliasMatch`,
  `ExecCGI`, `AddHandler cgi-script .py`, the `<Files "wsgi.py">` override
  and the CGI-only `SetEnv` lines are gone; `Alias /static/` and
  `WSGIScriptAlias /` remain. `install.sh` enables `wsgi` instead of
  `cgid`, and prints what to remove when an existing site file still has
  the CGI rules (see `docs/apache-deployment.md`, "Updating an existing
  server").

### Changed
- **Old `/<name>.py` URLs redirect in the app.** `web/old_urls.py` answers
  `/index.py`, `/sector.py?id=5`, `/search.py?...` and the rest with a 301
  to the page that replaced them, keeping their parameters, so old
  bookmarks still work once the vhost's `ScriptAliasMatch` is removed. Any
  other `.py` name is a 404.

## [7.0.1] - 2026-09-30

### Fixed
- **The 3D Galaxy Map drew a generated neighborhood as a lopsided
  half-sphere.** Each map tile returns at most 250 placed sectors, and past
  that it kept the lowest ids. Neighborhoods are generated in order
  outward from the core, so the lowest ids are the core-facing side: a
  100 ly sphere of ~2,770 sectors showed only its inner ~1,000, as a bowl
  with a flat face where the cap ran out. A full tile now returns every Nth
  sector by id, so the sample covers the whole neighborhood.

## [7.0.0] - 2026-09-30

### Changed
- **One sector standard: 4 parsecs.** Every sector, in the galaxy or
  standalone, is now 4 pc (about 13.05 ly) on a side instead of 11.5 ly. At
  the default Milky Way shape the galaxy then reaches about 50,000 ly, the
  real disk's radius. See `docs/design/galaxy-coordinate-system.md`,
  "Sector size".
- **The galaxy is a stack of aligned layers.** Ring `i` now holds
  `round(2π(i + ½))` slots (3, 9, 16, 22, …) instead of a multiple of 4,
  the same on every layer, so sectors line up in vertical columns. The
  skeleton is stored per layer: each layer, from the top of the galaxy to
  the bottom, runs from ring 0 out to the last ring that still expects a
  star per sector (`galaxy_layer`, replacing `galaxy_ring_band`).
- **Generation never lands outside the galaxy.** `generate.py plan` also
  stores each ring's column bound (`galaxy_column`: the highest and lowest
  layer it reaches), and every `generate.py galaxy` mode checks its address
  against the stored outline before generating anything. Explicit
  `--density` or `--num-systems` no longer skip that check; an address
  outside is refused with the reason, and a neighborhood near the edge
  leaves out the sectors past it. `generate.py galaxy` now refuses to run
  before `generate.py plan`.
- **Random starts are drawn from the real outline.** A random start picks a
  uniformly random sector inside the planned galaxy instead of from a fixed
  15,000 pc disk 2,000 pc tall, so the starting neighborhood can no longer
  land outside the galaxy. `--max-ring` defaults to the galaxy's own edge.
- The Galaxy Map's prisms use the generator's own one-star-per-sector
  threshold, so fully zoomed in their outline is exactly the galaxy's
  layers.
- `generate.py plan` no longer takes `--edge-ly` or `--empty-streak-to-stop`;
  the edge is always the standard, and the build takes a few milliseconds.

### Removed
- **Upgrading deletes every galaxy-placed sector again.** Schema v33's
  migration deletes each placed sector with its systems and phenomena, since
  nearly every address moves, and rebuilds the skeleton from the stored
  shape at 4 pc. Sectors that were never placed in the galaxy are kept.
  Regenerate the galaxy afterwards (`generate.py galaxy`), and take a
  backup first if you want the old data.

## [6.6.0] - 2026-09-30

### Changed
- **System and phenomenon pages moved to the new site.** A star system is
  now at `/system/<id>`, the phenomena list at `/phenomena` (paged with
  `?page=N`), and one phenomenon at `/phenomenon/<type>/<id>`. These are
  plain, bookmarkable URLs with no database name in them, and every link
  on them (neighbouring systems, sector, Navigate from/to here, the
  Wikitext/Markdown views with `?code=...`) is an ordinary link, so Back,
  reload and open-in-new-tab work. The pages use the new header and
  breadcrumbs, and the system map, body list, code views with the Copy
  button and phenomenon diagram work as before. Old `system.py`,
  `phenomena.py` and `phenomenon.py` links redirect to the new addresses.
- The admin "Upload to Wiki" form on a system page is now CSRF-protected
  and redirects back to the page afterwards with a fixed status message,
  so reloading never uploads twice. It also says so when the default
  admin credentials must be changed first.

## [6.5.0] - 2026-09-30

### Changed
- **The admin pages moved to the Flask app:** `/login`, `/logout`,
  `/account` (change username/password), `/admin` (API keys, a sector's
  wiki link; `?keys_page=N`) and `/admin/stats` (server and database
  stats; `?names_page=N`), replacing `login.py`, `logout.py`,
  `changecreds.py`, `admin.py` and `adminstats.py`, which now answer a
  `301` to their new URL. Sessions and the login rate limit work as
  before: the API's own session cookie is relayed unchanged.
- Visiting an admin page while logged out goes to `/login?next=<page>`
  and back there after login (only ever to a page on this site).
- `/admin/stats` and the wiki-link form use the site's configured
  database; the database picker and field are gone.

### Security
- Every admin form (log in, log out, change credentials, create or
  revoke a key, set a wiki link) is a POST with a CSRF token, answered
  with a redirect. A new API key reaches the page after the redirect in
  a signed, HttpOnly, SameSite=Strict one-time cookie, never the URL.
- Logging out is a POST: `GET /logout` only shows a "Log out" button,
  so a link or prefetch can't end a session.
- Admin pages are sent with `Cache-Control: no-store`.

## [6.4.0] - 2026-09-30

### Changed
- **The Galaxy Map moved to the Flask app at `/galaxy`.** It is a plain,
  bookmarkable GET URL with no database name in it
  (`/galaxy?quadrant=II&page=2` for one Quadrant's sector list), under
  the new header with the Galaxy section marked and breadcrumbs. The
  map's "View sector" button and every table row are real `<a href>`
  links. The map's script fetches its tiles from `/galaxy/tiles` (JSON,
  no database in the URL). Tiles are still cached on the server's disk,
  now by the WSGI process, and in the browser's `localStorage` under the
  same keys as before, so nothing already cached is lost. The map's
  viewport is larger on wide screens.
- `galaxy.py` and `galaxy_tiles.py` now answer 301 to `/galaxy` and
  `/galaxy/tiles`, keeping their parameters, so old links, bookmarks and
  open tabs still work. A tile cache directory set only with `SetEnv
  PLANETGEN_TILE_CACHE_DIR` in the vhost no longer applies (the WSGI
  daemon doesn't see `SetEnv`); set `tile_cache.dir` in `config.json`
  instead.

## [6.3.0] - 2026-09-30

### Added
- **One-off star systems from the admin site.** A new page,
  `/admin/generate/system` (linked from Generate), offers every
  `generate.py system` option, including a pasted system file and the
  debug narration, and shows the result as Markdown or wikitext with Copy,
  Download and a rendered preview. Nothing is saved to the database.
- **`generate.py system --output FILE`** (`-o`, `-` for stdout) writes
  the system's page instead of saving the system to the database.

## [6.2.0] - 2026-09-30

### Changed
- **The sector page moved to `/sector/<id>`.** It is a Flask page now, with
  the new header and breadcrumbs, bookmarkable Contents pages
  (`?contents_page=N`) and no database in the URL. The Sector Map's info
  panel buttons and its no-JavaScript list are plain links. The admin
  forms (wiki upload, generate neighborhood) carry a CSRF token and, once
  they succeed, redirect back to the page with a message, so reloading
  never repeats them. `sector.py` answers 301 to the new address.
- **The NAV page moved to `/nav`, and every step is a GET URL.** Endpoints
  read as `<kind>:<id>`: `/nav?from=system:12&to=nebula:3`; the pickers
  use `from_sector`/`to_sector`. The NAV Map's points and the route's
  stops are plain links, and a "Reverse course" link swaps the endpoints.
  `nav.py` answers 301 to the new address, translating its old
  parameters, and the old `from_id`/`from_kind`/`from_type` style
  redirects to the new one.

## [6.1.0] - 2026-09-30

### Added
- **Generate, plan and reset the galaxy from the web interface.** A new
  admin-only page, `/admin/generate` (the Generate link in the header),
  runs `generate.py plan`, `generate.py galaxy` (every mode: random start,
  whole ring, around a sector, one address) and `resetDb.py` as
  background jobs, plus a one-click "New galaxy" that resets, plans and
  generates a first neighborhood. Reset and New galaxy ask for the
  database name to be typed back. The running job shows its step, a
  progress bar, elapsed time and live output, and can be cancelled; the
  last 20 jobs keep their full output. New `jobs` section in
  `config.json` (`docs/config.md`).
- `generate.py` writes its progress to `$PLANETGEN_PROGRESS_FILE` when
  that is set (`stellarObjects/progressFile.py`).

## [6.0.0] - 2026-09-30

### Added
- **The Galaxy Map draws the sector grid itself.** Its prisms are the real
  cylindrical sector cells: one prism is one sector up close, and further
  out one prism stands for a block of whole sectors (3, 9, 27, ... a side),
  chosen so a block stays at least about 10 pixels across on screen and the
  view stays fast. Clicking a prism, or any empty spot, shows that sector or
  block: its address or ring and layer range, how many sectors it holds, its
  center in Cartesian, cylindrical and spherical coordinates, its size and
  its 8 corners, plus the command to generate a single sector.
- `GET /api/galaxy/cell?ring=&layer=&slot=` (or `?x=&y=&z=`) describes any
  sector cell in the galaxy, generated or not, the same way.

### Changed
- **Galaxy sectors now sit on a cylindrical grid instead of spherical
  shells.** Each galaxy-placed sector is one cell of rings 11.5 ly wide
  around the galactic axis, layers 11.5 ly tall (layer 0 centered on the
  galactic plane) and wedge-shaped slots about 11.5 ly across, so sectors
  follow the flat disk instead of a ball. A sector's address is
  `(ring, layer, slot)`; systems and phenomena are placed inside the real
  cell, with local axes pointing outward, along the ring and north. See
  `docs/design/galaxy-coordinate-system.md`, "Cylindrical sector grid".
- `generate.py galaxy` takes `--ring I [--layer J] [--slot K]` in place of
  `--shell K [--slot N]`, and `--max-ring` in place of `--max-shell`.
  `generate.py plan` builds one band of layers per ring and takes
  `--max-ring`; `--workers` and `--chunk-size` are gone, since the build
  now takes about half a second. The Galaxy Map's copied commands use the
  new flags.
- The Galaxy page's 100 ly radial groups are now called Zones.
- Designations encode ring, layer and slot, so every sector gets a new one.

### Removed
- **Upgrading deletes every galaxy-placed sector.** Schema v32's migration
  deletes each placed sector together with its systems and phenomena, since
  shell addresses have no matching cell, and rebuilds the skeleton from the
  stored shape. Sectors that were never placed in the galaxy are kept.
  Regenerate the galaxy afterwards (for example `generate.py plan`, then
  `generate.py galaxy`), and take a backup first if you want the old data.
- The `sector_vertices` and `galaxy_shell_band` tables.

## [5.59.0] - 2026-09-30

### Added
- **Browser checks for every Flask page.** A new test
  (`src/tests/test_web_a11y.py`, and its own `browser-a11y` CI job) loads
  each page in headless Chromium at phone and desktop widths, in light and
  dark, and fails on serious or critical axe-core (WCAG 2.1 AA)
  violations, horizontal page scroll, console or CSP errors, or a missing
  skip link or `aria-current` marker. The page list comes from the app's
  routes, so pages moved off CGI later are checked automatically.
  axe-core 4.13.0 is vendored under `src/tests/vendor/axe-core/`; the new
  `browser` extra installs Playwright.

### Fixed
- **The current section in the header was below WCAG AA contrast** in the
  light theme (4.45:1); it now uses the link colour.

## [5.58.0] - 2026-09-30

### Changed
- **`install.sh` now works on an externally managed Python (PEP 668),
  such as Ubuntu 26.04 LTS's.** Where pip used to fail with
  "externally-managed-environment", the installer now detects the
  `EXTERNALLY-MANAGED` marker and installs the libraries as apt packages
  instead, pip-installing only what the distribution lacks (or packages too
  old) into a venv at `/opt/planetgen/venv`. On an unmanaged Python it still
  uses pip exactly as before. It prints which path it took; see
  `docs/apache-deployment.md`'s "Managed Python".

## [5.57.0] - 2026-09-30

### Changed
- **The search page moved to `/search` (Flask), with bookmarkable GET
  URLs.** Every filter is a query parameter: `q` (the header search box)
  searches sector, system, star, planet and moon names at once; the
  per-object name fields (`sector_q`, ...), size ranges
  (`planet_min_radius_km`, ...), repeated tag facets
  (`spectral=G&spectral=K`) and each result panel's page
  (`stars_page=2`) follow it. Tags and "remove filter" chips are plain
  links, the per-object fields fold into a "Search by object and size"
  section, and results now appear above the tag browser. A submitted
  form's empty fields are dropped by a redirect to the short URL.
  `search.py` is now a shim that 301-redirects to `/search`, keeping every
  search parameter from an old link, bookmark or form post.

## [5.56.0] - 2026-09-30

### Added
- **Quasars.** A galaxy's nucleus can now be active: 10% of the time
  (`QUASAR_ACTIVE_NUCLEUS_CHANCE`) the first core sector gets a quasar at
  the exact galactic center, so a galaxy has at most one. Each has a
  supermassive black hole (1e8-1e10 solar masses), a luminosity set by
  its Eddington ratio, the matching accretion rate and broad-line-region
  size, and ~10% are radio-loud with jets. It shows on the Sector Map,
  in the sector's Contents table, in the phenomena list and on its own
  detail page, and `generate.py phenomenon --type quasar` makes one on
  demand (`--sector-id` must be a shell-0 sector). New `quasars` table,
  schema v31; run `migrateDb.py` (update.sh does).

### Fixed
- **A black hole's accretion-disk temperature could render one kelvin
  low.** Loading an anchored black hole back from the database truncated
  its fractional disk temperature, so the rendered page could read e.g.
  3,676,064 K instead of 3,676,065 K.

## [5.55.0] - 2026-09-27

### Changed
- **The 3D Galaxy Map shades predicted density with cylindrical segment
  prisms instead of spheres.** Space is cut into rings, wedges and layers
  on the same grid as cylindrical sectors, a power-of-two number of
  sector widths across so the prisms scale with the view (one prism is
  one sector at full zoom). Each prism is solid and lit, colored and
  sized inside its cell by its mean density. The density is computed in
  the browser from the galaxy's own shape (`static/galaxyprisms.js`), so
  the map no longer asks the server for a density point cloud.

## [5.54.0] - 2026-09-24

### Added
- **The home page is now served by the Flask app, with a new header.**
  `/` shows every sector and standalone system (each table paged on its
  own with `?sectors_page=N`/`?standalone_page=N`), and `/sectors` and
  `/systems` show one table each. These are plain, bookmarkable GET URLs
  with no database name in them (the database comes from `config.json`'s
  `mysql.database`). The new header has the sections (Galaxy, Sectors,
  Systems, Phenomena, Nav) with the current one marked, a search box,
  Login or Admin/Stats/Logout and the theme button. On phones these fold
  into a Menu that works without JavaScript. There is also a "Skip to
  content" link and breadcrumbs. The pages are Jinja2 templates in
  `src/html/web/` (autoescaped), served without a process start or an
  HTTP call back to the API. The other pages still run as CGI and move
  over in later releases; see `docs/html-interface.md`, "Flask pages".
- `secret_key` in `config.json` (or `PLANETGEN_SECRET_KEY`), used to sign
  CSRF tokens for the Flask pages' forms.

### Changed
- `index.py` and `browse.py` now answer 301 to `/`, keeping their page
  numbers, so old links and bookmarks still work.
- The example Apache vhost mounts the Flask app at `/` instead of `/api`,
  serves `/static/` with `Alias` and runs the remaining CGI pages through
  `ScriptAliasMatch`. **Existing servers need their vhost updated**: see
  `docs/apache-deployment.md`, "Updating an existing server".
- The API's app-wide default rate limit no longer counts the calls pages
  make in-process. Login and write limits still apply.
- Unknown URLs outside `/api` get an HTML 404 page instead of JSON.

## [5.53.2] - 2026-09-24

### Fixed
- **Airless planets and moons could keep an atmosphere they no longer
  had.** A body moved into an airless class after changing zones kept its
  old class's atmosphere density, molar density and scale height. Those
  values are now cleared when it's reclassified.
- **Bodies around very dim stars could be colder than space itself.**
  Surface temperatures are now floored at the cosmic microwave background
  (2.725 K).
- Schema v30 cleans up both problems in rows already saved: stale
  atmosphere values on airless bodies go back to empty, and temperatures
  below 2.725 K are raised to it. Run `migrateDb.py` (update.sh does).

## [5.53.1] - 2026-09-24

### Fixed
- Removed the old `/galaxy/view` Galaxy Map endpoint (`html/galaxy_view.py`, `GET /api/galaxy/view`). The map has used cube tiles since 5.47.0, but a browser holding a cached copy of the old map script kept calling it, and its unbounded query timed out and tied up the API. Such a request now gets a quick 404 instead.
- The Galaxy Map page now loads its script as `static/galaxymap3d.js?v=<version>`, so each release reaches browsers on their next page view instead of a stale cached copy.

## [5.53.0] - 2026-09-24

### Added
- **A light/dark/system theme switch** at the bottom of the side nav. The choice is
  remembered in the browser (`static/theme.js`) and applied before the page first
  paints, so there is no flash of the wrong theme. `style.css` gains the matching
  `:root[data-theme="dark"]` token block.
- **A favicon** (`static/favicon.svg`), a `<meta name="description">`, and
  `<meta name="theme-color">` for light and dark.

### Fixed
- **Phones got the desktop layout zoomed out.** The page shell now has
  `<meta name="viewport" content="width=device-width, initial-scale=1">`, so the
  existing narrow-screen rules apply. Those rules are also fixed: the side nav
  becomes a wrapping row of links instead of full-width buttons, the page gets an
  even gutter, and wide tables (including ones in wiki notes) scroll inside their
  own box instead of widening the page.

### Changed
- **Static files are versioned and cacheable.** Every CSS/JS/icon URL carries
  `?v=<release version>` (`fmt.static_url`), and the map scripts pass the same
  version on to three.js and `bodyRendering.js`. The Apache example now tells
  browsers to keep versioned `static/` files for a year, revalidates unversioned
  ones, and compresses HTML, CSS, JS, JSON and SVG with `mod_deflate` (three.js
  goes from about 740 KB to about 190 KB). Needs `a2enmod headers deflate`.

### Security
- **One source for the pages' security headers, and a stronger CSP.**
  `lib/page.py`'s `SECURITY_HEADERS` is the only place the HTML pages get them;
  the Apache example no longer repeats them vhost-wide (remove those three
  `Header always set` lines from an existing site config). The CSP adds
  `base-uri 'self'; form-action 'self'; frame-ancestors 'none'; object-src 'none'`
  to `default-src 'self'`.

## [5.52.0] - 2026-09-24

### Changed
- **System pages are rendered natively from the database.** The System
  panel is now an expandable list of the system's stars, planets (moons
  nested under each), asteroid belts and comets, each row showing class,
  terrestrial/gas giant, Habitable yes/no and Inhabited yes/no, and
  opening onto that body's own description. Wikitext and Markdown buttons
  show the full generated page in a code box with a Copy button.
- **Wiki page text is no longer stored (schema v29).**
  `star_systems.wikitext_content`/`markdown_content` are dropped; both
  formats are rendered on demand from the system's rows
  (`stellarObjects/systemRender.py`), so pages now follow renames, names
  made unique after generation, and orbit ticks. Wiki upload uses the
  fresh render. `GET /api/systems/<id>` no longer returns the two text
  fields; use the new `GET /api/systems/<id>/text?format=wikitext|markdown`
  and `GET /api/systems/<id>/sections`. **The migration deletes the stored
  copies: back up first** (`mysqldump`, or `src/checkRenderParity.py
  --export-dir`, which also compares the stored and rendered text on a
  database still at v28).

### Fixed
- **Binary stars could load with their primary and secondary swapped**
  when the second star generated was the heavier one, changing the
  "Barycenter Offset", planetary-limit and "Location" lines on reload.

## [5.51.0] - 2026-09-24

### Added
- **Every kind of stellar phenomenon now appears on the Sector Map.**
  Supernova remnants, rogue planets and interstellar comets had no galaxy
  position, so they were missing from the Sector Map, the sector's
  listing and NAV. They now have one (schema v28), drawn as a glowing
  shell, a dim world and an icy coma, and each is clickable like the
  others. A supernova remnant's leftover black hole or neutron star sits
  at the remnant's center. `migrateDb.py` gives existing ones a random
  spot inside the sector they were generated in.

### Changed
- **The sector page lists its phenomena alongside its systems.** One
  "Contents" table replaces the separate Systems and Nearby Exotic
  Phenomena tables, nearest the sector's center first, with a distance
  column, 50 rows a page. A phenomenon generated as part of a sector is always listed
  there, even if an older placement put it outside the cube.

## [5.50.0] - 2026-09-24

### Added
- **A debug log for the whole program.** `"debug": true` in `config.json`
  (off when missing) makes the generator CLI, the maintenance scripts, the
  web pages and the API all write to `/var/log/planetgen.log` (`log_file`
  to move it): every decision the generator makes and why, every random
  roll with the source line that asked for it and the probabilities it was
  compared against, every SQL statement, web request, API call and admin
  access check, and every error with its traceback, all timestamped to the
  millisecond. A seeded run generates the same result with it on or off.
- **Log rotation for it.** `install.sh`/`update.sh` create the log file
  when debug is on (writable by Apache and shell users alike) and install
  `/etc/logrotate.d/planetgen`: daily, or as soon as it passes 100 MB
  (checked hourly), keeping 7 compressed copies.

### Changed
- **Web pages no longer show tracebacks when debug is on.** They went to
  the page itself before; now they go to the debug log, and the 500 page
  says so.

## [5.49.2] - 2026-09-24

### Changed
- **Editing a sector no longer throws away the whole Galaxy Map cache.**
  The tile cache used to go stale on any change to the database. It now
  asks the new `GET /api/galaxy/changes`, which reads the schema-v27
  `modified_at` columns plus new sector and system ids, which cube tiles
  changed, and deletes just those, on the server's disk and in each
  visitor's browser. A rename refetches the dozen tiles holding that
  sector. Deleting a sector, re-planning the galaxy or a new release
  still refreshes everything, since a deleted row leaves nothing to
  locate its tiles by.

## [5.49.1] - 2026-09-24

### Changed
- **Every list on the site now pages 50 rows at a time with the same
  pager.** Browse's sectors and standalone systems, Phenomena, a sector's
  systems and nearby phenomena, a Galaxy Map Quadrant's sector list, each
  Search result panel, the admin API key list and the admin stats page's
  duplicate-names list all share one control
  (`src/html/lib/pagination.py`): a "Showing X-Y of Z" summary, First/Prev,
  numbered pages, Next/Last. Phenomena no longer stops at 500 rows and
  Search no longer stops at 300 matches per panel; every match is
  reachable a page at a time.
- `GET /api/search` pages each result panel: `limit` (default 300, as
  before) plus `sectors_offset`/`systems_offset`/`stars_offset`/
  `planets_offset`/`moons_offset`/`belts_offset`, and each panel now
  reports `total`, `limit` and `offset` alongside `truncated`.

## [5.49.0] - 2026-09-24

### Added
- **Admin stats page.** Logged-in admins get a new Stats page (sidenav
  and a link on the Admin page) showing server health (API version,
  uptime, load, memory, MySQL version/uptime/connections, galaxy tile
  cache usage and free disk) and stats about the current database: exact
  sector and system counts, size on disk, schema version, when rows were
  last created or modified, and per-table row estimates and sizes.
- **Names made unique.** The same page counts every name the uniqueness
  rules had to decorate (Alpha/Beta..., Little..., ...Kin) and lists each
  one with links to every sector and system carrying it; planets and
  moons link to their system.
- New admin-only endpoints `GET /api/admin/stats` and
  `GET /api/admin/duplicate-names` (see `docs/api.md`).

## [5.48.3] - 2026-09-24

### Fixed
- **Zoomed in on the Galaxy Map, unfilled sectors still drew as a ball.**
  Two things drew it. Unfilled (not-yet-generated) sector dots were
  clipped to a 20 pc ball around the view center while the view reached
  out to 200 pc, and near the galactic plane every slot qualifies, so the
  ball was solid. They now fill the whole view, shown while the view
  radius is 32 pc or less, and fade out toward its edge. Wider zoomed-in
  views get the density cloud instead, which had its own ball: a view a
  few hundred parsecs across kept almost none of its galaxy-wide samples
  and topped up with uniform points around the view. It now samples such
  views locally against the real density, so the cloud follows the disk.

## [5.48.2] - 2026-09-24

### Changed
- **The browse page's sector list now runs outward from the galactic
  core.** Sectors are sorted by distance from the core (nearest first),
  with sectors never placed in a galaxy listed after them by name, and a
  new "Distance from core" column shows each one's distance. Paging walks
  the same order.
- **A sector's systems table now runs outward from the sector's center.**
  Systems are sorted by distance from the sector's center (nearest first),
  with a new "From center" column. `GET /api/sectors` gains
  `galactic_radius_pc`/`galactic_radius_ly` and `GET /api/sectors/<id>`'s
  systems gain `center_distance_ly`.

## [5.48.1] - 2026-09-24

### Fixed
- **Galaxy Map clicks recentered the view far from where you clicked.**
  Clicking empty space picked a point on a sphere the camera itself sits
  on, so the view jumped thousands of parsecs away (often out of the
  galaxy entirely). It now lands on the spot under the cursor in the
  galactic plane (or facing the camera when the disk is seen edge-on),
  kept inside the galaxy.
- **Double-click zoomed somewhere other than where you double-clicked.**
  Each of its two clicks recentered again before the zoom; now only the
  first click recenters and the double-click zooms in on that point.
- **Clicking a sector centered on the edge of its marker, not the sector,**
  so it slid off-center and out of view as you zoomed in. It now centers
  on the sector's own position.
- **Scroll zoom ignored how far the wheel moved.** A trackpad's stream of
  tiny scroll events each took a full step, making zoom race; zoom now
  follows the scroll amount (about 1.28x per mouse-wheel notch, pinch
  supported).
- **Zooming all the way out didn't show the whole galaxy.** The farthest
  zoom now backs the camera off far enough for the galaxy's whole disk to
  fit in the map.

## [5.48.0] - 2026-09-24

### Added
- **Created and modified times on the main tables.** `sectors`,
  `star_systems` and the seven exotic-phenomenon tables now have
  `created_at` and `modified_at` columns, with an index on `modified_at`
  (schema v27). MySQL keeps `modified_at` current on every edit. Changing
  a planet or moon bumps its system's `modified_at` instead of the child
  row getting its own timestamp, and deleting a system bumps its sector's.
  Orbit updates from `updateOrbits.py` don't count as a modification. The
  migration adds the columns without rebuilding the tables where MySQL
  supports it, and builds the indexes online, so it's safe to run on a
  large database. Existing systems and sectors are backfilled from the
  systems' original creation times, in small batches; existing
  phenomena, which have no record of when they were made, get the
  migration's own time.

## [5.47.3] - 2026-09-24

### Added
- **`docs/server-checklist.md`**, a step-by-step check to run on the server
  after a deploy: confirms the code version, that `migrateDb.py` has brought
  every database to schema v26 (an Apache restart alone doesn't), the
  spatial indexes, `request-timeout=60`, the Galaxy Map tile cache, and that
  the Galaxy Map no longer times out or runs the API out of memory.

## [5.47.2] - 2026-09-24

### Changed
- **`install.sh` and `update.sh` now create the Galaxy Map's tile cache
  directory** (`/var/cache/planetgen/tiles`, or `tile_cache.dir` from
  `config.json`) and give it to Apache's user, so the cache doesn't fall
  back to a private `/tmp` folder that's cleared on every restart. The
  step is also available on its own as
  `sudo examples/apache/create-cache-dir.sh [dir]`.

## [5.47.1] - 2026-09-24

### Fixed
- **Galaxy Map: a generated sector could show up as a huge sphere
  filling the view.** Each placed sector's dot is kept a constant size
  on screen, but its two sprites were being sized from the camera's
  distance to the galactic center instead of to the sector, so a sector
  far from the core grew into a giant disc once zoomed in on it. The
  size now comes from the sector's own position, and the dots are
  smaller (3-7 px instead of 4-16 px).
- **System Map: orbit lines were missing.** The orbits layer's
  `z-index: -1` put it behind the map's own background, because the map
  box didn't form its own stacking context. The viewport now isolates
  its stacking, so orbits draw again, still hidden behind each body's
  sphere.
- **Star and atmosphere glows rendered as flat, opaque discs** on the
  System Map and Sector Map. The glow shader measured the rim the
  front-face way on a back-face-only shell, which clamps to full
  brightness everywhere. It now fades from full brightness at the body's
  limb to nothing at the shell's edge (`bodyRendering.makeGlowMaterial`
  takes the shell's scale for this), and the star and planet glow
  strengths were retuned for the new falloff.
- **The "nearest" systems on a system page weren't links.** They were
  parsed out of the stored `location` text, whose names are written
  before each neighbor's name is made unique and never follow a rename,
  so they rarely matched a real system. `GET /api/systems/<id>` now
  returns `nearest_neighbors` (id, current name, distance) computed from
  positions, and the page links every one of them.

### Added
- **`tests/test_skeleton_shape.py`** confirms the unfilled-sector
  skeleton's slot centers form a thin disk in galaxy-frame parsecs; the
  sphere the Galaxy Map used to draw came from its old 200 pc
  planned-tier radius cap, which the cube tiles in 5.47.0 replaced.

## [5.47.0] - 2026-09-24

### Fixed
- **The Galaxy Map could take the whole site down (Apache OOM-killed on
  2026-09-24).** One view request listed every not-yet-generated sector
  slot within 200 pc of the camera target (~770,000 of them) before
  keeping the nearest 4,000: ~400 MB and up to 140 s per request, several
  at once across the API's threads. Slot search now skips slots outside
  the view's azimuth (same results, ~100x faster), keeps only the nearest
  results in memory, and stops at 40 pc.

### Added
- **The Galaxy Map loads space one cube ("tile") at a time, like a web
  map.** New `GET /api/galaxy/tiles` returns fixed cubes of an octree
  (placed sectors capped per cube, not-yet-generated slots only for
  16 pc cubes) and `GET /api/galaxy/stamp` a token that changes when the
  galaxy's contents do. No request can scan an unbounded region.
- **Tiles are cached on the server's disk and in the browser.** The web
  layer (`html/lib/tilecache.py`, behind the new `galaxy_tiles.py`) only
  asks the API for tiles it hasn't cached, and the map keeps tiles in
  `localStorage`, so revisiting or reloading doesn't refetch them. New
  sectors show up within a minute. Configure with `tile_cache.dir`/
  `tile_cache.max_mb` or `PLANETGEN_TILE_CACHE_DIR`/
  `PLANETGEN_TILE_CACHE_MAX_MB`.

### Changed
- The example Apache config runs the API daemon with
  `request-timeout=60`, so a runaway request restarts the daemon instead
  of piling up until the OOM killer stops Apache.
- The planned-slot ball no longer appears around the galactic core at
  full-galaxy zoom; planned slots show within 20 pc of the camera target
  once zoomed in.

## [5.46.35] - 2026-09-24

### Changed
- **`docs/TODO.md` is now grouped and numbered in working order.** Bug
  fixes have their own group, led by the production crash; everything
  else is grouped by what it changes (performance, System Map, Sector
  Map, Galaxy Map, Web API) and numbered from the smallest change to the
  largest, with a short "what to do first" plan at the top.

## [5.46.34] - 2026-09-24

### Changed
- **`docs/TODO.md` now tracks open work only.** The long completed-work
  log is gone (`CHANGELOG.md` and git history are the record), and nine
  new open items were added with their symptoms and starting files: the
  production Apache OOM crash, a Galaxy Map rework, the unfilled-sector
  skeleton drawing as a central sphere, removing the leftover star marker
  from the Galaxy Map, a request/database cache, placing and drawing all
  seven phenomenon types on the Sector Map, the opaque star glow, missing
  System Map orbit lines, and a drawn, obstacle-avoiding Measure distance
  route.

## [5.46.33] - 2026-09-24

### Changed
- **PRs no longer bump the version themselves, so parallel PRs stop
  colliding on the same number** (as 5.41.0/5.41.1 and 5.46.15 did, each
  needing a hand-renumbering merge). A PR now adds one note file under
  `changes/` (see `changes/README.md`); after it merges, the new
  `.github/workflows/stamp-version.yml` runs `scripts/bump_version.py`,
  which gives each pending note its own version, writes it into
  `_version.py`, the README badge and `CHANGELOG.md` together, and commits
  `Release x.y.z` to `main`. The new `Release note` PR check fails a PR
  that edits the version itself or adds no note (unless it's labelled
  `no-release`).

## [5.46.32] - 2026-09-24

### Added
- **System Map: "Measure distance" -- click any two stars, planets, or
  moons in the same view for the real distance between them.** A new
  toggle button puts the map into selection mode; the two clicked bodies
  highlight and the info panel shows the real straight-line distance
  (from each body's own true, un-log-scaled km position, not its drawn
  pixel position -- the shared log radial scale that places markers on
  screen preserves real angle but not real distance). When that straight
  line would pass through the scene's own center body (the star, or --
  one level in, a moon scene -- the planet drilled into), a second
  "around it" figure is also shown: the exact shortest path that clears
  the obstacle (two tangent lines plus the arc between them), not just a
  flagged "blocked". Works for a binary's own two stars too, including a
  close pair's small real offset from their shared barycenter.

## [5.46.31] - 2026-09-24

### Added
- **Sector Map: clickable indicators for every immediately surrounding
  sector.** A small marker sits just past the scene's own edge in the
  real direction of each same-shell (lateral) Voronoi neighbor
  (`sectorGeometry.lateral_neighbor_slots`, exact) and the nearest
  inward/outward radial neighbor (`sectorGeometry.radial_neighbor_slot`,
  nearest-by-distance). An already-generated neighbor's indicator links
  straight to it; a not-yet-generated one shows its address and a
  copyable `generate.py galaxy --shell K --slot N` command, the same
  convention the Galaxy Map's own "planned" tier already uses.
  `queryDb.sector_neighbors` (also folded into `sector_detail`'s own
  `neighbors` key) drives this from the API side.

## [5.46.30] - 2026-09-24

### Changed
- **Sector Map: stars, nebulae, asteroid fields, black holes, and neutron
  stars are now real, textured, glowing 3D spheres** instead of flat
  camera-facing sprites, matching the System Map's own body rendering.
  Each body is a textured core mesh (star granulation, nebula/asteroid/
  compact-remnant textures reused unchanged as sphere surfaces) plus a
  fresnel rim-glow shell sized and colored per body kind (bright corona
  for stars/neutron stars/accreting black holes, a softer shell for
  nebulae/asteroid fields). The glow shader and star granulation texture
  are now shared with the System Map via a new `static/bodyRendering.js`
  module rather than duplicated between the two files.

## [5.46.29] - 2026-09-24

### Changed
- **Galaxy Map (3D): a placed (already-generated) sector's own marker is
  now colored by its real stellar density** (`system_count / edge_ly **
  3`, relative to `physical_constants.LOCAL_STELLAR_DENSITY_LY3` -- the
  real local-neighborhood average this whole generator already
  calibrates against), not just sized by raw system count. Marker size
  still scales with `system_count` as before; only the color (dim bronze
  at low density, bright gold at high) is new. The info panel gained a
  "Density" field showing the same ratio (e.g. "1.8x local average").
  `queryDb.galaxy_sectors_in_view`'s own returned shape gained `edge_ly`
  per sector to make this possible.

## [5.46.28] - 2026-09-24

### Fixed
- **Galaxy Map (3D): rapidly clicking/double-clicking or mashing the
  zoom buttons could fire off a live API request per click, faster than
  the server can service them.** `doFetch`'s own `activeAbort.abort()`
  only stops the *browser* from waiting on a superseded response -- it
  doesn't reliably stop the server from finishing a query it already
  started (Flask/WSGI doesn't check for a disconnected client mid-query
  unless specifically coded to), so rapid clicking still burned a real
  WSGI thread/DB-connection-pool slot per click even when every earlier
  response got thrown away client-side the instant the next one fired --
  a real contributor to the production connection-exhaustion pattern
  already fixed elsewhere in this and recent releases. Every interaction
  that requests an immediate fetch (click, double-click, the +/-/reset
  buttons) now funnels through one shared cap: at most 4 accepted
  immediate fetches per second. A click faster than that still moves the
  camera/selection instantly (never throttled), it just falls back to
  the existing debounced delay for its own data fetch instead of firing
  right away, so a rapid burst still settles on exactly one fetch
  shortly after it stops rather than either hammering the server once
  per click or never syncing to the final camera position at all.

## [5.46.27] - 2026-09-24

### Added
- **`browse.py`'s Sectors and Standalone Systems tables are now really
  paginated** (100 rows/page, independent Prev/Next controls per table)
  instead of a single page capped at 500 rows with a "try Search
  instead" hint and no way to ever reach anything past that cap. Each
  table's own `sector_offset`/`standalone_offset` page position is
  independent, so paginating one never resets the other back to page 1.

## [5.46.26] - 2026-09-24

### Fixed
- **Rogue planets and interstellar comets were completely invisible
  everywhere** -- generated at a non-trivial rate
  (`program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM`'s own
  `"rogue-planet": 0.1` and `"comet": 0.05`, roughly one rogue planet per
  ten star systems, far more common than a nebula) and saved to the
  `rogue_planets`/`interstellar_comets` tables the whole time, but no
  query function anywhere (`list_phenomena`, `count_phenomena`,
  `phenomenon_detail`) ever read either table, so they never appeared in
  the Phenomena listing or had a detail page of their own, despite
  existing in the database. Both are now wired up the same way
  `supernova_remnant` already was (no galaxy-frame placement columns of
  their own, so still absent from the Galaxy Map / Sector Map / NAV --
  only the flat listing and detail page gain them). `html/phenomenon.py`
  gained field specs for each type's own real columns (a rogue planet's
  `planet_type`/`mass_kg`/`composition`/etc., a comet's
  `nucleus_diameter_km`/`velocity_kms`/`is_active`/etc.).

  Nebulae, by contrast, were already fully wired up -- their apparent
  rarity is real, deliberately-researched astronomical calibration
  (`"nebula"` rate is a cited-literature `5e4 / 2e11` per star system,
  several orders of magnitude below a rogue planet's own rate), not a
  bug; expect one to actually appear only in a very large generated
  galaxy.

## [5.46.25] - 2026-09-24

### Fixed
- **`migrate_database` never actually applied 5.46.19's v26 migration**
  (the spatial indexes on `nebulae`/`asteroid_fields`/`black_holes`/
  `neutron_stars` that fix `sector.py`'s production timeout) -- the
  `_migrate_v25_to_v26` function existed but was never called from
  `migrate_database`'s own version cascade, so running `migrateDb.py`
  against an existing (pre-v26) database left it silently stuck at v25
  forever, `GET /api/health`'s `schema_current` reporting `false`
  indefinitely with no way to clear it short of applying the index by
  hand. Caught by running the full test suite against a real MySQL
  server rather than relying on syntax/logic review alone -- every
  `test_migrate_v*` regression test failed with `schema_version == 25`,
  not `SCHEMA_VERSION` (26). A fresh database (`_ensure_schema`, which
  reads the indexes straight from `schema.sql`'s own `CREATE TABLE`) was
  never affected -- only a database migrated from an earlier version.

## [5.46.24] - 2026-09-24

### Fixed
- **System Map: an orbit line drawn under a body's own live sphere
  visibly cut across it** instead of being occluded -- the sphere is
  drawn on a separate `<canvas>` layered by CSS z-index relative to the
  SVG, so within a single `<svg>` a body's own opaque sphere had no way
  to occlude a sibling orbit-line element painted in that same stacking
  context. Each scene is now built as two sibling `<svg>`s (an
  aria-hidden orbits-only layer, and the existing body-marker layer,
  both toggled together by `static/systemmap.js`'s `showScene`) stacked
  either side of the sphere canvas, so a sphere now actually covers the
  orbit line drawn under it. Belt rings (interactive markers, not
  decorative lines) stay in the body-marker layer.

### Changed
- **System Map: stars now render with a mottled granulation texture**
  (layered sine "turbulence" in both UV directions, tinted to the star's
  own spectral color) instead of a flat single-color sphere.
- **System Map: a star's own glow shell is now bigger and brighter than
  a planet's subtle atmosphere rim** (wider falloff, higher intensity,
  larger radius), reading clearly as a light source rather than the same
  faint haze a planet's atmosphere gets.

## [5.46.23] - 2026-09-24

### Changed
- **Sector Map's compass arrow (pointing toward the galactic center) now
  labels itself plain "N"**, matching a real map's compass-rose
  convention, instead of the more verbose "Galactic Center →" text.

## [5.46.22] - 2026-09-24

### Fixed
- **Galaxy Map (3D): a selected placed/planned sector could disappear
  entirely once zoomed in close to it.** The live re-fetch's own bounding
  box shrinks as the camera's orbit radius shrinks; a click/double-click
  that landed even slightly off a sector's own exact stored position
  (easy from far out, where its marker is only a handful of screen
  pixels) meant a later, smaller-radius re-fetch could legitimately no
  longer include it, and the client dropped anything missing from a
  fresh fetch. The selected entry is now pinned client-side and
  re-inserted into each fetch's own tier if the live query didn't happen
  to return it, so it stays in the scene for as long as it's selected.

## [5.46.21] - 2026-09-23

### Changed
- **Galaxy Map (3D) interaction model reworked:** left-click now only
  centers the view on the clicked dot/empty space and selects it
  (previously it also zoomed in, which punished an imprecise click by
  zooming into empty space nowhere near the intended target -- the
  likely cause of "I can't zoom into known space" once a dot was too
  small/far to click precisely from the full-galaxy starting view).
  Double-click now does what a single click used to (center, select,
  AND zoom in by one `clickZoomFactor` step). Right-click no longer does
  anything (previously zoomed out) -- the browser's own default context
  menu is left alone instead of being suppressed for nothing. The panel's
  own hint text/`aria-label` (`lib/galaxymap3d.py`) updated to match.
- **The illustrative density cloud (the "shows the spiral arms" layer)
  now renders as soft, translucent, additively-blended spheres
  (`THREE.InstancedMesh`) instead of tiny flat `THREE.Points` dots** --
  it read as a sparse scatter-plot rather than shaded spiral structure.
  Each sphere's size varies with its own `relative_density` (denser
  regions read as visibly bigger/brighter blobs) and with the camera's
  current orbit radius (so the cloud keeps a sensible relative size
  across zoom levels); additive blending lets overlapping spheres
  brighten rather than simply occlude, the cheap way many soft blobs
  merge into continuous-looking shading along a spiral arm.

## [5.46.20] - 2026-09-23

### Fixed
- **A `:80`-to-HTTPS redirect vhost with no exclusion for `/api` silently
  turns every CGI page's own internal API call into a real public
  round trip back into the same server, doubling Apache's connection
  load per page view.** Confirmed against a real production vhost:
  `html/lib/apiclient.py` defaults to `PLANETGEN_API_BASE_URL=http://
  127.0.0.1/api`, a plain-HTTP loopback call every page makes at least
  twice (its own data, plus `/api/auth/me`); a blanket `Redirect
  permanent /` (or an equivalent `RewriteRule`, including certbot's own
  default `--apache` rewrite) on the `:80` vhost catches that loopback
  call too, and `urllib` follows the redirect out through DNS/TLS/
  anything in front of the box and back in -- exactly the connection-
  exhaustion/timeout pattern (with refused connections that never reach
  the error log, since Apache never accepted them) reported alongside
  the `sector_detail`/`galaxy/view` timeouts this and recent releases
  already fixed. `examples/apache/planetgen.conf.example`'s own "HTTPS"
  section previously instructed copying the whole `:80` block (daemon
  process declaration included, which would conflict once duplicated)
  into `:443` and redirecting `:80` unconditionally -- it now excludes
  `/api` from that redirect and mounts `/api` locally on `:80` too,
  reusing one server-wide `WSGIDaemonProcess` declaration instead of two
  conflicting ones. `docs/api.md`'s "Deploying behind Apache" section
  cross-references the same warning.

## [5.46.19] - 2026-09-23

### Fixed
- **`GET /api/sectors/<id>` (`html/sector.py`'s page, and its Sector Map)
  and, under load, unrelated pages sharing the same single-process
  `planetgen-api` WSGIDaemonProcess (`html/system.py` included) were
  timing out / failing outright in production** ("Truncated or oversized
  response headers received from daemon process" and read-timeout errors
  in `planetgen.error.log`, clearing only after an Apache restart -- the
  same failure signature [5.46.16]'s `idx_sectors_center` fix addressed
  for the Galaxy Map). Root cause: `queryDb.phenomena_near_sector` called
  `_placed_phenomenon_rows` with no bounding box at all, so *every*
  `sector_detail` call did a genuine full-table scan across all four
  placed-phenomenon tables (`nebulae`/`asteroid_fields`/`black_holes`/
  `neutron_stars`), pulling every galaxy-placed phenomenon in the entire
  database into Python on every single sector page view. Once a database
  had a non-trivial number of placed phenomena, a handful of concurrent
  sector-page views were enough to hold every one of the API's 5 worker
  threads (and, in turn, its MySQL connection pool, which has no
  checkout timeout) in slow queries at once, starving every other
  request behind them until Apache was restarted.
  `_placed_phenomenon_rows` now takes an optional SQL bounding-box filter
  (the same `BETWEEN`-range-scan technique [5.46.16]'s fix used for
  `sectors`), and `phenomena_near_sector` uses it, padded by the widest
  currently-placed phenomenon radius (a cheap `MAX(radius_ly)` query, not
  a fixed assumption, so it stays exactly as correct for an arbitrarily
  large placed phenomenon as the old unconditional scan). New schema v26
  adds the matching spatial indexes
  (`idx_{nebulae,asteroid_fields,black_holes,neutron_stars}_center`) --
  **run `migrateDb.py` (or `update.sh`/`install.sh`) against any existing
  deployment's database for this fix to actually take effect**; `GET
  /api/health` (see [5.46.17]) will report `schema_current: false` in the
  meantime.

## [5.46.18] - 2026-09-23

### Fixed
- **`GET /api/health` was returning a bare 503 "error" for a database
  that's reachable but has never had `schema.sql`/`migrateDb.py` applied
  to it at all** (no `schema_migrations` table yet -- one step further
  back than "some migrations pending", which it already handled). Caught
  by CI: `test_health_ok` exercised exactly this case by accident (its
  `client` fixture's database starts completely empty) and failed after
  [5.46.17]'s health-reporting change merged. `/api/health` now reports
  this the same way it already reports a stale-but-present schema
  (`200`, `schema_current: false`, a `detail` naming the fix) rather than
  folding it into the "unreachable" 503 case, which is meant for an
  actually-unreachable server. `test_health_ok` now lays the schema down
  first (matching how a real deployment's database always already has
  one by the time its API is queried), and a new test covers the
  never-migrated-at-all case directly.

## [5.46.17] - 2026-09-23

### Fixed
- **[5.46.16]'s own changelog entry claimed a plain app restart applies a
  pending schema migration (e.g. that fix's own `idx_sectors_center`)
  "automatically" -- it doesn't.** `GET /api/health` now reports
  `schema_version`/`schema_current` (and a `detail` message naming the
  fix) by comparing the database's own `schema_migrations` table against
  the code's `SCHEMA_VERSION`, so a deployment that pulled in a
  schema-fixing code change but never actually ran `migrateDb.py` (or
  `update.sh`/`install.sh`) against its database is visible at `/api/health`
  instead of continuing to silently run the old, unmigrated schema -- the
  likely explanation if the same full-table-scan timeouts (and the site
  going down under load) kept happening after [5.46.16]'s code shipped
  but its database was never separately migrated.

## [5.46.16] - 2026-09-23

### Fixed
- **`GET /api/galaxy/view` (the interactive 3D Galaxy Map's live-viewport
  query, added in [5.46.13]) was a genuine full-table scan on `sectors`
  every single call -- confirmed in production as the site going down
  under load (`TimeoutError`/"Truncated or oversized response headers"
  from the WSGI daemon), the exact same failure mode already documented
  and fixed once before for the pre-v22 `/api/search` ([5.35.7], schema
  v22). `sectors.center_x/y/z_pc` had no index
  (`queryDb.galaxy_sectors_in_view`'s own bounding-box `WHERE` clause said
  so explicitly), and this endpoint is hit hard: once server-side on
  every galaxy-map page load (the zoomed-all-the-way-out starting view,
  spanning the whole galaxy) and repeatedly (debounced) as the 3D camera
  moves. New schema v25 adds a composite `idx_sectors_center` index
  (`_db._migrate_v24_to_v25`) so the query can range-scan instead of
  examining every row. **Run `migrateDb.py` (or `update.sh`/`install.sh`,
  which call it) against the database to pick this up** -- restarting the
  app alone does *not* apply it: the API's own connections are read-only
  (`ensure_schema=False`) and never run schema DDL at all, and even a
  read-write connection's `_ensure_schema` only runs `CREATE TABLE IF NOT
  EXISTS`, a no-op against a table that already exists (see [5.46.17]).

## [5.46.15] - 2026-09-23

### Changed
- **`galaxy.py` (the "Galaxy Map" page) now renders the 3D map directly**,
  in place of the flat SVG projection -- rather than living alongside it
  as a separate `galaxy3d.py` page (5.46.13's original approach). The
  flat SVG rendering code (`lib/galaxymap.py`'s former
  `render_galaxy_map_panel` and `static/galaxymap.js`) is removed
  entirely; `lib/galaxymap.py` keeps only the plain Quadrant/Ring
  classification math `galaxy.py`'s own data tables and
  `sector.py`/`browse.py`'s "Quadrant N" links still need.

### Fixed
- **The 3D map's illustrative density cloud was invisible from the
  starting full-galaxy view.** Its points were sized in world-space
  parsecs with perspective attenuation on -- correct for something meant
  to represent real physical size, but a 1-2 pc point shrinks to
  sub-pixel from thousands of parsecs away. Switched to a constant
  on-screen pixel size (`sizeAttenuation: false`) so the cloud stays
  visible at any zoom level.
- **The density cloud didn't read as a recognizable galaxy shape.**
  Points were drawn uniformly across the whole query volume -- at a wide
  view, almost all of that volume is near-empty halo, so only a sparse,
  shapeless scatter of points ever landed somewhere bright. Switched to
  importance sampling from a bulge+disk mixture shaped like the galaxy's
  own real mass distribution, so the cloud now visibly reads as a bright
  core plus a disk (spiral-arm contrast still comes through via each
  point's own real predicted density driving its color).
- **Placed/not-yet-generated sector markers were effectively invisible
  from a wide view** (the same world-space-sizing problem as the density
  cloud above) -- "I don't see the sectors we've generated anywhere."
  Both tiers now use the same constant-screen-size billboard technique
  (recomputed every frame from each marker's own live distance to the
  camera), so a generated sector stays a visible, clickable dot
  regardless of how far the camera currently is, without ever growing to
  dominate the view up close either.

### Added
- **`tests/galaxy_shape_visualizer_cli.py`** -- a diagnostic tool that
  renders the galaxy's real density model as an actual face-on/edge-on
  image (matplotlib), for directly eyeballing whether a set of shape
  parameters produces a recognizable spiral rather than only judging it
  through `relative_density` numbers. Optionally overlays every real,
  already-generated sector's own position when given `--mysql-*`
  connection args.

## [5.46.14] - 2026-09-23

### Changed
- **Every star/planet/moon marker on the System Map now renders its own
  live 3D sphere in place, not just a floating preview beside a click.**
  Previously, only the one planet/moon last clicked got a rotating 3D
  preview -- floated in a small box beside its flat marker (`#sysmap-
  preview`) rather than replacing it, and the whole-system view's star and
  every other unclicked body stayed flat 2D circles regardless. Every
  visible marker (star included) now gets its own sphere, sized and
  positioned to exactly cover -- and read as replacing -- its own flat
  circle, all drawn each frame through one shared WebGL canvas
  (`#sysmap-spheres-canvas`, `lib/systemmap.py` + `static/systemmap.js`)
  via a scissored sub-viewport per marker, rather than one `<canvas>`/
  context per body (which would risk exceeding a browser's cap on
  concurrent WebGL contexts in a crowded system). A star gets its own
  unlit, spectral-color-tinted sphere (`_star_color`, newly exposed as
  `data-color` same as a planet/moon's own class color) plus a matching
  glow shell; falls back to the plain flat marker for a browser that
  can't create a WebGL context at all.

## [5.46.13] - 2026-09-23

### Added
- **Interactive 3D Galaxy Map.** New "Galaxy Map (3D)" page
  (`galaxy3d.py`, linked from the existing flat Galaxy Map) -- a real
  perspective-camera WebGL scene (three.js, `lib/galaxymap3d.py` +
  `static/galaxymap3d.js`) a visitor can rotate, dolly, and click through,
  instead of only ever viewing the galaxy from directly above the disk.
  Because a real 3D camera scales sprite size with distance for free, this
  also fixes the flat map's own "star icon doesn't shrink as you zoom in"
  scaling problem, without any special-case code.
  - **Live viewport queries, not one whole-galaxy payload.** New
    `GET /api/galaxy/view` (`queryDb.galaxy_view`, backed by a new pure
    `planetgen.galaxy.viewport` module) returns, for whatever the
    camera's current view actually covers: real, already-generated
    sectors nearby; real, not-yet-generated sector addresses this
    galaxy's own density model predicts would qualify (exact, out to a
    200 pc cap); and, for the rest of a wider view, a coarse illustrative
    density point cloud. Fetched (debounced) by the page's own
    client-side JS directly from a new browser-facing proxy,
    `galaxy_view.py`, as the camera moves -- never baked into one page
    load the way the flat map's own dataset is.
  - **Logarithmic click-to-zoom.** Left-click zooms in on whatever's
    under the cursor (a sector, a real not-yet-generated address, or
    empty space), right-click zooms out -- both by a step size that
    shrinks the closer the camera already is (big multiplicative jumps
    while zoomed out over the whole galaxy, fine ones once close to a
    single sector), rather than a flat factor that's either too slow to
    cross the galaxy or too coarse to land on one sector.
  - **Sector designation/address, surfaced and copyable.** Clicking a
    real, not-yet-generated address now shows its provisional designation
    and a "Copy CLI command" button with the exact
    `generate.py galaxy --shell K --slot N` invocation to generate it.
  - **`generate.py galaxy --shell K --slot N`.** New single-address
    generation mode (on top of existing `--shell` batch and
    `--center-sector` neighborhood modes) -- generates exactly the one
    sector slot at that address via the existing lazy-generation entry
    point (`ensure_sector_generated`), the direct path from a designation
    copied out of the new 3D map into this script.

## [5.46.12] - 2026-09-23

### Added
- **Zoomable, real-scale diagram on every stellar phenomenon's page.**
  `phenomenon.py` gains a "Diagram" panel (new `lib/phenomenonmap.py` +
  `static/phenomenonmap.js`, reusing `static/mapzoom.js`'s shared zoom/pan)
  drawn directly to astronomical-unit scale: a nebula/asteroid field's real
  `radius_ly` becomes an actual to-scale circle, zoomable from about 1 AU
  up to 1 ly across. A black hole/neutron star (whose real size is
  negligible at this scale) instead shows a small fixed illustrative dot.
- **Supernova remnants are now a full first-class phenomenon type.**
  Previously `supernova_remnants` had no web page at all. `phenomena.py`'s
  listing and `phenomenon.py`'s detail page (morphology, progenitor type,
  age, radius, any compact remnant left behind, galactic orbit) now cover
  it, plus its own real-scale Diagram panel. It has no galaxy-frame
  placement columns of its own, though (unlike the other four phenomenon
  types), so it never appears on the Galaxy Map and can't be a NAV
  endpoint -- `phenomenon.py` shows a short note explaining this instead
  of offering "Navigate from/to here" buttons that would only fail.

### Fixed
- `queryDb.nav_between` would 500 with a raw "unknown column center_x_pc"
  SQL error if ever asked to resolve a supernova remnant as a NAV
  endpoint (its table genuinely has no such column). Now raises a clean
  `ValueError`, same as any other invalid NAV request.

## [5.46.11] - 2026-09-23

### Added
- **Real interactive zoom/pan on the Galaxy Map.** With enough placed
  sectors, the core cluster used to squash into what looked like a single
  dot no matter how many sectors actually existed -- the map was one fixed,
  non-interactive SVG scaled to fit the single farthest-placed sector, and
  `?quadrant=` only cropped that same squashed drawing. Scroll/wheel to
  zoom (centered on the cursor), drag to pan, and `+`/`-`/`Reset view`
  buttons, all the way in to about 100 ly across -- implemented as `viewBox`
  mutations on the existing server-drawn SVG (new shared
  `static/mapzoom.js`, reused as-is by a future phenomenon-diagram zoom;
  `static/galaxymap.js` wires it to the galaxy map's own scale readout).
  The `?quadrant=` crop is unchanged as the map's *starting* view; zoom/pan
  layers on top of it.

### Fixed
- **A marker/label click anywhere on the Galaxy Map stopped navigating**
  partway through implementing the above: capturing the pointer on
  `pointerdown` (needed so a drag that leaves the SVG mid-gesture keeps
  panning) retargeted the resulting `click` event to the `<svg>` itself
  per the Pointer Events spec, so `navform.js`'s delegated
  `closest("[data-nav-target]")` lookup never found the actual marker.
  Deferred `setPointerCapture` until a real drag is detected instead of
  calling it unconditionally on every pointerdown.

## [5.46.10] - 2026-09-23

### Added
- **"Navigate to here" (symmetric with the existing "Navigate from here"),
  and NAV support for standalone phenomena.** `system.py` now offers both
  directions; `phenomenon.py` gains both buttons too (nebulae, asteroid
  fields, black holes, and neutron stars can now be NAV origins/
  destinations, including a full optimal route via adjacent systems, not
  just a direct course). `nav.py`'s origin picker now accepts an
  already-known destination (from a "Navigate to here" link) and carries
  it through to the course instead of re-prompting. `GET /api/nav` gained
  `from_kind`/`to_kind`/`from_type`/`to_type` query parameters for this
  (see `docs/api.md`'s NAV section) -- existing system-to-system callers
  are unaffected.

### Fixed
- **`GET /api/nav` 500'd for any route involving a phenomenon endpoint.**
  `route.positions` could end up with both an int key (a system hop) and
  a string key (a phenomenon endpoint), and Flask's default JSON
  serialization sorts dict keys, which raises `TypeError` comparing an
  int to a string. Fixed by stringifying every key at the JSON boundary
  (`api/routes.py`), confirmed against a running instance and covered by
  a regression test.

## [5.46.9] - 2026-09-23

### Changed
- **Compacted the spread-out page header on System/Sector/Galaxy/Phenomenon
  pages into one line.** The breadcrumb, Octant/binary badges, "Navigate
  from here" button, and nearest-location text used to each be a
  separately stacked, full-width block -- on `system.py` alone that was 4
  lines of near-empty vertical space before the actual content started.
  Added a shared `.page-subhead` flex row (`static/style.css`) and wired
  it into `system.py`, `sector.py`, `galaxy.py`, and `phenomenon.py`.
  Verified in a browser: the header on a binary system's page shrank from
  roughly 650px of vertical space to about 150px.

## [5.46.8] - 2026-09-23

### Changed
- **System Map: better label collision avoidance, and the 3D body preview
  now lives inside the map itself.** A planet/moon label that collided
  with a neighbor used to try only "above"/"below" before giving up and
  hiding the label entirely; it now also tries "right"/"left", then a
  further-out "above"/"below" tier (connected back to its marker with a
  short leader line) before giving up -- in a stress test, this cut the
  hidden-label rate from ~40% to ~16% for a tightly packed cluster, with
  zero label-to-label overlaps either way. Separately, the rotating 3D
  sphere preview (`#sysmap-preview`) used to sit in a fixed sidebar box
  next to the map; it now floats inside the map viewport itself, next to
  whichever marker was just clicked.

## [5.46.7] - 2026-09-23

### Fixed
- **A planet's evolutionary/civilization narrative could contradict its own
  class description.** A planet's fixed per-class flavor text (e.g. Class
  G: "a rocky, barren world with simple life") and its evolutionary/tech
  milestone (`evolution.get_evolutionary_timeline`, up to "Technological
  Civilization") were generated by two completely independent code paths
  -- the planet's own class was never passed into the evolutionary-timeline
  calculation at all, so an old/fast-evolving star, or a forced
  `INTELLIGENT_LIFE=True`, could still report a full technological
  civilization for a planet whose own class text says life there tops out
  at "simple" or "bacterial." Added `program_constants.
  PLANET_CLASS_MAX_LIFE_STAGE`, a per-class ceiling on the highest
  milestone that class's description is consistent with (Class E capped at
  the most minimal stage, F/G at simple/bacterial, L at vegetation/
  multicellular; classes with no explicit life-complexity wording in their
  description stay uncapped), and threads the planet's class through so
  both the natural age-based roll and a forced `INTELLIGENT_LIFE=True` are
  capped by it. Includes a hard invariant assertion and dedicated
  regression tests.

## [5.46.6] - 2026-09-23

### Fixed
- **A binary system's own per-star property tables didn't render.** Each
  star's `###`/`===` section header was joined to its property table by a
  single newline instead of a blank line, so `html/lib/mdconvert.py`'s
  blank-line block splitter lumped the heading and table into one block --
  neither a valid single-line heading nor a valid table -- and rendered it
  as one escaped paragraph of literal `#`/`|` characters instead of a real
  heading plus table. Happened once per star, so every binary system showed
  two broken blocks. Fixed in `systemData.py`'s close- and wide-binary
  rendering paths.
- **Generated systems could place two asteroid belts overlapping each
  other, or a planet's orbit inside an asteroid belt.** Two independent
  causes: (1) a wide (S-type) binary's cross-star clearance check
  unconditionally skipped a trailing asteroid belt when finding a star's
  "outermost" object, so belt-vs-belt (or belt-vs-planet) overlap between
  the two stars' disks was never checked at all; (2) a forced
  habitable-world/explicit-class placement drew its distance uniformly
  within the target zone with no awareness of already-placed belts, and
  could land inside one, or otherwise leave the object list no longer
  sorted by distance -- which the existing overlap correction only ever
  checks between immediate list neighbors, so a resulting overlap with a
  non-adjacent belt went uncorrected. Fixed by making cross-star clearance
  belt-aware (using a belt's own outer edge and a fixed minimum-separation
  threshold, since a belt has no mass for the existing Hill-radius
  criterion to apply to) and by making zone-forced distance selection
  avoid already-placed belt spans. Added an independent, all-pairs overlap
  invariant check (not a re-derivation of the existing correction's own
  formula) to catch any future regression of this kind.
## [5.46.5] - 2026-09-20

### Fixed
- **`src/html/search.py` and `GET /api/search` crashed with a 500 error on any
  text search.** In `queryDb.py`, the SQL LIKE clauses across all five text search
  helpers (`_search_result_sectors`, `_search_result_systems`,
  `_search_result_stars`, `_search_result_planets`, `_search_result_moons`) used
  Python string literals `"ESCAPE '\\'"`. In Python string literals, `"\\"`
  resolves to a single backslash (`\`), passing `ESCAPE '\'` to MySQL. In
  MySQL/MariaDB, `\` is an escape character in string literals, so `\'` was
  parsed as an escaped single quote that left the string literal unclosed,
  causing MySQL syntax error 1064 and triggering a 500 Internal Server Error in
  the Flask API client (`planetGen API error (500): internal server error`).
  Updated all five queries to `"ESCAPE '\\\\'"` so MySQL receives `ESCAPE '\\'`
  and correctly evaluates the escape character as a single literal backslash.

## [5.46.4] - 2026-09-19

### Fixed
- **`generate.py sector`/`galaxy` crashed and discarded an entire sector's
  worth of already-generated systems once the sector ran out of physical
  room.** `SpaceSector.add_system`'s Hill-sphere-based random placement
  raises `ValueError` once a sector's cube has no space left for another
  system without overlapping an existing one's Hill sphere -- a real,
  physically expected outcome once `--num-systems`/`--density` (compounded
  by `--min-habitable` forcing extra large, larger-Hill-sphere stars) asks
  for more systems than a sector can hold at realistic stellar spacing
  (reproduced with `--num-systems 60` in a default 11.5 ly sector, which
  only fits ~35-44). `generate_sector` let that exception propagate,
  aborting the whole run and throwing away every system already
  generated -- including the expensive planet/moon generation behind each
  one. It now stops as soon as one system can't be placed, logs how many
  of the requested systems it actually placed, and returns the sector with
  whatever fit instead of crashing or generating further systems that were
  never going to fit either.

## [5.46.3] - 2026-09-19

### Fixed
- **`generate.py sector --debug` gave no way to see where a sector's generation
  time actually went.** Added `planetgen.util.log.timed_phase`, a debug-only
  context manager that logs `"<label>: <elapsed>ms"` (timestamped, like every
  other `--debug` line) around `generate_sector`'s and `StarSystem.__init__`'s
  major phases -- config building, each system's own generation and placement,
  phenomena generation, and the star/planet/comet/life-data passes inside each
  system -- so a `--debug` run's own output doubles as a per-phase profile with
  no separate profiling flag needed. Also added an opt-in benchmark
  (`PLANETGEN_RUN_PERF_BENCHMARK=1 pytest src/tests/test_sector_generation_perf.py`)
  that captures and aggregates those phase records across several generated
  sectors, cross-checked with a `cProfile` run -- this is what surfaced the
  sector-capacity crash fixed in the next release.

## [5.46.2] - 2026-09-19

### Fixed
- **`generate.py galaxy`'s random-start mode no longer shows a bogus
  single-sector progress bar.** The previous release ([5.46.1]) restored
  the "Sectors" progress bar for `galaxy`'s three modes, but random-start
  mode's own seed sector was wrapped in its own `"Sectors (random start)"`
  task with `total=1` -- a bar that goes straight from 0 to 100% in one
  `advance` call before it can render a meaningful rate or ETA, so in
  practice it was just a single completed bar sitting on screen, not a
  genuine progress display. That was a mistake in how [5.46.1] restored
  the bar: the whole point of the restore was the *galaxy*-level bar
  (sectors being filled across a shell/neighborhood/random-start's
  surrounding radius), not a bar for generating one sector. Random-start
  mode now generates its seed sector with no progress task of its own and
  falls straight through to `run_local_neighborhood`, whose own
  `"Sectors (local neighborhood)"` task (covering every sector generated
  around that seed) is the only bar shown -- matching `--shell`/
  `--center-sector` mode, which never had this problem.

## [5.46.1] - 2026-09-19

### Changed
- **Restored the progress bar for `generate.py galaxy` (`--shell`/
  `--center-sector`/random-start), but not for `generate.py sector`.** A
  prior release removed the progress bar from both commands entirely
  while chasing a flicker/scrolling bug; the bar itself (routed through
  `progress.console.print` so it stays pinned at the bottom with no
  flicker -- see `_generation_progress`'s own docstring) was fine and is
  genuinely useful for a `galaxy` run, which can mean thousands of
  sectors. `run_sector`'s own `--num-sectors` loop stays bar-free, since
  that run is normally short enough that a bar added more noise than it
  was worth.

## [5.46.0] - 2026-09-19

### Changed
- **A qualifying sector can no longer come out completely empty.** A
  sector's own system count and each exotic-phenomenon type's own count
  are independent Poisson draws, so all of them landing on zero at once
  is a real, expected outcome -- more likely the closer local density
  sits to the 1-star-per-sector qualification threshold (e.g. ~13% at a
  mean of 2 systems, ~37% right at the threshold itself, mean 1) -- but a
  sector with nothing in it at all isn't useful to anyone visiting it.
  `generate_sector` now force-adds exactly one system when a
  `--density`-driven sector's own draws (system count, every phenomenon
  type) all came back empty. Only applies when the count came from
  `--density` (explicit or `_BatchDensity`-resolved, as in `galaxy`
  mode); an explicit `--num-systems 0` (including on the plain `sector`
  subcommand) is a deliberate request this never second-guesses.

## [5.45.1] - 2026-09-19

### Changed
- **Cut database round trips for name-uniqueness bookkeeping** -- profiling
  a real save (`SHOW GLOBAL STATUS` query counters before/after) found
  `reserve_body_name`/`confirm_body_name` (run once per planet and once
  per moon -- hundreds of times in a single sector) issuing 4 separate
  SELECT/INSERT/UPDATE round trips each, the large majority of a saved
  sector's total query count. `reserve_body_name` now combines its two
  stateless cross-level collision checks (`sector_name_registry`,
  `system_name_registry`) into one query via `UNION ALL` instead of two
  sequential ones (its `body_name_registry` lookup keeps its own
  unchanged `FOR UPDATE` lock, still queried separately). `confirm_body_name`/
  `confirm_system_name`/`confirm_sector_name` now upsert
  (`INSERT ... ON DUPLICATE KEY UPDATE`, relying on each registry
  table's own `base_name UNIQUE` constraint) instead of a SELECT to
  decide between an INSERT and an UPDATE. Measured ~30-35% fewer total
  queries per planet/moon across repeated benchmark runs (same
  generation code, a throwaway database, and `SHOW GLOBAL STATUS`
  counters compared before/after this change) -- since DB writes
  dominate real generation time (a separate finding: ~70 systems/sec of
  pure in-memory generation vs. ~11-16 systems/sec once real database
  writes are included), this directly speeds up `sector`/`galaxy`
  generation. Every existing name-uniqueness test (collision handling,
  occurrence-count tracking, suffix/diminutive progression) still passes
  unchanged -- the upsert's `ON DUPLICATE KEY UPDATE` branch updates the
  exact same columns (`occurrence_count`, `suffix_index`/
  `diminutive_index`) the old SELECT-then-UPDATE branch did, and never
  touches `first_body_id`/`first_star_system_id`/`first_sector_id`,
  matching the old code exactly.

## [5.45.0] - 2026-09-19

### Added
- **`generate.py galaxy`'s random-start mode (no `--shell`/`--center-sector`)
  gained `--min-start-density`** -- requires the randomly chosen starting
  sector's own real `relative_density` (the same "expected" figure printed
  alongside each saved sector) to be at least the given value before
  accepting it, retried the same way an already-occupied or otherwise
  non-qualifying address already was. Lets an operator skip past the
  galaxy's own vast, sparse outskirts (a plain random start lands there
  most of the time, since a volume-weighted draw favors them) and start
  somewhere with real content to look at -- e.g. `--min-start-density 1.0`
  for at least as dense as the galaxy's own real local density. Only
  applies to random-start mode, and can't be combined with
  `--density`/`--num-systems` (those override every position's density
  uniformly, leaving nothing per-position to compare against).

## [5.44.1] - 2026-09-19

### Changed
- **Removed the "Sectors"/"Sectors (shell N)"/"Sectors (local
  neighborhood)"/"Sectors (random start)" progress bar entirely** from
  `generate.py sector`/`galaxy` -- the previous release only removed the
  nested per-sector bar and tried to fix the outer one's flicker by
  routing prints through it, but the outer bar itself was still visible
  and still wasn't what was wanted. `run_sector`/`run_shell_batch`/
  `run_local_neighborhood`/`run_random_start`/`run_galaxy` no longer take
  or build a `rich.progress.Progress` at all -- every status line is a
  plain `print` again, and `generate.py` no longer imports `rich.progress`.

## [5.44.0] - 2026-09-19

### Changed
- **Removed the per-sector "systems in this sector" progress bar** that
  `generate.py sector`/`galaxy` nested under the outer "Sectors" bar --
  most sectors, especially since the density-gating fix above, hold
  anywhere from zero to a handful of systems, and system generation
  itself is fast, so a bar that flashed on and off again within a single
  frame for nearly every sector added visual noise without conveying
  anything a viewer could actually track. `generate_sector`/
  `generate_and_save_sector_at` no longer take a `progress` parameter at
  all.
- **Fixed the remaining "Sectors" bar fighting with the status text
  printed alongside it, which is what actually caused the flicker/
  scrolling** -- `run_sector`/`run_shell_batch`/`run_local_neighborhood`/
  `run_random_start` now print every "Saved sector ..." status line (and
  the per-sector system/phenomena/density summary) via
  `progress.console.print(...)` instead of the builtin `print`, the
  correct way to write to the console alongside a live `rich.progress.
  Progress` display. Printing directly to stdout while `Progress`'s own
  `Live` region is active fights with its redraws -- each raw `print`
  forced the bar to erase itself, scroll up with the new text, and get
  redrawn at the bottom again, which is what showed up as flicker/
  scrolling on a real terminal. Routed through the shared console
  instead, rich prints each status line safely above the live region and
  leaves the bar itself pinned at the bottom, redrawn in place with no
  flicker -- verified against a real pseudo-terminal (`script`), and
  unchanged (a single plain line at the end, no ANSI live redraw) when
  stdout isn't a real terminal at all (piped to a file, a CI log, etc.),
  which `rich.Console` already detects and handles on its own.

## [5.43.1] - 2026-09-19

### Added
- **`test_galaxy_gen.py` now has two end-to-end tests that run against a
  *real* `generate.py plan` skeleton** (the actual `find_shell_bands`
  scan, not the file's existing `_seed_skeleton` shortcut) rather than an
  explicit `--num-systems`/`--density` that bypasses `_BatchDensity`'s own
  density/gating logic entirely -- every density-related test before this
  did one or the other, so none of them actually exercised the "run
  `plan`, then `galaxy` with neither flag given" workflow the previous
  release's empty-sectors regression slipped through.
  `test_random_start_neighborhood_matches_the_real_skeleton_plan` runs
  `galaxy`'s own default random-start mode -- pick a location, generate
  the nearest sectors out to `--radius-pc` (trimmed to 25 ly here, from
  the real default of 100 ly, to keep the test fast; `-planets` forced so
  each system skips its own planet/moon tree, since system *count* is
  what's under test) -- then independently recomputes, against the real
  stored skeleton, whether every candidate slot in that neighborhood
  should have been saved, and checks the aggregate system count generated
  is within a statistical band of what the plan's own density predicted.
  `test_shell_batch_generates_nothing_beyond_the_real_skeletons_outer_edge`
  is its deterministic companion: a shell chosen well past the real
  skeleton's own discovered edge must generate exactly zero sectors.
  Both fail against the pre-fix code (confirmed by hand, reverting
  `generate.py` locally and re-running).

## [5.43.0] - 2026-09-19

### Fixed
- **`generate.py galaxy` (`--shell`, `--center-sector`, and the default
  random-start mode) generated and saved a real, empty (0 systems, 0
  phenomena) sector row for every not-yet-occupied slot it visited, once
  `generate.py plan`'s skeleton existed to drive per-sector density.**
  Only the lazy, visit-triggered `ensure_sector_generated` was actually
  gating generation on `galaxy_shell_band`'s stored candidate bands and
  the exact `predicted_star_count >= 1.0` threshold; the three batch/
  neighborhood modes in `run_shell_batch`/`run_local_neighborhood`/
  `run_random_start` never consulted either, so every slot below that
  threshold -- most of a realistic galaxy's volume, off the spiral arms/
  disk plane -- still got a sector saved, almost always with 0 systems
  once its (correctly tiny) relative density was fed through the Poisson
  draw. A galaxy shell/neighborhood run could come back "full of empty
  sectors" as a result, especially the default random-start mode, whose
  volume-weighted starting-shell pick favors the sparse outer galaxy.
  `_BatchDensity.resolve` (`generate.py`) now applies the same
  band-then-exact-density gate `ensure_sector_generated` already used,
  returning `None` for a non-qualifying slot so every caller skips it
  instead of generating and saving it; the admin web UI's "generate more
  sectors around this one" action (`generate_sector_neighborhood`) picked
  up the same skeleton-driven density and gating, having previously used
  a flat `--num-systems 10` for every sector regardless of position at
  all.
- Every sector-generating command (`sector`, and `galaxy`'s `--shell`/
  `--center-sector`/random-start modes) now prints, right after saving a
  sector, how many star systems of each spectral class and phenomena of
  each type it actually holds, plus that sector's actual vs. expected
  star density (`1.0` = real local stellar density) -- enough to
  sanity-check a generation run from its own console output, without a
  separate database query.

## [5.42.0] - 2026-09-19

### Added
- **A 3D body-preview sphere on the System Map.** Clicking a planet or
  moon in `system.py`'s System Map now also redraws `#sysmap-preview`: a
  small rotating three.js sphere (reusing the same vendored build the
  Sector Map uses) shaded by the body's own class color, banded with a
  tilted ring for a gas giant, and wrapped in a fresnel-glow atmosphere
  shell -- tinted by surface temperature -- whenever the body actually has
  one. The info panel also gains "Atmosphere," "Surface composition," and
  "Surface temperature" fields (`planets`/`moons.atmosphere`/
  `composition`/`surface_temperature_k`, already generated and stored,
  just not previously surfaced here). The true-position SVG diagram itself
  is unchanged -- this is an appearance preview alongside it, not a
  replacement.

## [5.41.0] - 2026-09-19

### Changed
- **The Sector Map is now a real WebGL scene instead of a CSS 3D
  illusion.** `html/lib/starmap.py` no longer builds one `<div>` per star/
  outline edge/compass arrow positioned via CSS `transform-style:
  preserve-3d` -- it now serializes the same position/size/color/label
  data it always computed into a `<script type="application/json">`
  block, and a rewritten `html/static/sectormap.js` renders it with
  three.js (a real perspective camera, GPU-billboarded sprites for stars/
  nebulae/asteroid fields/black holes/neutron stars, and a wireframe
  outline for the sector's wedge or fallback cube) -- proper perspective/
  occlusion, a scale bar that now accounts for the panel's own responsive
  size instead of assuming a fixed 320px scene, and a `<noscript>` link
  list plus a hidden screen-reader-accessible button list (a canvas has
  no focusable children of its own the way the old per-star `<div
  role="button">`s were) so the sector's systems/phenomena stay reachable
  without JavaScript or with a keyboard/screen reader alike. three.js is
  vendored at `html/static/vendor/` (bundled and minified from the `three`
  npm package, not loaded from a CDN) so `html/lib/page.py`'s existing
  `Content-Security-Policy: default-src 'self'` needs no exception for it.
- **Every star/phenomenon marker on the Sector Map now navigates via
  `data-nav-target`/`data-nav-params` (`static/navform.js`) instead of a
  plain `href`**, catching the WebGL rewrite above up to the
  no-address-bar-params convention `5.40.0` (below) introduced for the
  rest of `html/` after this branch had already diverged from it.

## [5.40.0] - 2026-09-19

### Changed
- **Every navigational link in `html/` now posts its parameters as
  hidden form fields instead of putting them in a `<a href="page.py?
  db=...&id=...">`'s query string** -- `db`, a sector/system/phenomenon
  id, a search filter, a wiki-upload/admin-action field, and the like no
  longer show up in the browser's own address bar. `lib/fmt.py`'s new
  `post_link` builds a same-effect, no-JS-required `<form method="post">`
  submit button, styled (`static/style.css`'s `.link-btn`) to be visually
  indistinguishable from the plain link it replaces; `lib/page.py`'s new
  `nav_params`/`nav_multi_params` are what a page reads a followed link's
  params back with (a POST body when present, else the GET query string,
  so a bare `QUERY_STRING`-only smoke test still works). The one
  exception is a Galaxy Map/Sector Map/NAV Map marker plotted inside an
  `<svg>` (a `<form>` can't nest inside one) -- those still navigate via
  a real, focusable `<a>`, now carrying `data-nav-target`/
  `data-nav-params` (`fmt.data_nav_params`) instead of an `href` query
  string, intercepted by the new `static/navform.js` (loaded on every
  page) to post the same hidden form a click on any other link would.
  Every such marker still has a plain, no-JS-required row in a table
  below its own map, so a marker click is never the only way to reach
  something. `index.py` no longer redirects to `browse.py?db=...` for
  this deployment's one database (a redirect's `Location` URL would
  itself show `db` in the address bar) -- it calls `browse.handler`
  in-process instead and renders the result directly, guarded by
  `browse.py`'s own `if __name__ == "__main__":` so `browse.py` reached
  directly is unaffected. This makes every page un-bookmarkable/
  un-shareable by URL and, for a map marker specifically,
  JavaScript-dependent -- a deliberate trade-off (see `lib/page.py`'s
  module docstring) for keeping database names, record ids, and search
  terms out of browser history, address bars, and referrer headers.

## [5.39.1] - 2026-09-19

### Fixed
- **CI was red on every run.** `test_db_persistence.py` still asserted
  `_db.SCHEMA_VERSION == 22` (and `migrate_database(...) == 22`) in
  sixteen places, left over from before the v23 wiki-publishing and v24
  name-uniqueness migrations bumped `SCHEMA_VERSION` to 24; each
  assertion now just compares against `_db.SCHEMA_VERSION` instead of a
  stale literal. Fixing that surfaced a second, previously-masked bug in
  the same file: eight of those tests reset `schema_migrations` to an
  older version to replay later migrations, but never dropped the v22
  search-index migration's indexes first, so replaying it against a
  freshly-bootstrapped (already-v24) test database hit a "Duplicate key
  name" error -- a new `_drop_v22_search_indexes` helper (alongside the
  existing `_drop_v20_trajectory_columns`/`_drop_v21_phenomenon_columns`)
  fixes that.
- **Two `test_galaxy_gen.py` shell-batch tests raised `TypeError`.** Their
  `_fake_generate_sector` test doubles didn't accept the `progress`
  keyword the nested-progress-bar feature added to the real
  `generate_sector`, so `monkeypatch`ing it in broke as soon as
  `generate.py galaxy --shell` started passing one.
- **A MySQL connection-pool leak was exhausting the test server's
  `max_connections` partway through a full test run.** `_db.py` caches one
  `PooledDB` per distinct connection config in a module-level dict that's
  never evicted -- fine for the handful of long-lived databases a real
  deployment ever points at, but the test suite's own `mysql_config`
  fixture hands every single test a uniquely-named throwaway database, so
  each test left one more pool (and its `mincached` real connection)
  behind for the rest of the process's life. A new `_db.close_pool()`,
  called from that fixture's teardown once its database is dropped for
  good, closes and discards the pool immediately instead.

## [5.39.0] - 2026-09-18

### Added
- **Galaxy-wide name uniqueness.** `planetgen/names/uniqueness.py`
  tracks every sector/system/planet-or-moon base name ever generated
  (`sector_name_registry`/`system_name_registry`/`body_name_registry`,
  schema v24) and decorates a colliding name instead of letting two rows
  anywhere in the database share a display name -- Greek/Roman letters
  for a sector or system colliding with its own kind, a diminutive prefix
  for a system colliding with its sector, and a companion suffix for a
  planet/moon colliding with anything. `_db.py`'s `insert_sector`/
  `insert_star_system`/`insert_planet`/`insert_moon` all consult and
  update these registries now. `src/dedupeNames.py` is a new one-off
  script to decorate any duplicate names an existing database already
  has from before this feature existed.
- **A "generate more sectors around this one" admin action** on
  `html/sector.py`, for any already galaxy-placed sector -- fills in
  every not-yet-generated sector within a 100 ly sphere around it
  (`POST /api/sectors/<id>/generate-neighborhood`), the same
  local-neighborhood logic `generate.py galaxy --center-sector` already
  used from the CLI, now reachable from the web interface.
- **Live progress bars for `generate.py sector`/`galaxy`.** Both now show
  nested `rich.progress` bars (elapsed and estimated-remaining time) --
  an outer "Sectors" task for galaxy's shell-batch/local-neighborhood/
  random-start modes (and sector's own `--num-sectors` loop), and an
  inner "systems in this sector" task nested under it.

### Changed
- **Removed the separate write-capable MySQL config.** `config.json`'s
  `mysql_write` section (and the `PLANETGEN_MYSQL_WRITE_USER`/
  `_PASSWORD` env vars layered over it) is gone -- the Flask API's
  `WRITE_MYSQL_CONFIG` now simply reuses `MYSQL_CONFIG`, so there's a
  single account (`config.json`'s `mysql` section /
  `PLANETGEN_MYSQL_*`) for the generation CLIs, `install.sh`/
  `migrateDb.py`, and the API's reads and writes alike -- give that one
  account whatever grants the most demanding caller needs.

## [5.38.0] - 2026-09-18

### Fixed
- **A binary system's stored `markdown_content` mixed wikitext template
  blocks into it (and vice versa for `wikitext_content`)** -- only ever
  for the secondary star's own "Star Data" table/age sentence, and only
  for a binary system (`+binary_system`), never a single star.
  `StarSystem.__init__` gave the secondary star its own `copy.deepcopy`d
  `SystemConfig` (to force `LARGE_STAR` off on it without affecting the
  primary), which also forked `MARKDOWN` into its own disconnected copy.
  `_db.py.insert_star_system` renders *both* formats from one generated
  system by toggling `system_config.MARKDOWN` and calling `str(star_system)`
  again -- a toggle that only ever reached the primary/shared config, never
  the secondary's own deep copy, leaving its data table stuck rendering in
  whichever format was current at generation time regardless of which one
  was actually being requested afterward. The secondary star now shares
  the system's own `SystemConfig` object throughout (`LARGE_STAR` is still
  forced off for it, just transiently, restored right after).
- **A binary system's own name/title had "Binary System" literally baked
  into it** (e.g. "Sol Binary System" instead of "Sol") for a 'close'
  (P-type) pair -- `doubleStar.BinaryStarProxy.name` (the proxy's `.name`
  stands in for the whole system's own identity, stored as
  `star_systems.name`) now takes the primary star's own bare name, the
  same convention a 'wide' (S-type) pair's system name already used.
- **The Sector Map and System Map labeled a binary's two stars "Primary"/
  "Secondary"** instead of a real name -- both now read "&lt;name&gt; A"/
  "&lt;name&gt; B", matching the generator's own existing convention for
  the secondary star's *stored* name (`"<primary name> B"`).
- **A wide (S-type) binary's own System Map diagram never showed the
  secondary star's planets at all, and text/markers routinely rendered
  cut off or overlapping** -- both stars' planets, and the pair's own
  real separation (routinely tens to thousands of AU, per
  `wideBinary.py`), used to share one log-scaled radial pixel budget.
  Since a wide pair's real separation is so much larger than either
  star's own planetary system, that shared scale either crushed both
  stars' planets down near the frame's center to make room for the real
  separation, or pushed the stars themselves (and their planets' label
  text) toward the frame's outer edge and straight off the visible
  canvas. `lib/systemmap.py` now gives the primary its own full-budget
  "system" scene (its own planets only) with a companion marker for the
  secondary that swaps to the secondary's *own* full-budget scene when
  clicked -- the same "drill into it" pattern a planet with moons already
  used, one level up.
- **The Sector Map's default zoom could shrink every star dot into an
  illegible, barely-visible smear** for a galaxy-placed sector -- a
  sector's on-shell wedge wireframe is deliberately allowed to draw well
  past the fixed-size scene (so it can be dragged/zoomed into fully), but
  the *default* zoom used to fit that wedge's own extent alongside every
  star/cloud's, in one combined list. A wedge can be many times wider
  than the scene regardless of how tightly clustered the sector's own
  stars actually are, so a single oversized wedge could drag the whole
  default zoom down to its own floor, shrinking every star dot along with
  it -- confirmed by rendering a realistic sector this way and finding
  its dots reduced to a handful of barely-visible pixels, easily read as
  "almost nothing rendered". The default zoom now fits the real content
  (every star/cloud) on its own whenever there is any, falling back to
  fitting the wedge/cube outline only for a genuinely empty sector (the
  original problem that outline-fitting behavior was written to fix).
- **`GET /api/nav`'s cross-sector ("galaxy" scope) route could time out**
  once the generated galaxy grew large. `navGraph.build_knn_adjacency`
  built its routing graph with an O(n&sup2;) "compare every point to every
  other point" pass -- unnoticeable for one sector's own handful of
  systems, but this same function also runs over *every* system in *every*
  galaxy-placed sector generated so far for a cross-sector route, a set
  that only ever grows as more of the galaxy gets visited/generated. Now
  builds an in-memory 3D k-d tree and queries each point's true k nearest
  neighbors through it instead (O(n log n)) -- the exact same resulting
  graph, computed roughly 10x+ faster already at a couple thousand
  systems, comfortably under a second even at 50,000.
- The Sector page now also lists every nearby exotic phenomenon
  (`queryDb.phenomena_near_sector`, the same set the Sector Map's own
  clouds/points are drawn from) in its own table below the systems one,
  each row linking to that phenomenon's `phenomenon.py` detail page --
  previously only visible on the map itself, with no plain listing.

## [5.37.0] - 2026-09-18

### Added
- **Wiki publishing is wired up.** The standalone `wikiClient` library
  (`src/wikiClient/`, unified in [5.34.0]) is now actually called: an
  "Upload to Wiki" form on both `html/system.py` and `html/sector.py`
  (admin sessions only) publishes to whichever of Wiki.js/MediaWiki is
  configured deployment-wide, backed by new `POST /api/systems/<id>/wiki`/
  `POST /api/sectors/<id>/wiki` write routes. A system publishes its
  already-generated `markdown_content`/`wikitext_content`; a sector (which
  has no persisted page of its own) gets one built fresh at upload time
  from its own current detail. `config.json` gains a `wiki` section
  (`wikijs.base_url`/`.api_token`, `mediawiki.base_url`/`.username`/
  `.password`, each independently optional -- either, both, or neither
  backend may be configured at once, letting an upload choose "the wiki
  of their choice" when both are), read by the new `GET /api/wiki-config`
  endpoint the two forms use to know which backend(s) to offer.
  `star_systems.wikijs_url`/`mediawiki_url` (present in the schema since
  [5.34.0] but never populated) and a new `sectors.wiki_url` column
  (schema v23, `stellarObjects._db._migrate_v22_to_v23`) record where each
  page ends up; once set, `html/system.py`'s Description section is
  replaced by a link to the wiki page (opening in a new tab) instead of
  the locally rendered/source view, and `html/sector.py` shows the same
  kind of link. `html/admin.py` also gains a small form to manually set or
  clear a sector's `wiki_url` directly (`PATCH /api/sectors/<id>`), for a
  sector with a hand-written page from outside this app.

## [5.36.0] - 2026-09-18

### Added
- **A sector-then-system picker on the Nav page, reachable with no
  starting system already known.** Previously `nav.py` only ever worked
  when arriving via a specific system's own "Navigate from here" button
  (`?from=<id>` required); the sidenav's new "Nav" link now reaches it
  with nothing chosen yet, and a two-step `<select>` picker (every sector,
  then every system in the chosen one -- `GET /api/sectors` then
  `GET /api/sectors/<id>`) sets `from=` the same way arriving via
  `system.py` already did. The cross-sector half of the destination picker
  (choosing `to=` once an origin is known) gets the identical two-step
  sector-then-system cascade in place of its old plain numeric
  destination-system-id field.
- **A list and detail page for exotic phenomena** (nebula/asteroid
  field/black hole/neutron star) -- this project's first per-phenomenon
  pages; previously a phenomenon had no page of its own at all, only a
  hover tooltip on the Sector Map/Galaxy Map. `phenomena.py` lists every
  phenomenon across every sector, regardless of galaxy placement (`GET
  /api/phenomena`, paginated); each row links to `phenomenon.py`'s full
  detail view (`GET /api/phenomena/<type>/<id>`, each type's own real
  columns -- a nebula's `composition`/`formation_cause`, a black hole's
  `mass_solar`/`spin`/`has_accretion_disk`, etc. -- via new `queryDb.
  list_phenomena`/`phenomenon_detail`). Both pages are linked from the
  sidenav; the Sector Map's and Galaxy Map's own phenomenon markers
  (`lib/starmap.py`/`lib/galaxymap.py`) now click through to the same
  detail page instead of only showing a tooltip.

## [5.35.7] - 2026-09-18

### Fixed
- **The Search page timed out** (`TimeoutError`/`urllib.error.URLError`
  surfaced through `html/lib/apiclient.py`, rendered as an unexpected
  error page) once the database grew past a trivial size. `GET
  /api/search` always runs its full facet-count and autocomplete query set
  up front, on every visit, regardless of whether any filter is active
  (`queryDb.search`) -- ten `GROUP BY`/`SELECT DISTINCT ... ORDER BY`
  queries, none of them backed by an index on the column they group,
  filter, or sort by (`stars.yerkes_class`, `planets`/
  `moons`.`planet_class`/`body_type`/`life_chemical`,
  `asteroid_belts.density`, and every table's own `name`), so each one was
  a genuine full-table scan/sort. New schema v22
  (`_migrate_v21_to_v22`/`schema.sql`) adds the missing indexes; a name
  *term* search (`LIKE '%text%'`, a leading wildcard) isn't sped up by any
  of them -- that would need a FULLTEXT index, out of this fix's scope --
  but the facet counts and autocomplete lists that run unconditionally on
  every visit are. `apiclient.py`'s own request timeout also widened
  15s -> 30s as a second line of defense, not a replacement for the real
  fix.

## [5.35.6] - 2026-09-18

### Fixed
- **Every galaxy-generated sector's system count was flatly stuck at 10**,
  regardless of where it actually sits in the spiral galaxy -- a bulge
  sector and a sparse outer-disk sector generated the same way.
  `ensure_sector_generated` (the visit-triggered lazy-generation path)
  already correctly drove system count from the galaxy skeleton's real
  position-based `relative_density`, but `generate.py galaxy`'s own
  batch/local-neighborhood/random-start generation -- how every sector in
  a real deployment actually gets made -- never consulted it at all: every
  sector in a run shared one flat CLI value, defaulting to `num_systems =
  10` when neither `--density` nor `--num-systems` was given.
  `validate_shared_generation_args` now leaves both unset for `galaxy`
  mode specifically in that case (`sector` mode, which has no galaxy
  position to compute a density from, is unaffected); a new `_BatchDensity`
  helper resolves each sector's own `relative_density` from the stored
  skeleton (fetched once, reused for the whole run) and feeds it through
  exactly the way `ensure_sector_generated` already does, in
  `run_shell_batch`/`run_local_neighborhood`/`run_random_start` alike. An
  explicit `--density`/`--num-systems` still applies uniformly for the
  whole run, unchanged.
- Audited the actual system-placement code path (`SpaceSector.add_system`/
  `_random_position`, via `generate_sector`'s `for system, cfg in
  zip(...): sector.add_system(...)` loop) for whether it could silently
  place fewer systems than the (now real, skeleton-driven) requested
  count -- it can't: a sector too crowded to fit the next system's minimum
  Hill-sphere separation raises `ValueError` after
  `SECTOR_MAX_PLACEMENT_ATTEMPTS` tries rather than skipping it, so an
  under-delivered density would already be a loud failure, not a silent
  one. No code change needed for this part; noted here since it was the
  other half of what was reported.

## [5.35.5] - 2026-09-18

### Fixed
- **A galaxy-placed sector's Sector Map opened looking almost empty/broken**
  -- a couple of giant wireframe edges crossing the visible crop instead of
  a wedge shape, with any star near the outline's own edge invisible
  outside the fixed, non-panning viewport. The wedge wireframe (`lib/
  starmap.py`'s `_wedge_edges_px`) is deliberately allowed to extend well
  past the fixed 320px scene (the wedge's angular patch doesn't coincide
  with a cube's flat sides), but the map always *started* at `zoom = 1`
  regardless -- confirmed by rendering the real output in a browser and
  comparing that default against manually zooming all the way out, which
  showed the exact same content correctly. `render_map_panel` now computes
  a `_default_zoom` from the actual extent of everything being drawn
  (wedge/cube vertices, every star/cloud) and starts (and "Reset view"
  returns to) that fitted zoom instead of a flat default; `sectormap.js`'s
  own `MIN_ZOOM` floor widened to match.
- **The Galaxy Map rendered as a dense, unreadable smear of overlapping
  ring labels for any sector placed far from the core**, with its own dot
  sitting right at the visible circle's edge -- reproduced directly with a
  sector at `shell_index` ~1400, which implied 157 fixed-shell-width Rings,
  each drawn as its own guide circle + label, all crammed into the same
  480px panel. `lib/galaxymap.py`'s `_rings_to_show` (which picks the
  map's *scale*, so a far sector still fits) is now decoupled from how
  many ring guides `_ring_elements` actually *draws*: past
  `_MAX_RINGS_DRAWN` (10), it switches from one guide per literal
  fixed-shell-width Ring to 10 evenly-spaced distance markers spanning the
  same range -- still real, accurate distance labels, just no longer
  cluttering the map once there would be too many to read. The common
  near-core case (few real Rings) is unaffected -- confirmed with a
  regression render.

## [5.35.3] - 2026-09-18

### Fixed
- **Site `<title>`/browser-tab text always said "planetGen"**, ignoring
  `config.json`'s own `site_name` (e.g. "Molten Aether Starmap") that
  every other page-title path already honored. `lib/page.py`'s `render()`
  hardcoded the literal string instead of calling `load_config()` the way
  `api_base_url`/`base_url` already do; `index.py`'s own "no databases"
  title had the identical hardcoded string. Both now read `site_name` from
  config.
- **The "pick a database" landing page is gone.** This project deploys as
  one branded starmap per vhost now (`config.json`'s `site_name`/
  `api_base_url`), so a picker whose choice is realistically always length
  1 just added an extra click/page load in front of every visit.
  `index.py` now redirects straight to `browse.py` for the first database
  `GET /api/databases` returns, regardless of how many exist; every other
  page's breadcrumb (`browse.py`/`galaxy.py`/`sector.py`/`search.py`/
  `system.py`) drops its now-pointless leading "Databases" link, and the
  sidenav's own "Databases" item is removed (there is no longer a picker
  page for it to reach).

## [5.35.1] - 2026-09-18

### Fixed
- **Planet/moon/star infobox fields showed literal `<sup>7</sup>` markup
  instead of a superscript 7.** `tabledisplay.py`'s scientific-notation
  formatters (`format_body_distance`/`format_star_mass`/`format_star_radius`/
  `format_star_luminosity`) emit real HTML (`"5.3 × 10<sup>7</sup> km"`),
  correct for `system.py`'s static table cells (inserted unescaped on
  purpose) but wrong for `lib/systemmap.py`'s interactive System Map: it
  carries the same strings through `data-*` attributes that
  `static/systemmap.js` reads back with `.textContent` (deliberately never
  `innerHTML`, so database-derived values can never execute as markup) --
  which shows a `<sup>` tag as literal text instead of rendering it. Added
  `tabledisplay.to_plain_text`, converting the one `<sup>N</sup>` pattern
  into real Unicode superscript digits (`10⁷`), and applied it at every
  `data-*`-building call site in `systemmap.py` (distance, mass, radius,
  luminosity); `system.py`'s own raw-HTML table cells are untouched.

## [5.35.0] - 2026-09-17

### Changed
- **Unified `systemGen.py`/`sectorGen.py`/`galaxyGen.py`/`galaxyPlan.py`/
  `phenomenonGen.py` into a single `generate.py` script.** Those five
  root-level scripts are removed; every generator in this project is now
  reached through one program and one subcommand: `generate.py
  system|sector|galaxy|plan|phenomenon [options]`. Each subcommand
  accepts exactly the option surface its old standalone script offered
  and saves to the same database -- this is a pure consolidation, not a
  behavior change. The five scripts used to import each other
  (`sectorGen.py` called into `systemGen.py`, `galaxyGen.py` called into
  `sectorGen.py`, and so on); that logic now lives together in
  `generate.py`'s own sections (system -> sector -> galaxy -> galaxy
  skeleton -> exotic phenomena -> the unified CLI itself), calling each
  other directly instead of through cross-module imports.
  `setup.py`'s `py_modules`/console-script entry point were updated to
  match (`planetgen=generate:main`, replacing the old `systemgen`/
  `sectorgen` scripts), and `src/tests/test_examples.py`/
  `test_sector_gen.py`/`test_galaxy_gen.py` (the tests that imported the
  removed modules directly) now import `generate` instead.

## [5.34.0] - 2026-09-17

### Changed
- **Unified `wikijs`/MediaWiki publishing into a single `wikiClient`
  library.** `src/wikijs/` (the standalone Wiki.js GraphQL client) is
  replaced by `src/wikiClient/`, which exposes one `WikiClient` object
  (`backend="wikijs"` or `backend="mediawiki"`) whose `create_page` works
  the same way regardless of target -- both backends are create-only and
  share one exception hierarchy (`WikiClientAuthError`/
  `WikiClientPageExistsError`/`WikiClientRequestError`). The former
  `WikiJsClient` logic moves in unchanged as `wikijs.WikiJsBackend`; new
  alongside it is `mediawiki.MediaWikiBackend`, a from-scratch, stdlib-only
  MediaWiki Action API client (Bot Password login, CSRF token, `action=edit`
  with `createonly=1`) -- this project previously had no MediaWiki API
  client at all, only a wikitext *text format* option. Nothing in the app
  calls either backend yet (still an open item, see `docs/TODO.md`); this
  is purely the shared library those still-`# TODO` call sites
  (`routes.py`, `config.py`, `appconfig.py`, `system.py`) will build on.
  Tests renamed/moved to match (`test_wikiclient_wikijs(_integration).py`)
  and a `test_wikiclient_mediawiki(_integration).py` pair added, plus
  `test_wikiclient_client.py` for the new dispatch facade.

## [5.33.0] - 2026-09-17

### Added
- **Galaxy random-start mode.** `galaxyGen.py` run with neither `--shell`
  nor `--center-sector` now picks a uniformly random (by volume, not by
  shell index -- see new `_pick_random_shell_index`) not-yet-occupied
  sector address within a real Milky-Way-scale galaxy (`--max-shell`,
  defaulting to the shell nearest new `program_constants.GALAXY_RADIUS_PC`,
  15,000 pc), generates it, then falls straight through to
  local-neighborhood mode's own logic to generate every sector within
  `--radius-pc` of it too (defaulting to new `program_constants.
  RANDOM_START_NEIGHBORHOOD_RADIUS_LY`, 100 ly, in every direction) -- so
  a bare `galaxyGen.py` with no arguments at all creates a whole small
  starmap around a fresh, randomly chosen starting point in one run.
- **Science-based exotic phenomena as part of ordinary sector
  generation.** Every generated sector (`sectorGen.py` directly, or via
  `galaxyGen.py`) now also seeds a realistically sparse population of
  `phenomenonGen.py`'s own seven phenomenon types -- black holes, neutron
  stars, nebulae, supernova remnants, rogue planets, interstellar comets,
  standalone asteroid fields -- sampled independently per type via a
  Poisson draw whose mean is a cited real (or, where flagged, a
  deliberately conservative order-of-magnitude) astrophysical rate per
  star system (new `program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM`:
  Lamberts et al. 2018 for black holes, Sartore et al. 2010 for neutron
  stars, Frew & Parker 2010 for nebulae, Diehl et al. 2006 for supernova
  remnants, Sumi et al. 2011/Mroz et al. 2017 for rogue planets), scaled
  by however many star systems the sector actually ended up with -- a
  denser sector (`--density`) gets proportionally more, and most sectors,
  realistically, get none at all, exactly as real space this size usually
  is. New `sectorGen.generate_sector_phenomena`.
- **Hill-sphere-safe placement generalized to exotic phenomena.** A
  standalone black hole/neutron star is a real, stellar-mass gravitating
  body, so `spaceSector.py`'s existing star-to-star Hill-sphere placement
  logic (`hill_radius_ly`/`required_separation_ly`) now applies to it
  too, via new `SpaceSector.add_phenomenon`/`SectorPhenomenonEntry`:
  neither can ever land within a neighboring star system's or another
  compact remnant's own Hill sphere -- and, since the same neighbor list
  now includes already-placed phenomena, a star system placed afterward
  can't land inside a black hole's/neutron star's Hill sphere either.
  Every other phenomenon type has no comparable gravitational footprint
  at this generator's own scale and is placed at a random point in the
  sector's cube instead. Persisted via schema v21: `black_holes`/
  `neutron_stars` gain the same `sector_id`/`center_x/y/z_pc`/
  `galactic_radius_pc` galaxy-frame placement shape v18 gave `nebulae`/
  `asteroid_fields` (`_migrate_v20_to_v21`), with the new position
  converted directly from each phenomenon's own real in-sector offset
  rather than independently re-randomized; the other five phenomenon
  types' own `sector_id` column -- present since v16 but never actually
  populated -- finally gets wired through `_db.save_phenomenon`. Surfaced
  in the web interface everywhere nebulae/asteroid fields already were:
  `queryDb.py`'s phenomena queries, the Sector Map (a small glowing point
  instead of a translucent cloud, since a compact remnant's own physical
  radius is negligible at this scale), and the Galaxy Map.

### Fixed
- **`SectorPhenomenonEntry.distance_to` crashed** with `TypeError:
  'SectorPhenomenonEntry' object is not iterable`. `spaceSector.
  distance_between` only recognized `SectorSystemEntry` via `isinstance`,
  so calling `.distance_to()` on the new phenomenon-entry class (added
  alongside the Hill-sphere generalization above) tried to iterate the
  entry object itself instead of reading its `.position`. Generalized to
  duck-type on "has a `.position`" instead, covering both entry classes
  (and any future one) without an `isinstance` check needing to know
  about each concretely. Caught by this release's own new test suite,
  run against a real MariaDB server rather than skipped for lack of one
  (9,265 tests, 0 failed, 0 skipped) -- not by any pre-existing test.

## [5.32.0] - 2026-09-16

### Added
- **Proper two-body (barycentric) trajectories for binary stars, planets,
  and moons.** Every orbital pair previously modeled only "the lighter
  body orbits a fixed primary" -- a poor approximation for a binary
  secondary (sampled at 0.1-0.8x the primary's mass) and for a moon near
  this generator's own mass cap (up to 1/10 its parent planet's mass,
  approaching real "double planet" ratios like Pluto/Charon). Existing
  "relative position" columns (`planets`/`moons.position_x/y/z_km`,
  `star_systems.binary_mutual_position_x/y/z_km`) are unchanged -- still
  the true separation a large amount of existing physics (insolation,
  Hill sphere, tidal locking) depends on. New columns instead add the
  ORBITED body's own small "reflex offset"/"wobble" away from its
  nominal fixed point (new shared `utils.calculate_reflex_offset`
  helper): both binary members now visibly orbit their common barycenter
  (`star_systems.binary_primary_position_*_km`/
  `binary_secondary_position_*_km`, plus the constant
  `binary_secondary_mass_fraction` and, for a 'close' pair, the
  additional `binary_planetary_wobble_*_km` from its own circumbinary
  planets), a planet-hosting star gets its own wobble from the combined
  pull of its planets (`stars.reflex_offset_*_km`), and a moon-hosting
  planet gets its own wobble from its moons
  (`planets.reflex_offset_*_km`). `_db.advance_orbital_phases` recomputes
  all of these fresh on every run (a cheap derived value, no independent
  update-guard interval of its own), introducing a new correlated-
  subquery `UPDATE` technique to sum a parent's pull from multiple
  children in one set-based statement. Persisted via schema v20 (new
  columns on the three already-existing `star_systems`/`stars`/`planets`
  tables) and a `_migrate_v19_to_v20` migration step that backfills real
  values for every pre-existing row (unlike v17's migration, every value
  here is fully derivable from data already stored).
- **Standalone Wiki.js GraphQL client for page publishing.** New
  `wikijs` package (`src/wikijs/`): `client.py`'s `WikiJsClient` creates
  pages via Wiki.js's GraphQL API using Personal API Tokens, returning a
  `WikiPage` dataclass (`id`/`path`/`title`/a constructed `url`);
  `exceptions.py`'s hierarchy (`WikiJsError`/`WikiJsAuthError`/
  `WikiJsPageExistsError`/`WikiJsRequestError`) distinguishes auth
  failures (401/403, or a GraphQL-level auth error) from a duplicate-path
  rejection (message-text matching, since Wiki.js has no stable error
  code for it across versions) from any other failure. Stdlib-only
  (`urllib.request`/`urllib.error`/`json`, no external dependency) and
  create-only -- `create_page` never checks for an existing page first,
  letting Wiki.js reject a duplicate naturally rather than racing a
  check against it. Not yet wired into the app -- TODO comments mark
  where `apiclient.py`/`routes.py`/`system.py`/`appconfig.py`/
  `config.py` will eventually call it. (Superseded in `5.34.0` by the
  unified `wikiClient` library, which folds this client in unchanged as
  one of two backends.)

## [5.31.0] - 2026-09-12

### Added
- **Galaxy-frame placement for nebulae/asteroid fields, and web-interface
  map overhaul.** Nebulae/asteroid fields previously had no location at
  all -- `sector_id` existed but was explicitly "reserved for a future
  sector-context encounter" and always NULL. `phenomenonGen.py
  --sector-id` now actually places one: a galaxy-frame sphere
  (`center_x/y/z_pc`/`galactic_radius_pc`, schema v18) centered near an
  already galaxy-placed sector (`stellarObjects._db.
  compute_phenomenon_placement`), not a sector-relative offset -- a
  nebula can span up to 200 ly, far larger than a single ~11.5 ly sector,
  so it's a real sphere that may overlap several sectors' cubes, or none.
  `sector_id` itself becomes a real (if non-authoritative) "nearest
  sector" convenience link (`ON DELETE SET NULL`, not `CASCADE`). New
  `queryDb.phenomena_near_sector` (bounding-sphere-vs-sector-cube overlap
  test) and `galaxy_placed_phenomena`, exposed via `/api/sectors/<id>`'s
  new `phenomena` key and a new `GET /api/galaxy/phenomena` endpoint.
  `html/lib/starmap.py`'s Sector Map now draws a translucent spheroid
  cloud for every nebula/asteroid field near that sector (nebula color by
  type; a mottled tan/dark texture for an asteroid field); `html/
  lib/galaxymap.py`'s Galaxy Map plots every placed one as a small fixed-
  size dot instead -- a phenomenon's own real size is shown where it
  actually fits (the Sector Map), not misleadingly as a galaxy-scale dot.
- **System Map: real body positions, not a schematic.** `html/
  lib/systemmap.py`'s System Map previously always drew every planet due
  east of its star on one fixed concentric-orbit diagram, ignoring each
  body's own real orbital angle entirely. It's now a true top-down plot:
  every star/planet/moon/belt sits at its actual angle (from
  `planets`/`moons.position_x/y_km`) and a distance from its own anchor
  that's log-scaled into one shared pixel budget spanning the whole
  scene -- a log scale is what lets a close binary's ~0.05 AU separation
  and an outer planet's 30+ AU orbit coexist on the same diagram without
  either collapsing to a point or blowing out the frame. A binary pair's
  two stars are placed at their real mass-weighted offsets from the
  system's own barycenter (new `binary_mutual_position_x/y/z_km` in
  `queryDb.system_detail`'s response, split by each star's `mass_kg`) --
  for a 'wide' (S-type) pair, this also removes the old map's "primary
  only, see the tables below for the secondary" limitation: both stars
  and each one's own independent planets now draw in the same
  true-position scene. Asteroid belts are now full rings around their own
  anchor (the natural true-position shape) rather than a one-directional
  shaded band. A once-only pairwise-repulsion pass (`_relax_markers`)
  nudges apart any two markers real placement happened to put too close
  together -- real position first, decluttering only where legibility
  actually needs it.
- **Star-bound elliptical/parabolic comets.** Adds the two comet orbit
  types the existing `InterstellarComet` (always unbound/hyperbolic)
  couldn't represent: a periodic elliptical comet and a single-apparition
  parabolic one, both bound to their star. New `keplerMotion.py` solves
  the two-body Kepler/Barker equation so an eccentric or parabolic orbit
  moves with the correct non-uniform angular speed (Kepler's second law)
  instead of a linear phase advance; new `cometData.Comet` generates
  realistic orbital elements (Jupiter-family/Halley-type/long-period
  subtypes) and scales coma/tail activity chance by perihelion distance
  rather than a flat roll. Wired all the way through: `StarSystem` gets
  an independent `comets`/`secondary_comets` list (kept separate from
  planets/asteroid belts, since a comet's distance is a continuously
  varying current position, not a fixed orbital-slot distance) and a new
  `SystemConfig.COMETS` tri-state flag; `advance_comet_orbits`
  (`updateOrbits.py`) advances each comet's orbital anomaly and
  recomputes its position on every run; `queryDb.system_detail`/
  `html/system.py` gain a parallel comets table. Persisted via schema v19
  (`comets`/`comet_composition` tables, `_migrate_v18_to_v19`).
- **Five project-wide bugs, found in a general bug-check pass:**
  `StarSystem.validate_system`'s orbital-overlap correction
  double-counted `last_planet.distance`, roughly doubling the corrected
  distance for any overlap involving an asteroid belt; `utils.
  years_to_time_string` built its total from a 365.25-day year but
  decomposed it back with a plain 365-day divisor, leaking the
  discrepancy into every displayed orbital period's "days"/"hours" (e.g.
  1.0 year rendered as "1 year and 6 hours"); the admin session cookie
  (`html/api/auth.py`) was scoped to `Path=/api`, so it was never sent
  back to the CGI admin pages served from the deployment's own document
  root, locking every admin page immediately after a successful login;
  three latent bugs in `test_query_db_planets_moons.py` itself (a
  nonexistent `StarSystem.name`, `Planet.radius_km` vs. the real
  `radius`, a missing `body_type` guard before touching
  `AsteroidBelt.moons`); and `advance_comet_orbits` now batches its
  per-row updates into one `executemany` call instead of one `execute`
  per row.

### Fixed
- **CI: schema v18's null-together CHECK on `nebulae`/`asteroid_fields`
  broke `test_migrate_v*` on real MySQL 8.0** (`pymysql.err.
  OperationalError: (3959, "Check constraint 'nebulae_chk_2' uses column
  'center_x_pc', hence column cannot be dropped or renamed.")`) --
  MySQL 8.0 refuses `DROP COLUMN` on a column an anonymous CHECK still
  references, and the migration-rollback test helper had no name to drop
  it by. Both CHECKs are now named explicitly
  (`chk_nebulae_placement`/`chk_asteroid_fields_placement`) rather than
  left anonymous; `_migrate_v17_to_v18` now also adds the same named
  constraint (via the portable `ADD CONSTRAINT ... CHECK`/
  `DROP CONSTRAINT`, not MySQL-only `DROP CHECK`), which it had
  previously omitted -- a database migrated (not freshly created) would
  otherwise silently lack this invariant's enforcement. Caught locally by
  running the full suite against a real MariaDB server (which tolerated
  the original anonymous CHECK's `DROP COLUMN` fine, masking this) rather
  than against MySQL 8.0 as CI does. The workflow's own test matrix also
  gained `fail-fast: false`, since the default had been silently
  cancelling the 3.9 job the moment 3.12 failed, hiding whether a failure
  was version-specific.

## [5.30.0] - 2026-09-12

### Added
- **Standalone asteroid fields, the seventh exotic phenomenon.**
  `phenomenonGen.py --type asteroid-field` generates a field of asteroid
  debris drifting in open interstellar space (new `asteroidFieldData.
  AsteroidField`) -- physically the same object as an in-system
  `AsteroidBelt` (density + mineral composition), just without a host
  star/orbit, reusing `AsteroidBelt`'s own composition-generation logic
  (now extracted into shared `asteroidData.generate_asteroid_composition`/
  `format_composition_summary` functions rather than duplicated).
- **Galactic-orbital motion for every standalone exotic phenomenon.**
  A black hole/neutron star with no owning system, a nebula, a supernova
  remnant, a rogue planet, an interstellar comet, and a standalone
  asteroid field are all still gravitationally part of the galaxy even
  though none is bound to any specific star -- each now gets the same
  `galactic_orbital_speed_kms`/`_period_gy`/`_phase_deg`/
  `_min_update_interval_years` quartet a lone `Star` has (new shared
  `utils.generate_galactic_orbit_fields`/`format_galactic_orbit` helpers,
  also adopted by `Star`/`BinaryStarProxy`/`BlackHole`/`NeutronStar` for
  consistency), advanced over real elapsed time by `updateOrbits.py`/
  `_db.advance_orbital_phases` the identical way a star's already is.
  `advance_orbital_phases` now returns a name-keyed dict rather than a
  positional tuple, since the set of tables it advances keeps growing.
  Persisted via schema v17 (new columns on six pre-existing tables plus
  the new `asteroid_fields`/`asteroid_field_composition` tables) and a
  `_migrate_v16_to_v17` migration step.
- **Fixed: an anchored black hole/neutron star silently lost its identity
  on reload.** `_db.load_star_system`'s single-star branch always called
  `Star.from_dict`, with no dispatch on the owning `stars.yerkes_class`
  marker (`'BH'`/`'NS'`) and no query against the `black_holes`/
  `neutron_stars` satellite tables -- a system saved via `phenomenonGen.py
  --anchor-system` reloaded (via `queryDb.py` or the Flask API) as a
  generic `Star` carrying a nonsensical Yerkes class, missing every
  remnant-specific field, and rendered with `Star`'s own paragraph text
  instead of the remnant's. New `_db._load_single_star` now dispatches
  correctly, reusing the same satellite-row-plus-base-row combination
  `_binary_proxy_row_to_dict` already does for a close binary's merged
  proxy.
- **`planetgen/generation/phenomena_plausibility.py`, a statistical anomaly
  finder for all seven exotic phenomena** -- the same two-tier design as
  the existing planet-focused `plausibility.py` (analytically-derived hard
  invariants, e.g. an event horizon radius must match the Schwarzschild
  formula for its own mass, gated by `test_phenomena_plausibility.py`;
  Tukey's-fences statistical outliers plus category-frequency comparisons
  against each phenomenon's own configured chance, reported for human
  review via the new `src/tests/phenomena_plausibility_cli.py`, never
  asserted exactly).

## [5.29.0] - 2026-09-12

### Added
- **Exotic stellar phenomena, via a new, separate `phenomenonGen.py` CLI.**
  Six phenomena -- black holes, neutron stars, nebulae, supernova remnants,
  rogue planets, and interstellar comets -- can now be generated on demand,
  each grounded in real astrophysics (Schwarzschild radius for black
  holes; NICER-measured neutron star mass/radius ranges and ATNF-catalog
  pulsar spin/field populations; the four standard ISM nebula classes;
  Sedov-Taylor blast-wave expansion for supernova remnant age/size;
  'Oumuamua/Borisov-informed interstellar comet speed/composition).
  Deliberately **not** wired into `systemGen.py`/`sectorGen.py`'s normal
  per-slot generation odds -- `StarSystem._generate_planets` never
  produces one; they're reachable only through `phenomenonGen.py`'s own
  `--type` choice (uniformly random among all six when omitted). A black
  hole or neutron star (new `compactRemnant.py`, subclassing `Star` the
  same way `doubleStar.BinaryStarProxy` does) can optionally anchor a full
  `StarSystem` via `--anchor-system` -- real pulsar planets exist (PSR
  B1257+12) -- reusing all of `StarSystem`'s existing orbit-placement/
  rendering/serialization logic unchanged; its zero-or-near-zero
  luminosity naturally collapses the habitable zone to (0, 0) AU and the
  disk-physics planet-count ceiling to zero, matching the real rarity of
  confirmed planets around compact remnants, without any special-casing.
  Persisted via six new tables (schema v16: `black_holes`/`neutron_stars`
  as satellite tables extending a `stars` row when anchored,
  `nebulae`/`supernova_remnants`/`rogue_planets`/`interstellar_comets`
  always standalone) and a `_migrate_v15_to_v16` bookkeeping-only
  migration step (the six tables are brand new, so no existing table
  needed an `ALTER TABLE`).

## [5.28.0] - 2026-09-12

### Added
- **`queryDb.py`: `planets`/`moons` CLI subcommands.** Closes the gap
  `docs/TODO.md`'s "Open items" > "Search" tracked: the `systems`
  subcommand only ever filtered by star type/sector, with no way to ask
  this CLI "every Class D planet smaller than Earth" the way the web/API
  faceted search (`GET /api/search`, `queryDb.search`) already could.
  New `list_planets`/`list_moons` (mirroring `list_systems`'s shape) take
  an exact `--class`, a `--min-radius-km`/`--max-radius-km` range, and a
  `--sector-id`/`--system-id` scope, reusing the existing
  `_append_size_clause` helper the faceted-search result panels already
  share -- new `_body_filter_clause` factors the class/size/sector/system
  `WHERE` fragment the same way `_systems_filter_clause` already does for
  `systems`, so the two CLI-side query functions can't drift out of sync
  with each other. `moons` rows also report their parent planet's name,
  since a moon's own name alone doesn't say which planet it orbits.

### Fixed
- **README.md's version badge had drifted 3 releases stale** (5.24.0 while
  `_version.py`'s `__version__` was already 5.27.0), caught by hand during
  a deploy-readiness check with no CI signal at all. New
  `src/tests/test_version_sync.py` asserts the README badge and
  `CHANGELOG.md`'s own top entry both match `__version__`, so this can't
  silently recur.
- **CI's `dependency-audit` job had no explicit `setuptools` upgrade step**,
  leaving it exposed to whatever `setuptools` version happens to ship
  preinstalled on the runner's Python image -- caught locally via
  `pip-audit` flagging `PYSEC-2026-3447` against a preinstalled 79.0.1.
  `.github/workflows/ci.yml` now upgrades `setuptools` alongside `pip`
  before installing this project's own dependencies, the same way the
  `test` job's setup already keeps `pip` itself current.

## [5.27.0] - 2026-09-12

### Added
- **S-type (wide) binary star systems.** `+binary_system` previously only
  ever generated a P-type (close/circumbinary) pair -- the two stars merged
  into one effective star (`doubleStar.BinaryStarProxy`) for planet
  placement. A new `+wide_binary`/`-wide_binary` option (random if omitted)
  selects the other real binary configuration instead: an S-type pair,
  separated by tens to thousands of AU (log-uniformly sampled, matching the
  real, roughly log-normal spread of observed wide-binary separations),
  where each star keeps its own separate identity -- its own mass,
  luminosity, habitable zone -- and hosts its own independently-generated
  planets (new `doubleStar`-sibling module `wideBinary.py`'s
  `WideBinaryPair`). Each star's maximum stable planetary orbit is capped
  by Holman & Wiegert's (1999) empirical critical-semi-major-axis formula
  for the companion's long-term perturbation, and a Gladman (1993)
  mutual-Hill-radius check additionally prunes either star's outermost
  planet if the two stars' own disks would otherwise gravitationally
  encroach on each other -- a rare safety net for tight/eccentric pairs,
  not the common case. The pair's own orbital eccentricity is sampled from
  a realistic "thermal" distribution (unlike the close pair, a wide pair
  never tidally circularizes) and feeds both the stability formula and the
  reported periapsis/apoapsis separation, though (like every other orbit
  this generator tracks) the pair's live position/phase-advance tracking
  stays circular -- a deliberate, documented simplification consistent
  with the rest of the engine. Persisted via new `star_systems`/`stars`/
  `asteroid_belts` columns (schema v15) and a `_migrate_v14_to_v15`
  migration step; also fixes a latent bug the new columns exposed in
  `_db.advance_orbital_phases`, whose old single combined `UPDATE` could
  never have advanced a wide pair's mutual-orbit phase even after this
  release, had it not been caught -- now two independently-guarded
  `UPDATE`s.

### Fixed
- **Several tests that don't pin `BINARY_SYSTEM` could intermittently fail
  once merged against [5.26.0]'s real, spectral-class-dependent binary
  chance.** A companion star landing on a system those tests otherwise
  treat as single (or as a fixed planet count) could throw off an
  unrelated assertion -- most seriously, an S-type (wide) pair's own
  `a_crit_au` stability ceiling can make a forced `HABITABLE_WORLD`/
  `ASTEROID_BELT` guarantee geometrically impossible for an especially
  luminous host (an O-type supergiant's habitable zone can sit beyond any
  sampled companion's stability limit), which `test_full_matrix.py`'s
  full-star-type sweep and several of `test_systems.py`'s tri-state-flag
  tests surfaced; `test_disk_physics.py`, `test_db_persistence.py`, and
  `test_api.py` had similar exposure via an unexpected second star's own
  planet count/DB rows. Pinned `BINARY_SYSTEM=False` in each, following
  the same precedent already established when [5.26.0] itself pinned it
  in four other tests -- binary-vs-single was always incidental to what
  each of these was actually testing.
- **CI's `test (3.9)` job failed on every PR, unrelated to whatever the PR
  actually changed.** `test_planet_physics_fixes.py` called
  `statistics.correlation`, added in Python 3.10 -- this repo's CI matrix
  still runs a `3.9` job. Replaced with a small `_pearson_correlation`
  helper (matches `statistics.correlation` exactly; verified against it
  directly) used by both `test_atmospheric_pressure_correlates_positively_
  with_gravity` and `_spearman_correlation`'s own rank-based call.

## [5.26.0] - 2026-09-11

### Changed
- **`SystemConfig.BINARY_SYSTEM` now follows the same tri-state contract
  as every other flag (`HABITABLE_WORLD`, `ASTEROID_BELT`, etc.):
  `None` (the default) is no longer treated as "always single."** It now
  rolls real chance instead, from new
  `program_constants.BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS` --
  keyed by the primary star's own spectral letter, since real
  stellar-multiplicity surveys consistently find companionship rate
  rising with primary mass rather than sitting at one flat rate: ~26% for
  M dwarfs (Duchene & Kraus 2013) up through ~44% for solar-type F/G/K
  (anchored to Raghavan et al. 2010's 46%; Duchene & Kraus's own review
  groups F/G/K together at 44+/-2%) to ~90% for O-type primaries (Moe &
  Di Stefano 2017's 94+/-14%). New `StarSystem._should_generate_binary`
  (called from `__init__`, replacing the old flat `if self.system_config.
  BINARY_SYSTEM:` check) looks this up against `self.primary_star.type[0]`
  -- run *after* the primary star already exists, specifically so its
  real, already-rolled spectral type can drive the roll. `True`/`False`
  still force the outcome exactly as before; only `None`'s meaning
  changed, from "never" to "real chance for this star." `SystemConfig.
  BINARY_SYSTEM`'s own docstring updated to match.
- Four tests that generate systems without pinning `BINARY_SYSTEM`
  (two in `test_space_sector.py`'s name/position round-trip and
  recipe-fallback-reload coverage, one more in `test_space_sector.py`'s
  file save/load round trip, one in `test_serialization.py`'s
  single-star full-object-graph round trip) were relying on the old
  "`None` always means single" behavior to keep specific expected
  names/types deterministic, or (the `test_serialization.py` one) for
  `reloaded.star is reloaded.primary_star` to hold -- true only for a
  single star, since `__init__` never repoints `primary_star` at the
  `BinaryStarProxy` it reassigns `star` to for a real binary. All four
  now pin `BINARY_SYSTEM=False` explicitly, since binary-vs-single was
  always incidental to what each was actually testing.

## [5.25.0] - 2026-09-11

### Changed
- **`StarSystem.estimate_num_objects` now derives its planet/belt ceiling
  from real protoplanetary-disk physics instead of an arbitrary curve fit
  to stellar mass.** The old formula (`BASE_MAX_SYSTEM_OBJECTS *
  (1 + log10(solar_masses))`, base 15) had no grounding in orbital
  dynamics and no relationship at all to the mutual-Hill-radius spacing
  rule `validate_system` enforces ([5.24.0]) -- two disconnected dials
  governing "how many" and "how far apart," tuned independently by feel.
  The replacement, `StarSystem._estimate_max_objects_from_disk_physics`,
  walks outward from the same inner-edge distance `_generate_planets`
  itself seeds its first slot at, and at each step: computes the local
  *isolation mass* an oligarchic-growth embryo would reach there
  (Lissauer 1993; Kokubo & Ida 2000, 2002 -- new `utils.isolation_mass_kg`,
  closed-form-solved the same way `_mutual_min_distance_au` is, since the
  embryo's own Hill radius depends on its own still-unknown mass), from
  the Minimum Mass Solar Nebula's real solid surface-density profile
  (Hayashi 1981 -- new `utils.mmsn_surface_density_gcm2`,
  `physical_constants.MMSN_SOLID_SURFACE_DENSITY_SOL_GCM2 = 7.0 g/cm^2` at
  1 AU falling off as `distance^-1.5`, jumping `SNOW_LINE_ICE_BOOST_FACTOR
  = 30/7` beyond the snow line where ices condense), scaled for this
  star's own disk-mass budget (new `utils.disk_surface_density_scale`,
  `program_constants.DISK_MASS_STELLAR_MASS_EXPONENT = 1.8` --
  mm-continuum disk-demographics surveys, Andrews et al. 2013; Pascucci
  et al. 2016, find real disk dust mass scales roughly as
  `M_star^1.8-2.7`, not logarithmically), then advances by that same
  embryo's own mutual-Hill-radius feeding zone (the *same*
  `program_constants.MUTUAL_HILL_RADII_SEPARATION = 10` and kappa/clamp
  algebra `_mutual_min_distance_au` uses, so the count estimate and the
  spacing rule that will later constrain actual placement are provably
  consistent with each other) and counts a slot. The walk terminates at
  the disk's outer edge -- new `program_constants.
  DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER = 18` times the star's own snow
  line (new `utils.snow_line_au`, `physical_constants.
  SNOW_LINE_AU_AT_1_LSUN = 2.7`, the same `sqrt(luminosity)` shape
  `calculate_habitable_zone` already uses) -- deliberately *not*
  `star.system_perimeter` (that's the star's own galactic-tidal Hill
  sphere, tens to hundreds of thousands of AU; real disks are truncated
  far short of it by viscous spreading/photoevaporation) -- or
  `ABSOLUTE_MAX_SYSTEM_OBJECTS` isolation-mass slots, whichever comes
  first. Not every oligarch survives as a final planet: real N-body
  integrations of the subsequent giant-impact phase (Chambers 2001) show
  most merge or get ejected, so the raw slot count is scaled by new
  `program_constants.GIANT_IMPACT_SURVIVAL_FRACTION = 0.4` (tuned toward
  the middle of that literature's own range, empirically checked to keep
  a solar-mass star's typical resulting count in the same well-tested,
  playable range this generator already verified via repeated full-suite
  runs) before being returned. `estimate_num_objects`'s own override
  contract (`PLANETS`/`NUM_ORBITS`/`MAX_PLANETS`) is entirely unchanged --
  only what `max_objects` means changed. The resulting shape now tracks
  real demographics better than the old mass-scaling formula did: cool
  low-mass stars (whose smaller mutual-Hill spacing packs oligarchs more
  tightly per unit distance -- the real, observed TRAPPIST-1-style
  "compact multis favor small stars" pattern) come out *more*
  planet-rich on average than hot, luminous, high-mass stars (whose
  correspondingly larger isolation masses claim proportionally more of
  their own, larger disk per embryo, and whose short main-sequence
  lifetimes and intense UV output make real, confirmed planets around
  O/B-type stars genuinely rare) -- the reverse of the old formula's
  "bigger star, more objects" curve, and a better match to what's
  actually been observed. `BASE_MAX_SYSTEM_OBJECTS` removed (no longer
  referenced); `ABSOLUTE_MAX_SYSTEM_OBJECTS` unchanged, still the same
  hard safety cap.
- New `src/tests/test_disk_physics.py`: unit coverage for the new
  `utils` helpers directly (snow-line `sqrt(luminosity)` scaling,
  disk-density-scale monotonicity, MMSN falloff/snow-line jump,
  isolation mass matching the literature's own ~0.05-0.1 Earth-mass
  figure at 1 AU) plus `estimate_num_objects`'s override contract and the
  `ABSOLUTE_MAX_SYSTEM_OBJECTS` cap, exercised via `MAX_PLANETS=True`
  (deterministic, no `random.randint` draw) the same way test_systems.py
  already does.

## [5.24.0] - 2026-09-11

### Changed
- **Orbital spacing between adjacent planets now uses their *mutual* Hill
  radius, not either planet's own individual one.** The previous rule
  (`min_orbit_distance = 5 x this planet's own Hill radius`) was a
  reasonable approximation, but the real orbital-dynamics stability
  literature (the analytically rigorous two-planet minimum of `2*sqrt(3)`
  mutual Hill radii, Gladman 1993; a recommended ~8-10x margin for
  longer-term N-body stability, Chambers, Wetherill & Boslough 1996 and
  Smith & Lissauer 2009) expresses this in terms of the *pair's* combined
  mass and average distance instead. New `utils.mutual_hill_radius_m`
  (`((a1+a2)/2) * ((m1+m2)/(3*M_star))^(1/3)`) and
  `program_constants.MUTUAL_HILL_RADII_SEPARATION = 10` (the safer end of
  the literature's recommended range) back `StarSystem.
  _mutual_min_distance_au`, which `validate_system`'s planet-planet
  spacing check now calls instead of `max(planet.min_orbit_distance,
  last_planet.min_orbit_distance)`. Solved in closed form for the exact
  minimum distance rather than evaluated once at the pre-correction
  position and added on top: since the mutual radius depends on the
  *average* of both distances, that naive approach understates the
  requirement once the correction actually moves one of them -- confirmed
  by a real `assert_no_orbital_overlap` failure before the closed-form
  version replaced it. Reclassification (`planetPhysics.
  reconcile_zone_and_class`, [5.22.0]) can change a planet's own mass,
  which can in turn invalidate a spacing decision already made against
  its predecessor -- `validate_system` now retries the mutual-distance
  check (bounded to 3 iterations; converges in practice within 2) after
  any reclassification triggered by its own push. `Planet.
  min_orbit_distance` itself is unchanged and still single-body -- it
  remains the right tool for a *moon's* own orbital limit around its
  parent (`planetPhysics.generate_moons`), a different physical question
  from planet-to-planet spacing. `assert_no_orbital_overlap` (test_systems.py)
  updated to check the same mutual-radius formula, with a relative (not
  fixed-1e-9) tolerance -- the two independent computations of the same
  quantity can differ at the floating-point-noise level even when both
  are correct, and that noise scales with the (sometimes very large, for
  a massive star's own extreme systems) distances involved.
- **`html/search.py`/`GET /api/search` gained a min/max size filter for
  stars, planets, and moons.** `queryDb.search`'s new `sizes` argument
  (`{"star", "planet", "moon"} -> (min_km, max_km)`, either bound
  optional) filters/joins alongside the existing class/body/life-
  chemistry tags and per-entity name search, via a new shared
  `_append_size_clause` helper across `_search_result_stars`/
  `_search_result_planets`/`_search_result_moons` -- a continuous
  quantity like size has no discrete set of values to offer as a facet,
  so it's its own query-parameter pair
  (`star_min_radius_km`/`star_max_radius_km`, likewise `planet_`/
  `moon_`) rather than a tag. `GET /api/search` validates and parses
  these (`routes._parse_size_range`); `html/search.py` gained matching
  number-input fields, active-filter chips (e.g. "Planet size: 5,000–
  8,000 km"), and a Radius column on the Stars/Planets/Moons result
  tables. `queryDb.py`'s own CLI (`sectors`/`systems`/`near`
  subcommands) still has no `planets`/`moons` equivalent at all -- see
  `docs/TODO.md`'s "Open items" > "Search" for that narrower, separate,
  still-open gap.

### Fixed
- **`GET /api/health` could crash instead of returning a clean `503`.**
  `queryDb.open_readonly` raises a bare `SystemExit` (correct for its own
  CLI callers) when the configured MySQL server is unreachable --
  `SystemExit` is a `BaseException`, not an `Exception`, so left
  uncaught it would propagate straight through Flask's request dispatch
  (and the WSGI worker handling it) instead of becoming any HTTP
  response at all, from *every* route that opens a connection via
  `routes.get_db()`, not just `/health`'s own explicit check. `get_db()`
  now catches it and re-raises an `ApiError` (503 -- the same status
  `/health` already wanted to report), which the app's existing
  `ApiError` handler turns into the usual JSON error response for every
  other route, and which `/health`'s own `except Exception` catches
  directly. Regression test builds its own Flask app against a
  deliberately-unreachable config (a closed local port), so it runs
  without needing a real MySQL server the way every other API test does.
- **Sector Map star dots were plotted as if their own sector-local axes
  already ran parallel to the galaxy frame's.** The wedge outline and
  "Galactic Center" compass arrow (`_wedge_edges_px`/`_compass_html`)
  were always computed directly from galaxy-frame quantities
  (`sectors.center_x/y/z_pc`, `sector_wedge_vertices_pc`) and so were
  always correct on their own terms; star dots
  (`star_systems.position_x/y/z_mpc`) were plotted directly in the same
  scene without ever being rotated into that frame -- a design
  convention this project documents (`docs/design/
  galaxy-coordinate-system.md`'s "Cube orientation" section: radial-
  outward local `+Z`, projected-galactic-north local `+X`) but never
  actually applies at generation time. Rather than rotate stored
  positions, `starmap.py`'s new `_rotate_to_galaxy_frame` applies that
  convention at render time, computed fresh from the sector's own stored
  `center_x/y/z_pc` -- reusing `stellarObjects.sectorGeometry.
  cube_orientation`, the exact same basis that module already computes
  as the tangent-plane frame for this sector's own wedge vertices, so no
  new stored orientation column was needed. A sector with no galaxy
  placement (`center_pc=None`) keeps the previous unrotated behavior
  unchanged. New `src/tests/test_starmap.py` covers the rotation's
  length-preservation, the on-galactic-axis degeneracy case, and that
  `render_map_panel` actually renders a different on-screen position for
  a placed vs. unplaced sector.

## [5.23.0] - 2026-09-11

### Added
- **`config.json`: one unified deployment config file, replacing
  `webconfig.json`.** Every entry point in this project (generation CLIs,
  the Flask API, the `html/` CGI browser) used to read its own scattered
  `PLANETGEN_*` environment variables, each with its own hardcoded
  default -- fine per-variable, but it meant a deployment that just wants
  "one MySQL server, one account, one API base URL" still had to set half
  a dozen `SetEnv`/`EnvironmentFile` lines to get there. `webconfig.json`
  existed to solve exactly this for the web interface, but only ever
  covered `site_name`/`base_url` plus three `db_*` placeholders that
  predated the MySQL port and were never wired to anything.
  `planetgen.util.appconfig.load_config()` replaces it: a single
  `config.json` at the repo root, deep-merged onto built-in defaults, now
  covering the read-only and write-capable MySQL connections, the control
  schema name, the database-listing prefix, the API's rate limits, the
  admin cookie's `Secure` flag, the debug-page toggle, and the site's own
  name/base URL/API endpoint -- see `docs/config.md` for the full field
  list. Every `PLANETGEN_*` environment variable still works and still
  takes precedence over `config.json` (needed for, e.g.,
  `planetgen-orbits@.service`'s per-instance
  `PLANETGEN_MYSQL_DATABASE=%i`); `config.json` only adds a place to set
  the shared defaults once instead of repeating them everywhere.
  `config.json.example` (repo root) is the committed template;
  `config.json` itself is gitignored, next to the `webconfig.json` entry
  it replaces.

## [5.22.0] - 2026-09-11

### Fixed
- **`validate_system` could strand a planet outside the zone its class
  needs, and the sequential orbit-spacing loop had no outer bound.**
  Found via full `StarSystem` stress testing rather than the existing
  per-class plausibility tooling (`physical_plausibility_cli.py`), which
  only ever constructs isolated planets directly and never exercises the
  sequential placement loop or `validate_system` at all. A planet's class
  was chosen once, early, from its *initial* estimated distance -- but
  `validate_system`'s orbital-overlap correction could later push that
  distance arbitrarily far out (Hill-radius-based minimum spacing scales
  with a planet's own distance, so it compounds geometrically across a
  many-planet system) without ever re-checking whether the already-chosen
  class still made physical sense there. Measured before this fix: 10% of
  a 400-system sample had at least one ecosphere-class planet (M/K/N/etc.)
  stranded outside its own zone -- e.g. a "Class M, Earth-like world"
  ending up hundreds of thousands of AU out at ~30K -- concentrated almost
  entirely in O/B-type stars (92% of occurrences), none in G/K/M dwarfs.
  - **`planetPhysics.reconcile_zone_and_class(planet, primary_mass_kg,
    distance_override=None)`**: re-derives a body's zone from its
    *current* distance and, only if its existing class is no longer valid
    there, regenerates the class and everything derived from it
    (composition/radius/mass/density/atmosphere/period/gravity/orbital
    motion) -- the same thing real orbital migration does to a body's
    actual final conditions, not just its originally-assumed ones.
  - **`StarSystem._reconcile_moved_planet`** calls it every time
    `validate_system` moves a top-level planet's distance, for the planet
    itself and each of its moons: a moon's zone is always its parent's
    (`generate_moons`' `zone_override`), and a reclassified parent can
    come out with a different mass, which changes a moon's own
    period/position/rotation (Kepler's third law around the parent) even
    when the moon's own class didn't need to change.
  - **The sequential placement loop never bounded how far a system could
    grow.** `StarSystem._generate_planets` (split out of `__init__` so it
    can be retried -- see below) now stops adding slots once the next
    one would land beyond `star.system_perimeter` -- the star's own Hill
    sphere *relative to the galaxy*, already computed for an analogous
    purpose elsewhere (`spaceSector.py`, keeping neighboring systems'
    spheres of influence from overlapping) but never enforced during
    planet placement itself. A system that runs out of stable room this
    way ends up with fewer planets, the same outcome a real
    protoplanetary disk of finite extent would produce, rather than
    letting the geometric compounding above run unbounded.
  - **Reconciliation can occasionally reclassify away the specific body a
    requested `HABITABLE_WORLD`/`ASTEROID_BELT` guarantee was relying on**
    (previously silently masked by the same bug -- an invalid class
    sitting in the wrong zone still counted as satisfying it).
    `StarSystem.__init__` now retries the whole placement (same star,
    fresh positions and object count) up to the new
    `program_constants.MAX_SYSTEM_GENERATION_ATTEMPTS` (8) when a
    requested guarantee isn't met afterward, rather than silently
    dropping it. A smaller, complementary fix
    (`StarSystem._distance_within_zone_with_margin`) reserves a safety
    margin against the fixed `MIN_ASTEROID_BELT_SEPARATION` nudge when
    placing a forced-habitable or explicit-slot-class planet, reducing
    (not eliminating -- a large neighboring planet's own Hill-radius push
    has no fixed size to margin against, which is what the retry loop is
    for) how often the retry actually triggers.
  - Verified via a 500-system randomized stress test (mixed
    `HABITABLE_WORLD`/`ASTEROID_BELT`/`BINARY_SYSTEM` flags): zero orbital
    overlaps, zero misplaced ecosphere-class planets, zero guarantee
    failures, repeated across multiple runs.

## [5.21.0] - 2026-09-11

### Added
- **Ecosphere-zone classes now generate at a class-appropriate distance
  within the habitable zone, instead of every class sharing the same
  distance-blind draw.** Resolves `docs/TODO.md`'s "Class K (Mars analog)
  ... generated at the same zone-midpoint orbital distance as Class M"
  open item -- this generator previously picked a planet's *class* only
  after its *distance* was already fixed by unrelated orbital-spacing
  logic, so real-world position (Venus close-in, Mars farther out) had no
  influence on which class actually generated where. `program_constants.
  PLANET_CLASSES` gains a new per-class `zone_position_mode` (0.0-1.0,
  "how far through the zone's `[inner, outer]` AU range this class's real
  or reasoned analog sits") on every class with a single, fixed
  identity within the ecosphere zone: E/F/G/H/K/L/M/N/O/P/V. Class Q
  (eccentric orbit, extreme temperature swings) deliberately has none --
  it has no single fixed position by its own flavor. `planetPhysics.
  generate_planet_properties` reads it once a planet's class is settled
  and redraws `planet.distance` there via `utils.sample_bounded_bell`
  (the same bounded-bell-curve mechanism `size_mode` already uses for
  radius), for ordinary planets only -- explicitly skipped for moons,
  since a moon's own `distance` is its orbit around its *parent planet*,
  not an AU-scale position within the star's own zone `planet.
  habitable_zone` describes; redrawing it there would corrupt it, not
  correct it. `StarSystem.validate_system`'s existing orbital-overlap
  correction absorbs whatever reordering a class-biased redraw causes
  against already-placed neighbors, the same way it already absorbed
  `calculate_distance_for_class`'s explicit-slot distance nudging.
  `planetgen.generation.plausibility._extract_record` now also reports
  `distance`, letting a new `test_planets.py` regression suite verify
  the bias directly (per-class mean zone-fraction close to its declared
  mode; Class N < Class M < Class K in mean orbital distance; a moon's
  distance is provably untouched).

### Changed
- **Class K and Class N retuned now that they're placed at a real,
  class-appropriate distance instead of sharing Class M's midpoint.**
  Both carried an explicit "tuned to compensate for the wrong distance"
  comment (see above) -- with `zone_position_mode` now doing the
  distance part of the work for real, their climate ranges needed
  re-deriving via `climate_tuning_cli.py` rather than staying tuned
  against the old, distance-blind placement:
  - **K (Mars analog)**: `albedo_range` raised slightly (0.34-0.42 ->
    0.36-0.44) and `atm_density_range` widened/raised (0.012-0.025 ->
    0.022-0.042) to fit the new, colder starting point. Mean
    surface_temperature ~214K (real Mars ~210K, +1.9%, was +9.9% before
    this pass) and mean atmospheric_pressure ~611Pa (real Mars ~610Pa,
    +0.2%, was -11.6%) over a 1000-sample run -- the "as close as
    achievable without a zone change" caveat the old K tuning note
    carried is resolved by the zone change.
  - **N (Venus analog)**: `greenhouse_multiplier_range` cut roughly in
    half (370-420 -> 260-295) and `atm_density_range` raised (270-320 ->
    300-350) -- the old range was deliberately oversized specifically to
    compensate for N's too-cold midpoint placement, so keeping it at the
    same size now overshoots real Venus's temperature once N is
    correctly placed near the zone's hot inner edge. Mean
    surface_temperature ~737K (real Venus 737K, +0.0%) and mean
    atmospheric_pressure ~9.17MPa (real Venus ~9.2MPa, -0.3%, was -16.8%)
    over a 1000-sample run.
  - Every other affected class (E/F/G/H/L/O/P/V) keeps its existing
    albedo/atmosphere/greenhouse ranges unchanged -- none of them chase a
    single real-world numeric target the way K/N do, and each one's
    existing relative ordering (hotter/colder than its neighbors in the
    E->F->G->M/O progression, K/L, P) still holds with the new
    distance-aware placement. Their "Verified via climate_tuning_cli.py"
    comments are refreshed with new measured means reflecting the new
    placement.

## [5.20.0] - 2026-09-11

### Changed
- **`update.sh` skips the package reinstall when there's nothing new to
  install.** Previously it unconditionally re-ran the whole of
  `install.sh` (`pip install --force-reinstall`, an NLTK re-fetch,
  re-enabling Apache modules, a full permissions pass) every single
  invocation, even when `git pull` found no new commits at all -- pure
  wasted work for a caller like `examples/maintenance/planetgen-update.timer`
  that may run this monthly for years between real updates. Now, when the
  pull is a no-op, `update.sh` runs `src/migrateDb.py` directly instead
  (a cheap, idempotent no-op once the schema is already current) and
  skips the rest; `install.sh` (migration included, as its own step 2)
  still runs in full whenever the pull actually brought new commits, same
  as before.

### Added
- **`examples/maintenance/`: `update.sh` runs on the same schedule as the
  orbit update.** New `planetgen-update.service`/`.timer` run
  `sudo ./update.sh` (git pull + `install.sh`) monthly, at a fixed time
  30 minutes ahead of `planetgen-orbits@.timer`'s own now-fixed time (both
  timers dropped `RandomizedDelaySec` in favor of this deliberate,
  guaranteed ordering) -- update.sh can `pip install --force-reinstall` a
  new version of the very `stellarObjects` code `updateOrbits.py` imports,
  so the code update needs to land first, not run independently sometime
  in the same month. `planetgen-orbits@.service` also gained an
  `After=planetgen-update.service` ordering line for the case where both
  happen to be queued together. `install-maintenance-timer.sh` installs
  and enables both by default, sharing the same `/etc/planetgen/
  maintenance.env` credentials file (`update.sh`'s `install.sh` step needs
  DB credentials for `src/migrateDb.py` too); pass `--skip-update-timer`
  to opt out of unattended code updates and keep only the orbit timer, if
  this deployment's branch should only ever be updated by a human running
  `update.sh` deliberately.

## [5.18.0] - 2026-09-11

### Added
- **`examples/maintenance/`: systemd timer for `updateOrbits.py`.** An
  Ubuntu/Debian-native alternative to the raw crontab line
  `docs/database-schema.md` already documented for running the periodic
  orbital-motion update ("once a month or so"). `planetgen-orbits@.service`/
  `.timer` are a systemd *template* unit -- the instance name (e.g.
  `planetgen-orbits@planetgen.timer`) selects which database gets
  updated, so a deployment with more than one `PLANETGEN_MYSQL_DATABASE_PREFIX`
  schema enables one timer instance per database rather than needing a
  separate script per database. `install-maintenance-timer.sh` installs
  both units, writes `/etc/planetgen/maintenance.env` (mode 600) from
  `maintenance.env.example` for the shared read-write MySQL credentials
  (skipped if that file already exists, so it never clobbers credentials
  already set up), and enables the timer for each database name given on
  its command line (defaulting to `$PLANETGEN_MYSQL_DATABASE`, or
  "planetgen"). Output is captured by journald automatically, so there's
  no logfile/logrotate entry to maintain the way the crontab example
  needs.

## [5.17.0] - 2026-09-11

### Added
- **Binary mutual-orbit position.** `updateOrbits.py`/`_db.advance_orbital_phases`
  now recomputes `star_systems.binary_mutual_position_x_km`/`_y_km`/`_z_km`
  (the secondary star's Cartesian position relative to the primary) every
  time it advances `binary_mutual_orbital_phase_deg` -- the same "position
  has no independent update of its own, it just has to move whenever
  phase does" treatment `position_x/y/z_km` already gets for planets/
  moons. Derived via `utils.orbital_position_au` from
  `binary_separation_km` and the mutual orbit's own (non-near-ecliptic,
  full `[0, 180)`/`[0, 360)`-range) inclination/ascending-node/phase --
  the pair's true "direction of orbit" was already fully captured by
  those v13 orbital elements; this just keeps the derived position
  correctly in sync with them as time passes, rather than only at
  generation. Schema v14, with a migration that backfills real position
  values for existing binary rows from their own already-stored
  separation/orbital-element columns.

## [5.16.0] - 2026-09-11

### Added
- **Star motion: galactic orbit phase, plus binary mutual orbit.** Stars
  now get the same floating-point update guard planets/moons already
  have, and binary pairs get a real orbit around each other, both
  actively advanced over time.
  - Every star gains `galactic_orbital_phase_deg` (its current angular
    position around the galactic center, randomly rolled at generation --
    the same role `orbital_phase_deg` plays for a planet/moon) and
    `galactic_min_update_interval_years` (the same guard formula,
    `utils.minimum_update_interval_years` applied to
    `galactic_orbital_period_gy * 1e9` years). `_db.advance_orbital_phases`
    now advances this phase too, guarded by its own interval, reversing
    v12's "stars have no periodic update mechanism" scoping note -- they
    do now. Both stars of a binary pair, and `star_systems`' own mirrored
    `binary_galactic_orbital_phase_deg`/`binary_galactic_min_update_interval_years`,
    always carry the identical value: a binary's AU-scale separation is
    negligible next to its light-year-scale galactic orbit, so the pair
    moves around the galaxy together, not independently
    (`StarSystem.__init__` rolls the phase once and threads it to both
    stars and the proxy).
  - Binary pairs also get their own **mutual orbit** -- the two stars
    circling their common barycenter, entirely separate from (and vastly
    faster than) the galactic orbit above:
    `star_systems.binary_mutual_orbital_period_years`/`_speed_kms`
    (Kepler's third law / circular-orbit speed --
    `planetPhysics.calculate_orbital_period_years`/
    `utils.circular_orbital_speed_kms`, the same formulas a planet's own
    orbit already uses, applied to the pair's separation/combined mass),
    `_inclination_deg`/`_ascending_node_deg`/`_phase_deg` (the same
    `utils.orbital_position_au` orbital-element convention planets/moons
    use for "direction", but drawn from the full `[0, 180)`/`[0, 360)`
    range -- a binary's mutual orbital plane has no protoplanetary-disk
    reason to prefer any alignment), and its own
    `binary_mutual_min_update_interval_years` guard. `advance_orbital_phases`
    advances `binary_mutual_orbital_phase_deg` the same way.
  - Schema v13, with a migration that backfills every real-derivable
    value (both `*_min_update_interval_years` guards, plus the mutual
    orbit's period/speed) from existing rows' own already-stored data;
    the phase/orientation columns with no derivable "correct" value get
    an arbitrary `0` placeholder, the same treatment v9 already gives
    pre-existing planets'/moons' orbital orientation.

## [5.15.0] - 2026-09-10

### Added
- **Floating-point update guard.** Every planet and moon now stores
  `min_update_interval_years` -- the shortest `elapsed_years` worth
  calling `_db.advance_orbital_phases` for. Below this threshold, the
  phase delta `elapsed_years` would add is smaller than
  `orbital_phase_deg`'s own IEEE 754 double-precision resolution, so
  `MOD(orbital_phase_deg + delta, 360)` is guaranteed to round right back
  to the exact value already stored -- a wasted write that changes
  nothing. Derived purely from `period_years`
  (`utils.minimum_update_interval_years`: `period_years *
  math.ulp(360.0) / 360`, the coarsest representable step anywhere in
  `orbital_phase_deg`'s `[0, 360)` range). Not a narrative/display stat --
  purely a guard value: `advance_orbital_phases` now skips a row's
  `UPDATE` entirely (not just a no-op write, no attempt at all) when a
  call's `elapsed_years` is below it. In practice this floor sits many
  orders of magnitude below any realistic elapsed time
  (`updateOrbits.py` runs "once a month or so"), so it exists for
  correctness against a caller advancing time in much smaller steps, not
  because today's actual usage pattern comes close to tripping it.
  Scoped to `planets`/`moons` only -- `stars`' galactic-orbit values are
  fixed forever at generation time, with no periodic update mechanism to
  guard. Persisted as `planets`/`moons.min_update_interval_years` --
  schema v12, with a migration that backfills real derived values (not a
  placeholder) for existing rows.

## [5.14.0] - 2026-09-10

### Added
- **Planet/moon position.** Every generated planet and moon now gets an
  actual 3D Cartesian position (`Planet.position_x/y/z`, in AU), relative
  to its orbital anchor -- the star (or, for a binary system, the
  `BinaryStarProxy` standing in for the system's combined center) for a
  planet, the parent planet for a moon -- continuing the same "each body
  positioned relative to its immediate primary" hierarchy
  `docs/design/galaxy-coordinate-system.md` already uses one level up for
  sectors/systems relative to the galactic center. Derived from the
  existing orbital-motion elements (`orbital_inclination_deg`/
  `orbital_ascending_node_deg`/`orbital_phase_deg`, schema v9) via the
  standard circular-orbit-to-Cartesian transform
  (`utils.orbital_position_au`). Also adds `Planet.orbital_speed_kms`, a
  circular orbit's constant tangential speed (`utils.
  circular_orbital_speed_kms`, `v = 2*pi*r/T`) -- exact here, unlike the
  galactic orbit's rotation-curve model, since a planet's/moon's period is
  already known exactly from Kepler's third law. Surfaced as a new
  "Speed" row in every planet's/moon's rendered data table. Persisted as
  `planets`/`moons`.`position_x_km`/`_y_km`/`_z_km`/`orbital_speed_kms` --
  schema v11, with a migration that backfills real derived values (not a
  placeholder) for existing rows.
- `updateOrbits.py`/`_db.advance_orbital_phases` now recomputes position
  in lockstep with `orbital_phase_deg` as real time passes (a single
  set-based SQL `UPDATE` per table, matching phase advancement's own
  performance characteristics); `StarSystem.validate_system` recomputes
  position/speed too whenever it corrects a planet's `distance` post-hoc
  to resolve an orbital overlap, closing the same kind of staleness gap
  its own docstring already flagged for `period` before this.

## [5.13.0] - 2026-09-10

### Added
- **Galactic orbit.** Every generated star system now gets a circular
  orbital speed and period around the galactic center
  (`Star.galactic_orbital_speed_kms`/`galactic_orbital_period_gy`,
  `utils.calculate_galactic_orbit`), derived from the system's actual
  distance from the galactic center where known (a sector-placed system)
  or the same fixed `physical_constants.GALACTIC_CENTER_DISTANCE_LY`
  fallback `system_perimeter`/`heliosphere_radius` already use otherwise.
  Uses a simple rotation-curve model
  (`GALACTIC_ROTATION_FLAT_VELOCITY_KMS`/`GALACTIC_ROTATION_CORE_RADIUS_PC`)
  rather than a Keplerian point-mass orbit around the galaxy's total mass,
  which would overshoot Sol's real ~220-240 km/s orbital speed by roughly
  4x -- calibrated against Sol's own distance (~206 km/s, ~236 million
  years per orbit, both close to the real Sun's measured values). Surfaced
  as a new "Galactic Orbit" row in every star's/binary pair's rendered data
  table, and persisted as `stars.galactic_orbital_speed_kms`/
  `galactic_orbital_period_gy` (plus the `star_systems.binary_galactic_orbital_*`
  equivalents for binaries) -- schema v10, see `docs/database-schema.md`.

## [5.12.0] - 2026-09-10

_Originally developed and released as `5.9.0` on a separate branch; renumbered
on merge since `5.9.0` was independently already used below for the same-day
Class W removal. No functional difference from the original release._

### Changed
- **`../src/html/` is now a thin frontend for the Flask API instead of a
  direct MySQL client.** Every CGI page (`index.py`, `browse.py`,
  `sector.py`, `system.py`, `nav.py`, `search.py`, `galaxy.py`) fetches
  its data from `GET /api/...` via a new stdlib-only HTTP client
  (`html/lib/apiclient.py`, `PLANETGEN_API_BASE_URL`) rather than opening
  its own read-only MySQL connection -- the database-querying logic those
  pages used to duplicate now lives once in `queryDb.py`, shared with the
  API. `html/lib/dbutil.py` is gone; its non-DB formatting helpers
  (`esc`/`linkify_location`/`format_density`) moved to a new
  `html/lib/fmt.py`. The API itself gained what this needed: `?db=` on
  every route (multi-schema, matching the browser's own picker --
  `stellarObjects._db.list_databases`/`resolve_database`, moved out of
  `html/lib/dbutil.py` into the package itself), `GET /api/databases`,
  `GET /api/galaxy/sectors`, `GET /api/search` (the full faceted-search
  query layer, ported from `html/search.py`), `sector_id=none` on
  `GET /api/systems` (standalone systems), and richer `GET /api/sectors`/
  `GET /api/sectors/<id>`/`GET /api/systems/<id>` responses (display-
  ready fields -- ids, quadrant/location, star summaries, galaxy
  placement -- distinct from `stellarObjects._db.load_sector`/
  `load_star_system`'s *generation* object graph, still reachable the
  same way). Apache's example vhost
  (`examples/apache/planetgen.conf.example`) mounts the API at `/api/`
  (`WSGIScriptAlias`, in its own `WSGIDaemonProcess`) on the same vhost
  that serves `../src/html/` -- landed alongside `html/api/`'s own move
  into the same `DocumentRoot` (`[5.10.0]`), resolving the open question
  `TODO.md` left from the 5.3.2/5.3.3 file-system cleanup from two
  independent, orthogonal angles at once (where the API's source lives,
  and how `html/` gets its data) that happen to combine cleanly. See
  `docs/api.md` and `docs/html-interface.md`.
- System Map: a body with `life_chemical` set (habitable) now gets a
  small green badge on its marker, surfaced in the info panel's "Life
  Chemistry" field too (`html/lib/systemmap.py`'s `has_life`,
  `static/systemmap.js`) -- the one piece of state the existing
  circle+letter marker couldn't show at a glance. `docs/TODO.md`'s long-
  stale "sprite-based graphical system view" item is resolved: the System
  Map (an interactive scaled-orbit SVG diagram, not bitmap sprite art)
  already satisfied it, per that module's own docstring.

## [5.11.1] - 2026-09-10

### Changed
- **Moon tidal locking is now physics-based, not a flat probability.**
  [5.11.0] gave a moon a flat 75% (`MOON_TIDAL_LOCK_PROBABILITY`) chance
  of being tidally locked; replaced with an actual tidal-despinning
  timescale estimate (new `planetPhysics._tidal_locking_timescale_seconds`,
  the standard simplified Murray & Dermott formula: `t_lock = (2Q/15k2) *
  omega0 * a^6 * m_moon / (G * M_primary^2 * R_moon^3)`, with fixed
  representative `Q`/`k2` values for a rocky/icy body -- new
  `physical_constants.MOON_TIDAL_DISSIPATION_Q`/`_LOVE_NUMBER_K2`).
  Verified directly against three real systems spanning many orders of
  magnitude: ~47 million years for the real Earth-Moon system (real
  estimates: tens of millions of years), ~90 years for Mars/Deimos
  (consistent with Deimos being locked given its tiny size), and ~1
  billion years for Saturn/Iapetus (consistent with real estimates of a
  billion-year-plus despinning time for that unusually slow case). A
  candidate (pre-locking) rotation period is drawn first, same as any
  planet; a moon actually ends up locked only if that timescale is
  shorter than the system's age (`star.age`, the best available proxy --
  planets/moons don't carry an independent age of their own).
- **Moon orbital distance is now drawn log-uniformly, not
  linearly-uniformly.** Surfaced by the physics-based tidal locking
  above: `generate_moons`' distance range can span many orders of
  magnitude (its outer bound reaches 1/5 of the parent planet's own Hill
  radius -- tens to hundreds of millions of km for a large planet, far
  beyond where any real large moon actually orbits, e.g. our Moon at
  384,400 km), and a plain `random.uniform` over that range spends almost
  all its density in the single largest order of magnitude -- nearly
  every generated moon landed implausibly far out, so almost none had
  time to tidally lock (measured: ~1.3% of moons locked). Real moon
  systems are much closer to log-spaced (e.g. the Galilean moons run
  421,700 / 671,100 / 1,070,400 / 1,882,700 km, each roughly 1.5-1.6x the
  last), which log-uniform sampling matches far better while still
  allowing occasional genuinely distant/irregular moons. Measured effect:
  ~1.3% -> ~14% of moons locked over the same sample size -- still a
  minority overall (this generator's moons span a much wider population
  than just the handful of large, close, well-known real moons that
  dominate popular intuition about "most moons are locked"), but an order
  of magnitude more of them landing close enough to plausibly have
  locked.

## [5.11.0] - 2026-09-10

### Added
- **Orbital and rotational motion.** Every planet/moon now has a real 3D
  orbital orientation and a rotation period, plus a live position that a
  new, separately-run script advances over time -- resolves
  `docs/TODO.md`'s "Introducing realistic orbital paths and speeds..."
  entry.
  - `Planet` gains four new attributes (`planetPhysics.
    generate_orbital_motion_properties`): `orbital_inclination_deg`/
    `orbital_ascending_node_deg` (fixed at generation time -- together
    they orient this generator's circular-orbit model in 3D; planets draw
    from a tighter, real-solar-system-like range than moons, which can be
    tilted much further), `orbital_phase_deg` (this body's current
    position angle around its orbit -- the one field that changes after
    generation), and `rotation_period_hours` (axial "day length" -- a
    static descriptive stat; no rotational phase is tracked, by design,
    since nothing consumes "which side currently faces the primary"). A
    moon has a real chance (75%, `MOON_TIDAL_LOCK_PROBABILITY`) of being
    tidally locked (rotation period equal to orbital period) -- the norm
    for real large moons, not a rare special case; otherwise rotation
    period is drawn from a `body_type`-appropriate range (`physical_constants.
    ROTATION_PERIOD_RANGE_HOURS`).
  - **New `src/updateOrbits.py`.** Advances every planet's/moon's
    `orbital_phase_deg` in the configured database based on real elapsed
    time since the last run (`stellarObjects._db.advance_orbital_phases`,
    a single set-based `UPDATE` per table, not a per-row Python loop) --
    meant to be run periodically (e.g. cron, "once a month or so"), not
    on every generation run. A new `orbit_simulation_state` singleton row
    tracks when it last ran, measured server-side via `TIMESTAMPDIFF`
    rather than trusting the calling process' own clock to agree with the
    database server's (`get_orbit_update_elapsed_years`). The first run
    against a database just establishes that reference point (nothing to
    advance yet) rather than guessing a start time.
  - **Orbital period fix, needed for this feature to mean anything.**
    `Planet.period` used to be `sqrt(distance_au^3)` unconditionally --
    correct Kepler's-third-law shorthand only for a 1-solar-mass primary.
    For an ordinary planet orbiting a star of very different mass, and
    *especially* for a moon (whose real primary is its parent planet, not
    the grandparent star `self.star` still points at for other purposes
    like life chemistry) this was wrong by orders of magnitude -- a moon's
    period came out as if it orbited its host star directly at that same
    tiny distance. New `planetPhysics.calculate_orbital_period_years(distance_au,
    primary_mass_kg)` takes the actual primary's mass
    (`Planet.__init__`'s new `primary_mass_kg` parameter, threaded in
    from `generate_moons` as the parent planet's own `mass` for a moon,
    defaulting to `star.mass` for an ordinary planet -- transparently
    correct for a binary system too via `BinaryStarProxy.mass`'s existing
    effective-mass property). Also fixed a related latent staleness bug
    found while testing this: `StarSystem.validate_system` adjusts a
    planet's `distance` after generation to resolve orbital overlap and
    already recomputed atmospheric conditions to match, but never
    recomputed `period` -- now it does.
  - Schema v8 -> v9: `orbital_inclination_deg`/`orbital_ascending_node_deg`/
    `orbital_phase_deg`/`rotation_period_hours` on `planets`/`moons`, and
    the new `orbit_simulation_state` table. `stellarObjects._db.
    _migrate_v8_to_v9` is the first real per-version migration step of the
    MySQL era (every database before this one started fresh, already at
    the then-current schema) -- existing rows default to `0`/`24` (an
    arbitrary but harmless placeholder), since this generator never ran
    its actual random generation for them; every body generated from this
    point on gets real values. See `docs/database-schema.md`.

## [5.10.1] - 2026-09-10

### Fixed
- **Subdwarf (Yerkes VI) age modeling could produce a pre-Big-Bang star.**
  `Star._calculate_initial_star_age_and_lifespan` routed every non-main-
  sequence, non-white-dwarf Yerkes class (giants, subgiants, bright
  giants, supergiants, hypergiants, *and* subdwarfs) through the same
  model: derive main-sequence lifespan from the star's own already-
  generated mass, then draw age from the post-main-sequence window. That
  model's mass-sampling side already got a reject-and-resample guard
  against pre-Big-Bang progenitor masses ([5.3.1]) for every class it
  applies to -- except Yerkes VI, deliberately excluded because its
  *entire* allowed mass range (0.1-0.8 Msun) implies a main-sequence
  lifespan longer than the age of the universe, so every draw would have
  been rejected. That left subdwarfs still running through the age *model*
  itself, unguarded, which routinely produced ages of hundreds of billions
  of years -- because real subdwarfs (sdB/sdO) are thought to form via
  binary mass-stripping near the tip of a lower/intermediate-mass
  progenitor's red-giant branch, not single-star post-main-sequence
  evolution, so this star's own post-strip mass was never a valid stand-in
  for a progenitor's main-sequence lifespan in the first place -- the
  model was wrong for this class, not just missing a guard rail (flagged
  as an open design question in `docs/TODO.md`'s "Future ideas" ever
  since).
  - Yerkes VI now has its own age-generation branch, entirely independent
    of the star's own mass: age is drawn directly from a dedicated,
    old-population-biased range (new `program_constants.SUBDWARF_MIN_AGE_GY`
    /`SUBDWARF_MAX_AGE_GY`, 1.0-13.5 Gy, with the same young/old
    `SystemConfig.AGE` bias every other branch applies), reflecting that
    real subdwarf progenitors are typically old, low-mass population stars
    (a short-lived, higher-mass star wouldn't have had time to reach the
    RGB tip and get stripped). Lifespan is that age plus a short
    remaining-phase window (new `SUBDWARF_REMAINING_PHASE_MIN_GY`/
    `_MAX_GY`, 0.05-0.3 Gy) rather than a separately-derived value, since
    the real core-helium-burning subdwarf phase is short relative to the
    age itself.
  - New regression test `test_subdwarf_age_never_exceeds_universe_age`
    (`test_star_matrix.py`) locks this in across every spectral letter.

## [5.10.0] - 2026-09-10

### Changed
- **Moved the Flask API (`src/api/`) into `src/html/`.** Resolves the
  open question `docs/TODO.md` had carried since the 5.3.2/5.3.3
  file-system cleanup: the API now lives at `src/html/api/`, served from
  the same checkout/deployment tree as the interim CGI browser instead of
  a second, separately-tracked location. `src/wsgi.py` moved alongside it
  to `src/html/wsgi.py` (same "lives next to `api/`, so `import api` just
  works via `sys.path[0]`" property as before, plus an explicit
  `sys.path` entry for `src/` now that `queryDb`/`stellarObjects` are a
  directory farther away). `pytest.ini`'s `pythonpath` gained `src/html`
  so `test_api.py`'s `from api...` imports keep resolving unchanged.
  `examples/apache/planetgen.conf.example` gained a `WSGIScriptAlias /api`
  pointing at `html/wsgi.py`, plus deny-all `<Directory>`/`<Files>` blocks
  for `html/api/` and `html/wsgi.py` itself (mirroring the existing
  `html/lib/` block -- both hold source that's imported, never meant to
  be requested directly). Along the way, corrected a pre-existing
  inaccuracy in `docs/api.md`'s Apache deployment guidance: the vhost's
  `SetEnv` directives configure the CGI browser (`mod_cgi`/`mod_cgid`
  copies them into each script's real process environment) but never
  reach `os.environ` under `mod_wsgi` -- `PLANETGEN_MYSQL_*` for the API
  needs to come from the Apache service's own process environment instead
  (e.g. `/etc/apache2/envvars`).

## [5.9.1] - 2026-09-10

### Changed
- **Class P climate tuning.** Gave Class P the same
  `atm_molar_density_range`/`atm_density_range`/`greenhouse_multiplier_range`
  treatment the M/O/H/K/L/N/E/F/G/V pass ([5.3.7]) gave the other habitable
  classes; P had stopped at just its own `albedo_range` that pass. No
  single real-world analog for P, so this isn't chasing a target delta the
  way M/K/N are -- instead the new ranges make the class's own "cold,
  glaciated"/"thinning with age" flavor text physically real: molar
  density stays near Earth's real value (P's atmosphere text names
  oxygen/nitrogen/argon, not a heavier CO2-like mix), `atm_density` is set
  thin, and `greenhouse_multiplier` is set weak so the cold comes from
  genuine physics on top of the class's already-tuned high albedo, not
  albedo alone. Verified via `climate_tuning_cli.py --class P`: mean
  surface_temperature ~219K (well below freezing, clearly colder than
  Class M's ~286K), mean atmospheric_pressure ~10.4kPa (~0.1 atm) over a
  400-sample run across the full host-star grid.

## [5.9.0] - 2026-09-10

### Removed
- **Class W ("a tidally locked world with extreme temperature
  variations") removed entirely.** Its day/night-split identity can't be
  produced from a single global `surface_temperature` scalar under this
  generator's climate model -- no per-class range (albedo, molar density,
  greenhouse multiplier, or atmosphere density) reaches it without an
  actual dayside/nightside model this generator doesn't have, flagged as
  a known gap during the [5.3.7] climate-tuning pass and left open in
  `docs/TODO.md`'s "Investigate Further" section ever since. Rather than
  build a whole day/night thermal model for one class, cut it entirely --
  removed from `PLANET_CLASSES`, `PLANET_CLASS_PROBABILITIES` (its
  0.0001 weight just dropped; these are relative weights, not a
  normalized distribution, so nothing needed redistributing),
  `HABITABLE_PLANET_CLASSES`, and `MOON_BLACKLIST`
  (`program_constants.py`), plus its color entry in
  `html/lib/systemmap.py`'s `_CLASS_COLORS`.

## [5.8.3] - 2026-09-10

### Changed
- **Provisional sector designation is now a single packed hex number.**
  `provisional_sector_designation` (added in [5.8.2]) dropped its
  `"R<ring>-Q<quadrant>-<slot>"` letter/dash format for a plain
  bit-packed hex integer -- `ring` in the high bits, `quadrant - 1` in
  the next 2, the raw `shell_slot_index` in the low 32
  (`DESIGNATION_SLOT_BITS`/`DESIGNATION_QUADRANT_BITS`), e.g.
  `"1500002EE0"` (was `"R5-Q2-2EE0"`). Genuinely reversible back to
  `(ring, quadrant, shell_slot_index)` now too -- the fixed-width bit
  fields have no ambiguous boundary the way concatenating separately-
  sized hex numbers would.

## [5.8.2] - 2026-09-10

### Added
- **Provisional sector designations for un-generated/unvisited addresses.**
  `planetgen.galaxy.geometry.provisional_sector_designation(shell_index,
  shell_slot_index, edge_pc, edge_ly)` builds a short, human-readable label
  for a `(shell_index, shell_slot_index)` sector address --
  `R<ring>-Q<quadrant>-<slot>`, Ring and slot index in uppercase hex, Quadrant
  a plain 1-4 digit (e.g. `"R5-Q2-2EE0"`) -- the way a real astronomical
  catalog gives a not-yet-fully-characterized object a provisional name from
  its position rather than waiting for one. `sector_ring`/`sector_quadrant`
  (new, duplicating `html/lib/galaxymap.py`'s identically-named Ring/Quadrant
  concept in the dependency-free `galaxyGeometry` module so `galaxyGen.py`'s
  CLI doesn't need to import the web front-end layer) back it, and are both
  O(1) -- no scan over a shell's other slots, consistent with this codebase's
  galaxy-skeleton design principle of never doing per-sector work
  proportional to a shell's slot count. `galaxyGen.py`'s `--shell`/
  `--center-sector` generation modes now print this designation alongside
  each newly saved sector's procedural name and raw shell/slot numbers.

## [5.8.1] - 2026-09-10

### Added
- **NAV Map: a rendered plot of the NAV feature's origin/destination/route.**
  `queryDb.nav_between` now also returns `origin_position`/
  `destination_position` (the light-year positions its `direct` course was
  computed from) and, when a route was found, each hop's own position
  (`route["positions"]`) -- surfaced from `GET /api/nav`'s JSON response
  too (see `docs/api.md`'s updated "NAV" section). `src/html/lib/navmap.py`
  builds a new "NAV Map" panel on `html/nav.py`'s results page from that
  data: a flat, static, top-down SVG plot of the galactic X-Y plane --
  origin and destination as labeled, clickable points, a dashed line for
  the direct course, and a solid polyline through the optimal route's
  intermediate hops when one exists. Modeled on `galaxymap.py`'s flat 2D
  SVG rather than `starmap.py`'s rotatable 3D CSS scene: like a galaxy
  Quadrant, it's deliberately blind to altitude (the course panel's own
  Altitude figure already covers that axis), auto-scaled to whatever
  points it's given (no fixed sector size to normalize against) with one
  uniform light-years-per-pixel ratio on both axes so azimuth angles
  aren't visually distorted, plus a `+X` compass tick tying the plot's
  orientation to the course panel's own azimuth convention and a
  light-year scale-bar legend. Closes the rendered-image gap the original
  NAV feature (CHANGELOG [5.8.0]) left open in `docs/TODO.md` and
  `src/api/routes.py`'s `systems_near` TODO comment.

## [5.8.0] - 2026-09-10

### Added
- **NAV feature: course, distance, and optimal routing between two
  systems.** Three new pure/query modules plus an API endpoint and a web
  page:
  - `planetgen/galaxy/navigation.py` -- `course_between` (Euclidean
    distance plus galactic-plane-relative azimuth/altitude: azimuth from
    +X in the galactic X-Y plane, altitude as elevation above/below that
    plane) and `warp_travel_times` (`velocity_multiple_of_c = warp_factor
    ** (10/3)`, reported at warp 1/3/6/9, formatted via the existing
    `utils.years_to_time_string`).
  - `planetgen/galaxy/nav_graph.py` -- `build_knn_adjacency` (a symmetrized
    k-nearest-neighbor adjacency graph over a `{id: (x, y, z)}` position
    set) and `shortest_path` (Dijkstra) for the "optimal route via
    adjacent systems" half of NAV.
  - `queryDb.nav_between` -- resolves NAV availability between two
    systems (unavailable if either has no sector; same-sector always
    available; cross-sector only when both sectors have a galaxy
    placement), combining a sector's galaxy-frame center (parsecs) with a
    system's sector-local offset (milliparsecs) into one absolute
    position (no such combinator existed before this), then returns the
    direct course plus the optimal route.
  - `GET /api/nav?from=<id>&to=<id>` -- JSON endpoint over `nav_between`,
    `404` for an unknown system id, `400` (`NavUnavailable`) when NAV
    doesn't apply to the pair. See `docs/api.md`'s new "NAV" section.
  - `src/html/nav.py` -- a destination picker (a `<select>` of the
    origin's own sector-mates, plus a typed destination-id field when
    cross-sector NAV is available) and a results page (direct course,
    warp travel times, and the hop-by-hop optimal route, each hop linking
    to `system.py`). `system.py` links here ("Navigate from here")
    whenever a system has a sector.

## [5.7.0] - 2026-09-09

### Changed
- **Merged the parallel Phase 4 "galaxy density skeleton" work (schema
  v6-v8: `sector_vertices`, `galaxy_shape`, `galaxy_shell_band`,
  `galaxyGen.ensure_sector_generated` lazy generation -- see [5.4.7] and
  [5.4.8] below) into the MySQL backend from [5.5.0].** That work landed
  on `main` entirely against the pre-port SQLite backend while this
  branch's MySQL port was in flight from the same base commit, so every
  new table, `_db.py` function (`get_sector_id_at`, `save_galaxy_shape`,
  `get_galaxy_shape`, `replace_galaxy_shell_bands`,
  `get_galaxy_shell_bands`, the vertex-writing/-reading in
  `insert_sector`/`get_sector_galaxy_position`), and CLI script
  (`galaxyGen.py`'s rewrite, the new `galaxyPlan.py`) needed porting to
  MySQL/`pymysql`/`MySQLConfig` as part of reconciling the two lines of
  work, same conventions as the rest of the MySQL port: `?` placeholders
  (translated by the `Connection` wrapper), `--mysql-*`/`config=` instead
  of `--db-path`/`db_path=`, and MySQL's `INSERT ... ON DUPLICATE KEY
  UPDATE` instead of SQLite's `INSERT ... ON CONFLICT DO UPDATE` for
  `save_galaxy_shape`'s singleton-row upsert.

### Fixed
- **`planetgen/galaxy/geometry.py` called an undefined
  `_theta_for_index` in both `sector_position_pc` and
  `sector_wedge_vertices_pc`** -- a pre-existing bug on `main` (confirmed
  present there independent of this merge), only `_phi_for_index` had
  ever been defined despite the module's own `GOLDEN_RATIO` docstring
  describing the golden-angle azimuthal step it was supposed to compute.
  Every galaxy-placement test failed with `NameError` until this merge's
  full test run surfaced it. Added the missing function
  (`theta_i = (2*pi*i / GOLDEN_RATIO) mod 2*pi`), matching the exact
  formula `test_galaxy_geometry.py`'s own worked-example test already
  documented and asserted against.
- **`stellarObjects._db.Connection`'s `execute`/`executemany` could
  return a `tuple` instead of a `list` from `.fetchall()` when zero rows
  matched** -- confirmed by testing, `pymysql`'s own cursor returns `()`
  for no rows but a `list` when rows exist, unlike `sqlite3`'s cursor
  (always a `list`). Surfaced by a real MariaDB test run as a spurious
  `assert result == []` failure (`test_ensure_sector_generated_reports_no_content_outside_every_stored_band`)
  that never appeared without a live server. Added a `_Cursor` proxy
  wrapping every returned cursor to normalize `fetchall()` to always
  return a `list`, so no other call site (present or future) can trip
  over the same inconsistency.

## [5.6.0] - 2026-09-09

### Added
- **Rate limiting (Flask-Limiter), applied to the whole API out of the
  box.** Default limits (200/day, 50/hour per client IP -- Flask-Limiter's
  own quickstart example, overridable via `PLANETGEN_RATELIMIT_DEFAULT`)
  apply to every route except `/api/health`; every write endpoint (below)
  additionally layers a stricter 10/minute limit on top
  (`routes.WRITE_RATE_LIMIT`). Exceeding a limit returns `429
  {"error": "rate limit exceeded", "detail": "..."}` (never Flask-Limiter's
  own default plain-text body) with `Retry-After`/`X-RateLimit-*`
  headers. Storage backend defaults to in-memory
  (`PLANETGEN_RATELIMIT_STORAGE_URI`, correct for a single-process
  deployment; a multi-worker `mod_wsgi`/`gunicorn` deployment needs a
  shared backend, e.g. Redis, or each worker under-enforces the
  configured limit by tracking its own separate counters).
- **Write-stub endpoints**: `POST`/`PATCH`/`DELETE` on `/api/sectors` and
  `/api/systems`. Every one validates its JSON request body (sectors get
  a fully mapped-out `{"name", "edge_ly"}` schema; systems only check
  "is a JSON object" for now -- see docs/api.md for why that one's
  field-level schema is still an open design question) and applies
  `WRITE_RATE_LIMIT`, but always responds `501 {"error": "... is not
  implemented yet"}` -- no row is ever inserted, updated, or deleted.
  Exists now so the request/response contract is settled and testable
  before the real database logic lands; see docs/api.md's "Write
  endpoints" section for what filling them in for real will also need
  (a write-capable database account, and authentication/authorization --
  neither exists yet, and both are load-bearing on this API staying
  read-only in practice today).

## [5.5.0] - 2026-09-09

### Added
- **Flask API: pagination, input validation, health check, and JSON-only
  error handling.** `/api/sectors`/`/api/systems` now return a paginated
  `{"items", "total", "limit", "offset"}` envelope instead of a bare list
  (`limit` defaults to 100, clamped to 500; `offset` defaults to 0) --
  this project's own roadmap plans galaxy-scale generation, so an
  unbounded listing endpoint would eventually return an unbounded
  response. Every query parameter (`limit`/`offset`/`sector_id`/`radius`)
  is now validated and rejected with a `400 {"error": "..."}` rather than
  silently ignored or crashing. Added `GET /api/health` for
  liveness/readiness monitoring. Every error response -- 400, 404
  (including an unmatched route), 405, and 500 -- is now JSON, never
  Flask's default HTML error page, and a 500 never leaks exception detail
  to the client. See `docs/api.md`.
- **MySQL backend (`TODO.md` Phase 5), replacing SQLite entirely.**
  `stellarObjects/schema.sql` is now MySQL/InnoDB DDL (`BIGINT UNSIGNED`
  ids, `DOUBLE`/`VARCHAR`/`TEXT`/`LONGTEXT` typing, a real
  `schema_migrations` tracking table replacing `PRAGMA user_version`,
  every index/foreign key declared inline per table for
  `CREATE TABLE IF NOT EXISTS` idempotency). `stellarObjects/_db.py` now
  talks to MySQL via `pymysql` (pure-Python driver) through a small
  `Connection` wrapper that keeps every existing call site's
  `conn.execute(sql, params)`/`row["column"]` shape unchanged, backed by
  a real connection pool (`DBUtils.PooledDB`) per TODO.md's "add real
  connection pooling" note. Every tool that touches the database
  (`sectorGen.py`, `systemGen.py`, `galaxyGen.py`, `queryDb.py`,
  `migrateDb.py`, `src/api/`, `src/html/`) now takes `--mysql-*` flags/
  `PLANETGEN_MYSQL_*` environment variables (`stellarObjects._db.MySQLConfig`)
  instead of `--db-path`/`PLANETGEN_DB_PATH`; the CGI browser's database
  picker now lists MySQL schemas on the configured server (filtered by
  `PLANETGEN_MYSQL_DATABASE_PREFIX`) instead of `.db` files in a
  directory. A new one-time `src/migrateSqliteToMysql.py` script imports
  an existing pre-port SQLite database (already at schema v5) into
  MySQL. The SQLite-specific `v1`-`v5` in-place migration machinery
  (`migrate_database`'s per-version steps, gzip file backups,
  `BACKUP_MARKER`) is removed entirely, since every MySQL database this
  project creates now starts at the current schema directly. See
  `docs/database-schema.md`.

### Changed
- `setup.py`'s `install_requires` gained `pymysql`/`DBUtils` (core
  dependencies now, not just the `api` extra) -- every database-touching
  entry point needs them, not only the Flask API.

### Fixed
- **`test_mdconvert.py`'s own `sys.path` setup pointed at a
  nonexistent top-level `html/lib/` instead of `src/html/lib/`,
  silently masked by `test_db_migration.py` (collected first,
  alphabetically) inserting the correct path first.** Deleting
  `test_db_migration.py` (see above) exposed it; fixed the path
  computation directly.

## [5.4.8] - 2026-09-09

### Added
- **Galaxy-wide density "skeleton"** (schema v8: `galaxy_shape`,
  `galaxy_shell_band`; `planetgen/galaxy/skeleton.py`; `galaxyPlan.py`):
  precomputes and persists where the galaxy has any content at all,
  without storing a single sector's position/density/vertices -- those are
  pure deterministic functions of `(shell_index, shell_slot_index)` plus a
  handful of galaxy-wide shape parameters, so they're recomputed on demand
  instead. `galaxy_shape` holds the galaxy's shape parameters as one
  singleton row; `galaxy_shell_band` holds one row per contiguous
  *candidate* slot-index band per shell (almost always exactly one, found
  via an exact upper bound over spiral-arm azimuth, bisected to precision)
  -- a safe superset, not a per-sector list. Reduces the galaxy's
  structural storage from the ~1 TB a naive per-sector plan table would
  need down to ~350 KB at real Milky-Way scale (~4,076 rows), built in
  parallel across shells (`multiprocessing.Pool`) in under a second. See
  `docs/design/galaxy-coordinate-system.md` section 10 for the full
  analysis and measurements.
- **`galaxyGen.ensure_sector_generated(shell_index, shell_slot_index)`**:
  the lazy, visit-triggered generation entry point built on the skeleton
  above -- returns an already-generated sector if one exists; otherwise
  checks the stored band and exact density to decide whether the address
  holds anything, and if so generates and persists it on the spot, using
  that position's own `relative_density` as the actual system-count
  multiplier (so a bulge sector and a sparse outer-disk sector generate
  proportionally different counts, not a uniform default). A new
  `sectors.UNIQUE (shell_index, shell_slot_index)` constraint (schema v8)
  turns a concurrent visit race into a recoverable `IntegrityError`
  instead of a duplicate row.
- **Galaxy disk/spiral density model implemented** (`planetgen/galaxy/density.py`):
  exponential disk radial falloff x sech^2 vertical scale-height x
  logarithmic spiral-arm modulation, plus a spherical bulge, normalized so
  `relative_density == 1.0` at a calibration point. `build_galaxy_shape`
  auto-calibrates the normalization constant; `predicted_star_count` scales
  a caller-supplied baseline expected system count by relative density.
  Tested in `tests/test_galaxy_density.py`, and exercised end-to-end
  alongside sector geometry in a full (unfilled) small spiral galaxy
  simulation.

### Fixed
- **`galaxyDensity.relative_density` could raise `OverflowError` far off
  the galactic plane** relative to `disk_scale_height_pc`: its vertical
  falloff term computed `1.0 / math.cosh(x) ** 2` directly, which raises
  once `|x|` exceeds ~710 -- found by a new skeleton test using a
  toy-scale shape, reachable in production for any shape with a small
  scale height relative to its own radius. Fixed with an algebraically
  equivalent, overflow-safe rewrite (`_sech_squared`) that underflows to
  the correct `0.0` limit instead of crashing.
- **Two same-shell Voronoi tessellation bugs found by that end-to-end
  simulation, both specific to small/sparse shells** (`sectorGeometry.py`):
  (1) the same-shell candidate projection used the raw chord vector
  instead of a proper gnomonic (central) projection, understating real
  separation worse the farther away a candidate was; (2) the half-plane
  test applied after that projection used the flat-plane bisector formula,
  which only approximates the true spherical bisector in the small-angle
  limit and was too permissive on sparse shells, silently under-clipping
  cells. Both fixes are covered by new regression tests parametrized to
  include a small shell (shell 1), plus an exhaustive
  every-vertex-is-shared check for a fully-populated small shell. See
  `docs/design/galaxy-coordinate-system.md` section 9 for the full
  derivation. Re-running the simulation after both fixes: 0 unexplained
  gaps across 31,255 generated outer vertices (previously 4,798).

## [5.4.7] - 2026-09-09

### Added
- **Every galaxy-placed sector now has explicit vertices with genuinely
  zero gaps against its same-shell neighbors, stored as plain relational
  rows.** New `sector_vertices` table (schema v7) -- one row per vertex,
  no JSON or other serialized blob anywhere in the schema (matching this
  project's existing convention: variable-length structured lists always
  get their own child table, the same treatment `asteroid_belt_composition`/
  `planet_reflection_spectrum` already have). Each sector's vertices are
  built from an exact local spherical Voronoi tessellation among its
  same-shell neighbors (`stellarObjects/sectorGeometry.local_lateral_cell`)
  extruded radially between the shell's inner and outer bounding spheres
  (`prism_vertices`). Lateral sharing is exact, not approximate: each
  shared corner is the 3D circumcenter of a sector and two of its
  neighbors -- a plain geometric fact independent of which of the three
  computes it, so two real neighbors land on identical floating-point
  values (~1e-16 agreement, verified directly) rather than two nudged
  approximations. Vertex/face count varies per sector (typically 5-7, mean
  6.0) because it has to: a cube tiling of a sphere can't be gap-free in
  general (the same reason a soccer ball needs pentagons mixed with
  hexagons), so a fixed 8-vertex shape cannot exactly reconcile a sector
  with more real neighbors than it has faces -- confirmed by an earlier,
  never-released attempt at exactly that (fixed-cube corner averaging),
  which only ever reduced gaps (~35% aggregate), not eliminated them.
  Radially, adjacent shells match in total area covered, not
  vertex-for-vertex (a "non-conforming mesh interface," the same technique
  used where independently-meshed regions meet in finite-element/CFD
  meshing) -- shell k's outer bound and shell k+1's inner bound are the
  same sphere, and each shell's own sectors fully and independently tile
  it. Getting this fast required exploiting that the underlying placement
  is a Fibonacci lattice: true same-shell neighbors concentrate at
  Fibonacci-number index offsets (confirmed against real data: exact
  offsets of F_20 through F_23), which collapses same-shell neighbor
  search from ~0.4-0.5s/sector (the naive radius-search cost near a large
  outer shell's equator, exactly where the disk/spiral density model
  concentrates real generation) down to ~0.3-0.4ms/sector -- validated
  against a guaranteed-correct brute-force search across 840 cases
  spanning the full polar range and shell sizes from 3 to 211 million
  slots with zero mismatches, with small shells
  (`shell_sector_count(k) <= 2000`) falling back to unconditionally-correct
  brute force since the underlying asymptotic theory isn't reliable that
  close to the galactic core. `galaxyGen.py` computes and stores this for
  every sector it generates; `sectorGen.py`'s standalone (non-galaxy) CLI
  is unaffected, same as every other galaxy-frame data. New
  `_migrate_v6_to_v7` schema migration drops any v6 database's briefly-lived
  `vertices_pc` JSON column entirely rather than converting it into rows
  (a sector's vertices are cheap to recompute from its address if ever
  actually needed) -- `_migrate_v4_to_v5`/`_migrate_v5_to_v6` need no
  special-casing at all, since a v4 or v5 source's `sectors` table already
  has the same shape as the current one.

## [5.4.6] - 2026-09-07

### Fixed
- **System page table of contents was simply unavailable below a 90rem
  window width.** It lived exclusively in a fixed right-hand margin rail
  (`.toc` in `src/html/static/style.css`) that only existed at
  `min-width: 90rem`; narrower windows got `display: none` and no
  alternative. It's now a collapsible pulldown rendered inline above the
  description at any width, expanding on click, and only switches to the
  fixed sidebar (always shown open) once the window is wide enough for
  one. Built with a hidden checkbox + `<label>` rather than native
  `<details>`/`<summary>`: a closed `<details>`'s children turn out to
  stay out of layout/paint via internal browser state that isn't fully
  reachable through CSS overrides of `display`/`content-visibility`
  (confirmed by testing -- only toggling the element's `.open` property
  itself, not any stylesheet rule, restored it), which made "always open
  at the wide breakpoint" unreliable. A plain checkbox has no such
  internal state to fight.

## [5.4.5] - 2026-09-07

### Fixed
- **Sector Map: the 3D scene drew over surrounding page content instead
  of staying confined to its panel.** Zooming in (or rotating to an angle
  where the cube's diagonal grew past its nominal footprint) had nothing
  bounding where the map was visible, so it grew and drew over the rest
  of the page. Added `.starmap-viewport`, a fixed-size (`aspect-ratio: 1
  / 1`, matching the square footprint the old static image occupied),
  `overflow: hidden` window that the 3D scene now rotates/zooms/pans
  inside of -- clipped at its edges instead of spilling out. The 3D
  content inside renders/rotates/occludes exactly the same either way;
  this only bounds where it's visible from the outside.

## [5.4.4] - 2026-09-07

### Added
- **Sector Map is now interactive 3D: drag to rotate, scroll (or +/-
  buttons) to zoom.** Replaced the static server-baked isometric SVG
  projection in `src/html/lib/starmap.py` with a real CSS 3D scene
  (`transform-style: preserve-3d`, orthographic -- no `perspective`), so
  the browser's own compositor handles rotation and occlusion instead of
  a hand-rolled JS matrix routine; `src/html/static/sectormap.js` tracks
  two rotation angles and a zoom factor and feeds them straight to the
  scene's CSS transform. Each star dot is billboarded (counter-rotated
  every frame to keep facing the camera) so it stays a circle instead of
  going edge-on as the view turns, and dot-click detection is resolved by
  geometry (`getBoundingClientRect`) rather than native hit-testing,
  since the latter turns out unreliable for elements nested this deep in
  a rotated `preserve-3d` hierarchy.

### Changed
- **Sector Map star colors now come from the star's actual named
  spectral color** (`SPECTRAL_CLASS_COLORS`: Blue/Blue-White/White/
  Yellow-White/Yellow/Orange/Red, keyed off `star_type`'s leading letter)
  instead of a raw Kelvin-to-RGB blackbody approximation, so a "White
  Giant" reads white and a "Blue Giant" reads blue regardless of its
  exact temperature. Luminosity shades each color's vividness/lightness
  (brighter = more vivid, dimmer = more muted) and temperature nudges
  lightness within the star's own spectral band.
- **Binary systems draw two dots** (a larger primary and a smaller
  secondary, capped at 65% of the primary's radius and offset to its
  lower-right, overlapping) built from each component's own `stars` row,
  instead of one dot from the system-level `binary_type`/
  `binary_temperature_k` summary -- so each half of a binary is sized and
  colored from its own actual data.

See `docs/TODO.md` ("Near-term: interim `../src/html/` browser
enhancements") for where this started as a plan.

## [5.4.3] - 2026-09-07

### Changed
- **Web interface: system-page table of contents moved to a fixed
  right-margin rail.** The collapsed `<details>` TOC added in [5.4.2]
  still lived inline next to the description as a flex sibling, narrowing
  the prose column whenever it was open. `system.py`'s `_toc_html` now
  renders a plain, always-expanded `<nav>`, and `style.css` positions
  `.toc` fixed in the right margin (mirroring the left `.sidenav`),
  appearing only once the window is wide enough (`min-width: 90rem`) to
  hold it without crowding the main content -- narrower windows simply
  don't get one, rather than it floating over the page.

## [5.4.2] - 2026-09-07

### Removed
- **Class R ("an ejected, geologically active world") removed entirely.**
  It had `h`/`e`/`c` all `False` -- zero probability weight, unreachable
  outside a manual `zone_override`. A genuinely free-floating/rogue planet
  (no host star at all) is a real exoplanet category, but doesn't fit this
  generator's star-centric zone model -- "which zone" is the wrong
  question for an object with no star to be zoned relative to. Rather than
  leave it as permanent dead weight, cut from `PLANET_CLASSES`,
  `PLANET_CLASS_PROBABILITIES` (already 0.0000), and `MOON_BLACKLIST`
  (`program_constants.py`); the now-vacuous
  `test_known_issue_class_with_no_valid_zone_is_unreachable` regression
  test removed from `src/tests/test_planets.py`.

### Fixed
- **Naming: triple-consonant validation bug.** `utils.is_name_valid`
  enforced "no more than two consecutive vowels" but had no matching check
  for consonants, so a fixed denylist of specific clusters
  (`BAD_CONSONANTS`) was the only thing standing between a generated name
  and an arbitrary run of 3+ consonants -- empirically ~40% of generated
  star names had one. Added a `consonant_count` run-length counter
  symmetric to the existing `vowel_count` one. Verified 0/5000 across
  star/planet/moon/sector name generation post-fix.
- **Web interface: navbar had no way back to the database picker on a
  single-database deployment.** `lib/page.py`'s sidenav only linked to
  `index.py` when more than one `.db` file existed, since `index.py`
  itself auto-redirects past the picker for a single database -- but that
  meant the common single-database production deployment
  (`starmap.moltenaether.com`) had no Databases link at all. Fixed by
  always linking to `index.py?all=1`, and teaching `index.py` to render
  the full picker table (bypassing its own auto-redirect) when `?all=1`
  is present.
- **Web interface: system-page table of contents wasn't collapsible.**
  `system.py`'s `_toc_html` now builds a `<details>`/`<summary>` pair
  instead of a plain `<aside>`/`<h3>`, collapsed by default -- a native,
  no-JS toggle.
- Minor cleanups found in passing: an unused `program_constants` import in
  `sectorGen.py`, and a stray `f`-string prefix on a non-interpolated
  string in `evolution.py`.

## [5.4.1] - 2026-09-07

### Added
- **`sectorGen.py --density`: a controllable density value for sector
  generation.** Sectors previously always generated a flat count of systems
  (`--num-systems`, default 10) with no connection to sector volume -- there
  was no way to make one generated sector meaningfully denser or sparser
  than another. `--density` is a float multiplier on the real local stellar
  density this codebase already models (`physical_constants.LOCAL_STELLAR_DENSITY_LY3`
  via `SpaceSector.expected_system_count()`) -- 1.0 means a realistic sector
  this size, 2.0 twice as dense, 0.5 half. The actual per-sector count is
  drawn with the existing `_sample_poisson_count` helper (previously only
  used by `SpaceSector.grow_from_seed`), resolved fresh for each sector
  rather than once at parse time, so counts vary naturally sector to sector
  under `--num-sectors` and across `galaxyGen.py` runs, matching a real
  spatial Poisson process. Mutually exclusive with `--num-systems`; neither
  flag given keeps the original flat-10-systems default unchanged.
  `galaxyGen.py` gets the flag for free since it shares
  `sectorGen.add_shared_generation_options`/`validate_shared_generation_args`
  with no changes needed on its side.

## [5.4.0] - 2026-09-07

### Removed
- **Two brown-dwarf-scale gas-giant classes cut entirely.** Correcting
  their radius ranges to real sub-stellar physics in [5.3.9] left them as
  near-duplicates of each other (overlapping radius, differentiated only by
  an invisible density number no description text reflects), each carrying
  a vanishingly small (0.01%) generation weight, representing objects that
  aren't really planets in the first place (brown dwarfs are sub-stellar
  objects; this codebase already has a dedicated mechanism for stellar/
  sub-stellar companions via `BinaryStarProxy`/`BINARY_SYSTEM`). Their
  combined generation weight folded into the ordinary Jupiter/Saturn-class
  gas giant. All remaining references to the removed classes (and to the
  earlier-removed small-rocky-class pair) scrubbed from comments/docs,
  generically describing what they explain rather than naming
  now-nonexistent class letters.

### Added
- **Per-class `size_mode`: bell-curve size distributions for planets and
  moons.** Radius generation (`planetPhysics.generate_planet_properties`'s
  three radius-draw sites, and `generate_moons`' moon-radius draw) switched
  from a flat uniform draw across each class's declared radius range to a
  bounded Gaussian ("bell curve") draw peaking at a new per-class
  `size_mode` value (0.0-1.0, "what fraction through the available range is
  the statistically most common size") via new `utils.sample_bounded_bell`.
  For a moon, "available range" is its actual Hill-sphere/mass-capped
  window, not necessarily the class's full declared range -- `size_mode` is
  read relative to whatever range is actually being drawn from. Every
  surviving class's `size_mode` is set from a real-world single-body analog
  where one exists (Earth for M, Mars for K, Venus for N, Jupiter/Saturn
  for J, Uranus/Neptune for I/T, the real rocky-to-gaseous transition
  radius for V) or a reasoned default (small-body populations, both real
  asteroids/KBOs (Class C/D) and this generator's own hot-zone rocky
  classes (A/B), skew toward their smaller end, matching real
  size-frequency distributions). Verified empirically: generated radius
  means land within ~1% of each class's real-world anchor (e.g. Class M
  mean 6,337km vs Earth's 6,371km; Class K mean 3,410km vs Mars's 3,389.5km).
  The self-adjusting spread (`sample_bounded_bell`'s `spread_divisor`)
  keeps a legible bell shape even for a mode pinned near one edge of a
  range, via rejection sampling rather than clamping (which would pile
  spillover probability mass at the boundary).
- `PLANET_CLASS_PROBABILITIES`'s `J` weight increased 0.0529 -> 0.0531 to
  absorb the removed brown-dwarf-scale classes' combined weight.
- Atmospheric-pressure sanity bound (`test_full_matrix.py`/`test_planets.py`)
  lowered back down `1e9` -> `5e7` Pa now that the brown-dwarf-scale
  gravity/pressure extremes those classes produced are gone (observed max
  across the full star-type matrix is now ~1.1e7 Pa, from Class N).
- `test_planet_physics_fixes.py`'s gravity/pressure correlation test
  threshold recalibrated (Spearman `>0.5` -> `>0.1`) for the narrower
  gravity range the remaining gas giants span without the removed classes'
  two-orders-of-magnitude density spread; its density-blend-skip test
  rewritten to inject a temporary `density_range` onto an existing class
  via `monkeypatch` rather than depend on a specific class declaring one
  (no class currently does -- it's generic, reusable override
  infrastructure, same as every other per-class override this codebase has
  built up).

## [5.3.9] - 2026-09-07

### Removed
- **Classes X and Y removed, merged into B and A/B respectively.** Both
  were small, redundant variants of existing hot-zone rocky classes:
  - Class X ("a stripped core from a gas giant", no atmosphere, radius
    500-5000km) folded into Class B ("a small, molten world") as an
    alternate origin story ("occasionally the stripped core of a former
    gas giant") rather than a separate atmosphere-less class -- a
    Mercury-analog this close to its star already has only a negligible
    exosphere in reality, so B's existing thin atmosphere covers X's "no
    atmosphere" identity closely enough.
  - Class Y ("a 'demon' class world", radius 5000-7500km) folded into
    *both* Class A and Class B (both radius ceilings raised 5000->7500km
    to absorb Y's size range), its toxic/irradiated flavor folded into
    each class's description as a variant rather than kept as a third
    class. B's composition gained "and sulfur" for the shared
    volcanic/irradiated theme.
  - `MOON_BLACKLIST` and `PLANET_CLASS_PROBABILITIES` updated to drop X/Y
    (their combined generation weight folded into A/B rather than dropped).

### Changed
- **Gas-giant zones reworked using real exoplanet science.** Previously
  every gas/ice-giant class (I/J/S/T/U) was valid in zone `c` only --
  meaning no gas giant could ever appear close to its star or in the
  habitable zone, despite "hot Jupiters" and "warm Jupiters" being a
  standard, well-documented real classification (orbital period < 10 days
  / 10-365 days respectively) and Neptune-mass planets in or near a
  temperate zone being common too.
  - Class J now valid in `h`/`e`/`c` (hot/warm/cold Jupiter) -- a warm/cold
    Jupiter placed in zone `e` can also generate ordinary moons via the
    existing moon-generation path, including habitable-class ones (the
    "habitable exomoon around a giant planet" trope), verified working.
  - Class I now valid in `e`/`c` (warm/cold Neptune), deliberately **not**
    `h`: real close-in Neptune-mass planets are rare (the observed "hot
    Neptune desert") because stellar irradiation photoevaporates a
    Neptune-mass H/He envelope down to a bare rocky/metal core well before
    it could stay Class I -- that outcome is exactly Class B's newly-merged
    "stripped core" identity (see Removed, above).
  - Classes S/T/U left `c`-only: real directly-imaged super-Jovian/brown-dwarf
    companions are predominantly found at wide separations (formation and
    detection-bias reasons), and close-in high-mass companions, while known,
    are much rarer (the "brown dwarf desert").
  - Fixing `generate_moons` (planetPhysics.py) to respect
    `HABITABLE_WORLD=False` the same way direct planet generation already
    does -- unreachable before this change (no gas giant could ever be in
    zone `e`), but a warm/cold Jupiter placed there could otherwise roll a
    habitable-class moon even in a system explicitly configured to disallow
    habitable worlds.
- **Classes S ("supergiant") and U ("ultragiant") radius ranges corrected
  to real sub-stellar physics** -- previously 250,000-50,000,000km and
  25,000,000-60,000,000km respectively, i.e. up to ~86 solar radii, larger
  than most actual stars. Real brown dwarfs stay within ~15% of Jupiter's
  own radius (69,911km) across their *entire* 13-80 Jupiter-mass range
  (electron degeneracy pressure means more mass compresses them, R ~
  M^(-1/8)); even the smallest true hydrogen-fusing red dwarf stars are
  only ~0.1 solar radii (~69,600km, essentially Jupiter-sized). Corrected to
  S: 60,000-120,000km, U: 65,000-130,000km (deliberately overlapping --
  that overlap *is* the real physics, not an oversight). What actually
  differentiates a higher-mass sub-stellar object at essentially the same
  radius is **density**: real measured brown-dwarf densities run roughly
  10-200 g/cm^3 (~10-150x an ordinary gas giant's), so both classes gained
  a new `density_range` (S: 10.0-60.0, U: 60.0-150.0 g/cm^3). Class T ("gas
  dwarf") also corrected, 250,000-25,000,000km -> 15,000-55,000km -- its
  old floor matched S's own floor, letting a "dwarf" be exactly as large as
  a "supergiant"; now spans ice-giant-to-Saturn scale, meaningfully below
  both J and S.
- **New per-class `density_range` override** (`PLANET_CLASSES`, read by
  `planetPhysics.get_planet_mass_ranges`/`generate_planet_properties` and
  mirrored in `plausibility.theoretical_gravity_bounds_g`) -- same
  `.get(..., default)` pattern as `atm_molar_density_range` etc. Used by
  Classes S/U above.
- **Fixed a gas-giant density-blend bug that silently collapsed every gas
  giant's density to a near-zero, physically meaningless value**, previously
  flagged as necessary-but-deferred follow-up work in
  `test_gas_giant_sampled_densities_are_finite_positive_and_within_theoretical_bounds`'s
  own docstring. `generate_planet_properties`' core/envelope harmonic-mean
  density blend reused `planet.atm_density` (drawn from
  `ATMOSPHERE_DENSITY["g"]`) as the envelope term -- but that value is a
  thin, surface/pressure-layer density (a different physical layer, correct
  for the separate atmospheric-pressure/scale-height calculation), ~1000x
  lighter than a real gas-giant envelope's actual bulk density. A harmonic
  mean is dominated by whichever term is smaller almost regardless of mass
  fraction, so this collapsed density to ~0.001-0.003 g/cm^3 for every gas
  giant, independent of core density -- silently defeating the new S/U
  `density_range` work above (their much denser cores had almost no effect
  on the final blended density). Fixed by introducing
  `physical_constants.GAS_ENVELOPE_BULK_DENSITY` (0.06-0.3 g/cm^3, grounded
  in real measured "puffy" gas giants -- WASP-193b's ~0.06 g/cm^3 is the
  lowest confirmed bulk density known), used only for this blend. A class
  declaring its own `density_range` (S/U) now skips the blend entirely and
  uses that density directly -- real brown dwarfs don't have a meaningfully
  separate light envelope over a denser core the way an ordinary gas giant
  does, so blending toward a light "puffy" value would just dilute the
  elevated density right back down.
- Atmospheric-pressure sanity bound (`test_full_matrix.py`/`test_planets.py`)
  raised `5e7` -> `1e9` Pa -- Classes S/U's now-correct brown-dwarf-like
  gravity legitimately pushes pressure up to ~3.3e8 Pa across the full
  star-type matrix.
- `test_planet_physics_fixes.py`'s two gas-giant-blend tests updated to
  match the new formula (queuing an `envelope_density_gcm3` draw instead of
  reusing `atm_density`), plus a new
  `test_density_range_override_skips_the_blend` covering the S/U direct-
  density path.

## [5.3.8] - 2026-09-07

### Fixed
- **Habitable/life-bearing classes restricted to the ecosphere zone.** Class
  Q (`h`/`e`/`c` all `True`) and Class W (`h`/`e` `True`) both carry a
  `life_chemical`, but weren't restricted to zone `e` like every other
  life-bearing class already was — meaning a life-bearing Q or W world could
  be generated directly in the hot or cold zone. Both are now `e`-only
  (`h`/`c` `False`); Q's "eccentric orbit" flavor still holds fully confined
  to the ecosphere zone, and W's "tidally locked" identity is arguably more
  scientifically apt restricted to the ecosphere zone (real tidally-locked
  *habitable* worlds are a real, actively studied trope specifically because
  a cool star's habitable zone sits close enough in for tidal locking to be
  near-guaranteed — e.g. TRAPPIST-1's planets, Proxima b).
  New regression test `test_life_bearing_classes_are_ecosphere_only`
  (`src/tests/test_planets.py`) locks this in for every class with a
  `life_chemical`, not just the ones on the separately-curated
  `HABITABLE_PLANET_CLASSES` list.
- **Description/atmosphere-text cleanup across `PLANET_CLASSES`.** Several
  classes' `description` field repeated a word the render template
  (`planetData.to_paragraph_list`) already supplies via its own "with an
  atmosphere of {atmosphere}" or "with a composition of {composition}"
  clause, producing genuinely broken rendered sentences — e.g. Class B
  previously rendered "...a small, molten world **with a thin atmosphere
  with an atmosphere of** a mix of helium, sodium, and oxygen...". Fixed for
  Classes B, E, K, N, Q, W, X, and Y — descriptive detail that belonged on
  the atmosphere field (e.g. "thin", "dense, reducing") was moved there
  instead of dropped, and Y's atmosphere field (previously the only one not
  ending in a named gas mixture) reworded to "a turbulent, toxic, and
  irradiated mix of gases".

## [5.3.7] - 2026-09-07

### Changed
- **Greenhouse-formula fix and per-class climate tuning.** The prior
  `greenhouse_factor` formula (`planetPhysics.calculate_atmospheric_conditions`)
  used `atm_molar_density` as its only lever, scaled by a single shared
  `CO2_MAX_GREENHOUSE_FACTOR` cap (5) — real Earth air's own molar mass
  already produced `greenhouse_factor ≈ 3.33` under that formula, driving
  Class M's mean surface temperature to 362K instead of ~288K, and Mars and
  Venus (nearly identical real molar mass, ~100x different real greenhouse
  forcing) could never be told apart by molar mass alone.
  - `CO2_MAX_GREENHOUSE_FACTOR` (now 500) is a generous safety ceiling, not
    the calibration knob.
  - New per-class `PLANET_CLASSES` keys — `albedo_range` (extended beyond
    Class P, which introduced it), `atm_molar_density_range`,
    `atm_density_range`, and `greenhouse_multiplier_range` — give every
    tuned class independent control of composition (molar density),
    quantity (mass density), and potency (greenhouse multiplier), following
    the exact override pattern `albedo_range` established for Class P.
    Class N's old hardcoded `atm_density = 65` / `atm_molar_density = max`
    special case in `planetPhysics.generate_planet_properties` is folded
    into this same general mechanism.
  - **Tuned this pass** (via the new `climate_tuning_cli.py`, see below):
    M (Earth analog, ~286K/~99kPa), O (warm/wet ocean world, ~293K/~97kPa),
    H (hot/dry desert, ~325K/~42kPa), K (Mars analog, ~231K/~0.57kPa), L (K
    + vegetation, warmer/thicker than K, ~256K/~2.2kPa), N (Venus analog,
    ~740K/~9.4MPa), E/F/G (a young, cooling progression, ~373K -> ~329K ->
    ~292K, converging near M), and V (thick, hot Super-Earth, ~365K/~286kPa).
    Class P and W are unchanged this pass (P already had a working
    `albedo_range`; W's "extreme temperature variations" identity needs a
    day/night model this generator doesn't have, not just range tuning).
  - Text refinements: Class H's "and metals" -> "and mineral dust" (a real
    desert's atmosphere lofts particulate, not metal vapor); Class O's
    atmosphere text (previously byte-identical to Class M's) now mentions
    water vapor; Class E's vague "hydrogen compounds" now names a real
    Hadean/Archean-analog mix (water vapor, ammonia, methane); Class K's
    "carbon dioxide" -> "a thin mix of carbon dioxide and nitrogen" (names
    the real Mars-analog composition, not just class); Class V's atmosphere
    text now reflects its tuned CO2-retention (not primordial H/He) identity.
  - Class M's disabled gravity clamp (`planetPhysics.calculate_surface_gravity`)
    deleted outright (it had been commented out, inert, since the
    atmospheric-pressure formula fix; no longer worth keeping around "in
    case it needs restoring").
- **New: `src/tests/climate_tuning_cli.py`.** A human-driven tuning tool
  (mirrors `physical_plausibility_cli.py`'s batch-generate-and-report
  pattern, reusing `planetgen.generation.plausibility`'s engine): generates N
  bodies of one class, reports summary stats, and — for classes with a
  direct real-world analog (M/Earth, K/Mars, N/Venus) — a delta line against
  that reference. `--albedo`/`--molar-density`/`--density`/`--greenhouse`
  flags temporarily monkeypatch that class's `PLANET_CLASSES` entry for the
  run only, so candidate values can be iterated without editing source
  between runs.
- **New: `src/tests/test_climate_tuning.py`.** Regression suite locking in
  the tuning above via generously-toleranced bands (Class M/N/K within
  Earth/Venus/Mars-like ranges) and relative orderings (N hottest/
  highest-pressure of the tuned classes; K colder/thinner than L; L colder
  than M; H hotter/drier than O; E > F > G cooling progression converging
  near M; M < V < N) rather than brittle exact-value assertions, since
  generation is inherently stochastic.
- `test_full_matrix.py`/`test_planets.py`'s atmospheric-pressure sanity
  bound raised from `1e7` to `5e7` Pa — Class N now legitimately reaches
  ~9-15MPa depending on host star luminosity (previously ~2.98MPa mean, per
  the greenhouse-formula bug above), and the old bound was sized for the
  un-tuned, incorrectly-cold N. `test_planet_physics_fixes.py`'s
  `test_class_p_has_own_albedo_range_distinct_from_default` updated to
  check P's range differs from M's own new tuned range, rather than
  asserting M has no override at all.

## [5.3.6] - 2026-09-07

### Added
- **Galaxy-scale coordinate system (Track C), merged.** Sectors can now
  be placed on a galaxy-wide, radial shell/Fibonacci-sphere tiling
  (`docs/design/galaxy-coordinate-system.md` sections 0-8) instead of
  existing only in isolation:
  - **Schema v3 -> v4** (`stellarObjects/schema.sql`): `sectors` gains six
    nullable galaxy-frame columns (`center_x/y/z_pc`, `galactic_radius_pc`,
    `shell_index`, `shell_slot_index`), NULL together for a sector never
    placed in a galaxy. `_db.migrate_database`'s new `_migrate_v3_to_v4`
    handles the upgrade (existing `_migrate_v1_to_v2`/`_migrate_v2_to_v3`
    updated in lockstep, per that function's "each hop maps straight to
    the current schema" design); `docs/database-schema.md`'s schema
    history documents the change.
  - **`planetgen/galaxy/geometry.py`** (new): the shell/Fibonacci-sphere
    tiling primitives (`shell_sector_count`, `shell_radius_pc`,
    `sector_position_pc`) and `enumerate_sectors_within_radius` — an
    exact, two-prune generation-unit primitive that finds every sector
    address within a radius of an arbitrary galaxy-space point without
    ever scanning a whole shell's slot count (design doc section 8).
  - **`galaxyGen.py`** (new, repo root): a CLI generating many sectors as
    one galaxy, reusing `sectorGen.py`'s own per-sector generation/save
    path. `--shell K` batch-generates a whole radial shell (guarded by
    `LARGE_SHELL_WARNING_THRESHOLD`, needing `--limit`/`--yes` above it);
    `--center-sector ID --radius-pc R` generates a local neighborhood
    around an already galaxy-placed sector. Either mode skips slots a
    sector already occupies.
  - **`GALACTIC_CENTER_DISTANCE_LY` is now per-sector**, not a single
    fixed constant: `Star.calculate_system_perimeter`/
    `BinaryStarProxy._calculate_system_perimeter_static` accept a
    `galactic_center_dist_ly` override, threaded from `galaxyGen.py`
    through `StarSystem`/`Star`/`BinaryStarProxy`, falling back to the
    old fixed constant for unplaced/standalone sectors
    (`sectorGen.py`'s own CLI unaffected).
  - `stellarObjects/utils.py` gains `mpc_to_pc`/`pc_to_mpc` (exact) and
    `pc_to_ly`/`ly_to_pc` (display-string conversions), per the design
    doc's unit-choice section.
  - New end-to-end coverage: `src/tests/test_galaxy_gen.py` runs
    `galaxyGen.py`'s actual CLI entry point (both `--shell` batch mode
    and `--center-sector` local-neighborhood mode) against a real
    temporary database, asserting correct shell/slot addresses, correct
    stored positions, no duplicate slots, and that already-occupied
    slots are skipped on re-run — the gap `docs/TODO.md`'s Phase 4 entry
    flagged as not yet done when this was paused mid-session. Developed on
    a branch cut before Track A's completion and the TODO/FIXME-comment
    migration (5.3.5); merged into `main` after both, with no functional
    changes needed beyond a `test_db_migration.py` assertion that had to
    decompress the (now gzip-compressed, per Track B) v3->v4 migration
    backup the same way the existing v1->v2 backup test already did.

## [5.3.5] - 2026-09-07

### Fixed
- **Atmospheric pressure is no longer independent of gravity**
  (`planetgen/physics/planets.py`): the barometric-formula pressure
  calculation algebraically canceled gravity out entirely (`atmospheric_pressure
  = atm_density * R * T / atm_molar_density`), so a Neptune-gravity gas giant
  and a Jupiter-gravity one produced the same pressure. A new
  `_atmosphere_retention_factor(gravity_g)` (linear, normalized to 1.0 at
  Earth gravity) now scales an *effective* atmospheric density used only in
  the pressure calculation (not `planet.atm_density` itself, which also
  feeds the gas-giant density blend), reintroducing a real, tunable
  gravity/pressure relationship (Spearman correlation on a mixed
  terrestrial/gas-giant sample now > 0.5, vs. ~-0.11 before).
- **Class P ("cold, glaciated") is colder than Class M again**: gave Class P
  its own `albedo_range` (0.5-0.7, matching real ice/snow Bond albedo)
  instead of sharing the default rocky/Earth-like range (0.12-0.35) with
  every other terrestrial class. Previously P and M were statistically
  indistinguishable in temperature once the disabled clamp was removed (see
  `docs/analysis/habitability-atmosphere-sanity-review.md`); P's cold
  identity now comes from the unclamped physics instead of a post-hoc
  override.
- Completes Track A (see [5.3.4]'s gas-giant density/greenhouse fixes) --
  8682 tests passing, including 10 new regression tests in
  `src/tests/test_planet_physics_fixes.py`.

## [5.3.4] - 2026-09-07

### Fixed
- **Gas-giant density blend** (`planetgen/physics/planets.py`): the
  core/atmosphere blend used a mass fraction as an arithmetic-mean weight
  between two densities, which is dimensionally wrong and could produce
  gas giants as low as 0.026 g/cm^3. Replaced with the mass-weighted
  harmonic mean, the physically correct way to combine two component
  densities via a mass fraction. `plausibility.py`'s independently
  reimplemented copy of this formula (used to derive analytical
  hard-invariant gravity bounds) was updated in lockstep, including its
  docstring's justification for corner-evaluation (still valid: the new
  formula is monotonic in each argument, just no longer multilinear).
- **Inverted greenhouse factor** (`planetgen/physics/planets.py`): the
  formula rewarded an atmosphere for being *far* from CO2's molar density
  rather than for actually containing more CO2 — backwards from physical
  reality. Now scales directly with `atm_molar_density`, the only
  atmosphere-composition signal the data model has today.

### Changed
- **Database schema-migration backups are now gzip-compressed** and
  excluded from the web database picker and from a subsequent migration
  run (previously a plain `.db` copy that a naive `*.db` glob would both
  surface in the picker and silently re-migrate on the next run). See
  `docs/TODO.md`'s "File Management" section for detail.

### In progress, not yet merged (see `docs/TODO.md` for exact state)
- Atmospheric pressure is still algebraically independent of gravity, and
  Class M/Class P remain statistically indistinguishable — both scoped
  and partially started, paused mid-session in worktree
  `agent-a8acb02b98bed5b8d`.
- The galaxy-scale coordinate system's 8 open design questions were
  decided this session, and a schema v3->v4 migration,
  `GALACTIC_CENTER_DISTANCE_LY` fix, and a first `galaxyGen.py` were
  written but paused uncommitted in worktree `agent-a36f801e275fb2b71`
  before a full test pass — not part of this release.

## [5.3.3] - 2026-09-07

### Changed
- **Markdown consolidated into `docs/`, renamed for clarity.** Only
  `README.md`/`LICENSE.md`/`CHANGELOG.md` remain at the repo root; every
  other README and loose doc moved into `docs/` with a descriptive name:
  `TODO.md` -> `docs/TODO.md`, `src/html/README.md` -> `docs/html-interface.md`,
  `db/README.md` -> `docs/database-schema.md`,
  `src/api/README.md` -> `docs/api.md`, `apache/README.md` ->
  `docs/apache-deployment.md`, `examples/EXAMPLES.md` ->
  `docs/example-systems.md`, `examples/JSON.md` ->
  `docs/system-file-format.md`. Every cross-reference between them (and
  from code/scripts) was updated to match; a couple of pre-existing
  broken/mismatched links were caught and fixed along the way
  (`docs/system-file-format.md`'s "full command-line reference" link text
  didn't match its own target; `src/html/lib/dbutil.py`'s docstring still
  pointed at the pre-5.3.2 `stellarObjects/` path instead of
  `src/stellarObjects/`).
- **`wsgi.py`, `queryDb.py`, `migrateDb.py` moved into `src/`**, alongside
  `stellarObjects`/`api`/`tests`, so only the `*Gen.py` scripts
  (`sectorGen.py`/`systemGen.py`) are visible at the repo root as CLI
  entry points. Since these three now sit as direct siblings of
  `stellarObjects`/`api` under `src/`, Python's own sys.path[0] (the
  running script's own directory) already makes those packages
  importable -- the sys.path shims 5.3.2 added to them are gone, no
  longer needed (unlike `sectorGen.py`/`systemGen.py`, which stay one
  directory further away at the repo root and keep theirs).
  `migrateDb.py`'s `DEFAULT_DB_DIR` needed an extra `os.path.dirname()`
  level to still resolve to the repo-root `db/`, one directory deeper
  than before.
- **`physicalPlausibility.py` moved to `src/tests/physical_plausibility_cli.py`**,
  matching that directory's naming scheme, without a `test_` prefix so
  pytest doesn't try to collect it as a test module (it's a human-facing
  report generator, not a pass/fail check -- `src/tests/test_physical_plausibility.py`
  remains the actual automated test).
- **`apache/` moved into `examples/apache/`** (its `README.md` moved to
  `docs/apache-deployment.md` per the markdown-consolidation rule above);
  **`examples/*.json` moved into `examples/systems/`**, so `examples/`
  now holds two clearly-separated subfolders (`apache/` deployment
  config, `systems/` recipe files) instead of a flat mix of both kinds of
  example content. `.gitattributes`' LF-pinning rule for
  `apache/*.sh` was updated to `examples/apache/*.sh` to keep matching
  the actual file.
- `setup.py` and `pytest.ini` were **evaluated for a move into `src/` and
  kept at the repo root** -- both are genuinely not possible without
  breaking things, verified empirically rather than assumed:
  `pip install -e .` with `setup.py` moved silently built a bogus,
  empty `UNKNOWN-0.0.0` package instead of erroring (pip's PEP 517
  build only looks for `setup.py`/a full `pyproject.toml` project table
  at the invocation root, and this repo's `pyproject.toml` only declares
  a build backend, no project metadata of its own to fall back on); with
  `pytest.ini` moved, `pytest`'s own config-file discovery (which only
  searches the invocation directory and its parents, never a
  subdirectory) silently fell back to `pyproject.toml` and ignored
  `testpaths`/`pythonpath` entirely -- tests still happened to pass
  either way (coincidentally, via `pytest`'s own unrelated `__init__.py`-walkup
  sys.path behavior and an editable install's global registration of the
  root-level `py_modules`), which is exactly the kind of silent,
  environment-dependent fragility not worth introducing on purpose.

## [5.3.2] - 2026-09-07

### Changed
- **Repo layout: `stellarObjects`/`api`/`tests` moved under a new `src/`
  directory** (`src/stellarObjects/`, `src/api/`, `src/tests/`), so only
  the top-level CLI entry points (`sectorGen.py`, `systemGen.py`,
  `queryDb.py`, `migrateDb.py`, `physicalPlausibility.py`, `wsgi.py`) are
  visible at the repo root. `setup.py` now declares an explicit
  `package_dir` per discovered package rather than a blanket
  `package_dir={'': 'src'}`, since that would have also redirected the
  root-level `py_modules` lookups (`systemGen`/`sectorGen`) into `src/`,
  where they don't live. Every root entry script gained a small
  `sys.path` shim (inserting `src/` before its `stellarObjects`/`api`
  imports) so they keep working without requiring `pip install .` first,
  matching the no-install fallback `src/html/`'s CGI scripts already relied
  on -- those fallbacks (`src/html/lib/dbutil.py`, `src/html/sector.py`,
  `src/html/search.py`) were updated the same way. Two internal
  repo-root-relative path computations (`stellarObjects/webconfig.py`'s
  `_PROJECT_ROOT`, `stellarObjects/_db.py`'s `DEFAULT_DB_PATH`) needed an
  extra `os.path.dirname()` level to still resolve correctly one
  directory deeper; `pytest.ini` gained an explicit `pythonpath = . src`
  so both the entry scripts and the moved packages resolve during tests
  regardless of pytest's own import-mode heuristics.
- **`webconfig.json.example` moved into `src/html/`**; the real, gitignored
  `webconfig.json` stays at the repo root, outside Apache's `src/html/`
  `DocumentRoot`, for the same security reason `db/` already lives there
  (see `docs/webconfig.md`).
- **`WEBCONFIG.md` moved into `docs/`**, alongside this session's
  `docs/design/`/`docs/analysis/` additions, consolidating loose
  documentation in one place (`README.md`/`LICENSE.md`/`TODO.md`/
  `CHANGELOG.md` stay at the repo root, and per-directory READMEs
  `src/html/README.md`/`apache/README.md`/`db/README.md`/`src/api/README.md`
  stay next to the code they document).

## [5.3.1] - 2026-09-06

### Fixed
- **Evolved-star mass sampling could imply a pre-Big-Bang star.** An
  evolved-class star (`Yerkes != V`, e.g. a giant or supergiant) derives
  its required main-sequence lifespan from its own already-generated mass
  (`Star._calculate_initial_star_age_and_lifespan`'s evolved-star branch);
  for a sub-solar-mass progenitor (roughly under ~0.88 Msun), that implied
  lifespan alone already exceeded `UNIVERSE_AGE_GY` (13.8 Gy) -- meaning
  such a star couldn't actually have finished its main-sequence phase yet
  in the real universe. The [5.3.0] universe-age fix deliberately didn't
  paper over this by capping age below its own required floor, since that
  would produce a self-contradictory star (e.g. a red giant younger than
  its own progenitor's main-sequence lifespan); this is the deeper fix it
  called for. `Star.generate_star`'s evolved-star mass sampling now uses a
  new `_sample_evolved_star_mass_sol` helper (`planetgen/generation/star.py`)
  that rejects and resamples (not clamps, which would just pile an
  artificial spike at the cutoff) any candidate mass whose implied
  main-sequence lifespan would exceed `UNIVERSE_AGE_GY`, capped at
  `program_constants.EVOLVED_STAR_MASS_MAX_RESAMPLE_ATTEMPTS` (100)
  attempts before raising `ValueError` -- in practice a no-op resample for
  every class but III (Giant), whose 0.8-8 Msun range straddles the
  cutoff. Yerkes class VI (subdwarf) is deliberately excluded, since its
  entire allowed mass range (0.1-0.8 Msun) sits below the cutoff and would
  reject every draw; that's tracked as a separate, still-open modeling
  question in `TODO.md`.

## [5.3.0] - 2026-09-06

### Fixed
- **Atmospheric pressure was ~5 orders of magnitude too low for every
  planet class except M.** `planetPhysics.calculate_atmospheric_conditions`
  summed "shell" volumes derived from `planet.radius`/`scale_height`
  (stored in km) directly against `planet.atm_density` (kg/m³) without
  converting km³→m³, undercounting `atmospheric_mass` by ~10⁹×; an
  unexplained `* 7500` fudge factor only clawed back about 4 of those ~9
  orders of magnitude. This rendered as "0.0 kPa" for every atmosphere-
  bearing class but M (Class M never showed it, because a hardcoded clamp
  force-overrode its pressure into an Earth-like range regardless of what
  was computed). Replaced the whole shell-integration loop with the
  closed-form barometric formula for an isothermal, hydrostatic atmosphere
  (`P_surface = ρ_surface · g · H`), which needs no arbitrary zone count or
  fudge factor and lands within the right order of magnitude for both
  Earth-like and Mars-like test cases. `tests/test_planets.py`/
  `tests/test_full_matrix.py`'s pressure assertions were tightened from a
  no-op `>= 0` check into real physical sanity bounds (1 Pa – 10 MPa) so a
  regression like this can't silently pass again.
- **Stars could be reported as hundreds of billions of years old** (e.g.
  "918.77 Billion Years old"), which is impossible since the universe
  itself is only ~13.8 billion years old — even though the star's age never
  actually exceeded its own (very long, and realistically so: real M dwarfs
  are predicted to live trillions of years) lifespan. Added a
  `UNIVERSE_AGE_GY = 13.8` ceiling, applied in both the main-sequence and
  evolved-star branches of `Star._calculate_initial_star_age_and_lifespan`
  and in `adjust_age_for_planets`'s age-adjustment logic. (White dwarfs
  already had their own bounded 0.1–12.0 Gy cooling-age range and needed no
  change.) A sub-solar-mass evolved-star progenitor whose own
  main-sequence lifespan already exceeds the universe's age is a separate,
  deeper mass-sampling issue, noted in `TODO.md` rather than papered over
  here.
- **Rare `IndexError: string index out of range` crash in name
  generation.** `generate_phoneme_salad_name` (used for star, planet, and
  moon names) could crash when a base name containing a literal apostrophe
  (e.g. `PLANET_NAMES`'s "Hi'iaka") landed adjacent to a `UNIVERSAL_PHONEMES`
  chunk also ending in an apostrophe (e.g. "ch'"), producing a literal `''`
  in the assembled name; the final capitalization step split on `'` and
  indexed `part[0]` on every piece, crashing on the empty piece between the
  two apostrophes. Not specific to binary systems — any name generation
  call could hit it given enough attempts, which is why it only surfaced
  rarely across a large batch of generated systems.

### Changed
- Class M's (and Class P's) hardcoded gravity/pressure/temperature
  "forcing" clamps in `planetPhysics.py` are commented out, not deleted —
  the fixed atmospheric-pressure formula already lands close to realistic
  ranges on its own, so the band-aid that had been masking the pressure bug
  for Class M is no longer needed. May be restored later if Class M/P need
  tighter narrative guarantees again.
- `star_systems` gained a `location` column (schema bumped to v3): for any
  system placed in a sector, its sector's name plus distance (in
  light-years) to its 3 nearest neighboring systems, e.g. `"Voranthis
  Kelmoor -- nearest: Aldenar (4.2 ly), Brekthos (7.8 ly), Corvane (9.1
  ly)"`. Computed once at write time (mirroring how `quadrant` is already
  derived and persisted) via the existing `SpaceSector.nearest_neighbors`/
  `distance_between` helpers; `NULL` for systems never placed in a sector.
  Existing v1/v2 databases migrate automatically (backed up first, as
  usual) and have their `location` backfilled from already-stored position
  data.
- The web interface's Search link is now reachable from every page (it
  previously required detouring back through the database-browse page) —
  `src/html/lib/page.py`'s shared page header now includes it whenever a
  database is selected.
- The database-picker landing page (`src/html/index.py`) now redirects straight
  to browsing the one database present, if the configured database
  directory contains exactly one `.db` file, instead of always showing the
  picker.

### Added
- `webconfig.json` (repo root, gitignored — `webconfig.json.example` is the
  committed template): site-level configuration, currently `site_name` and
  `base_url`, plus unused placeholder fields (`db_username`, `db_password`,
  `db_name`) reserved for a possible future non-SQLite backend. Kept
  outside `src/html/`'s served document root, the same way `db/` already is.
  See [`docs/webconfig.md`](docs/webconfig.md) for full documentation.

## [5.2.5] - 2026-09-06

### Changed
- **Database schema bumped to v2**: moons split out of the shared
  `planets` table into their own `moons` table (`planet_id` FK to the
  planet they orbit, plus their own `moon_evolutionary_paragraphs`/
  `moon_reflection_spectrum` child tables), instead of self-referencing
  via `planets.parent_planet_id`/`is_moon`. This is what `src/html/search.py`'s
  attribute tags needed to actually tell a planet from a moon: previously
  a "Class D Planet" tag queried `planets` with no way to exclude
  `is_moon=1` rows of the same class, so it silently listed moons too.
  `stellarObjects/_db.py` gained a dedicated `insert_moon` (mirroring
  `insert_planet`, but writing to the new table); `insert_planet` no
  longer recurses into itself for moons.
- `src/html/search.py`'s tag facets and name search now follow the same
  split: "Planet Class"/"Planet Body Type"/"Planet Supported Life
  Chemistry" only ever match top-level planets, with an identically-
  shaped "Moon Class"/"Moon Body Type"/"Moon Supported Life Chemistry"
  set of tags (and a separate autocompleting "Moon name" field) for
  moons -- each still only rendering the tag buttons for values actually
  present in the chosen database. `src/html/system.py`'s planet/moon table
  now reads moons from the new table instead of a recursive
  self-join, and no longer needs to recurse (moons never generate their
  own moons -- confirmed by the existing
  `tests/test_moons.py::test_moons_cannot_themselves_have_moons`).

### Added
- `migrateDb.py` (repo root): converts every `*.db` file in a directory
  from schema v1 to v2, backing up each original first
  (`<name>.db.v1-backup-<timestamp>.db`) -- a no-op for a database
  that's already current. Backed by a new `stellarObjects._db.
  migrate_database`/`_migrate_v1_to_v2`, which builds the converted
  database at a temporary path and only atomically swaps it into place
  at the very end, so a crash or error partway through leaves the
  original file (and the backup already made) untouched. `install.sh`
  now runs this automatically (a new step, right after installing the
  Python package) over every database in `db/`, so `sudo ./update.sh`
  keeps an existing deployment's database working across the schema
  change with no manual step.
- `tests/test_db_migration.py` and `tests/test_db_persistence.py`: the
  first automated tests of `stellarObjects/_db.py`'s save/migrate paths
  (previously "manually smoke-tested via `sectorGen.py`" per `TODO.md`).
  Between the two, they cover the v1->v2 migration (including the
  moon-owned child-row split, which the real generated data used for
  manual smoke-testing happened not to exercise) and a genuine
  generate-then-save-then-query round trip through `insert_star_system` --
  the latter caught a real bug during development (`insert_moon`'s
  `INSERT` had one more column than value placeholder).

## [5.2.4] - 2026-09-06

### Added
- `src/html/search.py`: a faceted/name search page for the web interface.
  Two complementary ways to find an object in the chosen database:
  - Click-to-filter tag buttons for object type (star/planet/moon/
    asteroid belt), star spectral class and Yerkes luminosity class,
    planet class and body type, and supported life chemistry -- each
    group renders a button only for values actually present in that
    specific database (e.g. no Yerkes-luminosity buttons beyond "Main
    Sequence" if nothing else was generated), per `SELECT DISTINCT ...
    GROUP BY` queries against `stars`/`planets`/`asteroid_belts`, not a
    fixed enumeration. Multiple tags within one group OR together (any
    matching value); tags across groups AND together, except that
    selecting an explicit Object Type tag acts as a master filter (e.g.
    selecting only "Stars" hides the Planets panel even if a
    planet-class tag also happens to be selected). Clicking a tag
    toggles it via a plain link that rewrites the query string, so this
    works with JavaScript entirely disabled, same as every other page in
    `src/html/`.
  - A name search with one field each for sector, star system, star, and
    planet/moon names, each with its own HTML5 `<datalist>` for
    autocomplete (native browser suggestions, no JavaScript, populated
    from that entity's own distinct names). Asteroid belts have no name
    of their own (per `schema.sql`'s own note on this), so they're only
    reachable via the "Asteroid Belt" object-type tag.

  Every filter is a parameterized SQL query (`IN (...)`/`LIKE ... ESCAPE
  '\'`), and each result panel is capped (300 rows, with a "showing the
  first N" note) to keep a broad tag click from dumping an entire large
  database into one page. Linked from `browse.py`'s breadcrumb.

## [5.2.3] - 2026-09-06

### Changed
- `install.sh` no longer installs the Python package via the deprecated
  `python3 setup.py install`; it now runs a proper, build-isolated `pip
  install --upgrade --force-reinstall "$SCRIPT_DIR"` instead.
  `setup.py install` broke on a deployed Ubuntu 20.04/Python 3.8 host
  *twice* in a row, both times because setuptools' own vendoring shim
  (`extern`) prefers a real, already-installed copy of a dependency it
  vendors over its own newer bundled copy whenever one is importable —
  first `importlib_metadata` (`AttributeError: ... no attribute
  'EntryPoints'`), then, after that was patched around, `packaging`
  (`TypeError: canonicalize_version() got an unexpected keyword argument
  'strip_trailing_zero'`) — because that distribution's old apt-provided
  copies of each in turn shadowed the working vendored one. Chasing each
  shadowed dependency individually only fixes the one that broke that
  day, not the next one; `pip install`'s build isolation builds the
  package in a throwaway environment that can't see the system's
  site-packages at all, closing the whole bug class instead of patching
  it dependency-by-dependency, and means `install.sh` no longer needs to
  upgrade the system's global `setuptools`/`importlib_metadata` at all
  for this step. `--force-reinstall` is deliberate: a plain `pip install
  .` skips reinstalling when pip thinks the same version is already
  installed — true on every `update.sh` run between version bumps in
  `planetgen/_version.py` — which would otherwise silently leave the
  previous run's install in place instead of the source `update.sh` just
  pulled.
- `install.sh` and `update.sh` now re-`chmod +x` every `*.sh` file in the
  repo (not just `src/html/*.py` and one hardcoded `apache/set-permissions.sh`
  path), including themselves. A `core.fileMode=false` git config on the
  authoring machine drops the executable bit on any file type on
  checkout, not just `src/html/`'s `.py` scripts, and `update.sh` invokes
  `install.sh` directly (`"$SCRIPT_DIR/install.sh"`, not `bash
  install.sh`) — if a pull had dropped *its* executable bit, the shell
  would refuse to exec it before `install.sh`'s own permission fix ever
  got a chance to run. `update.sh` now fixes this itself, immediately
  before invoking `install.sh`, so a dropped bit on `install.sh` (or on
  `update.sh` for its own next run) no longer requires a manual `chmod
  +x` to recover from. (One unavoidable exception: the very first
  `update.sh` run to pull this fix is still executing the old script
  text it already had open when `git pull` replaced the file on disk out
  from under it, so that one migration may still need a manual `chmod +x
  install.sh` — every run after that starts from the fixed script.)

## [5.2.2] - 2026-09-06

### Added
- `src/html/system.py` now renders a system's description as actual HTML by
  default (`?view=rendered`), via a new small, purpose-built Markdown-to-
  HTML converter (`src/html/lib/mdconvert.py`) targeting exactly the narrow
  Markdown subset `StarSystem.__str__` generates (headers, pipe tables,
  paragraphs, `<sup>` exponents) — not a general-purpose parser, and no
  new dependency. The original raw wikitext/Markdown source is still one
  click away (`?view=source&format=...`) for copy-pasting into a wiki.
  Escapes every block in full before emitting markup, then narrowly
  re-enables only the one legitimate raw-HTML pattern generated content
  contains, so a mischievous `--name`/`--star-type` value can't inject
  live HTML into a rendered page (covered by new tests in
  `tests/test_mdconvert.py`).
- `src/html/static/style.css` rewritten as a small design system: CSS custom
  properties, automatic light/dark via `prefers-color-scheme`, card-style
  panels, badges, breadcrumbs, a proper type scale, hover states,
  focus-visible outlines, and a responsive breakpoint. Applied
  consistently across every page (`index.py`, `browse.py`, `sector.py`,
  `system.py`, and the shared error page in `lib/page.py`), not just the
  system page.
- `update.sh` (repo root, alongside `install.sh`): pulls the latest
  changes from git and re-runs `install.sh` so permissions stay correct
  afterward. Refuses to run over uncommitted local changes, and pulls
  with `--ff-only` (fails loudly rather than creating a surprise merge
  commit if history has diverged) instead of a plain `git pull`.

### Fixed
- `apache/set-permissions.sh` and `install.sh` both failed to make every
  `*.py` file under `src/html/` executable by `www-data` — `install.sh`'s own
  `chmod +x` step used `find -maxdepth 1`, silently skipping
  `src/html/lib/*.py`, and `set-permissions.sh` only `chgrp`'d (group
  ownership) rather than `chown`'d (user *and* group) the deployed
  directories. Both fixed: the `-maxdepth 1` restriction is gone, and
  `set-permissions.sh` now `chown -R`s to the detected Apache user:group
  and reports how many `.py` files it made executable, so a wrong path is
  obvious rather than silently matching nothing.

## [5.2.1] - 2026-09-06

### Added
- `src/html/`: a small, dependency-free web interface (plain Python CGI
  scripts, standard library only) for browsing the SQLite databases from
  `db/README.md` — pick a database, drill into its sectors and star
  systems, and view (or copy, via a `<textarea>`) the rendered
  wikitext/Markdown page saved for each one. See `src/html/README.md`.
- `apache/`: deployment tooling for the web interface —
  `planetgen.conf.example` (an example Apache2 virtual host) and
  `set-permissions.sh` (detects the user/group Apache2 actually runs as
  and sets ownership/permissions on the deployed `src/html/`/`db/`
  directories accordingly).
- `install.sh`: a one-shot Linux installer tying the above together —
  runs `setup.py install`, pre-fetches the NLTK `words` corpus into a
  shared, world-readable location, makes the CGI scripts executable,
  enables Apache's `cgid` module, and runs `apache/set-permissions.sh`;
  prints the one remaining manual step (copying/enabling the example
  vhost) rather than touching Apache's site configuration itself.
- `planetgen/names/wordlists.py` gained `UNIVERSAL_PHONEMES`: a pool of ~100
  short, ASCII-7-bit-printable phoneme chunks romanized from roughly 18
  language families (Romance, Germanic, Slavic, Arabic/Hebrew/Persian,
  Turkish, South Asian, Mandarin, Japanese, Korean, Vietnamese,
  Austronesian, Polynesian, Bantu, Mesoamerican, Andean, Celtic,
  Finno-Ugric, Caucasian). `generate_phoneme_salad_name` now splices one
  of these into every generated name — star, planet, moon, and sector
  alike — with `UNIVERSAL_PHONEME_CHANCE` (40%) odds, widening the
  cultural range generated names are drawn from beyond each type's own
  base name list.

### Changed
- Removed `SECTOR_DESIGNATORS` from `planetgen/names/wordlists.py` (dead code
  — defined but never referenced anywhere).
- `generate_phoneme_salad_name` gained an `allow_split` parameter
  (default `True`, unchanged for stars/planets/moons); see the sector
  name bug fix below.

### Fixed
- Sector names could come out as 3-4 words instead of the intended 2 (in
  a 500-name stress test, this happened 95% of the time).
  `generate_phoneme_salad_name` can split a long result into two
  space-separated words on its own (`split_long_word`), and
  `sectorGen.generate_sector_name` already joins two independent calls
  into one name — an internal split on either half silently produced
  3-4 words in the final result. `generate_sector_name` now passes the
  new `allow_split=False` for both halves.
- `planetgen/names/wordlists.py` unconditionally called
  `nltk.download('words', quiet=True)` at import time; `nltk`'s
  `download()` always targets the *current user's* default download
  directory and attempts to create it, regardless of whether the corpus
  already exists elsewhere on `nltk.data.path` — this broke the web
  interface entirely under Apache's locked-down `www-data` account
  (`PermissionError: ... '/var/www/nltk_data'`), since that account has
  no writable home directory. It had only ever "worked" for the CLI
  tools because those ran as `root`. Now checks
  `nltk.data.find('corpora/words')` first and only downloads on
  `LookupError`; `install.sh` pre-fetches the corpus into a shared,
  world-readable path so that check always succeeds once installed.
- `apache/set-permissions.sh` hard-failed when it couldn't detect
  Apache's user/group (e.g. because apache2 isn't started/enabled yet,
  which is exactly the case during a fresh `install.sh` run) — now warns
  and defaults to the standard Debian/Ubuntu `www-data:www-data` instead.
- CGI scripts could be deployed non-executable regardless of a local
  `chmod +x`, because a `core.fileMode=false` git config on the
  authoring machine silently drops the executable bit before it reaches
  a commit. `install.sh` now `chmod +x`s the CGI scripts and
  `apache/set-permissions.sh` directly on every install, independent of
  whatever mode git happened to store.

## [5.2.0] - 2026-09-06

### Added
- Full SQLite database persistence: `stellarObjects/schema.sql` defines the
  schema, and the new `stellarObjects/_db.py` writes an already-generated
  `SpaceSector` (every system, star, planet, moon, and asteroid belt,
  plus a rendered copy of the wiki page in both wikitext and Markdown)
  into it in a single transaction. Documented column-by-column in the new
  `db/README.md`. `sectorGen.py` now calls this automatically on every
  run, saving to `db/planetgen.db` by default (overridable via the new
  `--db-path` option); the database file itself is gitignored.
- Sector name generation: sectors now get a random two-word name (e.g.
  "Voranthis Kelmoor") drawn from a new sector-flavored name list
  (`SECTOR_NAMES`/`SECTOR_PREFIXES`/`SECTOR_SUFFIXES` in
  `planetgen/names/wordlists.py`, pulling from real galactic structures and
  well-known science-fiction sector names) instead of reusing star names
  with "Sector" appended.
- `Star`, `Planet`, `BinaryStarProxy`, and `AsteroidBelt` each gained a
  `get_table_properties()`/`get_composition_summary()` method that returns
  the same already-formatted values their `to_paragraph_list()` renders
  into text, so the new database layer can store exactly what was
  published without duplicating any formatting logic.

### Changed
- **Breaking:** `sectorGen.py`'s `-n` short flag now means `--name`
  (hard-sets the sector's own name) instead of `--num-systems`, which no
  longer has a short flag; the old `--sector-name` option was renamed to
  `--name`.
- Distances stored in the database use a two-tier unit convention
  (milliparsecs for sector-scale position/geometry, kilometers everywhere
  else); `stellarObjects/utils.py` gained `ly_to_milliparsecs`/
  `milliparsecs_to_ly` and `physical_constants.py` gained
  `AU_PER_PARSEC`/`AU_PER_MILLIPARSEC` to support the conversion, used only
  by the persistence layer.

### Fixed
- Atmospheric scale height (`planetPhysics.calculate_atmospheric_conditions`)
  was computed in meters but combined directly with `planet.radius` (in
  km) without converting units first, throwing off atmosphere thickness
  and volume for every planet with an atmosphere.
- `utils.calculate_object_mass` and `Planet.__init__` each independently
  (and inconsistently) converted radius-in-km to a volume, one of them
  mixing up its km/m factor; `Planet` now takes its volume solely from
  `calculate_object_mass`'s corrected calculation instead of recomputing
  it a second time.

## [5.1.0] - 2026-09-06

### Added
- `--version` option on both `systemGen.py` and `sectorGen.py`, printing the
  program's version, this repository's URL
  (https://github.com/dwhagar/planetGen), and a license summary, then
  exiting immediately.
- `planetgen/_version.py`: a single, dependency-free source of truth
  for the project's version number, shared by both CLI scripts and
  `setup.py`.

### Changed
- `setup.py` now reads its `version` from `planetgen/_version.py`
  instead of a hardcoded, never-updated placeholder.
- README updated with current version information and a link to this
  changelog.

## [5.0.1] - 2026-09-06

### Fixed
- Flavor text (both system-level and planet/moon-level) was being rolled —
  and shared `system_config` counters mutated — every time a system or planet
  was *rendered* (`__str__`/`to_paragraph_list()`) instead of once at
  generation time, so rendering the same object twice would silently re-roll
  and double-count it. It is now decided exactly once per object during
  generation (`StarSystem.__init__`), and rendering is a pure, idempotent
  read of the already-decided text.

### Changed
- Internal property cleanup across `Planet`, `Star`, and `AsteroidBelt`:
  removed five redundant star-property snapshots from `Planet` that
  duplicated data already reachable via its `star` reference (one of
  which — `star_radius` — had a latent double-unit-conversion bug whenever a
  moon was generated); renamed ambiguous or colliding attributes
  (`evolution` → `evolutionary_speed`, `hab` → `habitable_zone`, and
  `type` → `body_type` on `Planet`/`AsteroidBelt` so it no longer shares a
  name with the unrelated `Star.type`); filled in missing/stale docstring
  documentation.

### Added
- Regression tests proving rendering is idempotent (calling
  `to_paragraph_list()`/`__str__()` twice produces identical output and
  performs no further mutation).

## [5.0.0] - 2026-09-05 (evening)

### Added
- `spaceSector.py`: a new `SpaceSector`/`SectorSystemEntry` layer that places
  generated `StarSystem`s at `(x, y, z)` positions within a cubic sector, with
  save/load to JSON (`SystemConfig` gains `to_dict`/`from_dict` to support
  this).
- Poisson-disk-based sector growth for realistic inter-system spacing, plus
  named-location ("quadrant") formatting for a system's position within a
  sector.
- `TODO.md`: a long-term roadmap toward full database storage and a web
  interface.

### Changed
- Reworked sector growth to sample candidate positions and then fine-tune
  each one against every existing neighbor, rather than a single-pass
  placement.
- Sector density and the minimum allowed separation between two systems are
  now grounded in real astronomical data — density from real local stellar
  surveys, and minimum separation from each system's own gravitationally
  derived Hill-sphere radius — rather than arbitrary constants.

## [4.1.0] - 2026-09-05 (afternoon)

### Added
- JSON system-file support: a system's exact contents (star type, forced
  features, per-orbit slots) can now be specified via a JSON recipe file,
  with several example systems included under `examples/`.
- The project's first automated test suite (pytest), covering example
  systems, moon generation, and planet generation.
- `sectorGen.py` CLI entry point, ahead of the sector-generation logic it
  would call the next day.

### Changed
- Major standardization pass across CLI options and internal call
  signatures.
- Renamed the main script from `planetGen.py` to `systemGen.py`, reflecting
  that it now generates full star systems rather than a single planet.

## [4.0.0] - 2026-09-03

### Changed
- Major internal refactor: generation logic was split out of the `Planet`
  class into dedicated `planetPhysics.py` (physical/orbital generation) and
  `planetLife.py` (life chemistry and evolution) modules, so `Planet` itself
  holds state and presentation while generation logic lives alongside it as
  free functions.
- The monolithic `constants.py` was split into `physical_constants.py`
  (real-world physical constants) and `program_constants.py`
  (generation/tuning constants), for readability and maintainability.
- Import structure cleaned up throughout the package; install script
  (`setup.py`) updated.

## [3.0.0] - 2026-07-07

### Added
- Binary star systems: a new `doubleStar.py` module (`BinaryStarProxy`)
  represents a double star system as a single effective star (combined
  mass/luminosity/habitable zone) for the purposes of planet placement.

### Fixed
- Primary/secondary star identification and overall system age calculation
  for binary systems.
- An age/lifespan bug affecting planets.

### Changed
- Binary star output formatting refined; the large-star-forcing option
  updated to work correctly alongside binary generation; README updated.

## [2.2.3] - 2026-06-27

### Added
- More flavor text variants across different planet classes.

### Changed
- Refined repeat-prevention and selection logic for flavor text so the same
  text is less likely to recur in quick succession.

## [2.2.2] - 2026-06-26

### Fixed
- Flavor text selection options.

### Changed
- Switched random-number generation to Python's `secrets` module for
  higher-quality entropy, avoiding similar successive sequences.

## [2.2.1] - 2026-06-25

### Changed
- Wording and grammar adjustments throughout the generated text output.

## [2.2.0] - 2026-06-24

### Added
- Options to specify a system's age and to force (or forbid) intelligent
  life.
- Flavor text: a random chance of extra descriptive "sensor" text being
  appended to a system or a planet's description.

### Changed
- Centralized magic numbers into `constants.py` as named constants, making
  the physics/generation formulas easier to read.
- Unified the set of properties referenced across all object types onto a
  shared convention.
- README updated for the new CLI options.

## [2.1.1] - 2026-06-23

### Changed
- Asteroid belt data and logic split out into its own module
  (`asteroidData.py`).
- Expanded asteroid belt composition descriptions and grammar.

### Added
- Option to specify a custom name for the star system.

## [2.1.0] - 2026-06-21

### Added
- Stellar age is now generated and reported for every system.
- Evolutionary timeline narratives: systems now estimate the plausibility of
  life at every stage, from simple single cells through technological
  civilizations.

### Fixed
- Assorted output-formatting issues.

## [2.0.0] - 2026-06-20

### Added
- Life chemistry system: planets and star types now carry information about
  the chemical processes that could plausibly give rise to life, as a
  foundation for the evolutionary modeling that follows.

### Changed
- Command-line options reworked so multiple flags combine correctly together.

## [1.5.0] - 2026-06-12

### Added
- Option to generate a star system with no planets.
- Option to specify a star's exact spectral type from the command line.

## [1.4.0] - 2026-06-06

### Added
- Markdown export support, and the ability to write generated output to a
  file (in addition to the console).
- Heliosphere radius and stellar gravitational-influence ("system perimeter")
  calculations, describing the outer boundaries of a system.

### Fixed
- Name-generation issues.
- Scientific-notation formatting issues.

## [1.3.1] - 2026-06-05

### Changed
- Further cleanup of generated text output for clarity.
- Refinements to procedural naming.

## [1.3.0] - 2026-06-04

### Added
- Procedural name generation for stars, planets, and moons.
- Project README.

### Changed
- General code cleanup pass.

*(~23-month gap in development between 2024-07-15 and 2026-06-04.)*

## [1.2.2] - 2024-07-15

### Fixed
- Star age/temperature calculations corrected to consistently use Kelvin
  throughout.

## [1.2.1] - 2024-07-14

### Fixed
- Class M worlds are now reliably clamped to Earth-like gravity and
  atmospheric pressure.
- Class P (and other) planet classes now get appropriate temperature
  treatment for their class.

### Changed
- Asteroid belt sizing now makes better use of the habitable zone.

## [1.2.0] - 2024-07-13

### Added
- Options to force a habitable world, a large star, and/or an asteroid belt
  (forcing both a habitable world and an asteroid belt together
  automatically forces a large star).
- Options to control the overall size (object count) of a generated system.
- Expanded command-line help text.

### Fixed
- A long-standing infinite-loop bug in star generation.
- Several edge cases in generation logic.

### Changed
- Capped the maximum number of system objects at 500 to prevent runaway
  generations.
- Refined the system description text and radius output formatting.

## [1.1.1] - 2024-07-11

### Fixed
- Completed and debugged the moon system: moon orbital distances are now
  tracked and calculated correctly, and remaining orbital overlap between
  planets and asteroid belts was removed.

### Added
- Weighted probabilities for moon class selection.

### Changed
- Further wording, grammar, and formatting refinements.

## [1.1.0] - 2024-06-29

### Added
- First implementation of the moon-generation function (not yet wired into
  planet creation or tested at commit time).

## [1.0.2] - 2024-06-28

### Added
- Mass validation for a specified planet, laying the groundwork for moon
  generation.

## [1.0.1] - 2024-06-26

### Fixed
- Orbital placement now accounts for planetary position and Hill radius, so
  systems no longer generate overlapping or too-closely-spaced orbits.
- Asteroid belt inner/outer bounds no longer overlap neighboring planets.

### Changed
- Further readability tweaks to the generated text output.

## [1.0.0] - 2024-06-25

### Added
- First fully working end-to-end planet generator, validated against real
  Earth reference values ("I think it finally works!").
- Special-cased handling for Class N worlds (boosted atmospheric density and
  molar mass).

### Changed
- Multiple formatting and descriptive-text passes on the generated output.

## [0.2.0] - 2024-06-23

### Added
- Core stellar physics: working star generation and planet-count estimation.
- Initial atmospheric modeling: mass, gravity, and atmosphere calculations.

### Fixed
- Numerous early bugs in surface pressure and temperature calculations.

## [0.1.1] - 2024-06-22

### Fixed
- Continued debugging to get the initial prototype running end-to-end (moved
  development from VS Code to PyCharm along the way).

## [0.1.0] - 2024-06-21

### Added
- Initial project scaffolding: the `stellarObjects` package with the first
  `Planet`, `Star`, and `StarSystem` classes, and a `main.py` entry point.
- First rough pass at generating a single planet's basic properties.