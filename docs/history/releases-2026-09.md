# Release headlines, 2026-09

182 releases. Each headline is the single best feature that release brought (CHANGELOG.md, newest first). Releases without a headline predate the rule (2026-10-10).


## 2026-09-30

- **7.9.2**: Admin Generate jobs work on native Windows.
- **7.9.1**: Fix binary secondary outweighing its primary.
- **7.9.0**: Galaxy Map wedge lines, scale readout and one solid slice-able density block model.
- **7.8.0**: Warp times follow Boss's warp curve, plus dimensional fold travel times.
- **7.7.0**: pip installs only locked, hash-checked files.
- **7.6.1**: Fix star systems placed inside black hole Hill spheres.
- **7.6.0**: Deployment guides for every major platform and reverse proxy support.
- **7.5.0**: Random first admin password, rate limits and broad security hardening.
- **7.4.0**: Brute-force property-based tests across the whole codebase.
- **7.3.0**: Planets and moons named after their star, with API renaming.
- **7.2.1**: Fix asteroid belt overlapping the next planet.
- **7.2.0**: update.sh checks instead of reinstalling, with migrate-or-delete prompt and progress.
- **7.19.1**: Fix random binary star naming failure.
- **7.19.0**: Letter classes for nebulae, supernova remnants and asteroid fields.
- **7.18.1**: Moons now orbit within their planet's stable Hill sphere limits.
- **7.18.0**: Interstellar objects generated at real-world rates, plus brown dwarfs and runaway stars.
- **7.17.0**: Galaxy Map sizes sector blocks from the screen and zooms to single sectors sooner.
- **7.16.1**: Fix planet Hill spheres overlapping asteroid belts.
- **7.16.0**: Windows installers and macOS support for the bash installers.
- **7.15.0**: A supermassive black hole at the center of every galaxy.
- **7.14.1**: Failed logins locked out per username with growing delays.
- **7.14.0**: NAV courses read as bearing and mark on nested reference frames.
- **7.13.1**: Fix close binary planets orbiting inside the binary.
- **7.13.0**: Hybrid master-wedge sector slots that line up from center to edge.
- **7.12.0**: System page lists everything in orbit once, in order.
- **7.11.0**: Lighter top bar with a settings gear and responsive menu.
- **7.10.3**: Galaxy Map distances use the shared distance ladder.
- **7.10.2**: Upper bounds on admin generation inputs.
- **7.10.1**: Renames check only uniquely named objects for clashes.
- **7.10.0**: Every distance shown in its most meaningful unit, with exact constants.
- **7.1.0**: CGI pages removed, old URLs redirect in the app.
- **7.0.1**: Fix lopsided half-sphere neighborhood on the 3D Galaxy Map.
- **7.0.0**: One 4 parsec sector standard and a galaxy of aligned layers.
- **6.6.0**: System and phenomenon pages moved to the new site.
- **6.5.0**: Admin pages moved to the Flask app with CSRF protection.
- **6.4.0**: Galaxy Map moved to the Flask app at /galaxy.
- **6.3.0**: One-off star systems generated from the admin site.
- **6.2.0**: Sector and NAV pages moved to bookmarkable URLs.
- **6.1.0**: Generate, plan and reset the galaxy from the web interface.
- **6.0.0**: Galaxy Map draws the real cylindrical sector grid.
- **5.59.0**: Browser accessibility checks for every Flask page.
- **5.58.0**: install.sh works on externally managed Python.
- **5.57.0**: Search page moved to bookmarkable URLs with filters as parameters.
- **5.56.0**: Quasars at the galactic center.

## 2026-09-27

- **5.55.0**: Galaxy Map shades predicted density with lit cylindrical prisms computed in the browser.

## 2026-09-24

- **5.54.0**: Home page served by Flask with a new header, search box and mobile menu.
- **5.53.2**: Airless bodies lose leftover atmospheres and temperatures never drop below the cosmic background.
- **5.53.1**: Removed old Galaxy Map endpoint and versioned the map script against stale caches.
- **5.53.0**: Light, dark and system theme switch, plus a proper phone layout.
- **5.52.0**: System pages render natively from the database as an expandable list of bodies.
- **5.51.0**: Every stellar phenomenon type now appears on the Sector Map and in sector listings.
- **5.50.0**: Program-wide debug log recording every generator decision, random roll, query and request.
- **5.49.2**: Editing a sector refreshes only affected Galaxy Map tiles, not the whole cache.
- **5.49.1**: Every list on the site pages 50 rows at a time with one pager.
- **5.49.0**: Admin stats page showing server health, database statistics and duplicate-name report.
- **5.48.3**: Zoomed-in Galaxy Map no longer draws unfilled sectors and density cloud as a ball.
- **5.48.2**: Browse sector list and sector system tables sorted outward by distance.
- **5.48.1**: Galaxy Map clicks, double-clicks and scroll zoom now land where you point.
- **5.48.0**: Created and modified timestamps on sectors, systems and phenomenon tables.
- **5.47.3**: Server checklist document for verifying a deployment after an update.
- **5.47.2**: Install and update scripts create the Galaxy Map tile cache directory.
- **5.47.1**: Fixed giant sector spheres, missing orbit lines and flat star glows on the maps.
- **5.47.0**: Galaxy Map loads space as cached cube tiles, preventing server crashes.
- **5.46.35**: TODO list regrouped and numbered in working order.
- **5.46.34**: TODO list now tracks open work only, with new items added.
- **5.46.33**: Release notes in changes folder replace PR version bumps, with automatic version stamping.
- **5.46.32**: System Map measure-distance tool giving real distance and shortest path around the star.
- **5.46.31**: Sector Map shows clickable indicators for every neighboring sector, generated or not.
- **5.46.30**: Sector Map shows stars, nebulae and remnants as textured, glowing 3D spheres.
- **5.46.29**: Galaxy Map colors generated sectors by their real stellar density.
- **5.46.28**: Galaxy Map caps rapid-click server requests to protect the API.
- **5.46.27**: Browse sector and standalone system tables are truly paginated.
- **5.46.26**: Rogue planets and interstellar comets now appear in the Phenomena listing and detail pages.
- **5.46.25**: Database migration now correctly applies the spatial index schema upgrade.
- **5.46.24**: System Map stars get granulation texture and brighter glow; spheres hide orbit lines.
- **5.46.23**: Sector Map compass arrow now labeled simply N.
- **5.46.22**: Selected sector stays visible on the Galaxy Map when zoomed in close.

## 2026-09-23

- **5.46.9**: Compact one-line page header on system, sector, galaxy and phenomenon pages.
- **5.46.8**: Better System Map label collision avoidance and 3D body preview inside the map.
- **5.46.7**: Planet civilization narrative no longer contradicts its class description.
- **5.46.6**: Binary system star tables render correctly and generated belts no longer overlap.
- **5.46.21**: Reworked Galaxy Map click and double-click controls, with soft spherical density cloud.
- **5.46.20**: Apache example no longer doubles connection load through an HTTPS redirect loopback.
- **5.46.19**: Sector pages no longer time out the server thanks to spatial phenomenon indexes.
- **5.46.18**: Health endpoint reports a never-migrated database instead of failing.
- **5.46.17**: Health endpoint reports whether the database schema is current.
- **5.46.16**: Galaxy Map view query uses a new spatial index to stop site-wide timeouts.
- **5.46.15**: Galaxy Map page now renders the 3D map, with a recognizable galaxy density cloud.
- **5.46.14**: Every System Map star, planet and moon renders as its own live 3D sphere.
- **5.46.13**: Interactive 3D Galaxy Map with live viewport queries and copyable generate commands.
- **5.46.12**: Zoomable real-scale diagrams on every phenomenon page, with supernova remnants fully supported.
- **5.46.11**: Real zoom and pan on the Galaxy Map.
- **5.46.10**: Navigate-to-here and navigation routes to and from standalone phenomena.

## 2026-09-20

- **5.46.5**: Text search no longer crashes with a server error.

## 2026-09-19

- **5.46.4**: Sector generation keeps its systems when a sector runs out of room.
- **5.46.3**: Debug runs report per-phase generation timing for profiling.
- **5.46.2**: Galaxy random-start mode no longer shows a bogus single-sector progress bar.
- **5.46.1**: Progress bar restored for galaxy generation runs.
- **5.46.0**: Qualifying sectors can no longer come out completely empty.
- **5.45.1**: Name-uniqueness bookkeeping uses about a third fewer database queries, speeding generation.
- **5.45.0**: Galaxy random-start mode can require a minimum starting density.
- **5.44.1**: Removed the sector progress bars from generate command output.
- **5.44.0**: Removed per-sector progress bar and fixed progress bar flicker.
- **5.43.1**: End-to-end tests of galaxy generation against a real planned skeleton.
- **5.43.0**: Galaxy generation skips empty low-density sectors and reports per-sector star statistics.
- **5.42.0**: Rotating 3D body preview sphere with atmosphere glow on the System Map.
- **5.41.0**: Sector Map rebuilt as a real WebGL 3D scene with true perspective.
- **5.40.0**: Navigation parameters hidden from the address bar by posting links as forms.
- **5.39.1**: Continuous integration fixed and test connection-pool leak plugged.

## 2026-09-18

- **5.39.0**: Galaxy-wide unique names for sectors, systems, planets and moons.
- **5.38.0**: Wide binary System Map shows both stars' planets, with binary naming fixes.
- **5.37.0**: Publish systems and sectors to Wiki.js or MediaWiki from the web interface.
- **5.36.0**: Browsable list and detail pages for exotic phenomena.
- **5.35.7**: Search page no longer times out on large databases.
- **5.35.6**: Sector system counts now follow real spiral-galaxy density instead of a flat ten.
- **5.35.5**: Galaxy-placed Sector Map and far-from-core Galaxy Map open readable.
- **5.35.3**: Site title uses the configured site name and the database picker is gone.
- **5.35.1**: Superscript exponents display correctly in the System Map infobox.

## 2026-09-17

- **5.35.0**: All generators unified into a single generate command.
- **5.34.0**: One wiki client library supporting both Wiki.js and MediaWiki publishing.
- **5.33.0**: Science-based exotic phenomena seeded in every sector, plus galaxy random-start generation.

## 2026-09-16

- **5.32.0**: Two-body barycentric orbits so binary stars, planets and moons orbit a shared center.

## 2026-09-12

- **5.31.0**: True-position System Map and galaxy-placed nebulae and asteroid fields.
- **5.30.0**: Standalone asteroid fields and galactic orbital motion for all exotic phenomena.
- **5.29.0**: Exotic phenomena: black holes, neutron stars, nebulae, supernova remnants, rogue planets, comets.
- **5.28.0**: Command-line search for planets and moons by class and size.
- **5.27.0**: Wide S-type binary star systems with independent stars and planets.

## 2026-09-11

- **5.26.0**: Binary star chance now follows real multiplicity rates by spectral class.
- **5.25.0**: Planet counts derived from protoplanetary disk physics rather than a curve fit.
- **5.24.0**: Planet spacing uses mutual Hill radii, plus search by body size.
- **5.23.0**: One unified config.json file for all deployment settings.
- **5.22.0**: Planets no longer end up stranded outside the zone their class needs.
- **5.21.0**: Ecosphere planet classes now generate at class-appropriate distances within the habitable zone.
- **5.20.0**: Update script skips package reinstall when there is nothing new.
- **5.18.0**: Systemd timer for the periodic orbit update.
- **5.17.0**: Binary stars' mutual-orbit position advances over time.
- **5.16.0**: Stars orbit the galaxy and binary pairs orbit each other.

## 2026-09-10

- **5.9.1**: Class P climate tuned to a cold, glaciated world.
- **5.9.0**: Class W tidally locked world class removed.
- **5.8.3**: Provisional sector designations packed into a single reversible hex number.
- **5.8.2**: Provisional sector designations for ungenerated galaxy addresses.
- **5.8.1**: NAV Map plots the origin, destination, and route.
- **5.8.0**: NAV feature: course, distance, warp times, and optimal routes between systems.
- **5.15.0**: Floating-point guard skips useless orbital phase updates.
- **5.14.0**: Planets and moons get 3D positions and orbital speeds.
- **5.13.0**: Every star system gets a galactic orbital speed and period.
- **5.12.0**: Web interface becomes a thin frontend over the Flask API.
- **5.11.1**: Moon tidal locking now follows a physics-based despinning timescale.
- **5.11.0**: Planets and moons get real 3D orbits, rotation, and live orbital motion.
- **5.10.1**: Subdwarf stars no longer get ages older than the universe.
- **5.10.0**: Flask API moved into the web interface tree.

## 2026-09-09

- **5.7.0**: Galaxy density skeleton merged into the MySQL backend.
- **5.6.0**: API rate limiting and stubbed write endpoints.
- **5.5.0**: MySQL backend replaces SQLite, plus a paginated, validated API.
- **5.4.8**: Galaxy-wide density skeleton with lazy, visit-triggered sector generation.
- **5.4.7**: Galaxy sectors get gap-free vertices shared exactly with neighbors.

## 2026-09-07

- **5.4.6**: System page table of contents works at any window width.
- **5.4.5**: Sector Map 3D scene stays confined to its panel.
- **5.4.4**: Interactive 3D Sector Map with drag to rotate and zoom.
- **5.4.3**: System page table of contents moved to a right-margin rail.
- **5.4.2**: Removed the rogue-planet Class R and fixed triple-consonant names.
- **5.4.1**: Controllable sector density for generation.
- **5.4.0**: Bell-curve planet and moon sizes centered on real-world analogs.
- **5.3.9**: Gas giants can now appear in hot and habitable zones, with realistic brown-dwarf sizes.
- **5.3.8**: Life-bearing planet classes restricted to the ecosphere zone.
- **5.3.7**: Per-class climate tuning gives Earth, Mars and Venus analogs realistic temperatures and pressures.
- **5.3.6**: Galaxy-scale coordinate system with galaxyGen.py to generate many sectors as one galaxy.
- **5.3.5**: Atmospheric pressure now depends on gravity, and glaciated worlds run colder than Earth-like ones.
- **5.3.4**: Corrected gas-giant density and greenhouse physics, plus compressed migration backups.
- **5.3.3**: Repository tidied with all documentation consolidated into docs and scripts moved into src.
- **5.3.2**: Source code reorganized under a new src directory.

## 2026-09-06

- **5.3.1**: Evolved stars no longer sampled with impossibly old progenitor masses.
- **5.3.0**: Atmospheric pressure fixed for all planet classes, with star ages capped at universe age.
- **5.2.5**: Moons split into their own database table with a migration tool and separate moon search.
- **5.2.4**: Faceted tag and name search page for the web interface.
- **5.2.3**: More reliable installer using build-isolated pip install.
- **5.2.2**: System pages render as HTML with a redesigned light and dark web interface.
- **5.2.1**: Web interface for browsing saved databases, with Apache deployment and one-shot installer.
- **5.2.0**: Full SQLite database persistence of generated sectors, with named sectors.
- **5.1.0**: Version option on both command-line tools with a single version source.
- **5.0.1**: Flavor text decided once at generation so rendering is repeatable.

## 2026-09-03

- **4.0.0**: Generation logic refactored into separate planet physics and life modules.
