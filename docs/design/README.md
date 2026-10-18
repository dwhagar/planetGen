# Design notes

Every note in this folder, with what it is about and where it stands. The status words:

- **built**: the design is in the code.
- **partly built**: some of it is; the note's own status line says which.
- **plan**: decided, not built yet.
- **research**: findings and recommendations; read the note's status line for what was taken up.
- **reference**: facts and measurements, not a design.
- **history**: describes a design that was replaced; kept for the reasons.

A note's own status line, at its top, is the authority; this page only points at it. A test (`src/tests/test_docs_links.py`) fails when a note here is missing from this page. The architecture of the whole program is [architecture.md](architecture.md); the decisions behind it are in [design-decisions.md](design-decisions.md).

## The whole program

| Note | What it covers | Status |
|---|---|---|
| [architecture.md](architecture.md) | Every file and the main flows | reference (the code map is being rewritten, DOC.28) |
| [architecture-map.html](architecture-map.html) | The generated code map: areas, modules, functions and calls, with the main flows traced from the code (DOC.67) | generated page (open it in a browser; `scripts/build_architecture_map.py`) |
| [design-decisions.md](design-decisions.md) | The choices that shape the project, each with its status | reference |
| [later-phenomena-scatter-rules.md](later-phenomena-scatter-rules.md) | Which later phenomena are scattered, backfilled or fill-only (goal issue #1677): what each item stores, the magnetar scatter class, what the rule set and the profile gain, the columns for Core 1, and the changes to the TODO items | recommendation (Core 3, 2026-10-16; defaults stand, nothing built yet) |
| [todo-number-map.md](todo-number-map.md) | Old TODO numbers to category IDs, and the next free ID of each category | reference |
| [cache-warm-up.md](cache-warm-up.md) | The page cache shared through Redis, the data tags each answer reads and each change bumps, and the warm-up that fills it after an update (PERF.104) | built (Core 2, 2026-10-16) |
| [slow-page-answers.md](slow-page-answers.md) | A page request that is slow to build queues one job that stores the answer in the tagged page cache, with a cost rule, a status route, a node on the queue page, a cancel and a fallback (PERF.105) | built (Core 2, 2026-10-17; the search and phenomena pages use it, UX.118, #2089) |
| [library-migration.md](library-migration.md) | Moving to third-party libraries and into packages | built (queue being redone, PERF.86) |
| [math-and-library-efficiency.md](math-and-library-efficiency.md) | Where the heavy math runs; Python 3.12 and numba | partly built (PERF.88 open) |

## Galaxy, sectors and coordinates

| Note | What it covers | Status |
|---|---|---|
| [galaxy-coordinate-system.md](galaxy-coordinate-system.md) | The spherical shell, band and wedge sector grid | built (schema v91); the cylinder it replaced is in [archive/](archive/galaxy-coordinate-system-cylinder.md) |
| [galaxy-types.md](galaxy-types.md) | The 9.4 headline: one family per known galaxy type (spiral, lenticular, elliptical, irregular, dwarf), drawn from the seed through scaling relations, and the chain that leads to it | plan |
| [superbubble-density.md](superbubble-density.md) | Superbubbles shape the stars: the seeded bubble field becomes a density modulation on the young and intermediate populations, with bounds, sector counts and positions, and the chain that builds it | plan |
| [composition-model.md](composition-model.md) | The master note of the 9.0 composition headline: what every body, star and atmosphere is made of, the chain from element to derived value, the research lanes and where each item sits | checked design, for Boss's approval (GEN.437); the loader is built (GEN.472) |
| [atmosphere-colour.md](atmosphere-colour.md) | Atmosphere colour, the 9.0 feature: one gas and compound table, the mixture that sets the molar density, stored albedo and greenhouse values, ozone, haze, dust and cloud, then the absolute colour of the air, with the chain that builds it | checked design, for Boss's approval (GEN.409) |
| [bulk-composition.md](bulk-composition.md) | The bulk composition of every body (GEN.426): the bulk fractions, the equation of state and the reference solver, the draws, density, mass and gravity, moons and small bodies, the inverse and the migration | checked design, for Boss's approval (GEN.426) |
| [stellar-composition.md](stellar-composition.md) | The composition of every star and remnant (GEN.435): X, Y, Z and the elements from the stored [Fe/H], age and population, and what the star's iron, magnesium and silicon give the planets | checked design, for Boss's approval (GEN.435) |
| [radiation-model.md](radiation-model.md) | The radiation model (GEN.436): the bands, transmission and greenhouse from the infrared opacity, shielding by column, the ground's dose and radon, the surface ultraviolet | checked design, for Boss's approval (GEN.436) |
| [sky-colour-analysis.md](sky-colour-analysis.md) | The read-only study behind the atmosphere colour: what exists, the method, the effort and how well each source was checked | study |
| [atmosphere-sources-check.md](atmosphere-sources-check.md) | The check of the real-world sources behind the atmosphere colour chain: what was read, what is second-hand, what Boss must research | study |
| [bulk-composition-sources-check.md](bulk-composition-sources-check.md) | The check of the mass-radius and equation-of-state sources behind the bulk composition chain: constants, method order, the core fraction from the star and what Boss must research | study |
| [tunable-gaps.md](tunable-gaps.md) | Whether the unsourced composition weights can be seeded random draws tuned to published populations: the class of each of the 48 open atmosphere rows, the shared selection models, the 17 builder groups, the bracket corrections and what Boss must supply | study |
| [atmosphere-composition.md](atmosphere-composition.md) | The atmosphere composition research (GEN.431): the findings that change the design and the code, the one gas and compound registry and its consumers, the sources reached and refused, and what Boss must research | study |
| [atmosphere-composition-parts/](atmosphere-composition-parts/01-gas-table-thermodynamics.md) | The fifteen parts behind it: [gas table and thermodynamics](atmosphere-composition-parts/01-gas-table-thermodynamics.md), [optics](atmosphere-composition-parts/02-gas-optics-uv-visible.md), [infrared and greenhouse](atmosphere-composition-parts/03-infrared-greenhouse-thermal.md), [photochemistry, lifetimes and escape](atmosphere-composition-parts/04-photochemistry-lifetimes-escape.md), [mixing ratios by class](atmosphere-composition-parts/05-mixing-ratios-by-class.md), [condensates, clouds, hazes and dust](atmosphere-composition-parts/06-condensates-clouds-hazes-dust-minerals.md), [radiation shielding and ground sources](atmosphere-composition-parts/07-radiation-shielding-and-ground-sources.md), [verification of Boss's research note, A](atmosphere-composition-parts/08-verification-of-boss-added-research-A.md) and [B](atmosphere-composition-parts/08-verification-of-boss-added-research-B.md), [weighted composition tables](atmosphere-composition-parts/09-weighted-composition-tables.md) and the check of Boss's second note ([A](atmosphere-composition-parts/10-verification-second-note-A.md), [B](atmosphere-composition-parts/10-verification-second-note-B.md), [C](atmosphere-composition-parts/10-verification-second-note-C.md), [D](atmosphere-composition-parts/10-verification-second-note-D.md) and [E, the NH3 lifetime](atmosphere-composition-parts/10-verification-second-note-E-nh3-lifetime.md)) | study |
| [body-composition-research.md](body-composition-research.md) | The body composition research (GEN.432, GEN.433): what stars, rocky and icy planets, giants, small bodies, moons and the ground are made of, the chain from the star to the derived values, the findings against the code, the five decisions, the work clusters and the invariants a test can assert | study |
| [body-composition-parts/](body-composition-parts/A-stars.md) | The fifteen parts behind it: [A stars and remnants](body-composition-parts/A-stars.md), [B rocky planets](body-composition-parts/B-rocky.md), [C volatile-rich worlds and giants](body-composition-parts/C-volatile-and-giants.md), [D small bodies and moons](body-composition-parts/D-small-bodies-moons.md), [E ground, radiogenic heat and shielding](body-composition-parts/E-ground-radiogenic-shielding.md), [F check of Boss's research note](body-composition-parts/F-boss-note-check.md), [G1 comets and photolysis lifetimes](body-composition-parts/G1-comets.md), [G2 radiogenic ratios and Thorngren](body-composition-parts/G2-radiogenic-thorngren.md), [G3 magma ocean and lava worlds](body-composition-parts/G3-magma-ocean-lava.md), [G4 ices, brines and licences](body-composition-parts/G4-ices-brines-licences.md) | study |
| [spherical-sector-coordinates.md](spherical-sector-coordinates.md) | The plan for spherical shells, bands and wedges | plan, built (the switch merged to main on 2026-10-15, #1607; the Galaxy Map ladder is MAP.205 and MAP.206) |
| [spherical-shell-sectors.md](spherical-shell-sectors.md) | What spherical-shell sectors would mean (the question behind the plan) | research |
| [spherical-ladder-meridians.md](spherical-ladder-meridians.md) | Meridians for the spherical block ladder | research |
| [spherical-block-ladder.md](spherical-block-ladder.md) | The Galaxy Map's block ladder (243/81/27/9/3) and windows on shells | built (MAP.205, MAP.206) |
| [galaxy-disk-density.md](galaxy-disk-density.md) | The galaxy density model | built |
| [fill-order-curves-and-core.md](fill-order-curves-and-core.md) | Fill order, region selection and the galactic core | research |
| [realistic-galactic-core.md](realistic-galactic-core.md) | A realistic galactic core (GEN.279, reopens GEN.155): the nuclear star cluster and disc as density terms, the resolution cap, tiered mass and luminosity limits, planets and orbits at the centre | recommendation (Core 3, 2026-10-16; defaults stand, nothing built yet) |
| [sampling-backfill-and-resume.md](sampling-backfill-and-resume.md) | Sampling, backfill, void filling, resumable runs and directives | research (resume on hold, PERF.86) |
| [void-fill-census.md](void-fill-census.md) | GEN.102 again on the spherical grid: the sparse cells, the halo beyond the outline, the recommendation | research |
| [placement-by-population.md](placement-by-population.md) | GEN.103: where each kind is placed by population, and the three type-draw rules still to build | design |
| [mass-backfill.md](mass-backfill.md) | Four rings of heavier stars around generated sectors | built |
| [multiple-galaxies.md](multiple-galaxies.md) | Neighbour galaxies and a second galaxy | research |
| [globular-clusters.md](globular-clusters.md) | Globular clusters in generation | partly built (GEN.159: cluster table, King tables, generator; GEN.160: cluster density in the sector gate) |
| [intergalactic-and-halo-contents.md](intergalactic-and-halo-contents.md) | What fills the sphere off the disk: the stellar halo, satellites, ejected stars, gas, and how to give every cell contents | research (reach approved; build plan for TODO) |
| [galaxy-laws-research.md](galaxy-laws-research.md) | The sourced ranges and laws for the stellar halo, the population ages and dispersions, the thick disk, the dwarf and satellite sizes and the cluster system (GEN.406), with the fifteen open requests | research |
| [gap-filling-plan.md](gap-filling-plan.md) | How the open gaps are filled where the public record is silent: the 57 gaps with a derived, bracketed, theory or placeholder fill, the stand-in registry convention, the first Gemini batch (unverified) and the work order | research |
| [papers-to-fetch.md](papers-to-fetch.md) | The papers behind refused publishers that Boss downloads himself, and the arXiv papers a research thread reads instead, with Crossref citations (gap-filling plan section 8, GEN.531) | reference |
| [vizier-courteau-comeron.md](vizier-courteau-comeron.md) | The VizieR search for the Courteau 2007 scatter (GR5) and the Comeron thick-disc tables (GR3): the Comeron 2018 fit, the raw Vc/sigma0 scatter and the corrected references | research |
| [gemini-batch1-check.md](gemini-batch1-check.md) | The claim-by-claim check of Gemini batch 1 (51 claims) against read sources, the figures that can be used and the requests it answers | research |
| [boss-papers-2026-10-17.md](boss-papers-2026-10-17.md) | What Boss's 28 uploaded papers fill in the gap plan, what stays open, where they contradict our files and the printing faults found | research |
| [open-gaps-closest-laws-2026-10-17.md](open-gaps-closest-laws-2026-10-17.md) | Every open gap (64 lines) with the nearest law or model, where to look, the free tool, the stand-in quality and the deciding rule; the citation errors found | research |
| [ad1-methane-sources-check-2026-10-17.md](ad1-methane-sources-check-2026-10-17.md) | AD1: where the free methane absorption data is (Oxford k-tables from Karkoschka and Tomasko 2010, PDS volume gbat_0001), what is still open (ammonia, water) and the corrections to the pasted source note | research |
| [numeric-precision-audit.md](numeric-precision-audit.md) | GEN.372: whether 64-bit numbers are enough for positions at galactic scale (yes, if no code subtracts two large absolute coordinates), the numbers behind it and every precision concern found | research (Core 1, 2026-10-17; nothing built) |
| [spherical-region-queries.md](spherical-region-queries.md) | Selecting a sphere of the galaxy with a database query: the cell key, key runs, a covering index, measured against the box, the spatial index and the sector walk | research (2026-10-15; nothing built) |
| [phenomena-catalogue-and-fill.md](phenomena-catalogue-and-fill.md) | The master list of every phenomenon kind: observed rate, system odds, required contents, profile gate, scatter or fill, and the fill rules and audit | design (rules of 2026-10-15; items for TODO) |
| [multi-star-systems.md](multi-star-systems.md) | Multi-star and star-like systems at scatter time: the census, limits and probabilities, the stored member tree, motion, and the catalogue of about sixty observed and theoretical system types | design (plan of 2026-10-15; nothing built yet) |
| [disk-star-formation-history.md](disk-star-formation-history.md) | The thin disk's young, middle-aged and old star shares from its own gas history | built |
| [brown-dwarf-state.md](brown-dwarf-state.md) | A brown dwarf as a star: its L, T or Y type, temperature and radius from the cooling laws and the 20 pc census | built |
| [black-hole-formation-sources.md](black-hole-formation-sources.md) | The verified failed-supernova and direct-collapse shares and black hole remnant masses behind the black hole density (GEN.560, GEN.563) | research |
| [black-hole-density.md](black-hole-density.md) | The metallicity-dependent black hole density: the failed-supernova share f(Z) as two drawn anchors, the place factor, the rate sites, the tests and the seven open cards | plan |
| [unified-star-like-objects.md](unified-star-like-objects.md) | One recursive chain for anything that can be orbited: stars and star-like objects, planets, moons, rogue planets, free brown dwarfs and unbound comets, built through sector fill | plan |
| [unified-unit-quantity.md](unified-unit-quantity.md) | One registry of units, one Quantity value class and one declaration of every field's and column's quantity, copying the point-in-space class | plan |
| [unified-system-barycentre.md](unified-system-barycentre.md) | One System rule: a central object, members that can be systems themselves, the barycentre as the origin, derived offsets, and the hosting gate | plan |
| [mass-object-class.md](mass-object-class.md) | One master mass-object class: three family classes under it, state and parent as data, one core table plus three profile tables, a 10 m representability floor, the measured code savings | plan |
| [star-types-by-galactic-radius.md](star-types-by-galactic-radius.md) | How star types change with distance from the core | research (proposals for Boss) |
| [compact-remnant-regions.md](compact-remnant-regions.md) | Where neutron stars and black holes are observed | reference (feeds GEN.100) |

## Generation and reproducibility

| Note | What it covers | Status |
|---|---|---|
| [reproducible-galaxies.md](reproducible-galaxies.md) | Seeds, version keys and repair | partly built |
| [generation-determinism.md](generation-determinism.md) | What makes a seed rebuild the same galaxy; fingerprints and epochs | research |
| [object-ids.md](object-ids.md) | Interstellar object IDs | built |
| [object-id-options.md](object-id-options.md) | The options weighed for unique IDs and the decision | built |
| [interstellar-object-rates.md](interstellar-object-rates.md) | Rates of rogue planets, comets and other interstellar objects | built |
| [phenomenon-scatter-mass-cut.md](phenomenon-scatter-mass-cut.md) | The phenomenon scatter's lowest-mass cut | built |
| [rogue-planet-surface.md](rogue-planet-surface.md) | Rogue planet surface conditions | built |
| [comet-orbital-realism.md](comet-orbital-realism.md) | Elliptical, parabolic and interstellar comets | built |
| [nebula-and-asteroid-field-classes.md](nebula-and-asteroid-field-classes.md) | Nebula, remnant and asteroid field classes | built |
| [nebula-density-vs-reality.md](nebula-density-vs-reality.md) | Nebula density against the real Milky Way | partly built |
| [multistar-and-compact-systems.md](multistar-and-compact-systems.md) | Multi-star systems and compact-object primaries | research |
| [exotic-environments-planets-and-compact-binaries.md](exotic-environments-planets-and-compact-binaries.md) | Planets in nebulae, near pulsars, black holes and compact binaries | research |
| [anomalies.md](anomalies.md) | Anomalies for the starmap | research |
| [generation-performance-study.md](generation-performance-study.md) | Where generation time goes | research (recommendations built) |
| [generation-benchmark-report.md](generation-benchmark-report.md) | A measured benchmark of one build box | reference |
| [scatter-queue-feasibility.md](scatter-queue-feasibility.md) | A shuffled sector queue for scatters | research (recommendations built) |

## Habitability, species and population

| Note | What it covers | Status |
|---|---|---|
| [habitability-index.md](habitability-index.md) | The planetary habitability index (PHI-4) and the protection gear | partly built (scores built, GEN.84) |
| [atmospheres-retention-and-classes.md](atmospheres-retention-and-classes.md) | Atmosphere retention and class timing | research |
| [activity-magnetism-radiation-hydrosphere.md](activity-magnetism-radiation-hydrosphere.md) | Stellar activity, magnetism, radiation dose and hydrosphere | partly built |
| [population-and-politics.md](population-and-politics.md) | Species, civilizations, territories, tech levels and facilities | built |

## Orbits, navigation and courses

| Note | What it covers | Status |
|---|---|---|
| [orbital-updates.md](orbital-updates.md) | Moving things: positions, thresholds, gravity, collisions and spin | partly built |
| [orbital-solvers-and-integrators.md](orbital-solvers-and-integrators.md) | Kepler solvers, integrators, thresholds and the influence search | research |
| [galactic-potential.md](galactic-potential.md) | The galaxy's own gravity and galactic orbits | research |
| [collisions-and-mergers.md](collisions-and-mergers.md) | How often bodies collide and how to detect it | research |
| [navigation-frames.md](navigation-frames.md) | Nested navigation frames, bearing, mark and warp speed | built |
| [course-routing.md](course-routing.md) | Course routing and nearby search | partly built |
| [course-avoidance.md](course-avoidance.md) | Keep-out radii and path bending | research |
| [sky-view.md](sky-view.md) | The view from a planet | research |
| [gravity-map.md](gravity-map.md) | The field of a sector as a heat map | plan (the 9.1 headline) |

## Maps and the web interface

| Note | What it covers | Status |
|---|---|---|
| [galaxy-drilldown-navigation.md](galaxy-drilldown-navigation.md) | The Galaxy Map's drill-down | built (partly replaced by the fly-through) |
| [drilldown-region-sizes.md](drilldown-region-sizes.md) | Step sizes and region shapes for the drill-down | built (MAP.151) |
| [edit-everything.md](edit-everything.md) | The 8.6 headline: every admin change is a plan the admin approves; aspects, locks, regions, the rearranger | plan |
| [roles-and-requests.md](roles-and-requests.md) | The 10.0 headline: user, moderator and admin roles; edit requests are change plans that the role above approves | plan |
| [backups-and-repair.md](backups-and-repair.md) | One bz2 archive of the settings, the control database and the parity data, and restore and repair from it (settings 8.1, chain 8.9) | plan |
| [fly-through-view-distance.md](fly-through-view-distance.md) | The fly-through map: view distance, near field, free camera | partly built |
| [zoom-star-visibility.md](zoom-star-visibility.md) | Smoother star visibility while zooming | research (superseded in part by MAP.152) |
| [galaxy-map-wire-format.md](galaxy-map-wire-format.md) | What the browser downloads for the map | built |
| [nebula-map-visibility.md](nebula-map-visibility.md) | Making nebulae visible on the Galaxy Map | partly built |
| [map-ui-and-frontend-libraries.md](map-ui-and-frontend-libraries.md) | The map interface and front-end libraries | research |
| [3d-rendering-camera-research.md](3d-rendering-camera-research.md) | Camera proximity fading in dense 3D scenes | research (background) |
| [system-description.md](system-description.md) | One model of a system's text, many renderers | plan |
| [units-and-number-formatting.md](units-and-number-formatting.md) | The unit ladder and number format | plan (UX.23 open) |

## Server, queue and operations

| Note | What it covers | Status |
|---|---|---|
| [command-inventory.md](command-inventory.md) | Every command does one and only one thing and the order lives only in the wrapper's step list: Boss's principle (2026-10-17), and the inventory of every command and entry point with the item that reshapes each one (PERF.118, PERF.120 to PERF.128, OPS.60, OPS.61) | plan |
| [redis-and-the-work-queue.md](redis-and-the-work-queue.md) | How Redis is used and what to change (PERF.86) | research (decided; being built) |
| [work-queue-audit.md](work-queue-audit.md) | What starts work and what moves to the queue | built (replaced in part by PERF.86) |
| [performance-eta-queue-and-caching.md](performance-eta-queue-and-caching.md) | Progress estimates, the job queue and caching | research (queue under redo) |
| [slow-reads-and-timeouts.md](slow-reads-and-timeouts.md) | Keeping slow database reads inside the time limit | mostly built |
| [db-check-and-parity-repair.md](db-check-and-parity-repair.md) | Database check, parity repair and long migrations | partly built |
| [ops-scheduling-and-rotation.md](ops-scheduling-and-rotation.md) | Apache reload, daily maintenance and locks | partly built |
| [api-design-standards.md](api-design-standards.md) | API versions, key scopes, uploads and logging | research |
| [settings-seo-and-accounts.md](settings-seo-and-accounts.md) | Settings, the configuration page, SEO and accounts | research |
| [login-brute-force-protection.md](login-brute-force-protection.md) | Login lockouts and two-factor sign-in | built |
| [test-suite-rebuild.md](test-suite-rebuild.md) | The faster, modular test suite | partly built (stages 1 to 4) |

## Imported research

Boss's own research documents, kept as received. Each starts with a banner linking the notes above that use it. They say what Boss researched, not what planetGen does.

| File | Used by |
|---|---|
| [Atmospheric Toxicity](Atmospheric%20Toxicity.md), [Chemical Habitability](Chemical%20Habitability.md), [Exoplanetary Habitability](Exoplanetary%20Habitability.md), [Naturally Occuring Ionizing Radiation](Naturally%20Occuring%20Ionizing%20Radiation.md) | habitability and atmosphere notes |
| [Planetary Habitability Index](Planetary%20Habitability%20Index.md), [Mathematical and Algorithmic Implementation of the Planetary Habitability Index](Mathematical%20and%20Algorithmic%20Implementation%20of%20the%20Planetary%20Habitability%20Index.md), [Planetary Habitability and Speculative Xenobiology](Planetary%20Habitability%20and%20Speculative%20Xenobiology.md) | [habitability-index.md](habitability-index.md) |
| [Speculative Xenobiology Examples](Speculative%20Xenobiology%20Examples.md), [Speculative Xenobiology Extremes](Speculative%20Xenobiology%20Extremes.md) | [habitability-index.md](habitability-index.md) |
| [Environmental.Extremes.and.Protection.Research](Environmental.Extremes.and.Protection.Research.md), [Thermal.Respiratory.Protection.Research](Thermal.Respiratory.Protection.Research.md) | [habitability-index.md](habitability-index.md) |
| [Computational Astrodynamics](Computational%20Astrodynamics.md), [Observational Kinetics for Rotational Vectors](Observational%20Kinetics%20for%20Rotational%20Vectors.md) | [orbital-updates.md](orbital-updates.md) and related |
| [Orbital Position and Vector Mathematical Foundations](Orbital%20Position%20and%20Vector%20Mathematical%20Foundations.md), [Orbital Position and Vector Update Algorithms](Orbital%20Position%20and%20Vector%20Update%20Algorithms.md), [Orbital Update Collision Detection](Orbital%20Update%20Collision%20Detection.md), [Orbital Update Full Algorithm](Orbital%20Update%20Full%20Algorithm.md) | [orbital-updates.md](orbital-updates.md) |
| [Star Trek Anomalies and Science](Star%20Trek%20Anomalies%20and%20Science.md), [Star Trek Anomaly Probabilities](Star%20Trek%20Anomaly%20Probabilities.md) | [anomalies.md](anomalies.md) |
| [Technological Assessment](Technological%20Assessment.md) | [population-and-politics.md](population-and-politics.md) |
| [Library Migration Workflow](Library%20Migration%20Workflow.md), [Web UX Development Notes](Web%20UX%20Development%20Notes.md), [Web UX and Job Management Guide](Web%20UX%20and%20Job%20Management%20Guide.md) | [library-migration.md](library-migration.md) and related |

## Archive

[archive/](archive/) holds notes about designs that were replaced ([galaxy-coordinate-system-shells.md](archive/galaxy-coordinate-system-shells.md), [galaxy-disk-density-rev2.md](archive/galaxy-disk-density-rev2.md)). They are history.
