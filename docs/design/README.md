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
| [design-decisions.md](design-decisions.md) | The choices that shape the project, each with its status | reference |
| [todo-number-map.md](todo-number-map.md) | Old TODO numbers to category IDs, and the next free ID of each category | reference |
| [library-migration.md](library-migration.md) | Moving to third-party libraries and into packages | built (queue being redone, PERF.86) |
| [math-and-library-efficiency.md](math-and-library-efficiency.md) | Where the heavy math runs; Python 3.12 and numba | partly built (PERF.88 open) |

## Galaxy, sectors and coordinates

| Note | What it covers | Status |
|---|---|---|
| [galaxy-coordinate-system.md](galaxy-coordinate-system.md) | The spherical shell, band and wedge sector grid | built (schema v91); the cylinder it replaced is in [archive/](archive/galaxy-coordinate-system-cylinder.md) |
| [spherical-sector-coordinates.md](spherical-sector-coordinates.md) | The plan for spherical shells, bands and wedges | plan, built (the switch merged to main on 2026-10-15, #1607; the Galaxy Map ladder is MAP.205 and MAP.206) |
| [spherical-shell-sectors.md](spherical-shell-sectors.md) | What spherical-shell sectors would mean (the question behind the plan) | research |
| [spherical-ladder-meridians.md](spherical-ladder-meridians.md) | Meridians for the spherical block ladder | research |
| [spherical-block-ladder.md](spherical-block-ladder.md) | The Galaxy Map's block ladder (243/81/27/9/3) and windows on shells | built (MAP.205, MAP.206) |
| [galaxy-disk-density.md](galaxy-disk-density.md) | The galaxy density model | built |
| [fill-order-curves-and-core.md](fill-order-curves-and-core.md) | Fill order, region selection and the galactic core | research |
| [sampling-backfill-and-resume.md](sampling-backfill-and-resume.md) | Sampling, backfill, void filling, resumable runs and directives | research (resume on hold, PERF.86) |
| [mass-backfill.md](mass-backfill.md) | Four rings of heavier stars around generated sectors | built |
| [multiple-galaxies.md](multiple-galaxies.md) | Neighbour galaxies and a second galaxy | research |
| [globular-clusters.md](globular-clusters.md) | Globular clusters in generation | research (no model exists) |
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
| [gravity-map.md](gravity-map.md) | The field of a sector as a heat map | plan (the 9.0 headline) |

## Maps and the web interface

| Note | What it covers | Status |
|---|---|---|
| [galaxy-drilldown-navigation.md](galaxy-drilldown-navigation.md) | The Galaxy Map's drill-down | built (partly replaced by the fly-through) |
| [drilldown-region-sizes.md](drilldown-region-sizes.md) | Step sizes and region shapes for the drill-down | built (MAP.151) |
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
