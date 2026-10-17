# Version 9.0: Composition and atmosphere colour

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 9.0

## Goal

What every body, star and atmosphere is made of, worked out from real science, and what follows from it: the colour of the air and the radiation in the air and on the ground. Boss (2026-10-16, 11:25 PM EDT and 11:58 PM EDT): atmosphere composition and the composition of planets, moons, stars and everything is the joint headline of 9.0 with the atmosphere colour; every body has figures on its composition based on actual science, not made up where it can be avoided, and drawing the density first and filling in a fitting composition is as acceptable as deriving it. The atmosphere colour (sky colour and atmosphere colour are the same thing) is the absolute colour, not adapted to the star, and covers gas giants; it shows in the descriptions and replaces the planet colour on the maps. The atmospheric radiation (band transmission and the greenhouse effect from the mixture) and the surface radiation (the dose from the air's shielding and the ground's own composition) are derived the same way. Everything else is built before: the research, the designs and the data tables (the elements, the gases, the materials, the stars) are in 8.1 phases 1 and 2, the star composition function in 8.1 phase 3, and the generation and storage in 8.4 phases 1 to 4, so 9.0 holds the colour, the radiation and what shows them. A major version needs Boss's word. The master note is [composition-model.md](../design/composition-model.md); the atmosphere plan is [atmosphere-colour.md](../design/atmosphere-colour.md), based on [sky-colour-analysis.md](../design/sky-colour-analysis.md). The view from a planet and the gravity map are 9.1, and the backups and repair chain is 8.9 (Boss, 11:26 PM EDT).

## Phase 1: The colour model

The physics of the colour: the blackbody colour of the star, the Rayleigh optical depth of the mixture, the colour of a solid world's air and of a giant's cloud top, all as the absolute colour (not adapted to the star).

### Atmosphere colour

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.416 | Blackbody colour of a star: Planck to the CIE curves to sRGB, as the first part of the photometry library |  | Boss 11:24 PM EDT; VIEW.11 extends it. |
| GEN.417 | Rayleigh optical depth of the mixture: per-gas cross sections from the table and the air column from pressure and gravity | GEN.410, GEN.411 | Boss 11:24 PM EDT; pure function of the stored mixture. |
| GEN.418 | Atmosphere colour of a solid world: the absolute colour of the air with ozone, haze, dust, cloud and the ground's bounce | GEN.416, GEN.417, GEN.413, GEN.414, DB.67, DB.68, GEN.441, GEN.462 | Boss 11:24 PM EDT; absolute colour, not eye-adapted. |
| GEN.419 | Atmosphere colour of gas and ice giants, read at a reference level in the cloud tops | GEN.418 | Boss 11:24 PM EDT; reference level chosen in the design. |

## Phase 2: Atmospheric and surface radiation

The radiation derived from the stored composition: the band transmission and greenhouse of the atmosphere, the dose at the ground from the air's shielding and the ground's own content, and the habitability values that read them.

### Atmospheric and surface radiation

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.443 | Atmospheric radiation: transmission by band and the greenhouse effect from the stored mixture, the aerosols and the cloud | GEN.417, GEN.413, GEN.414, GEN.436 | Boss 11:58 PM EDT; band transmission and greenhouse from the mixture. |
| GEN.453 | The radiation constants are corrected: the dipole cutoff scales with moment over radius squared, the dose curve against the column and the field factor follow the measured and simulated tables (bug) | GEN.436 | Atmosphere research lane finding 4 (part 07 sections 2.1 and 4); generator epoch bump. Core 4 proposed. |
| GEN.444 | Surface radiation: the ground-level dose from the air's shielding and the ground's own radiogenic content | GEN.443, GEN.441, GEN.438, GEN.453, GEN.459, GEN.478 | Boss 11:58 PM EDT; replaces the proxies of physics/radiation.py; generator epoch bump. |
| GEN.445 | Habitability, toxicity and the UV index read the derived radiation instead of the proxies | GEN.444 | Boss 11:58 PM EDT; habitability reads the derived values. |

## Phase 3: Composition and colour in the product

The colour shown in the descriptions, on the maps in place of the planet colour, the composition and radiation text of stars, surfaces and small bodies, and the documents of the colour.

### Descriptions and maps

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.420 | The descriptions show the atmosphere colour: a swatch and a word in the planet and moon facts, in every output | GEN.418, GEN.419 | Boss 11:24 PM EDT; one place serves every output. |
| MAP.244 | The maps and pages use the atmosphere colour where they showed the planet colour | GEN.420 | Boss 11:24 PM EDT; replaces the temperature gesture. |
| GEN.446 | The descriptions state the composition of stars, surfaces and small bodies, and the radiation at the surface | GEN.440, GEN.441, GEN.442, GEN.444, GEN.445, GEN.221 | Boss 11:58 PM EDT; stars, surfaces, small bodies and radiation in the text. |
| DOC.64 | Docs: the atmosphere colour, its method and its sources | GEN.420, MAP.244 | Boss 11:24 PM EDT; docs only. |

### Stand-in replacements

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.473 | Replace the asteroid taxon shares and the zone kernel, measured in the Solar System belt, by a derived or general distribution | GEN.442 | Boss 1:44 AM EDT, Stand-in only: replace the measured Solar System taxon shares by a derived or general distribution. |
| GEN.474 | Replace the Earth rock-type library, measured on Earth, by ground types derived from the bulk composition and the melting of the mantle | GEN.441 | Boss 1:44 AM EDT, Stand-in only: replace the Earth rock-type library by derived ground types. |
| GEN.475 | Meteorite groups become outputs and not inputs: regions of composition space reachable from the drawn bulk, with the measured groups kept as calibration | GEN.473, GEN.474, GEN.451 | Boss 1:44 AM EDT, Stand-in only: R3, meteorite groups as outputs of the drawn chemistry. |
| GEN.476 | Aluminium-26 and iron-60 follow from the star's formation environment, not from the Solar System's own value and range | GEN.459, GEN.440 | Boss 1:44 AM EDT, Stand-in only: R4, Al-26 and Fe-60 from the formation environment. |

### Weight research (atmosphere)

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.480 | Source the occurrence weights of bare rock, hydrogen-helium envelope and water-rich worlds and of the envelope metal fraction, by radius, mass and insolation (AR1) | GEN.478 | Atmosphere lane research need R1 (rows b4, e3, e4, e5); replaces the seeded bracket of the DESIGN rows. |
| GEN.481 | Source the volatile delivery against formation distance: the water mass fraction by a_form and the snow-line law over time (AR2) | GEN.478 | Atmosphere lane research need R2 (rows b1 to b6); replaces the seeded bracket of the DESIGN rows. |
| GEN.482 | Source the magma-ocean fO2 prior against mass and the oxidation law of the later mantle (AR3) | GEN.478 | Atmosphere lane research need R3 (rows a1, a2, a3); replaces the seeded bracket of the DESIGN rows. |
| GEN.483 | Source the cosmic-shoreline constants and the pressure caps by retention regime (AR4) | GEN.478 | Atmosphere lane research need R4 (rows d6); replaces the seeded bracket of the DESIGN rows. |
| GEN.484 | Source the giant planets' C/O distribution and the fraction of carbon-rich hosts (AR5) | GEN.478 | Atmosphere lane research need R5 (rows f7); replaces the seeded bracket of the DESIGN rows. |
| GEN.485 | Source the hot-Jupiter thermochemistry above 2200 K, a Kzz distribution and the edges of the 1500 K and 2200 K regimes (AR6) | GEN.478 | Atmosphere lane research need R6 (rows h2, h3, h4); replaces the seeded bracket of the DESIGN rows. |
| GEN.486 | Source the lava-world vapour regimes and the fraction of lava worlds that keep volatiles (AR7) | GEN.478 | Atmosphere lane research need R7 (rows i6); replaces the seeded bracket of the DESIGN rows. |
| GEN.487 | Source the evolved and frost atmospheres: cold vapour pressures, the habitable-zone coefficients, the occurrence of Titan-, Pluto- and Venus-like bodies and the abiotic O2 frequency (AR8) | GEN.478 | Atmosphere lane research need R8 (rows j1 to j5); replaces the seeded bracket of the DESIGN rows. |
| GEN.488 | Source the comet coma classes and their volatile ratios (AR9) | GEN.478 | Atmosphere lane research need R9 (rows k1 to k5, k7); replaces the seeded bracket of the DESIGN rows. |
| GEN.489 | Source the white-dwarf and pulsar planet options (AR10) | GEN.478 | Atmosphere lane research need R10 (rows a4, l1, l2, l3); replaces the seeded bracket of the DESIGN rows. |
| GEN.490 | Source the trace-gas presence weights (AR11) | GEN.478 | Atmosphere lane research need R11 (rows m2 to m6, m11); replaces the seeded bracket of the DESIGN rows. |
| GEN.491 | Source the cloud and haze windows and yields (AR12) | GEN.478 | Atmosphere lane research need R12 (rows n6, n7); replaces the seeded bracket of the DESIGN rows. |
| GEN.492 | Source the ground radioisotope abundances by host and the potassium depletion by body class (AR13) | GEN.478 | Atmosphere lane research need R13 (rows o2, o3, o4, o6); replaces the seeded bracket of the DESIGN rows. |
| GEN.493 | Source the initial-rotation distribution and the saturation time against stellar mass for the host-star UV regime (AR14) | GEN.478 | Atmosphere lane research need R14 (rows p6); replaces the seeded bracket of the DESIGN rows. |
| GEN.494 | Source the brown-dwarf type against temperature and the T/Y boundary near 500 to 700 K (AR15) | GEN.478 | Atmosphere lane research need R15 (rows g1, g4); replaces the seeded bracket of the DESIGN rows. |
| GEN.501 | Source the activity factor of the NH3 photolysis rate from the Huebner and Mukherjee rate tables |  | Body lane reconciliation: NH3 lifetime 6.7e3 s quiet Sun, 5.6e3 s quiet/active mean; activity factor unknown. |

### Tuned draws and selection models

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.503 | The shared tuning harness: a target table, a simulate-and-compare function and a loader that refuses a tunable row without a target | GEN.472, GEN.478 | Tunable gaps report: The shared tuning harness: a target table, a simulate-and-compare func. |
| GEN.504 | The Kepler transit-plus-recovery selection model | GEN.503 | Tunable gaps report: The Kepler transit-plus-recovery selection model. |
| GEN.505 | The mass-radius sample selection model | GEN.503 | Tunable gaps report: The mass-radius sample selection model. |
| GEN.506 | The JWST atmosphere target-selection model | GEN.503 | Tunable gaps report: The JWST atmosphere target-selection model. |
| GEN.507 | The comet admission model by heliocentric distance | GEN.503 | Tunable gaps report: The comet admission model by heliocentric distance. |
| GEN.508 | The white-dwarf pollution detection and settling model | GEN.503 | Tunable gaps report: The white-dwarf pollution detection and settling model. |
| GEN.509 | Tuned draw builder B1: the radius valley and the envelope | GEN.503, GEN.504 | Tunable gaps report: Tuned draw builder B1: the radius valley and the envelope. |
| GEN.510 | Tuned draw builder B2: delivery and water worlds | GEN.503, GEN.505 | Tunable gaps report: Tuned draw builder B2: delivery and water worlds. |
| GEN.511 | Tuned draw builder B3: retention and lava worlds | GEN.503, GEN.506 | Tunable gaps report: Tuned draw builder B3: retention and lava worlds. |
| GEN.512 | Tuned draw builder B4: the iron branch | GEN.503, GEN.505 | Tunable gaps report: Tuned draw builder B4: the iron branch. |
| GEN.513 | Tuned draw builder B5: redox and white-dwarf debris | GEN.503, GEN.508 | Tunable gaps report: Tuned draw builder B5: redox and white-dwarf debris. |
| GEN.514 | Tuned draw builder B6: radioisotopes | GEN.503 | Tunable gaps report: Tuned draw builder B6: radioisotopes. |
| GEN.515 | Tuned draw builder B7: ground and petrology | GEN.503 | Tunable gaps report: Tuned draw builder B7: ground and petrology. |
| GEN.516 | Tuned draw builder B8: giant planets | GEN.503, GEN.506 | Tunable gaps report: Tuned draw builder B8: giant planets. |
| GEN.517 | Tuned draw builder B9: brown dwarfs | GEN.503 | Tunable gaps report: Tuned draw builder B9: brown dwarfs. |
| GEN.518 | Tuned draw builder B10: the chemistry module | GEN.503 | Tunable gaps report: Tuned draw builder B10: the chemistry module. |
| GEN.519 | Tuned draw builder B11: aerosol amplitude | GEN.503, GEN.506 | Tunable gaps report: Tuned draw builder B11: aerosol amplitude. |
| GEN.520 | Tuned draw builder B12: the comet sampler | GEN.503, GEN.507 | Tunable gaps report: Tuned draw builder B12: the comet sampler. |
| GEN.521 | Tuned draw builder B13: belts | GEN.503 | Tunable gaps report: Tuned draw builder B13: belts. |
| GEN.522 | Tuned draw builder B14: moons and Kuiper belt objects | GEN.503 | Tunable gaps report: Tuned draw builder B14: moons and Kuiper belt objects. |
| GEN.523 | Tuned draw builder B15: host-star activity | GEN.503 | Tunable gaps report: Tuned draw builder B15: host-star activity. |
| GEN.524 | Tuned draw builder B16: halogens and trace gases | GEN.503 | Tunable gaps report: Tuned draw builder B16: halogens and trace gases. |
| GEN.525 | Tuned draw builder B17: the design-flag registry | GEN.503 | Tunable gaps report: Tuned draw builder B17: the design-flag registry. |

## Phase 4: The proofs and the headline

The radiation sweep, the documents of the whole composition model and the headline.

### The headline

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.212 | A radiation sweep over every class: transmission, greenhouse, shielding and dose agree with the composition | GEN.443, GEN.444, GEN.445, GEN.446 | Boss 11:58 PM EDT; the radiation proof. |
| DOC.66 | Docs: the composition model of stars, bodies and atmospheres, the radiation it derives and the sources of every number | TEST.212, TEST.211, DOC.63, DOC.65 | Boss 11:58 PM EDT; docs only. |
| GEN.421 | Headline of 9.0: what every body and atmosphere is made of, and the atmosphere colour that follows from it | GEN.420, MAP.244, DOC.64, GEN.446, TEST.209, TEST.210, TEST.211, TEST.212, DOC.66 | Headline of 9.0 (Boss, 11:25 PM EDT). |
