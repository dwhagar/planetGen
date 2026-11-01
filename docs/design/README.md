# Design notes

Every note under `docs/design/`, with what it is about and where it stands. Notes are grouped in one folder per topic, and this page has one section per folder. The status words:

- **built**: the design is in the code.
- **partly built**: some of it is; the note's own status line says which.
- **plan**: decided, not built yet.
- **research**: findings and recommendations; read the note's status line for what was taken up.
- **reference**: facts and measurements, not a design.
- **history**: describes a design that was replaced; kept for the reasons.

A note's own status line, at its top, is the authority; this page only points at it. A test (`src/tests/test_docs_links.py`) fails when a note in a topic folder is missing from this page. A new note goes in the folder of its topic, with a row here. The architecture of the whole program is [architecture.md](architecture/architecture.md); the decisions behind it are in [design-decisions.md](architecture/design-decisions.md).

**Standing rule for downloaded data (Boss, 12:47 AM EDT on 2026-10-18):** Downloads of static files take place at setup (the step `update.sh` runs, into a cache folder outside the repository) unless the data needs to be constantly updated; any code that fetches data at run time must say in its design why the data needs constant updates (Boss, 12:47 AM EDT on 2026-10-18, project chat, as the coordinator relayed). The shared step is GEN.548; the rule and its offline test design are in [gap-filling-plan.md](research/gap-filling-plan.md) section 6.5.

## The whole program

`architecture/`: how the program fits together, and the shared data model. The generated code map and the number map are linked from here; they live in [generated/](../generated/) and [plan/](../plan/).

| Note | What it covers | Status |
|---|---|---|
| [architecture.md](architecture/architecture.md) | Every file and the main flows | reference (the code map is being rewritten, DOC.28) |
| [architecture-map.html](../generated/architecture-map.html) | The generated code map: areas, modules, functions and calls, with the main flows traced from the code (DOC.67) | generated page (open it in a browser; `scripts/build_architecture_map.py`) |
| [design-decisions.md](architecture/design-decisions.md) | The choices that shape the project, each with its status | reference |
| [todo-number-map.md](../plan/todo-number-map.md) | Old TODO numbers to category IDs, and the next free ID of each category | reference |
| [library-migration.md](architecture/library-migration.md) | Moving to third-party libraries and into packages | built (queue being redone, PERF.86) |
| [math-and-library-efficiency.md](architecture/math-and-library-efficiency.md) | Where the heavy math runs; Python 3.12 and numba | partly built (PERF.88 open) |
| [numeric-precision-audit.md](architecture/numeric-precision-audit.md) | GEN.372: whether 64-bit numbers are enough for positions at galactic scale (yes, if no code subtracts two large absolute coordinates), the numbers behind it and every precision concern found | research (Core 1, 2026-10-17; nothing built) |
| [unified-star-like-objects.md](architecture/unified-star-like-objects.md) | One recursive chain for anything that can be orbited: stars and star-like objects, planets, moons, rogue planets, free brown dwarfs and unbound comets, built through sector fill | plan |
| [unified-unit-quantity.md](architecture/unified-unit-quantity.md) | One registry of units, one Quantity value class and one declaration of every field's and column's quantity, copying the point-in-space class | plan |
| [unified-system-barycentre.md](architecture/unified-system-barycentre.md) | One System rule: a central object, members that can be systems themselves, the barycentre as the origin, derived offsets, and the hosting gate | plan |
| [mass-object-class.md](architecture/mass-object-class.md) | One master mass-object class: three family classes under it, state and parent as data, one core table plus three profile tables, a 10 m representability floor, the measured code savings | plan |
| [object-ids.md](architecture/object-ids.md) | Interstellar object IDs | built |
| [object-id-options.md](architecture/object-id-options.md) | The options weighed for unique IDs and the decision | built |
| [designation-history.md](architecture/designation-history.md) | Names and IDs may change; a retired ID is reserved forever and the history is kept (GEN.1121) | designed |
| [system-description.md](architecture/system-description.md) | One model of a system's text, many renderers | plan |
| [units-and-number-formatting.md](architecture/units-and-number-formatting.md) | The unit ladder and number format | plan (UX.23 open) |
| [command-inventory.md](architecture/command-inventory.md) | Every command does one and only one thing and the order lives only in the wrapper's step list: Boss's principle (2026-10-17), and the inventory of every command and entry point with the item that reshapes each one (PERF.118, PERF.120 to PERF.128, OPS.60, OPS.61) | plan |

## Galaxy, sectors and coordinates

`galaxy/`: the shape and contents of the galaxy, the spherical sector grid and the laws behind them.

| Note | What it covers | Status |
|---|---|---|
| [galaxy-coordinate-system.md](galaxy/galaxy-coordinate-system.md) | The spherical shell, band and wedge sector grid | built (schema v91); the cylinder it replaced is in [archive/](../archive/galaxy-coordinate-system-cylinder.md) |
| [galaxy-types.md](galaxy/galaxy-types.md) | The 9.4 headline: one family per known galaxy type (spiral, lenticular, elliptical, irregular, dwarf), drawn from the seed through scaling relations, and the chain that leads to it | plan |
| [superbubble-density.md](galaxy/superbubble-density.md) | Superbubbles shape the stars: the seeded bubble field becomes a density modulation on the young and intermediate populations, with bounds, sector counts and positions, and the chain that builds it | plan |
| [spherical-sector-coordinates.md](galaxy/spherical-sector-coordinates.md) | The plan for spherical shells, bands and wedges | plan, built (the switch merged to main on 2026-10-15, #1607; the Galaxy Map ladder is MAP.205 and MAP.206) |
| [spherical-shell-sectors.md](galaxy/spherical-shell-sectors.md) | What spherical-shell sectors would mean (the question behind the plan) | research |
| [spherical-ladder-meridians.md](galaxy/spherical-ladder-meridians.md) | Meridians for the spherical block ladder | research |
| [spherical-block-ladder.md](galaxy/spherical-block-ladder.md) | The Galaxy Map's block ladder (243/81/27/9/3) and windows on shells | built (MAP.205, MAP.206) |
| [galaxy-disk-density.md](galaxy/galaxy-disk-density.md) | The galaxy density model | built |
| [fill-order-curves-and-core.md](galaxy/fill-order-curves-and-core.md) | Fill order, region selection and the galactic core | research |
| [realistic-galactic-core.md](galaxy/realistic-galactic-core.md) | A realistic galactic core (GEN.279, reopens GEN.155): the nuclear star cluster and disc as density terms, the resolution cap, tiered mass and luminosity limits, planets and orbits at the centre | recommendation (Core 3, 2026-10-16; defaults stand, nothing built yet) |
| [named-regions.md](galaxy/named-regions.md) | Named regions: the 9-block cells of the block ladder with a name, anchor, framing sphere and extent, their stored table, and the uniqueness audit | built (GEN.825, GEN.827, DB.107, API.35); map, pages and search not built |
| [multiple-galaxies.md](galaxy/multiple-galaxies.md) | Neighbour galaxies and a second galaxy | research |
| [globular-clusters.md](galaxy/globular-clusters.md) | Globular clusters in generation | partly built (GEN.159: cluster table, King tables, generator; GEN.160: cluster density in the sector gate) |
| [intergalactic-and-halo-contents.md](galaxy/intergalactic-and-halo-contents.md) | What fills the sphere off the disk: the stellar halo, satellites, ejected stars, gas, and how to give every cell contents | research (reach approved; build plan for TODO) |
| [galaxy-laws-research.md](galaxy/galaxy-laws-research.md) | The sourced ranges and laws for the stellar halo, the population ages and dispersions, the thick disk, the dwarf and satellite sizes and the cluster system (GEN.406), with the fifteen open requests | research |
| [gen533-bar-sources.md](galaxy/gen533-bar-sources.md) | The sources read for the bar and the unbarred bulge (GEN.533): bar fraction against mass and stage, bar and bulge axis ratios, bar angle, with the conflicts and the rule that picked each draw | research |
| [spherical-region-queries.md](galaxy/spherical-region-queries.md) | Selecting a sphere of the galaxy with a database query: the cell key, key runs, a covering index, measured against the box, the spatial index and the sector walk | research (2026-10-15; nothing built) |
| [post-main-sequence-stars.md](galaxy/post-main-sequence-stars.md) | Stars after the main sequence (GEN.935): the MIST tracks rebuilt per mass at setup, the windows the generator's phase times place on them, Cepheids, the bright-star tables, what was removed | built (Star 1, 2026-10-29), part of the star mass function epoch (GEN.932) |
| [stellar-mass-function.md](galaxy/stellar-mass-function.md) | The stellar mass function (GEN.927 to GEN.931): the single-star curve, the derived primary curve, the white dwarf, neutron star and black hole remnant relations and what each consumer reads | built |
| [disk-star-formation-history.md](galaxy/disk-star-formation-history.md) | The thin disk's young, middle-aged and old star shares from its own gas history | built |
| [star-types-by-galactic-radius.md](galaxy/star-types-by-galactic-radius.md) | How star types change with distance from the core | research (proposals for Boss) |
| [galactic-potential.md](galaxy/galactic-potential.md) | The galaxy's own gravity and galactic orbits | research |

## Stars, remnants, phenomena and small bodies

`systems/`: what lives in the galaxy besides the disk and halo themselves.

| Note | What it covers | Status |
|---|---|---|
| [later-phenomena-scatter-rules.md](systems/later-phenomena-scatter-rules.md) | Which later phenomena are scattered, backfilled or fill-only (goal issue #1677): what each item stores, the magnetar scatter class, what the rule set and the profile gain, the columns for Core 1, and the changes to the TODO items | recommendation (Core 3, 2026-10-16; defaults stand, nothing built yet) |
| [phenomena-catalogue-and-fill.md](systems/phenomena-catalogue-and-fill.md) | The master list of every phenomenon kind: observed rate, system odds, required contents, profile gate, scatter or fill, and the fill rules and audit | design (rules of 2026-10-15; items for TODO) |
| [multi-star-systems.md](systems/multi-star-systems.md) | Multi-star and star-like systems at scatter time: the census, limits and probabilities, the stored member tree, motion, and the catalogue of about sixty observed and theoretical system types | design (plan of 2026-10-15; nothing built yet) |
| [overlap-free-systems.md](systems/overlap-free-systems.md) | Systems that cannot overlap by construction (GEN.1006): one frozen record per star that the tree, the build and the planets read, orbits drawn from allowed ranges, a clearance ledger for planets, moons and belts, an independent verifier, the redraw-rate test, and the build in ten parts | plan (checked against the papers; nothing built) |
| [lagrange-trojans.md](systems/lagrange-trojans.md) | What may sit at the L4 and L5 points of a pair (GEN.1149 to GEN.1151, MAP.354): the mass rule, libration, the swarm split, the map label, and the build in four parts | design (pure module built; wiring is GEN.1150) |
| [belt-evolution.md](systems/belt-evolution.md) | Asteroid belts that evolve: what formed, collisional grinding, giant clearing and resonance gaps, drift, star events, and a ledger of the individual bodies chosen by a mass-cut rule (every body above 1e19 kg, up to 90 percent of the belt mass, at most 300 per belt), with status words (planet, dwarf planet, minor body) and a sub-numbered designation for each member | plan (Boss said yes 4:27 PM EDT 2026-10-30; nothing built) |
| [belt-origins.md](systems/belt-origins.md) | Where belts come from (gap, unfinished tail, snow-line bump, stirring, instability depletion, collisional debris), the origin recorded on each belt, a scatter phase, and how many Class D and dwarf bodies a system needs | plan (Boss asked 4:19 PM EDT 2026-10-30; nothing built) |
| [belt-sources.md](systems/belt-sources.md) | What was read for the belt items (GEN.1060): Yarkovsky drift, the chaotic zone, resonance loss, early depletion, Soter and Margot planet numbers, round masses, the Kuiper belt size law, white dwarf pollution, debris disc incidence; with the registry rows BR1 to BR18 | read; numbers registered (section 6.95) |
| [rogue-planet-drift.md](systems/rogue-planet-drift.md) | Where ejected rogue planets drift (about 1 kpc per 0.5 Gyr along the galactic orbit) and where each is stored: the closed-form drift of the donor model; the virtual-donor storage is superseded by region-object-store.md (real donors, stored where they land) | plan (research 5:20 PM EDT 2026-10-30; the drift function is built, the storage is replaced) |
| [region-object-store.md](systems/region-object-store.md) | One galaxy-wide object store generated in region work units (whole sectors, bands of the galaxy), sectors as a navigation table with a chart ledger, rogues from real donors stored where they land, the three-stage rogue spread with the None/Few/Normal/Many level, and designation history (GEN.1118) | plan (study final 5:44 PM EDT 2026-10-30; defaults for tiles and bands until Boss answers; items GEN.1123 to GEN.1129, TEST.501) |
| [interstellar-asteroid-fields.md](systems/interstellar-asteroid-fields.md) | Asteroid fields in interstellar space from rogue-planet collisions and the better sources (supernova-unbound belts): the rate (about one collision per Galaxy life), outcomes, dispersal as a pure function of age, composition from two parents, members and designations | plan (items filed 4:51 PM EDT 2026-10-30; rate with Boss; nothing built) |
| [captured-asteroid-moons.md](systems/captured-asteroid-moons.md) | Irregular moons of the giants as captured belt and disc bodies: Hill-sphere stable zones, capture windows (gas era, instability, exchange), counts by grinding, orbits and families, storage by the belt rule with a 1e16 kg floor | plan (Boss said yes 4:27 PM EDT 2026-10-30, decided 4:34 PM EDT; nothing built) |
| [weird-worlds.md](systems/weird-worlds.md) | Physically plausible oddities (Venus-type spin, tilt-driven seasons, iron-rich cores, resonant chains, disintegrating planets, super-puffs, nightside collapse, bistable climates, abiotic oxygen): each a cause at an input of the chain from a gated tail rate, a weirdness dial, an Odd line on the page, and a plausibility harness | plan (Boss said yes 4:15 PM EDT 2026-10-30; nothing built) |
| [quadruple-stability-3p1.md](systems/quadruple-stability-3p1.md) | The 3+1 quadruple stability network of [R4575] ported to a plain numpy forward pass (GEN.1029), its accuracy check against N-body runs, and how GEN.1016 uses it | built (module merged; nothing calls it until GEN.1016) |
| [quadruple-stability-2p2.md](systems/quadruple-stability-2p2.md) | The 2+2 quadruple stability network of [R4575] ported to a plain numpy forward pass (GEN.1031), its accuracy against N-body runs (91 percent against 82 for two nested triples), and how GEN.1016 uses it | built (module merged; nothing calls it until GEN.1016) |
| [brown-dwarf-state.md](systems/brown-dwarf-state.md) | A brown dwarf as a star: its L, T or Y type, temperature and radius from the cooling laws and the 20 pc census | built |
| [white-dwarf-state.md](systems/white-dwarf-state.md) | A white dwarf as a star: temperature, luminosity and radius from the Bedard et al. cooling sequences at its cooling age and hydrogen layer (GEN.742) | built |
| [remnant-age.md](systems/remnant-age.md) | A neutron star's or stellar black hole's age: the population's birth age less the progenitor's time to collapse | built |
| [black-hole-formation-sources.md](systems/black-hole-formation-sources.md) | The verified failed-supernova and direct-collapse shares and black hole remnant masses behind the black hole density (GEN.560, GEN.563) | research |
| [wanderer-moon-sources.md](systems/wanderer-moon-sources.md) | The sources behind the moons of wanderers: the survival of moons when their planet is ejected and the disks of isolated planetary-mass objects, with the stated gaps (GEN.560, GEN.558) | research |
| [black-hole-density.md](systems/black-hole-density.md) | The metallicity-dependent black hole density: the failed-supernova share f(Z) as two drawn anchors, the place factor, the rate sites, the tests and the seven open cards | built (GEN.563) |
| [neutron-star-density.md](systems/neutron-star-density.md) | The neutron star density derived from the same massive-star budget as the black holes: the place factor with the survivors of the failed share, the unbound share, the numbers against 4.0e-3 | built (GEN.554 slice 1) |
| [neutron-star-structure.md](systems/neutron-star-structure.md) | A neutron star's radius, mass and crust from a TOV integration over the 21 equation-of-state tables, against NICER (GEN.968) | built (lock-free step; wiring GEN.969) |
| [stellar-spectra-library.md](systems/stellar-spectra-library.md) | A star's light by band from ATLAS9 and PHOENIX-ACES model atmospheres plus the corona's far ultraviolet from MUSCLES, against the Sun and eleven M and K stars (GEN.970) | built and wired (GEN.970, GEN.971) |
| [alpha-and-age-metallicity.md](systems/alpha-and-age-metallicity.md) | A high-alpha star read off a solar-scaled track at a higher [Fe/H], and a disc star's [Fe/H] tied to its age by branch (GEN.744) | built (lock-free step; wiring GEN.744) |
| [variability.md](systems/variability.md) | The kind of variable a star is (classical Cepheid, Mira or semiregular, T Tauri, Herbig Ae/Be), a derived row of the star sheet (GEN.944) | built |
| [carbon-stars.md](systems/carbon-stars.md) | Carbon stars and S stars of the asymptotic giant branch from a thermal-pulse count and the third dredge-up, a derived row of the star sheet (GEN.1026) | built |
| [chemically-peculiar.md](systems/chemically-peculiar.md) | Ap and Am stars: the shares, the slow rotation and the field of the chemically peculiar branch of the A-star draw (GEN.1027 part 1) | module built, not wired |
| [stellar-state.md](systems/stellar-state.md) | The state word of a star and the age since it stopped burning, from its Yerkes class and phase (GEN.785, GEN.786) | built and wired (GEN.785, GEN.786) |
| [explosion-islands.md](systems/explosion-islands.md) | Which massive stars leave a black hole: the failure windows of initial mass that hold the failed share, the progenitor draw over them, the direct collapse | built (GEN.554 slice 2) |
| [natal-kick.md](systems/natal-kick.md) | The share of a galaxy's neutron stars its natal kicks eject: the kick law against the galaxy's escape speed, replacing the 20% to 40% draw | built (GEN.554 slice 3) |
| [compact-remnant-regions.md](systems/compact-remnant-regions.md) | Where neutron stars and black holes are observed | reference (feeds GEN.100) |
| [interstellar-object-rates.md](systems/interstellar-object-rates.md) | Rates of rogue planets, comets and other interstellar objects | built |
| [phenomenon-scatter-mass-cut.md](systems/phenomenon-scatter-mass-cut.md) | The phenomenon scatter's lowest-mass cut | built |
| [rogue-planet-surface.md](systems/rogue-planet-surface.md) | Rogue planet surface conditions | built |
| [comet-orbital-realism.md](systems/comet-orbital-realism.md) | Elliptical, parabolic and interstellar comets | built |
| [nebula-and-asteroid-field-classes.md](systems/nebula-and-asteroid-field-classes.md) | Nebula, remnant and asteroid field classes | built |
| [nebula-density-vs-reality.md](systems/nebula-density-vs-reality.md) | Nebula density against the real Milky Way | partly built |
| [multistar-and-compact-systems.md](systems/multistar-and-compact-systems.md) | Multi-star systems and compact-object primaries | research |
| [exotic-environments-planets-and-compact-binaries.md](systems/exotic-environments-planets-and-compact-binaries.md) | Planets in nebulae, near pulsars, black holes and compact binaries | research |
| [anomalies.md](systems/anomalies.md) | Anomalies for the starmap | research |

## Generation and reproducibility

`generation/`: seeds, determinism, fill order, backfill and generation speed.

| Note | What it covers | Status |
|---|---|---|
| [sampling-backfill-and-resume.md](generation/sampling-backfill-and-resume.md) | Sampling, backfill, void filling, resumable runs and directives | research (resume on hold, PERF.86) |
| [void-fill-census.md](generation/void-fill-census.md) | GEN.102 again on the spherical grid: the sparse cells, the halo beyond the outline, the recommendation | research |
| [placement-by-population.md](generation/placement-by-population.md) | GEN.103: where each kind is placed by population, and the three type-draw rules still to build | design |
| [mass-backfill.md](generation/mass-backfill.md) | Four rings of heavier stars around generated sectors | built |
| [reproducible-galaxies.md](generation/reproducible-galaxies.md) | Seeds, version keys and repair | partly built |
| [generation-determinism.md](generation/generation-determinism.md) | What makes a seed rebuild the same galaxy; fingerprints and epochs | research |
| [generation-performance-study.md](generation/generation-performance-study.md) | Where generation time goes | research (recommendations built) |
| [generation-performance.md](generation/generation-performance.md) | Every generation step timed at four galaxy sizes with `planetgen benchmark --suite`, the fixes of the first round and a ranked list of what is left (PERF.146) | reference |
| [generation-benchmark-report.md](generation/generation-benchmark-report.md) | A measured benchmark of one build box | reference |
| [scatter-queue-feasibility.md](generation/scatter-queue-feasibility.md) | A shuffled sector queue for scatters | research (recommendations built) |

## Composition

`composition/`: what stars, planets, small bodies and atmospheres are made of, the colour of the sky, and the research parts behind them.

| Note | What it covers | Status |
|---|---|---|
| [composition-model.md](composition/composition-model.md) | The master note of the 9.0 composition headline: what every body, star and atmosphere is made of, the chain from element to derived value, the research lanes and where each item sits | checked design, for Boss's approval (GEN.437); the loader, element table, stellar relations and star composition are built (GEN.472, GEN.438 to GEN.440) |
| [surface-liquid-composition.md](composition/surface-liquid-composition.md) | What a body stores about its surface liquid: the acid, ammonia, salt and hydrocarbon makeup, a second liquid, the liquids and ices still missing, and how the life classes read them (GEN.1217) | draft, for Boss's answers to six questions (Tides 1, 2026-10-31); nothing built |
| [atmosphere-colour.md](composition/atmosphere-colour.md) | Atmosphere colour, the 9.0 feature: one gas and compound table, the mixture that sets the molar density, stored albedo and greenhouse values, ozone, haze, dust and cloud, then the absolute colour of the air, with the chain that builds it | approved 2026-10-17 (9:40 to 9:55 PM EDT) (GEN.409) |
| [bulk-composition.md](composition/bulk-composition.md) | The bulk composition of every body (GEN.426): the bulk fractions, the equation of state and the reference solver, the draws, density, mass and gravity, moons and small bodies, the inverse and the migration | approved 2026-10-17 (9:40 to 9:55 PM EDT) (GEN.426) |
| [planet-mass-law.md](composition/planet-mass-law.md) | How heavy a planet is once no class supplies the mass: isolation mass at the formation distance, giant-impact assembly, and the core-accretion runaway to a gas giant (GEN.650) | built as a pure law, not wired in (class lane, 2026-10-24); stand-ins to tune with GEN.509 and GEN.510 |
| [interior-heat.md](composition/interior-heat.md) | The heat of a rocky body's interior (GEN.699): radiogenic plus primordial heat from one thermal history of the mantle, the Urey ratio as an output, the two rheology draws, the stand-ins | built (lane H, 2026-10-24) |
| [magma-ocean-redox-model.md](composition/magma-ocean-redox-model.md) | The surface redox of a magma ocean from its base pressure through the ferric share (GEN.1005 part 1): the pure module, the stand-in rows MR1 to MR7, and what part 2 flips in the atmosphere table | module and registry rows built, nothing wired (Aux 3, 2026-10-30) |
| [dry-carbon-model.md](composition/dry-carbon-model.md) | The carbon a dry rocky body keeps in its mantle from the carbon of its building blocks, the share already differentiated and the core's share (GEN.1005 part 1 of b1): the pure module, the stand-in rows DC1 to DC3, and what the flip changes in b1. |
| [comet-species-model.md](composition/comet-species-model.md) | The N2, Ar, OCS, SO2, O2, Kr and Xe of a comet's ice (row k7) from the host's own ratios to carbon times a trapping share, with 67P as the calibration (GEN.1005): the pure module, the stand-in rows CK1 to CK6, and what remains a stand-in (GEN.1021). |
| [composition-spine-contracts.md](composition/composition-spine-contracts.md) | The interfaces of the composition spine: the bulk composition and mixture accessors, the structure kernels, the bulk draw, the boundaries with the other lanes | contracts of lane I, kept in step with the code (GEN.428, GEN.429, GEN.411) |
| [giant-interiors.md](composition/giant-interiors.md) | The interior of a giant planet, an ice giant and a brown dwarf (GEN.641, GEN.642): the layers, metallic hydrogen and helium rain, dilute and compact cores, with the state of every source |
| [giant-heavy-element-law.md](composition/giant-heavy-element-law.md) | The heavy-element mass of a giant: [R0710] read from the paper, against [R4352] (GEN.460's input) |
| [small-body-size-laws.md](composition/small-body-size-laws.md) | The size-frequency law of small bodies and the binary and satellite fraction by size (GEN.566, GEN.580), and the hosting gate's inputs |
| [small-bodies-composition.md](composition/small-bodies-composition.md) | The composition of belts, fields, comets, icy bodies and moons (GEN.457): which kinds draw composition first, the mixing rule and the compaction law with its hand-over finding, the taxon, comet, icy-body and moon-origin tables, what is stored, the parametrisations, a worked example and the test targets | plan, awaiting Boss's approval |
| [planetary-rings.md](composition/planetary-rings.md) | Planetary rings (GEN.732): kind, origin, edges against the Roche limit, mass, optical depth and ice share, every stand-in labelled in the gap plan registry (BC20 to BC23). |
| [stellar-composition.md](composition/stellar-composition.md) | The composition of every star and remnant (GEN.435): X, Y, Z and the elements from the stored [Fe/H], age and population, and what the star's iron, magnesium and silicon give the planets | checked design, for Boss's approval (GEN.435) |
| [radiation-model.md](composition/radiation-model.md) | The radiation model (GEN.436): the bands, transmission and greenhouse from the infrared opacity, shielding by column, the ground's dose and radon, the surface ultraviolet | checked design, for Boss's approval (GEN.436) |
| [atmosphere-state.md](composition/atmosphere-state.md) | The air's mean molar mass, density, scale height and greenhouse warming worked out from its mixture (GEN.423): the formulas, the greenhouse table of the radiative-convective model (GEN.850, `greenhouse-table.md`) and what changes when the mix changes |
| [greenhouse-table.md](composition/greenhouse-table.md) | The greenhouse the generator reads (GEN.850, GEN.851): a network fitted to the radiative-convective model on the airs the generator makes, the features, the fit and its stated error, the cloud and haze that go into it, the digest that ties it to the downloaded line data, and how to fit it again. |
| [climate-branches.md](composition/climate-branches.md) | The ice-albedo feedback as a loop with a snowball, a warm and an unstable fixed point (GEN.1046, with GEN.972): the pure module `climate_branches` and how the wiring takes a branch. |
| [cloud-albedo.md](composition/cloud-albedo.md) | The cloud term on the albedo of the surface-air fixed point: thinning with warming sized to the AR6 cloud feedback, and the deck of a locked ocean world (GEN.974): design and stand-in rows. |
| [line-opacity.md](composition/line-opacity.md) | The infrared opacity of an air from downloaded line data (GEN.715): ExoMol correlated-k tables of water, carbon dioxide and methane and HITRAN collision-induced pairs reduced at setup, the random-overlap band flux, the radiative-convective column built on it (GEN.717 carbon dioxide doubling, GEN.718 Titan's haze), the checks against Earth, Mars, Venus and Titan and the labelled stand-ins. |
| [surface-composition.md](composition/surface-composition.md) | The ground of a planet or moon (GEN.441): the surface and crust draw, where ice is stable, the ground colour by band, and the potassium, thorium and uranium reservoirs of ultramafic, basaltic and felsic crust. |
| [sky-colour-analysis.md](composition/sky-colour-analysis.md) | The read-only study behind the atmosphere colour: what exists, the method, the effort and how well each source was checked | study |
| [atmosphere-sources-check.md](composition/atmosphere-sources-check.md) | The check of the real-world sources behind the atmosphere colour chain: what was read, what is second-hand, what Boss must research | study |
| [bulk-composition-sources-check.md](composition/bulk-composition-sources-check.md) | The check of the mass-radius and equation-of-state sources behind the bulk composition chain: constants, method order, the core fraction from the star and what Boss must research | study |
| [tunable-gaps.md](composition/tunable-gaps.md) | Whether the unsourced composition weights can be seeded random draws tuned to published populations: the class of each of the 48 open atmosphere rows, the shared selection models, the 17 builder groups, the bracket corrections and what Boss must supply | study |
| [weight-research-dispositions.md](composition/weight-research-dispositions.md) | The disposition of every source of the weight research (GEN.480 to GEN.494): what each paper or data set gives, its label and what it can tune | research |
| [kepler-recovery.md](composition/kepler-recovery.md) | The Kepler transit-detection selection model (GEN.504): the chance a transiting planet is recovered, fitted to the DR25 injections | built |
| [mass-radius-sample.md](composition/mass-radius-sample.md) | The mass-radius sample selection model (GEN.505): the chance a transiting planet passes the mass and radius precision cut | built |
| [jwst-selection.md](composition/jwst-selection.md) | The JWST target-selection model (GEN.506): the chance a transiting planet has a published atmosphere spectrum | built |
| [comet-admission.md](composition/comet-admission.md) | The comet-admission selection model (GEN.507): the chance a known comet has a published coma composition, by group and perihelion | built |
| [wd-pollution-detection.md](composition/wd-pollution-detection.md) | The white-dwarf pollution-detection selection model (GEN.508): how an oxygen fugacity is reported and what the spread of the 22 values means | built |
| [retention-tuning.md](composition/retention-tuning.md) | The check of the retention and lava-world rows d1, i6 and j1 (GEN.511) | built |
| [wd-redox-tuning.md](composition/wd-redox-tuning.md) | The retune of the white-dwarf redox rows a3 and a4 with the detection model (GEN.513) | built |
| [radioisotope-tuning.md](composition/radioisotope-tuning.md) | The retune of the radioisotope rows o2, o4 and o6 (GEN.514) | built |
| [radius-valley-tuning.md](composition/radius-valley-tuning.md) | The tuning of the envelope rows e1 and e2 to the Kepler radius valley (GEN.509) | built |
| [planet-architecture.md](composition/planet-architecture.md) | The planet count, spacing and giant architecture model (GEN.728 to GEN.730, GEN.828, GEN.829, GEN.920, TEST.395): sourced occurrence bands, the Kepler multiplicity harness, the census and what wiring into the system will change | composition |
| [envelope-escape-wiring.md](composition/envelope-escape-wiring.md) | The hydro-based envelope loss wired to rows e1 and e2 and the radius valley, with the 1 Gyr onset as a tuned stand-in (GEN.822) | built |
| [water-delivery-tuning.md](composition/water-delivery-tuning.md) | The delivered-water rows b1 to b4 from the Bern population synthesis: the mass law of GEN.510 and the joint formation-distance law that replaced it (GEN.802) | built |
| [volatile-delivery-formation-distance.md](composition/volatile-delivery-formation-distance.md) | Water-rich share and delivery classes against formation distance from the Bern population synthesis, and how it relates to the GEN.510 mass law (GEN.481); wired into the rows by GEN.802 | research note |
| [magma-ocean-redox-against-mass.md](composition/magma-ocean-redox-against-mass.md) | Magma-ocean surface redox against mass, the later-mantle law and small-body building-block redox checked against first-hand sources; Na and K condensation temperatures (GEN.482) | research note |
| [impact-erosion-and-late-delivery.md](composition/impact-erosion-and-late-delivery.md) | Late accretion, the volatiles the impactors bring and the air they blow away, checked against the Earth, Venus and Mars (GEN.739) | research note |
| [volatile-inventory-against-planet-totals.md](composition/volatile-inventory-against-planet-totals.md) | The carbon and nitrogen of the rocky outgassing inventory against the totals of the Earth, Venus, Mars and Titan, and the two brackets that changed (GEN.813) | research note |
| [atmospheric-escape-model.md](composition/atmospheric-escape-model.md) | The thermosphere model, Jeans and diffusion-limited escape, the XUV history and the hydrodynamic mass loss of a hydrogen envelope, checked against the Earth, Venus, Mars, Titan and nine sub-Neptunes (GEN.738) | research note |
| [clathrate-brine-wiring.md](composition/clathrate-brine-wiring.md) | Hydrate cages in the ice skin and brine lids read by the generator (GEN.762) | built |
| [iron-branch-tuning.md](composition/iron-branch-tuning.md) | The iron branch of a rocky planet: star-matched or iron-enriched, tuned to the homogeneous samples (GEN.512) | built |
| [giant-planet-tuning.md](composition/giant-planet-tuning.md) | The retune of the giant enrichment and C/O rows f5 and f7 (GEN.516) | built |
| [brown-dwarf-types.md](composition/brown-dwarf-types.md) | The brown-dwarf type map from the [R2374] temperature polynomials, rows g1 to g5 (GEN.517) | built |
| [aerosol-tuning.md](composition/aerosol-tuning.md) | The titanium cold-trap step, row h4, and what waits for the chemistry module (GEN.519) | built |
| [design-flags.md](composition/design-flags.md) | The design-flag registry for the life-linked trace gases m3, m4 and m11 (GEN.525) | built |
| [comet-composition-tuning.md](composition/comet-composition-tuning.md) | The retune of the comet rows k1 to k5 and k7: coma class counts, CO2/H2O and CO/H2O, four single-comet species (GEN.520) | built |
| [host-activity-tuning.md](composition/host-activity-tuning.md) | The host-star rotation track, row p6: a uniform percentile through the three published tracks, and what waits for the Mors tracks (GEN.523) | built |
| [equilibrium-chemistry.md](composition/equilibrium-chemistry.md) | The gas-phase equilibrium module: Gibbs minimisation with the NASA polynomials, its check against an independent computation, and the pressure-dependent edges of the hot-Jupiter regimes (GEN.518) | built |
| [atmosphere-composition.md](composition/atmosphere-composition.md) | The atmosphere composition research (GEN.431): the findings that change the design and the code, the one gas and compound registry and its consumers, the sources reached and refused, and what Boss must research | study |
| [atmosphere-composition-parts/](composition/atmosphere-composition-parts/01-gas-table-thermodynamics.md) | The fifteen parts behind it: [gas table and thermodynamics](composition/atmosphere-composition-parts/01-gas-table-thermodynamics.md), [optics](composition/atmosphere-composition-parts/02-gas-optics-uv-visible.md), [infrared and greenhouse](composition/atmosphere-composition-parts/03-infrared-greenhouse-thermal.md), [photochemistry, lifetimes and escape](composition/atmosphere-composition-parts/04-photochemistry-lifetimes-escape.md), [mixing ratios by class](composition/atmosphere-composition-parts/05-mixing-ratios-by-class.md), [condensates, clouds, hazes and dust](composition/atmosphere-composition-parts/06-condensates-clouds-hazes-dust-minerals.md), [radiation shielding and ground sources](composition/atmosphere-composition-parts/07-radiation-shielding-and-ground-sources.md), [verification of Boss's research note, A](composition/atmosphere-composition-parts/08-verification-of-boss-added-research-A.md) and [B](composition/atmosphere-composition-parts/08-verification-of-boss-added-research-B.md), [weighted composition tables](composition/atmosphere-composition-parts/09-weighted-composition-tables.md) and the check of Boss's second note ([A](composition/atmosphere-composition-parts/10-verification-second-note-A.md), [B](composition/atmosphere-composition-parts/10-verification-second-note-B.md), [C](composition/atmosphere-composition-parts/10-verification-second-note-C.md), [D](composition/atmosphere-composition-parts/10-verification-second-note-D.md) and [E, the NH3 lifetime](composition/atmosphere-composition-parts/10-verification-second-note-E-nh3-lifetime.md)) | study |
| [mixture-draw.md](composition/mixture-draw.md) | How a body's atmosphere mixture is drawn from the regimes the possibility tables draw: the equilibrium solver over the NASA Glenn thermodynamics, the families (rocky, giant, evolved, ice, coma, lava), the stand-ins and the limits (GEN.479) | design |
| [surface-pressure-supply-law.md](composition/surface-pressure-supply-law.md) | How much gas a rocky body has delivered to its surface: the mantle inventory, the share degassed by the body's age from its heat release, the ocean and carbonate-silicate sinks, and the anchors (Earth, Venus) the tests keep (GEN.726) | design |
| [body-composition-research.md](composition/body-composition-research.md) | The body composition research (GEN.432, GEN.433): what stars, rocky and icy planets, giants, small bodies, moons and the ground are made of, the chain from the star to the derived values, the findings against the code, the five decisions, the work clusters and the invariants a test can assert | study |
| [body-composition-parts/](composition/body-composition-parts/A-stars.md) | The sixteen parts behind it: [A stars and remnants](composition/body-composition-parts/A-stars.md), [B rocky planets](composition/body-composition-parts/B-rocky.md), [C volatile-rich worlds and giants](composition/body-composition-parts/C-volatile-and-giants.md), [D small bodies and moons](composition/body-composition-parts/D-small-bodies-moons.md), [E ground, radiogenic heat and shielding](composition/body-composition-parts/E-ground-radiogenic-shielding.md), [F check of Boss's research note](composition/body-composition-parts/F-boss-note-check.md), [G1 comets and photolysis lifetimes](composition/body-composition-parts/G1-comets.md), [G2 radiogenic ratios and Thorngren](composition/body-composition-parts/G2-radiogenic-thorngren.md), [G3 magma ocean and lava worlds](composition/body-composition-parts/G3-magma-ocean-lava.md), [G4 ices, brines and licences](composition/body-composition-parts/G4-ices-brines-licences.md), [H the GEN.465 second extraction pass and the water-share result](composition/body-composition-parts/H-gen465-second-pass-tables.md) | study |

## Habitability, species and population

`habitability/`

| Note | What it covers | Status |
|---|---|---|
| [habitability-index.md](habitability/habitability-index.md) | The planetary habitability index (PHI-4) and the protection gear | partly built (scores built, GEN.84) |
| [survival-equipment-fields.md](habitability/survival-equipment-fields.md) | The survival equipment line: what a human needs on a world, as required and recommended items | built (GEN.655, DB.93; chips, regional temperatures, toxic rows and dust to come) |
| [day-and-night.md](habitability/day-and-night.md) | The share of the time a star is up on every planet and moon, the solar day and the hours of daylight, and the extra night a moon has in its planet's shadow | built (GEN.767) |
| [condition-based-planet-classes.md](habitability/condition-based-planet-classes.md) | A planet or moon's class worked out from its finished conditions (air, temperature, water, size, gravity) instead of drawn first: the rule file, the core classes, the shadow sweep and the cut-over (GEN.645 to GEN.648) | partly built (class lane, 2026-10-23; the classifier, sweep and class pages are built, the cut-over waits on the composition lanes) |
| [life-classes.md](habitability/life-classes.md) | The life class system: nine classes of carbon life (what each drinks, eats, breathes and gives off), the stage ladder, the gas ledger and the Life and air chart, scoring, the dominant class, storage and display (GEN.862, GEN.886 to GEN.894) | design, order final (Boss, 9:20 PM EDT 2026-10-28); nothing built yet |
| [life-timeline-and-oxygen.md](habitability/life-timeline-and-oxygen.md) | The life timeline (when each life stage is reached from conditions) and the free oxygen that photosynthesis puts into the air | built (GEN.750, GEN.855) |
| [life-class-rules.md](habitability/life-class-rules.md) | The rules the life class classifier reads (gates, windows, energy, supply per class, each number with its tier and tweak), generated from the rule file (GEN.886) | generated (reference) |
| [adaptation-traits.md](habitability/adaptation-traits.md) | Adaptation traits: how a body's life gets tougher against its own stress (radiation, UV, sulfide, acid, brine, heat or cold, CO2, low O2, toxins), the checked table, grades and animal caps, the draw, storage and tests (GEN.875, GEN.876) | design (lane Aux 2, 2026-10-28); nothing built yet |
| [atmospheres-retention-and-classes.md](habitability/atmospheres-retention-and-classes.md) | Atmosphere retention and class timing | research |
| [activity-magnetism-radiation-hydrosphere.md](habitability/activity-magnetism-radiation-hydrosphere.md) | Stellar activity, magnetism, radiation dose and hydrosphere | partly built |
| [population-and-politics.md](habitability/population-and-politics.md) | Species, civilizations, territories, tech levels and facilities | built |

## Orbits, navigation and courses

`orbits/`

| Note | What it covers | Status |
|---|---|---|
| [orbital-updates.md](orbits/orbital-updates.md) | Moving things: positions, thresholds, gravity, collisions and spin | partly built |
| [tidal-chain.md](orbits/tidal-chain.md) | The tide chain (GEN.842, GEN.954 to GEN.960): rheology from temperature and melt, spin torques of moons and air, moon orbit evolution and ice, in lock-free and wiring steps | partly built (stage 1a, GEN.954) |
| [ocean-tides.md](orbits/ocean-tides.md) | Ocean tide height and range at a world's surface from each moon and the star (GEN.1034): the equilibrium formula, the loading factor and its refusal for heavy liquids, the stand-ins | module built (physics/ocean_tides.py, lock-free); columns on train 7 |
| [nightside-collapse-model.md](composition/nightside-collapse-model.md) | Whether the air of a tidally locked rocky planet freezes onto its night side and how much stays (GEN.1044, oddity W19): the CO2 critical pressure from [R4767]'s fit, the thin-air night-side temperature, other gases by the same balance; pure module, wiring in GEN.1051. |
| [orbital-solvers-and-integrators.md](orbits/orbital-solvers-and-integrators.md) | Kepler solvers, integrators, thresholds and the influence search | research |
| [orbital-update-check-2026-10-18.md](orbits/orbital-update-check-2026-10-18.md) | Boss's four orbital motion update steps checked against the design, with the F1 to F5 answers | research |
| [collisions-and-mergers.md](orbits/collisions-and-mergers.md) | How often bodies collide and how to detect it | research |
| [navigation-frames.md](orbits/navigation-frames.md) | Nested navigation frames, bearing, mark and warp speed | built |
| [course-routing.md](orbits/course-routing.md) | Course routing and nearby search | partly built |
| [course-avoidance.md](orbits/course-avoidance.md) | Keep-out radii and path bending | research |
| [sky-view.md](orbits/sky-view.md) | The view from a planet | research |
| [gravity-map.md](orbits/gravity-map.md) | The field of a sector as a heat map | plan (the 9.1 headline) |

## Maps and the web interface

`interface/`

| Note | What it covers | Status |
|---|---|---|
| [admin-control-suite.md](interface/admin-control-suite.md) | The admin control suite: research on packing many options into a quiet interface, the Settings and Generate redesign, the shell, job strip and palette, Boss's 13 answers and the items ADM.101 to ADM.123 | designed (in 8.1 phase 1, goal #2350) |
| [galaxy-drilldown-navigation.md](interface/galaxy-drilldown-navigation.md) | The Galaxy Map's drill-down | built (partly replaced by the fly-through) |
| [drilldown-region-sizes.md](interface/drilldown-region-sizes.md) | Step sizes and region shapes for the drill-down | built (MAP.151) |
| [galaxy-map-performance.md](interface/galaxy-map-performance.md) | The Galaxy Map node path held to a budget: the MAP.249 benchmark over a synthetic 46,000-star galaxy, its numbers and where the time goes | built (MAP.249) |
| [text-style.md](interface/text-style.md) | The plain-words style guide for every sentence about a planet, moon, belt, star or system: short sentences, plain words, one name for a class, sentence-case labels, and the test that holds new text to it | built (DOC.83) |
| [references.md](interface/references.md) | Help > References, the one out-of-universe section: the address scheme and link helper, the word rule and page crawl, the shell, the generator and the page tree | R1 and R2 built (DOC.87, DOC.88), R3 to R8 planned |
| [body-rendering.md](interface/body-rendering.md) | Rendering bodies from composition: the look record, the shared painter, solid worlds, giants, stars and phenomena, items by version | approved (R1, R3, R4a built in 8.1) |
| [edit-everything.md](interface/edit-everything.md) | The 8.6 headline: every admin change is a plan the admin approves; aspects, locks, regions, the rearranger | plan |
| [roles-and-requests.md](interface/roles-and-requests.md) | The 10.0 headline: user, moderator and admin roles; edit requests are change plans that the role above approves | plan |
| [fly-through-view-distance.md](interface/fly-through-view-distance.md) | The fly-through map: view distance, near field, free camera | partly built |
| [object-pages.md](interface/object-pages.md) | Object pages (UX.168 to UX.178): the page of every star, planet, moon, belt and comet, the property sheet, the children API and the summary trees on the system and sector pages | built |
| [zoom-star-visibility.md](interface/zoom-star-visibility.md) | Smoother star visibility while zooming | research (superseded in part by MAP.152) |
| [galaxy-map-wire-format.md](interface/galaxy-map-wire-format.md) | What the browser downloads for the map | built |
| [nebula-map-visibility.md](interface/nebula-map-visibility.md) | Making nebulae visible on the Galaxy Map | partly built |
| [map-ui-and-frontend-libraries.md](interface/map-ui-and-frontend-libraries.md) | The map interface and front-end libraries | research |
| [3d-rendering-camera-research.md](interface/3d-rendering-camera-research.md) | Camera proximity fading in dense 3D scenes | research (background) |
| [settings-seo-and-accounts.md](interface/settings-seo-and-accounts.md) | Settings, the configuration page, SEO and accounts | research |

## Server, queue and operations

`operations/`

| Note | What it covers | Status |
|---|---|---|
| [cache-warm-up.md](operations/cache-warm-up.md) | The page cache shared through Redis, the data tags each answer reads and each change bumps, and the warm-up that fills it after an update (PERF.104) | built (Core 2, 2026-10-16) |
| [one-run-per-job.md](operations/one-run-per-job.md) | A Generate page or API job is one `run-job` process and one wrapper run over all the selected commands: the job spec, a follower list per verb, start checks in `prepare` hooks, job_runner following one run (PERF.203) | design note (Aux 2, 2026-10-30), not built |
| [logging-standard.md](operations/logging-standard.md) | One log record (time, RFC 5424 severity, component, message, fields), its file, JSON and display writings, the one wording of a stage's start, finish and failure, and every place the program logs today | partly built (the layer is built, UX.121; the places move in UX.122, UX.123, TEST.183) |
| [slow-page-answers.md](operations/slow-page-answers.md) | A page request that is slow to build queues one job that stores the answer in the tagged page cache, with a cost rule, a status route, a node on the queue page, a cancel and a fallback (PERF.105) | built (Core 2, 2026-10-17; the search and phenomena pages use it, UX.118, #2089) |
| [backups-and-repair.md](operations/backups-and-repair.md) | One bz2 archive of the settings, the control database and the parity data, and restore and repair from it (settings 8.1, chain 8.9) | plan |
| [redis-and-the-work-queue.md](operations/redis-and-the-work-queue.md) | How Redis is used and what to change (PERF.86) | research (decided; being built) |
| [work-queue-audit.md](operations/work-queue-audit.md) | What starts work and what moves to the queue | built (replaced in part by PERF.86) |
| [performance-eta-queue-and-caching.md](operations/performance-eta-queue-and-caching.md) | Progress estimates, the job queue and caching | research (queue under redo) |
| [slow-reads-and-timeouts.md](operations/slow-reads-and-timeouts.md) | Keeping slow database reads inside the time limit | mostly built |
| [db-check-and-parity-repair.md](operations/db-check-and-parity-repair.md) | Database check, parity repair and long migrations | partly built |
| [ops-scheduling-and-rotation.md](operations/ops-scheduling-and-rotation.md) | Apache reload, daily maintenance and locks | partly built |
| [api-design-standards.md](operations/api-design-standards.md) | API versions, key scopes, uploads and logging | research |
| [login-brute-force-protection.md](operations/login-brute-force-protection.md) | Login lockouts and two-factor sign-in | built |
| [test-suite-rebuild.md](operations/test-suite-rebuild.md) | The faster, modular test suite | partly built (stages 1 to 4) |

## Research, sources and gap filling

`research/`: source checks, the gap-filling plan and the papers to fetch.

| Note | What it covers | Status |
|---|---|---|
| [gap-filling-plan.md](research/gap-filling-plan.md) | How the open gaps are filled where the public record is silent: the 57 gaps with a derived, bracketed, theory or placeholder fill, the stand-in registry convention, the first Gemini batch (unverified) and the work order | research |
| [captured-moon-sources.md](research/captured-moon-sources.md) | The sourced head counts, size laws, capture probabilities and Kozai gap of the irregular moons of the four giants (GEN.1075), with the corrections to the captured asteroid moons design | built (report; registry rows CM1 to CM16) |
| [oddity-rates.md](research/oddity-rates.md) | The sourced rates and thresholds of the weird worlds oddities W3, W4, W5, W10, W14, W15 and W20 (GEN.1037), with the three that changed and the gates | built (table; rates read from papers) |
| [interstellar-field-sources.md](research/interstellar-field-sources.md) | What the sources say for the interstellar asteroid field modules (GEN.1098): collision speed window, unbound mass share, supernova belt share (registry rows IF1 to IF12) | built (read report) |
| [rogue-ejection-kick.md](research/rogue-ejection-kick.md) | Ejection speed of a rogue against scatterer mass and distance (about 3 km/s x sqrt(Jupiter masses / AU)), moon survival, direction statistics (GEN.1113) | built (read report; rows RK1, RK3 to RK8) |
| [rogue-planet-review.md](research/rogue-planet-review.md) | GEN.1109: what each rogue-planet quantity is (physics, sourced stand-in or flat rate), what the architecture ejects per star against microlensing, and the replacement items | research |
| [data-catalog.md](research/data-catalog.md) | Every dataset and table the project already holds, ranked by usefulness: the six setup downloads, the data held only in the project files, the repo tables and the stand-in registry, with the TODO items each could serve; required reading before hand-tuning a number | research |
| [wiring-register.md](research/wiring-register.md) | The wiring register of the guard GEN.754: every lock-free module, dataset reader, data file and registry row nothing reads yet, with the open TODO item that will wire it; `python scripts/wiring_guard.py --check` fails when one is missing | research |
| [burn-down-register.md](research/burn-down-register.md) | GEN.946 step 1: a disposition (replace with a derivation, source, keep as registered stand-in) for each L3 and L4 row of the stand-in registry, with the item that does it | research |
| [stand-in-values.md](research/stand-in-values.md) | Every stand-in value with its source tier (T1 to T9) and tweak level (L0 to L4), generated from the registry of the gap-filling plan (GEN.725) | generated (reference) |
| [papers-to-fetch.md](research/papers-to-fetch.md) | The papers behind refused publishers that Boss downloads himself, and the arXiv papers a research thread reads instead, with Crossref citations (gap-filling plan section 8, GEN.531) | reference |
| [vizier-courteau-comeron.md](research/vizier-courteau-comeron.md) | The VizieR search for the [R0908] scatter (GR5) and the Comeron thick-disc tables (GR3): the [R0865] fit, the raw Vc/sigma0 scatter and the corrected references | research |
| [gemini-batch1-check.md](research/gemini-batch1-check.md) | The claim-by-claim check of Gemini batch 1 (51 claims) against read sources, the figures that can be used and the requests it answers | research |
| [boss-papers-2026-10-17.md](research/boss-papers-2026-10-17.md) | What Boss's 28 uploaded papers fill in the gap plan, what stays open, where they contradict our files and the printing faults found | research |
| [open-gaps-closest-laws-2026-10-17.md](research/open-gaps-closest-laws-2026-10-17.md) | Every open gap (64 lines) with the nearest law or model, where to look, the free tool, the stand-in quality and the deciding rule; the citation errors found | research |
| [ad1-methane-sources-check-2026-10-17.md](research/ad1-methane-sources-check-2026-10-17.md) | AD1: where the free methane absorption data is (Oxford k-tables from [R2265], PDS volume gbat_0001), what is still open (ammonia, water) and the corrections to the pasted source note | research |
| [Atmospheric Toxicity](research/imported/atmospheric-toxicity.md), [Chemical Habitability](research/imported/chemical-habitability.md), [Exoplanetary Habitability](research/imported/exoplanetary-habitability.md), [Naturally Occuring Ionizing Radiation](research/imported/naturally-occuring-ionizing-radiation.md) | habitability and atmosphere notes |
| [Planetary Habitability Index](research/imported/planetary-habitability-index.md), [Mathematical and Algorithmic Implementation of the Planetary Habitability Index](research/imported/mathematical-and-algorithmic-implementation-of-the-planetary-habitability-index.md), [Planetary Habitability and Speculative Xenobiology](research/imported/planetary-habitability-and-speculative-xenobiology.md) | [habitability-index.md](habitability/habitability-index.md) |
| [Speculative Xenobiology Examples](research/imported/speculative-xenobiology-examples.md), [Speculative Xenobiology Extremes](research/imported/speculative-xenobiology-extremes.md) | [habitability-index.md](habitability/habitability-index.md) |
| [Environmental.Extremes.and.Protection.Research](research/imported/environmental-extremes-and-protection-research.md), [Thermal.Respiratory.Protection.Research](research/imported/thermal-respiratory-protection-research.md) | [habitability-index.md](habitability/habitability-index.md) |
| [Computational Astrodynamics](research/imported/computational-astrodynamics.md), [Observational Kinetics for Rotational Vectors](research/imported/observational-kinetics-for-rotational-vectors.md) | [orbital-updates.md](orbits/orbital-updates.md) and related |
| [Orbital Position and Vector Mathematical Foundations](research/imported/orbital-position-and-vector-mathematical-foundations.md), [Orbital Position and Vector Update Algorithms](research/imported/orbital-position-and-vector-update-algorithms.md), [Orbital Update Collision Detection](research/imported/orbital-update-collision-detection.md), [Orbital Update Full Algorithm](research/imported/orbital-update-full-algorithm.md) | [orbital-updates.md](orbits/orbital-updates.md) |
| [Star Trek Anomalies and Science](research/imported/star-trek-anomalies-and-science.md), [Star Trek Anomaly Probabilities](research/imported/star-trek-anomaly-probabilities.md) | [anomalies.md](systems/anomalies.md) |
| [Technological Assessment](research/imported/technological-assessment.md) | [population-and-politics.md](habitability/population-and-politics.md) |
| [Library Migration Workflow](research/imported/library-migration-workflow.md), [Web UX Development Notes](research/imported/web-ux-development-notes.md), [Web UX and Job Management Guide](research/imported/web-ux-and-job-management-guide.md) | [library-migration.md](architecture/library-migration.md) and related |

### Imported research

`research/imported/` holds Boss's own research documents, kept as received (content untouched; the file names are lower-case with hyphens). Each starts with a banner linking the notes above that use it. They say what Boss researched, not what planetGen does.

| File | Used by |
|---|---|
| [Atmospheric Toxicity](research/imported/atmospheric-toxicity.md), [Chemical Habitability](research/imported/chemical-habitability.md), [Exoplanetary Habitability](research/imported/exoplanetary-habitability.md), [Naturally Occuring Ionizing Radiation](research/imported/naturally-occuring-ionizing-radiation.md) | habitability and atmosphere notes |
| [Planetary Habitability Index](research/imported/planetary-habitability-index.md), [Mathematical and Algorithmic Implementation of the Planetary Habitability Index](research/imported/mathematical-and-algorithmic-implementation-of-the-planetary-habitability-index.md), [Planetary Habitability and Speculative Xenobiology](research/imported/planetary-habitability-and-speculative-xenobiology.md) | [habitability-index.md](habitability/habitability-index.md) |
| [Speculative Xenobiology Examples](research/imported/speculative-xenobiology-examples.md), [Speculative Xenobiology Extremes](research/imported/speculative-xenobiology-extremes.md) | [habitability-index.md](habitability/habitability-index.md) |
| [Environmental.Extremes.and.Protection.Research](research/imported/environmental-extremes-and-protection-research.md), [Thermal.Respiratory.Protection.Research](research/imported/thermal-respiratory-protection-research.md) | [habitability-index.md](habitability/habitability-index.md) |
| [Computational Astrodynamics](research/imported/computational-astrodynamics.md), [Observational Kinetics for Rotational Vectors](research/imported/observational-kinetics-for-rotational-vectors.md) | [orbital-updates.md](orbits/orbital-updates.md) and related |
| [Orbital Position and Vector Mathematical Foundations](research/imported/orbital-position-and-vector-mathematical-foundations.md), [Orbital Position and Vector Update Algorithms](research/imported/orbital-position-and-vector-update-algorithms.md), [Orbital Update Collision Detection](research/imported/orbital-update-collision-detection.md), [Orbital Update Full Algorithm](research/imported/orbital-update-full-algorithm.md) | [orbital-updates.md](orbits/orbital-updates.md) |
| [Star Trek Anomalies and Science](research/imported/star-trek-anomalies-and-science.md), [Star Trek Anomaly Probabilities](research/imported/star-trek-anomaly-probabilities.md) | [anomalies.md](systems/anomalies.md) |
| [Technological Assessment](research/imported/technological-assessment.md) | [population-and-politics.md](habitability/population-and-politics.md) |
| [Library Migration Workflow](research/imported/library-migration-workflow.md), [Web UX Development Notes](research/imported/web-ux-development-notes.md), [Web UX and Job Management Guide](research/imported/web-ux-and-job-management-guide.md) | [library-migration.md](architecture/library-migration.md) and related |

## Archive

Replaced designs and finished reviews are in [../archive/](../archive/): [galaxy-coordinate-system-shells.md](../archive/galaxy-coordinate-system-shells.md), [galaxy-disk-density-rev2.md](../archive/galaxy-disk-density-rev2.md), [galaxy-coordinate-system-cylinder.md](../archive/galaxy-coordinate-system-cylinder.md) and [habitability-atmosphere-sanity-review.md](../archive/habitability-atmosphere-sanity-review.md). They are history.
