# The composition model: what every body, star and atmosphere is made of

**Status (2026-10-17):** plan; nothing is built. This is the master note of Boss's directive (2026-10-16, 11:58 PM EDT). It joins [atmosphere-colour.md](atmosphere-colour.md) (the atmosphere and the colour), the bulk composition of section 5 of that note, and the new stellar and radiation parts. GEN.437 turns it into the checked design once the research lanes and the four designs below are done. Every section the research lanes write goes into the design documents named here.

## 1. What Boss decided

- Atmosphere composition and the composition of planets, moons, stars and everything is the joint headline of 9.0, together with the atmosphere colour.
- Two research lanes first, one for atmospheres and one for bodies (planets, moons, stars, remnants and the rest): the compounds, minerals and elements that make them up, what is there and all the statistics needed to use it in calculations. All notes go into the design documents.
- Every body gets figures on its composition from actual science, and not made up where that can be avoided. Drawing the density and the other figures that are drawn today and then filling in a composition that fits is acceptable, and so is building from the most basic composition up; either works, and the design picks the simplest consistent one per kind of body.
- 9.0 completes it, together with the atmosphere colour, and adds what becomes derivable: the atmospheric radiation (the atmosphere's transmission and greenhouse effect, since we know what is in the air) and the surface radiation (from the atmosphere and the ground composition).
- The building blocks go back from 9.0 to 8.1 (as early as their blockers allow), then a collision and consolidation sweep merges items that do almost the same thing and suit one worker, then the issue clusters and all documents are updated.
- Added 2026-10-17 at 1:40 AM EDT: every composition category is a gated table of possibilities with weighted probabilities that add up to 100 percent (section 4b).
- The standing rule stays: nothing fixed per class and nothing copied from the Solar System or the Milky Way; every input is a stored physical quantity, a derived value or a seeded draw, and the Sun, Earth and the rest are test targets.

## 2. The chain

One chain from the element to the body to the derived value:

1. The **element table** (GEN.438): atomic data and the radioactive isotopes. No abundances.
2. The **compound and material tables**: the gases and compounds of atmospheres (GEN.410) and the minerals, ices, metals and envelopes of bodies (GEN.427), each with its stats and sources; molar masses come from the elements.
3. The **stellar composition relations** (GEN.439) and their derived function (GEN.440): a star's composition from the stored [Fe/H], age and population, with no new column.
4. The **bodies**: bulk fractions that drive density, mass and gravity (GEN.428, GEN.429), the surface composition (GEN.441), the small bodies (GEN.442), and the atmosphere's mixture (GEN.411) with ozone and aerosols.
5. The **derived values**: the atmosphere colour (GEN.416 to GEN.419), the atmospheric radiation (GEN.443), the surface radiation (GEN.444), and what reads them (habitability, GEN.445; the text, GEN.446).

## 3. Where each item sits

Research, designs and tables are building blocks with no behaviour change and sit in 8.1. Generation and storage change stored values (generator epoch bumps and migrations) and sit in 8.4, where the planet work already is. The radiation, the colour and what shows them are the 9.0 features. Decided (Boss, 2026-10-17, 12:13 AM EDT, a card): the value-changing generation chain (composition fractions, the mixture, the albedo records, ozone and aerosols, the migrations and epoch bumps) stays in 8.4 and does not join the 8.1 stamp. Boss also tapped Core 4 as the owner of the star composition items (GEN.435, GEN.439, GEN.440, TEST.211) together with the other composition design and table items.

**Research (read-only lanes, drafts to /mnt/project-files/design-drafts/)**

| Where | Item | What |
|---|---|---|
| done | GEN.432 | Research: what rocky, icy and giant bodies and their small relatives are made of (merged, see section 4) |
| done | GEN.433 | Research: what stars and stellar remnants are made of (merged, see section 4) |
| done | GEN.431 | Research: what is in the atmospheres of every kind of world, and the numbers each compound needs |
| done | GEN.434 | Research: radiation through atmospheres and at the ground, and the numbers it needs (covered by the two research lanes and merged) |
| done | GEN.425 | Research: verify the mass-radius and equation-of-state sources and say what is missing ([bulk-composition-sources-check.md](bulk-composition-sources-check.md)) |

**Designs (docs/design)**

| Where | Item | What |
|---|---|---|
| 8.1 phase 1 | GEN.409 | The atmosphere colour design: the checked plan for the gas table, the mixture, the stored values, ozone, aerosols and the colour model |
| 8.1 phase 1 | GEN.426 | The bulk composition design: the fractions, the equation of state, the draws and the derived density, mass and gravity |
| 8.1 phase 1 | GEN.435 | The stellar composition design: the abundance relations, what is stored and what is derived, and the remnants |
| 8.1 phase 1 | GEN.436 | The radiation design: atmospheric transmission and greenhouse by band, shielding by column, and the ground's own dose |
| 8.1 phase 1 | GEN.457 | The small bodies and moons composition design: the taxon and comet draws, compaction, the moon origin types and the ice fraction |
| 8.1 phase 1 | GEN.465 | A second extraction pass over the sources the research could not reach and that are reachable now |
| 8.1 phase 1 | GEN.468 | Settle the water share beyond the snow line (the factor of two) |
| 8.1 phase 1 | GEN.437 | The composition model: how the element, compound, material and body tables link, and what is stored against derived |

**Data tables (no behaviour change)**

| Where | Item | What |
|---|---|---|
| 8.1 phase 2 | GEN.472 | One gated possibility-table format and loader: every composition category is a table of possibilities with a class gate and a weight |
| 8.1 phase 2 | GEN.500 | The licence audit: the licence column filled for every data source, and the permission list |
| 8.1 phase 2 | GEN.438 | One table of the elements and their isotopes: atomic data and the radioactive isotopes' half-lives and decay energies |
| 8.1 phase 2 | GEN.410 | One table of the gases and compounds that can be in an atmosphere, with every compound the descriptions mention |
| 8.1 phase 2 | GEN.427 | One table of the bulk materials with their equation-of-state constants and the source of each |
| 8.1 phase 2 | GEN.498 | One oxygen-fugacity scale for the generator (Hirschmann 2021), delta-IW stored on it |
| 8.1 phase 2 | GEN.497 | The cold-end vapour pressures: laboratory-corrected N2, CO, CH4, O2, NH3 and Span and Wagner CO2 |
| 8.1 phase 2 | GEN.496 | The photolysis lifetime table: species lifetimes at 1 AU, active-Sun factors, stellar-spectrum ratios |
| 8.1 phase 2 | GEN.470 | One saturation-pressure function p_sat(species, T) for the gases and ices (Fray and Schmitt polynomials, water by Feistel and Wagner) |
| 8.1 phase 2 | GEN.439 | The stellar composition relations: abundances by population, age and [Fe/H], and the remnants' layers |

**Star composition (derived, no column)**

| Where | Item | What |
|---|---|---|
| 8.1 phase 3 | GEN.440 | Star composition is derived from the stored [Fe/H], age and population, for every star, remnant and brown dwarf |
| 8.1 phase 1 | GEN.463 | The neutron star's radius is derived from sourced mass-radius relations (with the stellar tables) |
| 8.1 phase 1 | GEN.467 | The main-sequence lifetime is derived from sourced stellar relations |
| 8.1 phase 3 | GEN.464 | The disk solids law reads Z over Z0, not 10 to the power of [Fe/H] |
| 8.1 phase 4 | TEST.215 | The Solar System calibration rows live in test files only: a test fails when a generator table reads one |
| 8.1 phase 4 | TEST.214 | Every registered possibility table: the gated weights add up to 100 percent for every class, with a source on every row |
| 8.1 phase 4 | TEST.211 | A star composition sweep: every kind of star, remnant and brown dwarf, fractions, trends and determinism |

**Bodies: generation and storage**

| Where | Item | What |
|---|---|---|
| 8.4 phase 1 | GEN.458 | Each system draws its snow line and stores it once |
| 8.4 phase 1 | GEN.428 | Planets and moons draw their bulk composition fractions instead of carrying a fixed text |
| 8.4 phase 1 | GEN.429 | Bulk density, mass, radius and gravity are derived from the composition fractions with the mass-radius relation |
| 8.4 phase 1 | GEN.459 | One radiogenic module: heat from the isotope-wise decay, drawn from the host star's abundances |
| 8.4 phase 1 | GEN.460 | Giants and brown dwarfs: radius and heavy-element content from the checked laws |
| 8.4 phase 1 | DB.69 | Planets and moons store the composition fractions |
| 8.4 phase 1 | GEN.430 | The planet and moon descriptions state the composition from the stored fractions |
| 8.4 phase 1 | GEN.412 | Generation records the albedo and the greenhouse multiplier it drew, and the temperature is exactly what they give |
| 8.4 phase 2 | GEN.441 | Surface and crust composition: the minerals, ices and ocean on the ground, from the bulk fractions, the temperature and the air |
| 8.4 phase 2 | GEN.461 | The hydrosphere reads the stored water mass fraction (compressed layer depth, base pressure, ice phases) |
| 8.4 phase 2 | GEN.442 | Asteroids, comets, rings and dwarf planets draw their composition from the formation distance and the snow line |
| 8.4 phase 2 | GEN.466 | Moons draw their origin type, ice fraction and composition |
| 8.4 phase 3 | GEN.462 | Ground albedo by band, handed to the atmosphere colour once |
| 8.4 phase 3 | GEN.471 | Clathrate hydrates and eutectic brines for icy crusts and ocean moons |
| 8.4 phase 3 | DB.70 | Planets, moons and small bodies store the surface and small-body composition shares |
| 8.4 phase 2 | GEN.455 | Airless bodies and comets draw an exosphere or coma from their surface composition, temperature and the star's wind |

**Atmosphere: mixture, ozone, aerosols**

| Where | Item | What |
|---|---|---|
| 8.4 phase 2 | GEN.411 | The atmosphere's molar density comes from its per-gas partial pressures, so mixture, pressure and density always agree |
| 8.4 phase 2 | GEN.452 | The condensation table is replaced by the checked saturation curves (the hydrogen sulphide row is a copy of the sulphur dioxide row) |
| 8.4 phase 2 | GEN.454 | Giants, ice giants and brown dwarfs draw their mixture from the equilibrium chemistry sequence |
| 8.4 phase 2 | DB.71 | Planets and moons store the partial pressures of the added gases |
| 8.4 phase 3 | GEN.456 | Lava worlds and hot rocky planets draw a vapour atmosphere over the magma |
| 8.4 phase 2 | GEN.422 | Stop the recalculation redraw and order the steps so pressure, mixture and temperature agree after the cold-side freeze-out |
| 8.4 phase 2 | GEN.423 | Atmosphere density, scale height and the greenhouse multiplier follow from the mixture |
| 8.4 phase 2 | GEN.424 | The balance gas is a named gas with its own molar mass, and habitability reads and writes the mean molar mass of the stored mixture |
| 8.4 phase 2 | DB.67 | Planets and moons store the albedo, the greenhouse multiplier and the molar density they were built with |
| 8.4 phase 3 | GEN.413 | Ozone in the mixture: the column and partial pressure from the oxygen, the star's UV and the Chapman balance |
| 8.4 phase 3 | GEN.414 | Haze, dust and cloud: derived or drawn optical depths for every planet and moon |
| 8.4 phase 4 | DB.68 | Planets and moons store the ozone column and the haze, dust and cloud optical depths |
| 8.4 phase 4 | GEN.415 | Every compound the planet descriptions mention comes from the gas table, and the atmosphere text is written from the stored mixture |

**Proofs and docs of the generation half**

| Where | Item | What |
|---|---|---|
| 8.4 phase 4 | TEST.209 | One consistency sweep over every planet and moon class: composition, density, temperature, ozone, aerosols and text all agree |
| 8.4 phase 4 | TEST.210 | A bulk composition sweep over every class: fractions, density, mass, radius and gravity agree |
| 8.4 phase 4 | DOC.63 | Docs: the atmosphere data model, the gas table and how the values agree |
| 8.4 phase 4 | DOC.65 | Docs: bulk composition, the materials table and how density, mass and gravity follow |

**Colour model**

| Where | Item | What |
|---|---|---|
| 9.0 phase 1 | GEN.416 | Blackbody colour of a star: Planck to the CIE curves to sRGB, as the first part of the photometry library |
| 9.0 phase 1 | GEN.417 | Rayleigh optical depth of the mixture: per-gas cross sections from the table and the air column from pressure and gravity |
| 9.0 phase 1 | GEN.418 | Atmosphere colour of a solid world: the absolute colour of the air with ozone, haze, dust, cloud and the ground's bounce |
| 9.0 phase 1 | GEN.419 | Atmosphere colour of gas and ice giants, read at a reference level in the cloud tops |

**Radiation**

| Where | Item | What |
|---|---|---|
| 9.0 phase 2 | GEN.443 | Atmospheric radiation: transmission by band and the greenhouse effect from the stored mixture, the aerosols and the cloud |
| 9.0 phase 2 | GEN.453 | The radiation constants are corrected: cutoff with moment over radius squared, the dose curve and the field factor from the read tables |
| 9.0 phase 2 | GEN.444 | Surface radiation: the ground-level dose from the air's shielding and the ground's own radiogenic content |
| 9.0 phase 2 | GEN.445 | Habitability, toxicity and the UV index read the derived radiation instead of the proxies |

**In the product**

| Where | Item | What |
|---|---|---|
| 9.0 phase 3 | GEN.420 | The descriptions show the atmosphere colour: a swatch and a word in the planet and moon facts, in every output |
| 9.0 phase 3 | MAP.244 | The maps and pages use the atmosphere colour where they showed the planet colour |
| 9.0 phase 3 | GEN.446 | The descriptions state the composition of stars, surfaces and small bodies, and the radiation at the surface |
| 9.0 phase 3 | DOC.64 | Docs: the atmosphere colour, its method and its sources |

**Proofs, docs and the headline**

| Where | Item | What |
|---|---|---|
| 9.0 phase 4 | TEST.212 | A radiation sweep over every class: transmission, greenhouse, shielding and dose agree with the composition |
| 9.0 phase 3 | GEN.473 | Replace the asteroid taxon shares and the zone kernel by a derived or general distribution (R1) |
| 9.0 phase 3 | GEN.474 | Replace the Earth rock-type library by ground types derived from the bulk composition (R2) |
| 9.0 phase 3 | GEN.475 | Meteorite groups become outputs: regions of composition space reachable from the drawn bulk (R3) |
| 9.0 phase 3 | GEN.476 | Aluminium-26 and iron-60 follow from the star's formation environment (R4) |
| 9.0 phase 4 | DOC.66 | Docs: the composition model of stars, bodies and atmospheres, the radiation it derives and the sources of every number |
| 9.0 phase 4 | GEN.421 | Headline of 9.0: what every body and atmosphere is made of, and the atmosphere colour that follows from it |

## 4. The research lanes

- **Atmosphere research (2026-10-17):** GEN.431, with the re-check of the four sources the first check could not reach. Done and merged: [atmosphere-composition.md](atmosphere-composition.md) (the master: findings, design synthesis, sources, what Boss must research) and its seven parts in [atmosphere-composition-parts/](atmosphere-composition-parts/01-gas-table-thermodynamics.md) (gas table and thermodynamics, optics, infrared and greenhouse, photochemistry and escape, mixing ratios by class, condensates and aerosols, radiation shielding and ground sources). Its findings are GEN.452 to GEN.456, DB.71 and the notes on GEN.409 to GEN.411, GEN.413, GEN.414, GEN.417, GEN.422, GEN.423, GEN.436, GEN.444, GEN.445 and GEN.451.
- **Body composition research (2026-10-17):** GEN.432 (planets, moons, small bodies, giants) and GEN.433 (stars and remnants), done and merged: [body-composition-research.md](body-composition-research.md) (the master: what each body is made of, the chain, the findings against the code, the decisions, the clusters, the invariants) and its five parts in [body-composition-parts/](body-composition-parts/A-stars.md). GEN.425 (mass-radius and equation-of-state sources) is done, see the sources check.
- **Radiation research:** GEN.434 was covered by the two lanes above (the coordinator, 12:17 AM EDT) and is done with their merges.

Each lane writes its draft to /mnt/project-files/design-drafts/ and reports; the TODO thread merges it into the design documents (atmosphere-composition.md and its parts, the body-composition sections, stellar-composition.md, radiation-model.md, bulk-composition.md and atmosphere-colour.md) in a PR, with a dated comment on the goal issue. On 2026-10-17 at 12:02 AM EDT Boss widened the environment's network access, so the sources the first check could not reach can be re-read.

## 4a. Decided (Boss, 2026-10-17, as the coordinator relayed his card taps)

1. 1:21 AM EDT, "Composition first": new bodies get their composition drawn from the formation distance and the star and the bulk figures derived (route 2); route 1 stays as the inverse for editing and for an imported density (GEN.426, GEN.428, GEN.429).
2. 1:21 AM EDT, "Keep the switch": the water share of the solids beyond the snow line has both values selectable until a source settles it, the Marboeuf-derived 0.14 to 0.39 as the research draft's default and about 0.5 as the convention (GEN.428, GEN.458; the research is GEN.468).
3. 1:22 AM EDT, "Replace both": the neutron-star radius floor and the 10 Gyr main-sequence lifetime are derived from sourced relations (GEN.463, GEN.467, 8.1 with the stellar tables, each a generator epoch change joined to the one rescatter).
4. 1:24 AM EDT, "1 bar": the radius reference level (GEN.427, GEN.429).
5. 1:24 AM EDT, "Asplund 2021": the reference solar abundance set, with the spread among Asplund 2021, Lodders 2009 and Lodders 2019 as the uncertainty on ratios (GEN.435, GEN.439).
6. 1:43 AM EDT, "Weights, Dirichlet where unsure": the proportions of a drawn possibility use the table weights directly where a sourced kernel exists, and a Dirichlet draw around those weights, with the concentration from the known spread, where only a range is known (GEN.472, section 4b).
7. 1:44 AM EDT, "Stand-in only": measured Solar System shapes (the asteroid taxon shares and zone kernel, the Earth rock-type library, and the other content taken from measured Solar System objects) enter a table only as flagged temporary stand-ins, each with a research item that replaces it with a derived or general distribution (GEN.473 to GEN.476, TEST.215, section 4b).
8. 1:47 AM EDT, Boss's design point (typed in the body composition thread): anything an asteroid field or belt can be made of, a planet or moon can be made of too, so the materials and components are one shared table gated by class, and no object kind has a list of its own (GEN.427, GEN.472, GEN.477).
9. 1:51 AM EDT, "Yes" on "Are table weights the odds of picking a possibility, fractions normalised separately?": the weights of the gated rows are the odds of picking a possibility, the fractions or partial pressures inside the drawn row are normalised separately by the item that owns them, and both sums are tested (GEN.472, TEST.214).
10. 3:02 AM EDT, Boss's standing rule (typed in the project chat, relayed by the coordinator): everything verified and partially verified is usable as long as it is not contradictory; if a contradiction occurs the verified source wins. A partly verified value ships in a table as sourced with its stated limit (the validity window or caveat the check named) unless a verified source contradicts it; a contradicted or unsupported value does not ship; where two sources disagree the one that was verified is used.
11. 3:36 AM EDT, defaults of the Tunable gaps investigation ([tunable-gaps.md](tunable-gaps.md) section 6, the coordinator's relay; they stand as Decided unless Boss overrules them): D1 the water-rich share is a switch with the NGPPS logistic as the model default and the flagged bracket 0 to 0.87 (GEN.510); D2 the envelope loss is XUV energy-limited with variable efficiency, a core-mass distribution as a third knob and Lopez and Fortney 2014 as the retained-envelope radius law (GEN.509); D3 the white-dwarf delta-IW is the core-formation scale and the magma-ocean rows the surface scale, with a conversion test (GEN.498, GEN.513); D4 the life-linked gases are flagged design brackets, default off, never tuned (GEN.525); D5 Kzz and haze yield are nuisance parameters with the Brande scatter as prior, the LU[1e7, 1e11] Kzz draw replaced by a Teq threshold (GEN.519); D6 the comets take the Solar System population as the within-system shape plus a per-system offset knob bracketed by the two interstellar objects (GEN.520). The unsourced weights are seeded draws tuned to published populations through one harness (GEN.503, builders GEN.509 to GEN.525); the rows that no public data set defines stay labelled brackets.

## 4b. The gated possibility table (Boss, 2026-10-17, 1:40 AM EDT)

Boss (typed in the project chat): for each composition there is a list of possibilities with weighted probabilities, for all compositions, each category a table of possibilities and their associated values; the list items are gated, each table line lists the planet or object classes it supports, and the composition is drawn from the probabilities, which add up to 100 percent. The coordinator answered yes with two refinements, kept here: a weight may depend on the star and the formation distance (the route 2 decision of section 4a), and every weight carries its source (nothing invented; the standing no-fixed-content rule).

**The standing format** (one loader, GEN.472; the sum test is TEST.214):

| Column | Meaning |
|---|---|
| possibility | The row: a recipe, a taxon, an origin type, a population or a gas, named for what it is and never for a Solar System object |
| supported classes | The class gate: the planet, moon, star or object classes the row applies to |
| weight or weight function | A number, or a function of the star's stored [Fe/H], mass and age, the formation distance and the snow line |
| source and licence | Where the weight comes from and the licence of the data (the F5 rule: the loader refuses a row without them); `open` where none is found |
| notes | The ranges the row's drawn fractions or partial pressures take, the validity range, the check |

**The draw:**

1. Take the rows whose gate passes the body's class.
2. Evaluate each weight at the body's star and formation distance; drop rows whose weight is zero.
3. Normalise the rest to add up to 100 percent and pick one with the one random wrapper (GEN.56), in the table's fixed row order.
4. The drawn row supplies the ranges of its fractions or partial pressures; those are normalised by the item that owns them (the mixture in GEN.411, the bulk fractions in GEN.428). Two sums are therefore checked: the weights of the possibilities, and the fractions of the drawn one.
5. Proportions inside a drawn possibility (Decided 1:43 AM EDT, section 4a item 6): the table weights are used directly where the row names a sourced kernel, and as a Dirichlet draw around the weights where only a range is known, with the concentration from the spread, built from the one random wrapper.

**Rules that follow:**

- No row passes the gate, or every weight is zero: an error that names the class and the table. There is no silent fallback to a fixed row.
- A weight with no source is written `open`; the open rows of one gate share equally, the loader lists them in one report, and each is replaced when a source or a derivation is found. This is the one place the table holds a placeholder, and it is visible.
- Adding, removing or reweighting a row changes the draws of every seed in its classes, so each table has a version and a change is a generator epoch bump.
- A row built from a measured Solar System shape carries the marker `stand-in` (Decided 1:44 AM EDT, section 4a item 7), names its replacement research item (GEN.473 asteroid belt composition, GEN.474 rock types from the bulk oxides, GEN.475 meteorite groups as outputs, GEN.476 Al-26 and Fe-60 from the formation environment), and is listed with the `open` weights in the loader's report; the calibration rows themselves (Earth, Mars, Moon, the Galilean moons, Sirius B) live in test files only and no generator table reads them (TEST.215).
- The materials and components are one shared table (Decided 1:47 AM EDT, section 4a item 8): a belt's, a comet's, a moon's and a planet's components are rows of the table of GEN.427 selected by the class gate, and the loader refuses a component list that belongs to one kind alone; the belt and field draw of today, with its own list of 84 names and per-family lists, is replaced by it (GEN.477).
- One mixture storage and one text writer serve every body (the body composition research lane's section 14.6, proposed for the design GEN.426 to settle): a body's composition is a mixture of rows of the shared table stored as component mass fractions (DB.69 and DB.70 together, replacing `asteroid_belt_composition` and the planet composition text), the gate does the physics (condensation temperature at the formation distance, oxidation state, water and carbon inventory, star abundances) and not the object type, the rare-earth, uranium, platinum-group, gem and aromatic species are trace accessories scaled to the chemistry, and the descriptions come from the stored fractions through one thresholded writer (GEN.430, GEN.446), the class label becoming a summary of the mixture.
- The source column carries the verification state (Decided 3:02 AM EDT, section 4a item 10): `verified`, `partly verified` with the limit the check named (the row ships as sourced-with-limit), `open` (no source), or `stand-in`; a value the checks found contradicted or unsupported is not a row's source, and when a partly verified and a verified source disagree the verified one is used. The loader reports the count of each state per table.
- Tables register themselves with the loader, so TEST.214 walks them all and a table added later is covered without a change to the test.

**Where it applies:** the gas table (GEN.410, gate and source columns), the material table (GEN.427), the star abundance rows (GEN.439, GEN.440), the bulk recipes (GEN.428), the atmosphere recipes (GEN.411, GEN.454 to GEN.456), the small-body and moon tables (GEN.442, GEN.457, GEN.466) and the surface rows (GEN.441). Where an item is derived and not drawn (GEN.412, GEN.458 to GEN.462, GEN.464) it has no table and takes the drawn row or the weight arguments of the tables above. The element table (GEN.438) is a constants table with no weights.

### As built (GEN.472)

- **Where.** `planetgen/generation/possibility_table.py`: `Row`, `Table`, `Registry` (the module's `REGISTRY` is the generator's), `Subject`, `Condition`. The Dirichlet draw is `draw.dirichlet` (and `draw.gamma`) in `util/draw.py`, so it is built on the one random wrapper.
- **The gate evaluator.** A row passes when the subject's class is in its `classes` (or the row says `ANY_CLASS`) and every `Condition` holds. A condition reads one stored quantity of the subject (`minimum`, `maximum` or `equals`, a bound inclusive) and a quantity the subject lacks is an error, never a pass. The condition predicates of section 14.6 (the condensation temperature at the formation distance, the oxidation state, the water and carbon inventory, the star's abundances) are conditions on those stored or derived quantities, so no table needs class letters.
- **Open weights.** The sourced weights of the rows that pass are probabilities and the open rows of that gate share what is left over equally; with only open rows they share all of it. If the sourced weights already use up the whole, an open row gets nothing and stays in the report; if they add up to more than 1 while an open row passes, the draw is refused.
- **Fractions.** A row with `fractions` and no `concentration` gives them directly (scaled to 1); with a `concentration` they are a Dirichlet draw around them (`concentration_from_spread(mean, sd)` gives the concentration for a known spread). Component names are rows of the one shared table (`SHARED_COMPONENTS`, "materials", GEN.427): a table with fractions must name it, and `Registry.validate` refuses a component that table does not hold, so no kind has a list of its own.
- **Versions.** `Table.digest` hashes the version, every row's name, gate, weight (a function by its qualified name), state and fractions; `Registry.versions` lists `(version, digest)` per table. The first real tables (GEN.410, GEN.427, GEN.439) pin them with their epoch bump; the loader itself changes no stored value and has no epoch.
- **Not in the loader.** The sum test over every class and table (TEST.214), the target-table check of the tunable rows (GEN.503) and the tables themselves.

## 5. What Boss may research for the bodies

From the bulk composition sources check (section 5 of [bulk-composition-sources-check.md](bulk-composition-sources-check.md)): the low-pressure ice and silicate densities for small bodies and moons; the stellar Mg, Si, O and C against [Fe/H] (the project stores only [Fe/H]); a cited solar abundance set (Lodders 2003, Asplund 2009); the water inside the snow line; iron oxidation and the light elements in the core; the SCvH or CMS19 tables if any lane integrates an envelope itself.

The research list of the Tunable gaps investigation ([tunable-gaps.md](tunable-gaps.md) section 6) adds the exact papers and data for the six rows that need a named source (a1, a2, b4, b5, n7, k7: GEN.480, GEN.481, GEN.482, GEN.488, GEN.491): Putirka and Xu 2021 Table A1, the NGPPS formation-stage tables, Deng 2020 Fig. 3 data, Hirschmann 2021 full text and Dello Russo 2016, and the nine further lists of that section, with the licences to confirm (GEN.500).

## 6. Open points for the designs

- Which kinds of body draw the density first and fill in the composition, and which derive the density from the composition (GEN.426).
- Whether a star's composition changes the luminosity, lifetime and colour tables of GEN.407 (GEN.435).
- The order of the steps when the ground's colour is part of the albedo and the albedo sets the temperature (GEN.441, GEN.412, GEN.422).
- How a small body's mass follows from its composition and size (belts have a stated mass of zero today, NAV.62).
- Which radiation proxies of `physics/radiation.py` go and in which release (GEN.436, GEN.444).
- The licence of the data the ozone and band cross sections come from (open since the sources check).
