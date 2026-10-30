# The composition model: what every body, star and atmosphere is made of

**Status (2026-10-17):** checked design for Boss's approval (GEN.437, with the notes GEN.409, GEN.426, GEN.435 and GEN.436). The loader of section 4b (GEN.472, #2061), the element table (GEN.438), the stellar relations (GEN.439) and the star composition (GEN.440) are built (see "As built" in sections 4b and 7); the body and atmosphere tables and their draws are not. This is the master note of Boss's directive (2026-10-16, 11:58 PM EDT). It joins [atmosphere-colour.md](atmosphere-colour.md) (the atmosphere and the colour), [bulk-composition.md](bulk-composition.md) (the bulk of every body), [stellar-composition.md](stellar-composition.md) (the stars and remnants) and [radiation-model.md](radiation-model.md) (the radiation); section 7 joins the four. Every section the research lanes write goes into the design documents named here.

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
7. 1:43 AM EDT, "Stand-in only": measured Solar System shapes (the asteroid taxon shares and zone kernel, the Earth rock-type library, and the other content taken from measured Solar System objects) enter a table only as flagged temporary stand-ins, each with a research item that replaces it with a derived or general distribution (GEN.473 to GEN.476, TEST.215, section 4b).
8. 1:47 AM EDT, Boss's design point (typed in the body composition thread): anything an asteroid field or belt can be made of, a planet or moon can be made of too, so the materials and components are one shared table gated by class, and no object kind has a list of its own (GEN.427, GEN.472, GEN.477).
9. 1:51 AM EDT, "Yes" on "Are table weights the odds of picking a possibility, fractions normalised separately?": the weights of the gated rows are the odds of picking a possibility, the fractions or partial pressures inside the drawn row are normalised separately by the item that owns them, and both sums are tested (GEN.472, TEST.214).
10. 3:02 AM EDT, Boss's standing rule (typed in the project chat, relayed by the coordinator): everything verified and partially verified is usable as long as it is not contradictory; if a contradiction occurs the verified source wins. A partly verified value ships in a table as sourced with its stated limit (the validity window or caveat the check named) unless a verified source contradicts it; a contradicted or unsupported value does not ship; where two sources disagree the one that was verified is used.
11. 3:36 AM EDT, defaults of the Tunable gaps investigation ([tunable-gaps.md](tunable-gaps.md) section 6, the coordinator's relay; they stand as Decided unless Boss overrules them): D1 the water-rich share is a switch with the NGPPS logistic as the model default and the flagged bracket 0 to 0.87 (GEN.510); D2 the envelope loss is XUV energy-limited with variable efficiency, a core-mass distribution as a third knob and Lopez and Fortney 2014 as the retained-envelope radius law (GEN.509); D3 the white-dwarf delta-IW is the core-formation scale and the magma-ocean rows the surface scale, with a conversion test (GEN.498, GEN.513); D4 the life-linked gases are flagged design brackets, default off, never tuned (GEN.525); D5 Kzz and haze yield are nuisance parameters with the Brande scatter as prior, the LU[1e7, 1e11] Kzz draw replaced by a Teq threshold (GEN.519); D6 the comets take the Solar System population as the within-system shape plus a per-system offset knob bracketed by the two interstellar objects (GEN.520). The unsourced weights are seeded draws tuned to published populations through one harness (GEN.503, builders GEN.509 to GEN.525); the rows that no public data set defines stay labelled brackets.
12. 12:39 PM EDT, Boss's general rules (typed in the project chat at 12:39 PM EDT on 2026-10-17, relayed by the coordinator): when in doubt, verified wins; when two verified sources are in doubt, the newest wins; when two sources are both in question, the one that takes the most permission wins (his words are 'most permission take wins'; Boss confirmed by card at 1:12 PM EDT on 2026-10-17, relayed in words by the coordinator, that it means the widest range of the two sources). Boss answered a card at 12:40 PM EDT (relayed in words): his six earlier 'widen and keep both' answers stand (arm pitch against bulge fraction, nuclear cluster mass slope, blue to red cluster radius ratio, thick-to-thin height ratio, thick-to-thin mass ratio against Vc, cluster system radius), and the newest-wins rule applies to every conflict that comes after them. This extends item 10: it settles what item 10 left open when two verified sources disagree.

## 4b. The gated possibility table (Boss, 2026-10-17, 1:40 AM EDT)

Boss (typed in the project chat): for each composition there is a list of possibilities with weighted probabilities, for all compositions, each category a table of possibilities and their associated values; the list items are gated, each table line lists the planet or object classes it supports, and the composition is drawn from the probabilities, which add up to 100 percent. The coordinator answered yes with two refinements, kept here: a weight may depend on the star and the formation distance (the route 2 decision of section 4a), and every weight carries its source (nothing invented; the standing no-fixed-content rule).

**The standing format** (one loader, GEN.472; the sum test is TEST.214):

| Column | Meaning |
|---|---|
| possibility | The row: a recipe, a taxon, an origin type, a population or a gas, named for what it is and never for a Solar System object |
| supported classes | The class gate: the planet, moon, star or object classes the row applies to |
| weight or weight function | A number, or a function of the star's stored [Fe/H], mass and age, the formation distance and the snow line |
| source | Where the weight comes from (the F5 rule: the loader refuses a row without it). The licence is not a column: it is recorded in `DATA-LICENSE.md` in the repository root (Boss, 2026-10-22) |
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
- A row built from a measured Solar System shape carries the marker `stand-in` (Decided 1:43 AM EDT, section 4a item 7), names its replacement research item (GEN.473 asteroid belt composition, GEN.474 rock types from the bulk oxides, GEN.475 meteorite groups as outputs, GEN.476 Al-26 and Fe-60 from the formation environment), and is listed with the `open` weights in the loader's report; the calibration rows themselves (Earth, Mars, Moon, the Galilean moons, Sirius B) live in test files only and no generator table reads them (TEST.215).
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
- **Versions.** `Table.digest` hashes the version, every row's name, gate, weight (a function by its qualified name), state and fractions; `Registry.versions` lists `(version, digest)` per table. The first real tables (GEN.410, GEN.427, GEN.439) pin their table version; the 8.1 tables bump no epoch (Boss, 2026-10-17, 2:10 PM EDT, section 7.5), and the loader itself changes no stored value and has no epoch.
- **Not in the loader.** The sum test over every class and table (TEST.214), the target-table check of the tunable rows (GEN.503) and the tables themselves.

## 5. What Boss may research for the bodies

From the bulk composition sources check (section 5 of [bulk-composition-sources-check.md](bulk-composition-sources-check.md)): the low-pressure ice and silicate densities for small bodies and moons; the stellar Mg, Si, O and C against [Fe/H] (the project stores only [Fe/H]); a cited solar abundance set (Lodders 2003, Asplund 2009); the water inside the snow line; iron oxidation and the light elements in the core; the SCvH or CMS19 tables if any lane integrates an envelope itself.

The research list of the Tunable gaps investigation ([tunable-gaps.md](tunable-gaps.md) section 6) adds the exact papers and data for the six rows that need a named source (a1, a2, b4, b5, n7, k7: GEN.480, GEN.481, GEN.482, GEN.488, GEN.491): Putirka and Xu 2021 Table A1, the NGPPS formation-stage tables, Deng 2020 Fig. 3 data, Hirschmann 2021 full text and Dello Russo 2016, and the nine further lists of that section, (no licence gate, Boss 2:35 PM EDT).

## 6. Open points for the designs

- Decided (Boss, 2026-10-17, 1:21 AM EDT, "Composition first"; section 7.3): new bodies get their composition drawn first and the bulk figures derived (route 2), so this point is closed (GEN.426).
- Whether a star's composition changes the luminosity, lifetime and colour tables of GEN.407 (GEN.435).
- The order of the steps when the ground's colour is part of the albedo and the albedo sets the temperature (GEN.441, GEN.412, GEN.422).
- How a small body's mass follows from its composition and size (belts have a stated mass of zero today, NAV.62).
- Which radiation proxies of `physics/radiation.py` go and in which release (GEN.436, GEN.444).
- The licence of the data the ozone and band cross sections come from (open since the sources check).

## 7. The checked model (GEN.437)

The master note joins the four checked designs and the loader. This section is the joint: what links to what, in what units, which way each kind of body is drawn, what is stored, when the generator epoch moves, and which note an item reads. The numbers and the worked examples are in the four notes; this section repeats none of them.

| Note | Item | Settles |
|---|---|---|
| [atmosphere-colour.md](atmosphere-colour.md) | GEN.409 | The gas and compound table, the mixture, albedo and greenhouse (stored), ozone, aerosols and cloud, the giants' reference level, the colour model |
| [bulk-composition.md](bulk-composition.md) | GEN.426 | The bulk fractions, the equation of state and the reference solver, the draws, density, mass and gravity, moons and small bodies, the inverse and the migration |
| [stellar-composition.md](stellar-composition.md) | GEN.435 | X, Y, Z and the elements of every star from the stored [Fe/H], age and population, the remnants, what the star's iron, magnesium and silicon give the planets |
| [radiation-model.md](radiation-model.md) | GEN.436 | The bands, transmission and greenhouse from the infrared opacity, shielding by column, the ground's dose and radon, the surface ultraviolet |
| the loader, `generation/possibility_table.py` | GEN.472 | The one gated possibility-table format of section 4b |

### 7.1 The chain, as checked

1. **Elements** (GEN.438): atomic data and the radioactive isotopes. A constants table, no weights, no abundances.
   As built (GEN.438): `physics/elements.py` loads three data files in `physics/data/` (`elements.json`, `isotopes.json`, `decay_series.json`), built by `scripts/build_element_tables.py` from CIAAW's Atomic Weights 2021 (weights of the 118 elements, the abridged conventional value for the interval elements), NIST's isotope masses and natural compositions, NUBASE2020 (half-lives, decay modes and branchings, and the AME2020 mass excesses the chain steps' Q values come from), Lugaro et al. 2018 Table 4 (the energy released per decay of 23 radionuclides, the whole chain for U-238, U-235 and Th-232) and Dye 2012 Table 2 (the four long-lived parents' heat constants, kept as the cross-check). Each file states its sources and the day it was read, and the loader refuses a file that does not (the licence is not a column: it is recorded in `DATA-LICENSE.md`). `molar_mass(formula)` is the one place a compound's molar mass is computed (kg/mol); `heat_per_kg_of_isotope`, `specific_activity` and `isotopic_mass_fraction` give the constants the radiogenic heat (GEN.459) and the ground's dose read. Where the sources differ the newer verified one is used: NUBASE2020's half-lives (U-238 4.463 Gyr against Dye's 4.468) and Lugaro's energies (K-40 1.05e-6 erg against Dye's 1.10e-6, so the heat of K-40 is 2.2 percent below Dye's table; U and Th are within 0.3 percent). The test sums the molar mass of every stored atmosphere species from its elements; the gas and materials tables of GEN.410 and GEN.427 add their compounds to that check when they land. Finding for the owner of `habitability_world.MOLAR_KG`: its sulfur dioxide (0.06407) and methane (0.01604) are 1.9e-4 above the element sums (0.064058, 0.016043), from older weights.
   As built (GEN.410, first slice): `physics/gases.py` loads `physics/data/gases.json`, built by `scripts/build_gas_table.py` from the tables of research parts 01, 02 and 03 (`docs/design/composition/atmosphere-composition-parts/`): 78 gases and compounds (the 72 of part 01 T1 and the six it adds in section 5), each with formula, name (spelled as the code spells it), CAS number, molar mass computed from `physics.elements` (isotope masses for 3He, HD, HDO and radon; the research part's own figure is kept and the loader refuses a row more than 0.01 g/mol away), the infrared flag (monatomic none, homonuclear diatomics collision-induced only, the rest infrared), the volatility class, the photolysis threshold and channel, the near-ultraviolet flag, triple, boiling and critical points with the tag of where each came from, the acentric factor, ideal-gas heat capacity at 100, 300 and 1000 K, gamma, latent heats, and for 24 gases the Rayleigh cross section at 450, 550 and 650 nm with the King factor and the reference temperature its refractivity fit was made at (N2, Ar, CO2, CH4, N2O 288.15 K; O2, CO, H2, NH3 273.15 K; none for the polarizability route and for water vapor, whose refractivity has no known reference temperature, so it takes the generic value). CO takes the Allen/Cox form, CH4 the He 2021 dispersion and water vapor the vapor depolarization. Each row has `source`, `state` and, when not verified, the `limit` naming the block the research did not reach (17 rows have a phase-point, heat-capacity or identity block not reached; 13 Rayleigh entries are second-hand, assumed or on the polarizability route, good to -9 to +6 percent). It is a constants table with no weights and no class gate; the recipes of GEN.411 and GEN.454 to GEN.456 cite its rows by formula and add the class gate. The test `test_gas_table.py` fails when a gas named in `docs/design/research/imported/atmospheric-toxicity.md` is not a row, or when the text `atmosphere.describe` writes names a gas that is not a row (GEN.415: no planet class carries an atmosphere text any more, only `air`, whether it has air; every gas name in the generator's text is the gas table's, `atmosphere.gas_name` and `atmosphere.species_name`, one spelling, "water vapor"). Second slice, Henry constants and viscosity: 35 gases carry `henry` (the NIST WebBook median kH at 298.15 K in mol/(kg bar) with its range and count, the temperature coefficient d, the independent `thermo` fit for the audit, a `dissociates` flag for CO2, SO2, NH3, H2S, HCl and HCN whose constant is the physical solubility only, and `Henry.kh(T)`, which refuses a temperature outside 273.15 to 323.15 K where the research gives the law); 38 carry `viscosity` (Lennard-Jones sigma and epsilon/k, the NIST value at 300 K where one was read, and `Gas.viscosity_upa_s(T)` by Chapman-Enskog with the Neufeld integral, which refuses kT/epsilon outside 0.3 to 100, and refuses water vapor, which is polar, 10 percent off and has only its NIST 9.77 uPa s at 300 K). The coefficient part 01 printed (26.693) is micropoise: the code uses 2.6693 for uPa s, and a test reproduces all 37 non-polar viscosities of the research table to 0.01. A Henry entry is partly verified when it rests on one compiled value, the compiled values span more than a factor 3, or the NIST median and the thermo fit differ by more than a factor 2; a viscosity is partly verified when no NIST value was read or it is more than 5 percent from it. Third slice, infrared bands and ultraviolet bins: 37 gases carry `infrared`, the band centres of part 03 table 2.1 in cm-1 (`Band.wavelength_um` gives the wavelength; a band the research writes as a range, such as the CH3 modes of methanol, has `range_cm1` and no centre; `approximate` marks a value the research writes with "about", which it recalled and did not read), the rank of the gas in the Byrne and Goldblatt 2014b saturated-absorption ordering where the research states one (H2O 11, CO2 16, CH4 27, C2H2 6, SO2 8, CH3OH 1) and the research's own band and strength cells kept whole, because they carry the forcing statements, overlaps and broadening that are not numbers. A gas with no fundamentals (the noble gases, N2, O2, H2, the metals and oxides, HNO3) has an empty band list and the research's reason in the row; a row of recalled band positions is partly verified with that limit. 22 gases carry `uv`, the mean total photoabsorption cross section in cm2 per molecule over the 35 bins `gases.UV_BINS_NM` (10 to 1100 nm, part 02 section 2.3), `None` where the bin lies outside the data (the source package treats it as zero, and `UvAbsorption.at` returns 0.0), a flag where under half of the bin is covered, the wavelengths the data cover and where they come from (Heays 2017 or Huebner 2015 unless the research names another). The research's H2CO is the table's CH2O. The builder checks that every typed band number is in the research cell and that every bin is missing exactly where the coverage says it lies outside the data. The limits the research names are kept on the entry: the structured Schumann-Runge bands of O2 (175 to 200 nm), where a bin mean understates the transmission of a thick column by orders of magnitude and the generator needs a k-distribution or effective cross section, and the visible CH4 bins, which rest on Karkoschka 1994, a paper the research did not reach. `UvAbsorption.longest_wavelength_above(x)` gives the upper edge of the last bin at or above a cross section. The bin means are inputs for the photolysis and shielding items (GEN.413, GEN.436); nothing reads them yet, so there is no generator epoch. Still to build: the saturation-pressure curves (GEN.452, parts 01 sections 6 and 7, the composite of Wagner, Fray and Schmitt and the laboratory fits) and the class gate (GEN.411).
2. **Gases and materials** (GEN.410, GEN.427): the two tables every body reads. The materials table is the one shared components table of section 4a item 8; the gas table is the atmosphere's. Molar masses come from the elements. Every row is a row of the loader (section 4b) with a class gate, a source and a verification state. As built (GEN.427): `physics/materials.py` loads `physics/data/materials.json` (158 rows), built by `scripts/build_materials_table.py` from part E's materials table, part B's SLB phases and iron rows, part C's ice phases and the values read from Seager et al. 2007 Tables 1 to 3, the Handbook of Mineralogy and PubChem (`scripts/materials_data.py`). Every number has a source key and part E's label code, and a value with no source is `null` (the recalled cells of part E are dropped, not kept). Rows are gated `any` or `planetary`, equations of state are kept as each source gives them (Vinet, Birch-Murnaghan, SLB Mie-Grueneisen-Debye, the modified polytrope, tables), and the 84 asteroid components, the field families, the comet list and the class composition texts all resolve to rows (a test fails otherwise). The existing lists in `tuning` are not rewired to it yet: that rides GEN.428 and GEN.430, so no draw changes.
3. **Stars** (GEN.439, GEN.440): composition is a pure function of the stored [Fe/H], age and population (population from the age gate where a star has none). It adds no column to the 75 million stars.
   As built (GEN.439): `physics/stellar_abundances.py` holds the relations of docs/design/composition/stellar-composition.md as pure functions with no draw in them: the Asplund 2021 scale moved to proto-solar (`PROTO_MASS_FRACTION`, the derived helium-to-metal ratio 1.674, built from the element table's atomic weights), the age gate and the Bensby logistic for the branch (`BRANCH_TABLE`, a row pair of the GEN.472 loader that registers itself, the high-alpha row first so one fixed uniform number decides it and the branch flips at most once as [Fe/H] falls), the hinge fits of [X/Fe] with the Cholesky-correlated residuals and the carbon-enhanced share (`abundance_offsets`, which takes the normals and uniforms as arguments so GEN.440 can key them to the star), X, Y, Z and the element fractions (`composition`), the Fe/Si, Mg/Si, C/O and O/Si ratios, the first dredge-up of Boothroyd and Sackmann, the birth Eu, Th and U for GEN.459, and the remnants: the white dwarf's core type, the neutron star's constants, the black hole's quantities and the brown dwarf's limits and radius (GEN.590 removed the Chandrasekhar white dwarf table and the Chamel and Haensel outer-crust list, which the downloaded Bedard et al. sequences and Servais and Chamel crust files replaced). The fit constants are in `tuning.py` with their sources; `RELATIONS` lists every relation with its source and verification state, and the four that stay `open` or `stand-in` (neon, the neutron star's inner crust and core, and the brown dwarf's radius and hydrogen limit against metallicity) wait for a source or for an integrated equation of state. As built (GEN.590): `physics/stellar_composition.py` reads the downloaded files and requires them (Boss's 'Require the data'): a main-sequence star's MIST track and lifetime, a white dwarf's radius and layers at its temperature, a neutron star's outer crust; bulk generation refuses to start without them (gap-filling-plan.md 6.6). Nothing stored changes, so there is no epoch.
   As built (GEN.440): `physics/stellar_composition.py` is the derived function. `star_makeup(kind, fe_h, age_gy, mass_sol, key)` returns the birth and surface composition, the branch, the star's draws and what is particular to the kind (the dredge-up of a red giant, a brown dwarf's limits and radius, a white dwarf's core, atmosphere and radius, a neutron star's outer crust, a black hole's horizon quantities, a Wolf-Rayet star's surface fractions), and names what the design leaves open for that kind. The draws are a stream of their own, `abundance:{key}`, drawn in a fixed order and all always drawn, so the same star has the same composition at any worker count, an edited [Fe/H], age or mass re-derives it, and a group's members (a system's primary, a cluster's seed) share the branch and offsets and take their own phase. `kind_of` finds the kind from the stored Yerkes class (a Wolf-Rayet star's "W" included, GEN.831: its subtype "WN6", "WC4" or "WO2" picks the row of `tuning.STELLAR_WOLF_RAYET`), the phase and the mass. No column, no stored value and no edit to `evolve_star`: the test sweeps the kinds `evolve_star` makes and checks each against the relations.
4. **Bodies** (GEN.428 to GEN.430, GEN.441, GEN.442, GEN.457 to GEN.466): the bulk fractions, the surface and the small bodies. The formation distance, the snow line and the star's abundances decide the gate.
   As built (GEN.458): `formation.snow_line_au(luminosity)` is the 182 K blackbody position (`condensate_shares.ICE_CONDENSATION_K`, Lodders 2003; 2.34 AU for the Sun since GEN.950, 2.68 AU at the 170 K it first used; the fixed 2.7 AU constant is gone) and `formation.draw_snow_line_au` draws a star's own line around it, 0.2 dex in log10 cut at 3 sigma (`tuning.SNOW_LINE_SCATTER_DEX`, `SNOW_LINE_SCATTER_MAX_SIGMA`), from a stream of its own (`snow_line:<mass>:<age>:<galactic distance>`) so no other draw moves. Every star and compact remnant stores it (`Star.snow_line_au`, `stars.snow_line_km`, schema v132); a close pair's merged proxy takes the square root of the sum of the squares of its stars' lines, derived from the stored pair and not stored. The disk walk that counts planets and the disk's outer edge read the stored line. `Planet.born_beyond_snow_line` (`planets` and `moons` columns) is `distance >= the host's line`, set once the orbits are final, a moon taking its planet's; belts and comets read the stored line when GEN.442 lands. The generator epoch bump is the 8.4 chain's (Boss, fork 1), not this change's.
5. **Atmospheres** (GEN.411, GEN.412, GEN.454 to GEN.456): partial pressures first, then the molar mass, the pressure and the scale height from them. As built (GEN.478): `physics/atmosphere_recipes.py` loads `physics/data/atmosphere_recipes.json`, built by `scripts/build_atmosphere_recipes.py` from part 09 of the atmosphere research (its words) and `scripts/atmosphere_recipes_data.py` (the machine columns): the 104 rows of its 16 categories as the tables `atmosphere-a` to `atmosphere-p` of the loader, each row with its class gate, conditions on stored quantities, a weight expression, its source with the verification label, its source state and, for an open row, the research item that replaces it (R1 to R15 are GEN.480 to GEN.494). The rows that compete for 100 percent are 14 groups (a5 to a7, b1 to b4, c1 to c4, c5 and c6, d1 to d5, e1 to e4, f1 to f4, g, h, i, j, k, l, p); a group is drawn once, a level-2 row is gated on `drawn.<group>` so the joint probability is parent times conditional, and where every weight is zero the group's stored fallback row wins (c6, d1, e1, l2, or the explicit "none" of j and k). The other rows are presence probabilities (trace gases, clouds by temperature), continuous quantities as a parametric draw and as a binned table (normal, log-normal, uniform, log-uniform or a mixture), derived values, the cloud windows a profile crosses and one choice (the rotation track). A DESIGN bracket is a seeded constant drawn in the file's order (or the bracket's midpoint where nothing was drawn), never a fixed number. The expressions are Python arithmetic over named stored quantities (`INPUTS` lists them, with their units) and a few functions; the loader parses them with `ast` and refuses anything else, a name it does not know, a circle of `w()` or `p()` references, a probability above 1, a class it does not know, or a row without a source and a state. The corrections of the research are in: the shoreline of d is the Berta-Thompson 2025 logistic (the K window fitted to 13 Solar System bodies is gone and the four atmosphere regimes carry its probability by design edges), a3 and a4 follow the 22 polluted white dwarfs, f7 has the explicit core (centre 0.57, stopped at 0.8) and the planet tail of 0.043 to 0.087 (part 05's recipe put 13.9 percent of giants above C/O 0.8), the lava-world boundaries are the text edges (the caption list differs by 63 K), p4 takes the saturation time from the stellar mass, and j2 to j4 are outputs of the retention physics (Boss, fork 7), so 8 rows are SOURCED, 51 DERIVED and 45 DESIGN where the part had 8, 48 and 48. Of the 43 open rows two are stand-ins (a2, the law fitted to the Moon, Mars and Earth, and k7, one comet's species list). Nothing reads the tables yet, so no draw changes and there is no epoch of its own: GEN.411, GEN.413, GEN.414, GEN.451 and GEN.454 to GEN.456 take the one bump of the chain.
   As built (GEN.495): `physics/comet_population.py` loads `physics/data/comet_population.json` and draws a comet's coma as percent of H2O: the carbon oxides first, S = (CO + CO2)/H2O log-normal (median 18 percent, 0.2 dex since GEN.629) and f = CO/(CO + CO2) logit-normal around the class median (0.08 Jupiter-family and Halley-type, 0.5 long-period, logit scatter 1.5), CO = S f and CO2 = S (1 - f), redrawn together while either lies outside the observed bounds (CO 0.15 to 30, CO2 2 to 37; at most 1000 tries, then `PopulationError`, never a clipped value); then CH4, C2H6, NH3, HCN, CH3OH, H2CO and H2S as independent truncated log-normals with the geometric mean and scatter of part G1 section 2.3, one random number each by inverse transform; then the carbon-chain-depleted flag from the registered table `comet-carbon-chain` (a pair of gated rows per period class: 0.37, 0.25, 0.185). An interstellar comet draws CO/H2O log-uniform in 0.01 to 2 and CO2/H2O log-uniform in 0.02 to 10 (the object measured at 3.32 AU is far outside the Solar System range) and the other species from the population. Each number has its source, label and state in the file; 8 of the 9 species rows are verified; the partly verified ones (H2S from nine comets, the 0.2 dex width of S, the class medians of f, the flag fractions with the Halley-type one, the interstellar edges) name their limit. 67P, 103P/Hartley 2 and Hale-Bopp are test targets only (in `test_comet_population.py`: inside the bounds and where the draw has density), never a centre. The pooled CO and CO2 over the three classes have geometric means 2.5 and 12.4 percent and scatter 0.54 and 0.22 dex against the table's 2.6 and 11 and 0.64 and 0.27. **GEN.629:** the first build gave S a width of 0.1 dex, the +/- 4 of HP22's median converted, which is narrow against the spread of single comets (67P's bulk 7.8 percent was 3.6 sigma below the median). The width is now 0.2 dex, computed from the Robinson 2024 database (`df_comp_all.csv`, fetched 2026-10-23): the 9 comets with both CO and CO2 measured within 2.5 AU give 0.197 dex (their own median is 15 percent), the 8 within 2.0 AU give 0.138; the widest is kept (Boss's rule for sources in question), the table's median 18 stays, and 67P lies at 1.8 sigma. Nothing stored yet reads the draw (GEN.442, GEN.455, GEN.457 take it), so there is no generator epoch. It supersedes the comet rows k4 to k7 of `atmosphere_recipes`.
6. **Derived values** (GEN.413 to GEN.419, GEN.422, GEN.423, GEN.443 to GEN.446): the colour, the radiation, habitability and the text. Each is a pure function of stored values and the tables.

### 7.2 One set of units and one source field

- Storage is SI: kg, m, s, K, Pa, mol; kg/mol for molar mass (the `atm_molar_density` column keeps its name and unit, decided in the atmosphere note), and the units the existing columns already use stay. Doses are effective dose in mSv per year.
- A solid or a star is described by **mass fractions** (X, Y, Z for a star; component fractions of the shared materials table for a body); a gas by **partial pressures** (the mole fractions follow). Conversion between them uses the element table's molar masses and nothing else.
- Every number that comes from outside the repository sits in a table row or a named constant that carries its **source and its verification state** (`verified`, `partly verified` with the limit the check named, `open` with the seeded range used meanwhile, or `stand-in` with its replacing item). A number in a note has the same marks. The loader refuses a table row without them. The licence is not a column of any table (Boss, 2026-10-22): it is recorded in `DATA-LICENSE.md` in the repository root, and Boss manages it. Boss ruled at 2:35 PM EDT (GEN.500 removed) that every paper and dataset he acquired was obtained legally and is usable: the raw papers and downloaded datasets stay out of the public repository and are never redistributed or sold, the research results and derived values are fine to publish, and the program stays free forever.
- Where sources disagree the rules of section 4a item 12 decide, and each design note lists every case with the rule applied.

### 7.3 Drawn first, or the composition is the cause

Boss (11:58 PM EDT on 2026-10-16): either works, the design picks the simplest consistent one per kind of body. The choice, as the four notes make it:

| Kind | Composition is | Density, radius, gravity |
|---|---|---|
| Star (every class) | derived from the stored [Fe/H], age and population | luminosity, radius and the rest stay as drawn today; composition does not change the GEN.407 tables (the tip luminosity by metallicity is GEN.502) |
| Remnant (white dwarf, neutron star, black hole) | derived from the stored mass and the progenitor relation; a black hole carries none | the white-dwarf radius is integrated; the neutron-star radius and the main-sequence lifetime are derived (GEN.463, GEN.467) |
| Rocky, water and envelope planets, moons, small bodies | drawn first from the gated tables, at the formation distance, the snow line and the star's abundances (route 2) | derived from the fractions by the reference solver or its closed forms at 1 bar; the inverse (route 1) serves editing, an imported density and the migration |
| Giants and brown dwarfs | heavy-element mass from the mass law, the rest envelope | from the tabulated models; the reference level is the lower of 1 bar and the cloud top |
| Atmospheres | partial pressures first, from the gated recipes | pressure is their sum; molar mass and scale height follow; albedo and greenhouse are stored draws, the temperature follows |

### 7.4 Stored against derived

The four notes each hold their table. The rule they share: store a value only when it was drawn, cannot be derived, or is expensive and exact to recompute; derive everything else on demand and test that a cache equals the recomputation. A recalculation after an orbit validation or a zone change never redraws a stored composition (finding F10).

### 7.5 The generator epoch and the migrations

The designs and the loader change no stored value, so they bump no epoch. The build items do, in the groups the notes name:

- **8.1:** the loader (GEN.472) and the data tables (GEN.410, GEN.427, GEN.438, GEN.439, GEN.440) are building blocks with no behaviour change. GEN.463 and GEN.467 keep their own epochs, joined to the one rescatter. The epoch scope of GEN.427 is Decided (Boss, 2:10 PM EDT, card "Data table, no bump", relayed in words): it is a data table that bumps nothing, and the one bump is GEN.428 to GEN.430 plus DB.69 in 8.4 phase 1.
- **8.4:** one bump for the bulk chain (GEN.428 to GEN.430, DB.69 and the formation items GEN.458 to GEN.460, with the table versions); one for the atmosphere values (GEN.411, GEN.422, GEN.423, GEN.424, GEN.452, GEN.454, GEN.455 with DB.67 and DB.71); one for the added stored values (GEN.413, GEN.414, GEN.456 with DB.68). The revisions take the next free numbers at build time.
- **9.0:** the colour functions (GEN.416 to GEN.419) are pure functions of stored values and need no migration. The radiation items (GEN.443, GEN.444) replace the dose proxies with a bump and a migration of the dose columns, after GEN.453 has corrected the constants.

A change to a table row or weight is a change of that table's version and an epoch bump (section 4b); `Registry.versions` pins each table's version and digest.

### 7.6 How the items that follow read it

| Items | Read |
|---|---|
| GEN.438, GEN.410, GEN.427, GEN.472 | section 4b, then the gas table of the atmosphere note and the shared materials table of the bulk note |
| GEN.439, GEN.440, GEN.463, GEN.467, GEN.502 | the stellar note |
| GEN.428 to GEN.430, GEN.441, GEN.442, GEN.457 to GEN.466 | the bulk note |
| GEN.411 to GEN.419, GEN.422 to GEN.424, GEN.452 to GEN.456 | the atmosphere note |
| GEN.443 to GEN.446, GEN.453, GEN.459 | the radiation note, with the bulk note for the ground |
| TEST.211, TEST.212, TEST.214, TEST.215, GEN.503 | each note's test targets; the registry walk of TEST.214; the harness of GEN.503 |

### 7.7 What Boss rules on

Each note ends with its own list of decisions for Boss with a recommendation. The ones that cross notes:

1. **The reading of the third conflict rule.** Decided by Boss at 1:12 PM EDT: "most permission take wins" means the widest range of the two sources; both values are kept where a range cannot be widened. Nothing to approve.
2. **The age gate for the thick disc.** The stellar note moves the gate that replaces a stored population from 3.0 Gy to 8.0 Gy (range 6.5 to 9.5) on the verified thick-disc timing, which bears on the old-disc high-alpha share and on GEN.134.
3. **The epoch scope of GEN.427:** Decided (Boss, 2:10 PM EDT, card "Data table, no bump"). Nothing to approve.
4. **The ozone spread and the thin-air closure** (atmosphere and radiation notes), which the research leaves wider than a single seeded range.
5. **The giants' heavy-element law:** the newest (Chachan 2025) as the default with Thorngren selectable, by the newest-wins rule.
