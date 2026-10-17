# The composition model: what every body, star and atmosphere is made of

**Status (2026-10-17):** plan; nothing is built. This is the master note of Boss's directive (2026-10-16, 11:58 PM EDT). It joins [atmosphere-colour.md](atmosphere-colour.md) (the atmosphere and the colour), the bulk composition of section 5 of that note, and the new stellar and radiation parts. GEN.437 turns it into the checked design once the research lanes and the four designs below are done. Every section the research lanes write goes into the design documents named here.

## 1. What Boss decided

- Atmosphere composition and the composition of planets, moons, stars and everything is the joint headline of 9.0, together with the atmosphere colour.
- Two research lanes first, one for atmospheres and one for bodies (planets, moons, stars, remnants and the rest): the compounds, minerals and elements that make them up, what is there and all the statistics needed to use it in calculations. All notes go into the design documents.
- Every body gets figures on its composition from actual science, and not made up where that can be avoided. Drawing the density and the other figures that are drawn today and then filling in a composition that fits is acceptable, and so is building from the most basic composition up; either works, and the design picks the simplest consistent one per kind of body.
- 9.0 completes it, together with the atmosphere colour, and adds what becomes derivable: the atmospheric radiation (the atmosphere's transmission and greenhouse effect, since we know what is in the air) and the surface radiation (from the atmosphere and the ground composition).
- The building blocks go back from 9.0 to 8.1 (as early as their blockers allow), then a collision and consolidation sweep merges items that do almost the same thing and suit one worker, then the issue clusters and all documents are updated.
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
| 8.1 phase 1 | GEN.431 | Research: what is in the atmospheres of every kind of world, and the numbers each compound needs |
| 8.1 phase 1 | GEN.432 | Research: what rocky, icy and giant bodies and their small relatives are made of, and the numbers each material needs |
| 8.1 phase 1 | GEN.433 | Research: what stars and stellar remnants are made of, and the numbers each needs |
| 8.1 phase 1 | GEN.434 | Research: radiation through atmospheres and at the ground, and the numbers it needs |
| done | GEN.425 | Research: verify the mass-radius and equation-of-state sources and say what is missing ([bulk-composition-sources-check.md](bulk-composition-sources-check.md)) |

**Designs (docs/design)**

| Where | Item | What |
|---|---|---|
| 8.1 phase 1 | GEN.409 | The atmosphere colour design: the checked plan for the gas table, the mixture, the stored values, ozone, aerosols and the colour model |
| 8.1 phase 1 | GEN.426 | The bulk composition design: the fractions, the equation of state, the draws and the derived density, mass and gravity |
| 8.1 phase 1 | GEN.435 | The stellar composition design: the abundance relations, what is stored and what is derived, and the remnants |
| 8.1 phase 1 | GEN.436 | The radiation design: atmospheric transmission and greenhouse by band, shielding by column, and the ground's own dose |
| 8.1 phase 1 | GEN.437 | The composition model: how the element, compound, material and body tables link, and what is stored against derived |

**Data tables (no behaviour change)**

| Where | Item | What |
|---|---|---|
| 8.1 phase 2 | GEN.438 | One table of the elements and their isotopes: atomic data and the radioactive isotopes' half-lives and decay energies |
| 8.1 phase 2 | GEN.410 | One table of the gases and compounds that can be in an atmosphere, with every compound the descriptions mention |
| 8.1 phase 2 | GEN.427 | One table of the bulk materials with their equation-of-state constants and the source of each |
| 8.1 phase 2 | GEN.439 | The stellar composition relations: abundances by population, age and [Fe/H], and the remnants' layers |

**Star composition (derived, no column)**

| Where | Item | What |
|---|---|---|
| 8.1 phase 3 | GEN.440 | Star composition is derived from the stored [Fe/H], age and population, for every star, remnant and brown dwarf |
| 8.1 phase 4 | TEST.211 | A star composition sweep: every kind of star, remnant and brown dwarf, fractions, trends and determinism |

**Bodies: generation and storage**

| Where | Item | What |
|---|---|---|
| 8.4 phase 1 | GEN.428 | Planets and moons draw their bulk composition fractions instead of carrying a fixed text |
| 8.4 phase 1 | GEN.429 | Bulk density, mass, radius and gravity are derived from the composition fractions with the mass-radius relation |
| 8.4 phase 1 | DB.69 | Planets and moons store the composition fractions |
| 8.4 phase 1 | GEN.430 | The planet and moon descriptions state the composition from the stored fractions |
| 8.4 phase 1 | GEN.412 | Generation records the albedo and the greenhouse multiplier it drew, and the temperature is exactly what they give |
| 8.4 phase 2 | GEN.441 | Surface and crust composition: the minerals, ices and ocean on the ground, from the bulk fractions, the temperature and the air |
| 8.4 phase 2 | GEN.442 | Asteroids, comets, rings and dwarf planets draw their composition from the formation distance and the snow line |
| 8.4 phase 3 | DB.70 | Planets, moons and small bodies store the surface and small-body composition shares |

**Atmosphere: mixture, ozone, aerosols**

| Where | Item | What |
|---|---|---|
| 8.4 phase 2 | GEN.411 | The atmosphere's molar density comes from its per-gas partial pressures, so mixture, pressure and density always agree |
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
| 9.0 phase 4 | DOC.66 | Docs: the composition model of stars, bodies and atmospheres, the radiation it derives and the sources of every number |
| 9.0 phase 4 | GEN.421 | Headline of 9.0: what every body and atmosphere is made of, and the atmosphere colour that follows from it |

## 4. The research lanes

- **Atmosphere research (2026-10-17):** GEN.431, with the re-check of the four sources the first check could not reach.
- **Body composition research (2026-10-17):** GEN.432 (planets, moons, small bodies, giants), then GEN.433 (stars and remnants). GEN.425 (mass-radius and equation-of-state sources) is done, see the sources check.
- **Radiation research:** GEN.434 is covered by the two lanes above (the coordinator, 12:17 AM EDT); no separate thread unless Boss asks.

Each lane writes its draft to /mnt/project-files/design-drafts/ and reports; the TODO thread merges it into the design documents (atmosphere-composition and body-composition sections, stellar-composition.md, radiation-model.md, bulk-composition.md and atmosphere-colour.md) in a PR, with a dated comment on the goal issue. On 2026-10-17 at 12:02 AM EDT Boss widened the environment's network access, so the sources the first check could not reach can be re-read.

## 5. What Boss may research for the bodies

From the bulk composition sources check (section 5 of [bulk-composition-sources-check.md](bulk-composition-sources-check.md)): the low-pressure ice and silicate densities for small bodies and moons; the stellar Mg, Si, O and C against [Fe/H] (the project stores only [Fe/H]); a cited solar abundance set (Lodders 2003, Asplund 2009); the water inside the snow line; iron oxidation and the light elements in the core; the SCvH or CMS19 tables if any lane integrates an envelope itself.

## 6. Open points for the designs

- Which kinds of body draw the density first and fill in the composition, and which derive the density from the composition (GEN.426).
- Whether a star's composition changes the luminosity, lifetime and colour tables of GEN.407 (GEN.435).
- The order of the steps when the ground's colour is part of the albedo and the albedo sets the temperature (GEN.441, GEN.412, GEN.422).
- How a small body's mass follows from its composition and size (belts have a stated mass of zero today, NAV.62).
- Which radiation proxies of `physics/radiation.py` go and in which release (GEN.436, GEN.444).
- The licence of the data the ozone and band cross sections come from (open since the sources check).
