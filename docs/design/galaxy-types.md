# Galaxy types: one family per known galaxy type

**Status:** plan. Boss asked for it on 2026-10-16 at 8:48 PM EDT for version 9.4; nothing is built. GEN.378 turns this plan into the checked design with Boss's approval of its table before any build.

Boss (2026-10-16, 8:48 PM EDT): "For v9.4 or so, I want to add a feature that we have different galactic models for generation based on known galaxy types. TODO should break that apart into smaller chunks in a chain of required items, add those to our TODO plan, count its foundational steps back from 9.4 so they are ready for when we need to actually implement it."

## 1. What it means under the standing rules

- **No fixed content** (Boss, 2026-10-16, 4:11 and 4:14 PM EDT). A galaxy type is a parameter family, never a copy of a named galaxy: every number a family needs is calculated from the galaxy's stellar mass and a few seeded draws through published scaling relations, with a draw only where the science gives a distribution. The real galaxies (the Milky Way, Andromeda, M87 and so on) appear only as reference points for plausibility tests, never as a preset or an option.
- **One rule set.** The galaxy today is one spiral family with its numbers still copied from the Milky Way (`build_galaxy_shape` and the `tuning` defaults); GEN.368 plans the draw that replaces them. A type is the next level up: the draw first picks the family, then the family's laws fill the shape. The spiral family is GEN.368's draw, not a second implementation.
- **Reuse, not a second model.** The density model (`galaxy/density.py`), the potential (GEN.338, 8.3), the halo (GEN.234), the satellites (GEN.369), the clusters (GEN.297), the scatters (GEN.301) and the sector profile already take a galaxy shape. A family changes the shape's contents, so the chain below makes the shape carry its components and leaves the readers alone.
- **Multiple galaxies stays the existing "later" note.** A database still holds one galaxy; the family is a property of that galaxy. A second galaxy with its own family is the later feature (multiple-galaxies.md).

## 2. The families

| Family | What it has | What it leaves out |
|---|---|---|
| Spiral (unbarred and barred, early to late) | thin and thick disk, bulge, arms, bar by a draw, halo | |
| Lenticular | disk and bulge, no arms, old stars, almost no gas | arms, young stars |
| Elliptical | one Sérsic spheroid with an axis ratio, old stars, a halo | disk, arms, bar, young stars |
| Irregular and dwarf | a clumpy exponential or Sérsic body with young stars, or a Plummer dwarf spheroid | arms, bar, a defined bulge |

The draw (GEN.380) takes the galaxy's stellar mass from a Schechter function with the seed, then the family from the mass-dependent morphology fractions. No family is the default; the old shape defaults go with GEN.368.

## 3. The chain, counted back from 9.4

Each row is an item in TODO.md; a row sits in the earliest version that can hold it and be useful on its own, and everything that only makes sense with the families stays in 9.4. The build order rule applies: a need is earlier in its thread, in an earlier phase or in an earlier version.

| Step | Item | Version | Why there |
|---|---|---|---|
| 1 | GEN.378 the checked design: families, the scaling relation table with sources, the draw order, Boss's research requests | 8.1 phase 5 | needs only GEN.368's table; the shape-draw build that follows GEN.368 reads it, so it is not built twice |
| 2 | DB.66 the shape row carries the family key and the component amplitudes (migration only) | 8.3 phase 1 | one migration beside the other galactic-gravity work; Core 1 owns the numbering |
| 3 | GEN.379 the density model as a list of components, each with an amplitude that may be 0, and `component_masses` | 8.3 phase 1 | the potential (GEN.338) reads each component's mass from it, so it comes first; today's spiral reproduces today's density |
| 4 | GEN.380 the type draw; GEN.381 the Sérsic spheroid component; UX.135 one hand-drawn vector outline per family | 9.4 phase 1 | nothing calls them before the families exist; the outlines can be drawn any time after GEN.378 fixes the family list |
| 5 | GEN.382 spiral, GEN.383 elliptical, GEN.385 irregular and dwarf | 9.4 phase 2 | three independent threads |
| 6 | GEN.384 lenticular; GEN.386 populations; GEN.387 clusters, nuclear cluster and black hole; GEN.388 home position | 9.4 phase 3 | they read the families |
| 7 | GEN.389 creation draws the family and stores it | 9.4 phase 4 | needs every family |
| 8 | UX.134 the Generate page, MAP.242 the Galaxy Map, TEST.204 plausibility tests, DOC.59 docs | 9.4 phase 5 | four independent threads |

GEN.338 gains one need (GEN.379) and one test: a component of zero mass is legal, so an elliptical (no disk) needs no change in the potential.

## 4. The laws each family needs, to be checked in GEN.378

These are the relations the design looks up. The names are from memory and are not verified; GEN.378's job is to check each against the paper, state its range and scatter, and ask Boss for research where a relation cannot be found.

| Parameter | Relation | Source named from memory |
|---|---|---|
| Stellar mass | Schechter function | Baldry et al. 2012; Li and White 2009 |
| Family fractions by mass | morphology-mass relation | Conselice 2006; Kelvin et al. 2014 |
| Halo mass | stellar-to-halo mass relation | Behroozi et al. 2013; Moster et al. 2013 |
| Disk scale length | mass-size relation by family | Shen et al. 2003; van der Wel et al. 2014; Lange et al. 2015 |
| Disk scale height | thickness against length | Kregel et al. 2002 |
| Bulge-to-total, Sérsic index | Hubble stage and mass | Graham and Worley 2008; Graham and Driver 2005 |
| Bar fraction | mass and gas | Masters et al. 2011; Sheth et al. 2008 |
| Arm count and pitch angle | shear and mass | Seigar et al. 2006; Yu and Ho 2019 |
| Spheroid size, axis ratio | mass-size and triaxiality by mass | Shen et al. 2003; Weijmans et al. 2014 |
| Gas fraction, star-formation history | mass and family | Catinella et al. 2018 and the main-sequence relations |
| Central black hole | M-bulge and M-sigma | Kormendy and Ho 2013 |
| Cluster counts | specific frequency by family | Harris et al. 2013; Peng et al. 2008 |
| Dwarf profiles | size-mass, metallicity-mass | McConnachie 2012; Kirby et al. 2013 |

## 5. Defaults standing as decided

- The family is drawn from the seed; the Generate page may pin it ("auto", spiral, barred spiral, lenticular, elliptical, irregular). Pinning fixes the family, never a size or a count: the laws still draw those from the seed.
- Barred or unbarred is a draw inside the spiral family, not a family of its own.
- A satellite galaxy is made of the dwarf families (GEN.369 already draws profile by mass); a Plummer dwarf spheroid and an exponential dwarf irregular are the same objects the satellite laws draw.
- The root draw is the stellar mass; the halo mass follows from it by the stellar-to-halo relation, and GEN.368's table is read from that root.
- Each family has a simple hand-drawn static SVG outline in the icon library (UX.135, Boss, 2026-10-16, 9:00 PM EDT), shown beside the family's name on the Generate page and in the Galaxy Map key, so the user sees roughly what the model looks like. Static, kept in the repo, no generation step, works on the offline page.
- Interacting and merging galaxies are not part of it (idea S9.4-a in the roadmap).
