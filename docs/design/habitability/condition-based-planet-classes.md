# Condition-based planet and moon classes (GEN.645 to GEN.648)

Status: built (the class cut-over, GEN.646, 2026-10-28, the class lane). A body is made from its mass, where it formed and its star; its class letter is derived from what it comes out as, at the end of its build and again once its star's flux variation is known. Nothing draws a class first. Lane D's draft with the 13 places that read the class before the cut-over is `/mnt/project-files/design-drafts/condition-based-classes-2026-10-23.md`.

## Decisions (Boss, 2026-10-23, times EDT)

- 6:09 PM, GEN.644 option 3: a class is not drawn. Once a planet or moon is complete, its class is worked out from its air (pressure and mixture, none included), its mean surface temperature, its terrestrial or gas composition. The reason he gave for the old way: he did not know how to get the compositions right.
- 6:23 PM: classes are also classed by physical properties (size, surface gravity).
- 6:24 to 6:30 PM: keep the class letters as close to today's as possible, M above all (it is culturally identifiable); never retire a letter without asking him, and think of a replacement class that makes more sense first; letters already free stay in reserve for new designations. **Core classes: H, P, O, J, M, I, C, D and N.**
- 12:39 AM, 2026-10-24 (class card, option 1): keep all 20 letters A to V, redefined by conditions, with K, L, Q, S and V kept; R, U, W, X, Y and Z stay in reserve. This is what the rule file already does.
- 1:22 AM, 2026-10-24 (class B card, option 3; this replaces the 12:39 AM option 1 answer relayed earlier): class B never has hydrogen or helium, and a stripped giant core must be at least 5 Earth masses. At the class cut-over (C13) class B's text drops "helium" and the stripped-core line; the classifier rule for B (under 0.5 g, with air, above 210 K, no sealed ocean) does not change. Class I's note about hot-Neptune stripping is reworded to match.
- Habitable is not a class: it is the overall habitability colour Blue or Green (GEN.643).
- Testing mode: new epochs, database versions and API versions are fine.

## The rule file

`src/planetgen/physics/data/planet_classes.json`, read by `physics/planet_classes.py`. `classify(body)` is a pure function of stored values: the first row of the body's kind (terrestrial or giant) whose conditions all hold names the class. The last row of each kind has no conditions, so a body always gets one letter. 19 letters are in use (Q was retired by Boss, GEN.1135); the rest are reserved. Size and gravity are rule inputs, so there is no size suffix on the stored letter.

| Letter | Class | Rule (design thresholds, refitted by the sweep) |
|---|---|---|
| T | Gas dwarf | giant under 10 Earth masses |
| I | Ice giant | giant under 35 Earth masses |
| J | Gas giant | the rest of the giants |
| D | Small icy body | pressure under 10 Pa and water 1 percent of the mass or more |
| S | Airless super-Earth | pressure under 10 Pa and 1.5 Earth masses or more |
| C | Dead world | pressure under 10 Pa |
| A | Volcanic world | sulfur dioxide 3 percent of the air or more |
| T | Gas dwarf (terrestrial kind) | hydrogen 30 percent of the air or more and under 10 Earth masses (GEN.1142; the same T as the giant row) |
| N | Hot dense world | pressure of 10 bar or more, 400 K or more, hydrogen under 30 percent of the air (GEN.1142) |
| V | High-gravity world | gravity of 2 g or more (GEN.1136; 1.25 g before) |
| L | Marginal thick-air world | pressure 30 kPa to 3 bar, 273.15 to 320 K, carbon dioxide 1 percent or more, oxygen 3 to 16 percent, gravity 0.5 g or more (GEN.1142) |
| K | Thin-air world | pressure under 3 kPa, 273.15 to 320 K, oxygen under 10 percent (GEN.1142) |
| E | Deep-frozen world | air, mean surface under 233.15 K (-40 C) (GEN.1142; the letter of the old wet reducing world) |
| P | Glaciated world | mean surface 233.15 K to freezing (273.15 K) (GEN.1142) |
| B | Small hot world | gravity under 0.5 g |
| F | Reducing world | oxygen under 3 percent, methane and hydrogen 5 percent or more |
| G | Carbon dioxide world | carbon dioxide 20 percent of the air or more |
| O | Ocean world | water covers 80 percent or more of the surface (GEN.1142) |
| H | Desert world | water covers under 20 percent (GEN.1142) |
| M | Temperate world | the rest (air, water 20 to 80 percent) |

The rows are in the order of the file (giants, then D S C A T N V L K E P B F G O H M); the order decides when two rows hold. Each row carries a one-sentence reading and its `basis` ("design: ..."). A number the file marks design is a planetGen default, refitted when the compositions land; none is a measured constant and none is read from the Solar System.

**K, L, S and V** looked like letters to retire (zone-drawn guesses that PHI-6 now measures). Each instead became a condition that keeps its old meaning: K the cold thin-air world (below freezing, usually under ice), L the thin oxygen-bearing marginal world, S the airless super-Earth (size), V the high-gravity world (2 g or more since GEN.1136; at 1.25 g it took most temperate super-Earths and left M almost empty). **Q is retired** (Boss, 6:07 PM EDT 2026-10-30, GEN.1135): a class that was only an orbit. The swing of the starlight is the stored `flux_variation` (about 4 times the eccentricity) and the "At a glance" panel of a planet or moon flags it as very high from 0.75 (`describe.sheet.SWING_FLAG_MIN`, 4 percent of bodies). 19 classes are in use. Free letters (the reserve): Q, R, U, W, X, Y, Z.

**GEN.1142 (Boss, 6:57 to 7:12 PM EDT 2026-10-30): the lines move, conditions only, no zone.** The old wet reducing E (F with steam, no body in a sample of 1,788) is retired and its letter goes to a new class, the deep-frozen world: air-bearing bodies under 233.15 K (-40 C). P takes 233.15 K to freezing. K is the thin air (under 3 kPa) of 273.15 to about 320 K with little oxygen, L the 30 kPa to 3 bar air of 273.15 to about 320 K with carbon dioxide and 3 to 16 percent oxygen, O is 80 percent ocean or more, H under 20 percent, M the rest. A terrestrial-kind body under 10 Earth masses with 30 percent hydrogen or more in its air is a T (the gas dwarf, one letter in both kinds) and is tested ahead of N, which is a 10 bar air at 400 K or more that is not hydrogen-dominated. The 10 Earth-mass line is the bound chain's critical core mass (`GIANT_CORE_CRITICAL_MASS_EARTH`). Registry rows CL1 to CL7 (GEN.725) hold the lines and their tiers.

**Known weak spots.** The mass of an ice giant (10 to 60 Earth masses) overlaps a gas giant's (30 and up), so the mass cut at 35 misplaces some ice giants until the envelope share of GEN.428 is stored (Aux 17, GEN.822); the gas dwarf T stays rare and mostly merges into I until then; and a dense steam air (the solver of GEN.411 puts the water of an ice-rich world into the air; GEN.830 and GEN.837 of Aux 8 address it) makes class N a large share of the warm planets until they land. The class sweep (below) measures them.

## The class sweep

`python src/tests/class_sweep_cli.py` generates systems in memory, prints the share of every letter among planets and among moons, the mass bins of the planets, and how many bodies change letter when `classify` is run again on their stored conditions (none should). The thresholds of `planet_classes.json` and the mass law (`planet-mass-law.md`) are tuned against it and against Aux 16's check (`/mnt/project-files/reports/planet-class-distribution-check-2026-10-25.md`). The cut-over's own report is `/mnt/project-files/reports/class-cutover-sweep-2026-10-28.md`. Before the cut-over this section described a shadow sweep that compared `classify` with the drawn letter (98.3 percent of 32,487 bodies kept their letter, 2026-10-23); with no drawn letter there is nothing to compare, and the sweep now reports the distribution instead.

## What is built and what is next

Built (class lane): C1 the classifier and rule file, C2 the sweep, C4 the class pages from the rule rows (a rule and a core flag in place of the table of ranges), C5 a test that pins the modules that may assign a class (none of these changes a stored value). C3, life off the letter, is built in its own change with a generator epoch: a body carries a light-harvesting pigment when its stored conditions allow it (each `LIFE_CHEMICALS` entry's `requires`: free oxygen for the oxygenic pigments, none for the anoxygenic, an acidic ocean for the zinc one, air for every pigment but the radiotroph) and the star's spectrum favours it; the most advanced life stage is capped only by the habitability index (`life.stage_ceiling`), not by the letter; the flavour line reads whether the body carries life and the derived class A. Next: C13 to C15, the cut-over, when lanes I, J, K and L have landed the composition inputs: classify at the end of generation, the class column written only there, validation becomes "total and derived", the zone-drawn tables (`PLANET_CLASSES` generation fields, `PLANET_CLASS_PROBABILITIES`, `MOON_BLACKLIST`) deleted, no old-letter alias kept.

## Age needs from conditions (GEN.686, built)

A star's age is raised until it is old enough for its planets. That minimum no longer reads the class letter: `physics/age_needs.py` and `physics/data/age_needs.json` give each finished body a tier from its stored conditions (giants, airless bodies and small worlds need none; hot dense or reducing air is `molten`; heavy worlds, cold thin air and carbon-dioxide air are `settled`; temperate worlds with air are `watered`; cold worlds are `frozen`), and each tier names its minimum age on the fast, normal and slow evolutionary scales. The ages are the old `PLANET_CLASSES` `age_ranges` minimums (planetGen values with no source), so 99.4% of 4,612 sampled bodies need exactly the age they needed before. `age_ranges` is deleted.

## A sample that holds every class (GEN.791, planets built)

`SystemConfig.FORCE_CLASS` (a test and admin option, not saved with a system) gives a generated system at least one planet of the named class in the first orbital slot; `generation/class_sample.py` `fill(per_class, seed)` tops a sample up until every class has that many planets, for checks that need all twenty (the property sheet fill audit, UX.208). Random generation gives under 30 planets of ten classes in 400 systems. Moons are next, after Core 3's moon generation work (GEN.761) merges. At the class cut-over the switch becomes a condition request (an ocean world, PHI-6 Blue or Green) through the same path as `+habitable_world`.

