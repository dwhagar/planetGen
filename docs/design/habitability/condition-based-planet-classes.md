# Condition-based planet and moon classes (GEN.645 to GEN.648)

Status: partly built. The classifier, its rule file, the shadow sweep and the class pages are built (2026-10-23, the class lane); generation still draws a class first until the cut-over (GEN.646, GEN.647), which waits on the composition lanes. Lane D's draft with the 13 places that read the class today is `/mnt/project-files/design-drafts/condition-based-classes-2026-10-23.md`.

## Decisions (Boss, 2026-10-23, times EDT)

- 6:09 PM, GEN.644 option 3: a class is not drawn. Once a planet or moon is complete, its class is worked out from its air (pressure and mixture, none included), its mean surface temperature, its terrestrial or gas composition. The reason he gave for the old way: he did not know how to get the compositions right.
- 6:23 PM: classes are also classed by physical properties (size, surface gravity).
- 6:24 to 6:30 PM: keep the class letters as close to today's as possible, M above all (it is culturally identifiable); never retire a letter without asking him, and think of a replacement class that makes more sense first; letters already free stay in reserve for new designations. **Core classes: H, P, O, J, M, I, C, D and N.**
- 12:39 AM, 2026-10-24 (class card, option 1): keep all 20 letters A to V, redefined by conditions, with K, L, Q, S and V kept; R, U, W, X, Y and Z stay in reserve. This is what the rule file already does.
- 12:39 AM, 2026-10-24 (class B card, option 1): a stripped giant core keeps some hydrogen and helium. The class B template changes (cut-over, C13): its text "a thin mix of helium, sodium, and oxygen" gains hydrogen, and its "occasionally the stripped core of a former gas giant" origin is drawn into the body's mixture by the atmosphere chain (lanes J and L: GEN.411/412/454) as a retained H/He share, not a fixed string. Class I's note about hot-Neptune stripping points at the same B origin. The classifier rule for B (under 0.5 g, with air, above 210 K, no sealed ocean) does not change.
- Habitable is not a class: it is the overall habitability colour Blue or Green (GEN.643).
- Testing mode: new epochs, database versions and API versions are fine.

## The rule file

`src/planetgen/physics/data/planet_classes.json`, read by `physics/planet_classes.py`. `classify(body)` is a pure function of stored values: the first row of the body's kind (terrestrial or giant) whose conditions all hold names the class. The last row of each kind has no conditions, so a body always gets one letter. All 20 letters stay; none is retired. Size and gravity are rule inputs, so there is no size suffix on the stored letter.

| Letter | Class | Rule (design thresholds, refitted by the sweep) |
|---|---|---|
| T | Gas dwarf | giant under 10 Earth masses |
| I | Ice giant | giant under 35 Earth masses |
| J | Gas giant | the rest of the giants |
| D | Small icy body | pressure under 10 Pa and water 1 percent of the mass or more |
| S | Airless super-Earth | pressure under 10 Pa and 1.5 Earth masses or more |
| C | Dead world | pressure under 10 Pa |
| A | Volcanic world | sulfur dioxide 3 percent of the air or more |
| N | Hot dense world | pressure of 1,000,000 Pa (10 bar) or more |
| Q | Variable-climate world | starlight varies 25 percent or more over the orbit |
| V | High-gravity world | gravity of 1.25 g or more |
| L | Marginal thin-air world | pressure under 3 kPa, oxygen 10 percent of the air or more, gravity 0.5 g or more |
| K | Thin-air cold world | pressure under 3 kPa and below 273.15 K |
| P | Glaciated world | mean surface under 210 K |
| B | Small hot world | gravity under 0.5 g |
| E | Wet reducing world | oxygen under 3 percent, methane and hydrogen 5 percent or more, water vapour 15 percent or more |
| F | Reducing world | oxygen under 3 percent, methane and hydrogen 5 percent or more |
| G | Carbon dioxide world | carbon dioxide 20 percent of the air or more |
| O | Ocean world | water covers 90 percent or more of the surface |
| H | Desert world | water covers under 15 percent |
| M | Temperate world | the rest (air, water over 15 to 90 percent) |

The rows are in the order of the file (giants, then D S C A N Q V L K P B E F G O H M); the order decides when two rows hold. Each row carries a one-sentence reading and its `basis` ("design: ..."). A number the file marks design is a planetGen default, refitted when the compositions land; none is a measured constant and none is read from the Solar System.

**K, L, Q, S and V** looked like letters to retire (zone-drawn guesses that PHI-6 now measures). Each instead became a condition that keeps its old meaning: K the cold thin-air desert, L the thin oxygen-bearing marginal world, Q the world whose starlight swings by a quarter or more (an eccentric orbit), S the airless super-Earth (size), V the high-gravity world. Free letters (the reserve): R, U, W, X, Y, Z.

**Known weak spots.** The mass of an ice giant (10 to 60 Earth masses) overlaps a gas giant's (30 and up) in the old draw, so the mass cut at 35 misplaces about a quarter of the drawn ice giants until the envelope share of GEN.428 is stored; the gas dwarf T mostly merges into I; and a drawn Q cannot be reproduced because it was a draw with no stored cause, so the orbit now decides it. The shadow sweep (below) measures all three.

## The shadow sweep

`python src/tests/class_shadow_cli.py` generates systems (or `--per-class N` planets of every class and zone) in memory, classifies every planet and moon, and prints the share that keeps its drawn letter, the confusion table and the share of each letter before and after. First fit (2026-10-23): 98.3 percent of 32,487 bodies from 2,000 natural systems and 96.6 percent of 14,199 forced bodies keep their letter; M keeps 98 to 99 percent; the core classes H, P, O, C, D and N keep 95 to 100 percent. The report is `/mnt/project-files/reports/class-shadow-sweep-2026-10-23.md`.

## What is built and what is next

Built (class lane): C1 the classifier and rule file, C2 the sweep, C4 the class pages from the rule rows (a rule and a core flag in place of the table of ranges), C5 a test that pins the modules that may assign a class (none of these changes a stored value). C3, life off the letter, is built in its own change with a generator epoch: a body carries a light-harvesting pigment when its stored conditions allow it (each `LIFE_CHEMICALS` entry's `requires`: free oxygen for the oxygenic pigments, none for the anoxygenic, an acidic ocean for the zinc one, air for every pigment but the radiotroph) and the star's spectrum favours it; the most advanced life stage is capped only by the habitability index (`life.stage_ceiling`), not by the letter; the flavour line reads whether the body carries life and the derived class A. Next: C13 to C15, the cut-over, when lanes I, J, K and L have landed the composition inputs: classify at the end of generation, the class column written only there, validation becomes "total and derived", the zone-drawn tables (`PLANET_CLASSES` generation fields, `PLANET_CLASS_PROBABILITIES`, `MOON_BLACKLIST`) deleted, no old-letter alias kept.
