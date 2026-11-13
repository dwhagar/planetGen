# Condition-based planet and moon classes (GEN.645 to GEN.648)

Status: built (the class cut-over, GEN.646, 2026-10-28, the class lane). A body is made from its mass, where it formed and its star; its class letter is derived from what it comes out as, at the end of its build and again once its star's flux variation is known. Nothing draws a class first. Lane D's draft with the 13 places that read the class before the cut-over is `/mnt/project-files/design-drafts/condition-based-classes-2026-10-23.md`.

## Decisions (Boss, 2026-10-23, times EDT)

- 6:09 PM, GEN.644 option 3: a class is not drawn. Once a planet or moon is complete, its class is worked out from its air (pressure and mixture, none included), its mean surface temperature, its terrestrial or gas composition. The reason he gave for the old way: he did not know how to get the compositions right.
- 6:23 PM: classes are also classed by physical properties (size, surface gravity).
- 6:24 to 6:30 PM: keep the class letters as close to today's as possible, M above all (it is culturally identifiable); never retire a letter without asking him, and think of a replacement class that makes more sense first; letters already free stay in reserve for new designations. **Core classes: H, P, O, J, M, I, C, D and N.**
- 12:39 AM, 2026-10-24 (class card, option 1): keep all 20 letters A to V, redefined by conditions, with K, L, Q, S and V kept; R, U, W, X, Y and Z stay in reserve. This is what the rule file already does.
- 1:22 AM, 2026-10-24 (class B card, option 3; this replaces the 12:39 AM option 1 answer relayed earlier): class B never has hydrogen or helium, and a stripped giant core must be at least 5 Earth masses. At the class cut-over (C13) class B's text drops "helium" and the stripped-core line; the classifier rule for B (under 0.5 g, with air, above 210 K, no sealed ocean) does not change. Class I's note about hot-Neptune stripping is reworded to match.
- Habitable is not a class: it is the overall habitability colour Blue or Green (GEN.643).
- 10:21 AM EST on 2026-11-08 (class M rarity report, lever 9, option keep): the carbon dioxide limit for classes M, H, O, P and R stays at 1 kPa (judged by pressure, GEN.1321); class L's edge is unchanged. Raising it was not chosen: among 19 oxygen sea worlds a 2 kPa limit lets 5 pass and 10 kPa lets 8, but their temperature and ocean cover fail M anyway, so the lever is not where the rarity comes from.
- Testing mode: new epochs, database versions and API versions are fine.

## The rule file

`src/planetgen/physics/data/planet_classes.json`, read by `physics/planet_classes.py`. `classify(body)` is a pure function of stored values: since GEN.1207 the class is the best fit. Every class has its own conditions in `physics/data/planet_class_fit.json` (a condition just outside its edge lowers the match instead of ruling the class out) and the class with the highest percent match names the body, so no class relies on another having said no; a body that fits none is shown with its nearest class at a match of 0. The rows of `planet_classes.json` give only each class's name, phrase, sentence and basis. 25 letters are in use (Q is the one reserved letter; Q was retired by Boss, GEN.1135). Size and gravity are rule inputs, so there is no size suffix on the stored letter.

| Letter | Class | Rule (design thresholds, refitted by the sweep) |
|---|---|---|
| T | Gas dwarf | giant under 10 Earth masses |
| I | Ice giant | giant under 35 Earth masses |
| J | Gas giant | the rest of the giants |
| D | Small airless world | pressure under 10 Pa, round, and under 0.1 Earth masses (GEN.1160) |
| S | Airless super-Earth | pressure under 10 Pa and 1.5 Earth masses or more |
| C | Dead world | pressure under 10 Pa and 0.1 Earth masses or more |
| A | Volcanic world | sulfur dioxide 3 percent of the air or more |
| T | Gas dwarf (terrestrial kind) | hydrogen 30 percent of the air or more and under 10 Earth masses (GEN.1142; the same T as the giant row) |
| N | Hot dense world | pressure of 10 bar or more, 400 K or more, hydrogen under 30 percent of the air (GEN.1142) |
| V | High-gravity world | gravity of 2 g or more (GEN.1136; 1.25 g before) |
| L | Marginal thin-air world | pressure 30 kPa to 3 bar, 273.15 to 320 K, carbon dioxide 1 kPa or more (GEN.1321: a pressure too), oxygen 3 to 16 kPa (GEN.1321: judged by oxygen pressure; a share of 3 to 16 percent before), gravity 0.5 g or more (GEN.1142) |
| K | Thin-air world | pressure under 3 kPa, 273.15 to 320 K, any oxygen share (GEN.1142; the under-10-percent oxygen cap removed in GEN.1306, registry AW4) |
| E | Deep-frozen world | air, mean surface under 233.15 K (-40 C) (GEN.1142; the letter of the old wet reducing world) |
| P | Glaciated world | the Earth-like cold world: mean surface 240 to 270 K with a melting region in season, near-standard air with oxygen (GEN.1142, GEN.1221; the temperate family's too-cold member) |
| W | Frozen world | air, mean surface 233.15 K (-40 C) to freezing (273.15 K), under 2 g, not a hydrogen envelope and not a volcanic air; the cold air world class P does not hold (GEN.1306 part 4, registry AW8, Boss 11-06). P wins a tie, E is colder |
| Z | Hot thin-air world | no ocean, 0.5 to 2 g, 340 to 1500 K, air 1 to 10 kPa (GEN.1211; 1.4 g and 1100 K before GEN.1306, registry AW6 and AW7); the 1 and 10 kPa edges are the retention caps of the air |
| B | Small hot world | gravity under 0.5 g, air, and a mean surface above freezing (273.15 K); the name says hot but the class holds every temperature above freezing, up to 2,000 K in a census |
| F | Anoxic world | oxygen under 3 kPa (GEN.1321: a pressure; 3 percent of the air before) and carbon dioxide under 20 percent, 0.5 to 1.4 g (GEN.1281; the reducing world before, which needed methane and hydrogen 5 percent or more) |
| G | Carbon dioxide world | carbon dioxide 20 percent of the air or more, any oxygen share (GEN.1306, registry AW5) |
| O | Ocean world | liquid water covers 95 percent or more of the surface (GEN.1142; 80 percent before GEN.1570, the cover is the liquid sea since GEN.1567), at a mean 275 to 312 K: M's floor up to the old ceiling, so a cool global ocean is O and not left to the shoulders (GEN.1342) |
| H | Desert world | some water, but under 20 percent of the surface (GEN.1142), at a mean 275 to 315 K: M's floor up to the old ceiling; a world with no water at all is not H (GEN.1342) |
| M | Temperate world | the rest (air, liquid water 20 to 95 percent); M, H, O, P and R reach 2 g, V's line (GEN.1306, registry AW3; 1.25 g before); M, H, O, P and R judge oxygen (16 to 35.5 kPa) and carbon dioxide (under 1 kPa) by pressure, not share (GEN.1321, registry OP1) |
| R | Hothouse world | M's too hot (GEN.1221): mean 312 to 345 K, hottest region 315 K or more, coldest 295 K or more, wet bulb 31 C or more, liquid water phase, and the oxygen, air, gravity and carbon dioxide ranges of M |
| U | Methane-ethane world | methane-ethane liquid stable at the pressure and temperature (91 to 130 K, 50 to 500 kPa), methane 1 percent or more of the air, 0.05 to 2 g (GEN.1221) |
| Y | Sulfuric acid world | sulfuric acid liquid stable at the pressure and temperature (300 to 350 K), pressure 3 kPa or more, sulfur dioxide 3 percent or more of the air, 0.5 to 2 g (GEN.1221) |
| X | Ammonia-water world | no air; a liquid cover of 20 percent or more, ammonia 10 percent or more of the water's mass, 190 to 240 K, 0.05 to 2 g (GEN.1258) |

The table is in the order of the file. The order no longer decides anything (the best fit does; a tie goes to the class with more conditions, then to the letter), and R, U, X, Y and Z, which the rule file did not have when this table was written, are listed after M with their conditions as `planet_class_fit.json` holds them. Each row carries a one-sentence reading and its `basis` ("design: ..."). A number the file marks design is a planetGen default, refitted when the compositions land; none is a measured constant and none is read from the Solar System.

**K, L, S and V** looked like letters to retire (zone-drawn guesses that PHI-6 now measures). Each instead became a condition that keeps its old meaning: K the cold thin-air world (below freezing, usually under ice), L the thin oxygen-bearing marginal world, S the airless super-Earth (size), V the high-gravity world (2 g or more since GEN.1136; at 1.25 g it took most temperate super-Earths and left M almost empty). **Q is retired** (Boss, 6:07 PM EDT 2026-10-30, GEN.1135): a class that was only an orbit. The swing of the starlight is the stored `flux_variation` (about 4 times the eccentricity) and the "At a glance" panel of a planet or moon flags it as very high from 0.75 (`describe.sheet.SWING_FLAG_MIN`, 4 percent of bodies). 19 classes are in use. The only reserved letter is Q (R, U, W, X, Y and Z became classes in GEN.1221, GEN.1258 and GEN.1306).

**GEN.1142 (Boss, 6:57 to 7:12 PM EDT 2026-10-30): the lines move, conditions only, no zone.** The old wet reducing E (F with steam, no body in a sample of 1,788) is retired and its letter goes to a new class, the deep-frozen world: air-bearing bodies under 233.15 K (-40 C). P takes 233.15 K to freezing. K is the thin air (under 3 kPa) of 273.15 to about 320 K with little oxygen, L the 30 kPa to 3 bar air of 273.15 to about 320 K with carbon dioxide and 3 to 16 percent oxygen (GEN.1321 made this and the other oxygen edges pressures: 3 to 16 kPa), O is 80 percent ocean or more, H under 20 percent, M the rest. A terrestrial-kind body under 10 Earth masses with 30 percent hydrogen or more in its air is a T (the gas dwarf, one letter in both kinds) and is tested ahead of N, which is a 10 bar air at 400 K or more that is not hydrogen-dominated. The 10 Earth-mass line is the bound chain's critical core mass (`GIANT_CORE_CRITICAL_MASS_EARTH`). Registry rows CL1 to CL7 (GEN.725) hold the lines and their tiers.

**Known weak spots (as written in October; see the class description audit of 2026-11-12 for the text against code check).** The mass of an ice giant (10 to 60 Earth masses) overlaps a gas giant's (30 and up), so the mass cut at 35 misplaces some ice giants until the envelope share of GEN.428 is stored (Aux 17, GEN.822); the gas dwarf T stays rare and mostly merges into I until then; and a dense steam air (the solver of GEN.411 puts the water of an ice-rich world into the air; GEN.830 and GEN.837 of Aux 8 address it) makes class N a large share of the warm planets until they land. The class sweep (below) measures them.

## The class sweep

`python src/tests/class_sweep_cli.py` generates systems in memory, prints the share of every letter among planets and among moons, the mass bins of the planets, and how many bodies change letter when `classify` is run again on their stored conditions (none should). The thresholds of `planet_classes.json` and the mass law (`planet-mass-law.md`) are tuned against it and against Aux 16's check (`/mnt/project-files/reports/planet-class-distribution-check-2026-10-25.md`). The cut-over's own report is `/mnt/project-files/reports/class-cutover-sweep-2026-10-28.md`. Before the cut-over this section described a shadow sweep that compared `classify` with the drawn letter (98.3 percent of 32,487 bodies kept their letter, 2026-10-23); with no drawn letter there is nothing to compare, and the sweep now reports the distribution instead.

## How a body is placed when it fits no class (GEN.1314)

A body that meets every condition of a class takes that class and scores 1. Otherwise the fit looks for the class that needs the fewest values moved, and moves them only as far as the stand-ins allow. The width of each class edge is `physics/data/fudge_table.json` (GEN.1313): the root sum of squares of what the registered stand-ins move that quantity by, widened by the confidence of each (registry 6.261, FT1 to FT4). The hand shoulders that stood beside every edge are gone.

* A value may be moved by at most 0.25 width; if no class fits that way, by at most 1 width (registry CP1). The smaller reach is tried first.
* Among classes reachable at the same reach, the one with the fewest values moved wins (the temperatures of a body count as one value), then the one with the smallest total move, then the one with more conditions, then the letter.
* The match is `1 - total move / (values moved * 1 width)` (registry CP2); a body that fits fully scores 1.
* A body no class reaches by 1 width has its gravity and air-amount edges relaxed (digit 9), and if that still fails it takes the nearest class with a match of 0.
* The planet sheet says which values were moved and by how much (`class_fudge_values`, `class_fudge_amount`).
* What reads a class reads its match with it (GEN.1315, `physics/class_outcome`, registry FT5): a request for a class is met, and a Class M world is counted, only by a body whose `class_match` is at least `tuning.CLASS_MATCH_THRESHOLD` (0.5); the search takes a range on `class_match` (`planet_min_class_match=0.5`). The habitable count and the bar on habitable worlds never read a letter: they read the PHI-6 score and the life of the finished body (GEN.1387, GEN.1397).

## What is built and what is next

Built (class lane): C1 the classifier and rule file, C2 the sweep, C4 the class pages from the rule rows (a rule and a core flag in place of the table of ranges), C5 a test that pins the modules that may assign a class (none of these changes a stored value). C3, life off the letter, is built in its own change with a generator epoch: a body carries a light-harvesting pigment when its stored conditions allow it (each `LIFE_CHEMICALS` entry's `requires`: free oxygen for the oxygenic pigments, none for the anoxygenic, an acidic ocean for the zinc one, air for every pigment but the radiotroph) and the star's spectrum favours it; the most advanced life stage is capped only by the habitability index (`life.stage_ceiling`), not by the letter; the flavour line reads whether the body carries life and the derived class A. Next: C13 to C15, the cut-over, when lanes I, J, K and L have landed the composition inputs: classify at the end of generation, the class column written only there, validation becomes "total and derived", the zone-drawn tables (`PLANET_CLASSES` generation fields, `PLANET_CLASS_PROBABILITIES`, `MOON_BLACKLIST`) deleted, no old-letter alias kept.

## Age needs from conditions (GEN.686, removed by GEN.1546)

A star's age used to be raised until it was old enough for its planets (a tier of `age_needs` read from each finished body's conditions). That floor is gone: the star's age is drawn once, with the habitable-request tilt of `sample_living_star` (`LIFE_MIN_STAR_AGE_GY`) and the life timeline caps kept, and the planets are made on it. The class stays an output of the finished body, and nothing reads a class to set an age.

## A sample that holds every class (GEN.791, planets built)

`SystemConfig.FORCE_CLASS` (a test and admin option, not saved with a system) gives a generated system at least one planet of the named class in the first orbital slot; `generation/class_sample.py` `fill(per_class, seed)` tops a sample up until every class has that many planets, for checks that need all twenty (the property sheet fill audit, UX.208). Random generation gives under 30 planets of ten classes in 400 systems. Moons are next, after Core 3's moon generation work (GEN.761) merges. At the class cut-over the switch becomes a condition request (an ocean world, PHI-6 Blue or Green) through the same path as `+habitable_world`.

