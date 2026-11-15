# Class M abundance tiers

**Status:** being built (GEN.1704 the setting, GEN.1705 the meter, GEN.1706 the calibration, GEN.1721 the habitable family: the dry part, the desert seas and the quiet interiors are built, section 8; Boss approved the design by question cards, 11:17 to 11:19 PM EST 2026-11-13). Dial values are first-pass stand-ins until the census calibration that follows GEN.1694 to GEN.1698.

**Code:** `planetgen/generation/abundance_tier.py` (the tier and its dials), `tuning.ABUNDANCE_TIER_DIALS` (the values), the `class_m` row of `planetgen/galaxy/dials.py` (the control and its storage), `scripts/class_m_meter.py` (the meter, on the draw of `scripts/class_census.py`).

## 1. Why a tier exists

The generator made no Class M world and no world rated Blue in 2,750 plain systems and 455 Sun-like systems (`reports/class-m-blue-levers-2026-11-14.md`). Part of that was bugs, filed as GEN.1694 to GEN.1698. What is left is not a bug. With the bugs fixed, the generator's own laws give under one Class M world in ten thousand Sun-like systems, for three reasons that the science itself gives:

1. A rocky planet in the ecosphere is usually born far wetter than Earth: the Bern population synthesis (`physics/delivered_water.py`) gives most of them hundreds to thousands of times Earth's water, and sixteen of eighteen measured candidates were 100 percent ocean. Part of that was a unit bug (the law read against the snow line of GEN.1611, about 1.7 times too far out; sent to TODO on 2026-11-13), but the published models agree on the direction: Earth's water class (`low`) is 1 to 8 percent of habitable-zone planets (`reports/water-rich-hz-births-2026-11-14.md`).
2. Thirty-four of seventy-five right-mass candidates keep a hydrogen envelope, which is inside what the published work allows (`reports/hydrogen-envelopes-habitable-zone-2026-11-14.md`).
3. The band of sunlight that gives a Class M surface is narrow: an exact copy of Earth is Class M 23 percent of the time at 1.00 AU, 9 percent at 0.95 AU and 2 percent at 1.10 AU.

Nobody knows how common an Earth is. Earth may be one in a galaxy. So the rate is not a fact the generator can look up; it is a choice of the person who makes the galaxy, and the generator says plainly when that choice goes past the science.

## 2. The tiers

| Tier | Target, Class M worlds per Sun-like system | What it means |
|---|---|---|
| Science (the default) | none | Only the laws. Nothing is moved. |
| Sparse | 1 in 500 | The bottom rung (Boss, 9:28 PM EST 2026-11-14): the laws alone measured 1 in 2,800 on 2,800 Sun-like systems (`reports/train28-census-2026-11-14/meter2800.md`). |
| Scattered | 1 in 135 | |
| Plentiful | 1 in 37 | |
| Abundant | 1 in 10 | The top rung. The two between are split evenly on a log scale from Sparse to Abundant, each step about 3.7 times (Boss, 9:30 PM EST 2026-11-14). The first set was 1 in 1,000, 100, 30 and 10 (cards of 11:19 PM EST 2026-11-13). |

A Sun-like system is one whose primary is a main-sequence star of type F7 to K2 (the sample of the report and of the meter).

The target is what the dial values are calibrated to. The generator never checks a galaxy against it and never adds or removes a world to meet it.

## 3. The six dials

A tier turns six dials and nothing else. Each is a share from 0 to 1.

| Dial | What it moves | Where |
|---|---|---|
| `water_share` | The share of candidates whose water class is given by the tier and not drawn from the Bern law at the distance they formed at: Earth's (`low`, 0.025 to 0.16 percent of the planet's mass) or, for the part `dry_part` names, the dry one. The amount inside the class is still drawn. | `physics/bulk_draw.py`, `physics/delivered_water.py` |
| `dry_part` (GEN.1721) | The part of the planets the water dial reaches that take the dry class (`dry`, under 0.025 percent of the planet's mass: the water of a desert world) and not Earth's. It makes no draw of its own: it splits the water dial's one draw, so the same planets are reached whatever its value. | `physics/delivered_water.py` (`dialed_class`) |
| `quiet_share` (GEN.1721) | The share of candidates drawn with a quiet interior: their volcanic and metamorphic gases take a third of the oxygen a body of the same heat loses (`tuning.ABUNDANCE_QUIET_SINK_FACTOR`), which is the low end of what is published for Earth itself. Nothing else about the interior changes. | `generation/planet.py` (`Planet._draw_body`), `physics/biogenic_oxygen.py` (`reductant_flux`) |
| `cycle_share` (GEN.1720) | The share of candidates whose two carbon-cycle constants, the exponent of the heat flow in the volcanic outgassing (published 1 to 2) and the factor of the weathering law's thermodynamic limit (published a tenth to ten), are drawn from the part of those ranges that holds the carbon balance nearest the Earth's: the exponent from 1 to 1.25 and the factor from 0.1 to 1. The balance goes as the outgassing to about the fifth power, so over the whole ranges a planet with a heat flow of 0.6 or 1.5 times the Earth's holds under 2 Pa or over 100 kPa of CO2 and few hold the 2 Pa to 1 kPa the life classes ask for (`reports/gen1720-co2-attribution-2026-11-14.md`). It moves the draw of two constants and conditions nothing on the CO2 that follows. | `physics/continental_weathering.py` (`draw_cycle`), `physics/surface_air.py` (`Setting.cycle_share`) |
| `bare_share` | The share of candidates that are late cores: they finished growing after the gas disk went and captured no hydrogen envelope (the existing `Planet(late_core=True)`, GEN.1323). | `generation/planet.py`, `physics/envelope_capture.py` |
| `request_share` | The share of systems that are given the request for a habitable world, which multiplies the odds of the wanted draw by `tuning.REQUEST_ODDS` and never forces it (`util/request_odds.py`, GEN.1130). | `generation/system.py` (`StarSystem._draw_tier_request`) |

A **candidate** is a rocky planet, not a moon, whose orbit lies in its star's ecosphere (`abundance_tier.candidate`).

How each dial acts, as built:

- **Water.** One extra draw per candidate. A dialed planet takes the `low` class, or the `dry` class when the draw falls in the dry part, and its air takes the same class's recipe row (b2, or b1 for the dry class). A candidate that formed beyond its star's snow line is icy by origin and is not reached: it keeps its icy water range. A candidate the law gave `low` by itself keeps the air the law draws. The air takes the row only where the pick could take it for that body (its weight in the pick table is above zero); otherwise the pick is drawn as usual.
- **Zone edges.** The water and envelope dials read the zone a planet holds when its body is drawn. A planet that the system check later moves outward across a zone edge takes the new zone's name and nothing is drawn again, so some dialed planets end outside the ecosphere and some planets that end inside it were never candidates. This is not small: measured on 300 Sun-like systems at `abundant` (water share 1), of 276 candidates by final zone 96 were icy by origin (not reached), and of the 180 rocky ones 126 were dialed and 54 were not, because they were drawn in another zone and moved in afterwards. So a water share of 1 reaches about 70 percent of the rocky candidates and 46 percent of all candidates. The meter counts candidates by the final zone; the calibration (GEN.1706) sets the dials by the measured rate, which takes this in, and decides whether the dials should also act when a planet is moved into the ecosphere.
- **Desert worlds.** A planet whose draw falls in the dry part takes its water inside `tuning.ABUNDANCE_DESERT_WATER_RANGE`, the part of the dry class that leaves a few small seas (about 5 to 20 percent of the surface on the planets measured). On any other world so little water holds Earth's salt in so little water that it is saturated with rock salt (GEN.1730: a halite-saturated sea, water activity 0.753 where it is liquid above 252 K; `hydrosphere.halite_sea`). Such a sea of real depth (ten metres or more) is a `surface ocean`: life can start in it (salt-loving oxygen-making algae live at 0.75), and the oxygen law counts it as a sunlit sea of its own area. A tier makes no exception for it: every sea of so little water is the same law. A sea is what gives Class H its groundwater condition as well.
- **Carbon cycle.** One extra draw per candidate, made where the air's design constants are drawn (`surface_air.solve`); a candidate it reaches draws the two constants from the narrower ranges (`tuning.ABUNDANCE_CYCLE_OUTGASSING_RANGE`, `ABUNDANCE_CYCLE_LAMBDA_RANGE`) and records `cycle_nudged` in its design constants; one it does not reach draws over the published ranges. Science makes no draw.
- **Quiet interior.** One extra draw per candidate, made once when the planet is first drawn, beside the envelope draw. The body holds the decision as a field only when the draw says yes.
- **Envelope.** One extra draw per candidate, made once when the planet is first drawn; the decision stays if the orbit later settles into another zone.
- **Request.** One extra draw per system, made by the system itself after its stars are built and before any planet. It acts on the planets only. It is never written into the system's config, so the host star is drawn exactly as without a request (a request in the config thins out white dwarfs from 7.2 to 1.0 percent of systems), nothing is stored, and the system reports no request, met or unmet. A system whose user said yes or no to a habitable world keeps the user's answer.

What a tier may never do (Boss, card "Choose which dials the tiers may turn"):

- change a class limit (the conditions of Class M and the habitability ratings are the same on every tier);
- name a class letter or place a world;
- set or enforce a count.

The class of every body is still derived from its finished conditions. A candidate with Earth's water and no envelope can still freeze, boil, lose its air or keep too much carbon dioxide, and most do.

## 4. Science changes nothing

On Science every dial is 0. No dial makes a random draw, so a seeded system is byte for byte what it was before the tier existed: no generator epoch and no change to the golden systems. Each dial has a test that proves it.

On another tier a dial makes one extra draw per candidate (or per system, for the request) when its share is over 0.

## 5. Where the tier is kept

The tier is the `class_m` row of the galaxy lever table (`galaxy/dials.py`), in the Census dial beside the Rogues level, and it reuses that level's path end to end: the `--class-m` flag of `planetgen plan`, the field on the Generate page, the stored choice in `galaxy_shape_pins`, and the line on the galaxy page. No new table and no migration.

A system gets it through its config. The sector that makes a system sets `SystemConfig.ABUNDANCE_TIER` from the galaxy's shape (`FillContext.apply_population`), as it sets the population, so a system regenerated later is made under the same tier. The tier is not saved with the system's config: it belongs to the galaxy. A system made outside any galaxy is Science unless the caller sets the field. The single-system command takes `--class-m TIER` for that.

Two makers of systems inside a galaxy build their config without the galaxy's shape and so stay on Science: the central star of a planetary nebula (`generation/run_sector.py`) and the systems of a remnant's satellites (`generation/system_pass.py`). Neither holds a candidate in practice. There is no per-run override on the sectors command: a galaxy has one tier.

## 6. Dial values and calibration

| Tier | `water_share` | `bare_share` | `request_share` | `dry_part` | `quiet_share` | `cycle_share` |
|---|---|---|---|---|---|---|
| Sparse | 0.5 | 0 | 0 | 0.2 | 0 | 0.5 |
| Scattered | 1 | 1 | 0 | 0.2 | 0.5 | 1 |
| Plentiful | 1 | 1 | 0.5 | 0.2 | 0.5 | 1 |
| Abundant | 1 | 1 | 1 | 0.2 | 0.5 | 1 |

These are first-pass values from the estimates of the report (every candidate given Earth's water gives about 1 Class M world per 500 Sun-like systems; 1 per 100 needs the envelopes gone as well; the two top tiers need the request). They are registered as stand-ins of tier T9, tweak level L4 (`docs/design/research/gap-filling-plan.md` section 6.289).

Calibration is by the meter: `scripts/class_m_meter.py` draws a seeded sample of Sun-like systems under a tier and prints the funnel (candidates, right mass, ground, gravity, water cover, liquid water, Class M, Blue) and the rate of Class M and Blue worlds per 1,000 systems with its interval, beside the tier's target. The values are set again (GEN.1706) after GEN.1694 to GEN.1698 merge, because those fixes change the rate under every tier.

**A target may be out of reach.** The dials only move draws. If every candidate has Earth's water and no envelope and every system carries the request, and the meter still reads under 1 in 10, then Abundant gives what the dials can give and the page says the measured rate. The remedy is a decision for Boss (a wider dial), never a forced world.

## 7. What is not part of the tier

- The physics fixes of GEN.1694 to GEN.1698 and the weathering and oxygen work: they apply on every tier, Science included.
- The research on why habitable-zone planets are born water-rich and on hydrogen envelopes. If it shows the law itself is wrong, the law is fixed for Science and the dial values shrink.
- The habitable request's own conditions (`tuning.REQUEST_*`): GEN.1702 makes the request ask for land and pressure as well, for every request, tier or not.

## 8. The habitable family (GEN.1721)

Boss, 8:24 PM EST 2026-11-14: classes H (desert), O (ocean), P (glaciated) and R (hothouse) must be present beside M, with M the most common of the five, and none of them forced. The five classes share one air (oxygen 16 to 35.5 kPa, carbon dioxide under 1 kPa, pressure 50 to 250 kPa) and differ only by water cover and temperature (`physics/data/planet_class_fit.json`), so a planet that reaches the shared air with a different amount of water or sunlight is one of the other four.

What the laws give without help, on 800 Sun-like systems of main after Train 28 (`reports/train28-census-2026-11-14/five-classes-800.md`): planets with the shape of H 19, M 3, O 31, P 34, R 20, and classes H 0, M 1, O 1, P 0, R 0. Where they are lost:

- **O** has oxygen and is stopped by carbon dioxide (28 of 31 over 1 kPa). That is a bug of the carbon balance on a world with no land and is fixed for every tier (GEN.1720), not dialed.
- **H, P and R** hold no oxygen at all (0 of 73). The only life that makes oxygen needs a surface sea at the present day, which a desert, a glaciated world and a hothouse lack or hold too little of, and the oxygen law needs the production to be going at the present day.
- **H** also needs groundwater, which 9 of 19 lack.

The dials of the family, each on the tiers only and each off on Science:

| Dial | State | What it moves |
|---|---|---|
| `dry_part` | built (section 3) | Part of the dialed planets take the dry water class, so desert-shaped planets are made beside the temperate ones. |
| desert seas | built (section 3) | A desert world's water leaves a few small seas, and each is a halite-saturated sea of real depth that life can start in, not a seep. This also answers the groundwater lever Boss named: a world with a sea passes Class H's groundwater condition, and the dry planets that failed it had no surface water at all. |
| `quiet_share` | built (section 3) | Quiet volcanoes: a third of the oxygen sink on half of the candidates (Boss, 9:43 PM EST 2026-11-14). It helps every class with a living sea and does nothing for a world without one. |
| early warm phase | to build | A warm, wet early age of at least 0.8 Gyr in which oxygen-making life starts (the research of GEN.1722), before the planet settles to its present state. |
| hardier oxygen life | to build | How far past the sourced limits of real oxygen-making life the tiers let it live. The sourced limits themselves are GEN.1723, for every tier, Science included (Boss, 8:39 PM EST 2026-11-14: the minimal reasonable, educated guess). |

**Why desert worlds stay rare.** Boss chose "both, small seas first" (card, 9:49 PM EST 2026-11-14). The oxygen law counts the open sea alone and leaves out the land plants that make about half of Earth's production, so a world with under 20 percent sea holds 16 kPa of oxygen only when it is warm, bright and quiet inside. The land plants are GEN.1729, researched and then built for every galaxy. Measured on a prototype (Abundant, 400 Sun-like systems): when the small sea was left a brine, 31 of 44 desert worlds held no life at all, because a brine that strong (water activity under 0.6) is past what life is known to bear; that is why the tier makes the sea an ordinary one. With the sea an ordinary one (same sample, 2026-11-14) 41 of 53 desert worlds hold life, and still none is Class H: the four warm small-sea worlds hold 0.001 to 12 Pa of carbon dioxide, and oxygen-making life needs 2 Pa. The carbon dioxide a world holds goes as its outgassing to about the fifth power (GEN.1720), so only 5 of 30 Earth-shaped planets and 1 of 5 desert-sea planets land in the band of 2 Pa to 1 kPa that all five classes need. That is the next thing to decide, for every class at once.

The meter's family table (`scripts/class_m_meter.py report`) shows, for each of the five classes, the worlds per 1,000 systems and how many planets with the class's shape pass each of its air limits, so each dial is measured as it lands.
