# Class M abundance tiers

**Status:** being built (Boss approved the design by question cards, 11:17 to 11:19 PM EST 2026-11-13). Dial values are first-pass stand-ins until the census calibration that follows GEN.1694 to GEN.1698.

**Code:** `planetgen/generation/abundance_tier.py` (the tier and its dials), `tuning.ABUNDANCE_TIER_DIALS` (the values), the `class_m` row of `planetgen/galaxy/dials.py` (the control and its storage), `scripts/class_census.py` (the meter).

## 1. Why a tier exists

The generator made no Class M world and no world rated Blue in 2,750 plain systems and 455 Sun-like systems (`reports/class-m-blue-levers-2026-11-14.md`). Part of that was bugs, filed as GEN.1694 to GEN.1698. What is left is not a bug. With the bugs fixed, the generator's own laws give under one Class M world in ten thousand Sun-like systems, for three reasons that the science itself gives:

1. A rocky planet in the ecosphere forms at a median 0.85 snow lines, where the Bern population synthesis (`physics/delivered_water.py`) gives it thousands of times Earth's water. Sixteen of eighteen candidates are 100 percent ocean. Earth formed at 0.37 snow lines and is in the law's `low` class.
2. Thirty-four of seventy-five right-mass candidates keep a hydrogen envelope.
3. The band of sunlight that gives a Class M surface is narrow: an exact copy of Earth is Class M 23 percent of the time at 1.00 AU, 9 percent at 0.95 AU and 2 percent at 1.10 AU.

Nobody knows how common an Earth is. Earth may be one in a galaxy. So the rate is not a fact the generator can look up; it is a choice of the person who makes the galaxy, and the generator says plainly when that choice goes past the science.

## 2. The tiers

| Tier | Target, Class M worlds per Sun-like system | What it means |
|---|---|---|
| Science (the default) | none | Only the laws. Nothing is moved. |
| Sparse | 1 in 1,000 | |
| Scattered | 1 in 100 | |
| Plentiful | 1 in 30 | |
| Abundant | 1 in 10 | |

A Sun-like system is one whose primary is a main-sequence star of type F7 to K2 (the sample of the report and of the meter).

The target is what the dial values are calibrated to. The generator never checks a galaxy against it and never adds or removes a world to meet it.

## 3. The three dials

A tier turns three dials and nothing else. Each is a share from 0 to 1.

| Dial | What it moves | Where |
|---|---|---|
| `water_share` | The share of candidates whose water class is drawn as Earth's (`low`, 0.025 to 0.16 percent of the planet's mass) and not from the Bern law at the distance they formed at. The amount inside the class is still drawn. | `physics/bulk_draw.py`, `physics/delivered_water.py` |
| `bare_share` | The share of candidates that are late cores: they finished growing after the gas disk went and captured no hydrogen envelope (the existing `Planet(late_core=True)`, GEN.1323). | `generation/planet.py`, `physics/envelope_capture.py` |
| `request_share` | The share of systems that are given the request for a habitable world, which multiplies the odds of the wanted draw by `tuning.REQUEST_ODDS` and never forces it (`util/request_odds.py`, GEN.1130). | `generation/prevalence.py` |

A **candidate** is a rocky planet, not a moon, whose orbit lies in its star's ecosphere (`abundance_tier.candidate`).

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

A system gets it through its config. The sector that makes a system sets `SystemConfig.ABUNDANCE_TIER` from the galaxy's shape (`FillContext.apply_population`), as it sets the population, so a system regenerated later is made under the same tier. The tier is not saved with the system's config: it belongs to the galaxy. A system made outside any galaxy is Science unless the caller sets the field.

## 6. Dial values and calibration

| Tier | `water_share` | `bare_share` | `request_share` |
|---|---|---|---|
| Sparse | 0.5 | 0 | 0 |
| Scattered | 1 | 1 | 0 |
| Plentiful | 1 | 1 | 0.5 |
| Abundant | 1 | 1 | 1 |

These are first-pass values from the estimates of the report (every candidate given Earth's water gives about 1 Class M world per 500 Sun-like systems; 1 per 100 needs the envelopes gone as well; the two top tiers need the request). They are registered as stand-ins of tier T9, tweak level L4 (`docs/design/research/gap-filling-plan.md` section 6.288).

Calibration is by the meter: `scripts/class_census.py` draws a seeded sample of Sun-like systems under a tier and prints the funnel (candidates, right mass, ground, gravity, water cover, liquid water, Class M, Blue) and the rate of Class M and Blue worlds per 1,000 systems with its interval, beside the tier's target. The values are set again after GEN.1694 to GEN.1698 merge, because those fixes change the rate under every tier.

**A target may be out of reach.** The dials only move draws. If every candidate has Earth's water and no envelope and every system carries the request, and the meter still reads under 1 in 10, then Abundant gives what the dials can give and the page says the measured rate. The remedy is a decision for Boss (a wider dial), never a forced world.

## 7. What is not part of the tier

- The physics fixes of GEN.1694 to GEN.1698 and the weathering and oxygen work: they apply on every tier, Science included.
- The research on why habitable-zone planets are born water-rich and on hydrogen envelopes. If it shows the law itself is wrong, the law is fixed for Science and the dial values shrink.
- The habitable request's own conditions (`tuning.REQUEST_*`).
