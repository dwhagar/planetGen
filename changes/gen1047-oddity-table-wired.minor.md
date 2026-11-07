**Headline:** Every planet, moon and rogue planet now goes through the oddity table while it is built, so the rare outcomes of the weird-worlds design have a place to happen; the table holds no wired row yet, so no body differs (only the stored form of a planet gains an empty `oddities` list, which moves the golden structure and rosetta files), and a body that gets one will say so on its page and in the API.

### Changed
- `physics/oddities.py` gains `apply(body, levers)`: it asks the table for the rows whose gate holds for the body (`quantities_of`), rolls each on its own stream keyed by the star, orbit and mass (`body_key`, no draw from the run's stream), runs the effect the row names (`EFFECTS`, filled by `register`), and records the ids on `body.oddities`. The planet build calls it for the bulk inputs (lever L3) right after the composition is drawn; later wiring items call it at their own steps.
- `SystemConfig.FORCE_ODDITIES` forces rows whose gate holds (a body still cannot be given what its conditions forbid) and `SystemConfig.WEIRDNESS_DIAL` sets a run's own dial; the tuning dial applies when a run has none, and 0 turns every oddity off.
- `planets`, `moons` and `rogue_planets` gain `oddities` (VARCHAR(40), comma separated ids, NULL for none): schema v226, migration 0226 (number provisional), stored and loaded with the body.
- `GET /api/objects/<ref>` carries `oddities` (`[{id, name, text}]`) and the object page shows an "Odd" row in At a glance for each, name and then the sentence that says why.
- Oddity W10 (GEN.1048, iron-rich core) is wired: the planet architecture's merger count reaches the planet (`collision_count`), and a rocky planet under 10 Earth masses with a merger can have lost much of its mantle to a giant impact (about 1 in 20): it keeps its mass, its iron share rises (`physics/mantle_stripping.py`, `strip_mantle`), its water goes with the mantle, and its radius and density follow. The energy ratio of the impact is a new stand-in (MT8, tweak L3, burn-down line added).
- Not in this item: the dial in the settings page (it needs its own settings group), which follows as its own item.

generation-output: changed
