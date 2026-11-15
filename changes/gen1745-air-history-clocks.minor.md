**Headline:** A world's air history has its first part: the step schedule and the clocks (the star's light, the interior's heat, the oxygen law's reductant sink and the water lost) at every age, with no change to any generated world.

### Added

- GEN.1745 slice 1: `physics/air_history.py` gives the ages of a history (`schedule`: logarithmic from 10 kyr to 1 Gyr, then every 0.25 Gyr, cut at the ages asked and ending at the star's age) and `run(body)`, the steps with the star's luminosity over its present one, the interior's heat (`interior_heat.heat_at`), the reductant ratio and label of the oxygen law (`biogenic_oxygen.reductant_history_factor`, `profile_of`) and the water lost by that age. Reuse: every law is the existing one read at the age of the step. Nothing in the generator reads it yet, and nothing is stored.
- `runaway_history.light_at(star, age_gyr)` (the star's light on the history `loss_kg` follows, at one age) and a `trace` for `runaway_history.loss_kg` (the water lost at the end of each step; the answer is unchanged).
