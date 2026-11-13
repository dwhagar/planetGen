**Headline:** The air a world ends with is made in one step, with no change to any generated world.

### Changed
- **The steam-branch start temperature, the air oddities and the airless state are the solver's own step.** `planets.calculate_atmospheric_conditions` solved the air and then corrected the solver's input and output in three places of its own; `physics/air_final.py` now does all three (`make`) and returns one `FinalAir` that the planet stores without branching on how it came out (GEN.1657). The same draws come in the same order: 285 bodies of 30 seeds and 200 trials of the Earth analog are identical to main, so the golden does not move.
