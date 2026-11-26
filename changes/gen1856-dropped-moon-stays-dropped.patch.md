**Headline:** A moon that the size rule took out of its planet no longer comes back bare with no spin period, which stopped a sector save with "rotation_period_hours cannot be null".

### Fixed

- `planets.drop_unstable_moons` takes a moon out of `planet.moons` when it is too big for its planet, but the planet's born-moon list kept it, and `refinish_rocky` runs the moon history again from that list. `moon_history._settle` re-checked the orbit rules and not the size rule, so such a moon came back bare (deferred climate set, no spin), after `finish_moons` had already run (GEN.1856). `_settle` now holds the size limits of `planets.moon_size_limits` too, and `StarSystem._settle_to_final_age` finishes the moons again after `refinish_rocky`, so a moon the history gives back through a new orbit comes back finished. A system with no size-dropped moon is as before; one with such a moon now keeps it dropped.
- `test_dropped_moon_stays_dropped` builds the case (a moon the size rule dropped is not back; no moon is bare or without a spin at the end of the settle pass).

generation-output: changed
generation-affects: fill
