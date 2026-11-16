**Headline:** The air history now gives the surface temperature of a world's air at every step from the star's light and the interior's heat of that age, with no change to any generated world.

### Added

- GEN.1745 slice 1b: `air_history.run(body, cuts, air_input)` with the input of the body's final air solve gives each step a surface temperature. `aged_setting` is the air solver's input at an age (the star's light and the interior's heat of the age); the steps where the light or the heat has moved by 3 or 25 percent are solved in full from the draws the final air began with (at most 12), and the steps between are read between their two anchors in the starlight absorbed. Reuse: `surface_air.solve` and `Setting`, `runaway_history.light_at`, `interior_heat.heat_at` and `released_by`, `Air.entry_state`. Nothing in the generator reads it yet.
