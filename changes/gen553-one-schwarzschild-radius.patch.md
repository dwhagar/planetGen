**Headline:** The Schwarzschild radius and a star's birth age are each drawn or computed in one place (GEN.553 slices A and B1, no behaviour change beyond rounding in the black hole radius).

### Changed

- Schwarzschild radius (slice A): `multistar._schwarzschild_au` and `dark_light.disk_radius_km` called their own `2GM/c^2`; both now call `relations.schwarzschild_radius` (reused, no new function). The values are the same to rounding: about a third of the stored radii differ in the last binary digit (the multiplication order changes), nothing a draw or a limit can read.
- Birth age (slice B1): the draw of a birth age inside a window (uniform, or by the population's age law) was written out in `sample_star_age_gy`, in `sample_living_star`, in `remnant_progenitors.collapse_age_gy`; and the pick weighted by the age-metallicity relation in `sample_living_star` and twice in `star_population`. Now `stellar_evolution.draw_age_gy` and `stellar_evolution.birth_age_gy` hold them and every site calls them. The random numbers are drawn in the same order as before, so every star, remnant and system comes out the same.
- GEN.553's later slices (birth step on its own stream, the free brown dwarf leaving `rogue_planets` as GEN.1849) are not in this change.
