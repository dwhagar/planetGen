**Headline:** Three tests that failed on a clean main are fixed on the test side: a moved seed, a mixed composition and a stale plan key list.

### Changed
- **TEST.481: the Hill sphere moon test, the halo flux test and the component pass plan test pass on main.** (1) `test_shrinking_the_hill_sphere_drops_the_outer_moons_first` needs a planet with two or more moons, and how many of the six asked for keep a stable orbit is a draw (1 to 6 over seeds 3301 to 3320); the file's seed gave one after the star epoch moved the draws, so the test draws that planet from its own seed. (2) `test_the_halo_counts_match_the_evolutionary_flux_of_an_old_population` took the Sun's composition for the halo's mean luminosity and turn-off track while the count it checks is at the halo's reference [Fe/H] of -1.6; at -1.6, with the turn-off mass of the same 12.5 Gy lifetime (0.82 solar masses), the model's count is 1.06 to 1.09 times the flux's at 100, 500 and 1000 Lsun, inside the 35% the test allows. The mean of too few draws also moved it by a third (a 40,000 draw mean has a spread of 25%). (3) `test_a_large_components_phenomena_are_a_unit_per_slice_that_add_up_at_any_order` compared the whole plan with the shape's JSON, and the plan has said `outline_shells` and `min_mass_solar` since PERF.212. Nothing the generator makes changes.

generation-output: unchanged

TODO items resolved: TEST.481.
