**Headline:** The ice-species storage test now picks a system with two rocky planets, so its Pluto-like body is no longer overwritten by the CO2-frost one.

### Fixed
- `test_a_pluto_like_and_a_co2_frost_mars_like_body_keep_their_species_through_the_database` gave the first rocky planet a Pluto-like split and the last one a CO2-frost split, but took the first system of 60 seeds that had any rocky planet. When that system had one (the draw moved after the stored-composition and rogue work), both splits went to the same body, the second overwrote the first and the stored N2 amount read 0.0. Nothing was lost between the generated body and the stored row; the search now needs two rocky planets. No generator change.

generation-output: unchanged
