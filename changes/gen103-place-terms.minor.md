**Headline:** A star's chance of a close companion now follows its own metallicity, and a wide pair's widest axis follows the density and motion of the place it is in, in place of fixed numbers.

### Changed

- **Close companions follow [Fe/H], for every star mass.** The odds of a primary's companion inside 10 AU are multiplied by `F([Fe/H]) / F(0)`, the straight line through Moe, Kratter and Badenes (2019): 53, 40, 24 and 10% at [Fe/H] −3.0, −1.0, −0.2 and +0.5 (a G dwarf's chance goes 0.43, 0.36, 0.25, 0.12; the odds, not the chance, so no star's share passes 1 and an O star at −0.7 rises 29%, not 70%). The wide share does not move with [Fe/H], and the two together never pass 1.
- **A wide pair's widest axis is derived.** It is the smaller of the passing-star limit (Bahcall, Hut and Tremaine 1985, `a = 0.002 (M1 + M2) σ / (G ρ₂ t)`) and the galactic tide `(M / 3ρ)^(1/3)`, from the cell's density, its populations' velocity dispersion and the pair's age. The fixed 10,000 AU cap (`WIDE_BINARY_SEPARATION_MAX_AU`) is gone: the quiet disc keeps pairs out to about 29,000 AU, a bulge cell far fewer, and a place whose limit is under the 50 AU lower edge keeps no wide pair (the primary is alone, or a close pair when a companion is forced). A system drawn with no galaxy placement takes the Sun's neighbourhood.
- `StarSystem._companion` replaces `_should_generate_binary`: one roll decides none, close or wide, and the mass ratio is drawn with it. The same seed gives different companions than before.
- New `generation/place_terms.py` (the factor, the limits, `PlaceTerms` from a sector profile), `SystemConfig.PLACE_TERMS` (not saved), and the constants beside the binary mass ratio in `tuning.py`.

TODO items resolved: GEN.103 (the close-pair rule 3.1 and the wide-pair limit 3.3; the massive-multiples rule 3.2 is left to the multi-star items) / Added: none
