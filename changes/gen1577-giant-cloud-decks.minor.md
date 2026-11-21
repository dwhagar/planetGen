**Headline:** A giant planet now keeps where its hot-condensate clouds sit (chromium, manganese sulfide, sodium sulfide, zinc sulfide, potassium chloride) and the most each trace gas of its star's gas can hold, and its page shows them.

### Added

- GEN.1577: `physics/hot_cloud_deck.py` places each cloud on the giant's temperature-pressure profile (the grey profile of `hot_giant_profile`, 1e-6 to 1e3 bar, floored by a dry adiabat down from the 1 bar temperature its air is solved at) and gives the ceiling of each trace gas of `hot_condensates`, packed as JSON into the new `planets.cloud_decks` column (migration 0269, NULL for a body that is not a giant of a star and for a row made earlier). `hot_cloud_deck.generate` is called after the giant's air (`planets._calculate_atmospheric_conditions`); no draw is made, so no other value moves.
- The object sheet of a giant gains a row for the base of each cloud and one for each trace gas ceiling (`describe/deck_rows.py`).
- Stand-ins DA1 (the adiabat exponent, 2/7) and DA2 (the 1 bar anchor of the floor) of section 6.317; the wiring register drops the lines for `hot_giant_profile`, `hot_condensates`, `hot_cloud_deck` and the rows HGP1, HGP2 and HK1 to HK4, which are read now.

### Changed

- `hot_cloud_deck` drops the one-column-per-species plan (`column_names`, `stored_quantities`) it was built with; the packed column replaces it.
- `bulk_draw.composition_of(star)` names the element mixture of a star's gas that `bulk_draw` already built for its ratios.

generation-output: changed
generation-affects: fill
