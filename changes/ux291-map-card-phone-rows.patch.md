**Headline:** On a phone the Galaxy Map's info card for a star system fits in half the map again: each row's label now stands over its value instead of beside it in a column too narrow for it.

### Changed

- UX.291: at phone width (390 by 844) the star-system card was 338 px high against the half-window limit of 329 px, because the label column of its three rows ("Habitable worlds", "Uninhabitable worlds") broke into pieces of a few letters beside a value column only about 100 px wide. Below 40 rem the card's rows are now one column, the label over its value, so every row takes one or two lines and nothing is hidden.

generation-output: unchanged
