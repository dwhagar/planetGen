**Headline:** A body's roundness, strength and final-age ground no longer depend on its class letter.

### Changed
- **The strength pass and the final-age redraw no longer read the class letter (GEN.1398, class-last audit F1 and F2).** `body_strength._is_giant` now tests the body type and the hydrogen-helium envelope share, as the classifier does, the second strength pass in `apply_class` is gone, and the final-age water and ground redraw is seeded from the body's mass, radius, distance, age and moon flag only. This changes the generator output, so it takes the generator epoch of the train that carries it (the Train Manager assigns the number).
