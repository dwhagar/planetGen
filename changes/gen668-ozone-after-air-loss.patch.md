**Headline:** When a world loses part of its air, its ozone and its ultraviolet index can now follow it.

### Changed
- **GEN.668: a loss of air takes the ozone with it.** `radiation.apply_air_loss(body, retained_fraction)` scales a planet's or moon's stored ozone column by the share of the air that is kept and works the ultraviolet index out again from it (a world that loses half its air keeps half its ozone and has a higher index). No step of generation removes air after the ozone is drawn (escape and the water loss act on the air before it), so nothing is called yet and no generated value changes; the function is for the first loss event that does.
