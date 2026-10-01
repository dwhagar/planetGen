### Changed

- **Bright stars down to 100 solar luminosities by default.** `generate.py
  plan` now pre-places every star of 100 solar luminosities or more
  (about 220 million in a Milky Way, roughly an hour and a quarter plus
  the database load) instead of 500. `--bright-star-min-luminosity 500`
  still gives a quick test galaxy. The hottest white dwarfs, which sit
  right at 100, are never pre-placed and now always stay in a sector's
  own draw.
