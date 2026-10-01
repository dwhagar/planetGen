### Changed

- **Bright stars can go down to 100 solar luminosities.** The default
  stays at 500 (about 60 million stars in a Milky Way, about 10 GB);
  `generate.py plan --bright-star-min-luminosity 100` now works too
  (about 220 million stars, about 35 GB). White dwarfs, which are never
  pre-placed, always stay in a sector's own draw.
