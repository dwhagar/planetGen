**Headline:** White dwarfs now cool along real model sequences, so a white dwarf's temperature, brightness and size depend on its age and hydrogen layer.

### Changed
- **A white dwarf's temperature, luminosity and radius come from the Bedard et al. 2020 cooling sequences (GEN.742).** A young white dwarf is hotter than before (up to about 110,000 K), a 10 billion year old one of 0.6 solar masses is about 3,000 to 4,000 K, and the radius no longer depends on the mass alone (a thin hydrogen layer cools faster and the core crystallises late). A white dwarf draws one more number, for a hydrogen or a helium atmosphere (80 percent hydrogen, as before). The sequences are the ones setup already downloads. A white dwarf above 1.30 solar masses and one older than the last model are held at the sequence's edge and say so in the star's state. A run that disregards missing data keeps the old power law as its labelled fallback. Design note: `docs/design/systems/white-dwarf-state.md`.
- **A disc's solids follow the star's metal mass fraction over the Sun-like reference (GEN.464).** They followed 10 to the power [Fe/H], which is 1.7 to 3 times too low for a metal-poor star rich in alpha elements, so such stars' planets now get more solids. Metal-rich stars barely move.

generation-output: changed
