**Headline:** The Magellanic Clouds now hold young and massive stars, from a star formation history worked out for each satellite instead of a typed share.

### Added
- **A satellite's star formation history (GEN.526).** Each satellite galaxy's young and intermediate stars come from a history set by its stellar mass, size and mean stellar age: the rate rises or falls over 13 Gyr to give that mean age, the gas the present rate needs follows from the Schmidt-Kennicutt law, and a satellite whose gas is under the star formation threshold forms no stars now. The Large Cloud comes out at 1.4% young and 36% intermediate stars, with 0.36 solar masses a year and 4.4e8 solar masses of gas; the Small Cloud at 1.0% and 26%; Sagittarius, Fornax and the dwarf spheroidals form no stars now.
- Two stellar populations, `satellite_young` and `satellite_intermediate`, beside `satellite` (which is now the stars of 3 to 13 Gyr). The sector profile lists all three, and `bright_stars.population` takes them (migration 0107, which also widens the column to 24 characters).

### Changed
- The Clouds place massive stars: at the default floor of 9,000 solar luminosities the Clouds hold about 3.8e4 bright stars (none before) and 1.8e4 stars above the mass cut. A galaxy scattered before this needs `planetgen plan --redo-scatters luminosity mass phenomena`.
- A satellite place counts its young stars as young in the Type Ia share, the planetary nebula share and the molecular cloud share.

generation-output: changed

TODO items resolved: GEN.526. Added: none.
