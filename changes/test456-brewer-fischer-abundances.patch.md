**Headline:** A test now holds the generated star abundances to the Brewer and Fischer 2016 table of 1,615 dwarfs, and shows where they differ.

### Added
- **TEST.456.** `test_stellar_xfe_brewer_fischer.py` compares the generated [C/Fe], [O/Fe], [Mg/Fe], [Si/Fe] and [Ca/Fe] with the table by [Fe/H] bin (medians within 0.10 dex, scatter within a factor of three). The table is reduced to bands by `scripts/build_brewer_fischer_bands.py` (stars not kept). It finds two gaps, held as strict expected failures: calcium (and sulfur and argon, which follow it) has no scatter where the table has 0.04 to 0.10 dex, and the generated [Mg/Fe] sits 0.05 to 0.07 dex above the table from [Fe/H] -0.2 up.

### Changed
- **The wiring register** names the photolysis table against GEN.982 and the radiation belts against GEN.983.
