**Headline:** A star's page and sheet now show its chemistry: thin or thick disc, alpha enhancement, the metal mass fraction it formed with, and the helium a red giant's first dredge-up added.

### Added
- **Four derived rows in the Composition group of a star's sheet (UX.254).** `Disc population`, `Alpha enhancement [alpha/Fe]` (dex), `Metal mass fraction` and, for a red giant, `Helium added by the first dredge-up`, each with its source line. They are worked out when the sheet is made from the stored class, state, mass, [Fe/H] and age by `stellar_composition.star_makeup`; no column, no stored value. The draws are keyed as the planets' host ratios are (mass and the sector's distance from the centre), so the page and the planets agree. A star with no [Fe/H], or one whose downloaded data are missing, shows none.
- **`stellar_abundances.alpha_enhancement`**: the mean [X/Fe] of O, Mg, Si and Ca.
- **`src/tests/test_property_sheet.py`**: a sample of stars across both sequences and a red giant check the rows against `star_makeup` with their units and sources; a star with no [Fe/H] has none. `docs/design/composition/stellar-composition.md` section 14.2; the properties help page lists the four rows.

No generator output changes (the rows are derived at view time).

generation-output: unchanged

TODO items resolved: UX.254.
