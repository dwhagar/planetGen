**Headline:** A gas giant's carbon-to-oxygen ratio now comes from the disc it formed in (its snow lines and its star's carbon and oxygen), spread as the observed giants are, instead of a draw around 0.57 (GEN.1363 part C).

### Changed

- `physics/disk_condensation.envelope_c_over_o`: the C/O of an envelope made of the disc's gas plus the solids that froze where the giant formed, as many solids as its enrichment over the star asks for. A lightly enriched giant takes the gas C/O of its place; a heavily enriched one the oxygen-rich solids' (about 0.34 for the Sun's disc at 60 to 140 K).
- `mixture_subject` gives a giant that is not a moon its formation temperature (182 K at the snow line, falling as one over the square root of the distance) and its star's carbon and oxygen over the Sun's; `mixture_giant` multiplies the disc value by the recipe's drawn `f7` over 0.57 (its core's center), so the observed spread stays (stand-in DSC7, section 6.213). A giant with no disc (a brown dwarf, a rogue) keeps the recipe's value.
- The median C/O of the sampled G2V giants falls from 0.57 to about 0.34 (disc value 0.28 to 0.40 before the spread). The cause is that enrichments of 3 to 80 times the star's draw mostly solids, and the solids are oxygen-rich; the published median of 0.57 is not reproduced by the disc alone, which is reported to Boss in the part B report.

generation-output: changed
