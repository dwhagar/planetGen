**Headline:** The stars of a system are drawn as a tree, so systems of three or more stars now come out of the fill.

### Changed
- **A system draws its stars from the `star-count` table (GEN.631).** The one-companion draw and its per-class chances are gone: the primary's stars and the orbits between them are one tree (`draw_member_tree`), with the number of stars the table gives for its class (most O and B stars are not alone, and a G star has three or more stars in about a tenth of systems). The two stars of the tree's pair are the primary and secondary as before, and each further star is a member of the system, stored and shown with its place in the tree; members have no planets of their own yet. A brown dwarf, a neutron star or a black hole among a living star's companions is left out, since none can be built as a star of a pair. The widest orbit is still the place's, and `+binary_system`, `+wide_binary` and the binary prevalences ask for the same things. The system's seed gives different stars than before.
- **A stored tree of three or more stars is built.** A scattered star whose members are all ordinary stars is built with all of them, not the first alone, and its stored rows are spent.
- **The pair's kind follows its axis.** A pair is a `close` (circumbinary) pair under 0.25 AU and a `wide` pair otherwise, whatever the tree called the orbit; its node carries that kind.
- **The shares of binaries and wide pairs the Generate page starts from** are remeasured: 23.9% of systems have a pair and 81.0% of the pairs are wide.

- **The star count settings reach the draw (GEN.267).** `--stars N` forces the number of stars (a count the tree cannot give a star, such as a pair around a star at the hydrogen-burning limit, is reported as unmet), `--max-stars` caps the tree, and the `multi_star` prevalence scales the chance of five or more.

### Removed
- **The metallicity factor on close companions** (GEN.103) no longer reaches the fill, and its audit (TEST.196) is removed with the draw it checked; the tree does not take the factor yet.

generation-output: changed
