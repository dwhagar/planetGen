**Headline:** Every atmosphere recipe row that waited on research now says whether it is partly verified or a labelled stand-in, and which builder replaces it.

### Changed
- The 43 `open` rows of the atmosphere recipe tables (GEN.480 to GEN.494) are now 15 partly verified rows with a stated limit and 28 stand-ins that name the builder item (GEN.509 to GEN.525) that replaces them; no row is open. The recipes are not yet read by the generator, so no draw changes.
- The NH3 activity factor stays the flagged bracket 1.05 to 1.10 (GEN.501): the Huebner and Mukherjee 2015 rate table could not be reached.
- The disposition of every row and the reason are in `docs/design/weight-research-dispositions.md`.
