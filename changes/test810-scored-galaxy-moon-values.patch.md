**Headline:** The scored search sample keeps adding systems until enough moons hold every searched value, instead of stopping at six systems that may fall short.

### Fixed

- `tests/scored_site.py` stopped at `MAX_SYSTEMS` (6) systems even when some moon column (for example the period or ocean depth) was held by fewer than four moons, so `test_search_values.py` errored 27 times after the generator moved. The fixture now keeps scanning seeds 11 to 80 until four systems exist and every searched value is on four moons, and its failure message gives the real count per column. The requirement is unchanged (TEST.810).

generation-output: unchanged
