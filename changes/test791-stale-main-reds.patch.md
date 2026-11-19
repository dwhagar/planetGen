**Headline:** The cache warm-up leaves the same listing windows in the cache whatever size the queue cuts its units to, and five stale tests agree with the product again.

### Fixed
- `web/warmup._objects`: a slice of a warm-up unit asked the listing for exactly its own rows (`limit` = the slice, `offset` = its start), so a stage cut into four units left four windows no visitor ever asks for and a stage cut into one unit left a different one. A slice now asks for the listing's own windows (`PAGE_SIZE` rows from a multiple of `PAGE_SIZE`), which are the windows a visitor's listing page asks for, so any cut leaves the same ones and `warm-sectors`, `warm-systems` and `warm-phenomena` are cut-neutral as their flag says (`test_cut_neutral`).
- `tests/cut_neutral_recipes.py`: the stand-in settings for the node warm-up lacked `page_cache`, which the warm-up reads for its cap.
- `test_gen_fill_order`: the run also stores the sectors of the below-cut remnants it draws outside its wedge (GEN.1682); the test now picks the sector at the address it asked for.
- `test_unit_reports`: a piece of the map luminosity fill weighs its ids (OPS.243), so the unit of ids 1 to 10 reports 0 to 10, not 0 to 1.
- `test_generation_stats`: the console says a warning in lower case.
- `README.md`: the external-files section is regenerated from the catalogue (the Gaia query now carries its `g.` alias).

generation-output: unchanged
