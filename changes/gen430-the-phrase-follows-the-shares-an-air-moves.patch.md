**Headline:** A rocky body's composition phrase still names the level of each material after its air takes a sliver of the solids, so a share on a level's edge is not left one level too high.

### Changed
- **`envelope_share.fill` writes the composition phrase again after it scales the solids (GEN.430).** The phrase states the stored shares (`tests/audit/sweep_bulk.py` rule `text-materials`), but the air's share comes out of the solids after the phrase was made, so a body whose water share sat exactly on a level's edge kept "a large share of water" while its stored share was 0.19999999682. A promoted belt member (GEN.1510) takes the belt's round shares (0.8 rock, 0.2 water), so one of the 24 sweep systems (`system24/Eghintis ABC VI`) failed. Only a body that already holds a phrase gets it rewritten; the phrase changes only where the scaling crosses a level's edge.

generation-output: changed
generation-affects: fill
