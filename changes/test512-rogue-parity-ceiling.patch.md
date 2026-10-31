**Headline:** The rogue planet parity test is green again: the gap to a bound planet grew by twelve columns the planet gained, and the ceiling now says which.

### Changed
- TEST.512: `tests/test_rogue_planet_parity.py` was red on main (160 planet columns missing from a rogue, ceiling 148). No change to the rogue path opened the gap: GEN.474, GEN.768, GEN.1034, GEN.1035 and GEN.1097 each added planet columns (12 in all) while the rogue is still a generator of its own, and the rogue's own columns are what they were. `MISSING_CEILING` is 160 and `PLANET_COLUMNS_ADDED` names the twelve by item; a new test checks the ceiling rises by exactly those names and that each is still missing from the rogue, so the number cannot drift again without a name. The gap still falls to zero when the rogue is built through the planet code (GEN.1111 and the rogue train). Test-only.

generation-output: unchanged
