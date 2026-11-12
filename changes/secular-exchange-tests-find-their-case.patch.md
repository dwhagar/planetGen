**Headline:** The two belt-pull tests of the secular exchange find a system with a massive belt by search, so a change of the draw order no longer breaks them.

### Changed

- GEN.1667: `test_secular_exchange_generation.py` no longer relies on seed 1 and seed 2 making a belt with a mass; the belt-pull test searches seeds 1 to 60 for a massive belt whose rings move a planet's orbit, and the rings test takes the first system whose belt gives rings (as TEST.714 did for the others). Belts were not lost: the draw order of the train moved which seeds make one.

generation-output: unchanged
