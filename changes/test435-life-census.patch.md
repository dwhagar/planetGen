**Headline:** The prevalence suite now checks how far the life of generated worlds gets (microbes, plants, animals, sentient) against Boss's targets.

### Added
- **The life census (TEST.435).** `src/tests/audit/prevalence_life.py` generates 240 seeded systems, puts every terrestrial planet and moon that has life on the highest rung it reached, and checks the shares against the bands of `src/tests/data/prevalence_life_expectations.json` (the majority stop below animals, 15 to 35 percent reach animals, at most 2 percent are sentient, fewer than one in twenty animal worlds are sentient). The first run found 1,246 life worlds in 3,917, none past producers, so the animal band fails; it is tagged GEN.855 (no biogenic oxygen) and the test fails once the tag is removed and the band is still missed. It also counts the shadow life class (GEN.886) of each life world for the report: 794 of the 1,246 have no class at or above the line.

TODO items resolved: TEST.435 (the stored-class read waits on GEN.889). Added: none.
