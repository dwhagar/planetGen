**Headline:** The habitability browser tests have a seeded system that holds rated planets and moons again.

### Fixed
- TEST.754: `HABITABLE_SEED` in `tests/habitability_site.py` is 59 (it was 64). After the recent trains seed 64's system held three rated moons and no rated planet, so `test_the_lists_draw_each_bodys_stored_rating_as_the_symbol_of_its_tier[planet-planets]` failed ("the seeded system has no habitable planet"). Seed 59 holds five rated planets and nine rated moons; the 48 tests of the habitability, galaxy outlines, score filter and map icon browser files pass.
