**Headline:** The habitability browser tests have a new seeded system again, after the star redraw of Train 30.

### Fixed
- `tests/habitability_site.py`: `HABITABLE_SEED` is 65 (it was 64). After GEN.1623 and the star table fix redrew every star, seed 64's system held eleven gas dwarfs and two rated bodies, so all eight tests of `test_habitability_end_to_end_browser.py` stopped at setup; seed 65 holds fourteen rated planets and moons and the browser tests of that fixture (habitability, galaxy outlines, score filters, map icons) pass.
