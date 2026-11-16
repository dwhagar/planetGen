**Headline:** The tide-wears-the-orbit test now uses a body whose mantle melts, since a small body that conducts keeps its mantle solid at the distance the test used.

### Changed
- GEN.1770: `test_a_strong_tide_wears_the_orbit_down_until_it_spends_no_more_than_the_orbit_gives_up` took the 0.0013 Earth-mass sample-15 body at 0.0134 AU, which only melted because a body that small carried no heat out of its mantle. With conduction it stays at 1,082 K and makes 1e11 W, so the test's premise (a molten mantle) is gone; it takes the same body at 0.0104 AU, where the mantle melts (the circularisation time is 13 Gyr, the settled eccentricity 0.024), and a new test pins the solid case at 0.0134 AU.
- Tests only; no generation output changes.
