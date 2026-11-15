**Headline:** The life-air history API tests find a planet that can take the forced oxygen in the redrawn systems of Train 30.

### Fixed
- `tests/test_api_life_air.py`: the fixture asserted that `biogenic_oxygen.apply` took 2e19 mol of oxygen on the first airy rocky planet of each drawn system. After the star redraw of Train 30 the first such planet's air held the oxygen already (or could not take it), so two tests failed; a planet that cannot take it is now passed over for the next, and the fixture still skips the test only when none of 40 systems has one.
