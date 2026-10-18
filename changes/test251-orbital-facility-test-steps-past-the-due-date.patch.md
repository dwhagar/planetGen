**Headline:** The test that an orbital facility's phase advances no longer fails when the random giant it orbits is small enough that a quarter turn is not yet due.

### Fixed
- Tests: `test_orbital_facilities_advance_their_phase` stepped the orbit clock a quarter of the facility's period, but the update moves a facility only once it is due, which is when it has moved the planetary threshold (100,000 km). A facility on a close orbit round a small giant (about 56,000 km from the centre) moves 88,000 km in a quarter turn, so it is correctly not due and the update returned 0 (1 draw in 20). The test now steps whole turns past the due time plus the quarter (whole turns wrap away), so it still checks the 90 degrees whatever the giant. The update itself is right.

TODO items resolved: TEST.251 / Added: none
