**Headline:** A test of the least pressure a planet's air has no longer fails on a thin air stored a part in a billion under 1 Pa.

### Fixed
- **The thin-air floor in test_planets (TEST.485).** A thin regime caps the air at exactly 1 Pa (`surface_air.VACUUM_PRESSURE_PA`) and the solve then drops every gas whose partial is under a billionth of the total (`STORE_MIN_FRACTION`), so the stored sum came out as 0.99999999993 Pa and the check for at least 1.0 failed. About one planet in forty of a plain G2V sweep stores a pressure like that. The floor is now 1 Pa less that billionth share of each gas the registry knows, with a test that the floor allows the trim and no more. Nothing generated changes.

generation-output: unchanged

TODO items resolved: TEST.485.
