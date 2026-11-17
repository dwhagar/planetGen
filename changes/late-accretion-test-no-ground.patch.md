**Headline:** the two late-accretion air tests read the planets that have a ground again, so they pass on main.

### Fixed

- Since GEN.1710 (the ground cut reads `envelope_escape.BARE_BELOW`), an Earth-mass planet of a G star at 1 AU that keeps an envelope has no ground: its air is solved at the 1 bar level like a giant's and the late accretion is not drawn for it. `tests/test_late_accretion.py` took the last draw of every such planet and raised `IndexError` for seeds 1, 3 and 5. `_earth_like` now gives `None` for the percent of a planet with no draw, the replay test skips those and reads seeds 1 to 10 (still at least 3 checked), and the temperate-air search takes only planets with a draw. What the tests assert is unchanged.
