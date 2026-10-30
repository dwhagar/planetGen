**Headline:** The air of a tidally locked planet can now be tested for freezing onto its night side (not yet used by the generator).

### Changed

- GEN.1044: `physics/nightside_collapse.py` returns, for a locked planet's main gas, the pressure under which it freezes onto the dark side and the share that stays: CO2 by the fit of Wordsworth 2015 (0.14 bar at Earth's flux for an Earth-mass planet), other gases by the thin-air night-side balance. Registry rows NC1 to NC6 (section 6.75). Nothing the generator draws changes; the wiring is GEN.1051.
