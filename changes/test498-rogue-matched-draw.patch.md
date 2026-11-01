**Headline:** The rogue-planet parity test now draws a rogue and a bound planet of the same mass and compares them: same kind of body and size, and a strict expected failure for the one thing that differs today, the surface temperature far from light.

### Changed
- **A matched draw for the rogue parity contract (TEST.498).** `tests/test_rogue_planet_parity.py` builds a rogue planet of each mass bin and a bound planet of the same mass 1e5 AU from a star (until a starless host exists, GEN.1111) and checks that both are one kind of body, with a rocky radius inside the range its density gives. A strict expected failure records the difference found: far from light the bound chain returns the 2.725 K ambient floor, the rogue chain the temperature its own heat gives (270 K for a 5 Earth-mass body), because the bound chain's equilibrium temperature carries no internal heat. It comes off when the starless host is built through the planet code.

generation-output: unchanged

TODO items resolved: none (TEST.498 stays open until GEN.1111).
