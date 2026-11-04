**Headline:** A rogue planet can now be built through the planet code against a host with no light: the first part of making a rogue the same planet as a bound one.

### Changed
- **The starless host (GEN.1111, TEST.498).** `generation/phenomena/starless_host.py` holds `StarlessHost`, the donor star a body was ejected from with no luminosity and the rogue's own age, and `build_ejected_planet`, which makes a `Planet` of a given mass against it (the donor's mass and snow line, the formation distance of the loss, where `FAR_AU` puts it). Nothing reads it yet; the rogue rows from the architecture's ejections will, after the table rule below. Two strict expected failures in `test_starless_host.py` record why a whole body cannot be built yet: the atmosphere table's shoreline row (d1) is a logistic in the log of the insolation and of the host's luminosity, which are zero for a body with no light, so the weight is not a number. The rule for the table rows that read the star (a body with no host keeps its air as the ejection left it) is the next part and changes the table, so it goes in the lock wiring.

generation-output: unchanged

TODO items resolved: none (GEN.1111 stays open).
