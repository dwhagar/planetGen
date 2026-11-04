**Headline:** A rogue planet can now be built from an ejection event through the planet code, with the moons the ejection left it, and the sector it lands in can be found from its donor's place.

### Changed
- **The rogue body of an event (GEN.1111).** `starless_host.build_rogue_body` makes the planet an ejection event describes through the one planet chain: the event's mass and formation distance, no moons while it is made, then the Hill radius of its moon window set so that the outer edge is the one the ejection left (`EjectionEvent.moon_limit_km`), then its moons drawn inside that edge. A limit narrower than the planet's own clearance leaves no moons; a wider one lets them lie farther out.
- **Where the rogue is now (GEN.1111).** `generation/rogue_landing.py` takes the donor's place in the galaxy and the event's kick, orbital velocity and time since the ejection, and gives the rogue's position and velocity in the galaxy frame and the sector that holds it: the along-orbit offset is an angle round the galactic centre (`rogue_drift.ring_drift`), so an old ejection lands on the far side of the rotation rather than on a straight line, and a donor on the galaxy's axis uses the offset in the galaxy axes. The velocity is the circular speed along the rotation at the new place plus the motion against it.

Nothing reads either yet; the donor's save does in the lock wiring. The whole body of a starless rogue still waits on the atmosphere table's rule for a body with no light (two expected failures in `test_starless_host.py`); the moon bound is tested with a real star as host.

generation-output: unchanged

TODO items resolved: none (GEN.1111 stays open).
