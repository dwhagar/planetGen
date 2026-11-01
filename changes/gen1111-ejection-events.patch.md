**Headline:** The events a rogue planet is built from can now be drawn from a donor system's record: whether each thrown-out body left, when, how fast and in which direction, and how wide a moon system it kept.

### Changed
- **The ejection events (GEN.1111).** `generation/ejection_events.py` turns the bodies a system's planet architecture threw out (`StarSystem.ejected_losses`) into one event each: whether it ended ejected (by the share the architecture gives), the delay after formation, the age now and the time since, the kick of the scattering and the planet's orbital velocity at the ejection (both isotropic in direction, the orbital one at the circular speed of its distance), and the outer edge of the moon system it kept. Every number comes from the donor's own stream in a fixed order, so the set of rogues of a donor is the same whatever else is charted, and an event's index is the key a rogue row will carry. Nothing reads it yet; the donor's save does in the lock wiring.

generation-output: unchanged

TODO items resolved: none (GEN.1111 stays open).
