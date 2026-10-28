**Headline:** Changing a planet's class in the admin editor works again for a planet whose spacing is checked against an asteroid belt.

### Changed
- **Class change spacing check (TEST.423).** The admin editor's spacing stand-in (`admin/edits._StandIn`) had no eccentricity, which the belt-separation rule has read since migration 0163, so the class check of a planet next to an asteroid belt raised an error instead of answering. The stand-in is a circular orbit now (eccentricity 0). `test_admin_edits::test_recommended_class_change_saves_and_validates` passes again.

TODO items resolved: TEST.423. Added: none.
