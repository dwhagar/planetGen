**Headline:** A test now shows that a course between two systems of one sector also goes round an asteroid field.

### Changed
- **Asteroid-field avoidance is tested at the sector scale too.** `tests/test_nav_adjusted.py` places a field on the line between two systems of one sector and checks that the course bends round it (NAV.51; the between-sector case was already tested).
