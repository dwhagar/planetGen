**Headline:** A test now shows a course across the dense galactic core goes round what is made and counts the unresolved stars it crosses.

### Changed
- **Dense-core course test (NAV.76).** `tests/test_nav_dense_core.py` stores a galaxy with a made central hole and a made asteroid field on a course across the core, and checks that the course bends round both, that it carries the count of unresolved stars it crosses, that the count matches `UnresolvedField` along the course's polyline, and that the bent course crosses fewer than the straight line (it skirts the hole's sphere). A course far past the core carries no count.
