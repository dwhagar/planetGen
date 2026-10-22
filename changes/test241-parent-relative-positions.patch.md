**Headline:** A test now holds the rule that nothing inside a system stores an absolute position.

### Changed
- `db/point_columns.py` names the columns that tie a row to its system or body (`BODY_PARENT_COLUMNS`) and the four tables that keep a free galactic point despite one (`ABSOLUTE_WITH_A_PARENT`: bright_stars, black_holes, neutron_stars, facilities), each with its reason and its way out in the 8.4 revision.
- The rule (a body in a system is stored parent-relative, about the parent's central object, in the galactic axes, with its epoch) is written into the header of `physics/position.py` and into docs/design/unified-system-barycentre.md.

### Added
- `tests/test_parent_relative_positions.py` lists every table with a parent link and a point, fails if one is stored in another frame and is not allowed, fails for a table given an absolute position, and refuses an exception that is no longer needed.

TODO items resolved: TEST.241. Added: none.
