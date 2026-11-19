**Headline:** Kill test for a run canceled while its job is killed (OPS.246).

### Added

- A kill test with real processes: a run canceled while its job process is killed, in either order, ends with the unit canceled, never handed out again, and no piece done twice. The table of killed engine transitions names it.
