### Removed
- **Ten tests (2,819 parametrized cases) that no longer checked anything.**
  Assert-free "generates
  without error" tests whose very next test runs the same generation and
  asserts on it (planets, the full star-by-class matrix, every star type,
  the example files), three comet-orbit validation tests repeated word for
  word in `test_kepler_motion.py`, three diminutive-prefix tests covered by
  the fuzz walk over every prefix, a check that a constant equals its own
  definition, a `None` Markdown check already in `test_mdconvert.py`, and a
  one-off TODO-marker check `test_todo_tags.py` now covers. The guard that
  `spaceSector.py` never imports a root script checked for `systemGen`,
  which no longer exists; it now checks `generate` too.
