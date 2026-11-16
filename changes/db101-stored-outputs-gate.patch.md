**Headline:** A new number on a generated planet, moon or star that has no database column now fails a test, so a physics output cannot go unstored unnoticed.

### Added

- DB.101: `tests/test_stored_outputs.py` generates systems and fails when a generated planet, moon or star carries a public attribute that no field of its rows loads (`store.PLANET_FIELDS`, `MOON_FIELDS`, `STAR_FIELDS`) and that `HOMES` or `STAR_HOMES` does not place (another table, a derived value, or a reason it is not stored). A stale entry fails too. The audit of the 2026-10-24 wiring table is in `reports/db101-storage-audit-2026-11-16.md`: every built item stores its output; the unstored public attributes found are listed with their reasons.

generation-output: unchanged
generation-affects: none
