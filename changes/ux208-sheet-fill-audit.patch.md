**Headline:** The property sheet now states, for every column, the physical rule that makes it not apply to a body, and a test fails when an applicable column is empty without a filed item.

### Added
- **The fill audit rules (UX.208).** `describe.sheet.APPLIES` now says why a column does not apply: air columns only where the body has air pressure, the solar day and daylight only on a world that turns on its own, the shelter column and outside hours only where the radiation calls for a shelter, ice shell and high-pressure ice only on a body whose water is ice or frozen over, ring and aquifer columns only on a body that has them, the critical orbit only for a star with a companion, and the orbit period of a comet only on a closed orbit. `OPEN_GAPS` lists the columns that are empty where they apply, each excused by the filed item that closes it (GEN.768, 778 to 781, 784 to 789, 904, 905). `tests/test_sheet_fill.py` generates a sample and fails for any other empty applicable column. The report is `reports/sheet-fill-report-2026-10-29.md`.

TODO items resolved: UX.208. Added: GEN.904 (daylight columns stale after an orbit move), GEN.905 (class V mantle oxidation empty).
