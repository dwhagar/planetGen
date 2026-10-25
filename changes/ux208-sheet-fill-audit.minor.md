**Headline:** A report now shows which property-sheet columns are never filled for which kind of world.

### Changed
- `scripts/sheet_fill_report.py` generates sample systems and counts, for every kind, class and sheet column, how often it is filled, empty where it applies and empty where it does not. The rules that say where a column does not apply (giants have no surface, airless worlds have no air column, oceans need a liquid ocean, and so on) are in `describe/sheet.py` (`APPLIES`) with a sentence of reason each (UX.208).
