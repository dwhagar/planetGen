**Headline:** A redraw of a system keeps the effect ledger of its new bodies.

### Fixed
- **A redraw lost the `effect_ledger` rows of the new bodies and kept those of the old (DB.135).** `replace_star_system_content` moved the body rows to the kept system and left the ledger rows under the stand-in system, which deleting it removed, the way `body_history` and `lost_moons` were lost before DB.75b. The ledger rows now move with the bodies (the same `_SYSTEM_HISTORY_TABLES` step), and the old bodies' rows are deleted.

generation-output: unchanged
generation-affects: none
