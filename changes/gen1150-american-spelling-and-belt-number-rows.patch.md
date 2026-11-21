**Headline:** The trojan region files use American spelling, and the belt number test feeds the rows the query returns.

### Fixed
- GEN.1150: `centre` in the trojan overlay, the 3D system view, the system map, the style sheet and the schema reference is `center` (`scripts/american_spelling_sweep.py apply`; identifiers and test names follow).
- GEN.1150: the belt number test built its rows without the `lagrange_point` column that `belt_names.numbers_for_systems` selects; its rows now carry it and a trojan region among them, which takes no number and leaves no gap.

generation-output: unchanged
