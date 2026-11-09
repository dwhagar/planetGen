**Headline:** The effect ledger no longer stores a row for a body that has only one star sending it light, and an edit of a system clears the system's ledger rows. A default system now keeps about 2.4 light and ultraviolet rows on average, not 4 for every body.

### Changed

- A light or UV effect on a body with a single source writes no row: that source is the body's host star and the flux is on the body, so the rows begin where a second source matters (a binary or a multiple). The gravity effects keep the largest source and the summary row as before.
- `edits.save_system_edits` clears the system's ledger rows (an edited body's shares are no longer the ones generation noted; they stay empty until the system's ledger is refilled), and `edits.forget_ledger_rows` takes the rows of a standalone phenomenon that is deleted alone or with its sector, as receiver and as source.
- The module note of `physics/effect_ledger.py` says that light and UV sources are the stars of the body's own system, not the sky's.

generation-output: unchanged
generation-affects: none
