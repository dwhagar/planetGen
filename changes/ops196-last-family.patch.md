**Headline:** the clears, the drops and the reset now run as pieces the engine counts and can stop, so no command is left that works in one loop of its own.

### Changed

- OPS.196: `clear-bright-stars`, `empty-phenomenon-scatter`, `empty-nebula-field`, `drop-phenomenon-indexes`, `clear-scatter`, `clean-nebula-field` and `reset-galaxy` give the engine `pieces` and `do` instead of a `unit` function. A clear's statements are a plan (`store.clear_bright_stars_plan`, `clear_bright_star_pass_plan`, `clear_phenomenon_scatter_plan`, `clear_nebula_field_plan`), one piece to a statement; the reset is one piece to a table and one to tidy after; the same plans are what `clear_*` and `wipe` run directly, so there is one list of the statements. Their results are unchanged.
- `tests/test_commands_have_pieces.py` no longer has a list of commands still to convert: a command with a loop of its own fails. `tests/test_longest_silence.py` runs the commands of the last family under the unit beat and holds the longest stretch with no count, phase or statement moving to a budget.
- `tests/test_long_sql_steps.py`: the clear tests were red on main (a new galaxy has nothing to clear since OPS.254 and OPS.257); they now say the galaxy holds something, and the reset test seeds a sector.
