**Headline:** A sweep now walks every possibility table of the generator and checks that its odds add up, its gates are reachable, its draws follow its weights and its open rows are listed.

### Added
- **The possibility-table sweep (TEST.214).** `src/tests/audit/sweep_tables.py` walks all 33 registered tables (the materials table, the 16 atmosphere tables, the bulk, surface, moon, small-body, comet, star-count and satellite tables) with seeded subjects of every class each names. It checks that every row has a source and a state, that every class reaches a row, that the gated weights of every subject add up to 100 percent within 1e-12, that no subject fails but with an empty gate, that each row's drawn fractions add up to 1, that `Table.draw` follows the weights, that every row is reached for every class it names, that the open rows are the ones listed with a note in `src/tests/data/possibility_sweep_notes.json`, and that no table is left unregistered. It found four rows no subject reaches (atmosphere table e, rows e1 and e2 for class I and e3 and e4 for class T), each recorded in the notes file with the reason.

TODO items resolved: TEST.214. Added: none.
