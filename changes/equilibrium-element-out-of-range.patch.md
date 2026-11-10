**Headline:** A very hot giant no longer fails to generate when the gas table has no phosphorus gas above 6000 K; the gas is made of the other elements.

### Changed

- `gas_equilibrium.solve` raised "no species holds P at 7555 K" for an air with phosphorus (a metal of every giant's budget) at a chemistry temperature above 6000 K, where the table's two phosphorus gases (P2, PH3) end. A giant that hot (the example system of zeta Ophiuchi, an O9 star, draws one most runs) failed to generate, so `test_examples.py` was red on main for that file. The solve now leaves out an element the table has gases of but none in range at the temperature and solves the rest, for the default gases and for a named list alike; it still raises for an element the table has no gas of at all, and for a gas left with nothing to be made of (`NoSpeciesError`, a kind of `EquilibriumError`). Draws that used to fail now generate; nothing that generated before changes.
- Stand-in: the phosphorus left out above 6000 K is a gap in the gas table (atomic P and its ions), to be filled by research; until then it is simply not in the hottest gas.

generation-affects: none
