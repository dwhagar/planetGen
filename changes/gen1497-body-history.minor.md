**Headline:** A new table, `body_history`, keeps for every planet and moon the numbers its history leaves that cannot be worked out again from the stored values: the eccentricity and the spin (rate and tilt) it was born with and the two draws that decided its spin resonance captures. It changes no result.

### Added

- `body_history` table (migration 0241): one row a planet or moon that has any, filled when the system is saved from what the generator already keeps on the body (`_orbit_birth`, `_birth_spin`), with a `formula_version` so a later change to the formulas can tell the rows it can refill from the rows it must not. Nothing is drawn.
- An edit of a system keeps the rows of the bodies it keeps and deletes those of a planet or moon it deleted; deleting the system deletes its rows.
- `physics/body_history.py` names the columns and reads them off a body. The Kozai cycle numbers and the tidal lock age are not in it yet: each lane that builds them adds its own column.

generation-output: unchanged
generation-affects: none
