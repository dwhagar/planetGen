**Headline:** The precision table and the command touches know the new birth-number tables of GEN.1753.

### Changed

- GEN.1753 follow-up: `body_history.formation_distance_au`, `body_history.birth_distance_au` and `lost_moons.birth_distance_au` are numbers the generator drew and kept to run a body's history again, so `db/precision.py` lists them as `size` (no reach, no floor) and the tables that hold them are not tables with a position. `lost_moons` joins the tables `store.load_star_system` reads (`queue/touches.py`, `SYSTEM_RELOAD_READS`).

generation-output: unchanged
generation-affects: none
