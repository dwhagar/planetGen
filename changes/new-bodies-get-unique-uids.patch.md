### Changed
- GEN.173, GEN.174: a planet, moon or belt added by an admin edit now gets a uid (it was saved with none), and after a delete the new body's uid skips ones its siblings already carry (it failed with IntegrityError 1062 on `uq_planets_uid`).
