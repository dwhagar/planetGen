### Fixed

- `generate.py system` no longer saves a system that misses a forced body
  (GEN.49). A forced habitable world or asteroid belt that a star can't
  host (around an O5V star, a habitable world failed about two systems in
  five) is now tried on up to five whole systems; if none has it, the run
  exits with an error and saves nothing. Before, it saved the system anyway
  and reported success.
- `-planets` with `+asteroid_belt` is now refused, and for a single system
  `-planets` is also refused alongside the same keys, `num_orbits` or
  `slots` in a `--system-file` (GEN.50). Before, the belt was placed
  anyway.
