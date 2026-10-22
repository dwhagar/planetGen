**Headline:** The precision table's stability-ceiling test imports the Hill sphere from where it lives.

### Fixed
- `test_precision.py` imported `calculate_hill_sphere` from `generation/star.py`, which no longer has it since the orbital update; it takes it from `physics/orbits.py`, as the other tests do.
