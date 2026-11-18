**Headline:** A short N-body integration checks the generator's stability screens, with a report and a test; generation is unchanged.

### Added
- `planetgen/physics/nbody_check.py`: a leapfrog integration of a star's planets that says whether they collide, are ejected or have an encounter (GEN.1468).
- `scripts/nbody_spotcheck.py` runs it over generated systems and writes the report; `tests/test_nbody_spotcheck.py` fails when a system the spacing screen passes is broken by the integration.

generation-output: unchanged
