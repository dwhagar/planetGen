**Headline:** A test now fails when a new numeric constant is added to the generator's own modules without a registry row, a tag or a place in the constants module.

### Changed

- TEST.715 (the gate TEST.455 names, from the GEN.945 census): `scripts/constants_census.py` is Lane 21's census made importable; `test_constants_completeness.py` fails on a new unregistered constant in physics, galaxy, generation, population or `tuning.py`,. The committed list (`src/tests/baselines/unregistered_constants.txt`, 1,521 names today) may only shrink. Tests and a script only; nothing generated changes.

generation-output: unchanged
generation-affects: none
