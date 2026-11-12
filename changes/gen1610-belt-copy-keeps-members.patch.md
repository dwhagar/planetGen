**Headline:** A belt that is copied or sent to another process before its members were read now carries them.

### Changed
- GEN.1610 (bug): a belt draws its members on first use around the star it was built for, and that star is held apart from the belt in a weak map (`belt._STARS`) that a copy or a pickle does not carry, so a copied belt whose members had never been read returned none. `AsteroidBelt.__getstate__` draws the members first. No stored value changes: a belt that is stored has had its members drawn. `tests/test_belt_copy.py` pins `copy.copy`, `copy.deepcopy` and a pickle round trip (red before, green after).

generation-output: unchanged
generation-affects: none
