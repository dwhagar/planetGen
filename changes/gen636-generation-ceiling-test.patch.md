**Headline:** The test that ties a system's orbit ceiling to the validator's now says what the two rules are: the validator applies the stable orbit, and generation also stops at the place's survival limit.

### Changed

- `test_system_methods_still_delegate` failed on main since GEN.325 cut `StarSystem._orbit_ceiling_au` to the survival limit. It is replaced by a test that the generation ceiling is the stable orbit cut by the survival limit, and one that a body past the survival limit but inside the stable orbit still validates (the survival limit is a chance of outliving the passing systems, not a stability limit).
