**Headline:** The mixture subject of a body's air setting is a class and not a closure, so a setting can go to a process pool again

### Fixed

- TEST.805: `surface_air.setting_of` built `subject_of` as a local function that a process pool cannot pickle (`Can't get local object setting_of.<locals>.subject_of`), which failed `test_describe` (6 tests) and `test_insolation::test_a_lone_star_swings...` on main. It is now the module-level class `_BodySubject`, called the same way; no generated number changes.

generation-output: unchanged
generation-affects: none
