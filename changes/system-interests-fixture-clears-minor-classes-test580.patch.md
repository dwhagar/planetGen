**Headline:** The life-class tests of the system interests set the classes of their fixture system completely, so a minor class the generator drew there no longer adds a reason the test did not ask for.

### Changed

- `tests/test_system_interests.py` (TEST.580): `_set_life_classes` set `life_class` on the fixture system's planets and cleared it on its moons, but the other classes that live on a generated world are rows of `planet_life` and `moon_life` (GEN.1209), and the interest reasons read them too. After train 12 the fixture system (the first of 120 seeds with two planets) holds a Class 3a row there, so 4 tests found `lifeclass-3a` that they had not set (`test_the_world_query_leaves_the_life_classes_out...`, `test_a_v175_galaxy_gets_the_life_class_reasons...`, and two more). The helper now deletes those rows too, so the reasons are the ones the test sets.

generation-output: unchanged
