**Headline:** The tests of the object lists no longer fail now and then for want of a moon in their fixture galaxy.

### Fixed
- TEST.640: the galaxy that `test_object_list.py`, `test_object_list_values.py`, `test_systems_values.py` and `test_object_index.py` share was drawn afresh each run, and some draws gave no moons, so the tests that read moons failed at random. The three systems are now drawn from fixed streams, and the first stream that gives at least one moon is used, so the galaxy is the same in every run.
