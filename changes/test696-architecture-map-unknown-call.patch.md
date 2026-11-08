**Headline:** The architecture map test no longer depends on which methods the package happens to name.

### Fixed
- **Architecture map call test (TEST.696).** The test for "a method on an unknown value stays text" used `.enqueue`, which stopped being unknown once `queue/requires.py` gained the only function of that name, so the map resolved it as a probable call. The test now uses `connection.ping`, which no function in the package carries.

generation-output: unchanged

TODO items resolved: TEST.696 (architecture map).
