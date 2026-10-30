**Headline:** The test that checks a disk moon's ice against its disk's ice line draws enough systems to have moons on both sides of the line.

### Fixed
- **The disk-moon ice test (TEST.500).** `test_a_disk_moon_gains_ice_past_its_disks_ice_line` failed on main since train 4: the planet architecture makes fewer giants, the only hosts of disk moons, so the 40 seeded systems the test drew held three disk moons, all more than nine ice lines out, and none near the line. The fixture now draws at least 40 systems and goes on until it holds 30 disk moons (149 systems on this seed, with 8 near the line and 16 far). The generator is not wrong: the disk moons of a giant that has gas giants are placed by the same law as before. No generator output changes.

generation-output: unchanged

TODO items resolved: TEST.500.
