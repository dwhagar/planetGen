**Headline:** The limit on the nightly update's local offset is written down and pinned by a test.

### Changed

- TEST.200: the orbital-updates design note says how good the offset stored after a nightly rotation is (one unit in the last place of the absolute coordinate: 28 km at the Sun, 898 km of rounding at the 262 kpc reach, 3.4e5 to 1.1e7 times below the smallest stored move) and what to do if generating anywhere needs finer offsets. tests/test_local_offset_limit.py pins it at 8, 26 and 262 kpc over random centres, offsets, rotations up to a whole orbit and landing sectors; the worst error found is exactly one ulp. Tests and documentation only: the update and every stored value are unchanged.
