**Headline:** The air solve now has tests for what it promises about its start: the same answer for the same inputs, and one temperature for a giant whatever the start.

### Changed

- GEN.1555 (tests, built for Lane 14): `test_surface_air_start_independence.py` is rewritten after TODO's decision that real bistability stays: the strict-xfail "one fixed point from every start" test is replaced by (a) the same setting, start and draw state give the same answer twice, and (c) a giant (seeds 3009, 3010, 3000) ends at one temperature from every start (fails on the earlier one-step solver, passes with the giant fixed point); the seed 3032 test stays. The census script probes the pinned bodies under the production draw state as well as a fresh stream.

generation-output: unchanged
generation-affects: none
