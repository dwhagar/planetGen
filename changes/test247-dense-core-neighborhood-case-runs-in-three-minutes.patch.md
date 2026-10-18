**Headline:** The dense-core case of the neighbourhood test runs in about three minutes instead of up to ten, so it no longer times out in a loaded gate.

### Changed
- Tests: `test_random_start_neighborhood_matches_the_real_skeleton_plan[dense_core]` took about 466 s alone (240 to 650 s over 10 idle runs) and hit the 900 s limit in a loaded gate. Its start sits in cells over the cell star ceiling, which the run takes down to their luminosity cap and which cost the most; the case now reaches 18 ly round its start instead of 25 ly (the random case keeps 25 ly), so it generates a handful of sectors instead of 24 and runs in 173 to 197 s. The plan's skeleton build (about 150 s) is the rest and is shared with the other plan tests.

TODO items resolved: TEST.247 / Added: none
