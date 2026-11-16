**Headline:** The test that haze never cools the ground less than the clear air now replays the same world and skips airs that did not settle, so it tests something and is green again.

### Fixed

- TEST.778: `test_taking_the_haze_out_of_the_air_warms_the_ground_and_never_cools_it` was red on main since GEN.1754 (29fb3c4f1c) moved the water of rocky bodies and with it the first hazy body with ground, seed 98 (an earlier note) or 29 to seed 15, a 257 K ice-covered body whose air did not settle (`Air.converged` False, the knee of GEN.1344: the haze switches on and off between two answers and the stored air is their mean). The bound of 0.5 K is unchanged. Two causes, both of the test: the fixture takes an air that settles (`converged` and not `knee`), and the replay repeats every draw of the body's own solve from `Air.entry_state` (as `air_branches.other_branch` does); before it drew on a stream of its own and ended with a replayed haze of 0.0 in every body tried (seeds 15, 29, 140, 189), so the two airs it compared did not differ and the test passed on nothing. Now the first body is seed 29 (143.60 K with the haze, 144.09 K without: the haze shields the ground by 0.49 K, anti-greenhouse as GEN.851 says) and the test asserts the replay holds the haze (`tau` over 0). Seed 189 gives 117.27 K with and 132.14 K without.

generation-output: unchanged
generation-affects: none
