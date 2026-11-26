**Headline:** The giant-tilt test reads Sun-like stars, so it finds its four gas giants again.

### Fixed
- `test_a_generated_giant_tail_reaches_past_sixty_degrees_for_an_ice_giant_only` searched seeds 7100 to 7400 of random stars and ran out of seeds with one gas giant (a random star is mostly a red dwarf, and the star mix moved with the Train 62 mass function). The sweep now asks for a Sun-like star (G2V, no companion): 100 systems hold 10 gas giants and 29 giants (8 and 28 before Train 62), so the giant rate did not fall, and the test passes in 24 minutes under load instead of reading 300 systems for 50. The clause and the four giants of each kind it needs are unchanged.

generation-output: unchanged
