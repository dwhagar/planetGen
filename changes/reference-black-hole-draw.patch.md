**Headline:** The Milky Way-like reference galaxy keeps its 1e8 black holes now that the halo is lighter: its failed-supernova draw moves from the middle of the range to 0.55.

### Changed

- The reference structure (the shape built with no seed) takes its solar failed-supernova share from `tuning.REFERENCE_BLACK_HOLE_DRAWS` (0.55 of the range, the metal-free draw stays at the middle) through `black_hole_formation.reference_failed_shares`, which replaces `middle_failed_shares`. The lighter halo of the conditional chain cut the share from 0.175 to 0.161 (9.3e7 black holes in 5e10 solar masses of stars); it is 0.177 now (1.02e8). Galaxies with a seed draw their own shares as before. generation-output: changed
