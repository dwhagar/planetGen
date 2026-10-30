**Headline:** The command line and the run log now say what every number they print counts or measures.

### Changed
- **The skeleton line of a plan run** reads "Building skeleton: edge 12000 pc, stellar mass 5.2e10 solar masses, disk scale length 2600 pc, ... 4 spiral arms, pitch angle 12.3 degrees ..." instead of `edge_pc=` and `k_norm=` style pairs, and the three counts of the "Drew" line each end in "systems".
- **`planetgen near`** ends with "5 objects found within 10 pc; 3 more objects in range are not generated yet", and **`planetgen planets` and `planetgen moons`** end a page with "100 planets shown" or "100 moons shown" before the `--after` to carry on with.
- **The verbose log of a world's life** gives its complexity as "x Earth's net primary production" with the yield, rate factor and grade as "of 1" and "of 3", and its schedule as ages in Gyr.
