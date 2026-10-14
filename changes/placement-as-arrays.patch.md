### Changed
- Placing a system or compact remnant in a sector now tests the new position against every neighbour as one NumPy array instead of one distance call each; the same random numbers are drawn in the same order, so the placed positions and the stored galaxy are identical (PERF.87).
- The ring centerline length is cached like the ring and slot counts already were (PERF.90).
