**Headline:** A sphere in the galaxy can now be turned into at most 256 runs of cell keys, the first half of searching it by key.

### Added
- `planetgen.galaxy.cellkey.sphere_runs(center, r, budget=256)`: the ascending, merged runs of cell keys a sphere can hold (an octree cover in numpy and integers, no numba, so the web app can use it), for any stored reach, with a budget of runs. A sphere that reaches past the box gets the runs of the part inside it; one wholly outside gets none. The tests check every point of random spheres, the surface included, against the key of `cell_keys`, and the run list is ascending, not touching, and within the budget.

TODO items resolved: GEN.276 / Added: none
