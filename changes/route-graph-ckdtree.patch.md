### Changed
- **Route graph building is up to 18 times faster (NAV.52).** The nearest-neighbour search and the island joining behind NAV routes now use SciPy's cKDTree instead of a pure-Python k-d tree. On 106,529 systems in 2,000 clusters it takes 3.2 s instead of 57 s, and the graph is the same edge for edge and distance for distance (also checked at 3,000, 10,000 and 20,000 systems).
