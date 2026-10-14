### Changed
- **Saving a sector finds each system's nearest neighbours with one spatial tree (PERF.81).** It used to measure the distance from every system to every other one (1.8 million measurements for a sector of 1,916 systems, about 10 s); one cKDTree now answers all of them in 0.04 s, with the same neighbours, ties included.
