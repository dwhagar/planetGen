**Headline:** The plan estimate for a neighbourhood run is pinned by a test against the run that showed it wrong.

### Fixed
- **PERF.178: the plan estimate said 33 GB and 2 min 33 s for a run that took 4 min 28 s.** The pieces were already in: the galaxy-wide scatter tables are left out of the size per system (a few hundred megabytes for the 3,800 systems of that run, not 33 GB), the line no longer says "Refused ... Going ahead" (#2831), and a time worked out from the default speed reads as the range it lies within with "not timed yet on this server". A new test runs the logged case through the estimate and checks all three: the size is small, the range holds the 4 min 28 s the run took, and the summary says the time is not timed yet. No stored value changes.
