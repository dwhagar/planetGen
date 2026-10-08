### Changed
- A Galaxy Map view now has a stated cap on the stars it can fetch (70,000, about 8 MB before compression), built from the per-tile caps, with a test that fails if a budget is raised past it; a full view of the test galaxy is also timed (MAP.109).
