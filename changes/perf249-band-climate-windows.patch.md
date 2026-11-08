**Headline:** A world that spins very fast no longer makes the generator split its temperature record into hundreds of millions of empty windows; the numbers do not change (PERF.249).

### Changed

- `temperature_coverage.band_climate` cut the record into one window per solar day, so a planet that turns thousands of times per sample (a hypergiant's) was split into millions of empty windows on every call, costing minutes and gigabytes. The window count is now held at the number of samples (the extra windows were empty and dropped), and the windows are reduced in one pass. The result is bit-identical (tested against the old windowing), so no epoch.
