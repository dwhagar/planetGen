**Headline:** The Galaxy Map's node path has a benchmark that fails when a dive, a flight or a climb stalls.

### Added
- **A benchmark of the Galaxy Map over a large synthetic galaxy (MAP.249).** `large_site` serves 46,000 stars with
  every node carrying its real star budget; `tests/test_web_browser_map_perf.py` dives from the whole galaxy to a sector,
  flies across a level and climbs out, measuring each step's settle time, node requests and bytes, stars written to the
  buffer, frame times and longest task, and fails a step over its limit. The numbers, the limits and a profile of where
  the main thread's time goes are in `docs/design/interface/galaxy-map-performance.md`.
