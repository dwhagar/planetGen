### Changed
- Every random draw in generation goes through one wrapper (`planetgen/util/draw.py`) whose helpers are built on Python's `random()` alone and draw from the running unit's seeded stream, so a seed gives the same galaxy on any Python release, hash seed or locale, and changing the random source later means editing one file. This changes the numbers a given seed generates (GEN.56).
