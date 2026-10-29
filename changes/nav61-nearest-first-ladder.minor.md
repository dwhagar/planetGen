**Headline:** What's nearby no longer chokes on the dense core: a sphere that holds more than 100,000 systems is refused with a message that says to shrink the distance, and "the nearest few" is found by trying 2, 4, 8, 16 and 32 pc in turn.

### Added
- **Nearest first by a radius ladder (NAV.61).** `planetgen.db.near.nearest` and the API's `GET /api/near?...&nearest=K` (and `planetgen query near --nearest K`) return the nearest K rows by asking for 2, 4, 8, 16 and 32 pc and stopping at the first sphere that holds K, so the cost follows the neighbourhood: quick at the galactic core, where a 50 pc sphere holds millions of systems, as well as at the Sun. `distance` then only caps how far it looks.
- **Early-stop count and the dense-core guard.** Before a search reads any row it counts the sphere with `region.capped_count` (a `COUNT(*)` over the key-run statement stopped by `LIMIT`). A sphere over 100,000 systems is refused (a 400 on the API) with the text that says to shrink the distance or ask for the nearest few. Every answer carries `systems_in_range` (counted no further than 10,000) and `systems_in_range_capped`, and the What's nearby page shows "10,000+ star systems in range" when the sphere holds more.

TODO items resolved: NAV.61. Added: none.
