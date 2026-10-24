**Headline:** The object page walk test reads the bodies it opens from the database, not from a loop of API calls that could run out of the API's rate limit.

### Fixed
- **The `test_object_page_walk.py` fixtures errored with a system detail that had no `"stars"` (TEST.359).** The fixture asked the API about each system (two calls a system, up to every system of the sample galaxy); on a machine where the per-IP limit of 50 an hour was already spent, the 429 answer was read as a system and had no `stars`. The fixture now selects the first six rows of each kind (star, planet, moon, belt, comet, facility) from the tables and prints their IDs, makes no API call, and fails with a plain message if the galaxy holds no star or planet. The module passes in 57 s.
TODO items resolved: TEST.359. Added: none.
