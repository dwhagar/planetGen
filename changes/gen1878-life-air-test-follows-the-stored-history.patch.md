**Headline:** The life-air API test now says what the route documents: a lifeless world with a stored air history answers from that record.

### Fixed

- GEN.1878: `tests/test_api_life_air.py` expected `life_air` to be `null` for every planet with no living class. The route (`planet_life_air`, API.51) documents `null` only for a world with no living class and no stored history; since GEN.1745 slice 2 a lifeless world with air of 500 Pa or more keeps an air history record, and the route answers from it (its gases, shocks and temperatures, no classes and no stages of life). The API was right and the test stale; the test now pins both cases. No change to the route or to any world.

generation-output: unchanged
