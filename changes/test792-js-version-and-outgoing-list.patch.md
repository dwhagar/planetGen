**Headline:** The nebula layer of the Galaxy Map loads its sibling script with its own version, and the Gaia download is listed as an outgoing request that sends the project's user agent.

### Changed

- `galaxynebulae.js` imported `galaxyblocks.js` with a static import, which drops the `?v=` query so the sibling was cached unversioned (and could load twice); it now imports it the way every other script does, with `VERSION_QUERY`.
- `tests/test_user_agent.py` lists `datasets/gaia_nearby.py` among the outgoing request paths; it already sent `user_agent()` and no personal details, and the test now checks that it does.

generation-output: unchanged
