**Headline:** The route fuzz test no longer calls a nebula volume (bytes) or an emptied page parameter a failure, so it passes again.

### Fixed
- `test_fuzz_web_routes` (TEST.588) reads a response as bytes decoded with replacement, because the nebula volume route answers `application/octet-stream`; it was failing on a 200 with a binary body, not on a route error.
- An empty page parameter on a list page (`/systems?standalone_page=`) redirects to the clean address, as the search page does; the test follows the redirect and requires the clean page to answer 200. No route changes: no hostile or out-of-range page number gets an error page.
