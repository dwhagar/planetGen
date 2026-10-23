**Headline:** The uncharted black hole page test now says what the page answered when it fails.

### Fixed
- **`test_an_uncharted_black_hole_page_names_its_class` failed twice in loaded full gates and never alone (TEST.284).** The page was rendered without "Black Hole (Stellar)" and the failure printed nothing that said why. Run 150 times beside a six-worker load on the Bug Hunt lane it passed every time, so the cause is not found yet; the assertion now carries the status code, the headers and the start of the page, which tell a loading page (202), a redirect, a cached page of another row or a missing field apart on the next sighting.

TODO items resolved: none. Added: none.
