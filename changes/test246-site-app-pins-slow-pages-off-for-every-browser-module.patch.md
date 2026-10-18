**Headline:** Every browser test module that builds its pages from the shared site app now has the slow-page rule pinned off, so a loaded run no longer fails a module's setup with a 202 for /search.

### Fixed
- Tests: only `test_web_a11y.py` pinned `tuning.SLOW_PAGE_SECONDS` off (TEST.232); the layout, crumbs, edit menu and other browser modules import the same `site_app` and `page_targets` without that module's autouse fixture, so when a loaded machine made the first build of /search slower than 3 s the route was marked slow, the next request answered 202 with the loading state, and the module's `page_targets` fixture failed (`web.search: GET /search answered 202 (anonymous) / 202 (admin)`, Core 3's layout record). The pin now lives in `site_app` itself, for every module that uses it, and a test in the layout module pins it.

TODO items resolved: TEST.246 / Added: none
