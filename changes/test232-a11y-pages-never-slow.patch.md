**Headline:** The accessibility tests no longer fail on a loaded machine when the first build of a page takes longer than the slow-page limit.

### Changed
- `tests/test_web_a11y.py` pins the slow-page rule off for its module (`pages_are_never_slow`, an autouse module fixture setting `tuning.SLOW_PAGE_SECONDS` to infinity). On a loaded run the first build of `/search` or `/phenomena` took over 3 s, the cost rule marked the route slow, and the next request was queued and answered 202 with the loading state instead of the page, which failed the `page_targets` fixture for 199 tests in one gate. These tests check each page's layout and accessibility, not the slow-page answer, which `tests/test_slowpages.py` covers.

### Added
- Nothing.

TODO items resolved: TEST.232 / Added: none
