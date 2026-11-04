**Headline:** A browser test now runs the question the value filters were built for, planets of 0.8 to 1.2 g with a surface dose under 2 mSv a year, and finds the same planets on the Search page and in `GET /api/search`.

### Added
- `tests/test_web_search_values_browser.py` sets nine planets of a sample galaxy to known gravities and doses (on both edges of the gravity range, just outside it, over the dose limit, and with a value missing), asks the question by its own numbers, and checks that the page, its chips and the API's rows, `values` and `active_ranges` (with `g0` and `mSv/yr`) give exactly the three planets that fit.

generation-output: unchanged
TODO items resolved: UX.225 (the filters, sorts and units of `/api/search` shipped with UX.226 and are in `docs/reference/api.md`; this adds the one-query browser test)
