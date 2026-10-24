**Headline:** The browser tests' shared sample fixture no longer errors at setup on a fresh galaxy.

### Fixed
- **Twenty browser tests errored at setup with `KeyError: 'planets'` (TEST.349).** The shared `sample_params` fixture of `test_web_a11y.py` asked `/api/systems/<id>` for the detail of up to 100 systems to find one with planets. The API answers 50 requests an hour (`RATELIMIT_DEFAULT`), the sample galaxy holds about 375 systems in its sector and most have no planets, so after 50 calls every answer was `429 rate limit exceeded`, read as "no planets", and the next line indexed a body that was an error. The fixture now finds the system with one SQL query and makes one API call.
TODO items resolved: TEST.349. Added: none.
