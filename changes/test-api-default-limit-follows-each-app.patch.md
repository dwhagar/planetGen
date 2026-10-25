**Headline:** The browser tests no longer fail with "50 per 1 hour" rate-limit errors, and each app now follows its own API request limit.

### Fixed
- **The API's default request limit (200 per day, 50 per hour per address) was taken from the first app made in a process and kept for every app made after it** (`Limiter.init_app` keeps only the first default of a shared limiter). A test run builds many apps in one process, so a browser fixture that called the API from one address (`sample_params` makes about 7 calls per module) could not be given a different limit and ran into the hourly 50, which errored 259 browser tests at setup. The limit is now read from the current app's `RATELIMIT_DEFAULT` on each request. On a real server, which makes one app, nothing changes: the limit is still 200 per day and 50 per hour, and an empty setting still turns it off.
- The browser test sites (`habitability_site`, `map_site_support`, `test_web_a11y`, `test_web_browser_tables`, `test_web_settings_browser`) now set `RATELIMIT_DEFAULT = API_LIMIT_OFF`, so their fixtures cannot be throttled. A route's own limits (login, writes) still apply to them.

### Added
- `planetgen.api.limiter.API_LIMIT_OFF`, the `RATELIMIT_DEFAULT` that turns the default limit off for a test harness.
- Tests in `test_web_security_limits.py` that a direct API call is still limited by the default (429 with `rate limit exceeded`, per address), that `API_LIMIT_OFF` is never limited, that two apps in one process keep their own default in either order, and that the production default is unchanged.
TODO items resolved: TEST.385 (TODO closes it after reading this). Added: none.
