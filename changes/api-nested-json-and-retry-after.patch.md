### Fixed
- A 100,000-deep nested JSON body is a 400 ("nested too deeply") instead of a 500 (API.20), and `Retry-After` is only sent on a 429, no longer on every response Flask-Limiter counts (API.21).
