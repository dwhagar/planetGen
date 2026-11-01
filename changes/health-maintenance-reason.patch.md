**Headline:** The health check now says why the site is in maintenance mode.

### Fixed
- **`/api/health` left the reason for maintenance empty when the switch had none of its own.** The `maintenance` block now carries the line visitors see (the standard message when none was given), plus `eta` and `since`; existing fields are unchanged.
