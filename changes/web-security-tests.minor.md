### Added
- Tests for pages staying fresh after a command-line write (TEST.42) and
  the caches under real threads (TEST.54).
- Tests for who may call what across the whole API, generated from the
  route list (TEST.43), what an API key may do (TEST.44), more than one
  admin (TEST.45), trusted-device and two-factor edge cases (TEST.46),
  oversized requests (TEST.47) and security headers on every kind of
  response (TEST.48), the thinly tested API routes (TEST.49) and Galaxy
  Map URLs combined (TEST.50).

### Changed
- A request body over 2 MB is refused with a 413 (JSON under `/api`, the
  error page elsewhere) before the login checks or the database see it.
  Before, any size was read into memory.
- An API key can no longer make API keys, change credentials, set up or
  turn off two-factor sign-in, or log out: those answer 403 and need a
  browser sign-in, so a leaked key can't make itself a replacement.
- Turning two-factor sign-in off forgets every trusted device of that
  admin; the browser that turned it off gets a new one.

### Fixed
- An admin route whose control database can't be reached answers 503
  "database unavailable" instead of a 500.
- A phenomenon added from the command line (`generate.py phenomenon
  --sector-id N`) now shows up on cached sector pages and Galaxy Map
  tiles; before, it marked nothing as changed, so cached tiles never
  showed it.
- `/galaxy/locate` answers JSON when the API says "not found" (for
  example an unknown database), like `/galaxy/tiles` and `/galaxy/stage`
  already did, instead of the HTML 404 page.
