**Headline:** Every call to the API now leaves one line in the activity log with who made it, how it came in, the route and the HTTP response code.

### Changed
- **Every `/api/` request is logged (API.15).** `planetgen.log` gets an `API call` line when the answer is ready: the account (the key's owner, the signed-in user, `god` for a command-line tool, or `-`), `via=` (`key`, `session`, `anonymous`, `console`), method, route, path, HTTP status, milliseconds, bytes in and out, the user agent and the API version the client sent. Keys, tokens, passwords and bodies are never written.
- **`X-Request-Id`.** Each answer carries the request id that is in the log line; a client's own well-formed id is kept.
- **`api_log.internal`.** The calls the web pages make in-process (about six per page view) are marked `internal=1`; setting this to `off` leaves them out of the log. Default `all`.
- The command-line tools that lift a lockout, forget devices and reset two-factor write the same line with `user="god"` and `via=console`.
- `docs/reference/api.md` and `docs/reference/config.md` describe the log; `tests/test_api_call_log.py` makes a call each way and checks the line.

TODO items resolved: API.15. Added: none.
