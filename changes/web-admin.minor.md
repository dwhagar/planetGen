### Changed
- **The admin pages moved to the Flask app:** `/login`, `/logout`,
  `/account` (change username/password), `/admin` (API keys, a sector's
  wiki link; `?keys_page=N`) and `/admin/stats` (server and database
  stats; `?names_page=N`), replacing `login.py`, `logout.py`,
  `changecreds.py`, `admin.py` and `adminstats.py`, which now answer a
  `301` to their new URL. Sessions and the login rate limit work as
  before: the API's own session cookie is relayed unchanged.
- Visiting an admin page while logged out goes to `/login?next=<page>`
  and back there after login (only ever to a page on this site).
- `/admin/stats` and the wiki-link form use the site's configured
  database; the database picker and field are gone.

### Security
- Every admin form (log in, log out, change credentials, create or
  revoke a key, set a wiki link) is a POST with a CSRF token, answered
  with a redirect. A new API key reaches the page after the redirect in
  a signed, HttpOnly, SameSite=Strict one-time cookie, never the URL.
- Logging out is a POST: `GET /logout` only shows a "Log out" button,
  so a link or prefetch can't end a session.
- Admin pages are sent with `Cache-Control: no-store`.
