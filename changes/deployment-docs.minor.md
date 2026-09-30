### Added
- **Deployment guides for every major platform.** `docs/deployment/` has
  a comparison table and a guide for each setup: Apache2 + mod_wsgi (the
  existing guide, moved from `docs/apache-deployment.md`), nginx +
  gunicorn and Caddy + gunicorn on Linux, three native Windows setups
  (IIS + HttpPlatformHandler, Caddy, Apache Lounge, each with waitress)
  plus WSL2, Homebrew nginx + gunicorn under launchd on macOS (macOS
  Server is discontinued), and a VPS and platform-as-a-service note.
  Example configs are in `examples/nginx/`, `examples/caddy/`,
  `examples/systemd/`, `examples/windows/` and `examples/macos/`.
- **`proxy_fix` setting for running behind a reverse proxy.** Behind
  nginx, Caddy, IIS or Apache's `mod_proxy`, the app saw the proxy's
  address for every visitor (so everyone shared one rate-limit budget)
  and never sent `Strict-Transport-Security`. `config.json`'s
  `proxy_fix` (`x_for`, `x_proto`, `x_host`; environment variables
  `PLANETGEN_PROXY_FIX_X_FOR`, `_X_PROTO`, `_X_HOST`) sets how many
  proxies to trust for each `X-Forwarded-*` header, applied with
  werkzeug's `ProxyFix`. All 0 by default, so Apache + mod_wsgi is
  unchanged. A value that isn't a whole number 0 or above stops the app
  at startup.

### Changed
- `docs/server-checklist.md` works for every platform, and its schema
  and index checks match the current schema.
- `docs/api.md` no longer suggests a read-only database account for the
  web app: the app writes through its one account (admin pages, write
  endpoints), as `docs/config.md` already said.
- Docs and code comments that named the old `sectorGen.py`,
  `systemGen.py`, `phenomenonGen.py`, `galaxyGen.py` and `galaxyPlan.py`
  scripts now name the `generate.py` subcommands.
- `docs/TODO.md` item 54 tracks the admin Generate page's POSIX-only job
  handling, which doesn't work reliably on native Windows.
