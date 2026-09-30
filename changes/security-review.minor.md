### Security
- **No more published first login.** A new install's admin gets a random
  password, printed once by `migrateDb.py` (so by `install.sh` and
  `update.sh`); it must still be changed at first login. An existing
  install whose admin is still on the old `admin`/`password` login gets a
  random password the same way on the next update. See `docs/api.md`
  ("The first admin login", "Resetting the admin login").
- **Changing the admin username or password logs out every other
  browser** signed in as that admin. API keys keep working.
- **Login takes the same time whether or not the username exists.**
- **Rate limits on the web pages**, per IP: search 30 a minute, the Galaxy
  Map 60, its tiles 600, `/api/health` 60 and every other page 300,
  configurable as `ratelimit.pages` in `config.json` (an empty value
  turns one off).
- **CSRF tokens are tied to the login session**, and a new API key's
  one-time display cookie is only sent to `/admin`.
- **Installer permissions:** the code under `src/html` is now owned by
  root and read-only for Apache; only the runtime directories belong to
  Apache's user. `config.json` is set to `root:<apache group>`, mode 640.
  The debug log is mode 0660 (CLI users must be in Apache's group to
  append). The root-run cache-directory helper no longer imports code the
  web user could change.
- **The `/tmp` fallbacks** for the jobs and tile-cache directories are
  created private and refused when another user owns them.
- **A sector's wiki link must be an `http` or `https` URL**, and one
  stored earlier with another scheme is no longer shown as a link.
- **`/api/databases` no longer lists the control database**, and `?db=`
  can't select it.
- **Database errors no longer reach visitors**: `/api/health` and the
  read routes answer "database unavailable" and log the detail.
- **HSTS** is sent on HTTPS responses, and `%`/`_` in the `star_type`
  filter are matched literally.
