# html/api/limiter.py

"""
Rate limiting for the planetGen API (Flask-Limiter).

One shared `Limiter` instance, created here uninitialized (no `app` bound
yet) and wired up in `app.py`'s `create_app` via `limiter.init_app(app)`
-- the standard Flask-extension pattern, needed because `routes.py`'s
individual view functions want to decorate themselves with
`@limiter.limit(...)` for a stricter per-route limit, and importing a
single shared instance is how both modules end up talking about the same
limiter rather than each creating their own.

The global default limits and storage backend are configured through
`app.config` (`RATELIMIT_DEFAULT`/`RATELIMIT_STORAGE_URI`, see
`config.py`) rather than passed to the constructor here, since `Limiter`
reads those from the app at `init_app` time -- this module doesn't need
its own copy of that configuration.
"""

from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

limiter = Limiter(key_func=get_remote_address)
