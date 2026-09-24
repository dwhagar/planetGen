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

from flask import request
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

IN_PROCESS_ENVIRON_KEY = "planetgen.in_process"
"""str: WSGI environ key `web/transport.py` sets on the API requests the
Flask-served pages (`html/web/`) dispatch to this app in-process. Not an
`HTTP_*` key, so no client can set it from outside."""


def is_in_process_call():
    """
    True for an API request made in-process by one of the Flask-served
    pages (`web/transport.py`) rather than by a real client.

    Such a call is exempt from the *default* limits only: a page view
    already fetches its data this way (the CGI pages did the same over
    loopback HTTP), and counting each of those against the viewer's
    hourly budget would throttle ordinary browsing. A route's own
    explicit `@limiter.limit(...)` (login, every write) still applies,
    keyed by the real viewer's address, which `web/transport.py` passes
    through as `REMOTE_ADDR`.
    """
    return bool(request.environ.get(IN_PROCESS_ENVIRON_KEY))


limiter = Limiter(key_func=get_remote_address, default_limits_exempt_when=is_in_process_call)
