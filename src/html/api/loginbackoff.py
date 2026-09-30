# html/api/loginbackoff.py

"""
Per-username login backoff, on top of `POST /api/auth/login`'s per-IP
rate limit (`auth.LOGIN_RATE_LIMIT`).

The per-IP limit stops one address from guessing quickly, but an attacker
with many addresses could still guess one username's password at ten
tries a minute per address. This counts failed logins per username
instead: the first `FREE_FAILURES` cost nothing, then each further
failure locks that username for twice as long as the last (1 s, 2 s,
4 s, ...), up to `MAX_LOCK_SECONDS`. While a username is locked, a login
for it is refused with a 429 before the password is even checked, so a
correct guess during the lock gets nowhere. A successful login clears the
count, and a username with no failure for `FORGET_AFTER_SECONDS` starts
over.

Unknown usernames are counted exactly like real ones, so the lock never
reveals which usernames exist. The trade-off is that anyone can lock a
real admin out for up to `MAX_LOCK_SECONDS` at a time by failing on
purpose; the cap keeps that short.

Kept in memory, per process, like the per-IP limit's default storage:
a multi-worker server tracks each worker separately (the counts are
still bounded by `MAX_TRACKED_USERNAMES` per worker).
"""

import math
import threading
import time

FREE_FAILURES = 10
"""int: Failed logins a username gets before the first lock -- the same
as the per-IP limit's ten a minute, so one address guessing alone hits
that limit first."""

BASE_LOCK_SECONDS = 1.0
"""float: The first lock's length; each later failure doubles it."""

MAX_LOCK_SECONDS = 900.0
"""float: The longest lock (15 minutes)."""

FORGET_AFTER_SECONDS = 3600.0
"""float: A username whose last failure is this old (and not locked)
starts over at zero."""

MAX_TRACKED_USERNAMES = 10000
"""int: At most this many usernames are tracked; past it the entries
idle longest are dropped, so random usernames can't grow memory without
bound."""


def normalize(username):
    """The key a username is counted under: stripped and case-folded,
    so `Admin` and `admin ` share one count."""
    return username.strip().casefold()


class LoginBackoff:
    """
    Failed-login counts and locks, keyed by `normalize(username)`.
    Thread-safe; `clock` (seconds, monotonic) is injectable for tests.
    """

    def __init__(self, clock=time.monotonic):
        self._clock = clock
        self._lock = threading.Lock()
        # key -> [failures, locked_until, last_failure_at]
        self._entries = {}

    def retry_after(self, username):
        """Seconds until `username` may try again (rounded up to a whole
        second), or 0 when it isn't locked."""
        now = self._clock()
        with self._lock:
            entry = self._current(normalize(username), now)
            if entry is None or entry[1] <= now:
                return 0
            return max(1, math.ceil(entry[1] - now))

    def record_failure(self, username):
        """Counts one failed login; returns the new lock's length in
        seconds (0 while the failure is still free)."""
        key = normalize(username)
        now = self._clock()
        with self._lock:
            entry = self._current(key, now)
            if entry is None:
                entry = [0, 0.0, now]
                self._entries[key] = entry
                self._prune(now)
            entry[0] += 1
            entry[2] = now
            extra = entry[0] - FREE_FAILURES
            if extra <= 0:
                return 0.0
            seconds = min(MAX_LOCK_SECONDS, BASE_LOCK_SECONDS * 2 ** min(extra - 1, 32))
            entry[1] = now + seconds
            return seconds

    def record_success(self, username):
        """Clears `username`'s count after a successful login."""
        with self._lock:
            self._entries.pop(normalize(username), None)

    def clear(self):
        """Forgets every username (tests)."""
        with self._lock:
            self._entries.clear()

    def _current(self, key, now):
        entry = self._entries.get(key)
        if entry is not None and entry[1] <= now and now - entry[2] >= FORGET_AFTER_SECONDS:
            del self._entries[key]
            return None
        return entry

    def _prune(self, now):
        if len(self._entries) <= MAX_TRACKED_USERNAMES:
            return
        for key in [k for k, e in self._entries.items() if e[1] <= now and now - e[2] >= FORGET_AFTER_SECONDS]:
            del self._entries[key]
        overflow = len(self._entries) - MAX_TRACKED_USERNAMES
        if overflow > 0:
            # Unlocked usernames go first, oldest failure first, so a
            # flood of new names can't free a locked one.
            oldest = sorted(self._entries, key=lambda k: (self._entries[k][1] > now, self._entries[k][2]))
            for key in oldest[:overflow]:
                del self._entries[key]


backoff = LoginBackoff()
"""LoginBackoff: The one instance `auth.login` uses."""
