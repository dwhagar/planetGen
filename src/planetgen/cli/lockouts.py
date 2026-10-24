# planetgen.cli.lockouts

"""
Lists the login lockouts (SEC.1, SEC.21) from the command line, for an
admin locked out of the web interface itself:

    python3 -m planetgen.cli.lockouts       # list what is locked

It only lists. `planetgen.cli.lift_lockout` lifts a lockout,
`planetgen.cli.forget_devices` revokes an admin's trusted-device cookies
(SEC.22) and `planetgen.cli.reset_two_factor` turns off an admin's
two-factor sign-in (SEC.26); each does that one thing.

The lockouts live in Redis, on the rate-limit storage (`config.json`'s
`ratelimit.storage_uri`, else `redis.url`; SEC.30).
"""

import argparse
import sys
import time

import redis

from planetgen.admin import throttle
from planetgen.api.config import RATELIMIT_KEY_PREFIX, ratelimit_storage_uri
from planetgen.api.loginguard import REDIS_SCHEMES
from planetgen.util import log


def redis_store():
    """The lockouts' Redis store, or `None` after printing why there is
    none: the rate-limit storage keeps them in each web process's memory,
    and there is nothing to reach from here."""
    uri = ratelimit_storage_uri()
    if not uri.startswith(REDIS_SCHEMES):
        print(f"error: the rate-limit storage is {uri!r}, which keeps the lockouts in each web process's "
              f"memory; there is nothing to list or lift from here.", file=sys.stderr)
        return None
    return throttle.RedisStore.from_url(uri, RATELIMIT_KEY_PREFIX)


def main(argv=None):
    parser = argparse.ArgumentParser(description="List the login lockouts (per address and per username).")
    log.add_debug_arguments(parser)
    log.apply_debug_arguments(parser.parse_args(argv), parser)
    store = redis_store()
    if store is None:
        return 1
    try:
        rows = throttle.locked_subjects(store)
    except redis.RedisError as exc:
        print(f"error: could not use Redis at {ratelimit_storage_uri()} ({exc}).", file=sys.stderr)
        return 1
    if not rows:
        print("Nothing is locked.")
    for row in rows:
        until = time.strftime("%Y-%m-%d %H:%M:%SZ", time.gmtime(row["locked_until"]))
        print(f"{row['scope']:4}  {row['subject']:40}  until {until} ({row['retry_after']} s left)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
