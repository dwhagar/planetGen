# planetgen.cli.lift_lockout

"""
Lifts login lockouts (SEC.1, SEC.21) from the command line, for an admin
locked out of the web interface itself:

    python3 -m planetgen.cli.lift_lockout --ip 203.0.113.5    # one address
    python3 -m planetgen.cli.lift_lockout --user admin        # one username
    python3 -m planetgen.cli.lift_lockout --all               # every lockout

One of the three is needed. An IPv6 address is lifted by its /64, as it is
counted. Lifting also forgets the count (and an address's doubling level).
Each lift is written to the activity log as `DB lockout.lift`. The
lockouts live in Redis, on the rate-limit storage (`config.json`'s
`ratelimit.storage_uri`, else `redis.url`; SEC.30). `planetgen.cli.lockouts`
lists what is locked.
"""

import argparse
import sys

import redis

from planetgen.admin import activity_log, throttle
from planetgen.api.config import ratelimit_storage_uri
from planetgen.cli import lockouts
from planetgen.util import log


def main(argv=None):
    parser = argparse.ArgumentParser(description="Lift a login lockout (an address, a username, or every one).")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ip", help="Lift the lockout of this address (an IPv6 address: its /64).")
    group.add_argument("--user", help="Lift the lockout of this username.")
    group.add_argument("--all", action="store_true", help="Lift every lockout.")
    log.add_debug_arguments(parser)
    args = parser.parse_args(argv)
    log.apply_debug_arguments(args, parser)

    scope = subject = None
    if args.ip is not None:
        scope, subject = throttle.SCOPE_IP, throttle.ip_subject(args.ip)
        if subject is None:
            parser.error(f"{args.ip!r} is not an address that can be locked (invalid or loopback).")
    elif args.user:
        scope, subject = throttle.SCOPE_USER, throttle.normalize_username(args.user)
    store = lockouts.redis_store()
    if store is None:
        return 1
    try:
        lifted = store.lift(scope, subject)
    except redis.RedisError as exc:
        print(f"error: could not use Redis at {ratelimit_storage_uri()} ({exc}).", file=sys.stderr)
        return 1
    target = f"{scope}:{subject}" if scope else "all"
    activity_log.event("DB", "lockout.lift", user=activity_log.login_name(), target=target, lifted=lifted)
    print(f"Lifted {lifted} lockout{'' if lifted == 1 else 's'} ({target}).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
