# planetgen.cli.lockouts

"""
Lists and lifts login lockouts (SEC.1, SEC.21) from the command line, for
an admin locked out of the web interface itself:

    python3 -m planetgen.cli.lockouts                     # list what is locked
    python3 -m planetgen.cli.lockouts --ip 203.0.113.5    # lift one address
    python3 -m planetgen.cli.lockouts --user admin        # lift one username
    python3 -m planetgen.cli.lockouts --all               # lift every lockout
    python3 -m planetgen.cli.lockouts --forget-devices admin
                                    # that admin's browsers lose their
                                    # trusted-device cookies (SEC.22)
    python3 -m planetgen.cli.lockouts --reset-two-factor admin
                                    # turns off that admin's two-factor
                                    # sign-in (lost phone and recovery
                                    # codes; SEC.26)

An IPv6 address is lifted by its /64, as it is counted. Lifting also
forgets the count (and an address's doubling level). Each lift is written
to the activity log as `DB lockout.lift`. The lockouts live in Redis, on
the rate-limit storage (`config.json`'s `ratelimit.storage_uri`, else
`redis.url`; SEC.30); the device and two-factor options use the same MySQL
settings as every other script (`config.json`, `PLANETGEN_MYSQL_*`,
`--mysql-*`).
"""

import argparse
import getpass
import sys
import time

import pymysql
import redis

from planetgen.admin import activity_log, auth, throttle
from planetgen.api.config import RATELIMIT_KEY_PREFIX, ratelimit_storage_uri
from planetgen.api.loginguard import REDIS_SCHEMES
from planetgen.db.store import add_mysql_connection_args, control_mysql_config, get_control_connection, \
    mysql_config_from_args


def _who():
    try:
        return getpass.getuser()
    except Exception:  # noqa: BLE001
        return None


def _device_options(parser, args):
    """--forget-devices and --reset-two-factor: rows in the control
    database."""
    try:
        conn = get_control_connection(control_mysql_config(mysql_config_from_args(args)))
    except pymysql.MySQLError as exc:
        print(f"error: could not open the control database ({exc}).", file=sys.stderr)
        return 1
    try:
        if args.reset_two_factor:
            row = conn.execute("SELECT id FROM admin_users WHERE username = ?", (args.reset_two_factor,)).fetchone()
            if row is None:
                parser.error(f"no admin named {args.reset_two_factor!r}.")
            was_on = auth.disable_totp(conn, row["id"])
            activity_log.event("DB", "totp.reset", user=_who(), target=f"user:{args.reset_two_factor}")
            print(f"Two-factor sign-in for {args.reset_two_factor} is off"
                  f"{'' if was_on else ' (it was not set up)'}.")
            return 0
        row = conn.execute("SELECT id FROM admin_users WHERE username = ?", (args.forget_devices,)).fetchone()
        if row is None:
            parser.error(f"no admin named {args.forget_devices!r}.")
        revoked = auth.revoke_devices(conn, row["id"])
        activity_log.event("DB", "devices.revoke", user=_who(), target=f"user:{args.forget_devices}",
                          revoked=revoked)
        print(f"Revoked {revoked} trusted device{'' if revoked == 1 else 's'} of {args.forget_devices}.")
        return 0
    finally:
        conn.close()


def main(argv=None):
    parser = argparse.ArgumentParser(description="List or lift login lockouts (per address and per username).")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--ip", help="Lift the lockout of this address (an IPv6 address: its /64).")
    group.add_argument("--user", help="Lift the lockout of this username.")
    group.add_argument("--all", action="store_true", help="Lift every lockout.")
    group.add_argument("--forget-devices", metavar="USER",
                       help="Revoke every trusted-device cookie of this admin (SEC.22).")
    group.add_argument("--reset-two-factor", metavar="USER",
                       help="Turn off this admin's two-factor sign-in, for a lost phone (SEC.26).")
    add_mysql_connection_args(parser)
    args = parser.parse_args(argv)

    if args.reset_two_factor or args.forget_devices:
        return _device_options(parser, args)
    scope = subject = None
    if args.ip is not None:
        scope, subject = throttle.SCOPE_IP, throttle.ip_subject(args.ip)
        if subject is None:
            parser.error(f"{args.ip!r} is not an address that can be locked (invalid or loopback).")
    elif args.user:
        scope, subject = throttle.SCOPE_USER, throttle.normalize_username(args.user)
    uri = ratelimit_storage_uri()
    if not uri.startswith(REDIS_SCHEMES):
        print(f"error: the rate-limit storage is {uri!r}, which keeps the lockouts in each web process's "
              f"memory; there is nothing to list or lift from here.", file=sys.stderr)
        return 1
    store = throttle.RedisStore.from_url(uri, RATELIMIT_KEY_PREFIX)
    try:
        if scope is None and not args.all:
            rows = throttle.locked_subjects(store)
            if not rows:
                print("Nothing is locked.")
            for row in rows:
                until = time.strftime("%Y-%m-%d %H:%M:%SZ", time.gmtime(row["locked_until"]))
                print(f"{row['scope']:4}  {row['subject']:40}  until {until} ({row['retry_after']} s left)")
            return 0
        lifted = store.lift(scope, subject)
        target = f"{scope}:{subject}" if scope else "all"
    except redis.RedisError as exc:
        print(f"error: could not use Redis at {uri} ({exc}).", file=sys.stderr)
        return 1
    activity_log.event("DB", "lockout.lift", user=_who(), target=target, lifted=lifted)
    print(f"Lifted {lifted} lockout{'' if lifted == 1 else 's'} ({target}).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
