# src/loginLockouts.py

"""
Lists and lifts login lockouts (SEC.1, SEC.21) from the command line, for
an admin locked out of the web interface itself:

    python3 src/loginLockouts.py                     # list what is locked
    python3 src/loginLockouts.py --ip 203.0.113.5    # lift one address
    python3 src/loginLockouts.py --user admin        # lift one username
    python3 src/loginLockouts.py --all               # lift every lockout

An IPv6 address is lifted by its /64, as it is counted. Lifting also
forgets the count (and an address's doubling level). Each lift is written
to the activity log as `DB lockout.lift`. Uses the same MySQL settings as
every other script (`config.json`, `PLANETGEN_MYSQL_*`, `--mysql-*`).
"""

import argparse
import getpass
import sys
import time

from stellarObjects import activitylog, loginThrottle
from stellarObjects._db import add_mysql_connection_args, control_mysql_config, get_control_connection, \
    mysql_config_from_args


def _who():
    try:
        return getpass.getuser()
    except Exception:  # noqa: BLE001
        return None


def main(argv=None):
    parser = argparse.ArgumentParser(description="List or lift login lockouts (per address and per username).")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--ip", help="Lift the lockout of this address (an IPv6 address: its /64).")
    group.add_argument("--user", help="Lift the lockout of this username.")
    group.add_argument("--all", action="store_true", help="Lift every lockout.")
    add_mysql_connection_args(parser)
    args = parser.parse_args(argv)

    conn = get_control_connection(control_mysql_config(mysql_config_from_args(args)))
    try:
        store = loginThrottle.DbStore(conn)
        if args.all:
            lifted = store.lift()
            target = "all"
        elif args.ip:
            subject = loginThrottle.ip_subject(args.ip)
            if subject is None:
                parser.error(f"{args.ip!r} is not an address that can be locked (invalid or loopback).")
            lifted = store.lift(loginThrottle.SCOPE_IP, subject)
            target = f"ip:{subject}"
        elif args.user:
            subject = loginThrottle.normalize_username(args.user)
            lifted = store.lift(loginThrottle.SCOPE_USER, subject)
            target = f"user:{subject}"
        else:
            rows = loginThrottle.locked_subjects(store)
            if not rows:
                print("Nothing is locked.")
            for row in rows:
                until = time.strftime("%Y-%m-%d %H:%M:%SZ", time.gmtime(row["locked_until"]))
                print(f"{row['scope']:4}  {row['subject']:40}  until {until} ({row['retry_after']} s left)")
            return 0
    finally:
        conn.close()
    activitylog.event("DB", "lockout.lift", user=_who(), target=target, lifted=lifted)
    print(f"Lifted {lifted} lockout{'' if lifted == 1 else 's'} ({target}).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
