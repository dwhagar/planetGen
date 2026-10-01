# src/loginLockouts.py

"""
Lists and lifts login lockouts (SEC.1, SEC.21) from the command line, for
an admin locked out of the web interface itself:

    python3 src/loginLockouts.py                     # list what is locked
    python3 src/loginLockouts.py --ip 203.0.113.5    # lift one address
    python3 src/loginLockouts.py --user admin        # lift one username
    python3 src/loginLockouts.py --all               # lift every lockout
    python3 src/loginLockouts.py --forget-devices admin
                                    # that admin's browsers lose their
                                    # trusted-device cookies (SEC.22)
    python3 src/loginLockouts.py --reset-two-factor admin
                                    # turns off that admin's two-factor
                                    # sign-in (lost phone and recovery
                                    # codes; SEC.26)

An IPv6 address is lifted by its /64, as it is counted. Lifting also
forgets the count (and an address's doubling level). Each lift is written
to the activity log as `DB lockout.lift`. Uses the same MySQL settings as
every other script (`config.json`, `PLANETGEN_MYSQL_*`, `--mysql-*`).
"""

import argparse
import getpass
import sys
import time

import pymysql

from stellarObjects import activitylog, adminAuth, loginThrottle
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
    group.add_argument("--forget-devices", metavar="USER",
                       help="Revoke every trusted-device cookie of this admin (SEC.22).")
    group.add_argument("--reset-two-factor", metavar="USER",
                       help="Turn off this admin's two-factor sign-in, for a lost phone (SEC.26).")
    add_mysql_connection_args(parser)
    args = parser.parse_args(argv)

    try:
        conn = get_control_connection(control_mysql_config(mysql_config_from_args(args)))
    except pymysql.MySQLError as exc:
        print(f"error: could not open the control database ({exc}).", file=sys.stderr)
        return 1
    try:
        store = loginThrottle.DbStore(conn)
        if args.reset_two_factor:
            row = conn.execute("SELECT id FROM admin_users WHERE username = ?", (args.reset_two_factor,)).fetchone()
            if row is None:
                parser.error(f"no admin named {args.reset_two_factor!r}.")
            was_on = adminAuth.disable_totp(conn, row["id"])
            activitylog.event("DB", "totp.reset", user=_who(), target=f"user:{args.reset_two_factor}")
            print(f"Two-factor sign-in for {args.reset_two_factor} is off"
                  f"{'' if was_on else ' (it was not set up)'}.")
            return 0
        if args.forget_devices:
            row = conn.execute("SELECT id FROM admin_users WHERE username = ?", (args.forget_devices,)).fetchone()
            if row is None:
                parser.error(f"no admin named {args.forget_devices!r}.")
            revoked = adminAuth.revoke_devices(conn, row["id"])
            activitylog.event("DB", "devices.revoke", user=_who(), target=f"user:{args.forget_devices}",
                              revoked=revoked)
            print(f"Revoked {revoked} trusted device{'' if revoked == 1 else 's'} of {args.forget_devices}.")
            return 0
        if args.all:
            lifted = store.lift()
            target = "all"
        elif args.ip is not None:
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
