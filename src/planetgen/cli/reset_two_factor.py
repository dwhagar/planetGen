# planetgen.cli.reset_two_factor

"""
Turns off one admin's two-factor sign-in (SEC.26), for a lost phone and
lost recovery codes:

    python3 -m planetgen.cli.reset_two_factor admin

The admin signs in with the password alone and can set two-factor up
again. The setting is a row in the control database, which uses the same
MySQL settings as every other script (`config.json`, `PLANETGEN_MYSQL_*`,
`--mysql-*`). Written to the activity log as `DB totp.reset`.
"""

import argparse
import sys

import pymysql

from planetgen.admin import activity_log, auth
from planetgen.db.store import add_mysql_connection_args, control_mysql_config, get_control_connection, \
    mysql_config_from_args


def main(argv=None):
    parser = argparse.ArgumentParser(description="Turn off this admin's two-factor sign-in, for a lost phone.")
    parser.add_argument("user", metavar="USER", help="The admin's username.")
    add_mysql_connection_args(parser)
    args = parser.parse_args(argv)
    try:
        conn = get_control_connection(control_mysql_config(mysql_config_from_args(args)))
    except pymysql.MySQLError as exc:
        print(f"error: could not open the control database ({exc}).", file=sys.stderr)
        return 1
    try:
        row = conn.execute("SELECT id FROM admin_users WHERE username = ?", (args.user,)).fetchone()
        if row is None:
            parser.error(f"no admin named {args.user!r}.")
        was_on = auth.disable_totp(conn, row["id"])
        activity_log.event("DB", "totp.reset", user=activity_log.login_name(), target=f"user:{args.user}")
        print(f"Two-factor sign-in for {args.user} is off{'' if was_on else ' (it was not set up)'}.")
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
