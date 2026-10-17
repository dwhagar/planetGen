# planetgen.cli.forget_devices

"""
Revokes every trusted-device cookie of one admin (SEC.22), so that admin's
browsers are asked for the two-factor code again:

    python3 -m planetgen.cli.forget_devices admin

The devices are rows in the control database, which uses the same MySQL
settings as every other script (`config.json`, `PLANETGEN_MYSQL_*`,
`--mysql-*`). Written to the activity log as `DB devices.revoke`.
"""

import argparse
import sys

import pymysql

from planetgen.admin import activity_log, auth
from planetgen.db.store import add_mysql_connection_args, control_mysql_config, get_control_connection, \
    mysql_config_from_args


def main(argv=None):
    parser = argparse.ArgumentParser(description="Revoke every trusted-device cookie of this admin.")
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
        revoked = auth.revoke_devices(conn, row["id"])
        activity_log.event("DB", "devices.revoke", user=activity_log.login_name(), target=f"user:{args.user}",
                           revoked=revoked)
        print(f"Revoked {revoked} trusted device{'' if revoked == 1 else 's'} of {args.user}.")
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
