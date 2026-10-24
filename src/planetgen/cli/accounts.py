#!/usr/bin/env python3
# planetgen.cli.accounts

"""
The recovery path for the accounts (USR.31; the full accounts CLI is USR.26):

    python3 -m planetgen.cli.accounts --set-owner USER

Makes the account USER the one Owner, the old Owner becoming an admin. The Owner is the account that changes the
Owner-tier settings (the mail server, the site's address, the From address; ADM.86) after typing the password again.
The accounts are rows in the control database, which uses the same MySQL settings as every other script
(`config.json`, `PLANETGEN_MYSQL_*`, `--mysql-*`). Written to the activity log as `DB owner.set`.
"""

import argparse
import sys

import pymysql

from planetgen.admin import activity_log, auth
from planetgen.db.store import add_mysql_connection_args, control_mysql_config, get_control_connection, \
    mysql_config_from_args
from planetgen.util import log


def main(argv=None):
    parser = argparse.ArgumentParser(description="Recover the accounts: make an account the Owner.")
    parser.add_argument("--set-owner", metavar="USER", required=True, help="The username that becomes the Owner.")
    add_mysql_connection_args(parser)
    log.add_debug_arguments(parser)
    args = parser.parse_args(argv)
    log.apply_debug_arguments(args, parser)
    try:
        conn = get_control_connection(control_mysql_config(mysql_config_from_args(args)))
    except pymysql.MySQLError as exc:
        print(f"error: could not open the control database ({exc}).", file=sys.stderr)
        return 1
    try:
        try:
            before = auth.set_owner(conn, args.set_owner)
        except auth.AuthError as exc:
            parser.error(str(exc))
        activity_log.event("DB", "owner.set", user=activity_log.login_name(), target=f"user:{args.set_owner}")
        if before is None:
            print(f"{args.set_owner} is the Owner.")
        else:
            print(f"{args.set_owner} is the Owner now; {before} is an admin.")
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    sys.exit(main())
