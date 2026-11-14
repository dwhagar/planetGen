# planetgen.cli.admin_account_command

"""
The shape `forget_devices` and `reset_two_factor` share: a command line with one admin's username and the MySQL
connection and logging options, a control-database connection, the admin's row, an action on it written to the activity
log and said on the console.
"""

import argparse

import pymysql

from planetgen.admin import activity_log
from planetgen.db.store import add_mysql_connection_args, control_mysql_config, get_control_connection, \
    mysql_config_from_args
from planetgen.util import log


def parse_and_connect(parser, argv):
    """
    Adds the MySQL connection and debug options to `parser`, parses `argv` and opens the control database. Returns
    `(args, conn)`; `conn` is `None`, with the error logged, when the database cannot be opened (the command returns 1).
    """
    add_mysql_connection_args(parser)
    log.add_debug_arguments(parser)
    args = parser.parse_args(argv)
    log.apply_debug_arguments(args, parser, errors_to_stderr=True)
    try:
        return args, get_control_connection(control_mysql_config(mysql_config_from_args(args)))
    except pymysql.MySQLError as exc:
        log.error(f"The control database could not be opened ({exc}).")
        return args, None


def run(argv, description, console_name, act):
    """
    Parses `argv` (`USER` plus the connection and debug options), opens the control database (1 and an error line when it
    cannot be opened), finds the admin (a usage error for an unknown name) and calls `act(conn, user_id, username)`, which
    does the work, writes its own `activity_log.event` and returns the sentence to print; the call is written as
    `console_name` with the user as target. Returns 0.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("user", metavar="USER", help="The admin's username.")
    args, conn = parse_and_connect(parser, argv)
    if conn is None:
        return 1
    try:
        row = conn.execute("SELECT id FROM admin_users WHERE username = ?", (args.user,)).fetchone()
        if row is None:
            parser.error(f"no admin named {args.user!r}.")
        said = act(conn, row["id"], args.user)
        activity_log.console_call(console_name, target=f"user:{args.user}")
        print(said)
        return 0
    finally:
        conn.close()
