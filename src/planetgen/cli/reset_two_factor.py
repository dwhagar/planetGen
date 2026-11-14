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

import sys

from planetgen.admin import activity_log, auth
from planetgen.cli import admin_account_command


def main(argv=None):
    return admin_account_command.run(argv, "Turn off this admin's two-factor sign-in, for a lost phone.", "cli/reset-two-factor", _reset)


def _reset(conn, user_id, username):
    was_on = auth.disable_totp(conn, user_id)
    activity_log.event("DB", "totp.reset", user=activity_log.login_name(), target=f"user:{username}")
    return f"Two-factor sign-in for {username} is off{'' if was_on else ' (it was not set up)'}."


if __name__ == "__main__":
    sys.exit(main())
