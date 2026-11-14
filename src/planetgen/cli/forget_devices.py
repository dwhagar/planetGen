# planetgen.cli.forget_devices

"""
Revokes every trusted-device cookie of one admin (SEC.22), so that admin's
browsers are asked for the two-factor code again:

    python3 -m planetgen.cli.forget_devices admin

The devices are rows in the control database, which uses the same MySQL
settings as every other script (`config.json`, `PLANETGEN_MYSQL_*`,
`--mysql-*`). Written to the activity log as `DB devices.revoke`.
"""

import sys

from planetgen.admin import activity_log, auth
from planetgen.cli import admin_account_command


def main(argv=None):
    return admin_account_command.run(argv, "Revoke every trusted-device cookie of this admin.", "cli/forget-devices", _revoke)


def _revoke(conn, user_id, username):
    revoked = auth.revoke_devices(conn, user_id)
    activity_log.event("DB", "devices.revoke", user=activity_log.login_name(), target=f"user:{username}", revoked=revoked)
    return f"Revoked {revoked} trusted device{'' if revoked == 1 else 's'} of {username}."


if __name__ == "__main__":
    sys.exit(main())
