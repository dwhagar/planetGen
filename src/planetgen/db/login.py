# planetgen/db/login.py

"""
Logging in to MySQL (OPS.49)
============================

`connect(...)` is `pymysql.connect(...)` that explains a failed login: the
account's password is wrong, or its authentication plugin can't finish with
this pymysql. Every place that opens a database connection calls it, so
`update.sh` and the commands print
"Couldn't log in to the database server at 127.0.0.1:3306 as 'planetgen' ..."
and not a bare "'NoneType' object has no attribute 'is_auth_switch_request'"
(pymysql 1.2.1 to 1.2.3 raise that for a caching_sha2_password login without
TLS).

**The caching_sha2_password step.** MySQL 8's default accounts use that
plugin. pymysql 1.2.3 asks for TLS when the server offers it, and then logs
in. When the server has no TLS (`have_ssl` DISABLED, or `skip-ssl`) and the
password isn't cached yet (the first login after a server start, or after
`FLUSH PRIVILEGES`), pymysql sends the RSA-encrypted password and reads the
server's OK, but its `caching_sha2_password_auth` ends without `return`ing
that packet, and the login loop that 1.2.1 added then calls
`.is_auth_switch_request()` on `None`. PyMySQL merged the missing `return`
(PR #1275, 2026-09-18) after 1.2.3 and has not released it. `connect` hands
pymysql's own handler through its documented `auth_plugin_map` so the loop
takes the packet or the end of the exchange (`None`) from it: the login is
still pymysql's, with the same checks. The hook stays in place after the fix
is released (Boss, 2026-10-15): pymysql calls the handler in place of its own
built-in call to the same step, so the step still runs exactly once, and a
release that returns the packet passes it on unchanged. The requirement
`pymysql>=1.2.3` takes either. tests/test_db_login_server.py logs in to a
fake MySQL 8 server with both versions of the step.
"""

import pymysql

LOGIN_ERROR_CODES = (1045, 1251, 1698, 2059, 2061)
"""tuple: The MySQL and client error codes of a refused or unfinished login (access denied; the client doesn't support
the authentication protocol; access denied by plugin; authentication plugin can't be loaded; SSL required)."""


class DatabaseLoginError(pymysql.err.OperationalError):
    """The login to the database server failed (`connect`): the message names the server, the account and what to
    check. An `OperationalError` like the refusal it explains (code 1045, or 2059 when the handshake broke), so code
    that handles a database that can't be reached handles this too."""

    def __str__(self):
        return self.args[1]


class _CachingSha2Login:
    """pymysql's own caching_sha2_password step, registered through `auth_plugin_map` (see the module docstring)."""

    def __init__(self, connection):
        self.connection = connection

    def authenticate(self, auth_packet):
        from pymysql import _auth

        return _auth.caching_sha2_password_auth(self.connection, auth_packet)


def connect(**kwargs):
    """`pymysql.connect(**kwargs)`, raising `DatabaseLoginError` when the login itself fails."""
    kwargs.setdefault("auth_plugin_map", {"caching_sha2_password": _CachingSha2Login})
    try:
        return pymysql.connect(**kwargs)
    except pymysql.err.OperationalError as exc:
        if not exc.args or exc.args[0] not in LOGIN_ERROR_CODES:
            raise
        code, problem = exc.args[0], f"{exc.args[1] if len(exc.args) > 1 else exc.args[0]}"
    except (OSError, pymysql.err.Error):
        raise
    except Exception as exc:  # noqa: BLE001 -- pymysql broke inside the login handshake itself
        code, problem = 2059, f"{type(exc).__name__}: {exc}"
    raise DatabaseLoginError(
        code,
        f"Couldn't log in to the database server at {kwargs.get('host')}:{kwargs.get('port')} as "
        f"'{kwargs.get('user')}' ({problem}). This is the login step, before any query: check mysql.user and "
        f"mysql.password in config.json, that the account exists for this host, and that its authentication plugin "
        f"works with pymysql {pymysql.VERSION_STRING} (a MySQL 8 caching_sha2_password account needs the "
        f"cryptography package)."
    )
