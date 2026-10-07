# html/api/loginguard.py

"""
The checks around every password check (SEC.1, SEC.20, SEC.21): the
per-address lockout and the per-username backoff
(`planetgen/admin/throttle.py`), and the activity-log line and audit
row each refused sign-in writes.

A view builds one `LoginGuard` for the username being tried and the
request's client address, then:

    refused = guard.refusal()       # a 429 response while either is locked
    ... check the password ...
    guard.failed("login.failed")    # or guard.succeeded()

The counters live in the control database's `login_throttle` table. If
that can't be used (an install where `update.sh` hasn't created it yet,
or the database is down), the guard falls back to this process's memory
(`memory_store`) and says so once on stderr, rather than failing every
login or dropping the protection.

`LOGIN_BACKOFF_ENABLED = False` in the app config turns both counters off
(a test harness that fails many logins on purpose); the activity log
and audit rows are still written.
"""

import sys

from flask import current_app, jsonify, request

from planetgen.admin import activity_log, auth, throttle
from planetgen.util import log

from .common import get_control_db

memory_store = throttle.MemoryStore()
"""MemoryStore: The fallback while `login_throttle` can't be used (and
what tests clear between runs)."""

_warned = set()


def _warn_once(key, message):
    if key not in _warned:
        _warned.add(key)
        print(message, file=sys.stderr)
        log.error(message)


def allowlist():
    """The configured `login_allowlist` as networks; a bad entry is
    skipped with one warning."""
    networks, bad = throttle.parse_allowlist(current_app.config.get("LOGIN_ALLOWLIST") or ())
    for entry in bad:
        _warn_once(f"allowlist:{entry}", f"planetgen: ignoring login_allowlist entry {entry!r}: "
                                         f"not an address or network.")
    return networks


def with_store(operation):
    """Runs `operation(store)` against the database, or against
    `memory_store` when the database can't be used."""
    try:
        return operation(throttle.DbStore(get_control_db()))
    except Exception as exc:  # noqa: BLE001 -- never fail a login because of the counters
        _warn_once("store", f"planetgen: the login_throttle table can't be used ({exc}); counting failed "
                            f"logins in memory per process until it can. Run update.sh to create it.")
        return operation(memory_store)


def wait_text(seconds):
    """`seconds` as words for a message: `1 second`, `40 seconds`,
    `5 minutes`, `3 hours`."""
    seconds = max(1, int(seconds))
    if seconds < 120:
        return f"{seconds} second{'' if seconds == 1 else 's'}"
    if seconds < 2 * 3600:
        minutes = -(-seconds // 60)
        return f"{minutes} minutes"
    hours = -(-seconds // 3600)
    return f"{hours} hours"


def record_failure_row(action, username, admin_user_id=None):
    """One `admin_audit_log` row for a refused sign-in; a database error is
    logged, never raised."""
    try:
        auth.record_login_failure(get_control_db(), action, username,
                                       ip=activity_log.clean_ip(request.remote_addr), admin_user_id=admin_user_id)
    except Exception as exc:  # noqa: BLE001
        log.error(f"Could not record {action} for {username!r} in admin_audit_log: {exc}")


class LoginGuard:
    """The lockout checks for one sign-in attempt by `username` from this
    request's address."""

    def __init__(self, username, admin_user_id=None, trusted_device=False):
        self.username = username
        self.admin_user_id = admin_user_id
        self.enabled = current_app.config.get("LOGIN_BACKOFF_ENABLED", True)
        self.address = request.remote_addr
        self.ip = throttle.ip_subject(self.address, allowlist()) if self.enabled else None
        # A browser holding this admin's device cookie (SEC.22) is neither
        # held back nor counted by the per-username lock, so failures on
        # purpose from elsewhere can't lock the real admin out; the
        # per-address lock still applies to it.
        self.user = throttle.normalize_username(username) if self.enabled and not trusted_device else None

    def refusal(self):
        """
        A 429 response while this address or username is locked (checked
        before the password, so a guess made during a lock learns
        nothing), after logging it as `login.locked`; `None` otherwise.
        """
        for scope, subject, what in ((throttle.SCOPE_IP, self.ip, "from this address"),
                                     (throttle.SCOPE_USER, self.user, "for this username")):
            if not subject:
                continue
            wait = with_store(lambda store: throttle.check(store, scope, subject))
            if wait:
                activity_log.event("AUTH", "login.locked", user=self.username, scope=scope, retry_after=wait)
                record_failure_row("login.locked", self.username, self.admin_user_id)
                resp = jsonify({"error": f"too many failed logins {what}; try again in {wait_text(wait)}",
                                "retry_after": wait, "scope": scope})
                resp.status_code = 429
                resp.headers["Retry-After"] = str(wait)
                return resp
        return None

    def failed(self, action="login.failed"):
        """Counts a wrong password against the address and the username,
        and logs it (plus a `lockout.start` line for each lock it starts)."""
        locks = {}
        for scope, subject in ((throttle.SCOPE_IP, self.ip), (throttle.SCOPE_USER, self.user)):
            if subject:
                locks[scope] = with_store(lambda store: throttle.record_failure(store, scope, subject))
        activity_log.event("AUTH", action, user=self.username)
        record_failure_row(action, self.username, self.admin_user_id)
        for scope, seconds in locks.items():
            if seconds:
                activity_log.event("AUTH", "lockout.start", user=self.username, scope=scope,
                                  subject=self.ip if scope == throttle.SCOPE_IP else self.user,
                                  seconds=int(seconds))
        if locks.get(throttle.SCOPE_IP) and throttle.is_private_address(self.address) \
                and not (current_app.config.get("PROXY_FIX") or {}).get("x_for"):
            _warn_once("proxy", f"planetgen: locked out the private address {self.address} after failed logins. "
                                f"If the site is behind a reverse proxy, every visitor has that address: set "
                                f"proxy_fix.x_for in config.json (docs/config.md), or add the proxy to "
                                f"login_allowlist.")

    def succeeded(self):
        """Clears the username's count and the address's failure count
        (its doubling level stays and decays)."""
        for scope, subject in ((throttle.SCOPE_IP, self.ip), (throttle.SCOPE_USER, self.user)):
            if subject:
                with_store(lambda store: throttle.record_success(store, scope, subject))
