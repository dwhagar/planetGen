# planetgen/admin/throttle.py

"""
Login lockouts per client address (SEC.1) and per username (SEC.21)
===================================================================

Two counters guard every password check:

- **Per address** (`IP_POLICY`, Boss's numbers): 3 failed logins from
  one address lock it for 5 minutes; each further lockout doubles that,
  up to 1 day. A successful login clears the address's failure count but
  not its doubling level, which halves for every day without a lockout,
  so one right guess can't wipe it. An IPv6 address is counted by its
  /64 (one machine often holds a whole /64). Loopback addresses and the
  configured allowlist (`config.json`'s `"login_allowlist"`) are never
  locked, so the server's own admin can always get in.
- **Per username** (`USER_POLICY`, the numbers `api/loginbackoff.py` had
  since PR #131): 10 free failures, then 1 s doubling to 15 minutes; an
  hour without failures starts over; a successful login clears it.
  Unknown usernames are counted like real ones, so a lock never reveals
  which usernames exist.

While either is locked a login is refused with 429 before the password
is checked. Both live in the control database's `login_throttle` table
(`DbStore`), so every worker process shares them and a restart keeps
them. `MemoryStore` holds the same state in memory: the fallback while
that table can't be used (before `update.sh` has created it), and for
tests.

The rules are pure functions on a small state dict (`fail`, `succeed`,
`retry_after`), so both stores apply exactly the same ones.
"""

import ipaddress
import math
import threading
import time
from collections import namedtuple

SCOPE_IP = "ip"
SCOPE_USER = "user"
SCOPES = (SCOPE_IP, SCOPE_USER)

Policy = namedtuple("Policy", "free_failures base_lock_seconds max_lock_seconds forget_after_seconds "
                              "per_lockout decay_seconds")
"""A counter's rules. `per_lockout` False (the username rule): every
failure past the free ones doubles the lock. `per_lockout` True (the
address rule): reaching `free_failures` locks and starts the count over,
and each lockout doubles the next, the level halving every
`decay_seconds` without one."""

IP_POLICY = Policy(free_failures=3, base_lock_seconds=300.0, max_lock_seconds=86400.0,
                   forget_after_seconds=86400.0, per_lockout=True, decay_seconds=86400.0)
"""Policy: Per address: 3 failures, then 5 minutes doubling to 1 day."""

USER_POLICY = Policy(free_failures=10, base_lock_seconds=1.0, max_lock_seconds=900.0,
                     forget_after_seconds=3600.0, per_lockout=False, decay_seconds=None)
"""Policy: Per username: 10 free failures, then 1 s doubling to 15 minutes."""

POLICIES = {SCOPE_IP: IP_POLICY, SCOPE_USER: USER_POLICY}

PRUNE_AFTER_SECONDS = 7 * 86400.0
"""float: Rows idle (unlocked, no failure) this long are deleted; by then
an address's doubling level has halved seven times."""

MAX_SUBJECT_LENGTH = 128
"""int: `login_throttle.subject` is VARCHAR(128)."""

MAX_MEMORY_ENTRIES = 10000
"""int: `MemoryStore` keeps at most this many entries; past it the idle
and then the oldest unlocked ones go first."""


# ---------------------------------------------------------------------
# Subjects
# ---------------------------------------------------------------------

def normalize_username(username):
    """The key a username is counted under: stripped, case-folded and
    cut to `MAX_SUBJECT_LENGTH`, so `Admin` and `admin ` share a count."""
    return (username or "").strip().casefold()[:MAX_SUBJECT_LENGTH]


def parse_allowlist(entries):
    """
    `ipaddress` networks for the allowlist's entries (single addresses or
    CIDR ranges, as strings). Returns `(networks, bad_entries)`.
    """
    networks, bad = [], []
    if isinstance(entries, str):
        entries = [part for part in entries.replace(",", " ").split()]
    for entry in entries or ():
        try:
            networks.append(ipaddress.ip_network(str(entry).strip(), strict=False))
        except ValueError:
            bad.append(entry)
    return networks, bad


def _address(address):
    try:
        ip = ipaddress.ip_address(str(address).strip())
    except ValueError:
        return None
    if ip.version == 6 and ip.ipv4_mapped is not None:
        ip = ip.ipv4_mapped
    return ip


def ip_subject(address, allowlist=()):
    """
    The key an address is counted under, or `None` when it is never
    locked: not a valid address, loopback, or in `allowlist` (networks
    from `parse_allowlist`). IPv4 addresses (and IPv4-mapped IPv6 ones)
    count alone; IPv6 addresses count by their /64.
    """
    ip = _address(address)
    if ip is None or ip.is_loopback or any(ip in network for network in allowlist):
        return None
    if ip.version == 6:
        return str(ipaddress.ip_network(f"{ip}/64", strict=False))
    return str(ip)


def is_private_address(address):
    """Whether `address` is on a private network (10/8, 192.168/16,
    fc00::/7, ...): what a reverse proxy in front of the app usually has."""
    ip = _address(address)
    return ip is not None and ip.is_private and not ip.is_loopback


# ---------------------------------------------------------------------
# The rules
# ---------------------------------------------------------------------

def new_state():
    """A subject with no history."""
    return {"failures": 0, "level": 0, "locked_until": 0.0, "last_failure_at": 0.0, "last_lockout_at": 0.0}


def retry_after(state, now):
    """Whole seconds until a locked subject may try again (at least 1),
    or 0 when it isn't locked."""
    if state is None or state["locked_until"] <= now:
        return 0
    return max(1, math.ceil(state["locked_until"] - now))


def current_level(policy, state, now):
    """The doubling level after decay: halved for every `decay_seconds`
    since the last lockout."""
    level = state["level"]
    if policy.decay_seconds and state["last_lockout_at"] and level:
        halvings = int((now - state["last_lockout_at"]) // policy.decay_seconds)
        level >>= min(max(halvings, 0), 63)
    return level


def fail(policy, state, now):
    """
    Counts one failed login against `state` (changed in place).

    Returns:
        float: The lock this failure started, in seconds; 0 when none.
    """
    if state["locked_until"] <= now and state["last_failure_at"] \
            and now - state["last_failure_at"] >= policy.forget_after_seconds:
        state["failures"] = 0
    state["failures"] += 1
    state["last_failure_at"] = now
    if policy.per_lockout:
        if state["failures"] < policy.free_failures:
            return 0.0
        level = current_level(policy, state, now)
        seconds = min(policy.max_lock_seconds, policy.base_lock_seconds * 2 ** min(level, 32))
        state.update(failures=0, level=level + 1, locked_until=now + seconds, last_lockout_at=now)
        return seconds
    extra = state["failures"] - policy.free_failures
    if extra <= 0:
        return 0.0
    seconds = min(policy.max_lock_seconds, policy.base_lock_seconds * 2 ** min(extra - 1, 32))
    state["locked_until"] = now + seconds
    state["last_lockout_at"] = now
    return seconds


def succeed(policy, state, now):
    """
    A successful login: returns the state to keep, or `None` to forget the
    subject. A username starts over; an address keeps its (decayed)
    doubling level.
    """
    if not policy.per_lockout:
        return None
    if not current_level(policy, state, now):
        return None
    state["failures"] = 0
    return state


def is_idle(state, now):
    """Unlocked, and no failure for `PRUNE_AFTER_SECONDS`."""
    return state["locked_until"] <= now and now - max(state["last_failure_at"], state["last_lockout_at"]) \
        >= PRUNE_AFTER_SECONDS


# ---------------------------------------------------------------------
# Stores
# ---------------------------------------------------------------------

_COLUMNS = ("failures", "level", "locked_until", "last_failure_at", "last_lockout_at")


class MemoryStore:
    """The counters in this process's memory. Thread-safe."""

    def __init__(self):
        self._lock = threading.Lock()
        self._entries = {}

    def get(self, scope, subject):
        with self._lock:
            state = self._entries.get((scope, subject))
            return dict(state) if state is not None else None

    def update(self, scope, subject, change):
        """Runs `change(state) -> (new state or None, result)` atomically
        and returns `result`."""
        with self._lock:
            state = dict(self._entries.get((scope, subject)) or new_state())
            new, result = change(state)
            if new is None:
                self._entries.pop((scope, subject), None)
            else:
                self._entries[(scope, subject)] = new
                self._prune()
            return result

    def locked(self, now):
        with self._lock:
            return [{"scope": scope, "subject": subject, **state}
                    for (scope, subject), state in self._entries.items() if state["locked_until"] > now]

    def lift(self, scope=None, subject=None):
        with self._lock:
            keys = [key for key in self._entries
                    if (scope is None or key[0] == scope) and (subject is None or key[1] == subject)]
            for key in keys:
                del self._entries[key]
            return len(keys)

    def clear(self):
        with self._lock:
            self._entries.clear()

    def __len__(self):
        return len(self._entries)

    def _prune(self):
        if len(self._entries) <= MAX_MEMORY_ENTRIES:
            return
        now = time.time()
        for key in [k for k, s in self._entries.items() if is_idle(s, now)]:
            del self._entries[key]
        overflow = len(self._entries) - MAX_MEMORY_ENTRIES
        if overflow > 0:
            # Unlocked subjects go first, oldest failure first, so a flood
            # of new names or addresses can't free a locked one.
            oldest = sorted(self._entries, key=lambda k: (self._entries[k]["locked_until"] > now,
                                                          self._entries[k]["last_failure_at"]))
            for key in oldest[:overflow]:
                del self._entries[key]


class DbStore:
    """The counters in the control database's `login_throttle` table,
    through an open control-schema `Connection`."""

    _last_prune = 0.0

    def __init__(self, conn):
        self.conn = conn

    def get(self, scope, subject):
        row = self.conn.execute(
            "SELECT failures, level, locked_until, last_failure_at, last_lockout_at FROM login_throttle "
            "WHERE scope = ? AND subject = ?", (scope, subject),
        ).fetchone()
        return _state_from_row(row) if row is not None else None

    def update(self, scope, subject, change):
        """As `MemoryStore.update`, in one transaction holding the row's
        lock, so two workers counting the same subject can't lose a
        failure."""
        conn = self.conn
        try:
            conn.execute("INSERT IGNORE INTO login_throttle (scope, subject) VALUES (?, ?)", (scope, subject))
            row = conn.execute(
                "SELECT failures, level, locked_until, last_failure_at, last_lockout_at FROM login_throttle "
                "WHERE scope = ? AND subject = ? FOR UPDATE", (scope, subject),
            ).fetchone()
            new, result = change(_state_from_row(row) if row is not None else new_state())
            if new is None:
                conn.execute("DELETE FROM login_throttle WHERE scope = ? AND subject = ?", (scope, subject))
            else:
                conn.execute(
                    "UPDATE login_throttle SET failures = ?, level = ?, locked_until = ?, last_failure_at = ?, "
                    "last_lockout_at = ? WHERE scope = ? AND subject = ?",
                    tuple(new[c] for c in _COLUMNS) + (scope, subject),
                )
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        self._maybe_prune()
        return result

    def locked(self, now):
        rows = self.conn.execute(
            "SELECT scope, subject, failures, level, locked_until, last_failure_at, last_lockout_at "
            "FROM login_throttle WHERE locked_until > ? ORDER BY locked_until DESC", (now,),
        ).fetchall()
        return [{"scope": row["scope"], "subject": row["subject"], **_state_from_row(row)} for row in rows]

    def lift(self, scope=None, subject=None):
        where, params = [], []
        if scope is not None:
            where.append("scope = ?")
            params.append(scope)
        if subject is not None:
            where.append("subject = ?")
            params.append(subject)
        sql = "DELETE FROM login_throttle" + (" WHERE " + " AND ".join(where) if where else "")
        cur = self.conn.execute(sql, tuple(params))
        self.conn.commit()
        return cur.rowcount

    def _maybe_prune(self):
        """Deletes idle rows, at most once an hour per process."""
        now = time.time()
        if now - DbStore._last_prune < 3600:
            return
        DbStore._last_prune = now
        cutoff = now - PRUNE_AFTER_SECONDS
        self.conn.execute(
            "DELETE FROM login_throttle WHERE locked_until <= ? AND last_failure_at < ? AND last_lockout_at < ?",
            (now, cutoff, cutoff),
        )
        self.conn.commit()


def _state_from_row(row):
    return {
        "failures": int(row["failures"]), "level": int(row["level"]),
        "locked_until": float(row["locked_until"]), "last_failure_at": float(row["last_failure_at"]),
        "last_lockout_at": float(row["last_lockout_at"]),
    }


# ---------------------------------------------------------------------
# Operations on a store
# ---------------------------------------------------------------------

def check(store, scope, subject, now=None):
    """Seconds `subject` must still wait (0 when it may try)."""
    now = time.time() if now is None else now
    return retry_after(store.get(scope, subject), now)


def record_failure(store, scope, subject, now=None):
    """Counts a failed login; returns the lock it started (seconds, 0 for none)."""
    now = time.time() if now is None else now
    policy = POLICIES[scope]

    def change(state):
        seconds = fail(policy, state, now)
        return state, seconds
    return store.update(scope, subject, change)


def record_success(store, scope, subject, now=None):
    """A successful login for `subject`."""
    now = time.time() if now is None else now
    policy = POLICIES[scope]
    if store.get(scope, subject) is None:
        return
    store.update(scope, subject, lambda state: (succeed(policy, state, now), None))


def locked_subjects(store, now=None):
    """Every subject locked right now, with `retry_after` filled in."""
    now = time.time() if now is None else now
    return [{**row, "retry_after": retry_after(row, now)} for row in store.locked(now)]
