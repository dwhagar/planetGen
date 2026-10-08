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
is checked. Both live in Redis (`RedisStore`), the same server as
Flask-Limiter's request counts, so every worker process shares them and a
restart keeps them. `MemoryStore` holds the same state in memory: the
fallback while Redis can't be reached, and for tests.

The rules are pure functions on a small state dict (`fail`, `succeed`,
`retry_after`), so both stores apply exactly the same ones.
"""

import ipaddress
import json
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
"""float: A subject idle (unlocked, no failure) this long is forgotten (its
Redis key expires); by then an address's doubling level has halved seven
times."""

MAX_SUBJECT_LENGTH = 128
"""int: The longest subject (address or username) counted."""

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


class RedisStore:
    """The counters in Redis (SEC.30), next to Flask-Limiter's own request
    counts, so every worker process shares them and a restart keeps them.
    One JSON value per subject, written under WATCH so two workers
    counting the same subject can't lose a failure; an idle subject's key
    expires by itself after `PRUNE_AFTER_SECONDS`."""

    def __init__(self, client, prefix):
        self.client = client
        self.key_prefix = f"{prefix}:login:"

    @classmethod
    def from_url(cls, url, prefix):
        """A store on the Redis server at `url`, its keys under `prefix`
        (the app's `RATELIMIT_KEY_PREFIX`, so one setting keeps a test's
        or an install's counts apart)."""
        import redis
        return cls(redis.Redis.from_url(url, decode_responses=True), prefix)

    def _key(self, scope, subject):
        return f"{self.key_prefix}{scope}:{subject}"

    def get(self, scope, subject):
        raw = self.client.get(self._key(scope, subject))
        return json.loads(raw) if raw else None

    def update(self, scope, subject, change):
        """As `MemoryStore.update`, as one optimistic transaction."""
        key = self._key(scope, subject)

        def run(pipe):
            raw = pipe.get(key)
            new, result = change(json.loads(raw) if raw else new_state())
            pipe.multi()
            if new is None:
                pipe.delete(key)
            else:
                pipe.set(key, json.dumps(new), ex=int(PRUNE_AFTER_SECONDS))
            return result
        return self.client.transaction(run, key, value_from_callable=True)

    def _keys(self, scope=None, subject=None):
        if scope is not None and subject is not None:
            return [self._key(scope, subject)]
        pattern = f"{self.key_prefix}{scope}:*" if scope is not None else f"{self.key_prefix}*"
        return list(self.client.scan_iter(match=pattern, count=500))

    def locked(self, now):
        keys = self._keys()
        rows = []
        for key, raw in zip(keys, self.client.mget(keys) if keys else ()):
            state = json.loads(raw) if raw else None
            if state and state["locked_until"] > now:
                scope, _, subject = key[len(self.key_prefix):].partition(":")
                rows.append({"scope": scope, "subject": subject, **state})
        return sorted(rows, key=lambda row: -row["locked_until"])

    def lift(self, scope=None, subject=None):
        keys = self._keys(scope, subject)
        return self.client.delete(*keys) if keys else 0

    def clear(self):
        self.lift()


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
