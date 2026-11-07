#!/usr/bin/env python3
"""
Warns about Redis server settings that don't suit planetGen (OPS.48), for
scripts/deploy-common.sh (`check_redis_settings`, run by install.sh and
update.sh).

    check_redis_settings.py REDIS_URL

Reads the running server's settings with CONFIG GET and prints one warning
per setting that is off, each with the exact line to add to the Redis
config (examples/redis/planetgen.conf has them all). Never changes a
setting or a file, and always exits 0: this only warns. Needs the `redis`
package from planetGen's virtual environment.
"""
import sys

WANTED = {
    "appendonly": ("yes", "the append-only file is off, so a crash or restart loses every queued job and lockout",
                   "appendonly yes\nappendfsync everysec"),
    "maxmemory-policy": ("noeviction", "Redis may throw away queue entries or lockouts when memory fills",
                         "maxmemory-policy noeviction"),
}
LOOPBACK = {"127.0.0.1", "::1", "localhost"}


def listens_beyond_loopback(bind):
    """Whether a `bind` value (space-separated addresses, `-` marking optional ones) reaches past this machine."""
    addresses = [a.lstrip("-") for a in bind.split()]
    return not addresses or any(a not in LOOPBACK and not a.startswith("127.") for a in addresses)

ACL_LINES = ("user default off\n"
             "user planetgen on >A-LONG-RANDOM-PASSWORD ~* &* +@all -@dangerous +config|get +info +client|list\n"
             "and in config.json: \"redis\": {\"url\": \"redis://planetgen:A-LONG-RANDOM-PASSWORD@127.0.0.1:6379/0\"}")
"""str: What the warning about a Redis that anyone on the machine may use prints (`examples/redis/planetgen-acl.conf` has it)."""


def open_access_problem(whoami, default_flags):
    """
    `(what is wrong, lines to add)` when this connection is the server's passwordless default user (SEC.39): any local user may
    then queue a job that a planetGen worker runs. `None` when the connection logged in as a user of its own, or the default user
    is off or has a password.
    """
    flags = {str(flag).lower() for flag in (default_flags or ())}
    if str(whoami).lower() == "default" and "nopass" in flags and "on" in flags:
        return ("any user on this machine can use Redis without a password, and a planetGen worker runs the function a queued job "
                "names", ACL_LINES)
    return None


MIN_MAXMEMORY = 2 * 1024 ** 3
"""int: Below this, the Redis memory limit is called too small."""


def problems(config):
    """`[(what is wrong, line to add)]` for a dict of CONFIG GET values."""
    found = []
    for key, (wanted, why, line) in WANTED.items():
        if str(config.get(key, "")).lower() != wanted:
            found.append((f"{key} is {config.get(key)!r}, not {wanted!r}: {why}", line))
    if config.get("appendonly", "").lower() == "yes" and config.get("appendfsync", "everysec").lower() == "no":
        found.append(("appendfsync is 'no': the operating system decides when data reaches the disk",
                      "appendfsync everysec"))
    if str(config.get("maxmemory", "0")) in ("0", ""):
        found.append(("maxmemory is 0 (no limit): a runaway queue could fill the host's memory",
                      "maxmemory 4gb"))
    elif int(config.get("maxmemory") or 0) < MIN_MAXMEMORY:
        found.append((f"maxmemory is {config.get('maxmemory')} bytes: the page cache and the run records share it, "
                      "and a full Redis stops a galaxy run", "maxmemory 4gb"))
    if listens_beyond_loopback(config.get("bind", "")):
        found.append((f"bind is {config.get('bind', '')!r}, so it listens beyond this machine; only planetGen should "
                      "reach it", "bind 127.0.0.1 -::1\nprotected-mode yes"))
    return found


def main(argv):
    import redis

    try:
        client = redis.Redis.from_url(argv[1], socket_connect_timeout=3, socket_timeout=5, decode_responses=True)
        config = {}
        for pattern in ("appendonly", "appendfsync", "maxmemory", "maxmemory-policy", "bind"):
            config.update(client.config_get(pattern))
    except Exception as exc:  # noqa: BLE001 -- warn only
        print(f"warning: couldn't read the Redis settings ({type(exc).__name__}); not checked "
              "(CONFIG may be disabled).", file=sys.stderr)
        return 0
    found = problems(config)
    try:
        opened = open_access_problem(client.execute_command("ACL", "WHOAMI"), client.acl_getuser("default").get("flags"))
    except Exception:  # noqa: BLE001 -- Redis before 6, or a user not allowed ACL: nothing to say
        opened = None
    if opened:
        found.append(opened)
    if not found:
        print("Redis settings: append-only file on, noeviction, memory limit set, loopback only, a user of its own")
        return 0
    print(f"warning: the Redis server at {argv[1]} has settings planetGen doesn't want:", file=sys.stderr)
    for what, line in found:
        print(f"  - {what}.", file=sys.stderr)
        print("    Add to your Redis config:", file=sys.stderr)
        for text in line.splitlines():
            print(f"      {text}", file=sys.stderr)
    print("  Then restart Redis. /etc/redis/redis.conf is the usual file; examples/redis/planetgen.conf has every "
          "line. update.sh never edits your Redis config.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
