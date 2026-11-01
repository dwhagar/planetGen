# planetgen/util/secret_key.py

"""
The key that signs the pages' forms (`secret_key` / `PLANETGEN_SECRET_KEY`).

An explicit `secret_key` in `config.json` or the environment always wins. With
none set, one is generated and kept in `secret_key` inside the web-owned state
folder (`settings.state_dir()`), mode 0600, so every worker and every restart
reads the same key. Several processes may start at once: the key is written to
a private temporary file and hard-linked into place, which only one of them
can do; the others read what the winner wrote.

The key itself is never printed, logged or put in an error: `ensure` and the
command line (`python -m planetgen.util.secret_key`, run by `update.sh`) say
only that a key was missing and one is being set.
"""

import os
import secrets
import sys

from planetgen.util.settings import get_settings, state_dir

FILE_NAME = "secret_key"
MISSING_NOTICE = "No secret_key is set (config.json or PLANETGEN_SECRET_KEY): setting one now and keeping it in {path} (mode 0600)."


def key_path():
    return os.path.join(state_dir(), FILE_NAME)


def _read(path):
    try:
        with open(path, "r", encoding="ascii") as f:
            return f.read().strip()
    except (OSError, UnicodeDecodeError):
        return ""


def _create(path):
    """Write a new key at `path` unless another process got there first; returns the key in force."""
    folder = os.path.dirname(path)
    os.makedirs(folder, mode=0o700, exist_ok=True)
    temp = os.path.join(folder, f".{FILE_NAME}.{os.getpid()}.{secrets.token_hex(4)}.tmp")
    fd = os.open(temp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        with os.fdopen(fd, "w", encoding="ascii") as f:
            f.write(secrets.token_hex(32) + "\n")
        try:
            os.link(temp, path)
        except FileExistsError:
            pass
    finally:
        os.unlink(temp)
    return _read(path)


def ensure(configured="", announce=None, path=None):
    """
    `(key, kind)`; kind is "configured", "existing", "created" or "ephemeral" (the state folder
    can't be written, so this process alone gets a random key). `announce(message)` is called
    only when a key is missing, with a message that never holds the key.
    """
    if configured:
        return configured, "configured"
    path = path or key_path()
    existing = _read(path)
    if existing:
        return existing, "existing"
    say = announce or (lambda message: None)
    say(MISSING_NOTICE.format(path=path))
    try:
        key = _create(path)
    except OSError as exc:
        say(f"The key could not be kept ({type(exc).__name__}): this process uses a random one, and forms left open across a restart must be resubmitted.")
        return secrets.token_hex(32), "ephemeral"
    return (key, "created") if key else (secrets.token_hex(32), "ephemeral")


def main():
    key, kind = ensure(get_settings().secret_key, announce=lambda message: print(message, file=sys.stderr))
    if kind == "created":
        print("The secret key is set.")
    elif kind == "ephemeral":
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
