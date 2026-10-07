# planetgen/admin/activity_log.py

"""
The activity log: an always-on record of who did what (SEC.28)
================================================================

Separate from the debug log (`planetgen.util.log`), which is written only
while `config.json`'s `"debug"` is on and records every SQL statement and
random roll, this file is written all the time, at one line per event,
and holds only what an administrator (or a tool such as fail2ban) needs
to see:

- sign-ins: successful, failed and locked logins, logouts, credential
  changes and wrong current passwords (`AUTH`);
- refused requests: a missing, expired or unknown session or API key, an
  admin route refused, a failed CSRF check (`AUTHZ`);
- database changes: every write the web interface or API makes, every
  schema migration, and each generation run or job as one start and one
  finish line (`DB`, `GEN`).

Where it goes (`appconfig.activity_log_path`): `/var/log/planetgen/
planetgen.log` on Linux, `/Library/Logs/planetgen/planetgen.log` on
macOS, and `logs\\planetgen.log` under the checkout on Windows;
`config.json`'s `"log_dir"` or `PLANETGEN_LOG_DIR` moves it.

Every line has the same shape, documented in `docs/config.md`, so tools
can match it::

    2026-10-01T08:00:00Z planetgen[1234]: AUTH login.failed ip=203.0.113.5 user="admin"

the time in UTC, the process id, the category, the event, then `ip=` (the
client address, or `-` when there is none, such as a command-line run)
and `user=` (always quoted), then any further `key=value` fields. The
address always comes before any text a visitor typed, and that text is
always quoted with `"` and `\\` escaped and control characters written
as `\\xNN`, so a crafted username can't start a new line or fake a
field. Passwords, tokens and keys are never passed in.

Rotation (`appconfig.log_rotation_mode`): where install.sh/update.sh
installed `/etc/logrotate.d/planetgen-log` (Linux) or
`/etc/newsyslog.d/planetgen-log.conf` (macOS), the file is opened with a
`WatchedFileHandler`, which reopens it after the system tool moves it.
Otherwise (Windows, or a setup without root) the program rotates it
itself: 100 MB per file, 30 old copies.

Writing to it never stops the program: if the file can't be opened, one
warning goes to stderr (only when its folder exists, so a development
checkout without the folder stays quiet) and later events are dropped.
With the debug log on, every line is copied there too.
"""

import ipaddress
import logging
import logging.handlers
import os
import re
import sys
import threading
import time

from planetgen.util import appconfig

ROTATE_MAX_BYTES = 100 * 1024 * 1024
"""int: The size at which the program rotates the file itself (when no
system rotation is installed)."""

ROTATE_BACKUP_COUNT = 30
"""int: Old copies kept by the program's own rotation."""

FILE_MODE = 0o660
"""int: Mode for a file this module creates: the web server's user and
group can write it (a command-line user in that group can append), other
users can't read it."""

MAX_VALUE_LENGTH = 200
"""int: Quoted values (usernames, paths, details) are cut to this many
characters, so one request can't write an arbitrarily long line."""

CATEGORIES = ("AUTH", "AUTHZ", "DB", "GEN")
"""tuple: The categories a line may carry."""

_logger = logging.getLogger("planetgen.activity")
_logger.propagate = False
_logger.setLevel(logging.INFO)

_lock = threading.Lock()
_handler = None
_handler_path = None
_disabled_path = None
_warned = set()

_PLAIN_VALUE = re.compile(r"^[A-Za-z0-9_.:/@+-]{1,200}$")
_ACTION = re.compile(r"^[a-z][a-z0-9_.-]{0,63}$")
_KEY = re.compile(r"^[a-z][a-z0-9_]{0,31}$")


class _CreateModeMixin:
    """Opens the log file with `FILE_MODE` when it has to be created, rather
    than the process umask's usual world-readable 0644."""

    def _open(self):
        fd = os.open(self.baseFilename, os.O_WRONLY | os.O_APPEND | os.O_CREAT, FILE_MODE)
        return os.fdopen(fd, "a", encoding=self.encoding or "utf-8")


class _WatchedHandler(_CreateModeMixin, logging.handlers.WatchedFileHandler):
    pass


class _RotatingHandler(_CreateModeMixin, logging.handlers.RotatingFileHandler):
    pass


def quote(value):
    """
    `value` as a quoted field: `"..."`, with `"` and `\\` escaped, every
    control character (and the Unicode line and paragraph separators)
    written as `\\xNN`/`\\uNNNN`, and cut to `MAX_VALUE_LENGTH` characters.
    """
    text = "" if value is None else str(value)
    if len(text) > MAX_VALUE_LENGTH:
        text = text[:MAX_VALUE_LENGTH] + "..."
    out = []
    for char in text:
        code = ord(char)
        if char in ('"', "\\"):
            out.append("\\" + char)
        elif code < 0x20 or code == 0x7F:
            out.append(f"\\x{code:02x}")
        elif 0x80 <= code < 0xA0 or char in (" ", " "):
            out.append(f"\\u{code:04x}")
        else:
            out.append(char)
    return '"' + "".join(out) + '"'


def _value(value):
    """A field value: bare when it is a number or a plain token, else quoted."""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.3f}".rstrip("0").rstrip(".")
    text = "" if value is None else str(value)
    return text if _PLAIN_VALUE.match(text) else quote(text)


def clean_ip(address):
    """
    `address` when it is a real IPv4 or IPv6 address, else `-`. The
    address can come from an `X-Forwarded-For` header (`proxy_fix`), so it
    is checked before it goes into the line unquoted.
    """
    if not address:
        return "-"
    try:
        return str(ipaddress.ip_address(str(address).strip()))
    except ValueError:
        return "-"


def format_line(category, action, ip=None, user=None, fields=None, when=None, pid=None):
    """
    One activity-log line, without the trailing newline.

    Args:
        category (str): One of `CATEGORIES`.
        action (str): The event, lower case with dots (`login.failed`).
        ip (str, optional): The client address; `-` when missing or not
            a valid address.
        user (str, optional): Who: the admin, or the username as typed;
            `-` when nobody.
        fields (dict, optional): Further `key=value` pairs, in order.
        when (float, optional): Seconds since the epoch (default now).
        pid (int, optional): Process id (default this process).
    """
    if category not in CATEGORIES:
        raise ValueError(f"unknown activity log category {category!r}")
    if not _ACTION.match(action):
        raise ValueError(f"bad activity log action {action!r}")
    stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(time.time() if when is None else when))
    parts = [f"{stamp} planetgen[{os.getpid() if pid is None else pid}]: {category} {action}",
             f"ip={clean_ip(ip)}", f"user={quote('-' if user is None else user)}"]
    for key, value in (fields or {}).items():
        if not _KEY.match(key) or key in ("ip", "user"):
            raise ValueError(f"bad activity log field name {key!r}")
        if value is None:
            continue
        parts.append(f"{key}={_value(value)}")
    return " ".join(parts)


def _current_handler():
    """The open file handler, opening it on first use; `None` when the
    file can't be opened (warned about once)."""
    global _handler, _handler_path, _disabled_path
    if _handler is not None:
        return _handler
    try:
        config = appconfig.load_config()
        path = appconfig.activity_log_path(config)
        mode = appconfig.log_rotation_mode(config)
    except Exception as exc:  # noqa: BLE001 -- logging must never stop the program
        _warn_once("config", f"planetgen: activity log disabled, could not read config.json: {exc}")
        return None
    if path == _disabled_path:
        return None
    folder = os.path.dirname(path)
    try:
        if not os.path.isdir(folder):
            os.makedirs(folder, mode=0o770, exist_ok=True)
        if mode == "app":
            handler = _RotatingHandler(path, maxBytes=ROTATE_MAX_BYTES, backupCount=ROTATE_BACKUP_COUNT,
                                       encoding="utf-8")
        else:
            handler = _WatchedHandler(path, encoding="utf-8")
    except OSError as exc:
        _disabled_path = path
        if os.path.isdir(folder):
            _warn_once(path, f"planetgen: the activity log {path} can't be opened ({exc}); run "
                             f"install.sh/update.sh as root to set it up, or set \"log_dir\" in config.json.")
        return None
    handler.setFormatter(logging.Formatter("%(message)s"))
    _logger.addHandler(handler)
    _handler, _handler_path = handler, path
    return handler


def _warn_once(key, message):
    if key not in _warned:
        _warned.add(key)
        print(message, file=sys.stderr)


def _request_ip():
    """The current Flask request's client address, if there is one."""
    try:
        from flask import has_request_context, request
    except ImportError:
        return None
    if not has_request_context():
        return None
    return request.remote_addr


def event(category, action, user=None, ip=None, **fields):
    """
    Writes one line (see `format_line`). Inside a Flask request the
    client address is taken from it unless `ip` is given. Never raises.

    Returns:
        str | None: The line written (also when only the debug log got
            it), or `None` if it couldn't be formatted.
    """
    try:
        line = format_line(category, action, ip=ip if ip is not None else _request_ip(), user=user, fields=fields)
    except Exception as exc:  # noqa: BLE001
        _warn_once(f"format:{category}:{action}", f"planetgen: activity log line not written: {exc}")
        return None
    with _lock:
        handler = _current_handler()
    if handler is not None:
        try:
            _logger.info(line)
        except Exception:  # noqa: BLE001 -- logging's own handleError already reported it
            pass
    from planetgen.util import log
    log.trace("activity: %s", line)
    return line


def log_path():
    """The file this process writes to, or `None` while it isn't open."""
    return _handler_path


def reset():
    """Closes the file and forgets the path, so the next event re-reads
    the configuration (tests, and a process whose `log_dir` changed)."""
    global _handler, _handler_path, _disabled_path
    with _lock:
        if _handler is not None:
            _logger.removeHandler(_handler)
            _handler.close()
        _handler = _handler_path = _disabled_path = None
        _warned.clear()
