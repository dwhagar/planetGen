# planetgen/util/logrecord.py

"""
The logging layer (UX.121)
==========================

The one record every log line the program writes for a person or a tool to read is made from, and the three ways it is
written, so a file and a screen can never disagree. The standard, and every place the program logs today, are in
`docs/design/operations/logging-standard.md`.

A `Record` is a time (whole milliseconds since the epoch), a severity (the RFC 5424 levels `debug`, `info`, `notice`, `warning`,
`error`, `critical`), a short component name, a message and ordered `key=value` fields for numbers and names a tool
may want (a count, an elapsed time). It is written as:

* a **file line** (`file_line`), logfmt, one record per line, the time local with the offset of its zone (OPS.97)::

      ts=2026-10-24T14:40:00.123-0400 level=info component=generate msg="Stage 4 of 12: Scatter the massive stars" stage=4 stages=12

* a **JSON line** (`json_line`), the same record as JSON Lines, for machines;
* a **display line** (`display_line`), for the command line and the web log: the time in local time and labelled with its
  zone, the severity, the component, the message and the fields, with a stage's start, finish and failure drawn as one
  divider (`=== Stage 4 of 12: label ===`). `plain_line` is the command line's untimed text.

`parse_file_line` and `parse_json_line` read a line back to the record it came from, so what a file holds and what the
screen shows are the same records. This module is the record and its lines, and nothing that writes: `planetgen.util.log` is the one
place a record is made from a `logging` call and written (`log.get(component)` is a standard logger that takes fields as keyword
arguments, and `log.LineFormatter` writes what it logs in every form). Standard library only.
"""

import datetime
import json
import logging
import math
import re
import time

from planetgen.util import stagemark

SEVERITIES = ("debug", "info", "notice", "warning", "error", "critical")
"""tuple: The severities, least severe first: the RFC 5424 levels, named as the syslog standard names them."""

SYSLOG_CODES = {"debug": 7, "info": 6, "notice": 5, "warning": 4, "error": 3, "critical": 2}
"""dict: The RFC 5424 severity number of each severity (a lower number is more severe)."""

NOTICE = 25
"""int: The `logging` level of `notice`, between `info` (20) and `warning` (30), named in `logging` by `planetgen.util.log`."""

LOGGING_LEVELS = {"debug": logging.DEBUG, "info": logging.INFO, "notice": NOTICE, "warning": logging.WARNING,
                  "error": logging.ERROR, "critical": logging.CRITICAL}
"""dict: The `logging` level of each severity."""

RESERVED_FIELDS = ("ts", "level", "component", "msg")
"""tuple: The names the record itself uses in a file line and a JSON line; a field cannot take one."""

COMPONENTS = {
    "update": "`update.sh`, the update of an installed site",
    "install": "`install.sh`, the first install",
    "warm": "`planetgen warm`: the cache warm-up, its progress and its log",
    "tiles": "`planetgen tiles`: the nightly whole-map build",
    "orbits": "`planetgen orbits`: the nightly orbit update",
    "generate": "the generator: `planetgen generate` and its stages, the physics and the galaxy code a run calls",
    "queue": "the work queue: runs, units, workers and the machine's lease",
    "web": "the web app and the API: requests, pages, the admin pages",
    "job": "the Generate page's job runner, which runs a job's steps as processes",
    "activity": "the always-on activity log (`AUTH`, `AUTHZ`, `DB`, `GEN` events)",
    "command": "any other command line tool: `migrate`, `accounts`, `backup`, `config`, `query` and the rest",
    "lib": "another library's own log lines (Flask, urllib3, ...) that the debug log picks up",
    "logs": "`examples/apache/setup-debug-log.sh`: the log files and their rotation",
    "permissions": "`examples/apache/set-permissions.sh`: ownership and permissions",
    "cache-dir": "`examples/apache/create-cache-dir.sh`: the cache, jobs and state folders",
    "swap": "`scripts/enable-swap.sh`: the swap file",
    "python-deps": "`scripts/install-python-deps.sh`: the virtual environment and its libraries",
}
"""dict: The component names, each with what writes under it (OPS.103; docs/design/operations/logging-standard.md, "Components"). A record
takes its component from this table and from nowhere else; `Record` refuses a name that is not in it, `component_of` finds the one a logger or a
source file belongs to, and `src/tests/test_no_own_log_formatting.py` fails on a literal that is not here."""


COMPONENT_LOGS = {
    "warm": "warmup.log",
    "tiles": "tiles.log",
    "orbits": "orbits.log",
    "update": "update.log",
    "install": "console",
    "web": "apache",
    "activity": "planetgen.log",
}
"""dict: The log a component's lines are also written to (OPS.116): `error.log` names it on each error (`log=warmup.log`) so the line can be
traced back. A component not here is a command that prints to its terminal (`console`); the web app's stderr is the web server's (`apache`)."""


def origin_log(component):
    """The name of the log (or `console`, or `apache`) a `component`'s lines also go to; see `COMPONENT_LOGS`."""
    return COMPONENT_LOGS.get(component, "console")


FIELD_PATTERN = re.compile(r"^[a-z][a-z0-9_]{0,31}$")

COUNTING_WORDS = ("Stage", "Part", "Phase")
"""tuple: The three words that say where a script or a run is (OPS.103): `Stage N of M` for the top-level steps of a script or a run (update,
install, generate), `Part N of M` for the runs inside one (the warm-up's two), `Phase N of M` for a named group inside a run. No other word
(`step`, `pass`) counts, and a family of the warm-up is named, not numbered; `chunk` is only ever a count (`28 of 180 chunks`)."""

STAGE_START, STAGE_FINISH, STAGE_FAILED = "stage.start", "stage.finish", "stage.failed"
STAGE_SKIPPED = "stage.skipped"
"""str: The `event` field of the records that mark a stage (it starts, it finishes, it fails, it is skipped); the display
draws each as a divider."""

STAGE_EVENTS = (STAGE_START, STAGE_FINISH, STAGE_FAILED, STAGE_SKIPPED)

STDERR_SEVERITIES = ("warning", "error", "critical")
"""tuple: The severities a command writes to standard error (the rest go to standard output)."""

DIVIDER = "==="
"""str: What a stage divider is drawn with (`=== Stage 4 of 12: label ===`)."""

_NUMBER_LIKE = re.compile(r"^-?(?:\d+(?:\.\d*)?(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?)$")
_BARE_VALUE = re.compile(r'^[^\s"\\=\x00-\x1f\x7f\u0085\u2028\u2029]+$')
_EPOCH = datetime.datetime(1970, 1, 1, tzinfo=datetime.timezone.utc)

logging.addLevelName(NOTICE, "NOTICE")


class Record:
    """
    One log record. `ms` is the time in whole milliseconds since the Unix epoch, `severity` one of `SEVERITIES`,
    `component` a short lower-case name (`generate`, `queue`), `message` the text and `fields` an ordered dict of
    names to a string, an integer, a finite number or a boolean.

    Raises:
        ValueError: A severity, component or field name that is not allowed, a field named like one of the record's
            own (`RESERVED_FIELDS`), a field value that is not a string, an integer, a finite number or a boolean.
    """

    __slots__ = ("ms", "severity", "component", "message", "fields")

    def __init__(self, ms, severity, component, message, fields=None):
        if severity not in SYSLOG_CODES:
            raise ValueError(f"no severity {severity!r}; the severities are {', '.join(SEVERITIES)}")
        if component not in COMPONENTS:
            raise ValueError(f"component {component!r} is not in the component table ({', '.join(COMPONENTS)})")
        fields = dict(fields or {})
        for name, value in fields.items():
            if not FIELD_PATTERN.match(name) or name in RESERVED_FIELDS:
                raise ValueError(f"{name!r} is not a field name (lower-case a-z, 0-9, _; not {', '.join(RESERVED_FIELDS)})")
            if not isinstance(value, (str, int, float, bool)) or (isinstance(value, float) and not math.isfinite(value)):
                raise ValueError(f"field {name}={value!r}: a field holds a string, an integer, a finite number or a boolean")
        self.ms = int(ms)
        self.severity = severity
        self.component = component
        self.message = str(message)
        self.fields = fields

    @classmethod
    def at(cls, epoch, severity, component, message, /, **fields):
        """A record for the time `epoch` (seconds since the Unix epoch; rounded to the millisecond)."""
        return cls(round(epoch * 1000), severity, component, message, fields)

    @classmethod
    def now(cls, severity, component, message, /, **fields):
        return cls.at(time.time(), severity, component, message, **fields)

    @property
    def time(self):
        return self.ms / 1000

    @property
    def event(self):
        """The `event` field (a stage's start, finish or failure), or `None`."""
        return self.fields.get("event")

    def _key(self):
        return (self.ms, self.severity, self.component, self.message, tuple(self.fields.items()))

    def __eq__(self, other):
        return isinstance(other, Record) and self._key() == other._key()

    def __hash__(self):
        return hash(self._key())

    def __repr__(self):
        return f"Record({self.ms}, {self.severity!r}, {self.component!r}, {self.message!r}, {self.fields!r})"


def _stamp(ms):
    """`"2026-10-24T14:40:00.123-0400"` for `ms` milliseconds since the epoch: the computer's local time with the offset of its zone (OPS.97)."""
    moment = (_EPOCH + datetime.timedelta(milliseconds=ms)).astimezone()
    return moment.strftime("%Y-%m-%dT%H:%M:%S") + f".{moment.microsecond // 1000:03d}" + moment.strftime("%z")


def _ms_of(stamp):
    """The milliseconds since the epoch of a `_stamp` text.

    Raises:
        ValueError: Not a time in that form.
    """
    match = re.match(r"^(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d)\.(\d{3})([+-]\d{4})$", stamp)
    if not match:
        raise ValueError(f"{stamp!r} is not a time like 2026-10-24T14:40:00.123-0400")
    moment = datetime.datetime.strptime(match.group(1) + match.group(3), "%Y-%m-%dT%H:%M:%S%z")
    return int((moment - _EPOCH) / datetime.timedelta(milliseconds=1)) + int(match.group(2))


def local_text(seconds):
    """`"2026-10-24 14:40:00 EDT"` for `seconds` since the epoch: the local time labelled with its zone, as the display line writes it."""
    moment = (_EPOCH + datetime.timedelta(seconds=seconds)).astimezone()
    return f"{moment.strftime('%Y-%m-%d %H:%M:%S')} {_zone_name(moment)}"


_ESCAPES = {"\\": "\\\\", '"': '\\"', "\n": "\\n", "\r": "\\r", "\t": "\\t"}
_UNESCAPES = {"\\": "\\", '"': '"', "n": "\n", "r": "\r", "t": "\t"}


def _quote(text):
    out = []
    for char in text:
        if char in _ESCAPES:
            out.append(_ESCAPES[char])
        elif ord(char) < 0x20 or ord(char) == 0x7F:
            out.append(f"\\x{ord(char):02x}")
        elif char in "\u0085\u2028\u2029":
            out.append(f"\\u{ord(char):04x}")
        else:
            out.append(char)
    return '"' + "".join(out) + '"'


def _value_text(value):
    """A field value as the file line writes it: a number or a boolean bare, a string quoted when it needs it (or when
    it would read back as something else)."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        text = format(value, ".15g") if float(format(value, ".15g")) == value else repr(value)
        return text if re.search(r"[.eE]", text) else text + ".0"
    if _BARE_VALUE.match(value) and not _NUMBER_LIKE.match(value) and value not in ("true", "false"):
        return value
    return _quote(value)


def file_line(record):
    """The record as a file line (logfmt, the time local with its offset, no line break)."""
    parts = [f"ts={_stamp(record.ms)}", f"level={record.severity}", f"component={record.component}",
             f"msg={_quote(record.message)}"]
    parts.extend(f"{name}={_value_text(value)}" for name, value in record.fields.items())
    return " ".join(parts)


def json_line(record):
    """The record as a JSON Lines line: the same names as the file line, the fields beside them."""
    body = {"ts": _stamp(record.ms), "level": record.severity, "component": record.component, "msg": record.message}
    body.update(record.fields)
    text = json.dumps(body, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    return text.replace("\u0085", "\\u0085").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")


_PAIR = re.compile(r'([a-z][a-z0-9_]*)=("(?:[^"\\]|\\.)*"|[^\s"]*)(?:\s+|$)')


def _unquote(text):
    out, index = [], 1
    while index < len(text) - 1:
        char = text[index]
        if char != "\\":
            out.append(char)
            index += 1
            continue
        code = text[index + 1]
        if code in _UNESCAPES:
            out.append(_UNESCAPES[code])
            index += 2
        elif code == "x":
            out.append(chr(int(text[index + 2:index + 4], 16)))
            index += 4
        elif code == "u":
            out.append(chr(int(text[index + 2:index + 6], 16)))
            index += 6
        else:
            raise ValueError(f"unknown escape \\{code} in {text!r}")
    return "".join(out)


def _bare_value(text):
    if text == "true":
        return True
    if text == "false":
        return False
    if _NUMBER_LIKE.match(text):
        return float(text) if re.search(r"[.eE]", text) else int(text)
    return text


def parse_file_line(line):
    """The record a file line was written from.

    Raises:
        ValueError: The line is not a record in the form `file_line` writes.
    """
    line = line.rstrip("\n")
    pairs, position = [], 0
    while position < len(line):
        match = _PAIR.match(line, position)
        if not match:
            raise ValueError(f"not a log line at column {position + 1}: {line!r}")
        raw = match.group(2)
        pairs.append((match.group(1), _unquote(raw) if raw.startswith('"') else _bare_value(raw)))
        position = match.end()
    names = [name for name, _value in pairs[:4]]
    if names != list(RESERVED_FIELDS):
        raise ValueError(f"a log line starts with {', '.join(RESERVED_FIELDS)}: {line!r}")
    (_ts, stamp), (_level, severity), (_component, component), (_msg, message) = pairs[:4]
    return Record(_ms_of(stamp), severity, component, message, dict(pairs[4:]))


def parse_json_line(line):
    """The record a JSON line was written from.

    Raises:
        ValueError: The line is not a record in the form `json_line` writes.
    """
    body = json.loads(line)
    if not isinstance(body, dict) or any(name not in body for name in RESERVED_FIELDS):
        raise ValueError(f"not a log record: {line!r}")
    fields = {name: value for name, value in body.items() if name not in RESERVED_FIELDS}
    return Record(_ms_of(body["ts"]), body["level"], body["component"], body["msg"], fields)


def stage_start(first, total, label, count=1):
    """The message of a stage's start record: `"Stage 4 of 12: label"`, or `"Stages 4 to 6 of 12: label"` for a run of them."""
    return f"{stagemark.marker(first, total, count)}: {label}"


def stage_finish(first, total, seconds, count=1, peak_mb=0.0):
    """The message of a stage's finish record: `"Stage 4 of 12 finished after 31 s"`, and, when a unit reported its memory (PERF.226),
    `"Stage 4 of 12 finished after 31 s (largest worker 640 MB)"`."""
    memory = f" (largest worker {peak_mb:.0f} MB)" if peak_mb else ""
    return f"{stagemark.marker(first, total, count)} finished after {seconds:.0f} s{memory}"


def stage_skipped(first, total, label, reason, count=1):
    """The message of a skipped stage's record: `"Stage 4 of 12: label -- skipped: reason"`."""
    return f"{stagemark.marker(first, total, count)}: {label} -- skipped: {reason}"


def stage_failed(first, total, seconds, reason, count=1):
    """The message of a stage's failure record: `"Stage 4 of 12 failed after 31 s: reason"`."""
    return f"{stagemark.marker(first, total, count)} failed after {seconds:.0f} s: {reason}"


def stage_record(event, component, first, total, label=None, seconds=0.0, reason=None, count=1, now=None):
    """
    The record of a stage's start (`STAGE_START`, with its `label`), finish (`STAGE_FINISH`, after `seconds`), failure
    (`STAGE_FAILED`, after `seconds`, with the `reason`) or skip (`STAGE_SKIPPED`, with its `label` and the `reason`): the
    one wording of each, the severity `info`, `info`, `error` and `info`, and the fields `event`, `stage` (the first stage), `stages` (how many there are in all), `count` (when more
    than one) and `seconds` (finish and failure).

    Raises:
        ValueError: An `event` that is not one of `STAGE_EVENTS`.
    """
    if event == STAGE_START:
        severity, message = "info", stage_start(first, total, label, count)
    elif event == STAGE_FINISH:
        severity, message = "info", stage_finish(first, total, seconds, count)
    elif event == STAGE_FAILED:
        severity, message = "error", stage_failed(first, total, seconds, reason, count)
    elif event == STAGE_SKIPPED:
        severity, message = "info", stage_skipped(first, total, label, reason, count)
    else:
        raise ValueError(f"no stage event {event!r}")
    fields = {"event": event, "stage": first, "stages": total}
    if count != 1:
        fields["count"] = count
    if event in (STAGE_FINISH, STAGE_FAILED):
        fields["seconds"] = round(seconds, 1)
    return Record.at(time.time() if now is None else now, severity, component, message, **fields)


_COMPONENT_OF_PATH = (
    ("/planetgen/cli/warm", "warm"), ("/planetgen/web/warmup", "warm"),
    ("/planetgen/cli/tiles", "tiles"), ("/planetgen/web/tilebuild", "tiles"),
    ("/planetgen/cli/orbits", "orbits"),
    ("/planetgen/web/job_runner", "job"),
    ("/planetgen/queue/", "queue"),
    ("/planetgen/web/", "web"), ("/planetgen/api/", "web"), ("/planetgen/admin/", "web"),
    ("/planetgen/cli/generate", "generate"),
    ("/planetgen/cli/", "command"),
    ("/planetgen/util/joblog", "command"),
)
"""tuple: Source-file fragments and the component of what is logged from there, first match wins; anything else of planetGen's own is `generate`."""


def component_of(logger_name, pathname=""):
    """
    The component of a `logging` record (OPS.103): the logger's name after `planetgen.` when that is a component (`planetgen.queue`,
    `planetgen.activity`); else the one the source file belongs to (`_COMPONENT_OF_PATH`); another library's logger is `lib`.
    """
    if logger_name == "planetgen_app":
        return "web"            # Flask's own logger for the web app (`web/app.py`), outside the program's logger tree
    if logger_name != "planetgen" and not logger_name.startswith("planetgen."):
        return "lib"
    name = logger_name.removeprefix("planetgen.").lower().split(".")[0]      # `planetgen.queue.units` is the queue's
    if name in COMPONENTS:
        return name
    path = pathname.replace("\\", "/")
    for fragment, component in _COMPONENT_OF_PATH:
        if fragment in path:
            return component
    return "generate"


def counting_text(word, number, total):
    """`Part 2 of 5` (OPS.103): where a script or a run is, in one of the three `COUNTING_WORDS`."""
    if word not in COUNTING_WORDS:
        raise ValueError(f"no counting word {word!r}; the words are {', '.join(COUNTING_WORDS)}")
    return f"{word} {number} of {total}"


def _zone_name(moment):
    return moment.tzname() or moment.strftime("%z")


def display_line(record, tz=None, times=True, component=True, fields=True):
    """
    The record as the screen shows it. The time is local time, labelled with its zone (`2026-10-24 14:40:00 EDT`; `tz`
    is a `datetime.tzinfo`, the computer's own zone when `None`), then the severity, the component and the message,
    then the fields as `key=value` in the form of the file line. A stage's start, finish and failure are one divider,
    `=== Stage 4 of 12: label ===`, with no fields (the message says it all). `times`, `component` and `fields` leave
    their part out.
    """
    parts = []
    if times:
        moment = (_EPOCH + datetime.timedelta(milliseconds=record.ms)).astimezone(tz)
        parts.append(f"{moment.strftime('%Y-%m-%d %H:%M:%S')} {_zone_name(moment)}")
    parts.append(f"{record.severity.upper():<8}")
    if component:
        parts.append(f"{record.component}:")
    divider = record.event in STAGE_EVENTS
    message = record.message
    if message.lower().startswith(f"{record.severity}: "):
        message = message[len(record.severity) + 2:]         # the severity column says it: not "ERROR  migrate: error: ..."
    parts.append(f"{DIVIDER} {message} {DIVIDER}" if divider else message)
    if fields and not divider and record.fields:
        parts.append("(" + " ".join(f"{name}={_value_text(value)}" for name, value in record.fields.items()) + ")")
    return " ".join(parts)


_DISPLAY_SEVERITY = re.compile(r"^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d \S+ (DEBUG|INFO|NOTICE|WARNING|ERROR|CRITICAL) ")


def display_severity(line):
    """The severity (one of `SEVERITIES`) of a line `display_line` drew with its time, or `None` for any other line (the
    output of a command, a continuation). The admin pages colour a log by it (UX.284)."""
    found = _DISPLAY_SEVERITY.match(line)
    return found.group(1).lower() if found else None


def plain_line(record):
    """The command line's text without times: a divider for a stage's marks, the message for `debug`, `info` and
    `notice`, and the message after `warning: `, `error: ` or `critical: ` for the rest (unless the message already
    begins with that word, as `Error: ...` does)."""
    if record.event in STAGE_EVENTS:
        return f"{DIVIDER} {record.message} {DIVIDER}"
    if SYSLOG_CODES[record.severity] >= SYSLOG_CODES["notice"] or record.message.lower().startswith(record.severity):
        return record.message
    return f"{record.severity}: {record.message}"


SHELL_BEGIN = "# BEGIN generated vocabulary (scripts/build_log_vocabulary.py from planetgen/util/logrecord.py; do not edit)"
SHELL_END = "# END generated vocabulary"


def shell_vocabulary():
    """
    The vocabulary `scripts/logrecord.sh` is generated with (OPS.104): the severities, the ones that go to standard error, the word of a
    stage and the divider, as shell variable assignments between `SHELL_BEGIN` and `SHELL_END`, so the shell's side cannot name a word the
    module does not (`tests/test_log_shell.py` also runs both against the same records).
    """
    stage = stagemark.run(1).split(" ")[0]
    return "\n".join((SHELL_BEGIN, f'LOG_SEVERITIES="{" ".join(SEVERITIES)}"', f'LOG_STDERR_SEVERITIES="{" ".join(STDERR_SEVERITIES)}"',
                      f'LOG_STAGE_WORD="{stage}"', f'LOG_DIVIDER="{DIVIDER}"', f'LOG_COUNTING_WORDS="{" ".join(COUNTING_WORDS)}"',
                      f'LOG_COMPONENTS="{" ".join(COMPONENTS)}"', SHELL_END))
