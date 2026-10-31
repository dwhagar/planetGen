# planetgen/util/logrecord.py

"""
The logging layer (UX.121)
==========================

The one record every log line the program writes for a person or a tool to read is made from, and the three ways it is
written, so a file and a screen can never disagree. The standard, and every place the program logs today, are in
`docs/design/operations/logging-standard.md`.

A `Record` is a time (whole milliseconds, UTC), a severity (the RFC 5424 levels `debug`, `info`, `notice`, `warning`,
`error`, `critical`), a short component name, a message and ordered `key=value` fields for numbers and names a tool
may want (a count, an elapsed time). It is written as:

* a **file line** (`file_line`), logfmt, one record per line, the time always UTC::

      ts=2026-10-24T18:40:00.123Z level=info component=generate msg="Stage 4 of 12: Scatter the massive stars" stage=4 stages=12

* a **JSON line** (`json_line`), the same record as JSON Lines, for machines;
* a **display line** (`display_line`), for the command line and the web log: the time in local time and labelled with its
  zone, the severity, the component, the message and the fields, with a stage's start, finish and failure drawn as one
  divider (`=== Stage 4 of 12: label ===`). `plain_line` is the command line's untimed text.

`parse_file_line` and `parse_json_line` read a line back to the record it came from, so what a file holds and what the
screen shows are the same records. Python's `logging` stays the entry point: `get_logger(component)` is a standard
logger that takes fields as keyword arguments, and `FileFormatter`, `JsonFormatter` and `DisplayFormatter` write what
it logs. Standard library only.
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
"""int: The `logging` level of `notice`, between `info` (20) and `warning` (30), named in `logging` by `get_logger`."""

LOGGING_LEVELS = {"debug": logging.DEBUG, "info": logging.INFO, "notice": NOTICE, "warning": logging.WARNING,
                  "error": logging.ERROR, "critical": logging.CRITICAL}
"""dict: The `logging` level of each severity."""

RESERVED_FIELDS = ("ts", "level", "component", "msg")
"""tuple: The names the record itself uses in a file line and a JSON line; a field cannot take one."""

COMPONENT_PATTERN = re.compile(r"^[a-z][a-z0-9_.-]{0,23}$")
FIELD_PATTERN = re.compile(r"^[a-z][a-z0-9_]{0,31}$")

STAGE_START, STAGE_FINISH, STAGE_FAILED = "stage.start", "stage.finish", "stage.failed"
STAGE_SKIPPED = "stage.skipped"
"""str: The `event` field of the records that mark a stage (it starts, it finishes, it fails, it is skipped); the display
draws each as a divider."""

STAGE_EVENTS = (STAGE_START, STAGE_FINISH, STAGE_FAILED, STAGE_SKIPPED)

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
        if not COMPONENT_PATTERN.match(component):
            raise ValueError(f"component {component!r} is not a short lower-case name (a-z, 0-9, _ . -)")
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


def _utc(ms):
    """`"2026-10-24T18:40:00.123Z"` for `ms` milliseconds since the epoch."""
    moment = _EPOCH + datetime.timedelta(milliseconds=ms)
    return moment.strftime("%Y-%m-%dT%H:%M:%S") + f".{moment.microsecond // 1000:03d}Z"


def _ms_of(stamp):
    """The milliseconds since the epoch of a `_utc` text.

    Raises:
        ValueError: Not a time in that form.
    """
    match = re.match(r"^(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d)\.(\d{3})Z$", stamp)
    if not match:
        raise ValueError(f"{stamp!r} is not a time like 2026-10-24T18:40:00.123Z")
    moment = datetime.datetime.strptime(match.group(1), "%Y-%m-%dT%H:%M:%S").replace(tzinfo=datetime.timezone.utc)
    return int((moment - _EPOCH) / datetime.timedelta(milliseconds=1)) + int(match.group(2))


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
    """The record as a file line (logfmt, the time in UTC, no line break)."""
    parts = [f"ts={_utc(record.ms)}", f"level={record.severity}", f"component={record.component}",
             f"msg={_quote(record.message)}"]
    parts.extend(f"{name}={_value_text(value)}" for name, value in record.fields.items())
    return " ".join(parts)


def json_line(record):
    """The record as a JSON Lines line: the same names as the file line, the fields beside them."""
    body = {"ts": _utc(record.ms), "level": record.severity, "component": record.component, "msg": record.message}
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
    parts.append(f"=== {record.message} ===" if divider else record.message)
    if fields and not divider and record.fields:
        parts.append("(" + " ".join(f"{name}={_value_text(value)}" for name, value in record.fields.items()) + ")")
    return " ".join(parts)


def plain_line(record):
    """The command line's text without times: a divider for a stage's marks, the message for `debug`, `info` and
    `notice`, and the message after `warning: `, `error: ` or `critical: ` for the rest (unless the message already
    begins with that word, as `Error: ...` does)."""
    if record.event in STAGE_EVENTS:
        return f"=== {record.message} ==="
    if SYSLOG_CODES[record.severity] >= SYSLOG_CODES["notice"] or record.message.lower().startswith(record.severity):
        return record.message
    return f"{record.severity}: {record.message}"


def record_of(log_record):
    """The `Record` of a standard `logging.LogRecord`: its level as a severity, the logger's name (less a leading
    `planetgen.`; for the program's own logger the module that logged) as the component, the fields `get_logger` attached, and its time."""
    level = log_record.levelno
    severity = "debug"
    for name in SEVERITIES:
        if level >= LOGGING_LEVELS[name]:
            severity = name
    component = log_record.name.removeprefix("planetgen.").lower() or "planetgen"
    if component == "planetgen" and getattr(log_record, "module", None):
        component = log_record.module.lower()          # the program's own logger: the module that logged is the component
    component = re.sub(r"[^a-z0-9_.-]", "_", component)[:24]
    if not COMPONENT_PATTERN.match(component):
        component = "planetgen"
    return Record.at(log_record.created, severity, component, log_record.getMessage(),
                     **getattr(log_record, "fields", {}))


class FileFormatter(logging.Formatter):
    """Writes what a logger logs as file lines (`file_line`)."""

    def format(self, log_record):
        return file_line(record_of(log_record))


class JsonFormatter(logging.Formatter):
    """Writes what a logger logs as JSON lines (`json_line`)."""

    def format(self, log_record):
        return json_line(record_of(log_record))


class DisplayFormatter(logging.Formatter):
    """Writes what a logger logs as display lines (`display_line`, or `plain_line` with `plain=True`)."""

    def __init__(self, tz=None, plain=False, **parts):
        super().__init__()
        self.tz, self.plain, self.parts = tz, plain, parts

    def format(self, log_record):
        record = record_of(log_record)
        return plain_line(record) if self.plain else display_line(record, self.tz, **self.parts)


class Logger:
    """
    A standard logger for one component that takes fields as keyword arguments: `log.info("Placed sectors", count=12,
    seconds=1.5)`. Its stage methods log the start, finish and failure of a stage in the one wording.
    """

    def __init__(self, component):
        self.component = component
        self.logger = logging.getLogger(f"planetgen.{component}")

    def log(self, severity, message, /, **fields):
        if self.logger.isEnabledFor(LOGGING_LEVELS[severity]):
            Record.now(severity, self.component, str(message), **fields)       # a bad field is the caller's bug, now
            self.logger.log(LOGGING_LEVELS[severity], message, extra={"fields": fields})

    def debug(self, message, /, **fields):
        self.log("debug", message, **fields)

    def info(self, message, /, **fields):
        self.log("info", message, **fields)

    def notice(self, message, /, **fields):
        self.log("notice", message, **fields)

    def warning(self, message, /, **fields):
        self.log("warning", message, **fields)

    def error(self, message, /, **fields):
        self.log("error", message, **fields)

    def critical(self, message, /, **fields):
        self.log("critical", message, **fields)

    def stage(self, event, first, total, label=None, seconds=0.0, reason=None, count=1):
        """Logs the start, finish or failure of a stage (`stage_record`)."""
        record = stage_record(event, self.component, first, total, label, seconds, reason, count)
        self.log(record.severity, record.message, **record.fields)


def get_logger(component):
    """The `Logger` of `component`.

    Raises:
        ValueError: A component name that is not allowed.
    """
    if not COMPONENT_PATTERN.match(component):
        raise ValueError(f"component {component!r} is not a short lower-case name (a-z, 0-9, _ . -)")
    return Logger(component)
