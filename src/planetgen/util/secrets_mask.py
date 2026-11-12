"""
Masks keys and credentials in every log line (SEC.48): the first 4 characters, the last 2 and X's between (`abcdXXXXXXyz`); a value of 6
characters or fewer is all X's. One helper for the whole log layer: `logrecord.Record` (every field and message), the standard-logging
formatters and handlers (`log`), the job logs (`joblog`) and the `error.log` and Queue stream paths all call it, so no line reaches a screen,
a file or a stream with a secret in it.

A secret is found three ways:
- by the name beside it: `api_key=...`, `"password": "..."`, `Authorization: Bearer ...`, `--mysql-password ...`, and a log field with such a name;
- by its place in an address: `scheme://user:password@host` (the database and Redis URLs);
- by its value: every secret option of `config.json` / the environment (`x-secret`), every environment variable whose name says secret, and any
  value given to `register`, wherever it turns up.

Masking is idempotent: a masked value masks to itself.
"""

import os
import re
import threading

KEEP_FIRST = 4
KEEP_LAST = 2
FILL = "X"
MIN_KNOWN = 4
"""int: A known value shorter than this is not searched for (it would mask ordinary words); it is still masked by name."""

SECRET_WORDS = ("password", "passwd", "passphrase", "secret", "token", "api_key", "apikey", "api-key", "credential", "authorization",
                "private_key", "access_key", "session_id", "sessionid", "cookie", "csrf", "bearer")
"""tuple: A name containing one of these (any case, `-` or `_`) holds a secret."""

_WORD = "|".join(re.escape(word).replace("_", "[_-]?") for word in SECRET_WORDS if word != "bearer")
_NAME = rf"[\w.-]*(?:{_WORD}|(?<![a-z])key(?![a-z])|reset[_-]?(?:code|key|link))[\w.-]*"
_SENSITIVE_NAME = re.compile(rf"^{_NAME}$", re.IGNORECASE)
_VALUE = r"""(?:"[^"\r\n]*"|'[^'\r\n]*'|[^\s,;&"'<>)\]}]+)"""
_ASSIGNED = re.compile(rf"""(?P<name>(?<![\w.-])(?:--)?{_NAME})(?P<sep>["']?\s*(?:=|:)\s*|(?<=[a-z0-9])\s+(?=[^\s-]))(?P<value>{_VALUE})""", re.IGNORECASE)
_BEARER = re.compile(r"(?P<head>\b(?:Bearer|Basic)\s+)(?P<value>[A-Za-z0-9._~+/=-]+)", re.IGNORECASE)
_URL_CREDENTIALS = re.compile(r"(?P<head>\b[a-z][a-z0-9+.-]*://[^\s/:@\"']*:)(?P<value>[^\s/@\"']+)(?=@)", re.IGNORECASE)
_URL_USERONLY = re.compile(r"(?P<head>\b(?:redis|rediss|unix)://:)(?P<value>[^\s/@\"']+)(?=@)", re.IGNORECASE)

_lock = threading.Lock()
_registered = set()
_settings_seen = {"id": None, "values": frozenset()}
_ENV_NAME = re.compile(r"(?:PASSWORD|PASSWD|SECRET|TOKEN|API_?KEY|PRIVATE_KEY|ACCESS_KEY)", re.IGNORECASE)


def mask_value(value):
    """`value` masked: its first 4 characters, its last 2 and X's between; all X's when it has 6 characters or fewer. A value of `""` stays empty."""
    text = str(value)
    if len(text) <= KEEP_FIRST + KEEP_LAST:
        return FILL * len(text)
    return text[:KEEP_FIRST] + FILL * (len(text) - KEEP_FIRST - KEEP_LAST) + text[-KEEP_LAST:]


def is_secret_name(name):
    """Whether a field, option or variable called `name` holds a secret."""
    return bool(_SENSITIVE_NAME.match(str(name)))


def register(value):
    """Makes `value` a known secret: it is masked wherever it appears in a line, whatever stands beside it."""
    text = str(value) if value is not None else ""
    if len(text) >= MIN_KNOWN:
        with _lock:
            _registered.add(text)


def forget_registered():
    """Forgets the values given to `register` (tests)."""
    with _lock:
        _registered.clear()
        _settings_seen.update(id=None, values=frozenset())


def _settings_secrets():
    """The values of the secret options of the settings now in force (`x-secret`), found again when the settings object changes."""
    try:
        from planetgen.util import settings
        current = settings.get_settings()
    except Exception:  # noqa: BLE001 -- a log line never fails because the settings do
        return frozenset()
    if _settings_seen["id"] == id(current):
        return _settings_seen["values"]
    found = set()
    try:
        for path, _info, meta in settings.iter_fields():
            if meta.get("x-secret"):
                node = current
                for part in path.split("."):
                    node = getattr(node, part)
                if isinstance(node, str) and len(node) >= MIN_KNOWN:
                    found.add(node)
    except Exception:  # noqa: BLE001 -- as above
        pass
    values = frozenset(found)
    _settings_seen.update(id=id(current), values=values)
    return values


def known_values():
    """Every secret value to search for, longest first (so one that holds another is masked whole)."""
    values = set(_registered) | _settings_secrets()
    for name, value in os.environ.items():
        if len(value) >= MIN_KNOWN and _ENV_NAME.search(name):
            values.add(value)
    return sorted(values, key=len, reverse=True)


def _looks_secret(text):
    """Whether a word that only a space follows a secret name with looks like a key: 12 characters or more with a letter and a digit."""
    return len(text) >= 12 and any(c.isdigit() for c in text) and any(c.isalpha() for c in text)


def _assigned(match):
    value = match.group("value")
    quote = value[0] if value[:1] in "\"'" and value[-1:] == value[:1] and len(value) >= 2 else ""
    inner = value[1:-1] if quote else value
    if not inner:
        return match.group(0)
    if inner in ("<withheld>", "(set)"):
        return match.group(0)
    if match.group("sep").isspace() and not match.group("name").startswith("--") and not _looks_secret(inner):
        return match.group(0)            # "token ended" is a sentence; "token 9f8e7d6c5b4a" is a secret
    if inner.lower() in ("bearer", "basic"):
        return match.group(0)
    return f"{match.group('name')}{match.group('sep')}{quote}{mask_value(inner)}{quote}"


def mask_text(text):
    """`text` with every secret in it masked (see the module note); text without one comes back unchanged."""
    if not text or not isinstance(text, str):
        return text
    text = _URL_CREDENTIALS.sub(lambda m: m.group("head") + mask_value(m.group("value")), text)
    text = _URL_USERONLY.sub(lambda m: m.group("head") + mask_value(m.group("value")), text)
    text = _BEARER.sub(lambda m: m.group("head") + mask_value(m.group("value")), text)
    text = _ASSIGNED.sub(_assigned, text)
    for secret in known_values():
        if secret in text:
            text = text.replace(secret, mask_value(secret))
    return text


def mask_fields(fields):
    """A copy of the log fields `fields` with a secret-named field's value masked and every other string value scanned by `mask_text`."""
    if not fields:
        return fields
    shown = {}
    for name, value in fields.items():
        if isinstance(value, str):
            shown[name] = mask_value(value) if is_secret_name(name) and value not in ("<withheld>", "(set)") else mask_text(value)
        elif is_secret_name(name) and not isinstance(value, bool) and isinstance(value, (int, float)):
            shown[name] = mask_value(value)
        else:
            shown[name] = value
    return shown


__all__ = ["FILL", "KEEP_FIRST", "KEEP_LAST", "forget_registered", "is_secret_name", "known_values", "mask_fields", "mask_text", "mask_value", "register"]
