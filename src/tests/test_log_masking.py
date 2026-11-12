"""SEC.48: keys, tokens, passwords and credentialed addresses are masked in every log line (first 4 characters, last 2, X's between)."""

import io
import json
import logging

import pytest

from planetgen.util import joblog, log, logrecord, secrets_mask

KEY = "ak_live_0123456789abcdefXYZ"
TOKEN = "tok-9f8e7d6c5b4a39281706"
PASSWORD = "Hunter2-correct-horse"
RESET = "7c1d5e90a2b34f68"
SESSION = "sess-aabbccddeeff00112233"
KNOWN = "plainly-known-value-77"
SECRETS = (KEY, TOKEN, PASSWORD, RESET, SESSION, KNOWN)
DB_URL = f"mysql+pymysql://planet:{PASSWORD}@127.0.0.1:3306/planetgen"
REDIS_URL = f"redis://:{TOKEN}@127.0.0.1:6379/0"


def leaks(text):
    return [s for s in SECRETS if s in text]


@pytest.fixture(autouse=True)
def known_value():
    secrets_mask.register(KNOWN)
    yield
    secrets_mask.forget_registered()


def test_mask_value_shape():
    assert secrets_mask.mask_value("abcdefghij") == "abcdXXXXij"
    assert secrets_mask.mask_value("abcdefg") == "abcdXfg"[:4] + "X" + "fg"
    assert secrets_mask.mask_value("abcdef") == "XXXXXX"
    assert secrets_mask.mask_value("ab") == "XX"
    assert secrets_mask.mask_value("") == ""


def test_masking_is_idempotent():
    line = f"api_key={KEY} {DB_URL} password: {PASSWORD}"
    once = secrets_mask.mask_text(line)
    assert secrets_mask.mask_text(once) == once


def test_text_by_name_by_address_and_by_value():
    line = (f"api_key={KEY} token: {TOKEN} {{\"password\": \"{PASSWORD}\"}} reset_token={RESET} --mysql-password {PASSWORD} "
            f"Authorization: Bearer {SESSION} db {DB_URL} redis {REDIS_URL} loose {KNOWN}")
    assert leaks(secrets_mask.mask_text(line)) == []
    assert "ak_lXXXXXXXXXXXXXXXXXXXYZ"[:4] in secrets_mask.mask_text(f"api_key={KEY}")


def test_ordinary_words_are_left_alone():
    line = "Placed sectors count=12 session ended; keyboard monkey basalt_fraction=0.3 the password is wrong"
    assert secrets_mask.mask_text(line) == line


def test_fields_by_name():
    shown = secrets_mask.mask_fields({"api_key": KEY, "session_id": SESSION, "count": 3, "note": f"x {DB_URL}"})
    assert leaks(json.dumps(shown)) == []
    assert shown["count"] == 3


def test_record_holds_no_secret():
    record = logrecord.Record.now("error", "web", f"login failed for {DB_URL}", api_key=KEY, token=TOKEN)
    assert leaks(repr(record)) == []
    assert leaks(logrecord.json_line(record)) == [] and leaks(logrecord.file_line(record)) == []
    assert leaks(logrecord.display_line(record)) == []


@pytest.mark.parametrize("form", ["display", "file", "json", "plain", "raw", "debug"])
def test_every_formatter_masks(form):
    handler = logging.Handler()
    handler.setFormatter(log.LineFormatter(form, traceback=True))
    try:
        raise RuntimeError(f"connect failed {REDIS_URL} password={PASSWORD}")
    except RuntimeError:
        import sys
        record = logging.getLogger("planetgen.web").makeRecord("planetgen.web", logging.ERROR, __file__, 1,
                                                               f"api_key=%s {DB_URL}", (KEY,), sys.exc_info())
    record.origin = "web[1]"
    record.fields = {"token": TOKEN, "note": RESET and f"reset_token={RESET}"}
    assert leaks(handler.format(record)) == []


def test_every_logger_path_is_masked(caplog):
    stream = io.StringIO()
    handler = logging.StreamHandler(stream)                      # a plain handler a library might add, no planetgen formatter
    logger = logging.getLogger("some.library")
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    try:
        try:
            raise ValueError(f"bad {DB_URL} token={TOKEN}")
        except ValueError:
            logger.exception("failed with api_key=%s", KEY)
        log.warning(f"password={PASSWORD}")
        log.get("web").error("login", token=TOKEN, session_id=SESSION)
    finally:
        logger.removeHandler(handler)
    assert leaks(stream.getvalue()) == []
    assert leaks(caplog.text) == []


def test_error_log_job_log_callback_and_relay_are_masked(tmp_path):
    seen = []
    callback = log.callback_handler(lambda message, level: seen.append(message), logging.DEBUG)
    log.attach(callback)
    path = tmp_path / "job.log"
    try:
        with joblog.JobLog("warm", "warm", "Masking test", path=str(path)) as job:
            job.say(f"using api_key={KEY} and {DB_URL}")
            log.error(f"token {TOKEN}", stderr=False)
    finally:
        log.detach(callback)
    assert leaks(path.read_text()) == []
    assert leaks("\n".join(seen)) == []

    relay_lines = []
    handler = log._RelayHandler.__new__(log._RelayHandler)       # the relay's own formatting and masking, without a Redis
    logging.Handler.__init__(handler)
    handler._lines, handler._guard = relay_lines, __import__("threading").Lock()
    handler.flush = lambda: None
    handler.setFormatter(log.LineFormatter("file"))
    handler.emit(logging.LogRecord("planetgen.queue", logging.INFO, "", 0, f"x api_key={KEY}", None, None))
    assert leaks("\n".join(relay_lines)) == []


def test_write_line_masks(capsys):
    log.write_line(f"another process said password={PASSWORD}")
    assert leaks(capsys.readouterr().out) == []
