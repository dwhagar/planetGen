# tests/test_settings.py

"""
The settings model (ADM.42, `planetgen.util.settings`): how the layers
combine, what each file may hold, and the drift tests that keep the model,
the generated files, the log-path module and the code's environment
variables in step.

Run with: pytest tests/test_settings.py
"""

import json
import os
import re

import pytest

from planetgen.cli import config as config_cli
from planetgen.util import logpaths, settings

SRC = os.path.join(logpaths.PROJECT_ROOT, "src", "planetgen")


@pytest.fixture
def files(tmp_path, monkeypatch):
    """Points the loader at a temporary config.json and settings.json."""
    config = tmp_path / "config.json"
    overlay = tmp_path / "settings.json"
    monkeypatch.setattr(logpaths, "CONFIG_PATH", str(config))
    monkeypatch.setenv("PLANETGEN_SETTINGS_FILE", str(overlay))
    settings.reset_cache()
    yield config, overlay
    settings.reset_cache()


def _load(config=None, overlay=None, environ=None):
    return settings.build(config, overlay, environ if environ is not None else {})


# --- Layers ---------------------------------------------------------------

def test_defaults_when_there_is_nothing():
    result = _load()
    assert result.site_name == "planetGen" and result.mysql.port == 3306
    assert result.proxy_fix.model_dump() == {"x_for": 0, "x_proto": 0, "x_host": 0}


def test_a_partial_section_keeps_the_rest_of_its_defaults():
    result = _load({"mysql": {"host": "db.example.com", "password": "secret"}})
    assert (result.mysql.host, result.mysql.password) == ("db.example.com", "secret")
    assert result.mysql.user == settings.Settings().mysql.user and result.debug is False


def test_the_layers_stack_in_order_defaults_config_overlay_environment():
    config = {"site_name": "From config", "tile_cache": {"max_mb": 10}, "jobs": {"keep": 5}}
    overlay = {"site_name": "From overlay", "tile_cache": {"max_mb": 20}}
    result = _load(config, overlay, {"PLANETGEN_TILE_CACHE_MAX_MB": "30"})
    assert result.site_name == "From overlay"
    assert result.tile_cache.max_mb == 30 and result.jobs.keep == 5


def test_an_empty_variable_counts_as_unset_except_where_the_model_says_otherwise():
    config = {"jobs": {"dir": "/from/config"}, "mysql": {"password": "kept?"}}
    result = _load(config, environ={"PLANETGEN_JOBS_DIR": "", "PLANETGEN_MYSQL_PASSWORD": ""})
    assert result.jobs.dir == "/from/config"
    assert result.mysql.password == ""      # a blank password set in the environment is a password


@pytest.mark.parametrize("raw, expected", [("off", False), ("0", False), ("false", False), ("1", True), ("yes", True)])
def test_a_boolean_reads_text_the_way_debug_always_did(raw, expected):
    assert _load({"debug": raw}).debug is expected
    assert _load(environ={"PLANETGEN_DEBUG": raw}).debug is expected


def test_the_allowlist_reads_comma_or_space_separated_text():
    result = _load(environ={"PLANETGEN_LOGIN_ALLOWLIST": "1.2.3.4, 10.0.0.0/8  ::1"})
    assert result.login_allowlist == ["1.2.3.4", "10.0.0.0/8", "::1"]


@pytest.mark.parametrize("bad", [-1, "one", True, 1.5])
def test_a_proxy_count_must_be_a_whole_number_0_or_above(bad):
    with pytest.raises(settings.SettingsError, match="proxy_fix.x_proto"):
        _load({"proxy_fix": {"x_proto": bad}})


def test_errors_name_every_field():
    with pytest.raises(settings.SettingsError) as raised:
        _load({"mysql": {"port": 99999}, "log_rotation": "sometimes", "tile_cache": {"max_mb": "lots"}})
    text = str(raised.value)
    assert "mysql.port" in text and "log_rotation" in text and "tile_cache.max_mb" in text


# --- The files ------------------------------------------------------------

def test_the_loader_rereads_when_a_file_changes(files):
    config, _overlay = files
    assert settings.get_settings().site_name == "planetGen"
    config.write_text(json.dumps({"site_name": "One"}))
    assert settings.get_settings().site_name == "One"
    config.write_text(json.dumps({"site_name": "Two!"}))
    assert settings.get_settings().site_name == "Two!"


def test_a_broken_or_non_object_config_file_fails_loudly(files):
    config, _overlay = files
    for raw in ("{not json", "[1, 2]", '"text"', "null"):
        config.write_text(raw)
        settings.reset_cache()
        with pytest.raises(settings.SettingsError, match="config.json"):
            settings.get_settings()


def test_the_overlay_may_set_only_editable_options(files):
    config, overlay = files
    overlay.write_text(json.dumps({"site_name": "Web edited", "tile_cache": {"max_mb": 5}}))
    result = settings.get_settings()
    assert result.site_name == "Web edited" and result.tile_cache.max_mb == 5
    for key in ({"mysql": {"host": "evil"}}, {"jobs": {"python": "/bin/sh"}}, {"secret_key": "x"},
                {"not_an_option": 1}):
        overlay.write_text(json.dumps(key))
        settings.reset_cache()
        with pytest.raises(settings.SettingsError, match="settings.json"):
            settings.get_settings()


def test_strict_mode_names_unknown_keys_in_config_json():
    _load({"ratelimt": {"default": "1 per day"}})                       # loading ignores them
    with pytest.raises(settings.SettingsError, match="ratelimt is not a setting"):
        settings.build({"ratelimt": {}, "mysql": {"hots": "x"}}, None, {}, strict=True)


def test_a_newer_file_is_refused_and_an_unversioned_one_is_version_1():
    assert _load({}).config_version == 1
    with pytest.raises(settings.SettingsError, match="newer"):
        _load({"config_version": settings.CONFIG_VERSION + 1})


# --- Drift ----------------------------------------------------------------

def test_the_generated_files_match_the_model():
    """`python -m planetgen.cli.config docs` rewrites them."""
    for path, content in config_cli.generated_files().items():
        with open(path, encoding="utf-8") as f:
            assert f.read() == content, f"{os.path.relpath(path, logpaths.PROJECT_ROOT)} is out of date; run " \
                                        "python3 -m planetgen.cli.config docs"


def test_the_example_file_holds_the_defaults():
    with open(config_cli.EXAMPLE_PATH, encoding="utf-8") as f:
        assert json.load(f) == settings.Settings().model_dump(mode="json")


def test_the_log_options_default_to_what_logpaths_reads():
    model = settings.Settings()
    assert (model.log_file, model.log_dir, model.log_rotation) == (
        logpaths.DEFAULTS["log_file"], logpaths.DEFAULTS["log_dir"], logpaths.DEFAULTS["log_rotation"])
    assert set(logpaths.DEFAULTS) <= {path for path, _i, _m in settings.iter_fields()}
    assert tuple(settings.Settings.model_fields["log_rotation"].annotation.__args__) == logpaths.LOG_ROTATION_MODES


def test_every_option_is_documented_and_has_complete_metadata():
    with open(config_cli.DOCS_PATH, encoding="utf-8") as f:
        docs = f.read()
    for path, _info, meta in settings.iter_fields():
        assert meta["description"], path
        assert f"`{path}`" in docs, path
        assert meta["x-min_role"] in ("admin", "owner"), path
        if meta.get("x-env"):
            assert f"`{meta['x-env']}`" in docs, path
    assert all(meta["x-restart"] or True for _p, _i, meta in settings.iter_fields())


def test_secrets_are_never_web_editable_by_admins_and_database_options_never_editable():
    for path, _info, meta in settings.iter_fields():
        if path.startswith("mysql.") and path != "mysql.statement_timeout_seconds":
            assert not meta["x-editable"], path
        if meta["x-secret"] and path.startswith(("mysql.", "secret_key")):
            assert not meta["x-editable"], path
    assert not dict((p, m) for p, _i, m in settings.iter_fields())["control_database"]["x-editable"]


# Variables the program reads that are not `config.json` options: switches
# for one process, tools and tests.
NON_OPTION_VARIABLES = {
    "PLANETGEN_WORKERS": "worker processes of one run",
    "PLANETGEN_PROGRESS_FILE": "a worker's progress channel",
    "PLANETGEN_GENERATION_STATS": "off in the test suite",
    "PLANETGEN_MIGRATIONS_DIR": "the migration scripts' folder",
    "PLANETGEN_SETTINGS_DIR": "the galaxy settings files' folder (ADM.18), not config.json",
    "PLANETGEN_SETTINGS_FILE": "the web-owned settings.json",
    "PLANETGEN_MYSQL_SQL_MODE": "the test suite's SQL mode",
    "PLANETGEN_TEST_REDIS_URL": "the test suite's Redis",
    "PLANETGEN_WORK_PARENT": "set by a job for its child processes",
    "PLANETGEN_WORKER_HOOKS": "set by a work queue for its workers",
}


def test_the_code_reads_no_variable_the_model_does_not_name():
    """A new `PLANETGEN_*` variable in code must be an option's `x-env` or
    be listed above with its reason."""
    known = set(settings.env_names()) | set(NON_OPTION_VARIABLES)
    stray = {}
    for root, _dirs, names in os.walk(SRC):
        for name in names:
            if not name.endswith(".py"):
                continue
            path = os.path.join(root, name)
            with open(path, encoding="utf-8") as f:
                text = f.read()
            for variable in re.findall(r"""["'](PLANETGEN_[A-Z0-9_]+)["']""", text):
                if variable not in known and not variable.endswith("_"):
                    stray.setdefault(variable, os.path.relpath(path, SRC))
    assert not stray, f"environment variables read in code but not in the settings model: {stray}"


def test_no_module_reads_the_file_directly_or_the_old_loader():
    offenders = []
    for root, _dirs, names in os.walk(SRC):
        for name in names:
            if name.endswith(".py") and name not in ("settings.py", "logpaths.py") \
                    and os.path.join(root, name) != os.path.join(SRC, "cli", "config.py"):
                with open(os.path.join(root, name), encoding="utf-8") as f:
                    text = f.read()
                if re.search(r"appconfig|load_config|DEFAULT_CONFIG|CONFIG_PATH", text):
                    offenders.append(name)
    assert not offenders, offenders
