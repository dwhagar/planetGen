# planetgen/util/logpaths.py

"""
Where planetGen's logs go (the stdlib-only part of the configuration)
=====================================================================

The debug switch and the two log locations, read straight from
`config.json` (or the environment) with the standard library alone. They
live apart from `planetgen.util.settings`, which holds every other option
as one validated Pydantic model, because `examples/apache/
log-locations.py` runs as root under the system's `python3 -I`, before any
virtual environment exists, and loads this file by its path: it can't
import Pydantic. These options are the settings model's `x-editable:
false` ones (file-only, restart to apply), so reading them here and
nowhere else never disagrees with `settings.json`; `settings.py` takes its
defaults for them from the constants below, and a test fails if the two
drift.

Precedence: the matching `PLANETGEN_*` environment variable, then
`config.json`, then the default.
"""

import json
import os
import sys

# This file is src/planetgen/util/logpaths.py: four levels up is the repo
# root.
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

CONFIG_PATH = os.path.join(PROJECT_ROOT, "config.json")
"""str: Absolute path to the real config file, at the repo root -- a
sibling of `html/`, `db/`, and `src/`, deliberately outside `html/`'s
Apache `DocumentRoot`. See `docs/config.md` for why the location
matters."""

DEFAULT_LOG_FILE = "/var/log/planetgen.log"
"""str: The debug log when `log_file` and `PLANETGEN_LOG_FILE` are unset."""

DEFAULT_LOG_ROTATION = "auto"
"""str: `log_rotation`'s default."""

DEFAULTS = {"log_file": DEFAULT_LOG_FILE, "log_dir": "", "log_rotation": DEFAULT_LOG_ROTATION}
"""dict: Each log option's default, for a tool that reports where a value
came from (`examples/apache/log-locations.py`)."""

LOG_ROTATION_MODES = ("auto", "system", "app")
"""tuple: What `log_rotation` accepts."""


def read_config_file(path=None):
    """
    `config.json` as a dict (empty when the file doesn't exist).

    Raises:
        ValueError: The file isn't valid JSON, or isn't a JSON object.
    """
    path = path or CONFIG_PATH
    if not os.path.isfile(path):
        return {}
    with open(path, "r", encoding="utf-8") as f:
        config = json.load(f)
    if not isinstance(config, dict):
        # Same failure mode as malformed JSON (json.JSONDecodeError is a
        # ValueError), with a message that says what's actually wrong.
        raise ValueError(f"{path} must contain a JSON object, not {type(config).__name__}")
    return config


_FALSE_STRINGS = ("", "0", "false", "no", "off")


def debug_enabled(config=None):
    """
    Whether this deployment's debug mode is on: the `PLANETGEN_DEBUG`
    environment variable when set (`0`/`false`/`no`/`off`/empty mean off,
    anything else on), else `config.json`'s `"debug"`, which defaults to
    off when missing (a string there is read the same way as the variable). Debug mode turns on the verbose debug log (see
    `planetgen.util.log`) and the web interface's traceback-in-page 500
    responses.

    Args:
        config (dict, optional): An already-loaded `read_config_file()` result,
            to skip re-reading the file.

    Returns:
        bool
    """
    env = os.environ.get("PLANETGEN_DEBUG")
    if env is not None:
        return env.strip().lower() not in _FALSE_STRINGS
    if config is None:
        config = read_config_file()
    value = config.get("debug")
    if isinstance(value, str):
        # "false"/"0"/"off"/"no" in config.json mean off, exactly as they
        # do for PLANETGEN_DEBUG -- bool("false") would be True.
        return value.strip().lower() not in _FALSE_STRINGS
    return bool(value)


def log_file_path(config=None):
    """
    Where the debug log goes: `PLANETGEN_LOG_FILE` when set, else
    `config.json`'s `"log_file"`, else `/var/log/planetgen.log`.

    Args:
        config (dict, optional): An already-loaded `read_config_file()` result.

    Returns:
        str
    """
    env = os.environ.get("PLANETGEN_LOG_FILE")
    if env:
        return env
    if config is None:
        config = read_config_file()
    return config.get("log_file") or DEFAULT_LOG_FILE


ACTIVITY_LOG_NAME = "planetgen.log"
"""str: The activity log's file name inside `log_dir_path()`."""

SYSTEM_ROTATION_FILES = ("/etc/logrotate.d/planetgen-log", "/etc/newsyslog.d/planetgen-log.conf")
"""tuple: The rotation configs `examples/apache/setup-debug-log.sh`
installs for the activity log (Linux, macOS). While one exists, `"auto"`
rotation leaves the file to the system tool."""


def default_log_dir(platform=None):
    """
    The platform's standard folder for the activity log: `/var/log/
    planetgen` on Linux (and other Unix systems), `/Library/Logs/planetgen`
    on macOS, and `logs` under the checkout on Windows.

    Args:
        platform (str, optional): A `sys.platform` value (default this
            one's).
    """
    platform = platform or sys.platform
    if platform.startswith("win"):
        return os.path.join(PROJECT_ROOT, "logs")
    if platform == "darwin":
        return "/Library/Logs/planetgen"
    return "/var/log/planetgen"


def log_dir_path(config=None):
    """
    The activity log's folder: `PLANETGEN_LOG_DIR` when set, else
    `config.json`'s `"log_dir"`, else `default_log_dir()`.
    """
    env = os.environ.get("PLANETGEN_LOG_DIR")
    if env:
        return env
    if config is None:
        config = read_config_file()
    value = config.get("log_dir")
    if value and not isinstance(value, str):
        raise TypeError(f'"log_dir" must be a string, not {type(value).__name__}')
    return value or default_log_dir()


def activity_log_path(config=None):
    """
    The activity log file: `planetgen.log` in `log_dir_path()`, or
    `planetgen-activity.log` there when that would be the very file the
    debug log uses (`log_file_path()`), so the two never share a file.
    """
    if config is None:
        config = read_config_file()
    folder = os.path.abspath(log_dir_path(config))
    path = os.path.join(folder, ACTIVITY_LOG_NAME)
    if os.path.normcase(path) == os.path.normcase(os.path.abspath(log_file_path(config))):
        path = os.path.join(folder, "planetgen-activity.log")
    return path


def log_rotation_mode(config=None):
    """
    How the activity log is rotated: `"system"` (logrotate or newsyslog
    moves the file; the program reopens it) or `"app"` (the program
    rotates it itself). `config.json`'s `"log_rotation"` picks one;
    `"auto"` (the default) means `"system"` where one of
    `SYSTEM_ROTATION_FILES` exists, else `"app"` (always on Windows).
    """
    if config is None:
        config = read_config_file()
    value = config.get("log_rotation")
    value = value.strip().lower() if isinstance(value, str) else "auto"
    if value in ("system", "app"):
        return value
    if sys.platform.startswith("win"):
        return "app"
    return "system" if any(os.path.exists(p) for p in SYSTEM_ROTATION_FILES) else "app"
