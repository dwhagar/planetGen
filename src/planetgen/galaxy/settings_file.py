# planetgen/galaxy/settings_file.py

"""
The galaxy's creation settings as a JSON file (ADM.18;
docs/design/reproducible-galaxies.md, section 7).

When a galaxy is planned, a file holding only what reproduction needs is
written: every plan option (defaults included), the 128-bit seed, the version
key with its parts spelled out, the naming key (GEN.70), the SHA-256 of
`requirements.lock`, and the word lists the name generator checks against
(gzip-compressed and base64-encoded, each with its SHA-256). It is named
`<32-hex seed>-<22-hex version key>-<YYYYMMDD>-<HHMMSS>Z.json` in UTC, and lives in `galaxy-settings` inside the Generate page's jobs
directory (`jobs.dir` in `config.json`, or `PLANETGEN_SETTINGS_DIR`).

The file is written once for a plan. Planning again with other settings
writes a new file under the new date-time and leaves the old one as a dated
backup, so the newest file for the galaxy's seed is the current one.
"""

import base64
import datetime
import gzip
import hashlib
import json
import os
import re
import sys

from planetgen import tuning
from planetgen._version import __version__
from planetgen.galaxy import seed as galaxy_seed, version_history, version_key

FORMAT = 1
"""int: The file's layout version."""

DIR_ENV_VAR = "PLANETGEN_SETTINGS_DIR"
DIR_NAME = "galaxy-settings"

PLAN_SETTINGS = (
    "disk_scale_length_pc", "disk_scale_height_pc", "bulge_scale_radius_pc", "bulge_amplitude", "arm_count",
    "pitch_angle_deg", "arm_density", "interarm_density", "core_density", "calibration_radius_pc", "max_ring", "bright_star_min_luminosity",
    "no_bright_stars", "phenomenon_min_mass", "compact_min_mass",
)
"""tuple: The `planetgen plan` options that shape the galaxy. The prevalence
options belong to the generate runs that follow (kept in the run history), not
to the plan."""

_NAME = re.compile(r"^([0-9A-F]{32})-([0-9A-F]{22})-(\d{8})-(\d{6})Z\.json$")


def settings_dir():
    """Where the files live."""
    configured = os.environ.get(DIR_ENV_VAR)
    if configured:
        return configured
    from planetgen.web import jobs  # not at the top: the web layer sits above this one
    return os.path.join(jobs.configured_jobs_dir(), DIR_NAME)


def file_name(seed, key, when):
    """`<seed>-<key>-<YYYYMMDD>-<HHMMSS>Z.json` for a UTC `when`."""
    return f"{bytes(seed).hex().upper()}-{key}-{when:%Y%m%d-%H%M%S}Z.json"


def encode_words(words):
    """A word set as `{"count", "sha256", "gzip_base64"}`: the sorted words,
    one per line, hashed, then gzip-compressed (no timestamp) and base64-encoded."""
    text = "\n".join(sorted(words)).encode("utf-8")
    return {
        "count": len(words),
        "sha256": hashlib.sha256(text).hexdigest(),
        "gzip_base64": base64.b64encode(gzip.compress(text, mtime=0)).decode("ascii"),
    }


def decode_words(entry):
    """The word set `encode_words` made (its hash is checked)."""
    text = gzip.decompress(base64.b64decode(entry["gzip_base64"]))
    if hashlib.sha256(text).hexdigest() != entry["sha256"]:
        raise ValueError("the stored word list does not match its SHA-256")
    return set(text.decode("utf-8").split("\n")) if text else set()


def version_parts():
    """The version key with its parts spelled out."""
    major, revision, build = (int(part) for part in __version__.split("."))
    return {
        "version_key": version_key.version_key(),
        "planetgen": {"major": major, "revision": revision, "build": build, "version": __version__},
        "python": {"version": version_key.python_version(), "major": sys.version_info[0],
                   "minor": sys.version_info[1], "micro": sys.version_info[2]},
        "platform": version_key.platform_name(),
        "os_digit": version_key.os_digit(),
        "architecture_digit": version_key.arch_digit(),
    }


def build(settings, seed, naming_key=None, shape=None, words=None):
    """
    The file's content, as a dict (no timestamp: the name carries it).

    Args:
        settings (dict): The plan options (`PLAN_SETTINGS`), defaults included.
        seed (bytes): The 16-byte galaxy seed.
        naming_key (str, optional): The galaxy's naming key (GEN.70).
        shape (dict, optional): The resolved skeleton numbers (edge, outer ring...).
        words (dict, optional): `{"dictionary": set, "offensive": set}`; the
            name generator's lists by default.
    """
    if words is None:
        from planetgen.names import wordlists
        words = {"dictionary": wordlists.DICTIONARY_WORDS, "offensive": wordlists.NSFW_WORDS}
    return {
        "format": FORMAT,
        "seed": galaxy_seed.format_seed(seed),
        "version": version_parts(),
        "settings": {name: settings.get(name) for name in PLAN_SETTINGS},
        "skeleton": shape or {},
        "naming_key": naming_key,
        "hashes": {"requirements_lock_sha256": version_history.requirements_sha256()},
        "sector_edge_pc": float(tuning.DEFAULT_SECTOR_EDGE_PC),
        "word_lists": {name: encode_words(found) for name, found in words.items()},
    }


def list_files(directory=None):
    """`[{"name", "seed", "key", "when", "path", "size"}]` for every settings
    file in the directory, newest first."""
    directory = directory or settings_dir()
    found = []
    try:
        names = os.listdir(directory)
    except OSError:
        return []
    for name in names:
        match = _NAME.match(name)
        if match is None:
            continue
        when = datetime.datetime.strptime(match.group(3) + match.group(4), "%Y%m%d%H%M%S")
        path = os.path.join(directory, name)
        found.append({"name": name, "seed": match.group(1), "key": match.group(2), "when": when,
                      "path": path, "size": os.path.getsize(path)})
    found.sort(key=lambda entry: (entry["when"], entry["name"]), reverse=True)
    return found


def current_file(seed, directory=None):
    """The newest file for `seed` (bytes), or `None`."""
    wanted = bytes(seed).hex().upper()
    for entry in list_files(directory):
        if entry["seed"] == wanted:
            return entry
    return None


def _same_content(document, path):
    try:
        with open(path, "r", encoding="utf-8") as handle:
            return json.load(handle) == document
    except (OSError, ValueError):
        return False


def write(document, seed, directory=None, now=None):
    """
    Writes `document` as the newest file for `seed`, unless the newest one
    already holds exactly this (then its entry is returned and nothing is
    written). An earlier file stays as a dated backup. Written atomically,
    readable by the web server's user.

    Returns:
        dict: The file's `list_files` entry.
    """
    directory = directory or settings_dir()
    current = current_file(seed, directory)
    if current is not None and _same_content(document, current["path"]):
        return current
    os.makedirs(directory, mode=0o750, exist_ok=True)
    now = (now or datetime.datetime.now(datetime.timezone.utc)).replace(tzinfo=None, microsecond=0)
    path = os.path.join(directory, file_name(seed, document["version"]["version_key"], now))
    temporary = path + ".tmp"
    with open(temporary, "w", encoding="utf-8") as handle:
        json.dump(document, handle, indent=1, sort_keys=True)
        handle.write("\n")
    os.chmod(temporary, 0o644)
    os.replace(temporary, path)
    return next(entry for entry in list_files(directory) if entry["path"] == path)
