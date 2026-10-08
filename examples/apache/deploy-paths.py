"""
examples/apache/deploy-paths.py

Prints the two runtime directories create-cache-dir.sh makes for Apache,
one per line: the Galaxy Map tile cache (empty when the disk cache is off)
and the Generate page's jobs directory.

create-cache-dir.sh runs this as root, so it deliberately imports nothing
from the repo: it is run with `python3 -I` (no script directory, current
directory, user site-packages or PYTHON* variables on sys.path) and reads
config.json with the standard library alone. It mirrors
`tilecache.configured_cache_dir()` (src/planetgen/web/lib/tilecache.py) and
`jobs.jobs_dir()`'s configured directory (src/planetgen/web/jobs.py); keep the
three in step.

Usage:
    python3 -I examples/apache/deploy-paths.py <repo-dir>
    python3 -I examples/apache/deploy-paths.py --move-old-jobs <jobs-dir>

The second form (OPS.19) moves Generate jobs from the old default,
`OLD_JOBS_DIR`, into `<jobs-dir>` when that is the new default, and says
what it did.
"""

import json
import os
import shutil
import sys

DEFAULT_CACHE_DIR = "/var/cache/planetgen/tiles"
DEFAULT_JOBS_DIR = "/var/lib/planetGen/jobs"
OLD_JOBS_DIR = "/var/lib/planetgen/jobs"
"""The jobs default before OPS.19: `/var/lib/planetgen` differs from the
checkout's `/var/lib/planetGen` only by case, so a default install had
two folders."""
LOCK_NAME = "active"
"""The running job's lock in the jobs directory
(`planetgen.web.jobs.LOCK_NAME`): its first line is the job's id."""


def load_config(repo_dir):
    path = os.path.join(repo_dir, "config.json")
    if not os.path.isfile(path):
        return {}
    with open(path, "r", encoding="utf-8") as f:
        config = json.load(f)
    if not isinstance(config, dict):
        raise ValueError(f"{path} must contain a JSON object, not {type(config).__name__}")
    return config


def _section(config, name):
    section = config.get(name)
    return section if isinstance(section, dict) else {}


def tile_cache_dir(config):
    tile_cache = _section(config, "tile_cache")
    raw = os.environ.get("PLANETGEN_TILE_CACHE_MAX_MB")
    if raw is None:
        raw = tile_cache.get("max_mb", 200)
    try:
        max_bytes = max(0, int(float(raw) * 1024 * 1024))
    except (TypeError, ValueError):
        max_bytes = 200 * 1024 * 1024
    if max_bytes == 0:
        return ""
    return os.environ.get("PLANETGEN_TILE_CACHE_DIR") or tile_cache.get("dir") or DEFAULT_CACHE_DIR


def jobs_dir(config):
    return os.environ.get("PLANETGEN_JOBS_DIR") or _section(config, "jobs").get("dir") or DEFAULT_JOBS_DIR


def _finished(job_dir):
    """Whether the job in `job_dir` says it finished (`state.json`)."""
    try:
        with open(os.path.join(job_dir, "state.json"), "r", encoding="utf-8") as f:
            return json.load(f).get("status") in ("succeeded", "failed", "cancelled")
    except (OSError, ValueError, AttributeError):
        return False


def move_old_jobs(new_dir, old_dir=OLD_JOBS_DIR):
    """
    Moves each job folder from `old_dir` into `new_dir` (OPS.19), keeping
    its history, when `new_dir` is the default (a configured `jobs.dir`
    is left alone). A running job (named by `old_dir`'s lock) stays where
    it is, with its lock, and so does a job `new_dir` already has; the
    old folder, and its parent, are removed once empty.

    Returns:
        list[str]: What it did, one line each (empty when nothing).
    """
    if os.path.normpath(new_dir) != DEFAULT_JOBS_DIR or not os.path.isdir(old_dir) \
            or os.path.realpath(old_dir) == os.path.realpath(new_dir):
        return []
    running = None
    try:
        with open(os.path.join(old_dir, LOCK_NAME), "r", encoding="utf-8") as f:
            running = f.readline().strip() or None
    except OSError:
        pass
    if running and _finished(os.path.join(old_dir, running)):
        running = None  # a stale lock
    os.makedirs(new_dir, mode=0o750, exist_ok=True)
    moved, kept = [], []
    for name in sorted(os.listdir(old_dir)):
        source = os.path.join(old_dir, name)
        if not os.path.isdir(source) or os.path.islink(source):
            continue
        if name == running or os.path.exists(os.path.join(new_dir, name)):
            kept.append(name)
            continue
        shutil.move(source, os.path.join(new_dir, name))
        moved.append(name)
    lines = []
    if moved:
        lines.append(f"Moved {len(moved)} Generate job(s) from {old_dir} to {new_dir}.")
    if kept:
        lines.append(f"Left in {old_dir}: {', '.join(kept)} (running, or already in {new_dir}); "
                     f"move by hand once finished.")
    else:
        try:
            os.remove(os.path.join(old_dir, LOCK_NAME))  # a stale lock: no job left to hold it
        except OSError:
            pass
        try:
            os.rmdir(old_dir)
            lines.append(f"Removed the old {old_dir}.")
            os.rmdir(os.path.dirname(old_dir))
            lines.append(f"Removed the old {os.path.dirname(old_dir)}.")
        except OSError:
            pass
        if os.path.isdir(old_dir):
            lines.append(f"{old_dir} still holds other files; left in place.")
    return lines


def main(argv):
    if len(argv) == 3 and argv[1] == "--move-old-jobs":
        for line in move_old_jobs(argv[2]):
            print(line)
        return 0
    if len(argv) != 2:
        print("usage: python3 -I deploy-paths.py <repo-dir>", file=sys.stderr)
        return 2
    config = load_config(argv[1])
    print(tile_cache_dir(config))
    print(jobs_dir(config))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
