"""
examples/apache/deploy-paths.py

Prints the two runtime directories create-cache-dir.sh makes for Apache,
one per line: the Galaxy Map tile cache (empty when the disk cache is off)
and the Generate page's jobs directory.

create-cache-dir.sh runs this as root, so it deliberately imports nothing
from the repo: it is run with `python3 -I` (no script directory, current
directory, user site-packages or PYTHON* variables on sys.path) and reads
config.json with the standard library alone. It mirrors
`tilecache.configured_cache_dir()` (src/html/lib/tilecache.py) and
`jobs.jobs_dir()`'s configured directory (src/html/web/jobs.py); keep the
three in step.

Usage:
    python3 -I examples/apache/deploy-paths.py <repo-dir>
"""

import json
import os
import sys

DEFAULT_CACHE_DIR = "/var/cache/planetgen/tiles"
DEFAULT_JOBS_DIR = "/var/lib/planetgen/jobs"


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


def main(argv):
    if len(argv) != 2:
        print("usage: python3 -I deploy-paths.py <repo-dir>", file=sys.stderr)
        return 2
    config = load_config(argv[1])
    print(tile_cache_dir(config))
    print(jobs_dir(config))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
