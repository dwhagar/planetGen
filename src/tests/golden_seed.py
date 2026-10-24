# tests/golden_seed.py

"""
The golden-seed records (TEST.77): what a fixed seed makes, in the three forms the golden test pins.

- `battery_record()`: the database-free battery (OPS.57) as a digest per case (tier A);
- `rosetta_text()`: three of the battery's smallest cases as readable canonical text, one value a line, so a failing
  digest can be read as a diff;
- `galaxy_record(conn)`: the fingerprint (GEN.58) of a small stored galaxy, a digest per sector, per named section of
  each sector, and for the region (tier B).

They are kept under `src/tests/golden/` as `battery.json`, `galaxy.json` and `rosetta.txt`. A change of what the
generator makes changes them, and the generator epoch goes up with them (`tests/golden_update.py`, OPS.56). With
`PLANETGEN_WRITE_GOLDEN=1` the golden test writes the files instead of comparing; `golden_update` sets it, and nothing
sets it under CI.
"""

import json
import os
import platform
import sys
from pathlib import Path

from planetgen.db import fingerprint
from planetgen.generation import generator_battery as battery
from planetgen.util import fingerprint_encoding as encoding

GOLDEN_DIR = Path(__file__).resolve().parent / "golden"
"""Path: Where the golden files are kept."""

BATTERY_FILE = GOLDEN_DIR / "battery.json"
GALAXY_FILE = GOLDEN_DIR / "galaxy.json"
ROSETTA_FILE = GOLDEN_DIR / "rosetta.txt"

WRITE_ENVIRONMENT = "PLANETGEN_WRITE_GOLDEN"
"""str: The environment variable that makes the golden test write the files instead of comparing."""

GALAXY_SEED = bytes.fromhex("5EED00000000A1B2C3D4E5F607182930")
"""bytes: The 128-bit seed of the golden galaxy."""

GALAXY_ARGUMENTS = ("--shells", "1", "--num-systems", "4")
"""tuple: The `planetgen galaxy` arguments (before `--workers`) that build the golden galaxy: shell 1 with four
systems a sector, then the run's own scatter of the bright stars and the phenomena, its backfill and its settle."""

ROSETTA_CASES = ("phenomenon", "nebula-shape", "unit-seed")
"""tuple: The battery cases written out in full."""


def writing():
    """Whether this run rewrites the golden files; never under CI."""
    return bool(os.environ.get(WRITE_ENVIRONMENT)) and not os.environ.get("CI")


def environment():
    """What the files were made on, for the reader only: the digests are the same on every leg."""
    found = {"python": platform.python_version(), "platform": platform.platform()}
    for name in ("numpy", "astropy"):
        try:
            found[name] = __import__(name).__version__
        except ImportError:
            found[name] = None
    return found


def battery_record(found=None):
    """The battery's digests as the golden file keeps them: the overall one and, per case, `[full, structural]`."""
    found = battery.generator_battery_digest() if found is None else found
    return {
        "overall": found.overall,
        "cases": {name: [found.cases[name].full, found.cases[name].structural] for name in battery.CASES},
    }


def rosetta_text():
    """The text of `ROSETTA_CASES` in canonical form, one value a line (a break after each comma), under a heading
    for each case."""
    parts = []
    for name in ROSETTA_CASES:
        text = encoding.encode(battery.case_value(name))
        parts.append(f"== {name} ==\n" + text.replace(",", ",\n") + "\n")
    return "\n".join(parts)


def galaxy_record(conn):
    """The golden galaxy's fingerprint: `sectors` maps a sector's address to its digest and to
    `{section: [full, structural]}` (for a sector that has systems; the others are told apart by the digest alone), and `region` is the digest over the plan and every sector."""
    found = fingerprint.region_fingerprint(conn)
    sectors = {}
    for label, digest in found.sectors:
        address = tuple(int(part) for part in label.split())
        sections = fingerprint.sector_sections(conn, address=address)
        kept = sections if "system" in sections else {}   # only a sector with systems is told apart by section
        sectors[label] = {"digest": digest,
                          "sections": {name: [value.full, value.structural] for name, value in kept.items()}}
    return {"plan": found.plan, "region": found.region, "sectors": sectors}


def read(path):
    """The golden file at `path`: parsed for a `.json` one, else its text."""
    text = Path(path).read_text(encoding="utf-8")
    return json.loads(text) if str(path).endswith(".json") else text


def write(path, value):
    """Writes `value` (JSON for a `.json` file, else text) as the golden file at `path`."""
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    text = json.dumps(value, indent=2, sort_keys=True) + "\n" if str(path).endswith(".json") else value
    Path(path).write_text(text, encoding="utf-8")


def written_battery():
    """The `battery.json` content: the record plus where it was made."""
    return {**battery_record(), "made_on": environment()}


def describe_difference(expected, actual):
    """One line per case whose digest differs: a structural difference is a change of the code, a difference in
    the floats alone a change of the platform or its libraries."""
    lines = []
    for name in battery.CASES:
        want, got = expected["cases"].get(name), actual["cases"].get(name)
        if want == got:
            continue
        kind = "structure (the code changed)" if want and got and want[1] != got[1] else "numbers only (the platform?)"
        lines.append(f"case {name}: {kind}: expected {want and want[0][:12]}, got {got and got[0][:12]}")
    return lines


if __name__ == "__main__":   # pragma: no cover -- a quick look at the records
    print(json.dumps(written_battery(), indent=2), file=sys.stdout)
