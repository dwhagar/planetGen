"""
TEST.700: the candidate numbers `tests/data/civilized_pins.json` pins, so the tests that need a system with a technological civilization do not hunt for one.

    python -m tests.civilized_pins refresh            keep the numbers the last test run found (the temp file) in the pins file
    python -m tests.civilized_pins refresh --hunt 4   make candidates from 0 up until 6 hold a civilization, and pin those

Run it after a change that alters generation (the same moments as `python -m tests.golden_update`): a pin that no longer holds a civilization is only a wasted
candidate, but `test_the_pinned_civilized_seeds_still_qualify` fails to say the file is old.
"""

import argparse
import json
import pathlib
import sys

from tests import test_population as population

COUNT = 4
"""int: How many numbers the file keeps; the first two are the ones a test takes, the rest cover a number that stops qualifying."""


def _write(seeds):
    path = pathlib.Path(__file__).parent / "data" / population.PINS_FILE
    path.write_text(json.dumps({"note": "TEST.700: candidate numbers of test_population._candidate that hold a technological civilization; written by "
                                        "python -m tests.civilized_pins refresh", "seeds": sorted(seeds)[:COUNT]}, indent=2) + "\n")
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("command", choices=["refresh"])
    parser.add_argument("--hunt", type=int, default=0, help="make candidates until this many hold a civilization")
    args = parser.parse_args(argv)
    if args.hunt:
        found, number = [], 0
        while len(found) < min(args.hunt, COUNT):
            system, _cfg = population._candidate(number)
            print(f"tried {number}", flush=True)
            if population._has_civilization(system):
                found.append(number)
                print(f"candidate {number} holds a civilization ({len(found)} of {args.hunt})", flush=True)
            number += 1
    else:
        try:
            found = [int(number) for number in json.loads(population._seed_file().read_text())["seeds"]]
        except (OSError, ValueError, KeyError, TypeError):
            print("no test run has found a civilization with this generator yet: use --hunt", file=sys.stderr)
            return 1
    print(f"wrote {_write(found)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
