"""
TEST.700: the candidate numbers `tests/data/civilized_pins.json` pins, so the tests that need a system with a technological civilization do not hunt for one.

    python -m tests.civilized_pins refresh            keep the numbers the last test run found (the temp file) in the pins file
    python -m tests.civilized_pins refresh --hunt 4   make candidates from 0 up until 4 hold a civilization that gets a species, and pin those
    python -m tests.civilized_pins refresh --hunt 4 --from 200   the same, from candidate 200 (run several at once on different ranges, then pin by hand)

Run it after a change that alters generation (the same moments as `python -m tests.golden_update`): a pin that no longer holds a civilization is only a wasted
candidate, but `test_the_pinned_civilized_seeds_still_qualify` fails to say the file is old.
"""

import argparse
import contextlib
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


@contextlib.contextmanager
def _database():
    """A test database, as the `mysql_config` fixture leases one: the check the galaxy fixture makes of a candidate saves it."""
    from tests import dbpool, publicids
    from tests.conftest import _test_server_kwargs
    worker = dbpool.WorkerDatabase(_test_server_kwargs())
    config = worker.lease()
    publicids.use(config)
    try:
        yield config
    finally:
        publicids.use(None)
        worker.release()
        worker.close()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("command", choices=["refresh"])
    parser.add_argument("--hunt", type=int, default=0, help="make candidates until this many hold a civilization that gets a species")
    parser.add_argument("--from", dest="first", type=int, default=0, help="the candidate to start the hunt from")
    args = parser.parse_args(argv)
    if args.hunt:
        found, number = [], args.first
        with _database() as mysql_config:
            while len(found) < min(args.hunt, COUNT):
                system, cfg = population._candidate(number)
                civilized = population._has_civilization(system)
                qualifies = civilized and population.qualifies(mysql_config, system, cfg) is not None
                print(f"tried {number}" + (" (a civilization, but no species)" if civilized and not qualifies else ""), flush=True)
                if qualifies:
                    found.append(number)
                    print(f"candidate {number} holds a civilization and gets a species ({len(found)} of {args.hunt})", flush=True)
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
