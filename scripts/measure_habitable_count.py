#!/usr/bin/env python3
"""
The habitable bodies per system, measured on the current generator (GEN.1245).

    python scripts/measure_habitable_count.py [--systems 4000] [--seed 5000000] [--jobs 4] [--chunk 100] [--out FILE]

Generates `--systems` systems from the seeds `--seed` upward with the default settings (`StarSystem(SystemConfig())`) and counts
the bodies of a life class in each (`System.habitable_class_count`: planets and moons whose class is in
`tuning.HABITABLE_PLANET_CLASSES`). It writes the histogram to `src/planetgen/physics/data/habitable_count.json`, which
`generation.directives` reads for `HABITABLE_COUNT_PMF`: the forced fill of a `habitable` directive (GEN.180) and the "impossible"
check size themselves from it.

It saves after every chunk and carries on from where the file stops (same `--seed`), so a long run can be interrupted. Re-run it
(a fresh file, `--fresh`) when a generator epoch changes how many bodies are of a life class; the file keeps the epoch it was
measured under and `tests/test_habitable_count_table.py` shows the table to be consistent with itself.
"""

import argparse
import json
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from planetgen.generation.config import SystemConfig  # noqa: E402
from planetgen.generation.system import StarSystem  # noqa: E402
from planetgen.physics import elements  # noqa: E402
from planetgen.util import draw, generator_epoch  # noqa: E402

DEFAULT_FILE = os.path.join(elements.DATA_DIRECTORY, "habitable_count.json")


def count(seed):
    """The habitable bodies of the one system drawn from `seed`."""
    with draw.bound(seed):
        return StarSystem(system_config=SystemConfig()).habitable_class_count


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--systems", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=5_000_000)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--chunk", type=int, default=100)
    parser.add_argument("--out", default=DEFAULT_FILE)
    parser.add_argument("--fresh", action="store_true", help="Ignore a file already there.")
    args = parser.parse_args()

    epoch = generator_epoch.generator_epoch()
    state = {"counts": [], "systems": 0}
    if not args.fresh and os.path.exists(args.out):
        with open(args.out, encoding="utf-8") as handle:
            state = json.load(handle)
        if state.get("epoch") != epoch or state.get("seed") != args.seed:
            sys.exit(f"{args.out} was measured under epoch {state.get('epoch')} from seed {state.get('seed')}; use --fresh")
    counts = list(state["counts"])
    done = state["systems"]
    with Pool(args.jobs) as pool:
        while done < args.systems:
            seeds = range(args.seed + done, args.seed + min(args.systems, done + args.chunk))
            for found in pool.map(count, seeds, chunksize=2):
                counts.extend([0] * (found + 1 - len(counts)))
                counts[found] += 1
            done += len(seeds)
            with open(args.out, "w", encoding="utf-8") as handle:
                json.dump({"epoch": epoch, "seed": args.seed, "systems": done, "counts": counts}, handle)
                handle.write("\n")
            print(f"{done} systems, mean {sum(i * n for i, n in enumerate(counts)) / done:.3f}, none {counts[0] / done:.3f}", flush=True)


if __name__ == "__main__":
    main()
