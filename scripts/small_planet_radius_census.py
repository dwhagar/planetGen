"""GEN.1845: the radii of the small planets the generator makes, for comparison with the deconvolved CKS occurrence ([R1478: Table 4]).

Makes `--systems` G2V systems in blocks (each on a run seed of its own, as `tests/test_radius_valley_sweep.py` does) and writes one JSON line per block of
10 systems (its rocky planets: days, radius, mass, envelope fraction, density) to `--out`, resumable by block.

    PYTHONPATH=src python3 scripts/small_planet_radius_census.py --out /tmp/x.jsonl --systems 1200 --workers 4
"""
import argparse
import concurrent.futures
import json
import multiprocessing
import os

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants
from planetgen.util import draw

EARTH = constants.EARTH_MASS_TO_KG
PER_BLOCK = 10
SEED = 20261125


def block_rows(block):
    draw.set_run_seed(f"{SEED}:{block}")
    rows = []
    for _ in range(PER_BLOCK):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        system = StarSystem(system_config=cfg)
        for planet in system.all_planets:
            if getattr(planet, "body_type", None) != "t" or planet.is_moon:
                continue
            rows.append({"days": planet.period * 365.25, "radius": planet.radius / constants.EARTH_RADIUS_KM,
                         "mass": planet.mass / EARTH, "envelope": planet.envelope_mass_fraction or 0.0, "density": planet.density})
    return block, {"block": block, "systems": PER_BLOCK, "planets": rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--systems", type=int, default=1200)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--median-start", type=float, default=None, help="experiment: override envelope_capture.MEDIAN_START")
    args = parser.parse_args()
    blocks = args.systems // PER_BLOCK
    done = set()
    if os.path.exists(args.out):
        with open(args.out) as handle:
            for line in handle:
                try:
                    row = json.loads(line)
                except ValueError:
                    continue
                done.add(row["block"])
    todo = [b for b in range(blocks) if b not in done]
    if args.median_start is not None:
        from planetgen.physics import envelope_capture
        envelope_capture.MEDIAN_START = args.median_start
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context("fork")) as pool, open(args.out, "a") as out:
        for block, row in pool.map(block_rows, todo):
            out.write(json.dumps(row) + "\n")        # one line a block: a block of an interrupted run is redone whole
            out.flush()
            print("block", block, "done", flush=True)


if __name__ == "__main__":
    main()
