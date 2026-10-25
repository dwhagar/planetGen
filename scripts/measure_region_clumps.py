#!/usr/bin/env python3
"""
How often do the generated sectors of one region cell form two far-apart clumps (GEN.827)?

    python scripts/measure_region_clumps.py [--trials 400] [--seed 5]

A model, not the galaxy: the cell is 9 sectors of 4 pc on a side (the 36 pc of a 9-block, the grid taken as cubic), and k fills
of radius r (sectors whose centre is within r of a random point of the cell) are made in it. For each (k, r) the script prints
the share of cells whose sectors form two or more face-connected clumps, and the share where the bounding sphere of all the
sectors is more than 1.5 times the sphere of the largest clump (the gap the single sphere would cover). The bounding sphere is
the one of `named_regions.bounding_sphere` (centres plus a sector's circumradius). A standard fill at an anchor (GEN.826) makes
one clump; the cases here are the ad hoc fills a user makes by hand.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from planetgen import tuning  # noqa: E402
from planetgen.galaxy import named_regions  # noqa: E402

EDGE = float(tuning.DEFAULT_SECTOR_EDGE_PC)
BOX = 9 * EDGE
FACES = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))


def clumps(cells):
    """The face-connected groups of integer cells."""
    left, groups = set(cells), []
    while left:
        stack = [left.pop()]
        group = list(stack)
        while stack:
            x, y, z = stack.pop()
            for dx, dy, dz in FACES:
                other = (x + dx, y + dy, z + dz)
                if other in left:
                    left.remove(other)
                    stack.append(other)
                    group.append(other)
        groups.append(group)
    return groups


def fill(rng, radius):
    """The cells whose centre is within `radius` of a random point of the cell."""
    centre = rng.uniform(0, BOX, 3)
    reach = int(radius // EDGE) + 1
    base = np.floor(centre / EDGE).astype(int)
    out = []
    for dx in range(-reach, reach + 1):
        for dy in range(-reach, reach + 1):
            for dz in range(-reach, reach + 1):
                cell = base + np.array([dx, dy, dz])
                if np.linalg.norm((cell + 0.5) * EDGE - centre) <= radius:
                    out.append(tuple(int(v) for v in cell))
    return out


def sphere_radius(cells):
    points = (np.array(list(cells), float) + 0.5) * EDGE
    return named_regions.bounding_sphere(points, EDGE)[1]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--trials", type=int, default=400)
    parser.add_argument("--seed", type=int, default=5)
    args = parser.parse_args(argv)
    rng = np.random.default_rng(args.seed)
    print(f"{'fills':>5} {'radius pc':>10} {'2+ clumps':>10} {'sphere >1.5x clump':>19} {'median sectors':>15}")
    for k, radius in ((2, 12.0), (3, 12.0), (2, 8.0), (3, 8.0), (2, 16.0)):
        multi = wide = 0
        counts = []
        for _ in range(args.trials):
            cells = set()
            for _fill in range(k):
                cells |= set(fill(rng, radius))
            if not cells:
                continue
            groups = clumps(cells)
            counts.append(len(cells))
            if len(groups) > 1:
                multi += 1
                if sphere_radius(cells) > 1.5 * sphere_radius(max(groups, key=len)):
                    wide += 1
        print(f"{k:>5} {radius:>10.0f} {multi / args.trials:>10.1%} {wide / args.trials:>19.1%} {int(np.median(counts)):>15}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
