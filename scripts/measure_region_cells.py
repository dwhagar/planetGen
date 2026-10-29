#!/usr/bin/env python3
"""
How many 9-blocks of the sector grid the grid squeezes (GEN.825): the cells whose inscribed radius
(`named_regions.inscribed_radius_pc`) is under the standard region radius (12 pc) or under 8 pc.

    python scripts/measure_region_cells.py [--edge-pc 4] [--max-shell 9227]

A cell's inscribed radius depends on its shell and its row only, so each (shell, row) is measured once and counted as many
times as the row has wedges. The block shells are measured one in `step` beyond the first `exact` ones, each standing for the
shells up to the next; the grid is scale-free there, so the share barely moves between the sampled shells. Prints a table by
block-shell range and the totals; the numbers go in the design draft.
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from planetgen import tuning  # noqa: E402
from planetgen.galaxy import drill, named_regions, spherical  # noqa: E402


def measure_shell(shell, edge_pc):
    """`(cells, under 12 pc, under 8 pc, smallest inscribed radius)` of the 9-blocks of block shell `shell`."""
    rows = drill.drill_band_count(named_regions.REGION_LEVEL, shell)
    cells = under_standard = under_eight = 0
    smallest = math.inf
    for band in range(-(rows // 2), rows // 2 + 1):
        wedges = drill.drill_wedge_count(named_regions.REGION_LEVEL, shell, band)
        radius = named_regions.inscribed_radius_pc(drill.DrillBlock(named_regions.REGION_LEVEL, shell, band, 0), edge_pc)
        cells += wedges
        smallest = min(smallest, radius)
        if radius < named_regions.REGION_RADIUS_PC:
            under_standard += wedges
        if radius < 8.0:
            under_eight += wedges
    return cells, under_standard, under_eight, smallest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--edge-pc", type=float, default=tuning.DEFAULT_SECTOR_EDGE_PC)
    parser.add_argument("--max-shell", type=int, default=spherical.MAX_ADDRESSABLE_SHELL // named_regions.REGION_LEVEL)
    parser.add_argument("--exact", type=int, default=40, help="the first block shells measured one by one")
    parser.add_argument("--step", type=float, default=1.08, help="the growth of the gap between sampled shells beyond them")
    args = parser.parse_args(argv)
    shells = list(range(min(args.exact, args.max_shell + 1)))
    gap = 1.0
    shell = shells[-1] + 1
    while shell <= args.max_shell:
        shells.append(shell)
        gap *= args.step
        shell += max(1, int(gap))
    totals = [0.0, 0.0, 0.0]
    smallest = math.inf
    print(f"{'block shells':>16} {'cells':>16} {'under 12 pc':>14} {'share':>8} {'under 8 pc':>12} {'share':>8}")
    ranges = [(0, 0), (1, 9), (10, 99), (100, 999), (1000, args.max_shell)]
    by_range = {r: [0.0, 0.0, 0.0] for r in ranges}
    for index, shell in enumerate(shells):
        weight = (shells[index + 1] if index + 1 < len(shells) else args.max_shell + 1) - shell
        cells, under_standard, under_eight, low = measure_shell(shell, args.edge_pc)
        smallest = min(smallest, low)
        for r in ranges:
            if r[0] <= shell <= r[1]:
                for i, v in enumerate((cells, under_standard, under_eight)):
                    by_range[r][i] += v * weight
        for i, v in enumerate((cells, under_standard, under_eight)):
            totals[i] += v * weight
    for r, (cells, a, b) in by_range.items():
        if cells:
            print(f"{r[0]:>7} to {r[1]:<6} {cells:16.0f} {a:14.0f} {a / cells:8.2%} {b:12.0f} {b / cells:8.2%}")
    print(f"{'all':>16} {totals[0]:16.0f} {totals[1]:14.0f} {totals[1] / totals[0]:8.2%} {totals[2]:12.0f} "
          f"{totals[2] / totals[0]:8.2%}")
    print(f"smallest inscribed radius seen: {smallest:.3f} pc ({len(shells)} block shells measured)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
