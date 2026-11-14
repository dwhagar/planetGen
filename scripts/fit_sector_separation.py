#!/usr/bin/env python3
"""
GEN.1534: fits the share of the galactic Hill spheres that random placement keeps between two systems (`tuning.SECTOR_HILL_SEPARATION_SHARE`, stand-in SEP1) to the spacing of the real
stars, the Gaia DR3 sample within 25 pc (`planetgen.datasets.gaia_nearby`).

The statistic is the distance of a system to its nearest neighbour (stars within 0.1 pc are one system; only systems within 15 pc ask, so the neighbour is always inside the complete
sample). The model is the galaxy's own placement rule (`galaxy/sector.py`): systems are put at random positions in a box at the density the sample gives, each rejected when it falls
nearer to one already placed than `share * (h_new + h_old)`, with Hill radii `h` log-normal about the median of the audit galaxy at the solar circle (1.3 pc, scatter `HILL_SIGMA`
in the natural log). The fit is the Kolmogorov-Smirnov distance between the model's nearest-neighbour distribution and the sample's, for each share; the shares inside the 95 percent
band of the sample's size are the ones the data allow.

    python scripts/fit_sector_separation.py [--box 70] [--seed 7]
"""

import argparse
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

HILL_MEDIAN_PC = 1.3
"""float: The median galactic Hill radius of a system at the solar circle, parsecs (the audit galaxy of GEN.1284, `galaxy/sector.py`)."""

HILL_SIGMA = 0.2
"""float: The scatter of the Hill radius in the natural logarithm; it follows the star mass to the one third (`Star.calculate_system_perimeter`), so 0.2 is a spread of masses of about a factor 2."""

SHARES = (0.0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5)
"""tuple: The shares the fit tries."""

GRID_PC = np.linspace(0.0, 3.5, 36)
"""The distances at which the two distributions are compared, parsecs."""


def simulate(share, density_per_pc3, box_pc=70.0, seed=7, median_pc=HILL_MEDIAN_PC, sigma=HILL_SIGMA):
    """
    The nearest-neighbour distances (parsecs, sorted) of the systems in the middle of a box of side `box_pc` filled by the galaxy's rejection placement at `density_per_pc3`
    (`share` of the sum of the two Hill radii kept between two systems). Also the share of the wanted systems that found room.
    """
    from scipy.spatial import cKDTree

    rng = np.random.default_rng(seed)
    wanted = int(density_per_pc3 * box_pc ** 3)
    points, hill, placed = np.empty((wanted, 3)), np.empty(wanted), 0
    attempts = 0
    while placed < wanted and attempts < wanted * 50:
        attempts += 1
        candidate, radius = rng.uniform(0.0, box_pc, 3), median_pc * math.exp(rng.normal(0.0, sigma))
        if placed and (np.linalg.norm(points[:placed] - candidate, axis=1) < share * (radius + hill[:placed])).any():
            continue
        points[placed], hill[placed] = candidate, radius
        placed += 1
    points = points[:placed]
    middle = np.all((points > 0.2 * box_pc) & (points < 0.8 * box_pc), axis=1)
    distances, _ = cKDTree(points).query(points[middle], k=2)
    return np.sort(distances[:, 1]), placed / wanted


def cdf(distances, grid=GRID_PC):
    """The share of `distances` (sorted array) under each distance of `grid`."""
    return np.searchsorted(distances, grid, side="left") / len(distances)


def ks_distance(model, sample):
    """The largest gap between the two distributions' cumulative shares over `GRID_PC`."""
    return float(np.max(np.abs(cdf(model) - cdf(sample))))


def critical_ks(count):
    """The Kolmogorov-Smirnov distance past which a model differs from a sample of `count` at the 95 percent level."""
    return 1.36 / math.sqrt(count)


def fit(shares=SHARES, box_pc=70.0, seed=7):
    """`[(share, ks distance, share of systems placed)]` against the Gaia sample, and the sample's size and density."""
    from planetgen.datasets import gaia_nearby

    centers = gaia_nearby.systems()
    sample = gaia_nearby.nearest_neighbour_distances_pc(centers)
    density = gaia_nearby.system_density_per_pc3(centers)
    rows = []
    for share in shares:
        model, kept = simulate(share, density, box_pc, seed)
        rows.append((share, ks_distance(model, sample), kept))
    return rows, len(sample), density


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--box", type=float, default=70.0)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args(argv)
    rows, count, density = fit(box_pc=args.box, seed=args.seed)
    print(f"Gaia sample: {count} systems within 15 pc ask, {density:.4f} systems per pc3 within 25 pc; the 95 percent KS band is {critical_ks(count):.3f}")
    print("share  KS distance  placed")
    for share, ks, kept in rows:
        print(f"{share:5.3f}  {ks:11.3f}  {kept:6.2f}  {'inside the band' if ks <= critical_ks(count) else ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
