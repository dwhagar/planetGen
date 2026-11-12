# planetgen/physics/disk_dispersal.py

"""
How long a young star keeps its disk (GEN.1551, built ahead by Lane 07 for Lane 11)
=====================================================================================

Each star gets its own disk dispersal time, drawn once from a stream of its own, and holds its class 0, I or II disk while it is younger than that time
(`reports/young-stars-disks-no-bodies-2026-11-11.md`, recommendation 2). The same time replaces the fixed 10 Myr planet cut (`PLANET_MIN_STAR_AGE_GY`, GEN.1552): bodies
appear once the star's own disk is gone.

* The time follows the survey curve of the share of stars that still hold a dusty disk, read as a survival curve: 1 at birth, `DISK_SURVIVAL_POINTS` at 1, 3 and 5.5 Myr
  and none by 10 Myr (stand-in DD1). One uniform draw is turned into a time by walking the curve backwards (straight lines between the points).
* A close companion inside `COMPANION_GAS_CUT_AU` (47 au, GEN.1475) disperses the gas early: two thirds of such stars lose the disk within 1 Myr (stand-in DD2). The test draws
  one more uniform for the companion case, in a fixed order, so a star with a companion and the same star without one draw the same first number.
* The draw is keyed on the star (its mass and a birth key the caller gives, such as its galactic position), never on the age text or the final age, so remaking the star
  at another age keeps its time (GEN.1547). A star's mass does not move the curve (the report holds no source for it; a mass slope is the first thing a later source adds).

Pure functions of their arguments and their stream; nothing in the old order calls them. Lane 11 wires `disk_dispersal_time_yr` into the star's disk stage (`pre_main_sequence.draw_disk`
reads `disk_present` instead of the fixed class ages) and the planet cut (GEN.1552), and deletes the wiring-register line.
"""

import math

from planetgen.util import draw

DISK_SURVIVAL_POINTS = ((0.0, 1.0), (1.0e6, 0.80), (3.0e6, 0.50), (5.5e6, 0.10), (10.0e6, 0.0))
"""tuple: `(age in years, share of stars still holding a disk)`, joined by straight lines. Stand-in DD1 (the table of the design note holds 0.8 at 1 Myr and 0.1 at 5 to 6 Myr;
0.5 at 3 Myr is the usual survey reading, Haisch, Lada and Lada 2001; Mamajek 2009 fits an e-folding of about 2.5 Myr), tier Confident Unverified, tweak L0."""

COMPANION_GAS_CUT_AU = 47.0
"""float: A companion inside this distance, au, cuts the gas of the disk early (GEN.1475)."""

COMPANION_EARLY_SHARE = 2.0 / 3.0
"""float: The share of stars with such a companion that lose their disk inside `COMPANION_EARLY_TIME_YR`. Stand-in DD2, tier Marginal Unverified, tweak L0."""

COMPANION_EARLY_TIME_YR = 1.0e6
"""float: The time inside which an early-cut disk is gone, years."""


def survival_to_time_yr(survival, points=DISK_SURVIVAL_POINTS):
    """The age, years, at which the share of disks still present falls to `survival` (0 to 1), by the lines between `points`; the inverse of `disk_share`."""
    survival = min(1.0, max(0.0, survival))
    for (t0, s0), (t1, s1) in zip(points, points[1:]):
        if s1 <= survival <= s0:
            return t0 if s0 == s1 else t0 + (s0 - survival) / (s0 - s1) * (t1 - t0)
    return points[-1][0]


def disk_share(age_yr, points=DISK_SURVIVAL_POINTS):
    """The share of stars of that age that still hold a disk (the survey curve itself): 1 at birth, 0 from the last point."""
    if age_yr <= points[0][0]:
        return points[0][1]
    for (t0, s0), (t1, s1) in zip(points, points[1:]):
        if age_yr <= t1:
            return s0 + (age_yr - t0) / (t1 - t0) * (s1 - s0)
    return points[-1][1]


def dispersal_stream(mass_sol, birth_key):
    """The star's own stream: its mass and the birth key the caller gives (its position, say), never its age."""
    return draw.Stream(f"disk_dispersal:{mass_sol!r}:{birth_key!r}")


def disk_dispersal_time_yr(mass_sol, birth_key, companion_au=None):
    """
    The age, years, at which this star's disk is gone. `companion_au` is the distance of the nearest companion (None for a single star): inside `COMPANION_GAS_CUT_AU`
    two thirds of the stars lose the disk inside a million years. The draws are made in a fixed order (the survey time, then the companion case) so a companion never moves
    the survey time.
    """
    stream = dispersal_stream(mass_sol, birth_key)
    survey_time = survival_to_time_yr(stream.random())
    early_roll, early_time = stream.random(), stream.random() * COMPANION_EARLY_TIME_YR
    if companion_au is not None and companion_au < COMPANION_GAS_CUT_AU and early_roll < COMPANION_EARLY_SHARE:
        return min(survey_time, early_time)
    return survey_time


def disk_present(age_yr, mass_sol, birth_key, companion_au=None):
    """Whether the star still holds its disk at `age_yr`: younger than its own dispersal time."""
    return age_yr < disk_dispersal_time_yr(mass_sol, birth_key, companion_au)
