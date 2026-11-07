# planetgen/generation/disk_ejection.py

"""
The bodies an outer disk's giant threw out, for the rogue builders (GEN.1308, docs/design/systems/belt-origins.md section 6)
==============================================================================================================================

The mass a giant cleared from an outer disk and scattered is shared among the places the scatter step names (`physics/swept_hand_on.py`); the share that is *ejected* leaves the system. Almost all of it is
small bodies, which are interstellar debris and not planets. The part of it that is heavy enough to be a rogue planet (`tuning.ROGUE_PLANET_MASS_BIN_EDGES_EARTH`, the lightest bin's lower edge) is counted from the
disk's size law (`belt_ledger.ICY_LAW`, a disk beyond the snow line is icy) and handed to the rogue builders as one record of kind `disk`: the same `planet_architecture.Loss` the thrown planets leave
(`StarSystem.ejected_losses`), with the giant as its scatterer, so a rogue of a disk is built, kicked and placed like any other (`ejection_events`, `rogue_yield`).

    losses = disk_losses(bodies, star)      # [planet_architecture.Loss] for each outer disk with an ejected share that holds a body of rogue mass; [] for none
"""

import math

from planetgen import tuning
from planetgen.physics import belt_ledger, planet_architecture, scatter, unit_registry as units
from planetgen.util import draw

DISK_EJECTION_TIME_YR = scatter.DISK_LIFETIME_YR
"""float: When the bodies of a disk are thrown out, years after the system formed: the lifetime of the disk the scatter step uses for the bodies a planet scatters (`scatter.DISK_LIFETIME_YR`, row EJ2)."""

_EARTH_KG = units.convert(1.0, "Mearth", "kg")
_ROGUE_MIN_KG = tuning.ROGUE_PLANET_MASS_BIN_EDGES_EARTH["terrestrial"][0] * _EARTH_KG


def _mean_over_minimum():
    """The mean mass of the bodies above a mass in the steep stretch of the icy law (cumulative diameter exponent `q`, so `q / 3` in mass): `a / (a - 1)` times the minimum."""
    a = belt_ledger.ICY_LAW.exponents[-1] / 3.0
    return a / (a - 1.0)


def rogue_bodies(ejected_kg, stream):
    """`(count, total_kg)` of the bodies of rogue mass among `ejected_kg` of a disk thrown out: the expected number above the lightest rogue (`belt_ledger.count_above`), rounded by a draw of `stream`, and their mass
    (the count times the mean mass above the minimum, never more than the ejected mass)."""
    if not ejected_kg > 0.0:
        return 0, 0.0
    expected = belt_ledger.count_above(belt_ledger.ICY_LAW, ejected_kg, _ROGUE_MIN_KG)
    count = int(math.floor(expected)) + (1 if stream.random() < expected - math.floor(expected) else 0)
    return count, min(count * _mean_over_minimum() * _ROGUE_MIN_KG, ejected_kg)


def disk_losses(bodies, star):
    """
    The `Loss` records (kind `disk`) of the outer disks among `bodies` (the belts of `star` that `belt_birth.hand_on` handed their mass on): one for each disk whose ejected share holds at least one body of rogue mass.
    """
    star_solar = units.convert(star.mass, "kg", "Msun")
    found = []
    for belt in bodies:
        swept = getattr(belt, "swept_kg", None)
        if getattr(belt, "body_type", None) != 'a' or not swept or belt.scatter_giant_kg is None or belt.scatter_giant_au is None:
            continue
        count, total_kg = rogue_bodies(swept["ejected"], draw.Stream(f"disk-ejection:{belt.distance!r}:{belt.lower_limit!r}:{star.mass!r}"))
        if not count:
            continue
        each_earth = total_kg / count / _EARTH_KG
        theta = scatter.safronov_number(units.convert(belt.scatter_giant_kg, "kg", "Mearth"), units.convert(belt.scatter_giant_radius_km, "km", "Rearth"), belt.scatter_giant_au, star_solar)
        found.append(planet_architecture.Loss("disk", count, total_kg, total_kg, belt.distance, belt.distance, DISK_EJECTION_TIME_YR, belt.scatter_giant_au, belt.scatter_giant_kg, theta,
                                              scatter.hill_radius_au(each_earth, belt.distance, star_solar), belt.distance >= star.snow_line_au, 1.0, 0.0, 0.0))
    return found
