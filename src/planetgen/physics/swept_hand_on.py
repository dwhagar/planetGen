# planetgen/physics/swept_hand_on.py

"""
Where the mass an outer disk lost went (GEN.1308, docs/design/systems/belt-origins.md section 6)
===================================================================================================

An outer disk that a giant acted on has lost mass in the ledger of `belt_origin`: what the giant cleared around its orbit (`cleared_earth`) and what its move through the disk scattered
(`scattered_earth`), the *gone* mass of the disk. The scatter step (`scatter.scatter_shares`, GEN.1089) says where swept mass goes, by two rules that do not depend on where the swept band lies: the share
thrown out is `1 - exp(-Theta / THETA_EJECT)` of the swept mass (`Theta` the giant's Safronov number) and the rest is split by `SWEPT_SPLIT`. This module applies those two rules to the gone mass, so the disk's
own mass is never changed again (a second subtraction would count the same sweep twice) and the shares of the disk that stay in it (`scatter_shares`) are left as they are.

    parts = hand_on(gone_kg, theta)    # {"ejected", "implanted", "captured", "other"}, kilograms, adding up to gone_kg

* `ejected`: thrown out of the system; the rogue builders take the bodies of it that are planets (`generation/disk_ejection.py`), the rest is the small interstellar debris, counted and not drawn.
* `implanted`: the share put on orbits of an inner belt (`scatter.SWEPT_SPLIT`, `implanted_inner_belt`); the nearest inner belt of the star receives it (`generation/belt_birth.hand_on`).
* `captured`: the share the planets caught as irregular moons (`captured_moons`); the irregular moons of the star's planets are caught from it first (`irregular_moons.take_from_disks`).
* `other`: the scattered disk, the centaurs, the detached bodies and the trojans: they stay in the star's system as small bodies that no object here holds.

A disk with no giant has no Safronov number, and the gone mass has no places to go (`hand_on` returns `None`).
"""

import math

from planetgen.physics import scatter

_SPLIT = dict(scatter.SWEPT_SPLIT)


def ejected_share(theta):
    """The share of swept mass a giant of Safronov number `theta` throws out of the system (the rule of `scatter.scatter_shares`)."""
    return 1.0 - math.exp(-max(theta, 0.0) / scatter.THETA_EJECT)


def hand_on(gone_kg, theta):
    """
    The `gone_kg` of a disk (what was cleared from it and scattered) shared among `ejected`, `implanted`, `captured` and `other`, kilograms, adding up to `gone_kg`: the ejected share by `ejected_share(theta)`,
    the rest by `scatter.SWEPT_SPLIT`. `None` for a sweep of nothing or a giant with no Safronov number.
    """
    if not gone_kg > 0.0 or theta is None or not theta > 0.0:
        return None
    ejected = gone_kg * ejected_share(theta)
    bound = gone_kg - ejected
    parts = {"ejected": ejected, "implanted": bound * _SPLIT["implanted_inner_belt"], "captured": bound * _SPLIT["captured_moons"]}
    parts["other"] = gone_kg - sum(parts.values())
    return parts
