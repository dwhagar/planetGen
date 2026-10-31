# planetgen/physics/trojans.py

"""
Trojans: bodies at the Lagrange points L4 and L5 of a pair (GEN.725 rows TJ1 to TJ5, docs/design/systems/lagrange-trojans.md, reports/lagrange-trojans-2026-10-30.md)
===========================================================================================================================================================

A secondary (a giant, a brown dwarf, a large moon) of mass `m2` orbits a host of mass `m1`. A third body `m3` can keep the orbit of the secondary
60 degrees ahead (L4) or behind (L5) when the pair is lopsided enough. L1, L2 and L3 are saddles: no natural body stays there, and this module has none.

* The rule: L4 and L5 are linearly stable when `m2 / (m1 + m2)` is under `MU_CRIT = (1 - sqrt(23 / 27)) / 2 = 0.0385` (Gascheau 1843, Routh 1875; the text
  of Schwarz et al. 2016). With a finite trojan the secondary plus the trojan must stay under that share of the three (`stable`, `max_trojan_mass`).
  Two stars never pass; a star and a companion under about 4 percent of it do; a Class J or a Class M world passes at a star-and-giant when the
  sum of the giant and the trojan stays under that share of the star (a Class M world always does beside a giant).
* The motion: a tadpole librates about the point; its slow period is `2 pi / s` with `s^2 = (n^2 / 2) (1 - sqrt(1 - 27 mu (1 - mu)))` (`libration_period`;
  Jupiter's 12.5 orbits, about 150 years). Past a libration amplitude of `TADPOLE_LIMIT_DEG` the orbit is a horseshoe (`motion`).
* The swarm: the scatter keeps a share of an outer disc's swept mass as `trojans` (`physics/scatter.py`); `swarm_split` puts it into the leading (L4) and
  trailing (L5) regions.
* The label the map shows at the centre of a point (`point_label`): the secondary's name and the point, such as `Jupiter-L4`, with the pair (host, secondary) in the text.

Every number is in `physics/data/trojans.json` with its source tier and tweak level. Nothing here sets a rate; it says what is allowed.

Use:

    stable(1.0, 0.001, 3e-6)               # True: a Sun, a Jupiter and an Earth-mass trojan
    max_trojan_mass(1.0, 0.001)            # 0.0395...: host masses in whatever unit
    libration_period(11.86, 9.5e-4)        # about 148 years
    swarm_split(trojan_mass, 0.58)         # (leading, trailing)
    point_label("Jupiter", "L4")           # "Jupiter-L4"
"""

import json
import math
import os

from planetgen.physics import elements


class TrojanError(ValueError):
    """A rule file that is refused."""


def _load():
    with open(os.path.join(elements.DATA_DIRECTORY, "trojans.json"), encoding="utf-8") as handle:
        document = json.load(handle)
    for name, entry in document["numbers"].items():
        if entry.get("tier") not in tuple(f"T{n}" for n in range(1, 10)) or entry.get("tweak") not in tuple(f"L{n}" for n in range(5)) \
                or not entry.get("basis") or "value" not in entry:
            raise TrojanError(f"{name}: every number needs a value, a tier T1..T9, a tweak L0..L4 and a basis, got {entry!r}")
    return document


DOCUMENT = _load()
NUMBERS = DOCUMENT["numbers"]

MU_CRIT = (1.0 - math.sqrt(23.0 / 27.0)) / 2.0
"""float: The largest m2 / (m1 + m2) for which L4 and L5 are linearly stable. [Stand-in TJ1, tier T1, tweak L0.]"""

TADPOLE_LIMIT_DEG = NUMBERS["tadpole_limit_deg"]["value"]
"""float: The libration amplitude (degrees) past which a tadpole is a horseshoe. [Stand-in TJ3, tier T2, tweak L0.]"""

LEADING_SHARE_RANGE = tuple(NUMBERS["leading_share_range"]["value"])
"""tuple: The range the share of a swarm in L4 is drawn from. [Stand-in TJ4, tier T4, tweak L4.]"""

ECC_MAX = NUMBERS["ecc_max"]["value"]
"""float: The largest orbital eccentricity a trojan may have. [Stand-in TJ5, tier T4, tweak L4.]"""

INCLINATION_MAX_DEG = NUMBERS["inclination_max_deg"]["value"]
"""float: The largest inclination to the secondary's orbit (degrees). [Stand-in TJ5, tier T4, tweak L4.]"""

POINTS = ("L4", "L5")
"""tuple: The only points that hold natural bodies."""


def mass_parameter(host_mass, secondary_mass):
    """`m2 / (m1 + m2)` of the pair, in any one mass unit."""
    if host_mass <= 0.0 or secondary_mass < 0.0:
        raise TrojanError(f"masses must be a positive host and a secondary of at least 0, got {host_mass!r} and {secondary_mass!r}")
    return secondary_mass / (host_mass + secondary_mass)


def stable(host_mass, secondary_mass, trojan_mass=0.0):
    """True when L4 and L5 of the pair hold a trojan: the secondary plus the trojan under `MU_CRIT` of the three masses (TJ1, TJ2)."""
    if trojan_mass < 0.0:
        raise TrojanError(f"a trojan's mass is at least 0, got {trojan_mass!r}")
    return mass_parameter(host_mass, secondary_mass) < MU_CRIT and (secondary_mass + trojan_mass) / (host_mass + secondary_mass + trojan_mass) < MU_CRIT


def max_trojan_mass(host_mass, secondary_mass):
    """The heaviest trojan the pair can hold (the unit of the arguments); 0 when the pair holds none."""
    room = MU_CRIT * host_mass / (1.0 - MU_CRIT) - secondary_mass
    return max(room, 0.0)


def libration_period(orbital_period, mu):
    """The slow libration period of a tadpole of small amplitude, in the unit of `orbital_period`; `None` for a pair at or past `MU_CRIT` (TJ3)."""
    if not 0.0 < mu < MU_CRIT:
        return None
    slow_squared = 0.5 * (1.0 - math.sqrt(1.0 - 27.0 * mu * (1.0 - mu)))
    return orbital_period / math.sqrt(slow_squared)


def motion(amplitude_deg, mu):
    """`"tadpole"` below `TADPOLE_LIMIT_DEG` of libration, `"horseshoe"` above it up to 180 degrees; `None` for a pair that holds nothing (TJ3)."""
    if not 0.0 < mu < MU_CRIT or not 0.0 <= amplitude_deg <= 180.0:
        return None
    return "tadpole" if amplitude_deg < TADPOLE_LIMIT_DEG else "horseshoe"


def longitude_offset_deg(point):
    """Degrees of mean longitude from the secondary: +60 for L4 (ahead), -60 for L5 (behind)."""
    if point not in POINTS:
        raise TrojanError(f"only {POINTS} hold natural bodies, got {point!r}")
    return 60.0 if point == "L4" else -60.0


def orbit_allowed(eccentricity, inclination_deg):
    """True when a trojan's own orbit is inside `ECC_MAX` and `INCLINATION_MAX_DEG` of the secondary's plane (TJ5)."""
    return 0.0 <= eccentricity <= ECC_MAX and abs(inclination_deg) <= INCLINATION_MAX_DEG


def allows(host_mass, secondary_mass, trojan_mass, eccentricity=0.0, inclination_deg=0.0):
    """A body of `trojan_mass` may sit at L4 or L5 of this pair: the mass rule and the orbit limits both hold."""
    return stable(host_mass, secondary_mass, trojan_mass) and orbit_allowed(eccentricity, inclination_deg)


def leading_share(u):
    """The share of a swarm in L4 for `u` in 0..1 drawn by the caller (TJ4)."""
    low, high = LEADING_SHARE_RANGE
    return low + (high - low) * u


def swarm_split(trojan_mass, leading):
    """`(L4 mass, L5 mass)` of a swarm of `trojan_mass` with `leading` of it in L4."""
    if not 0.0 <= leading <= 1.0 or trojan_mass < 0.0:
        raise TrojanError(f"a share from 0 to 1 and a mass of at least 0 are needed, got {leading!r} and {trojan_mass!r}")
    return trojan_mass * leading, trojan_mass * (1.0 - leading)


def point_label(secondary_name, point):
    """The map label at the centre of a point: `"Jupiter-L4"` (the secondary's name, then the point)."""
    longitude_offset_deg(point)
    return f"{secondary_name}-{point}"


def pair_text(host_name, secondary_name, point):
    """The full words for a label: `"L4 of Sun and Jupiter"`."""
    longitude_offset_deg(point)
    return f"{point} of {host_name} and {secondary_name}"
