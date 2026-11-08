# planetgen/physics/mantle_reservoirs.py

"""
The water and carbon a rocky mantle stores (GEN.1365)
=====================================================

A mantle holds water in its minerals, and how much depends on which minerals the pressure makes. The data are
`physics/data/mantle_reservoirs.json`; three zones of a mantle, by pressure:

* **the upper mantle**, olivine, pyroxene and garnet, to 13 GPa (the olivine to wadsleyite change; 410 km on Earth): the upper mantle of Earth holds
  50 to 200 ppm of water [R5142: Abstract] and olivine can hold under 0.9 wt% [R5139: Introduction];
* **the transition zone**, wadsleyite and ringwoodite, from 13 to 23 GPa (410 to 660 km on Earth): up to about 3 wt% [R5139: Introduction], 1.5 to 2 wt%
  for ringwoodite [R5140: Abstract], 1.43 wt% in the diamond-hosted ringwoodite [R5144], about 1.0 wt% by several lines of evidence [R5141: Introduction];
* **the lower mantle**, bridgmanite and ferropericlase, from 23 GPa: at most 0.1 wt% [R5141: Introduction] and under 20 ppm [R5142: Abstract].

`zone_shares` puts the mantle's mass in those zones from the body's mass, radius and rock fraction: the mass between two pressures is
`4 pi r^2 dP / g`, the depth of each pressure comes from `bulk_structure.mantle_density` under the surface gravity, and the walk down stops where
the mantle's mass has run out (the core-mantle boundary). For the Earth it gives the core-mantle boundary at 132 GPa (136 in the PREM model) and a
transition zone of 10.9 percent of the mantle's mass. `water_capacity_fraction` is the sum of each zone's share times its capacity, as a mass fraction
of the body: what `hydrosphere.surface_water_fraction` takes as the most water the mantle keeps (the rest stays at the surface).

The capacity of each zone is one number, `CAPACITY_WT_PERCENT` (stand-ins MR1 to MR3): the sources print no dependence of the capacity on temperature,
so there is none here, and a hot or cold mantle of the same minerals stores the same.

The carbon of a mantle is in carbonates in an oxidized mantle and in diamond or graphite in a reduced one; `carbon_form` splits them at `CARBON_REDOX_DIW`,
about 3 log units over the iron-wuestite buffer, where diamond and carbonate meet at 15 to 45 GPa (stand-in MR4).

Use:

    zone_shares(mass_kg, radius_m, rock_fraction)    # (upper mantle, transition zone, lower mantle) shares of the mantle's mass, the core-mantle boundary Pa
    water_capacity_fraction(mass_kg, radius_m, rock_fraction)   # the mantle's capacity as a mass fraction of the body
    carbon_form(delta_iw)                            # "carbonate" or "reduced"
"""

import json
import math
import os
from collections import namedtuple

from planetgen.physics import bulk_structure, constants, elements
from planetgen.physics.elements import DataError

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "mantle_reservoirs.json")
"""str: The data file."""

ZONES = ("upper mantle", "transition zone", "lower mantle")
"""tuple: The three zones, from the surface down."""

PRESSURE_STEP_PA = 5.0e7
"""float: The pressure step of the walk down the mantle, Pa (0.05 GPa)."""

MANTLE_X_FE = 0.11
"""float: The Fe/(Fe+Mg) of the mantle the density is read at (`bulk_structure`'s default, the Earth's)."""

Zones = namedtuple("Zones", "shares cmb_pressure_pa")
"""The shares of the mantle's mass in `ZONES`, and the pressure where its mass runs out, Pa."""


def load(path=DATA_FILE):
    """Reads and checks the file at `path`: `(top of the transition zone Pa, bottom Pa, capacities (mass fractions) by zone, carbon redox dIW)`."""
    try:
        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, ValueError) as error:
        raise DataError(f"cannot read the mantle reservoir table {path}: {error}") from error

    def positive(value, what):
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value) or value <= 0:
            raise DataError(f"{what} must be a positive number, not {value!r}")
        return float(value)

    edges = document.get("boundaries_gpa") or {}
    top, bottom = positive(edges.get("transition_top"), "the top of the transition zone"), positive(edges.get("transition_bottom"), "its bottom")
    if not top < bottom:
        raise DataError("the transition zone must start above where it ends")
    raw = document.get("capacity_wt_percent") or {}
    capacity = tuple(positive(raw.get(name.replace(" ", "_")), f"the {name} capacity") / 100.0 for name in ZONES)
    if capacity[1] <= capacity[0] or capacity[1] <= capacity[2]:
        raise DataError("the transition zone holds the most water of the three")
    redox = document.get("carbon_redox_diw")
    if not isinstance(redox, (int, float)) or isinstance(redox, bool) or not math.isfinite(redox):
        raise DataError(f"the carbon redox limit must be a number, not {redox!r}")
    return top * 1e9, bottom * 1e9, capacity, float(redox)


TRANSITION_TOP_PA, TRANSITION_BOTTOM_PA, CAPACITY_FRACTION, CARBON_REDOX_DIW = load()
"""The pressures bounding the transition zone (Pa), the water each zone holds at most as a mass fraction of its rock, and the dIW above which carbon is carbonate."""

CAPACITY_WT_PERCENT = {name: fraction * 100.0 for name, fraction in zip(ZONES, CAPACITY_FRACTION)}
"""dict: The capacity of each zone, wt% water (stand-ins MR1 to MR3)."""


def zone_shares(mass_kg, radius_m, rock_fraction):
    """
    The share of a body's mantle mass in each zone of `ZONES`, and the pressure at the core-mantle boundary (Pa).

    The mantle's mass is `rock_fraction` of `mass_kg`; the walk goes down in steps of `PRESSURE_STEP_PA` from the surface, each step taking
    `4 pi r^2 dP / g` of the mantle with `g` the surface gravity and `r` the radius at the step's depth, and stops where the mantle's mass
    is used. A body whose mantle ends above the transition zone has all of it in the upper mantle.

    Raises:
        ValueError: A mass, radius or rock fraction that is not positive, or a rock fraction over 1.
    """
    if not (mass_kg > 0.0 and radius_m > 0.0 and 0.0 < rock_fraction <= 1.0):
        raise ValueError("the mass, the radius and the rock fraction must be positive and the fraction at most 1")
    gravity = constants.G * mass_kg / radius_m ** 2
    mantle = rock_fraction * mass_kg
    shares = [0.0, 0.0, 0.0]
    used, depth, pressure = 0.0, 0.0, 0.0
    while used < mantle:
        radius = radius_m - depth
        if radius <= 0.0:
            break
        mid = pressure + PRESSURE_STEP_PA / 2.0
        taken = min(4.0 * math.pi * radius * radius * PRESSURE_STEP_PA / gravity, mantle - used)
        shares[0 if mid < TRANSITION_TOP_PA else (1 if mid < TRANSITION_BOTTOM_PA else 2)] += taken
        used += taken
        depth += PRESSURE_STEP_PA / (bulk_structure.mantle_density(mid, MANTLE_X_FE) * gravity)
        pressure += PRESSURE_STEP_PA
    return Zones(tuple(share / mantle for share in shares), pressure)


def water_capacity_fraction(mass_kg, radius_m, rock_fraction):
    """The most water the mantle of a body holds, as a mass fraction of the whole body: the zones' shares of the mantle's mass times their capacities,
    times the mantle's share of the body."""
    zones = zone_shares(mass_kg, radius_m, rock_fraction)
    return rock_fraction * sum(share * capacity for share, capacity in zip(zones.shares, CAPACITY_FRACTION))


def carbon_form(delta_iw):
    """`"carbonate"` for a mantle whose oxygen fugacity is `delta_iw` log units over the iron-wuestite buffer or more above `CARBON_REDOX_DIW`, and
    `"reduced"` (diamond, graphite and iron carbide) below it. Raises `ValueError` for a value that is not a number."""
    if not isinstance(delta_iw, (int, float)) or isinstance(delta_iw, bool) or math.isnan(delta_iw):
        raise ValueError(f"the oxygen fugacity must be a number, not {delta_iw!r}")
    return "carbonate" if delta_iw > CARBON_REDOX_DIW else "reduced"
