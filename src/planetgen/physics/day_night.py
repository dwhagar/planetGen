# planetgen/physics/day_night.py

"""
The day-and-night state of an object, and its counter (GEN.1097; Boss, 4:47 PM EDT 2026-10-30)
=========================================================================================

Every object has a day and night, or says plainly why it has none. This module is the pure part: from the spin, the orbit and the
polar circle that `daylight` derives it names the state and counts the days in a year. It draws nothing and has no stand-in but two
thresholds (`NON_ROTATING_TURNS`, `POLAR_CIRCLE_MAX_DEG`, registry row DN1).

The state is one word, the first that holds:

1. **`no night`**: a star (it shines, it does not stand under one).
2. **`no solid surface`**: a giant. Its cloud tops still turn under the star (the solar day is stored), but there is no ground to stand on.
3. **`tidally locked`**: one turn per orbit round the star (`spin.spin_state`): a permanent day side, a terminator and a permanent night side.
4. **`spin-orbit resonance`**: 3/2, 2 or 5/2 turns per orbit (Mercury: the star returns every two years).
5. **`non-rotating`**: under `NON_ROTATING_TURNS` of a turn per orbit, so the star goes round once a year: the solar day is the year.
6. **`polar day and night`**: the polar circle is within `POLAR_CIRCLE_MAX_DEG` of the equator (a tilt of 60 degrees or more), so half the surface or more has
   a day and a night as long as a season (Uranus at 98 degrees).
7. **`day and night`**: the rest, Earth's kind.

The counter is `days_per_year`: the year over the solar day (365.24 for Earth; 0.5 for Mercury; 1 for a non-rotating world); `None` when there
is no solar day (locked, a star, a giant with no period). `night_hours` is the solar day less the hours of daylight at the equator.
"""

from planetgen.physics import spin

NO_NIGHT = "no night"
NO_SURFACE = "no solid surface"
LOCKED = "tidally locked"
RESONANT = "spin-orbit resonance"
NON_ROTATING = "non-rotating"
POLAR = "polar day and night"
CYCLE = "day and night"
STATES = (NO_NIGHT, NO_SURFACE, LOCKED, RESONANT, NON_ROTATING, POLAR, CYCLE)

NON_ROTATING_TURNS = 0.01
"""float: A body that turns fewer times than this in one orbit is taken as not rotating (the star goes round once a year to within a percent). [Stand-in DN1, tier T9, tweak L0.]"""

POLAR_CIRCLE_MAX_DEG = 30.0
"""float: A body whose polar circle is within this latitude has half its surface or more in the polar zone (the area beyond latitude L is `1 - sin L` of a sphere). [Stand-in DN1, tier T9, tweak L0.]"""


def state(rotation_period_hours, period_years, polar_circle_deg=None, solid=True, star=False):
    """
    The state word of an object that spins `rotation_period_hours` and orbits its star in `period_years` (the parent's, for a moon), with the latitude of its
    `polar_circle_deg`. `solid` is False for a giant; `star` is True for an object that is a star. `None` for a planet or moon whose spin is not known.
    """
    if star:
        return NO_NIGHT
    if not solid:
        return NO_SURFACE
    if not rotation_period_hours or not period_years:
        return None
    kind, turns = spin.spin_state(rotation_period_hours, period_years)
    if kind == "locked":
        return LOCKED
    if kind == "resonant":
        return RESONANT
    if abs(turns) < NON_ROTATING_TURNS:
        return NON_ROTATING
    if polar_circle_deg is not None and polar_circle_deg <= POLAR_CIRCLE_MAX_DEG:
        return POLAR
    return CYCLE


def days_per_year(period_hours, solar_day_hours):
    """How many solar days one orbit holds; `None` where there is no solar day."""
    if not period_hours or not solar_day_hours:
        return None
    return period_hours / solar_day_hours


def night_hours(solar_day_hours, daylight_hours):
    """The hours of night in a solar day at the equator; `None` where either is missing."""
    if solar_day_hours is None or daylight_hours is None:
        return None
    return solar_day_hours - daylight_hours
