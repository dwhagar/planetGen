# planetgen/util/ladders.py

"""
Unit Ladders
============

The data behind every measurement the site shows (UX.23): for each quantity
a ladder of units, smallest first, and the rule for which one a value is
shown in. `planetgen.util.format.format_quantity` is the one algorithm that
reads these tables in Python; `scripts/build_unit_ladders.py` writes the
same tables to `html/static/unitladders.js` for `html/static/unitladder.js`,
the one algorithm in the browser. A test fails when the committed file is
out of date, so the two copies cannot drift apart
(docs/design/units-and-number-formatting.md).

Each ladder's base is the SI unit of its quantity, and a rung's `factor` is
how many base units one of it is. A value is shown in the largest rung it
is at least `start` of (a hair of tolerance, so a value computed as exactly
1 unit isn't shown as 1,000 of the one below it). A value that would show
under `TINY_NUMBER` of a rung is shown in scientific notation instead of a
row of zeros (UX.81).
"""

from dataclasses import dataclass
from typing import Optional, Tuple

from planetgen.physics import constants as physical_constants

TINY_NUMBER = 1e-3
"""float: A number below this in the unit shown is written as scientific
notation ("2.50 × 10⁻¹³ Pa"), not "0.00000000000025 Pa" (UX.81)."""


@dataclass(frozen=True)
class Rung:
    """One unit on a ladder."""

    symbol: str
    """str: Shown after the number ("AU", "km/s")."""

    factor: float
    """float: Base units in one of this unit."""

    name: str
    """str: Spelled out, for screen readers and tooltips ("astronomical units")."""

    start: Optional[float] = None
    """float, optional: Smallest size, in base units, shown in this rung; the
    `factor` when unset. Zero for the lowest rung of a ladder that has
    nothing smaller to fall back to."""

    singular: Optional[str] = None
    """str, optional: Shown instead of `symbol` when the number is exactly 1
    ("year", "minute")."""

    whole_from: Optional[float] = None
    """float, optional: From this size, in base units, the number is shown as
    a whole number, not to three significant figures ("384,400 km")."""

    with_aside: bool = False
    """bool: The text carries one more familiar unit in parentheses,
    chosen by the ladder's `asides` ("4.2 pc (13.7 ly)")."""


@dataclass(frozen=True)
class Companion:
    """A second unit shown beside the primary one, in parentheses ("288 K
    (15 °C, 59 °F)"), to help a reader judge the size of a value."""

    symbol: str
    """str: Shown after the number."""

    name: str
    """str: Spelled out, for screen readers."""

    factor: float
    """float: Base units in one of this unit."""

    offset: float = 0.0
    """float: The base value that is zero in this unit (273.15 K for °C), so
    the number shown is (value - offset) / factor."""

    whole: bool = False
    """bool: Shown to whole units, or one decimal under 10 in size ("15",
    "5.9"), not to three significant figures."""

    min_size: Optional[float] = None
    """float, optional: Only shown when the value's size, in base units, is
    at least this."""

    max_size: Optional[float] = None
    """float, optional: Only shown when the value's size, in base units, is
    under this."""


@dataclass(frozen=True)
class Ladder:
    """The units one quantity is shown in."""

    quantity: str
    """str: What is measured ("length", "speed", "time")."""

    base: str
    """str: The SI unit the values are in ("m", "m/s", "s")."""

    rungs: Tuple[Rung, ...]
    """Tuple[Rung, ...]: Smallest first."""

    nonfinite: str
    """str: The symbol after "inf" or "-inf" for an infinite value."""

    zero: Optional[str] = None
    """str, optional: The symbol zero is shown in; the lowest rung when unset."""

    companions: Tuple[Companion, ...] = ()
    """Tuple[Companion, ...]: Units shown beside the primary one, in this order."""

    asides: Tuple[Tuple[float, str], ...] = ()
    """Tuple[Tuple[float, str], ...]: For a rung `with_aside`, the first
    (smallest size in base units, rung symbol) that the value reaches, in
    order high to low; the last entry's size is 0 so one always matches."""

    def rung(self, symbol):
        """The rung with `symbol`."""
        for candidate in self.rungs:
            if candidate.symbol == symbol:
                return candidate
        raise KeyError(symbol)


LENGTH = Ladder(
    quantity="length",
    base="m",
    nonfinite="km",
    rungs=(
        # Whole kilometers read better than three figures ("384,400 km").
        Rung("km", physical_constants.KM_M, "kilometers", start=0.0, whole_from=1e6),
        # From a million km (0.0067 AU), so Mercury is "0.387 AU", not "5.79 × 10⁷ km" (UX.81).
        Rung("AU", physical_constants.AU_M, "astronomical units", start=1e6 * physical_constants.KM_M),
        Rung("mpc", physical_constants.MILLIPARSEC_M, "milliparsecs", with_aside=True),
        Rung("cpc", physical_constants.CENTIPARSEC_M, "centiparsecs", with_aside=True),
        Rung("ly", physical_constants.LIGHTYEAR_M, "light-years"),
        Rung("pc", physical_constants.PARSEC_M, "parsecs", with_aside=True),
        Rung("kpc", physical_constants.KILOPARSEC_M, "kiloparsecs", with_aside=True),
        Rung("Mpc", physical_constants.MEGAPARSEC_M, "megaparsecs", with_aside=True),
        Rung("Gpc", physical_constants.GIGAPARSEC_M, "gigaparsecs", with_aside=True),
    ),
    # A parsec value's parenthetical is in ly when the distance is at least
    # 0.01 ly, else AU when at least 0.01 AU, else km (Boss, 2026-09-30:
    # "4.2 pc (13.7 ly)", "2.4 mpc (495 AU)").
    asides=(
        (0.01 * physical_constants.LIGHTYEAR_M, "ly"),
        (0.01 * physical_constants.AU_M, "AU"),
        (0.0, "km"),
    ),
)

# Boss, 2026-10-01 (UX.13): "km/h on the low speed end to Mm/s on the high
# speed end", then multiples of c from a tenth of light speed.
SPEED = Ladder(
    quantity="speed",
    base="m/s",
    nonfinite="km/s",
    rungs=(
        Rung("km/h", physical_constants.KM_M / 3600, "kilometers per hour", start=0.0),
        Rung("km/s", physical_constants.KM_M, "kilometers per second"),
        Rung("Mm/s", 1e6, "megameters per second"),
        Rung("c", physical_constants.SPEED_OF_LIGHT_M_S, "times the speed of light",
             start=0.1 * physical_constants.SPEED_OF_LIGHT_M_S),
    ),
)

# A year is the Julian 365.25 days (UX.14).
TIME = Ladder(
    quantity="time",
    base="s",
    nonfinite="years",
    zero="s",
    rungs=(
        Rung("µs", 1e-6, "microseconds", start=0.0),
        Rung("ms", 1e-3, "milliseconds"),
        Rung("s", 1.0, "seconds"),
        Rung("minutes", 60.0, "minutes", singular="minute"),
        Rung("hours", 3600.0, "hours", singular="hour"),
        Rung("days", 86400.0, "days", singular="day"),
        Rung("years", physical_constants.SECONDS_PER_YEAR, "years", singular="year"),
        # Not "Gy": that is the gray in SI (UX.81).
        Rung("kyr", physical_constants.SECONDS_PER_YEAR * 1e3, "thousand years"),
        Rung("Myr", physical_constants.SECONDS_PER_YEAR * 1e6, "million years"),
        Rung("Gyr", physical_constants.SECONDS_PER_YEAR * 1e9, "billion years"),
    ),
)

# A surface temperature (Boss, 2026-10-01: "temperatures should be reported
# in K, C, and F for surface conditions"); a star's effective temperature
# is shown in K alone (`companions=False`). Under 1 K the usual SI prefixes.
TEMPERATURE = Ladder(
    quantity="temperature",
    base="K",
    nonfinite="K",
    zero="K",
    rungs=(
        Rung("nK", 1e-9, "nanokelvins", start=0.0),
        Rung("µK", 1e-6, "microkelvins"),
        Rung("mK", 1e-3, "millikelvins"),
        Rung("K", 1.0, "kelvins"),
    ),
    companions=(
        Companion("°C", "degrees Celsius", 1.0, offset=physical_constants.CELSIUS_ZERO_K, whole=True),
        Companion("°F", "degrees Fahrenheit", 5 / 9, offset=physical_constants.CELSIUS_ZERO_K - 32 * 5 / 9, whole=True),
    ),
)

# A pressure on the SI ladder with customary units beside it (Boss,
# 2026-10-01: "Atmospheric pressure and surface conditions should show
# customary units as well as a secondary to help contextualize the metric
# values given"): atm and psi always, bar from 10 atm up to 10,000 atm (Venus's
# 92 bar, a gas giant's depths), mbar from 1 Pa up to 0.01 atm (Mars's 6 mbar).
PRESSURE = Ladder(
    quantity="pressure",
    base="Pa",
    nonfinite="Pa",
    zero="Pa",
    rungs=(
        Rung("nPa", 1e-9, "nanopascals", start=0.0),
        Rung("µPa", 1e-6, "micropascals"),
        Rung("mPa", 1e-3, "millipascals"),
        Rung("Pa", 1.0, "pascals"),
        Rung("kPa", 1e3, "kilopascals"),
        Rung("MPa", 1e6, "megapascals"),
        Rung("GPa", 1e9, "gigapascals"),
        Rung("TPa", 1e12, "terapascals"),
    ),
    companions=(
        Companion("atm", "standard atmospheres", physical_constants.STANDARD_ATMOSPHERE_PA),
        Companion("psi", "pounds per square inch", physical_constants.PSI_PA),
        Companion("bar", "bars", 1e5, min_size=10 * physical_constants.STANDARD_ATMOSPHERE_PA,
                  max_size=1e4 * physical_constants.STANDARD_ATMOSPHERE_PA),
        Companion("mbar", "millibars", 1e2, min_size=1.0, max_size=0.01 * physical_constants.STANDARD_ATMOSPHERE_PA),
    ),
)

# Gravity: m/s² with standard gravity beside it ("3.71 m/s² (0.378 g)").
GRAVITY = Ladder(
    quantity="gravity",
    base="m/s²",
    nonfinite="m/s²",
    rungs=(Rung("m/s²", 1.0, "meters per second squared", start=0.0),),
    companions=(Companion("g", "times standard gravity", physical_constants.EARTH_GRAVITY),),
)

LADDERS = (LENGTH, SPEED, TIME, TEMPERATURE, PRESSURE, GRAVITY)
"""Tuple[Ladder, ...]: Every ladder, in the order `unitladders.js` lists them."""
