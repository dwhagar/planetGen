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
(docs/design/architecture/units-and-number-formatting.md).

Each ladder's base is the SI unit of its quantity, and a rung's `factor` is
how many base units one of it is. A rung or companion is built from its row of the
unit registry (`rung`, `companion`; GEN.568): the factor, offset and name are the
registry's, stated once, and what the ladder adds is display data (the symbol where it
shows another, `start`, `singular`, `systems`). A value is shown in the largest rung it
is at least `start` of (a hair of tolerance, so a value computed as exactly
1 unit isn't shown as 1,000 of the one below it). A value that would show
under `TINY_NUMBER` of a rung is shown in scientific notation instead of a
row of zeros (UX.81).
"""

import functools
from dataclasses import dataclass
from typing import Optional, Tuple

from planetgen.physics import constants as physical_constants
from planetgen.physics import unit_registry
from planetgen.physics import unit_registry as units

_mearth_to_kg = units.converter("Mearth", "kg")
_mjup_to_kg = units.converter("Mjup", "kg")
_msun_to_kg = units.converter("Msun", "kg")

TINY_NUMBER = 1e-3
"""float: A number below this in the unit shown is written as scientific
notation ("2.50 × 10⁻¹³ Pa"), not "0.00000000000025 Pa" (UX.81)."""


AUTOMATIC, METRIC, CUSTOMARY = "auto", "metric", "customary"
SYSTEMS = (AUTOMATIC, METRIC, CUSTOMARY)
"""tuple: The unit preferences (UX.78): Automatic is the ladders as Boss set
them, Metric only drops the astronomical and customary units for SI ones, and
Customary puts feet, miles, pounds, °F and psi first."""


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

    offset: float = 0.0
    """float: The base value that is zero in this unit (255.37 K for °F), so
    the number shown is (value - offset) / factor."""

    systems: Tuple[str, ...] = SYSTEMS
    """Tuple[str, ...]: The unit preferences this rung belongs to."""


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

    systems: Tuple[str, ...] = SYSTEMS
    """Tuple[str, ...]: The unit preferences this companion belongs to."""


def rung(key, symbol=None, name=None, **kwargs):
    """
    The rung for the registry unit `key` (`planetgen.physics.unit_registry`): its factor, offset and name come from the
    registry, and its symbol too unless the ladder shows another ("lunar masses" for M☾). `kwargs` are the display data
    (`start`, `singular`, `whole_from`, `with_aside`, `systems`).
    """
    row = unit_registry.unit(key)
    return Rung(row.symbol if symbol is None else symbol, row.factor, row.name if name is None else name,
                offset=row.offset, **kwargs)


def companion(key, symbol=None, name=None, **kwargs):
    """The companion for the registry unit `key`, as `rung` takes a rung from it; `kwargs` are `whole`, `min_size`,
    `max_size` and `systems`."""
    row = unit_registry.unit(key)
    return Companion(row.symbol if symbol is None else symbol, row.name if name is None else name, row.factor,
                     offset=row.offset, **kwargs)


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

    @functools.lru_cache(maxsize=None)
    def for_system(self, system):
        """
        This ladder as it reads under a unit preference (UX.78): only the
        rungs and companions of `system`, the asides that still have a rung
        (the last one taking over as the fallback), and the zero rung only
        if it is still there.
        """
        if system not in SYSTEMS:
            raise ValueError(f"unknown unit system {system!r}")
        rungs = tuple(rung for rung in self.rungs if system in rung.systems)
        symbols = {rung.symbol for rung in rungs}
        asides = tuple((size, symbol) for size, symbol in self.asides if symbol in symbols)
        if not any(rung.with_aside for rung in rungs):
            asides = ()
        if asides and asides[-1][0] != 0.0:
            asides = asides[:-1] + ((0.0, asides[-1][1]),)
        return Ladder(
            quantity=self.quantity, base=self.base, rungs=rungs, nonfinite=self.nonfinite,
            zero=self.zero if self.zero in symbols else None,
            companions=tuple(companion for companion in self.companions if system in companion.systems),
            asides=asides)


LENGTH = Ladder(
    quantity="length",
    base="m",
    nonfinite="km",
    zero="km",
    rungs=(
        # Customary: inches to a foot, feet to a mile, miles to the astronomical units.
        rung("in", start=0.0, systems=(CUSTOMARY,)),
        rung("ft", systems=(CUSTOMARY,)),
        rung("mi", whole_from=1e6, systems=(CUSTOMARY,)),
        # A thin ice shell reads "1.09 m", not "0.00109 km" (GEN.769): millimeters, centimeters, meters, kilometers.
        rung("mm", start=0.0, systems=(AUTOMATIC, METRIC)),
        rung("cm", systems=(AUTOMATIC, METRIC)),
        rung("m", systems=(AUTOMATIC, METRIC)),
        # Whole kilometers read better than three figures ("384,400 km").
        rung("km", whole_from=1e6, systems=(AUTOMATIC, METRIC)),
        # Metric only: SI prefixes in place of the astronomical units.
        rung("Gm", systems=(METRIC,)),
        rung("Tm", systems=(METRIC,)),
        rung("Pm", systems=(METRIC,)),
        rung("Em", systems=(METRIC,)),
        rung("Zm", systems=(METRIC,)),
        rung("Ym", systems=(METRIC,)),
        # From a million km (0.0067 AU), so Mercury is "0.387 AU", not "5.79 × 10⁷ km" (UX.81).
        rung("au", start=1e6 * physical_constants.KM_M, systems=(AUTOMATIC, CUSTOMARY)),
        rung("mpc", with_aside=True, systems=(AUTOMATIC, CUSTOMARY)),
        rung("cpc", with_aside=True, systems=(AUTOMATIC, CUSTOMARY)),
        rung("ly", systems=(AUTOMATIC, CUSTOMARY)),
        rung("pc", with_aside=True, systems=(AUTOMATIC, CUSTOMARY)),
        rung("kpc", with_aside=True, systems=(AUTOMATIC, CUSTOMARY)),
        rung("Mpc", with_aside=True, systems=(AUTOMATIC, CUSTOMARY)),
        rung("Gpc", with_aside=True, systems=(AUTOMATIC, CUSTOMARY)),
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
        rung("mph", start=0.0, systems=(CUSTOMARY,)),
        rung("km_h", start=0.0, systems=(AUTOMATIC, METRIC)),
        rung("mi_s", systems=(CUSTOMARY,)),
        rung("km_s", systems=(AUTOMATIC, METRIC)),
        rung("Mm_s", systems=(AUTOMATIC, METRIC)),
        rung("c", start=0.1 * physical_constants.SPEED_OF_LIGHT_M_S),
    ),
)

# A year is the Julian 365.25 days (UX.14).
TIME = Ladder(
    quantity="time",
    base="s",
    nonfinite="years",
    zero="s",
    rungs=(
        rung("us", start=0.0),
        rung("ms"),
        rung("s"),
        rung("min", symbol="minutes", singular="minute"),
        rung("h", symbol="hours", singular="hour"),
        rung("d", symbol="days", singular="day"),
        rung("yr", symbol="years", singular="year"),
        # Not "Gy": that is the gray in SI (UX.81).
        rung("kyr"),
        rung("Myr"),
        rung("Gyr"),
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
        rung("nK", start=0.0, systems=(AUTOMATIC, METRIC)),
        rung("uK", systems=(AUTOMATIC, METRIC)),
        rung("mK", systems=(AUTOMATIC, METRIC)),
        rung("K", systems=(AUTOMATIC, METRIC)),
        # Customary: °F all the way down, with °C and K beside it.
        rung("degF", start=0.0, systems=(CUSTOMARY,)),
    ),
    companions=(
        companion("degC", whole=True),
        companion("degF", whole=True, systems=(AUTOMATIC,)),
        companion("K", systems=(CUSTOMARY,)),
    ),
)

# A temperature difference (GEN.770 follow-up, Boss 9:13 PM EDT 2026-10-24: the day and night swing of 9.07 K read "-264 °C"):
# a size, not a place on a scale, so it has no offsets: 9.07 K is 9.07 °C and 16.3 °F.
TEMPERATURE_DIFFERENCE = Ladder(
    quantity="temperature difference",
    base="dK",
    nonfinite="K",
    zero="dK",
    rungs=(
        rung("dK", start=0.0, systems=(AUTOMATIC, METRIC)),
        rung("ddegF", start=0.0, systems=(CUSTOMARY,)),
    ),
    companions=(
        companion("ddegC"),
        companion("ddegF", systems=(AUTOMATIC,)),
        companion("dK", systems=(CUSTOMARY,)),
    ),
)

# A share of a whole (GEN.770): "water share of the mass 400 ppm", "surface under ocean 71 %". A percent from 1 percent,
# parts per million below it, so a trace reads as 400 ppm and not as 0.0004.
FRACTION = Ladder(
    quantity="fraction",
    base="frac",
    nonfinite="%",
    zero="%",
    rungs=(
        rung("ppm", start=0.0),
        rung("pct"),
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
        rung("nPa", start=0.0, systems=(AUTOMATIC, METRIC)),
        rung("uPa", systems=(AUTOMATIC, METRIC)),
        rung("mPa", systems=(AUTOMATIC, METRIC)),
        rung("Pa", systems=(AUTOMATIC, METRIC)),
        rung("kPa", systems=(AUTOMATIC, METRIC)),
        rung("MPa", systems=(AUTOMATIC, METRIC)),
        rung("GPa", systems=(AUTOMATIC, METRIC)),
        rung("TPa", systems=(AUTOMATIC, METRIC)),
        # Customary: psi all the way down (zero reads "0 psi"), with atm beside it.
        rung("psi", start=0.0, systems=(CUSTOMARY,)),
    ),
    companions=(
        companion("atm", systems=(AUTOMATIC, CUSTOMARY)),
        companion("psi", systems=(AUTOMATIC,)),
        companion("bar", min_size=10 * physical_constants.STANDARD_ATMOSPHERE_PA, max_size=1e4 * physical_constants.STANDARD_ATMOSPHERE_PA, systems=(AUTOMATIC, METRIC)),
        companion("mbar", min_size=1.0, max_size=0.01 * physical_constants.STANDARD_ATMOSPHERE_PA, systems=(AUTOMATIC, METRIC)),
    ),
)

# Gravity: m/s² with standard gravity beside it ("3.71 m/s² (0.378 g)").
GRAVITY = Ladder(
    quantity="gravity",
    base="m/s²",
    nonfinite="m/s²",
    rungs=(
        rung("m_s2", start=0.0, systems=(AUTOMATIC, METRIC)),
        rung("ft_s2", start=0.0, systems=(CUSTOMARY,)),
    ),
    companions=(companion("g0", symbol="g", name="times standard gravity", systems=(AUTOMATIC, CUSTOMARY)),),
)

# Mass (UX.22; Boss, 2026-10-10: lunar-mass rung wanted): kilograms (scientific
# by the number rule) under 0.01 lunar mass, then lunar masses to 0.1 Earth
# mass, Earth masses to 0.1 Jupiter mass, Jupiter masses to the hydrogen-burning
# limit of 0.075 solar mass, and solar masses above.
MASS = Ladder(
    quantity="mass",
    base="kg",
    nonfinite="kg",
    zero="kg",
    rungs=(
        rung("kg", start=0.0, systems=(AUTOMATIC, METRIC)),
        rung("lb", start=0.0, systems=(CUSTOMARY,)),
        rung("Mmoon", symbol="lunar masses", singular="lunar mass", start=0.01 * physical_constants.MOON_MASS_KG, systems=(AUTOMATIC,)),
        rung("Mearth", symbol="Earth masses", singular="Earth mass", start=_mearth_to_kg(0.1), systems=(AUTOMATIC,)),
        rung("Mjup", symbol="Jupiter masses", singular="Jupiter mass", start=_mjup_to_kg(0.1), systems=(AUTOMATIC,)),
        rung("Msun", symbol="solar masses", singular="solar mass", start=_msun_to_kg(0.075), systems=(AUTOMATIC,)),
    ),
)

# Luminosity and power (UX.22): watts with the SI prefixes, then solar
# luminosities from a millionth of the Sun's (3.8e20 W), where a star's
# "1.00 L☉" reads better than "383 EW". Below one watt (Hawking luminosities)
# scientific text, by the tiny-number rule.
POWER = Ladder(
    quantity="power",
    base="W",
    nonfinite="W",
    zero="W",
    rungs=(
        rung("W", start=0.0),
        rung("kW"),
        rung("MW"),
        rung("GW"),
        rung("TW"),
        rung("PW"),
        rung("EW"),
        # Metric only: more SI prefixes in place of the solar luminosity.
        rung("ZW", systems=(METRIC,)),
        rung("YW", systems=(METRIC,)),
        rung("RW", systems=(METRIC,)),
        rung("QW", systems=(METRIC,)),
        rung("Lsun", start=1e-6 * physical_constants.SOLAR_LUMINOSITY, systems=(AUTOMATIC, CUSTOMARY)),
    ),
)

# Heat flow and flux (UX.22): W/m² with SI prefixes, so Earth's 0.087 W/m²
# reads "87 mW/m²" as in the literature.
FLUX = Ladder(
    quantity="flux",
    base="W/m²",
    nonfinite="W/m²",
    zero="W/m²",
    rungs=(
        rung("nW_m2", start=0.0),
        rung("uW_m2"),
        rung("mW_m2"),
        rung("W_m2"),
        rung("kW_m2"),
        rung("MW_m2"),
    ),
)

# Gas and dust in clouds are always quoted per cubic centimeter (a nebula runs
# from a fraction of a particle to millions), so the one rung stands in for
# the ladder, as the design note says for number densities.
NUMBER_DENSITY = Ladder(
    quantity="number density",
    base="particles/m³",
    nonfinite="particles/cm³",
    zero="particles/cm³",
    rungs=(rung("cm-3", start=0.0),),
)

# A magnetic field (UX.22): tesla with the SI prefixes for planets, then gauss
# from 10,000 T (a hundred million gauss), where the literature quotes the
# fields of pulsars and magnetars ("1.00 × 10¹² G").
MAGNETIC_FIELD = Ladder(
    quantity="magnetic field",
    base="T",
    nonfinite="G",
    zero="T",
    rungs=(
        rung("nT", start=0.0),
        rung("uT"),
        rung("mT"),
        rung("T"),
        rung("G", start=1e4, systems=(AUTOMATIC, CUSTOMARY)),
    ),
)

LADDERS = (LENGTH, SPEED, TIME, TEMPERATURE, TEMPERATURE_DIFFERENCE, FRACTION, PRESSURE, GRAVITY, MASS, POWER, FLUX, NUMBER_DENSITY, MAGNETIC_FIELD)
"""Tuple[Ladder, ...]: Every ladder, in the order `unitladders.js` lists them."""

_BY_KIND = dict(zip(
    ("length", "speed", "duration", "temperature", "delta_temperature", "fraction", "pressure", "acceleration", "mass", "power", "flux", "number_density",
     "magnetic_field"), LADDERS))
_BY_KIND_OF_LADDER = {ladder.quantity: kind for kind, ladder in _BY_KIND.items()}
LADDER_OF_KIND = {kind: ladder.quantity for kind, ladder in _BY_KIND.items()}
"""dict: The name of the ladder each registry kind is shown on (`unitladders.js` `LADDER_OF_KIND`)."""


def for_kind(kind):
    """The ladder a registry kind (`unit_registry.KINDS`) is shown on; `ValueError` for a kind with none (an angle, an
    instant, a velocity)."""
    try:
        return _BY_KIND[kind]
    except KeyError:
        raise ValueError(f"no ladder shows a {kind!r}") from None
