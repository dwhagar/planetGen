# planetgen/physics/unit_registry.py

"""
One registry of every unit the project uses (GEN.568; docs/design/architecture/unified-unit-quantity.md)
============================================================================================

Every unit the code, the database, the display ladders and the browser use is one row here, with its kind (what it
measures), the factor that takes it to the kind's base unit, the offset of its zero, a display symbol, an ASCII key
and a name. Nothing else in the repository states a factor twice once the sweeps of GEN.571 have run; this module
touches no existing file, so the duplicate names in `constants.py` stay until then.

The rules:

* **SI base units.** A kind's base is the SI unit of its dimension (m, kg, s, K, rad, and the compounds of them:
  m/s, W, Pa, ...). A unit's `factor` is how many base units one of it is, and its `offset` the base value that is
  zero in it (255.372 K for a degree Fahrenheit), so `base = value * factor + offset` and `value = (base - offset) /
  factor`. Only the temperature scales have an offset.
* **Exact numbers come from `physics/constants.py`.** The registry reads them; the only numbers it states itself
  are definitions (the SI prefixes, 60 s in a minute, 360 degrees in a turn, the 0.3048 m foot is read from
  `constants`). The tests check each factor against its definition written out in full.
* **A dimension is five integer exponents** of (length, mass, time, temperature, angle). A unit declared with
  `parts` (km/s is km over s) has the dimension and factor its parts give, and the tests check both against its
  kind's dimension and its own factor.
* **A kind says what values are allowed** (`Bounds`): finite always; a mass, duration, power, ... at least 0; a speed
  below c (`check`). The checks are for the edges, where a value is stored, loaded or shown; the numeric kernels stay
  on plain floats (the Quantity class of GEN.569 builds on this).
* **A temperature difference is its own kind** (`delta_temperature`, no offsets): 10 K warmer is 10 degrees Celsius
  warmer, not 283 degrees Celsius.
* **Relations are data** (`RELATIONS`): a pair of kinds and the two functions between them, in base units
  (the mass and its gravitational parameter, a radius and its diameter, a period and its frequency, a rotation
  period and its angular rate).
* **Keys are ASCII and unique** across the registry (`km_s`, `Msun`, `degC`, `uW`); symbols are what the display
  shows (`km/s`, `M☉`, `°C`, `µW`) and are unique within a kind.
* **Unit systems and constants** (GEN.585): a `UnitSystem` is four registry units (length, mass, time, temperature); SI
  is the kernels' own and `au-yr-msun` is the orbit work's. A physical constant's number in a system (`constant_in`)
  is worked out from its dimension and the system's base units, never typed. The typed `4 pi squared` of
  `kepler.AU3_PER_YR2_PER_SOLAR_MASS` is 39.4784; G for this solar mass in the Julian year is 39.4769, 3.8e-5 less
  (the typed one is the Gaussian year's); the sweep of GEN.587 decides which stands.
* **No astropy.** The registry has no dependency beyond `constants`.

The functions are pure.
"""

import math
from dataclasses import dataclass
from typing import Callable, Optional, Tuple

from planetgen.physics import constants

DIMENSION_BASES = ("length", "mass", "time", "temperature", "angle")
"""tuple: What the five exponents of a dimension count, in order (m, kg, s, K, rad)."""


def dimension(length=0, mass=0, time=0, temperature=0, angle=0):
    """The dimension with these exponents of (length, mass, time, temperature, angle)."""
    return (length, mass, time, temperature, angle)


def multiply_dimensions(first, second, power=1):
    """`first * second ** power`, the exponents added."""
    return tuple(a + power * b for a, b in zip(first, second))


@dataclass(frozen=True)
class Bounds:
    """The values a kind allows, in its base unit. `low` and `high` are inclusive unless open."""

    low: Optional[float] = None
    high: Optional[float] = None
    low_open: bool = False
    high_open: bool = False

    def holds(self, value):
        """Whether `value` is inside."""
        if self.low is not None and (value <= self.low if self.low_open else value < self.low):
            return False
        if self.high is not None and (value >= self.high if self.high_open else value > self.high):
            return False
        return True

    def describe(self):
        """The bounds as text, for an error: "at least 0", "below 299792458", "any value"."""
        parts = []
        if self.low is not None:
            parts.append(f"{'above' if self.low_open else 'at least'} {self.low:g}")
        if self.high is not None:
            parts.append(f"{'below' if self.high_open else 'at most'} {self.high:g}")
        return " and ".join(parts) or "any value"


ANY = Bounds()
NON_NEGATIVE = Bounds(low=0.0)


@dataclass(frozen=True)
class Kind:
    """What a unit measures."""

    name: str
    """str: "length", "mass", "speed", ..."""
    dimension: Tuple[int, ...]
    """tuple: The five exponents (`DIMENSION_BASES`)."""
    base: str
    """str: The key of the unit whose factor is 1 and offset 0."""
    bounds: Bounds = ANY
    """Bounds: The values allowed, in the base unit."""


@dataclass(frozen=True)
class Unit:
    """One row of the registry."""

    key: str
    """str: ASCII, unique in the registry ("km_s", "Msun", "degC", "uW")."""
    kind: str
    """str: The `Kind.name` it measures."""
    symbol: str
    """str: Shown after a number ("km/s", "M☉", "°C"); unique within its kind."""
    name: str
    """str: Spelled out, plural ("kilometers per second")."""
    factor: float
    """float: Base units in one of it."""
    offset: float = 0.0
    """float: The base value that is zero in it."""
    singular: Optional[str] = None
    """str, optional: The name for exactly one ("kilometer per second"); `name` when unset."""
    parts: Tuple[Tuple[str, int], ...] = ()
    """tuple: `((unit key, exponent), ...)` this unit is made of (km/s: km 1, s -1); the tests check its factor and
    its kind's dimension against them."""


@dataclass(frozen=True)
class Relation:
    """A declared pair of kinds with the functions between them, both in base units (`apply`)."""

    name: str
    """str: "mass-mu", "radius-diameter", ..."""
    first_role: str
    first_kind: str
    second_role: str
    second_kind: str
    forward: Callable
    """callable: The first quantity's base value to the second's."""
    inverse: Callable
    """callable: The second's to the first's."""


KINDS = {}
"""dict: `Kind.name` -> `Kind`."""
UNITS = {}
"""dict: `Unit.key` -> `Unit`."""
RELATIONS = {}
"""dict: `Relation.name` -> `Relation`."""

_BY_KIND = {}

SI_PREFIXES = (
    ("n", "n", "nano", -9), ("u", "µ", "micro", -6), ("m", "m", "milli", -3), ("", "", "", 0),
    ("k", "k", "kilo", 3), ("M", "M", "mega", 6), ("G", "G", "giga", 9), ("T", "T", "tera", 12),
    ("P", "P", "peta", 15), ("E", "E", "exa", 18), ("Z", "Z", "zetta", 21), ("Y", "Y", "yotta", 24),
    ("R", "R", "ronna", 27), ("Q", "Q", "quetta", 30),
)
"""tuple: `(ASCII key prefix, symbol prefix, name prefix, power of ten)`: the SI prefixes the ladders use."""
PARSEC_PREFIXES = (("m", "m", "milli", -3), ("c", "c", "centi", -2), ("", "", "", 0), ("k", "k", "kilo", 3),
                   ("M", "M", "mega", 6), ("G", "G", "giga", 9))
"""tuple: The prefixes of the parsec (the ladders' mpc, cpc, pc, kpc, Mpc, Gpc)."""


def _power_of_ten(exponent):
    """`10 ** exponent` correctly rounded, so a prefix factor is the same float as its literal (1e-9, not 1.0000000000000001e-9)."""
    return float(f"1e{exponent}")


def _kind(name, dim, base, bounds=ANY):
    KINDS[name] = Kind(name, dim, base, bounds)
    _BY_KIND[name] = []


def _unit(key, kind, symbol, name, factor, offset=0.0, singular=None, parts=()):
    if key in UNITS:
        raise ValueError(f"unit key {key!r} is declared twice")
    UNITS[key] = Unit(key, kind, symbol, name, factor, offset, singular, tuple(parts))
    _BY_KIND[kind].append(key)


def _prefixed(kind, key, symbol, name, factor, prefixes=SI_PREFIXES, singular=None, parts=(), name_suffix=""):
    """One unit per prefix: `key` ("W"), `symbol` ("W"), `name` ("watts") with the prefix put before it; `parts` of the
    unprefixed unit are not carried (a prefixed unit has the base's factor times its prefix)."""
    for ascii_prefix, symbol_prefix, name_prefix, power in prefixes:
        _unit(ascii_prefix + key, kind, symbol_prefix + symbol, name_prefix + name + name_suffix,
              factor * _power_of_ten(power) if power else factor,
              singular=None if singular is None else name_prefix + singular + name_suffix,
              parts=parts if not power else ())


# --- Kinds ---------------------------------------------------------------------------------------------------------

_LENGTH, _MASS, _TIME = dimension(length=1), dimension(mass=1), dimension(time=1)
_SPEED = dimension(length=1, time=-1)
_ACCELERATION = dimension(length=1, time=-2)
_POWER = dimension(length=2, mass=1, time=-3)
_ENERGY = dimension(length=2, mass=1, time=-2)
_PRESSURE = dimension(length=-1, mass=1, time=-2)

_kind("length", _LENGTH, "m", NON_NEGATIVE)
_kind("mass", _MASS, "kg", NON_NEGATIVE)
_kind("duration", _TIME, "s", NON_NEGATIVE)
_kind("instant", _TIME, "unix_s")
_kind("frequency", dimension(time=-1), "Hz", NON_NEGATIVE)
_kind("speed", _SPEED, "m_s", Bounds(low=0.0, high=float(constants.SPEED_OF_LIGHT_M_S), high_open=True))
_kind("acceleration", _ACCELERATION, "m_s2")
_kind("temperature", dimension(temperature=1), "K", NON_NEGATIVE)
_kind("delta_temperature", dimension(temperature=1), "dK")
_kind("fraction", dimension(), "frac", NON_NEGATIVE)
_kind("angle", dimension(angle=1), "rad")
_kind("angular_rate", dimension(time=-1, angle=1), "rad_s")
_kind("angular_momentum", dimension(length=2, mass=1, time=-1), "kg_m2_s")
_kind("power", _POWER, "W", NON_NEGATIVE)
_kind("energy", _ENERGY, "J", NON_NEGATIVE)
_kind("force", dimension(length=1, mass=1, time=-2), "N")
_kind("mass_flow", dimension(mass=1, time=-1), "kg_s", NON_NEGATIVE)
_kind("pressure", _PRESSURE, "Pa", NON_NEGATIVE)
_kind("density", dimension(length=-3, mass=1), "kg_m3", NON_NEGATIVE)
_kind("column_mass", dimension(length=-2, mass=1), "kg_m2", NON_NEGATIVE)
_kind("area", dimension(length=2), "m2", NON_NEGATIVE)
_kind("surface_number_density", dimension(length=-2), "m-2", NON_NEGATIVE)
_kind("volume", dimension(length=3), "m3", NON_NEGATIVE)
_kind("number_density", dimension(length=-3), "m-3", NON_NEGATIVE)
_kind("flux", dimension(mass=1, time=-3), "W_m2", NON_NEGATIVE)
_kind("magnetic_field", dimension(mass=1, time=-2), "T", NON_NEGATIVE)
_kind("gravitational_parameter", dimension(length=3, time=-2), "m3_s2", NON_NEGATIVE)


# --- Units ---------------------------------------------------------------------------------------------------------

_MINUTE_S, _HOUR_S, _DAY_S = 60.0, 3600.0, 86400.0
_JULIAN_YEAR_S = constants.SECONDS_PER_YEAR


def _lengths():
    _prefixed("length", "m", "m", "meters", 1.0, singular="meter")
    for key, symbol, name, factor in (
            ("au", "AU", "astronomical units", constants.AU_M),
            ("ly", "ly", "light-years", constants.LIGHTYEAR_M),
            ("ld", "ld", "light-days", constants.SPEED_OF_LIGHT_M_S * _DAY_S),
            ("cm", "cm", "centimeters", 1e-2),
            ("in", "in", "inches", constants.INCH_M),
            ("ft", "ft", "feet", constants.FOOT_M),
            ("mi", "mi", "miles", constants.MILE_M),
            ("Rsun", "R☉", "solar radii", constants.SOLAR_RADIUS_M),
            ("Rearth", "R⊕", "Earth radii", constants.EARTH_RADIUS_KM * constants.KM_M),
            ("Rjup", "R♃", "Jupiter radii", constants.JUPITER_RADIUS_KM * constants.KM_M)):
        _unit(key, "length", symbol, name, factor)
    _prefixed("length", "pc", "pc", "parsecs", constants.PARSEC_M, PARSEC_PREFIXES, singular="parsec")


def _masses():
    _unit("kg", "mass", "kg", "kilograms", 1.0, singular="kilogram")
    _unit("g", "mass", "g", "grams", 1e-3, singular="gram")
    for key, symbol, name, factor in (
            ("lb", "lb", "pounds", constants.POUND_KG),
            ("Mmoon", "M☾", "lunar masses", constants.MOON_MASS_KG),
            ("Mearth", "M⊕", "Earth masses", constants.EARTH_MASS_TO_KG),
            ("Mjup", "M♃", "Jupiter masses", constants.JUPITER_MASS_TO_KG),
            ("Msun", "M☉", "solar masses", constants.SOLAR_MASS_TO_KG)):
        _unit(key, "mass", symbol, name, factor)


def _times():
    _unit("s", "duration", "s", "seconds", 1.0, singular="second")
    for ascii_prefix, symbol_prefix, name_prefix, power in SI_PREFIXES:
        if ascii_prefix in ("n", "u", "m"):
            _unit(ascii_prefix + "s", "duration", symbol_prefix + "s", name_prefix + "seconds", _power_of_ten(power),
                  singular=name_prefix + "second")
    for key, symbol, name, singular, factor in (
            ("min", "min", "minutes", "minute", _MINUTE_S), ("h", "h", "hours", "hour", _HOUR_S),
            ("d", "d", "days", "day", _DAY_S), ("yr", "yr", "years", "year", _JULIAN_YEAR_S),
            ("kyr", "kyr", "thousand years", None, _JULIAN_YEAR_S * 1e3),
            ("Myr", "Myr", "million years", None, _JULIAN_YEAR_S * 1e6),
            ("Gyr", "Gyr", "billion years", None, _JULIAN_YEAR_S * 1e9)):
        _unit(key, "duration", symbol, name, factor, singular=singular)
    _unit("unix_s", "instant", "s (Unix)", "seconds since the Unix epoch", 1.0)
    _unit("unix_ms", "instant", "ms (Unix)", "milliseconds since the Unix epoch", 1e-3)
    _unit("unix_yr", "instant", "yr (Unix)", "years since the Unix epoch", _JULIAN_YEAR_S)
    _prefixed("frequency", "Hz", "Hz", "hertz", 1.0, SI_PREFIXES[2:9])
    _unit("1/s", "frequency", "1/s", "per second", 1.0, parts=(("s", -1),))
    _unit("1/yr", "frequency", "1/yr", "per year", 1.0 / _JULIAN_YEAR_S, parts=(("yr", -1),))
    _unit("1/Gy", "frequency", "1/Gy", "per billion years", 1.0 / (_JULIAN_YEAR_S * 1e9), parts=(("Gyr", -1),))


def _kinematics():
    _unit("m_s", "speed", "m/s", "meters per second", 1.0, singular="meter per second", parts=(("m", 1), ("s", -1)))
    for key, symbol, name, factor, parts in (
            ("km_s", "km/s", "kilometers per second", 1e3, (("km", 1), ("s", -1))),
            ("Mm_s", "Mm/s", "megameters per second", 1e6, (("Mm", 1), ("s", -1))),
            ("au_d", "AU/d", "astronomical units per day", constants.AU_M / _DAY_S, (("au", 1), ("d", -1))),
            ("au_yr", "AU/yr", "astronomical units per year", constants.AU_M / _JULIAN_YEAR_S,
             (("au", 1), ("yr", -1))),
            ("pc_Myr", "pc/Myr", "parsecs per million years", constants.PARSEC_M / (_JULIAN_YEAR_S * 1e6),
             (("pc", 1), ("Myr", -1))),
            ("km_h", "km/h", "kilometers per hour", 1e3 / _HOUR_S, (("km", 1), ("h", -1))),
            ("mph", "mph", "miles per hour", constants.MILE_M / _HOUR_S, (("mi", 1), ("h", -1))),
            ("mi_s", "mi/s", "miles per second", constants.MILE_M, (("mi", 1), ("s", -1))),
            ("c", "c", "times the speed of light", constants.SPEED_OF_LIGHT_M_S, ())):
        _unit(key, "speed", symbol, name, factor, parts=parts)
    _unit("m_s2", "acceleration", "m/s²", "meters per second squared", 1.0, parts=(("m", 1), ("s", -2)))
    _unit("ft_s2", "acceleration", "ft/s²", "feet per second squared", constants.FOOT_M,
          parts=(("ft", 1), ("s", -2)))
    _unit("g0", "acceleration", "g₀", "standard gravities", constants.EARTH_GRAVITY)


def _thermal():
    _prefixed("temperature", "K", "K", "kelvins", 1.0, SI_PREFIXES[:4], singular="kelvin")
    _unit("degC", "temperature", "°C", "degrees Celsius", 1.0, offset=constants.CELSIUS_ZERO_K, singular="degree Celsius")
    _unit("degF", "temperature", "°F", "degrees Fahrenheit", 5.0 / 9.0, offset=constants.CELSIUS_ZERO_K - 32 * 5.0 / 9.0,
          singular="degree Fahrenheit")
    _unit("dK", "delta_temperature", "K", "kelvins", 1.0, singular="kelvin")
    _unit("ddegC", "delta_temperature", "°C", "degrees Celsius", 1.0, singular="degree Celsius")
    _unit("ddegF", "delta_temperature", "°F", "degrees Fahrenheit", 5.0 / 9.0, singular="degree Fahrenheit")


def _angular():
    _unit("rad", "angle", "rad", "radians", 1.0, singular="radian")
    _unit("deg", "angle", "°", "degrees", math.pi / 180.0, singular="degree")
    _unit("arcmin", "angle", "′", "arcminutes", math.pi / 180.0 / 60.0, singular="arcminute")
    _unit("arcsec", "angle", "″", "arcseconds", math.pi / 180.0 / 3600.0, singular="arcsecond")
    _unit("mas", "angle", "mas", "milliarcseconds", math.pi / 180.0 / 3.6e6, singular="milliarcsecond")
    _unit("turn", "angle", "turn", "turns", 2.0 * math.pi, singular="turn")
    _unit("rad_s", "angular_rate", "rad/s", "radians per second", 1.0, parts=(("rad", 1), ("s", -1)))
    _unit("deg_d", "angular_rate", "°/d", "degrees per day", math.pi / 180.0 / _DAY_S, parts=(("deg", 1), ("d", -1)))
    _unit("turn_d", "angular_rate", "turn/d", "turns per day", 2.0 * math.pi / _DAY_S, parts=(("turn", 1), ("d", -1)))
    _unit("kg_m2_s", "angular_momentum", "kg m²/s", "kilogram meters squared per second", 1.0,
          parts=(("kg", 1), ("m", 2), ("s", -1)))


def _energetic():
    _prefixed("power", "W", "W", "watts", 1.0, singular="watt")
    _unit("Lsun", "power", "L☉", "solar luminosities", constants.SOLAR_LUMINOSITY)
    _unit("erg_s", "power", "erg/s", "ergs per second", 1e-7, parts=(("erg", 1), ("s", -1)))
    _prefixed("energy", "J", "J", "joules", 1.0, SI_PREFIXES[2:9], singular="joule")
    _unit("erg", "energy", "erg", "ergs", 1e-7)
    _unit("eV", "energy", "eV", "electronvolts", constants.ELECTRON_VOLT_J, singular="electronvolt")
    _unit("keV", "energy", "keV", "kiloelectronvolts", constants.ELECTRON_VOLT_J * 1e3, singular="kiloelectronvolt")
    _unit("MeV", "energy", "MeV", "megaelectronvolts", constants.ELECTRON_VOLT_J * 1e6, singular="megaelectronvolt")
    _prefixed("force", "N", "N", "newtons", 1.0, SI_PREFIXES[2:9], singular="newton")
    _unit("dyn", "force", "dyn", "dynes", 1e-5)
    _unit("kg_s", "mass_flow", "kg/s", "kilograms per second", 1.0, parts=(("kg", 1), ("s", -1)))
    _unit("g_s", "mass_flow", "g/s", "grams per second", 1e-3, parts=(("g", 1), ("s", -1)))
    _unit("Msun_yr", "mass_flow", "M☉/yr", "solar masses per year", constants.SOLAR_MASS_TO_KG / _JULIAN_YEAR_S,
          parts=(("Msun", 1), ("yr", -1)))
    _unit("Mearth_yr", "mass_flow", "M⊕/yr", "Earth masses per year", constants.EARTH_MASS_TO_KG / _JULIAN_YEAR_S,
          parts=(("Mearth", 1), ("yr", -1)))
    _prefixed("pressure", "Pa", "Pa", "pascals", 1.0, SI_PREFIXES[:12], singular="pascal")
    _unit("bar", "pressure", "bar", "bars", 1e5)
    _unit("mbar", "pressure", "mbar", "millibars", 1e2)
    _unit("atm", "pressure", "atm", "standard atmospheres", constants.STANDARD_ATMOSPHERE_PA)
    _unit("psi", "pressure", "psi", "pounds per square inch", constants.PSI_PA)
    _prefixed("flux", "W_m2", "W/m²", "watts per square meter", 1.0, SI_PREFIXES[:12],
              singular="watt per square meter")
    _prefixed("magnetic_field", "T", "T", "teslas", 1.0, SI_PREFIXES[:4], singular="tesla")
    _unit("G", "magnetic_field", "G", "gauss", 1e-4)


def _fields():
    _unit("m2", "area", "m²", "square meters", 1.0, parts=(("m", 2),))
    _unit("km2", "area", "km²", "square kilometers", 1e6, parts=(("km", 2),))
    _unit("cm2", "area", "cm²", "square centimeters", 1e-4)
    _unit("m3", "volume", "m³", "cubic meters", 1.0, parts=(("m", 3),))
    _unit("km3", "volume", "km³", "cubic kilometers", 1e9, parts=(("km", 3),))
    _unit("cm3", "volume", "cm³", "cubic centimeters", 1e-6)
    _unit("L", "volume", "L", "liters", 1e-3, singular="liter")
    _unit("kg_m3", "density", "kg/m³", "kilograms per cubic meter", 1.0, parts=(("kg", 1), ("m", -3)))
    _unit("g_cm3", "density", "g/cm³", "grams per cubic centimeter", 1e3)
    _unit("kg_m2", "column_mass", "kg/m²", "kilograms per square meter", 1.0, parts=(("kg", 1), ("m", -2)))
    _unit("g_cm2", "column_mass", "g/cm²", "grams per square centimeter", 10.0)
    _unit("kg_km2", "column_mass", "kg/km²", "kilograms per square kilometer", 1e-6, parts=(("kg", 1), ("km", -2)))
    _unit("kg_au2", "column_mass", "kg/AU²", "kilograms per square astronomical unit", constants.AU_M ** -2, parts=(("kg", 1), ("au", -2)))
    _unit("Mearth_au2", "column_mass", "M⊕/AU²", "Earth masses per square astronomical unit", constants.EARTH_MASS_TO_KG / constants.AU_M ** 2,
          parts=(("Mearth", 1), ("au", -2)))
    _unit("m-2", "surface_number_density", "objects/m²", "objects per square meter", 1.0)
    _unit("km-2", "surface_number_density", "objects/km²", "objects per square kilometer", 1e-6, parts=(("km", -2),))
    _unit("au-2", "surface_number_density", "objects/AU²", "objects per square astronomical unit", constants.AU_M ** -2, parts=(("au", -2),))
    _unit("Msun_pc3", "density", "M☉/pc³", "solar masses per cubic parsec",
          constants.SOLAR_MASS_TO_KG / constants.PARSEC_M ** 3, parts=(("Msun", 1), ("pc", -3)))
    _unit("m-3", "number_density", "particles/m³", "particles per cubic meter", 1.0)
    _unit("cm-3", "number_density", "particles/cm³", "particles per cubic centimeter", constants.CM3_TO_M3)
    _unit("pc-3", "number_density", "particles/pc³", "particles per cubic parsec", constants.PARSEC_M ** -3,
          parts=(("pc", -3),))
    _unit("m3_s2", "gravitational_parameter", "m³/s²", "cubic meters per second squared", 1.0,
          parts=(("m", 3), ("s", -2)))
    _unit("km3_s2", "gravitational_parameter", "km³/s²", "cubic kilometers per second squared", 1e9,
          parts=(("km", 3), ("s", -2)))


def _fractions():
    """A share of a whole (GEN.770): the base is the plain fraction 0 to 1; the page writes it as a percent from 1 percent
    and as parts per million below."""
    _unit("frac", "fraction", "(fraction)", "fractions of the whole", 1.0)
    _unit("pct", "fraction", "%", "percent", 1e-2)
    _unit("ppm", "fraction", "ppm", "parts per million", 1e-6)


_lengths()
_masses()
_fractions()
_times()
_kinematics()
_thermal()
_angular()
_energetic()
_fields()


# --- Relations -----------------------------------------------------------------------------------------------------

def _relation(name, first_role, first_kind, second_role, second_kind, forward, inverse):
    RELATIONS[name] = Relation(name, first_role, first_kind, second_role, second_kind, forward, inverse)


_relation("mass-mu", "mass", "mass", "mu", "gravitational_parameter",
          lambda mass: constants.G * mass, lambda mu: mu / constants.G)
_relation("radius-diameter", "radius", "length", "diameter", "length", lambda radius: 2.0 * radius,
          lambda diameter: 0.5 * diameter)
_relation("period-frequency", "period", "duration", "frequency", "frequency", lambda period: 1.0 / period,
          lambda frequency: 1.0 / frequency)
_relation("rotation_period-angular_rate", "rotation_period", "duration", "angular_rate", "angular_rate",
          lambda period: 2.0 * math.pi / period, lambda rate: 2.0 * math.pi / rate)


# --- Unit systems and constants ------------------------------------------------------------------------------------

@dataclass(frozen=True)
class UnitSystem:
    """A coherent set of base units: one unit of length, mass, time and temperature (keys of the registry). The kernels
    work in SI; another system (AU, year and solar mass for orbits) is generated from the registry where a kernel wants
    one, never typed."""

    name: str
    length: str
    mass: str
    time: str
    temperature: str

    def base_factors(self):
        """`(length, mass, time, temperature)`, the SI value of one of each base unit."""
        return (UNITS[self.length].factor, UNITS[self.mass].factor, UNITS[self.time].factor,
                UNITS[self.temperature].factor)


UNIT_SYSTEMS = {}
"""dict: `UnitSystem.name` -> `UnitSystem`."""


def _system(name, length, mass, time, temperature="K"):
    for key, expected in ((length, "length"), (mass, "mass"), (time, "duration"), (temperature, "temperature")):
        if UNITS[key].kind != expected or UNITS[key].offset:
            raise ValueError(f"unit system {name!r}: {key!r} is not a {expected} unit without an offset")
    UNIT_SYSTEMS[name] = UnitSystem(name, length, mass, time, temperature)


_system("si", "m", "kg", "s")
_system("au-yr-msun", "au", "Msun", "yr")

PHYSICAL_CONSTANTS = {
    "G": (constants.G, dimension(length=3, mass=-1, time=-2)),
    "c": (constants.SPEED_OF_LIGHT_M_S, dimension(length=1, time=-1)),
    "sigma": (constants.STEFAN_BOLTZMANN_CONSTANT, dimension(mass=1, time=-3, temperature=-4)),
    "k_B": (constants.BOLTZMANN, dimension(length=2, mass=1, time=-2, temperature=-1)),
    "hbar": (constants.REDUCED_PLANCK, dimension(length=2, mass=1, time=-1)),
}
"""dict: Name -> `(SI value, dimension)` of the constants the kernels use; a constant's number in another system is
worked out from its dimension (`constant_in`)."""


def system_factor(dimension_vector, system):
    """How many SI units one unit of dimension `dimension_vector` is in `system` (`UnitSystem` or its name): the product
    of each base unit's SI value to the dimension's exponent; the angle exponent counts for nothing (an angle is
    dimensionless here)."""
    system = UNIT_SYSTEMS[system] if isinstance(system, str) else system
    factor = 1.0
    for base, exponent in zip(system.base_factors(), (dimension_vector[0], dimension_vector[1], dimension_vector[2],
                                                    dimension_vector[3])):
        if exponent:
            factor *= base ** exponent
    return factor


def in_system(value_si, dimension_vector, system):
    """`value_si`, a value of dimension `dimension_vector` in SI, as its number in `system`."""
    return value_si / system_factor(dimension_vector, system)


def constant_in(name, system="si"):
    """The number of the physical constant `name` ("G", "c", "sigma", "k_B", "hbar") in `system`: G is 4 pi squared
    (about) in AU, year and solar mass, worked out here from its dimension and the three units, not typed."""
    try:
        value_si, vector = PHYSICAL_CONSTANTS[name]
    except KeyError:
        raise ValueError(f"unknown constant {name!r} ({', '.join(PHYSICAL_CONSTANTS)})") from None
    return in_system(value_si, vector, system)


# --- Lookups ---------------------------------------------------------------------------------------------------------

def kind(name):
    """The `Kind` called `name`."""
    try:
        return KINDS[name]
    except KeyError:
        raise ValueError(f"unknown kind {name!r} (the kinds are {', '.join(KINDS)})") from None


def unit(key):
    """The `Unit` with this key."""
    try:
        return UNITS[key]
    except KeyError:
        raise ValueError(f"unknown unit {key!r}") from None


def units_of(kind_name):
    """The units of a kind, in the order they were declared."""
    kind(kind_name)
    return tuple(UNITS[key] for key in _BY_KIND[kind_name])


def unit_dimension(key):
    """The dimension a unit's kind gives it."""
    return kind(unit(key).kind).dimension


def parts_dimension(key):
    """The dimension a compound unit's `parts` give it (`None` for a unit with no parts)."""
    parts = unit(key).parts
    if not parts:
        return None
    result = dimension()
    for part_key, exponent in parts:
        result = multiply_dimensions(result, unit_dimension(part_key), exponent)
    return result


def parts_factor(key):
    """The factor a compound unit's `parts` give it (`None` for a unit with no parts)."""
    parts = unit(key).parts
    if not parts:
        return None
    result = 1.0
    for part_key, exponent in parts:
        result *= unit(part_key).factor ** exponent
    return result


# --- Conversion and checks -------------------------------------------------------------------------------------------

def to_base(value, key):
    """`value` of the unit `key` in its kind's base unit."""
    row = unit(key)
    return value * row.factor + row.offset if row.offset else value * row.factor


def from_base(base_value, key):
    """A base value of the kind of unit `key`, in that unit."""
    row = unit(key)
    return (base_value - row.offset) / row.factor if row.offset else base_value / row.factor


_TEN_POWERS = tuple(float(f"1e{exponent}") for exponent in range(23))
"""tuple: `10.0 ** n` for n = 0 to 22, the powers of ten a float holds exactly."""


def scale(value, source_factor, target_factor):
    """
    `value * source_factor / target_factor`, exact where the ratio of the factors is a power of ten: then `value` is
    multiplied or divided by that power, one correctly rounded step (km to m, 1.5 -> 1500.0; m to km, 1500 -> 1.5),
    where the two-step product would round twice (`v * 1e3 / 1e-6` differs from `v * 1e9` for a third of all `v`).
    """
    ratio = source_factor / target_factor
    if ratio > 0.0 and math.isfinite(ratio):
        exponent = round(math.log10(ratio))
        if -22 <= exponent <= 22 and abs(ratio - float(f"1e{exponent}")) <= 1e-14 * ratio:
            return value * _TEN_POWERS[exponent] if exponent >= 0 else value / _TEN_POWERS[-exponent]
    return value * source_factor / target_factor


def convert(value, from_key, to_key):
    """`value` of the unit `from_key` in the unit `to_key`; a unit of another kind raises `ValueError`."""
    source, target = unit(from_key), unit(to_key)
    if source.kind != target.kind:
        raise ValueError(f"cannot convert {source.symbol} ({source.kind}) to {target.symbol} ({target.kind})")
    if source.factor == target.factor and source.offset == target.offset:
        return value
    if source.offset == target.offset:
        return scale(value, source.factor, target.factor)
    return (value * source.factor + source.offset - target.offset) / target.factor


def converter(from_key, to_key):
    """
    A function that turns a number of the unit `from_key` into the unit `to_key`, for a module-level name a kernel calls
    (`km_to_m = converter("km", "m")`, then `km_to_m(radius_km)`). The two units are looked up once, so a unit of
    another kind raises `ValueError` when the module loads, not when the kernel runs. The factor between them is
    worked out once and applied in one step, as the hand-written sites did: a multiply by the ratio of the larger
    unit to the smaller (`AU_TO_KM`, 1000.0) or a divide by the ratio of the smaller to the larger (`/ AU_TO_KM`,
    `/ 1000.0`), a power of ten taken as the exact power. Two units of equal factor return the value unchanged.
    """
    source, target = unit(from_key), unit(to_key)
    if source.kind != target.kind:
        raise ValueError(f"cannot convert {source.symbol} ({source.kind}) to {target.symbol} ({target.kind})")
    if source.offset != target.offset:
        return lambda value: (value * source.factor + source.offset - target.offset) / target.factor
    if source.factor == target.factor:
        return lambda value: value
    if source.factor > target.factor:
        factor = _step(source.factor / target.factor)
        return lambda value: value * factor
    factor = _step(target.factor / source.factor)
    return lambda value: value / factor


def _step(ratio):
    """`ratio`, or the power of ten it is within rounding of."""
    exponent = round(math.log10(ratio))
    if -22 <= exponent <= 22 and abs(ratio - float(f"1e{exponent}")) <= 1e-14 * ratio:
        return _TEN_POWERS[exponent]
    return ratio


def check(value, key):
    """`value` of the unit `key`, returned when it is finite and inside the kind's bounds; `ValueError` otherwise."""
    row = unit(key)
    if not math.isfinite(value):
        raise ValueError(f"a {row.kind} must be finite, not {value!r}")
    base_value = to_base(value, key)
    bounds = KINDS[row.kind].bounds
    if not bounds.holds(base_value):
        raise ValueError(f"a {row.kind} must be {bounds.describe()} (in {KINDS[row.kind].base}), "
                         f"not {value:g} {row.symbol}")
    return value


def apply_relation(name, base_value, from_role):
    """The partner of `base_value` (in base units) in the relation `name`, the one `from_role` names being given."""
    try:
        relation = RELATIONS[name]
    except KeyError:
        raise ValueError(f"unknown relation {name!r} ({', '.join(RELATIONS)})") from None
    if from_role == relation.first_role:
        return relation.forward(base_value)
    if from_role == relation.second_role:
        return relation.inverse(base_value)
    raise ValueError(f"{name} relates {relation.first_role} and {relation.second_role}, not {from_role!r}")
