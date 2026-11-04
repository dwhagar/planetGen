# planetgen/physics/relations.py

"""
Named physical laws, declared once (GEN.586; docs/design/architecture/unified-unit-quantity.md, section 11)
==============================================================================================

Each law is one function here: its name, the kind of each input, the kind of its answer and a formula written with the
ordinary operators (`+ - * / **`, `sqrt`, `log`). Two entry points come from that one declaration:

* **The kernel** is the declared function itself: `relations.surface_gravity(mass_kg, radius_m)`. It takes plain floats
  in SI base units (kilograms, meters, seconds, kelvin, radians) and answers a plain float in SI, so a generation loop
  pays for nothing but the formula. A kernel never sees a `Quantity`, and a value in AU or solar masses is converted by
  the caller before the call.
* **The Quantity form** is for the edges (describe, the API, one-off values, tests):
  `relations.quantities.surface_gravity(mass_quantity, radius_quantity)` takes `Quantity` values of the declared kinds
  in any of their units and answers a `Quantity` of the declared kind in its base unit, every other unit one read away.
  It runs the very same formula on the inputs' SI values, so the two forms agree bit for bit.

The unit rules the formulas follow:

* A formula is written in SI. A constant comes from its `c` argument (`c.G`, `c.c`, `c.sigma`, `c.k_B`, `c.hbar`,
  `c.radian`, `c.turn`, a full turn in radians), never as a number of its own; `c` is the SI table for a kernel and a table of
  dimension probes for the dimension test, both read from `unit_registry.PHYSICAL_CONSTANTS`. A number that is a plain
  coefficient (2, 4/3, the Roche coefficient) is written in the formula; a dimensionless input (an albedo) has the kind
  `None`.
* A law in AU, years and solar masses is not a second formula: `unit_registry.constant_in("G", "au-yr-msun")` is 4 pi
  squared worked out from G's dimension, and a caller that wants an orbit in those units converts its inputs and its
  answer with `unit_registry.in_system`.
* **The dimension test** (`check_dimensions`, run by `mathcheck`'s invariant group) runs each formula once on probe
  objects that carry a dimension and checks the answer has the declared output kind's dimension. A length added to a
  mass, a missing square, a missing constant or a wrong root fails it. A factor is not a dimension, so a km where an m
  was meant is not found here: that is what the Quantity form is for, and what the kernels' SI rule prevents. Nothing
  is checked at run time.

The vector laws (the obliquity of an axis against an orbit normal, a galactic orbit from a position and a velocity) are
not here: the formulas take plain scalars, so GEN.584 declares them beside the columns that read them.
"""

import math
from dataclasses import dataclass
from fractions import Fraction
from types import SimpleNamespace
from typing import Callable, Tuple

from planetgen.physics import unit_registry as units
from planetgen.physics.quantity import Quantity

ROCHE_FLUID_COEFFICIENT = 2.44
"""float: The fluid-body Roche limit, d = 2.44 R (rho_M / rho_m)^(1/3): a self-gravitating satellite
with no strength comes apart inside it. A rigid body survives closer (1.26), so the fluid limit is the cautious one."""


# --- Dimension probes ----------------------------------------------------------------------------------------------

class DimensionError(ValueError):
    """A formula added or subtracted quantities of different dimensions, or took a function of one that has some."""


class Probe:
    """A stand-in for a number that carries a dimension (a `unit_registry.dimension` vector). The ordinary operators
    work on it as on a float: a product adds the exponents, a power multiplies them, a sum or a difference refuses a
    different dimension. A plain number in a formula counts as dimensionless."""

    __slots__ = ("vector",)

    def __init__(self, vector=(0, 0, 0, 0, 0)):
        self.vector = tuple(Fraction(exponent) for exponent in vector)

    @staticmethod
    def _of(value):
        if isinstance(value, Probe):
            return value.vector
        if isinstance(value, (int, float)):
            return (Fraction(0),) * len(units.DIMENSION_BASES)
        raise TypeError(f"a formula may combine numbers and probes, not {type(value).__name__}")

    def __mul__(self, other):
        return Probe(tuple(a + b for a, b in zip(self.vector, self._of(other))))

    __rmul__ = __mul__

    def __truediv__(self, other):
        return Probe(tuple(a - b for a, b in zip(self.vector, self._of(other))))

    def __rtruediv__(self, other):
        return Probe(tuple(b - a for a, b in zip(self.vector, self._of(other))))

    def __pow__(self, exponent):
        if isinstance(exponent, Probe):
            raise DimensionError("an exponent has no dimension")
        return Probe(tuple(a * Fraction(exponent).limit_denominator(1000) for a in self.vector))

    def __add__(self, other):
        if self._of(other) != self.vector:
            raise DimensionError(f"cannot add {_describe(self.vector)} to {_describe(self._of(other))}")
        return Probe(self.vector)

    __radd__ = __sub__ = __rsub__ = __add__

    def __neg__(self):
        return Probe(self.vector)

    __abs__ = __pos__ = __neg__

    def __repr__(self):
        return f"Probe({_describe(self.vector)})"


def _describe(vector):
    parts = [f"{base}^{exponent}" for base, exponent in zip(units.DIMENSION_BASES, vector) if exponent]
    return " ".join(parts) or "dimensionless"


def sqrt(value):
    """The square root of a float, or of a probe (the dimension halved)."""
    return value ** 0.5 if isinstance(value, Probe) else math.sqrt(value)


def log(value):
    """The natural logarithm of a float, or of a probe, which must be dimensionless."""
    if isinstance(value, Probe):
        if any(value.vector):
            raise DimensionError(f"the logarithm of {_describe(value.vector)}")
        return Probe()
    return math.log(value)


# --- Constants -----------------------------------------------------------------------------------------------------

_TURN = 2.0 * math.pi
_ANGLE = units.dimension(angle=1)

SI = SimpleNamespace(**{name: value for name, (value, _vector) in units.PHYSICAL_CONSTANTS.items()}, radian=1.0,
                     turn=_TURN)
"""SimpleNamespace: The constants a kernel reads, in SI: G, c, sigma, k_B, hbar, `radian` (1, the angle unit a rate in
radians per second carries) and `turn`, a full turn in radians."""

PROBES = SimpleNamespace(**{name: Probe(vector) for name, (_value, vector) in units.PHYSICAL_CONSTANTS.items()},
                         radian=Probe(_ANGLE), turn=Probe(_ANGLE))
"""SimpleNamespace: The same names as `SI`, each a `Probe` of the constant's dimension (`turn` of angle)."""


# --- Declarations --------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Law:
    """One declared law."""

    name: str
    """str: The function's name, also its key in `LAWS`."""
    inputs: Tuple[Tuple[str, object], ...]
    """tuple: `(argument name, kind)` in call order; a kind of `unit_registry.KINDS`, or None for a plain number."""
    output: str
    """str: The kind of the answer, a key of `unit_registry.KINDS`."""
    formula: Callable
    """callable: The declared function, `formula(*inputs, c=SI)`: floats in, a float out, all in SI."""

    def input_dimension(self, kind):
        return (0,) * len(units.DIMENSION_BASES) if kind is None else units.KINDS[kind].dimension

    def check_dimensions(self):
        """Runs the formula once on probes of the declared input kinds. Returns the dimension found when it is the
        declared output kind's; raises `DimensionError` when it is not, or when the formula combines unlike dimensions."""
        probes = [Probe(self.input_dimension(kind)) for _argument, kind in self.inputs]
        answer = self.formula(*probes, c=PROBES)
        expected = tuple(Fraction(exponent) for exponent in units.KINDS[self.output].dimension)
        found = answer.vector if isinstance(answer, Probe) else Probe._of(answer)
        if found != expected:
            raise DimensionError(f"{self.name} gives {_describe(found)}, a {self.output} is {_describe(expected)}")
        return found

    def of_quantities(self, *given, **named):
        """The law on `Quantity` values (a plain number for a dimensionless input): the answer is a `Quantity` of the
        output kind in its base unit. The inputs' kinds are checked and each is read in its base unit."""
        arguments = dict(zip((argument for argument, _kind in self.inputs), given))
        for key, value in named.items():
            if key in arguments or key not in dict(self.inputs):
                raise TypeError(f"{self.name}() got an unexpected or repeated argument {key!r}")
            arguments[key] = value
        missing = [argument for argument, _kind in self.inputs if argument not in arguments]
        if missing or len(given) > len(self.inputs):
            raise TypeError(f"{self.name}() takes {', '.join(argument for argument, _kind in self.inputs)}")
        bases = []
        for argument, kind in self.inputs:
            value = arguments[argument]
            if kind is None:
                bases.append(float(value))
                continue
            if not isinstance(value, Quantity) or value.kind != kind:
                raise ValueError(f"{self.name}: {argument} must be a Quantity of kind {kind}, not {value!r}")
            bases.append(value.base())
        return Quantity(self.output, self.formula(*bases), units.KINDS[self.output].base)


LAWS = {}
"""dict: `Law.name` -> `Law`, every declared law."""

quantities = SimpleNamespace()
"""SimpleNamespace: The Quantity form of each law, by name: `quantities.surface_gravity(mass, radius)`."""


def law(output, **inputs):
    """Declares a law: the decorated function `name(inputs..., *, c=SI)` is the kernel and is returned unchanged; its
    Quantity form is set on `quantities`. `inputs` name each argument in order with its kind (None for a plain number)."""

    def declare(formula):
        name = formula.__name__
        if name in LAWS:
            raise ValueError(f"law {name!r} is declared twice")
        argument_names = formula.__code__.co_varnames[:formula.__code__.co_argcount]
        if tuple(inputs) != argument_names:
            raise ValueError(f"law {name!r}: kinds are named for {tuple(inputs)}, the function takes {argument_names}")
        for kind in (output, *inputs.values()):
            if kind is not None:
                units.kind(kind)
        declared = LAWS[name] = Law(name, tuple(inputs.items()), output, formula)
        setattr(quantities, name, declared.of_quantities)
        return formula

    return declare


def check_dimensions():
    """The names of the laws whose formula does not give its declared kind's dimension; an empty list when every law
    is dimensionally right."""
    wrong = []
    for declared in LAWS.values():
        try:
            declared.check_dimensions()
        except DimensionError:
            wrong.append(declared.name)
    return wrong


# --- The laws ------------------------------------------------------------------------------------------------------

@law("length", mass="mass")
def schwarzschild_radius(mass, *, c=SI):
    """The Schwarzschild radius, 2GM/c^2."""
    return 2 * c.G * mass / c.c ** 2


@law("acceleration", mass="mass", radius="length")
def surface_gravity(mass, radius, *, c=SI):
    """The gravity at a body's surface, GM/r^2."""
    return c.G * mass / radius ** 2


@law("speed", mass="mass", radius="length")
def escape_speed(mass, radius, *, c=SI):
    """The speed that escapes from a distance `radius` of `mass`, sqrt(2GM/r)."""
    return sqrt(2 * c.G * mass / radius)


@law("speed", mass="mass", radius="length")
def circular_speed(mass, radius, *, c=SI):
    """The speed of a circular orbit at `radius` round `mass`, sqrt(GM/r); at a body's own radius, its break-up speed."""
    return sqrt(c.G * mass / radius)


@law("angular_rate", mass="mass", radius="length")
def breakup_rate(mass, radius, *, c=SI):
    """The spin rate at which a body's equator orbits, sqrt(GM/r^3); faster and it sheds mass."""
    return c.radian * sqrt(c.G * mass / radius ** 3)


@law("speed", distance="length", semi_major_axis="length", mass="mass")
def vis_viva_speed(distance, semi_major_axis, mass, *, c=SI):
    """The speed of an orbit of semi-major axis `a` at `distance` from `mass`, sqrt(GM (2/r - 1/a))."""
    return sqrt(c.G * mass * (2 / distance - 1 / semi_major_axis))


@law("duration", semi_major_axis="length", mass="mass")
def orbital_period(semi_major_axis, mass, *, c=SI):
    """Kepler's third law, 2 pi sqrt(a^3 / GM), the mass being the pair's total."""
    return 2 * math.pi * sqrt(semi_major_axis ** 3 / (c.G * mass))


@law("length", period="duration", mass="mass")
def semi_major_axis(period, mass, *, c=SI):
    """Kepler's third law the other way, (GM T^2 / 4 pi^2)^(1/3)."""
    return (c.G * mass * period ** 2 / (4 * math.pi ** 2)) ** (1 / 3)


@law("angular_rate", semi_major_axis="length", mass="mass")
def mean_motion(semi_major_axis, mass, *, c=SI):
    """The mean motion of an orbit, sqrt(GM / a^3), in radians per second."""
    return c.radian * sqrt(c.G * mass / semi_major_axis ** 3)


@law("angular_rate", rotation_period="duration")
def rotation_rate(rotation_period, *, c=SI):
    """The angular rate of a spin of this period, a turn over the period."""
    return c.turn / rotation_period


@law("duration", angular_rate="angular_rate")
def rotation_period(angular_rate, *, c=SI):
    """The period of a spin at this angular rate, a turn over the rate."""
    return c.turn / angular_rate


@law("length", distance="length", planet_mass="mass", star_mass="mass")
def hill_radius(distance, planet_mass, star_mass, *, c=SI):
    """The Hill sphere of a planet at `distance` from a star, d (m / 3 (M + m))^(1/3)."""
    return distance * (planet_mass / (3 * (star_mass + planet_mass))) ** (1 / 3)


@law("length", primary_radius="length", primary_density="density", satellite_density="density")
def roche_limit(primary_radius, primary_density, satellite_density, *, c=SI):
    """The fluid Roche limit, 2.44 R (rho_M / rho_m)^(1/3)."""
    return ROCHE_FLUID_COEFFICIENT * primary_radius * (primary_density / satellite_density) ** (1 / 3)


@law("volume", radius="length")
def sphere_volume(radius, *, c=SI):
    """The volume of a sphere, 4/3 pi r^3."""
    return 4 / 3 * math.pi * radius ** 3


@law("density", mass="mass", radius="length")
def mean_density(mass, radius, *, c=SI):
    """The mean density of a sphere of this mass and radius."""
    return mass / (4 / 3 * math.pi * radius ** 3)


@law("length", mass="mass", density="density")
def sphere_radius(mass, density, *, c=SI):
    """The radius of a sphere of this mass and mean density."""
    return (3 * mass / (4 * math.pi * density)) ** (1 / 3)


@law("power", radius="length", temperature="temperature")
def luminosity(radius, temperature, *, c=SI):
    """The power a black body radiates, 4 pi R^2 sigma T^4."""
    return 4 * math.pi * radius ** 2 * c.sigma * temperature ** 4


@law("temperature", luminosity="power", radius="length")
def effective_temperature(luminosity, radius, *, c=SI):
    """The surface temperature of a black body of this luminosity and radius, (L / 4 pi R^2 sigma)^(1/4)."""
    return (luminosity / (4 * math.pi * radius ** 2 * c.sigma)) ** (1 / 4)


@law("flux", luminosity="power", distance="length")
def inverse_square_flux(luminosity, distance, *, c=SI):
    """The flux of a source at `distance`, L / 4 pi d^2."""
    return luminosity / (4 * math.pi * distance ** 2)


@law("temperature", luminosity="power", distance="length", albedo=None)
def equilibrium_temperature(luminosity, distance, albedo, *, c=SI):
    """The temperature of a fast-rotating body with no atmosphere at `distance` from a star, `((1 - A) L / 16 pi sigma
    d^2)^(1/4)`."""
    return ((1 - albedo) * luminosity / (16 * math.pi * c.sigma * distance ** 2)) ** (1 / 4)


@law("mass_flow", mass="mass", duration="duration")
def mass_flow(mass, duration, *, c=SI):
    """A mass lost or gained over a duration, at the mean rate."""
    return mass / duration
