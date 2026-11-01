# planetgen/physics/mixture_base.py

"""
What every family of the atmosphere mixture draw shares (GEN.479; docs/design/composition/mixture-draw.md)
==============================================================================================

`MixtureDraw`, the error, the seeded design constant, the mapping of the solver's species to the gas table, and the
step that closes a draw: take an equilibrium, drop what the gas table does not hold (refusing when that is too much), and
normalise. The families (`mixture_rocky`, `mixture_giant`, `mixture_frost`, `mixture_lava`) each decide the element
budget, the temperature, the pressure and the oxygen fugacity; `mixture_draw.draw_mixture` picks the family.
"""

import math
from dataclasses import dataclass, field

from planetgen.physics import gases, species_thermo
from planetgen.util import draw

UNLISTED_LIMIT = 1e-3
"""float: The most mole fraction that species the gas table does not hold may add up to before the draw says so (`MixtureError`)."""

MIN_FRACTION = 1e-30
"""float: Species below this mole fraction are dropped from the result."""


class MixtureError(RuntimeError):
    """A mixture that cannot be drawn: a quantity the tables read that the subject lacks, or a gas outside the gas table."""


class NotApplicable(MixtureError):
    """A regime the recipe drew whose physics does not hold at this body's temperature or pressure (an ice atmosphere on a surface above the
    ices' critical points): the dispatcher leaves the body to its outgassed air."""


@dataclass(frozen=True)
class MixtureDraw:
    """
    A drawn mixture.

    Attributes:
        mole_fractions (dict): Gas-table formula -> mole fraction, summing to 1 (species below `MIN_FRACTION` dropped).
        family (str): Which chemistry made it (`rocky`, `giant`, `brown-dwarf`, `hot-jupiter`, `envelope`, `frost`, `lava`,
            `coma`).
        drawn (dict): Recipe group -> the row drawn (`{"a": "a6", "b": "b2", "c": "c3"}`).
        redox_diw (float or None): log10 fO2 above the iron-wuestite buffer of the mantle the gas is in contact with
            ([R1951] scale); `None` for a gas that has no mantle (a giant).
        temperature_k (float): The temperature of the equilibrium the mixture is the result of.
        pressure_bar (float): The pressure of the equilibrium.
        inventory (dict): The element budget the chemistry read (atoms per 100 g of body for a rocky body, per hydrogen atom for a
            giant), by element symbol.
        constants (dict): The design constants this draw drew, name -> value.
        condensed (dict): What left the gas above the layer (clouds and ices): name -> a condensation temperature or the moles taken.
        recipe (RecipeDraw or None): The recipe draw the mixture followed (set by `mixture_draw.draw_mixture`).
        unlisted (dict): Species the equilibrium holds that the gas table does not (their mole fractions before the remainder was
            renormalised): atoms and radicals of a hot gas.
        moles (float or None): The moles of gas the element budget makes, per 100 g of mantle (the scale of `inventory`: the elements it
            holds as atoms per 100 g); `None` for a family whose pressure is not a supply (a giant's reference level, the vapour over a
            melt or an ice). After the loss of the light gases (`mixture_escape`) it is what stays.
    """

    mole_fractions: dict
    family: str
    drawn: dict
    redox_diw: float | None
    temperature_k: float
    pressure_bar: float
    inventory: dict
    constants: dict = field(default_factory=dict)
    condensed: dict = field(default_factory=dict)
    unlisted: dict = field(default_factory=dict)
    recipe: object = None
    moles: float | None = None


def bracket(constants, name, spec):
    """Draws one design constant of `spec` = (low, high, shape) ("log" or "uniform") and records it under `name`."""
    low, high, shape = spec
    if low == high:
        value = low
    elif shape == "log":
        value = 10.0 ** draw.uniform(math.log10(low), math.log10(high))
    else:
        value = draw.uniform(low, high)
    constants[name] = value
    return value


def require(subject, names):
    """Raises `MixtureError` naming every quantity of `names` the subject lacks."""
    missing = [name for name in names if subject.get(name) is None]
    if missing:
        raise MixtureError(f"the subject of class {subject.body_class!r} lacks {', '.join(missing)}")


def atoms_of(equilibrium):
    """The moles of atoms of every element in the gas of `equilibrium`, oxygen included."""
    atoms = {}
    for name, moles in equilibrium.moles.items():
        for symbol, count in species_thermo.SPECIES[name].elements.items():
            atoms[symbol] = atoms.get(symbol, 0.0) + count * moles
    return atoms


def to_gas_table(fractions):
    """`(listed, unlisted)`: the mole fractions of species names as the gas table's formulas, and the ones it does not hold."""
    listed, unlisted = {}, {}
    for name, value in fractions.items():
        formula = species_thermo.SPECIES[name].formula
        target = listed if gases.has(formula) else unlisted
        target[formula] = target.get(formula, 0.0) + value
    return listed, unlisted


def finish(equilibrium, *, family, drawn, redox_diw, inventory, constants, condensed=None, supplied=False):
    """
    Closes a draw: the equilibrium's species as gas-table formulas, normalised, largest first. `supplied` says the equilibrium's
    abundances are an element supply per 100 g of mantle, so its total moles are the draw's `moles`.

    Raises:
        MixtureError: More than `UNLISTED_LIMIT` of the gas is species the gas table lacks.
    """
    listed, unlisted = to_gas_table({name: value for name, value in equilibrium.mole_fractions.items() if value > MIN_FRACTION})
    if sum(unlisted.values()) > UNLISTED_LIMIT:
        raise MixtureError(f"{sum(unlisted.values()):.3g} of the gas is species the gas table lacks: {sorted(unlisted)}")
    total = sum(listed.values())
    return MixtureDraw(
        mole_fractions={formula: value / total for formula, value in sorted(listed.items(), key=lambda item: -item[1])},
        family=family,
        drawn=drawn,
        redox_diw=redox_diw,
        temperature_k=equilibrium.temperature_k,
        pressure_bar=equilibrium.pressure_bar,
        inventory=inventory,
        constants=constants,
        condensed=condensed or {},
        unlisted=unlisted,
        moles=equilibrium.total_moles if supplied else None,
    )

