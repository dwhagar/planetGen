# planetgen/physics/air_column.py

"""
The air above one place, as the light meets it (GEN.417, GEN.443; docs/design/composition/radiation-model.md section 3).

An `AirColumn` is the mixture, the gravity and the total pressure: from them follow the column mass `X = P / g`
(g/cm2), the mean molar mass and the number of molecules of each gas above a square centimetre. Every optical depth
of the atmosphere (Rayleigh scattering, ultraviolet absorption, the greenhouse) is a sum over these column numbers, so
this is the one place that turns partial pressures into columns and the one place that reads the stored partial
pressures of a body.

The mixture is a plain mapping of gas-table formula (`physics/gases.py`, `"CO2"`, `"Ar"`) to partial pressure in Pa. The
functions in `physics/rayleigh.py`, `physics/aerosols.py` and `physics/atmospheric_radiation.py` take an `AirColumn`
and never read a planet's columns, so a change of how the mixture is stored (GEN.411, GEN.424) changes `column_of`
and nothing else.

Units: Pa, m/s2, kg/mol, g/cm2 and cm-2 (the unit the optics tables are published in).
"""

from dataclasses import dataclass

from planetgen.physics import atmosphere, constants, gases

AVOGADRO = 6.02214076e23
"""float: Molecules per mole (the 2019 SI value)."""


class EmptyColumn(ValueError):
    """The mixture has no pressure: there is no air above the ground to ask about."""


@dataclass(frozen=True)
class AirColumn:
    """
    The air above a place at the surface.

    Attributes:
        partials_pa (dict): Gas-table formula to partial pressure at the surface, Pa; zero entries are dropped.
        gravity_m_s2 (float): Surface gravity, m/s2.
        pressure_pa (float): The sum of the partial pressures, Pa.
        mean_molar_mass_kg_mol (float): The mole-fraction-weighted molar mass, kg/mol.
        column_g_cm2 (float): The air mass above one square centimetre, `X = P / g`, g/cm2.
    """

    partials_pa: dict
    gravity_m_s2: float
    pressure_pa: float
    mean_molar_mass_kg_mol: float
    column_g_cm2: float

    def fraction(self, formula):
        """The mole fraction of `formula` (0 when it is not in the mixture)."""
        return self.partials_pa.get(formula, 0.0) / self.pressure_pa

    def number_column_cm2(self, formula):
        """The molecules of `formula` above one square centimetre, cm-2: `p_i N_A / (m g)` with `m` the mean molar mass."""
        pressure = self.partials_pa.get(formula, 0.0)
        return pressure * AVOGADRO / (self.mean_molar_mass_kg_mol * self.gravity_m_s2) * 1e-4

    def total_number_column_cm2(self):
        """All the molecules above one square centimetre, cm-2."""
        return self.pressure_pa * AVOGADRO / (self.mean_molar_mass_kg_mol * self.gravity_m_s2) * 1e-4


def from_partials(partials_pa, gravity_m_s2):
    """
    An `AirColumn` from partial pressures (Pa, keyed by gas-table formula) and surface gravity (m/s2).

    Raises:
        EmptyColumn: The partial pressures sum to zero or gravity is not above zero.
        gases.UnknownGas: A formula is not in the gas table.
    """
    kept = {formula: float(value) for formula, value in partials_pa.items() if value and value > 0.0}
    total = sum(kept.values())
    if total <= 0.0 or not gravity_m_s2 or gravity_m_s2 <= 0.0:
        raise EmptyColumn("the air column needs a pressure above zero and a gravity above zero")
    mean_molar = sum(value * gases.gas(formula).molar_mass for formula, value in kept.items()) / total
    return AirColumn(kept, float(gravity_m_s2), total, mean_molar, total / gravity_m_s2 * 0.1)


def partials_pa(body):
    """
    A body's stored mixture as formula -> Pa, the balance gases included (`atmosphere.gas_partials_pa`), less any gas with no
    Rayleigh cross section in the gas table: sodium (class B's second balance gas) is an exosphere, not a column.
    """
    return {formula: value for formula, value in atmosphere.gas_partials_pa(body).items() if gases.gas(formula).rayleigh is not None}


def column_of(body):
    """
    The `AirColumn` of a planet or moon, or `None` for a body with no air or no gravity (an airless body has no column).
    """
    gravity = (getattr(body, "gravity", None) or 0.0) * constants.EARTH_GRAVITY
    mixture = partials_pa(body)
    if gravity <= 0.0 or sum(mixture.values()) <= 0.0:
        return None
    return from_partials(mixture, gravity)
