# planetgen/physics/surface_dose.py

"""
Surface dose kernels (GEN.444, docs/design/composition/radiation-model.md sections 5 and 6)
==========================================================================================

Pure pieces of the ground-level dose that `physics/radiation.py` reads: the Band form of the stellar-particle spectrum
(`solar_particles` builds its measured events from it) and the nuclear-interaction shielding of a mixed column. The share of
an event that gets past a column and a magnetic cutoff is `solar_particles.dose_share` with `proton_range`
(NIST PSTAR, mixtures by Bragg's rule); the first slice's own copy of it, its rigidity function and its ICRP neutron weight
were retired (GEN.444): the share duplicated the production code, and the cosmic-ray table is already a dose weighted for
its neutrons, so a second weight would count them twice.

Shielding (section 5). A column of several gases is as thick for a primary as `sum X_k (lambda_air / lambda_k)` g/cm2 of air.
The nuclear interaction length of a gas follows from its elements, `1 / lambda = sum w_i / lambda_i` over the mass
fractions `w_i` of its elements (`ELEMENT_INTERACTION_LENGTH_G_CM2`, the Particle Data Group [R4926]); the rule gives the
Particle Data Group's own molecules to 0.2 percent (CO2 88.95 against 88.9, H2O 83.35 against 83.3, CH4 73.7 against 73.8).
The reference, air, is the same rule over dry Earth air, so an N2/O2 sky counts as air to the gram. An element not in the
table is an error, not a guess. The protons of a stellar event are stopped by ionization, and that stopping goes by Bragg's
rule over the same elements (`stopping_mixture`, `proton_range.energy_mev` takes the mixture).
"""

import functools
import math

from planetgen.physics import elements
from planetgen.physics.elements import DataError

# --- Stellar particles ------------------------------------------------------------------------------------------


class BandSpectrum:
    """
    The integral rigidity spectrum of a solar particle event, `J(>R) = J0 R^-g1 exp(-R/R0)` up to `R1 = (g2 - g1) R0`
    and `J0 R1^-g1 exp(-R1/R0) (R/R1)^-g2` above, continuous at `R1` (the published form, not the
    discontinuous one of the research note). Rigidity in GV; `J0` sets the unit of the flux and cancels in a share.
    """

    def __init__(self, j0, gamma1, gamma2, r0_gv):
        if j0 <= 0 or r0_gv <= 0 or gamma2 <= gamma1:
            raise DataError("a Band spectrum needs J0 > 0, R0 > 0 and gamma2 > gamma1")
        self.j0, self.gamma1, self.gamma2, self.r0_gv = j0, gamma1, gamma2, r0_gv
        self.r1_gv = (gamma2 - gamma1) * r0_gv

    def flux(self, rigidity_gv):
        """`J(>R)`."""
        r = max(rigidity_gv, 1e-9)
        if r <= self.r1_gv:
            return self.j0 * r ** -self.gamma1 * math.exp(-r / self.r0_gv)
        return (self.j0 * self.r1_gv ** -self.gamma1 * math.exp(-self.r1_gv / self.r0_gv)
                * (r / self.r1_gv) ** -self.gamma2)


# --- Shielding by a mixed column --------------------------------------------------------------------------------

ELEMENT_INTERACTION_LENGTH_G_CM2 = {
    "H": 52.0, "D": 71.8, "He": 71.0, "Li": 71.3, "C": 85.8, "N": 89.7, "O": 90.2, "F": 97.4, "Ne": 99.0, "Na": 102.6, "Mg": 104.1,
    "Al": 107.2, "Si": 108.4, "P": 111.4, "S": 112.4, "Cl": 115.7, "Ar": 119.7, "K": 119.0, "Ca": 119.8, "Ti": 126.2, "V": 128.5,
    "Cr": 129.3, "Fe": 132.1, "Xe": 172.1, "Rn": 204.3,
}
"""dict: The nuclear interaction length, g/cm2, of each element the gas table holds, as the Particle Data Group's atomic and nuclear properties
tables print it [R4926]: the gas where the element is one at room temperature (H2, N2, O2, F2, Cl2, the noble gases), the solid otherwise.
Krypton and bromine are not in those tables and come from `_INTERPOLATED` below."""

_INTERPOLATED = ("Kr", "Br")
"""tuple: The elements whose length is the power law in atomic weight through iron and xenon (stand-in SD1)."""


def _interpolate_missing():
    fe, xe = elements.element("Fe").atomic_weight, elements.element("Xe").atomic_weight
    slope = math.log(ELEMENT_INTERACTION_LENGTH_G_CM2["Xe"] / ELEMENT_INTERACTION_LENGTH_G_CM2["Fe"]) / math.log(xe / fe)
    for symbol in _INTERPOLATED:
        ELEMENT_INTERACTION_LENGTH_G_CM2[symbol] = ELEMENT_INTERACTION_LENGTH_G_CM2["Fe"] * (elements.element(symbol).atomic_weight / fe) ** slope


_interpolate_missing()


DEUTERIUM_WEIGHT = 2.014
"""float: The atomic weight of deuterium, g/mol, as published (the element table has hydrogen's natural weight only)."""


def _atomic_weight(symbol):
    if symbol == "D":
        return DEUTERIUM_WEIGHT
    element = elements.element(symbol)
    return element.atomic_weight if element.atomic_weight is not None else float(element.longest_lived_mass_number)


def element_mass_fractions(columns):
    """The mass fraction of each element of a mapping of a gas formula (or `"air"`) to an amount (a column, g/cm2, or a mass): the amount of a
    gas is split by the mass of its atoms. Deuterium is its own element (HD, HDO); a leading mass number is dropped (`3He`)."""
    totals = {}
    for gas, amount in columns.items():
        counts = _atoms(gas)
        mass = sum(_atomic_weight(symbol) * count for symbol, count in counts.items())
        for symbol, count in counts.items():
            totals[symbol] = totals.get(symbol, 0.0) + amount * _atomic_weight(symbol) * count / mass
    return totals


@functools.lru_cache(maxsize=None)
def _atoms(gas):
    """The atoms of a gas formula as `{symbol: count}`."""
    if gas == "air":
        return _AIR_ATOMS
    text = gas.lstrip("0123456789")
    deuterium = text.count("D")
    counts = dict(elements.parse_formula(text.replace("D", ""))) if text.replace("D", "") else {}
    if deuterium:
        counts["D"] = deuterium
    return counts


DRY_AIR_MOLE_FRACTIONS = {"N2": 0.78084, "O2": 0.20946, "Ar": 0.00934}
"""dict: Dry Earth air at sea level, by mole fraction, as far as it matters for a shield (the trace gases weigh under 0.05 percent)."""


def interaction_length_g_cm2(mass_fractions):
    """The nuclear interaction length, g/cm2, of a material of element mass fractions (any total): `1 / lambda = sum w_i / lambda_i`."""
    total = sum(mass_fractions.values())
    if total <= 0:
        raise DataError("a material needs a positive mass")
    inverse = 0.0
    for symbol, weight in mass_fractions.items():
        try:
            inverse += weight / total / ELEMENT_INTERACTION_LENGTH_G_CM2[symbol]
        except KeyError:
            raise DataError("no nuclear interaction length for the element %r; the elements are %s"
                            % (symbol, ", ".join(ELEMENT_INTERACTION_LENGTH_G_CM2))) from None
    return 1.0 / inverse


def _air_atoms():
    atoms = {}
    for gas, fraction in DRY_AIR_MOLE_FRACTIONS.items():
        for symbol, count in elements.parse_formula(gas).items():
            atoms[symbol] = atoms.get(symbol, 0.0) + fraction * count
    return atoms


_AIR_ATOMS = _air_atoms()
AIR_INTERACTION_LENGTH_G_CM2 = interaction_length_g_cm2(element_mass_fractions({"air": 1.0}))
"""float: The nuclear interaction length of dry Earth air by the same rule as every other gas, g/cm2 (the Particle Data Group prints 90.1 for it;
the reference of the shield is the rule's own number so that Earth's air counts as air exactly)."""


def effective_air_column_g_cm2(columns_g_cm2):
    """The column, g/cm2 of air, that stops as many primaries as `columns_g_cm2`, a mapping of a gas formula (or
    "air") to its own column in g/cm2: `sum X_k (lambda_air / lambda_k)`, `lambda_k` from the gas's elements."""
    total = 0.0
    for gas, column in columns_g_cm2.items():
        if column:
            total += column * AIR_INTERACTION_LENGTH_G_CM2 / interaction_length_g_cm2(element_mass_fractions({gas: 1.0}))
    return total


def _is_element(name):
    try:
        elements.element(name)
    except elements.UnknownElement:
        return False
    return name == elements.element(name).symbol


def stopping_mixture(columns, rows, atomic_numbers=None):
    """The `{row: mass fraction}` of a mixture of gas columns for Bragg's rule in `proton_range` (`rows`, the rows the range table holds, of which the elements count): an
    element the table does not hold takes the row of the nearest atomic number, the lighter on a tie, and deuterium takes hydrogen (stand-in 6.223 SD2)."""
    symbols = [row for row in rows if _is_element(row)]
    by_z = {elements.element(symbol).z: symbol for symbol in symbols}
    mix = {}
    for symbol, weight in element_mass_fractions(columns).items():
        if symbol in by_z.values():
            row = symbol
        elif symbol == "D":
            row = "H"
        else:
            z = elements.element(symbol).z
            row = by_z[min(by_z, key=lambda other: (abs(other - z), other))]
        mix[row] = mix.get(row, 0.0) + weight
    return mix
