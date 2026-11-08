# planetgen/physics/disc_condensation.py

"""
What a protoplanetary disc condenses at a temperature, species by species (GEN.1363 part A; docs/design/composition/disc-condensation.md)
==========================================================================================================================================

`condensate_shares` has three ices (water 182 K, ammonia 131 K, methane 41 K) and puts all the carbon in methane. A disc measured by its
gas and ice is not like that [R3321, R3322]: about half of the carbon is refractory (organic and graphite-like grains that survive to about
500 K); the rest is mostly CO and CO2, with a little CH4; most of the nitrogen is N2; sulfur is mostly in troilite with some H2S; and the noble gases
come out as clathrates in the water ice. This module is the budget and the lines of that disc, and nothing draws from it yet (part B wires the
ice shares of a body to it and part C the C/O of a giant).

1. **The budget.** Per 10^4 hydrogen atoms (the unit of [R3322: Table 1]) the Sun has C 2.88, N 0.68, O 4.90, S 0.13, Ar 0.024, Kr 1.3e-5 and Xe 1.7e-6
   ([R0168: Table 2], photospheric log eps; `tuning.STELLAR_SOLAR_LOG_EPS` and `NOBLE_LOG_EPS`). A star's ratios (its element over its rock-forming
   metals, over the Sun's: 1 is the Sun, as `condensate_shares` takes them) scale each. Rock holds `condensate_shares.ROCK_OXYGEN_SHARE` of the Sun's oxygen.
2. **The carbon split.** `REFRACTORY_CARBON_SHARE` of the carbon is grains. The rest goes `CH4_SHARE_OF_VOLATILE_CARBON` to methane and the remainder to CO and CO2 in the
   5 to 1 of [R3322: Table 1] (CO 1.5 and CO2 0.3 per 10^4 hydrogen). The oxygen that CO and CO2 hold is out of the water (a gas takes oxygen whether or not it is
   frozen), so the water is the oxygen the rock, the CO and the CO2 leave.
3. **The nitrogen and the sulfur.** `NITROGEN_AS_N2` of the nitrogen is N2, `NITROGEN_AS_SALTS` ammonium salts and the rest ammonia hydrate; `SULFUR_AS_H2S`
   of the sulfur is H2S ice (the rest is troilite, in the rock, from 704 K).
4. **The lines.** Each species freezes below its line (`LINES_K`). A line is the temperature of the table the repo already uses for water, ammonia and
   methane (the equilibrium at 10^-4 bar, [R2726: Table 10]); CO and CO2 are the disc-density values of [R3322: Table 1], N2 and H2S are recalled
   (registered so). Noble-gas clathrate lines are [R2726: Table 8].
5. **The C/O across a line.** The gas C/O rises outward past each snow line (water takes oxygen first, then CO2, then CO); the solids' C/O is the
   other half. A giant that gathers gas between the water and CO lines is carbon-rich, one that gathers solids is oxygen-rich [R3322].

Limits. Equilibrium supply: every species freezes whole below its line. A body that formed warm and trapped a species in water ice (a comet) keeps a tiny
share of it (`comet_ice_formation`, `comet_ice_budget.TRAPPED_FRACTION`), and that is for the body, not the disc. The lines do not move as the disc cools. The H2S and the
ammonium-salt shares of the nitrogen and sulfur are stand-ins (registry rows DC4 to DC6).
"""

import dataclasses
import math

from planetgen import tuning
from planetgen.physics import condensate_shares, constants as physical_constants, unit_registry as units
from planetgen.util.checks import finite_domain

NOBLE_LOG_EPS = {"Kr": 3.12, "Xe": 2.22}
"""dict: The Sun's photospheric log eps (H = 12) of krypton and xenon ([R0168: Table 2]; argon is in `tuning.STELLAR_SOLAR_LOG_EPS`). [Stand-in DC1, tier T1, tweak L0.]"""

LINES_K = {
    "H2O": condensate_shares.ICE_CONDENSATION_K["H2O"],
    "NH3": condensate_shares.ICE_CONDENSATION_K["NH3"],
    "H2S": 80.0,
    "Xe": 69.0,
    "Kr": 53.0,
    "Ar": 48.0,
    "CO2": 47.0,
    "CH4": condensate_shares.ICE_CONDENSATION_K["CH4"],
    "N2": 22.0,
    "CO": 20.0,
}
"""dict: The temperature (K) below which each ice freezes, warmest first. Water, ammonia and methane are the repo's equilibrium values at 10^-4 bar; the noble-gas
clathrates are [R2726: Table 8] (Xe 69, Kr 53, Ar 48); CO 20 and CO2 47 are [R3322: Table 1] (CO2 42 to 52 K there, 70 K in other discs); N2 22 K is recalled
([R3321] section 1.2 puts CH4, CO and N2 near 20 to 40 K) and H2S 80 K is recalled (ice sublimates near 50 to 90 K by the matrix). [Stand-in DC2, tier T1 for the read
lines and T3 for N2 and H2S, tweak L0.]"""

REFRACTORY_LINE_K = 500.0
"""float: The temperature below which the refractory carbon is a grain ([R3322: Table 1] adopt more than 150 K and 500 K in the model). [Stand-in DC3, tier T2, tweak L0.]"""

REFRACTORY_CARBON_SHARE = 0.5
"""float: The share of the carbon that is refractory grains: about half ([R3321] section 1.2, from the interstellar medium; [R3322: Table 1] has 0.6 of 2.4 of the carbon
they account for in the disc). [Stand-in DC3, tier T2, tweak L1.]"""

CO_PER_CO2 = 1.5 / 0.3
"""float: The CO over the CO2 of the volatile carbon that is not methane ([R3322: Table 1], CO 1.5 and CO2 0.3 molecules per 10^4 hydrogen). [Stand-in DC3, tier T2, tweak L0.]"""

CH4_SHARE_OF_VOLATILE_CARBON = 0.34 / (0.34 + 3.1 + 4.7)
"""float: The methane share of the volatile carbon: comet 67P's CH4 over CH4 plus CO plus CO2, 0.34, 3.1 and 4.7 percent of the water ([R3835: Table 2]), 4.2 percent.
One comet, taken for the disc's methane. [Stand-in DC4, tier T3, tweak L1.]"""

NITROGEN_AS_N2 = 0.6
NITROGEN_AS_SALTS = 0.1
"""float: The shares of the nitrogen that are N2 and ammonium salts, the rest being ammonia hydrate: most nitrogen is N2 and unobserved ([R3321] section 1.2) and the
salts hold a minority ([R3321], citing the Rosetta result second-hand; the split is recalled). [Stand-ins DC5, tier T3, tweak L3.]"""

SULFUR_AS_H2S = 0.3
"""float: The share of the sulfur that is H2S ice: 67P's H2S is 1.1 percent of the water ([R3835: Table 2]) against 3.5 percent for all the sulfur over the water (S 0.13 over the
3.8 oxygen the rock leaves). The rest is troilite in the rock. [Stand-in DC6, tier T3, tweak L2.]"""

MOLAR_MASS = {"H2O": 18.015, "NH3": 17.031, "NH4": 17.031, "H2S": 34.08, "Xe": 131.29, "Kr": 83.80, "Ar": 39.948, "CO2": 44.009, "CH4": 16.043, "N2": 28.014,
              "CO": 28.010, "C": 12.011}
"""dict: The molar mass (g/mol) of each species ('NH4' is the ammonium salts, taken as one ammonia of nitrogen: the anion's mass is not counted)."""

_m_to_au = units.converter("m", "au")

ORDER = tuple(LINES_K)
"""tuple: The ices warmest first."""

PERCENT_PER_H2O_MOLE = condensate_shares.TABLE_11_MASS_PERCENT["H2O"] / (
    10.0 ** (tuning.STELLAR_SOLAR_LOG_EPS["O"] - 12.0) * 1e4 * (1.0 - condensate_shares.ROCK_OXYGEN_SHARE))
"""float: The percent of the nebula one mole of water per 10^4 hydrogen makes, calibrated on [R2726: Table 11] (water 0.571 percent for the 3.8 oxygen the rock leaves); the percent of
any species is its moles times its molar mass over water's times this."""


def solar_per_1e4_h(element):
    """The Sun's atoms of `element` per 10^4 hydrogen [R0168]; `element` is a symbol in `tuning.STELLAR_SOLAR_LOG_EPS` or `NOBLE_LOG_EPS`."""
    log_eps = tuning.STELLAR_SOLAR_LOG_EPS.get(element, NOBLE_LOG_EPS.get(element))
    if log_eps is None:
        raise ValueError(f"no solar abundance for {element!r}")
    return 1e4 * 10.0 ** (log_eps - 12.0)


@dataclasses.dataclass(frozen=True)
class Inventory:
    """
    The molecules per 10^4 hydrogen that a disc has to condense, by species, and the oxygen and carbon the rock and the grains hold.

    Attributes:
        moles (dict): Species (`ORDER`) to its molecules, or atoms for a noble gas, per 10^4 hydrogen.
        rock_oxygen (float): The oxygen atoms in the rock (a fixed share of the Sun's).
        refractory_carbon (float): The carbon atoms in grains.
        oxygen, carbon (float): The star's total atoms of each per 10^4 hydrogen.
    """

    moles: dict
    rock_oxygen: float
    refractory_carbon: float
    oxygen: float
    carbon: float


@finite_domain()
def inventory(carbon=1.0, oxygen=1.0, nitrogen=1.0, sulfur=1.0, noble=1.0):
    """
    The `Inventory` of a disc whose star has these element ratios (each over the rock-forming metals, over the Sun's; 1 is the Sun).

    The oxygen of the CO and CO2 comes out of the water; with too little oxygen the water is zero and the CO and CO2 keep their carbon (the oxygen is then short, and `oxygen_short`
    of the result is what is missing).
    """
    for name, value in (("carbon", carbon), ("oxygen", oxygen), ("nitrogen", nitrogen), ("sulfur", sulfur), ("noble", noble)):
        if not value >= 0.0:
            raise ValueError(f"the {name} ratio must be zero or more, not {value!r}")
    total_c = solar_per_1e4_h("C") * carbon
    total_o = solar_per_1e4_h("O") * oxygen
    total_n = solar_per_1e4_h("N") * nitrogen
    total_s = solar_per_1e4_h("S") * sulfur
    rock_o = condensate_shares.ROCK_OXYGEN_SHARE * solar_per_1e4_h("O")
    refractory = REFRACTORY_CARBON_SHARE * total_c
    volatile = total_c - refractory
    ch4 = CH4_SHARE_OF_VOLATILE_CARBON * volatile
    co2 = (volatile - ch4) / (CO_PER_CO2 + 1.0)
    co = volatile - ch4 - co2
    water = max(total_o - rock_o - co - 2.0 * co2, 0.0)
    moles = {
        "H2O": water,
        "NH3": total_n * (1.0 - NITROGEN_AS_N2 - NITROGEN_AS_SALTS),
        "H2S": total_s * SULFUR_AS_H2S,
        "Xe": solar_per_1e4_h("Xe") * noble,
        "Kr": solar_per_1e4_h("Kr") * noble,
        "Ar": solar_per_1e4_h("Ar") * noble,
        "CO2": co2,
        "CH4": ch4,
        "N2": 0.5 * total_n * NITROGEN_AS_N2,
        "CO": co,
    }
    return Inventory(moles=moles, rock_oxygen=rock_o, refractory_carbon=refractory, oxygen=total_o, carbon=total_c)


def ammonium_salts(carbon=1.0, oxygen=1.0, nitrogen=1.0, sulfur=1.0, noble=1.0):
    """The ammonium salts, nitrogen atoms per 10^4 hydrogen; they freeze with the ammonia hydrate (`LINES_K['NH3']`)."""
    return solar_per_1e4_h("N") * nitrogen * NITROGEN_AS_SALTS


@finite_domain()
def condensed(temperature_k, carbon=1.0, oxygen=1.0, nitrogen=1.0, sulfur=1.0, noble=1.0):
    """
    The species frozen at `temperature_k`, as `{species: molecules per 10^4 hydrogen}` ('NH4' for the ammonium salts, 'C' for the refractory carbon grains). A species is in the result when
    the temperature is below its line (`LINES_K`, `REFRACTORY_LINE_K`) and it has any.
    """
    if not temperature_k > 0.0:
        raise ValueError("the temperature must be positive")
    held = inventory(carbon, oxygen, nitrogen, sulfur, noble)
    found = {}
    for species, line in LINES_K.items():
        if temperature_k < line and held.moles[species] > 0.0:
            found[species] = held.moles[species]
    if temperature_k < LINES_K["NH3"]:
        salts = ammonium_salts(carbon, oxygen, nitrogen, sulfur, noble)
        if salts > 0.0:
            found["NH4"] = salts
    if temperature_k < REFRACTORY_LINE_K and held.refractory_carbon > 0.0:
        found["C"] = held.refractory_carbon
    return found


def masses_percent(temperature_k, carbon=1.0, oxygen=1.0, nitrogen=1.0, sulfur=1.0, noble=1.0):
    """
    The mass each condensate makes at `temperature_k`, in percent of the nebula, as `condensate_shares.condensed_masses` gives it: `rock` (the silicates, oxides, metal and
    troilite of [R2726: Table 11], whole below the troilite line), the grains `C`, and each ice in `condensed`.
    """
    masses = {condensate_shares.ROCK: condensate_shares.TABLE_11_MASS_PERCENT[condensate_shares.ROCK]}
    for species, moles in condensed(temperature_k, carbon, oxygen, nitrogen, sulfur, noble).items():
        masses[species] = moles * MOLAR_MASS[species] / MOLAR_MASS["H2O"] * PERCENT_PER_H2O_MOLE
    return masses


@finite_domain()
def gas_and_solid_c_over_o(temperature_k, carbon=1.0, oxygen=1.0, nitrogen=1.0, sulfur=1.0, noble=1.0):
    """
    The carbon over oxygen atoms of the gas and of the solids at `temperature_k`, as `(gas, solid)` ([R3322]: the gas C/O rises outward past the water, CO2 and CO lines). The gas is what
    is not frozen: the carbon not in grains and not in a frozen CO, CO2 or CH4, over the oxygen not in rock and not in frozen water, CO or CO2. The solids are the rest. `None` for a
    side that holds no oxygen.
    """
    held = inventory(carbon, oxygen, nitrogen, sulfur, noble)
    frozen = condensed(temperature_k, carbon, oxygen, nitrogen, sulfur, noble)
    carbon_solid = frozen.get("C", 0.0) + frozen.get("CO", 0.0) + frozen.get("CO2", 0.0) + frozen.get("CH4", 0.0)
    oxygen_solid = held.rock_oxygen + frozen.get("H2O", 0.0) + frozen.get("CO", 0.0) + 2.0 * frozen.get("CO2", 0.0)
    carbon_gas = held.carbon - carbon_solid
    oxygen_gas = held.oxygen - oxygen_solid
    gas = carbon_gas / oxygen_gas if oxygen_gas > 1e-12 else None
    solid = carbon_solid / oxygen_solid if oxygen_solid > 1e-12 else None
    return gas, solid


@finite_domain()
def line_au(species, luminosity_w):
    """The distance, AU, at which a blackbody in the light of a star of `luminosity_w` is as cold as the line of `species` (`LINES_K`)."""
    if species not in LINES_K:
        raise ValueError(f"no line for {species!r}")
    distance_m = math.sqrt(luminosity_w / (16.0 * math.pi * physical_constants.STEFAN_BOLTZMANN_CONSTANT)) / LINES_K[species] ** 2
    return _m_to_au(distance_m)
