# planetgen/physics/disk_condensation.py

"""
What a protoplanetary disk condenses at a temperature, species by species (GEN.1363 part A; docs/design/composition/disk-condensation.md)
==========================================================================================================================================

`condensate_shares` has three ices (water 182 K, ammonia 131 K, methane 41 K) and puts all the carbon in methane. A disk measured by its
gas and ice is not like that [R3321, R3322]: about half of the carbon is refractory (organic and graphite-like grains that survive to about
500 K); the rest is mostly CO and CO2, with a little CH4; most of the nitrogen is N2; sulfur is mostly in troilite with some H2S; and the noble gases
come out as clathrates in the water ice. This module is the budget and the lines of that disk, and nothing draws from it yet (part B wires the
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
   methane (the equilibrium at 10^-4 bar, [R2726: Table 10]); CO and CO2 are the disk-density values of [R3322: Table 1], N2 and H2S are recalled
   (registered so). Noble-gas clathrate lines are [R2726: Table 8].
5. **The C/O across a line.** The gas C/O rises outward past each snow line (water takes oxygen first, then CO2, then CO); the solids' C/O is the
   other half. A giant that gathers gas between the water and CO lines is carbon-rich, one that gathers solids is oxygen-rich [R3322].

Limits. Equilibrium supply: every species freezes whole below its line. A body that formed warm and trapped a species in water ice (a comet) keeps a tiny
share of it (`comet_ice_formation`, `comet_ice_budget.TRAPPED_FRACTION`), and that is for the body, not the disk. The lines do not move as the disk cools. The H2S and the
ammonium-salt shares of the nitrogen and sulfur are stand-ins (registry rows DSC4 to DSC6).
"""

import dataclasses
import math

from planetgen import tuning
from planetgen.physics import condensate_shares, constants as physical_constants, unit_registry as units
from planetgen.util.checks import finite_domain

NOBLE_LOG_EPS = {"Kr": 3.12, "Xe": 2.22}
"""dict: The Sun's photospheric log eps (H = 12) of krypton and xenon ([R0168: Table 2]; argon is in `tuning.STELLAR_SOLAR_LOG_EPS`). [Stand-in DSC1, tier T1, tweak L0.]"""

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
clathrates are [R2726: Table 8] (Xe 69, Kr 53, Ar 48); CO 20 and CO2 47 are [R3322: Table 1] (CO2 42 to 52 K there, 70 K in other disks); N2 22 K is recalled
([R3321] section 1.2 puts CH4, CO and N2 near 20 to 40 K) and H2S 80 K is recalled (ice sublimates near 50 to 90 K by the matrix). [Stand-in DSC2, tier T1 for the read
lines and T3 for N2 and H2S, tweak L0.]"""

REFRACTORY_LINE_K = 500.0
"""float: The temperature below which the refractory carbon is a grain ([R3322: Table 1] adopt more than 150 K and 500 K in the model). [Stand-in DSC3, tier T2, tweak L0.]"""

REFRACTORY_CARBON_SHARE = 0.5
"""float: The share of the carbon that is refractory grains: about half ([R3321] section 1.2, from the interstellar medium; [R3322: Table 1] has 0.6 of 2.4 of the carbon
they account for in the disk). [Stand-in DSC3, tier T2, tweak L1.]"""

CO_PER_CO2 = 1.5 / 0.3
"""float: The CO over the CO2 of the volatile carbon that is not methane ([R3322: Table 1], CO 1.5 and CO2 0.3 molecules per 10^4 hydrogen). [Stand-in DSC3, tier T2, tweak L0.]"""

CH4_SHARE_OF_VOLATILE_CARBON = 0.34 / (0.34 + 3.1 + 4.7)
"""float: The methane share of the volatile carbon: comet 67P's CH4 over CH4 plus CO plus CO2, 0.34, 3.1 and 4.7 percent of the water ([R3835: Table 2]), 4.2 percent.
One comet, taken for the disk's methane. [Stand-in DSC4, tier T3, tweak L1.]"""

NITROGEN_AS_N2 = 0.6
NITROGEN_AS_SALTS = 0.1
"""float: The shares of the nitrogen that are N2 and ammonium salts, the rest being ammonia hydrate: most nitrogen is N2 and unobserved ([R3321] section 1.2) and the
salts hold a minority ([R3321], citing the Rosetta result second-hand; the split is recalled). [Stand-ins DSC5, tier T3, tweak L3.]"""

SULFUR_AS_H2S = 0.3
"""float: The share of the sulfur that is H2S ice: 67P's H2S is 1.1 percent of the water ([R3835: Table 2]) against 3.5 percent for all the sulfur over the water (S 0.13 over the
3.8 oxygen the rock leaves). The rest is troilite in the rock. [Stand-in DSC6, tier T3, tweak L2.]"""

RECIPE_C_OVER_O_CENTER = 0.57
"""float: The center of the core of the giant recipe's C/O row `f7` (0.57, the median of LExACoM's 69 retrieved giants, pinned by tests/test_atmosphere_recipes.py)."""

RECIPE_C_OVER_O_CORE_SIGMA_DEX = 0.21
RECIPE_C_OVER_O_CORE_RANGE = (0.2, 0.8)
RECIPE_C_OVER_O_TAIL_RANGE = (0.8, 1.5)
RECIPE_C_OVER_O_TAIL_SHARE = math.sqrt(0.043 * 0.087)
"""The rest of the shape of `f7` (the log-normal core cut to 0.2 to 0.8, the uniform tail to 1.5, the tail share log-uniform from 0.043 to 0.087; tests/test_atmosphere_recipes.py pins them
against the recipe file): `recipe_c_over_o_median` works out the median of the draw from them."""


def recipe_c_over_o_median():
    """
    The median of the drawn C/O row `f7`: the core is cut at 0.8, above its center's 0.57 by only 0.7 sigma, so the median of the draw is 0.48 and not the core's 0.57 (GEN.1471; GEN.1363
    part C had divided by 0.57 and so lowered every giant by a sixth). Found by bisection on the mixture's distribution function at the geometric middle of the tail share.
    """
    center = math.log10(RECIPE_C_OVER_O_CENTER)
    low, high = (math.log10(x) for x in RECIPE_C_OVER_O_CORE_RANGE)

    def normal(x):
        return 0.5 * (1.0 + math.erf((x - center) / (RECIPE_C_OVER_O_CORE_SIGMA_DEX * math.sqrt(2.0))))

    inside = normal(high) - normal(low)
    tail_low, tail_high = RECIPE_C_OVER_O_TAIL_RANGE

    def cumulative(x):
        core = (normal(min(max(math.log10(x), low), high)) - normal(low)) / inside
        tail = min(max((x - tail_low) / (tail_high - tail_low), 0.0), 1.0)
        return (1.0 - RECIPE_C_OVER_O_TAIL_SHARE) * core + RECIPE_C_OVER_O_TAIL_SHARE * tail

    lo, hi = RECIPE_C_OVER_O_CORE_RANGE[0], tail_high
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if cumulative(mid) < 0.5 else (lo, mid)
    return 0.5 * (lo + hi)


RECIPE_C_OVER_O_MEDIAN = recipe_c_over_o_median()
"""float: The median of the drawn `f7`, 0.477: the drawn `f7` over this is the spread of a giant's C/O around its disk's value, which the disk alone (one place, one star) does not give, without
moving the center (GEN.1363 part C, stand-in DSC7, tier T3, tweak L2; the center was 0.57 before GEN.1471)."""

RECIPE_C_OVER_O_SPREAD_EXPONENT = 0.5
"""float: The power of that ratio that multiplies the disk's C/O. The recipe's spread (the history of what a giant swallowed) and the disk's own (the star's C/O, the place of formation)
overlap, and the full ratio counted both (GEN.1471: 16 percent of the giants under 0.35 against the retrieved 17 with the half power; the full power puts 20 to 26 under it at the same median). [Stand-in DSC8, tier T6, tweak L2.]"""

ENVELOPE_METAL_SHARE = 0.06
"""float: The share of a giant's bulk heavy-element enrichment (`f5`, the planet's whole heavy-element mass over its star's, [R4352]) that its envelope carries when its C/O is set. The bulk
metallicity is an upper limit on the atmosphere's ([R4350: Abstract]) and the water of gas giants is substellar ([R4663: Abstract]); the whole bulk enrichment, 18 times the star's for the
median giant, is 90 percent solids whose C/O is 0.44 whatever their amount. At 0.06 (a median envelope enrichment of 1.1) the generated median is 0.57 and 16 percent of the giants sit under 0.35,
the LExACoM figures (0.57 and 17). [Stand-in DSC9, tier T6, tweak L2.]"""

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
    The molecules per 10^4 hydrogen that a disk has to condense, by species, and the oxygen and carbon the rock and the grains hold.

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
    The `Inventory` of a disk whose star has these element ratios (each over the rock-forming metals, over the Sun's; 1 is the Sun).

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
def water_share_of_solids(temperature_k, carbon=1.0, oxygen=1.0, nitrogen=1.0, sulfur=1.0, noble=1.0):
    """
    Water ice over everything solid (`masses_percent`: rock, the carbon grains and every ice) at `temperature_k`; 0 above the water line. The oxygen that CO and CO2 hold is out of the
    water, so the Sun's water share is 0.34 against the 0.54 of Lodders's equilibrium with all the spare oxygen as water (`condensate_shares`, before GEN.1363 part B).
    """
    masses = masses_percent(temperature_k, carbon, oxygen, nitrogen, sulfur, noble)
    return masses.get("H2O", 0.0) / sum(masses.values())


SNOW_LINE_SOLID_JUMP_K = 0.5 * (LINES_K["H2O"] + LINES_K["NH3"])
"""float: The temperature at which the jump of the solid density at the snow line is read, K: the middle of the water and ammonia lines, where the water is frozen and the ammonia is not."""


@finite_domain()
def snow_line_solid_jump(carbon=1.0, oxygen=1.0, nitrogen=1.0, sulfur=1.0, noble=1.0):
    """
    How many times the rock's surface density the solids make just beyond the water line (GEN.1396): everything condensed at
    `SNOW_LINE_SOLID_JUMP_K` (the rock, the refractory carbon grains and the water ice, from `masses_percent`) over the rock alone. The divisor is
    the rock because the density scale of the solids is the rock's: the minimum-mass nebula is built from the rock of the planets ([R1854] puts the
    jump at 30 / 7), and the planets of the inner disk hold far less carbon than the grains would bring (meteorites and the bulk silicate Earth have
    a C/Si well under the Sun's, so the grains are destroyed there [R0357]), so the grains are a part of the solids only where they survive.
    Lodders's equilibrium table gave 2.17 (all the spare oxygen as water); the disk holds part of that oxygen in CO and CO2 gas, which is not frozen
    here, and the grains add 0.145 percent of the nebula, so the Sun's value is 1.97.

    Args:
        carbon, oxygen, nitrogen, sulfur, noble (float): The star's ratios of the element over the rock-forming metals, over the Sun's (1 is the Sun).

    Returns:
        float: The jump, 1.97 for the Sun.
    """
    masses = masses_percent(SNOW_LINE_SOLID_JUMP_K, carbon, oxygen, nitrogen, sulfur, noble)
    return sum(masses.values()) / masses[condensate_shares.ROCK]


SNOW_LINE_SOLID_JUMP = snow_line_solid_jump()
"""float: The jump of the solid surface density at the snow line for the Sun's composition (1.97); `formation.mmsn_surface_density_gcm2` multiplies the disk's solids by it beyond a star's own snow line."""


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


def _gas_and_solid_atoms(temperature_k, carbon, oxygen, nitrogen, sulfur, noble, ice_weights=None):
    """
    `(carbon gas, oxygen gas, carbon solid, oxygen solid)` in atoms per 10^4 hydrogen at `temperature_k` (the rock's oxygen is solid at any temperature below its line).

    `ice_weights` (`{species: E}`, from `disk_drift.ice_weights`; GEN.1473) multiplies the vapor of each ice that evaporated on the way to this place: its atoms in the gas are `E` times
    the static ones, the solids stay the static ones. `None` and weights of 1 are the static budget (GEN.1363) exactly.
    """
    held = inventory(carbon, oxygen, nitrogen, sulfur, noble)
    frozen = condensed(temperature_k, carbon, oxygen, nitrogen, sulfur, noble)
    carbon_solid = frozen.get("C", 0.0) + frozen.get("CO", 0.0) + frozen.get("CO2", 0.0) + frozen.get("CH4", 0.0)
    oxygen_solid = held.rock_oxygen + frozen.get("H2O", 0.0) + frozen.get("CO", 0.0) + 2.0 * frozen.get("CO2", 0.0)
    carbon_gas = held.carbon - carbon_solid
    oxygen_gas = held.oxygen - oxygen_solid
    for species, weight in (ice_weights or {}).items():
        if species not in held.moles or temperature_k < LINES_K[species]:
            continue                                  # a frozen ice is not in the gas, whatever its weight
        extra = (weight - 1.0) * held.moles[species]
        carbon_gas += extra * ICE_CARBON.get(species, 0)
        oxygen_gas += extra * ICE_OXYGEN.get(species, 0)
    return carbon_gas, oxygen_gas, carbon_solid, oxygen_solid


ICE_CARBON = {"CO": 1, "CO2": 1, "CH4": 1}
"""dict: The carbon atoms in a molecule of each ice (`_gas_and_solid_atoms` weights the vapor by these)."""

ICE_OXYGEN = {"H2O": 1, "CO": 1, "CO2": 2}
"""dict: The oxygen atoms in a molecule of each ice."""


def _planetesimal_share(enrichment, carbon_gas, oxygen_gas, carbon_solid, oxygen_solid, metals_star=None):
    """`sigma` of `envelope_c_over_o`: the multiple of the frozen solids the envelope holds, that makes its carbon and oxygen `enrichment` times the star's (`metals_star`, the budget's carbon and oxygen; the sum of the four when not given)."""
    if metals_star is None:
        metals_star = carbon_gas + oxygen_gas + carbon_solid + oxygen_solid
    metals_solid = carbon_solid + oxygen_solid
    return max(enrichment * metals_star - (carbon_gas + oxygen_gas), 0.0) / metals_solid if metals_solid > 0.0 else 0.0


@finite_domain()
def envelope_c_over_o(temperature_k, enrichment, carbon=1.0, oxygen=1.0, ice_weights=None):
    """
    The carbon over oxygen atoms of a giant's envelope that formed where the disk was at `temperature_k` and holds `enrichment` times its star's heavy elements (GEN.1363 part C).

    The envelope is the disk's gas plus `sigma` times the solids that froze there, `sigma` being what makes the carbon and oxygen of the envelope `enrichment` times the star's
    (`sigma` is 1 at enrichment 1: gas and solids in the disk's own proportion, the star's C/O; it is 0, gas alone, when the enrichment is under the gas's share of the star's carbon
    and oxygen). A cold or inner-disk giant with little enrichment takes the gas C/O of its place ([R3322]: it rises outward past the water line), a heavily enriched one the
    solids' (oxygen-rich) C/O. `None` when the disk has no oxygen left in the gas and the envelope holds none of the solids' either.

    Args:
        temperature_k (float): The disk temperature where the giant formed.
        enrichment (float): The envelope's heavy elements over its star's (`f5` of the giant recipe); zero or more.
        carbon, oxygen (float): The star's C and O over its rock-forming metals, over the Sun's (`bulk_draw.ice_forming_ratios`).
        ice_weights (dict): `{species: E}` of the ices that evaporated on the way to the giant's place at the age it gathered its gas (`disk_drift.giant_ice_weights`, GEN.1473: the pebbles
            carry ice inward faster than the gas); `None` for the static disk, which is every `E` equal to 1.
    """
    if not temperature_k > 0.0:
        raise ValueError("the temperature must be positive")
    if not enrichment >= 0.0:
        raise ValueError("the enrichment must be zero or more")
    carbon_gas, oxygen_gas, carbon_solid, oxygen_solid = _gas_and_solid_atoms(temperature_k, carbon, oxygen, 1.0, 1.0, 1.0, ice_weights)
    held = inventory(carbon, oxygen)
    sigma = _planetesimal_share(enrichment, carbon_gas, oxygen_gas, carbon_solid, oxygen_solid, held.carbon + held.oxygen)
    oxygen_total = oxygen_gas + sigma * oxygen_solid
    if oxygen_total <= 1e-12:
        return None
    return (carbon_gas + sigma * carbon_solid) / oxygen_total


ENRICHED_ELEMENTS = ("N", "S", "Ar", "Kr", "Xe")
"""tuple: The elements `envelope_enrichments` gives their own factor beyond C and O (GEN.1372). P is not among them: its carriers (phosphides, apatite) condense above 1000 K, so it
is rock at every formation temperature of a giant and follows the rock-forming metals' scale. Ne stays gas and is not enriched by solids."""


@finite_domain()
def envelope_enrichments(temperature_k, enrichment, carbon=1.0, oxygen=1.0, nitrogen=1.0, sulfur=1.0, noble=1.0, ice_weights=None):
    """
    The nitrogen, sulfur and noble gases of a giant's envelope over its star's, element by element (GEN.1372), for a giant that formed where the disk was at `temperature_k` and holds
    `enrichment` times its star's heavy elements.

    The envelope is the disk's gas plus `sigma` times what froze there (`envelope_c_over_o` gives `sigma`: the share that makes the carbon and oxygen `enrichment` times the star's),
    so an element's factor is `(gas + sigma * solid) / (gas + solid)` of its atoms: one for an element that stays gas (the noble gases of a giant that formed warmer than their
    lines, unless the pebbles carried them in: `ice_weights`, GEN.1516), `sigma` for one that froze whole as an ice, the rock's own scale `enrichment` for the troilite of the sulfur, and
    between for the nitrogen (N2 gas, the ammonia and the salts frozen below the ammonia line) and the rest of the sulfur (H2S gas above its line).

    Args:
        temperature_k (float): The disk temperature where the giant formed.
        enrichment (float): The envelope's heavy elements over its star's; zero or more.
        carbon, oxygen, nitrogen, sulfur, noble (float): The star's ratios over the rock-forming metals, over the Sun's (`inventory`); `noble` is the argon's, taken for krypton and xenon.
        ice_weights (dict): `{species: E}` of the ices that evaporated on the way to the giant's place (`disk_drift.giant_ice_weights`, GEN.1473; GEN.1516 gives it to the nitrogen,
            sulfur and noble gases as `envelope_c_over_o` gives it to the carbon and oxygen): the gas holds `E` times the static vapor of each, so a noble gas the pebbles carried in is enriched
            past what the solids alone give. `None` is the static disk, every `E` equal to 1.

    Returns:
        dict: `ENRICHED_ELEMENTS` symbol to the envelope's atoms per hydrogen over the star's.
    """
    if not temperature_k > 0.0:
        raise ValueError("the temperature must be positive")
    if not enrichment >= 0.0:
        raise ValueError("the enrichment must be zero or more")
    held = inventory(carbon, oxygen, nitrogen, sulfur, noble)
    frozen = condensed(temperature_k, carbon, oxygen, nitrogen, sulfur, noble)
    if ice_weights:
        sigma = _planetesimal_share(enrichment, *_gas_and_solid_atoms(temperature_k, carbon, oxygen, nitrogen, sulfur, noble, ice_weights), held.carbon + held.oxygen)
    else:
        sigma = _planetesimal_share(enrichment, *_gas_and_solid_atoms(temperature_k, carbon, oxygen, nitrogen, sulfur, noble))
    extra = {species: (weight - 1.0) * held.moles[species] for species, weight in (ice_weights or {}).items() if species in held.moles and temperature_k >= LINES_K[species]}
    nitrogen_total = 2.0 * held.moles["N2"] + held.moles["NH3"] + ammonium_salts(carbon, oxygen, nitrogen, sulfur, noble)
    nitrogen_solid = 2.0 * frozen.get("N2", 0.0) + frozen.get("NH3", 0.0) + frozen.get("NH4", 0.0)
    nitrogen_vapor = 2.0 * extra.get("N2", 0.0) + extra.get("NH3", 0.0)
    sulfur_total = solar_per_1e4_h("S") * sulfur
    sulfur_ice = frozen.get("H2S", 0.0)
    troilite = sulfur_total * (1.0 - SULFUR_AS_H2S)
    solids = {"N": (nitrogen_total, nitrogen_solid, nitrogen_vapor), "Ar": (held.moles["Ar"], frozen.get("Ar", 0.0), extra.get("Ar", 0.0)),
              "Kr": (held.moles["Kr"], frozen.get("Kr", 0.0), extra.get("Kr", 0.0)), "Xe": (held.moles["Xe"], frozen.get("Xe", 0.0), extra.get("Xe", 0.0))}
    factors = {symbol: (total - solid + vapor + sigma * solid) / total for symbol, (total, solid, vapor) in solids.items()}
    factors["S"] = (sulfur_total - troilite - sulfur_ice + extra.get("H2S", 0.0) + enrichment * troilite + sigma * sulfur_ice) / sulfur_total      # the troilite is rock: the rock's scale, not the ices'
    return factors
