# planetgen/physics/sulfur_balance.py

"""
Where the degassed sulfur goes: the sulfate and sulfide balance (GEN.815; docs/design/composition/mixture-draw.md)
=================================================================================================================

The supply law (`physics/surface_pressure.py`) says how much sulfur the volcanoes have brought to the surface. Nothing in the air
holds it there: sulfur has two solid sinks that take it out of the gas as fast as it arrives, and the air keeps only what the sinks
have not yet taken and what they cannot take.

* **Sulfate.** Oxidized sulfur (SO2, SO3) is taken by calcium as anhydrite or gypsum, CaSO4. With liquid water the sulfur dioxide
  is scavenged by cloud and rain and oxidized to sulfate within days (the Earth's tropospheric lifetime is one to a few days);
  on a dry surface the gas reacts with the calcium-bearing rock, the sulfation of calcite and silicates that [R1343]
  put at about two million years for the SO2 of Venus.
* **Sulfide.** Reduced sulfur (H2S, S2, COS) is taken by iron as pyrite or troilite, FeS2 and FeS, on the same two lifetimes.
* **The crust that holds them.** Calcium and iron are supplied by the crust the same volcanism built: the mass of crust is the
  mass of mantle the interior has processed (`Theta` of `surface_pressure`) times the crust made per unit of processed mantle, the
  Earth twin's share. A sulfate sink holds at most one calcium per sulfur and a sulfide sink one iron per sulfur. Sulfur beyond
  that capacity has nowhere to go and stays in the air whole.

Which sink takes a mole of sulfur is the redox of the gas it arrived in: `oxidized_share` is the share of the sulfur that sits in
SO2 and SO3 in the equilibrium gas of the outgassing stage (`oxidized_share_of`), so a reducing mantle fills the sulfide sink and an
oxidizing one the sulfate sink.

The balance is a steady state, as the carbonate-silicate balance is. The sulfur that goes to the sinks arrived at the mean rate
`S / age`; with a removal lifetime `tau` the air holds `S tau / age` of it (all of it when `tau` exceeds the age). So the sulfur of the
air is the capacity excess plus that burden:

    air = (S - bound) + bound * min(1, tau / age)        with  bound = min(f S, C_sulfate) + min((1 - f) S, C_sulfide)

The balance conserves sulfur: air, sulfate and sulfide sum to `S`.

What this does not do: the speciation of the sulfur that stays in the air is the equilibrium chemistry's (the caller puts the sulfur
atoms into the closed system); the sulfur the ocean dissolves is `mixture_sinks` (Henry's law) and acts afterwards; calcium that the
carbonate balance has taken is not taken from the sulfate sink (sulfate and carbonate share the calcium, so the sulfate capacity is an
upper bound where a carbonate sink is full); the sulfur that stays in the melt is the degassing efficiency of `surface_pressure`.

Use:

    life = draw_lifetimes(constants)                                  # drawn, recorded in constants
    crust = crust_column_kg_m2(released_j, radius_m)
    result = balance(degassed_mol_m2, age_s=age_gy * 3.156e16, ocean_present=True, oxidized_share=0.9, crust_kg_m2=crust, lifetimes=life)
    result.air_mol_m2                                                 # the sulfur to put into the air's element budget

Every number that is not read from a source is a design constant with its source tier and tweak level in
`docs/design/research/gap-filling-plan.md` (GEN.725 rule).
"""

import math
from dataclasses import dataclass

from planetgen.physics import elements
from planetgen.physics.mixture_base import bracket
from planetgen.physics.surface_pressure import PROCESSING_ENERGY_J_PER_KG, THETA_EARTH

DAY_S = 86400.0
"""float: One day, s."""

YEAR_S = 3.15576e7
"""float: One Julian year, s."""

OXIDIZED_SULFUR = ("SO2", "SO3", "SO")
"""tuple: The gas-table species whose sulfur is bonded to oxygen (sulfur at +2 or more): the supply of the sulfate sink."""

REDUCED_SULFUR = ("H2S", "S2", "S", "OCS", "CS2")
"""tuple: The gas-table species whose sulfur is not bonded to oxygen: the supply of the sulfide sink."""

WET_LIFETIME_DAYS = (0.5, 7.0, "log")
"""tuple: The lifetime of atmospheric sulfur gas against rain-out and oxidation to sulfate over liquid water, days: the Earth's tropospheric
SO2 and H2S live one to a few days, the sulfate aerosol up to a week (read from memory, not re-checked: tier 3, tweak level 0)."""

DRY_LIFETIME_YEARS = (2.0e6, 2.0e7, "log")
"""tuple: The lifetime of atmospheric SO2 against sulfation of the surface rock when there is no liquid water, years, at `SULFATION_REFERENCE_K` (Venus's surface): about 1.9 million years from the
sulfation rate of [R1343] and 20 million from a later estimate, both second-hand through [R0157] (GEN.1562 N1: the range was 0.5 to 5 Myr, which
sat below the second figure). Confident unverified, tier 2, tweak level 0."""

SULFATION_ACTIVATION_K = 15248.0
"""float: The activation temperature of the sulfation of rock by SO2, K: the rate `10^19.64 exp(-15248 / T)` molecules per square centimetre per second of [R1343] (126.8 +/- 5.4 kJ/mol),
second-hand and checked only for consistency (15248 K x 8.314 = 126.8 kJ/mol; 4.4e10 at 735 K against the quoted 4.6e10 mean). The lab runs were 873 to 1123 K; below that the law is an
extrapolation, and below about 450 K another route (water films on the grains) would take over, so a cool dry surface is taken to sulfate slowly, not faster (GEN.1562 N2). Confident
unverified, tier 2, tweak level 0."""

SULFATION_REFERENCE_K = 735.0
"""float: The surface temperature of Venus, K (measured), where `DRY_LIFETIME_YEARS` is the lifetime: the lifetime at a temperature `T` is that times `exp(SULFATION_ACTIVATION_K (1 / T - 1 / 735))`."""

WATER_LAYER_FULL_M = 2.75
"""float: The depth, m, of a global layer of water (1e-3 of the Earth's ocean, 2.75 km deep) at and above which the sulfur is taken to be removed by water as fast as the wet lifetime says, ocean or not
([R2735], whose models need more than that to hold a sustained sulfuric acid haze) (GEN.1562 N6). Under it the water removes at most one mole of sulfur per mole of water (SO2 + H2O to
H2SO4 is one water per sulfur) and the rest meets only the dry sulfation lifetime. Confident unverified, tier 2, tweak level 1."""

WATER_MOLAR_MASS_KG = 0.01801528
"""float: Kilograms per mole of water (measured: the atomic masses of H and O)."""

EARTH_MANTLE_KG = 4.0e24
"""float: The mass of the Earth's mantle, kg (the mass `surface_pressure.PROCESSING_ENERGY_J_PER_KG` is calibrated on)."""

EARTH_CRUST_KG = 2.8e22
"""float: The mass of the Earth's crust, oceanic and continental, kg (about 0.47 percent of the Earth's mass; read from memory, not re-checked)."""

CRUST_PER_PROCESSED_MANTLE = EARTH_CRUST_KG / (EARTH_MANTLE_KG * THETA_EARTH)
"""float: Kilograms of crust made per kilogram of mantle processed through melting: the Earth twin makes `EARTH_CRUST_KG` out of
`THETA_EARTH` mantle masses (a single calibration point: tier 6, tweak level 0)."""

CRUST_CALCIUM_WT = 0.085
"""float: The calcium mass fraction of the crust the volcanism builds: the basalt of mid-ocean ridges, 11.9 wt% CaO ([R4687], read
from memory, not re-checked: tier 3, tweak level 0)."""

CRUST_IRON_WT = 0.074
"""float: The iron mass fraction of the crust the volcanism builds: the same basalt, 9.5 wt% FeO ([R4687], from memory: tier 3)."""


@dataclass(frozen=True)
class Lifetimes:
    """
    The removal lifetimes of sulfur gas in the air.

    Attributes:
        wet_s (float): With liquid water at the surface, s.
        dry_s (float): With none, s.
    """

    wet_s: float
    dry_s: float


@dataclass(frozen=True)
class SulfurBalance:
    """
    Where the degassed sulfur is, mol per square metre of surface. The three sum to the sulfur that was degassed.

    Attributes:
        air_mol_m2 (float): In the air: the capacity excess plus the steady-state burden.
        sulfate_mol_m2 (float): Bound as sulfate in the crust.
        sulfide_mol_m2 (float): Bound as sulfide in the crust.
        lifetime_s (float): The removal lifetime that was used.
    """

    air_mol_m2: float
    sulfate_mol_m2: float
    sulfide_mol_m2: float
    lifetime_s: float


def sulfation_lifetime_s(reference_s, temperature_k):
    """
    The sulfation lifetime of SO2 at a surface temperature (GEN.1562 N2): the lifetime at `SULFATION_REFERENCE_K` times `exp(SULFATION_ACTIVATION_K (1 / T - 1 / T_ref))`, so a cool dry surface
    keeps its SO2 far longer than Venus does (a factor 5e5 at 450 K) and a hotter one shorter. The exponent is held under 60 (the lifetime is far past any age by then). `temperature_k` of
    `None` or under 1 K gives `reference_s`.
    """
    if temperature_k is None or temperature_k < 1.0:
        return reference_s
    exponent = SULFATION_ACTIVATION_K * (1.0 / temperature_k - 1.0 / SULFATION_REFERENCE_K)
    return reference_s * math.exp(min(exponent, 60.0))


def draw_lifetimes(constants):
    """The two lifetimes, drawn and recorded in `constants` (`sulfur_lifetime_wet_days`, `sulfur_lifetime_dry_years`)."""
    wet = bracket(constants, "sulfur_lifetime_wet_days", WET_LIFETIME_DAYS) * DAY_S
    dry = bracket(constants, "sulfur_lifetime_dry_years", DRY_LIFETIME_YEARS) * YEAR_S
    return Lifetimes(wet_s=wet, dry_s=dry)


def oxidized_share_of(mole_fractions):
    """
    The share of the sulfur of a gas that is bonded to oxygen.

    Args:
        mole_fractions (dict): Gas-table formula -> mole fraction or moles; the sulfur gases of `OXIDIZED_SULFUR` and `REDUCED_SULFUR`
            are read (SO2 holds one sulfur, S2 and CS2 their own count).

    Returns:
        float: Between 0 and 1; 1 when the gas holds no sulfur at all (nothing to sink either way).
    """
    atoms = {"SO2": 1, "SO3": 1, "SO": 1, "H2S": 1, "S2": 2, "S": 1, "OCS": 1, "CS2": 1}
    oxidized = math.fsum(mole_fractions.get(name, 0.0) * atoms[name] for name in OXIDIZED_SULFUR)
    reduced = math.fsum(mole_fractions.get(name, 0.0) * atoms[name] for name in REDUCED_SULFUR)
    total = oxidized + reduced
    return oxidized / total if total > 0.0 else 1.0


def crust_column_kg_m2(released_j, radius_m):
    """
    The crust the interior has built, as a mass per square metre of the surface.

    Args:
        released_j (float): The heat the interior has carried out (`interior_heat.Heat.released_j`).
        radius_m (float): The body's radius.

    Returns:
        float: `CRUST_PER_PROCESSED_MANTLE * Theta * mantle / (4 pi R^2)` with `Theta = released / (mantle e_proc)`, which is
        `CRUST_PER_PROCESSED_MANTLE * released / (e_proc 4 pi R^2)` (the mantle's mass cancels), kg/m2; 0 for no heat.
    """
    if released_j <= 0.0 or radius_m <= 0.0:
        return 0.0
    processed_kg = released_j / PROCESSING_ENERGY_J_PER_KG
    return CRUST_PER_PROCESSED_MANTLE * processed_kg / (4.0 * math.pi * radius_m ** 2)


def sink_capacities_mol_m2(crust_kg_m2):
    """
    What a crust column can bind, mol of sulfur per square metre.

    Returns:
        tuple: `(sulfate, sulfide)`: one calcium per sulfur and one iron per sulfur.
    """
    return (crust_kg_m2 * CRUST_CALCIUM_WT / elements.molar_mass("Ca"), crust_kg_m2 * CRUST_IRON_WT / elements.molar_mass("Fe"))


def balance(degassed_mol_m2, *, age_s, ocean_present, oxidized_share, crust_kg_m2, lifetimes, temperature_k=None, water_layer_m=None):
    """
    Splits the degassed sulfur between the air, the sulfate sink and the sulfide sink.

    Args:
        degassed_mol_m2 (float): The sulfur the volcanoes have brought to the surface, mol per square metre.
        age_s (float): The body's age, s; the sulfur is taken to have arrived at its mean rate over it.
        ocean_present (bool): Whether liquid water is at the surface (picks the wet or the dry lifetime).
        oxidized_share (float): The share of the sulfur that is oxidized (`oxidized_share_of`), 0 to 1.
        crust_kg_m2 (float): The crust that holds the sinks (`crust_column_kg_m2`).
        lifetimes (Lifetimes): `draw_lifetimes`.
        temperature_k (float or None): The surface temperature, K, that scales the dry sulfation lifetime (`sulfation_lifetime_s`, GEN.1562 N2); `None` keeps the lifetime at Venus's.
        water_layer_m (float or None): The depth, m, of the global layer the body's surface water would make (GEN.1562 N6). Without an ocean, water under `WATER_LAYER_FULL_M` removes at most one mole of
            oxidized sulfur per mole of water over the wet lifetime and the rest meets the dry lifetime; at or over it the wet lifetime holds for all of it. `None` keeps the old rule (the ocean
            alone picks the lifetime).

    Returns:
        SulfurBalance: The three reservoirs; they sum to `degassed_mol_m2`.
    """
    total = max(degassed_mol_m2, 0.0)
    share = min(max(oxidized_share, 0.0), 1.0)
    capacity_sulfate, capacity_sulfide = sink_capacities_mol_m2(max(crust_kg_m2, 0.0))
    sulfate = min(share * total, capacity_sulfate)
    sulfide = min((1.0 - share) * total, capacity_sulfide)
    dry = sulfation_lifetime_s(lifetimes.dry_s, temperature_k)
    wet_full = ocean_present or (water_layer_m is not None and water_layer_m >= WATER_LAYER_FULL_M)
    if wet_full:
        wet_share = 1.0
    elif water_layer_m is None or water_layer_m <= 0.0 or sulfate <= 0.0:
        wet_share = 0.0
    else:       # GEN.1562 N6: acid is made from the water there is, one mole of oxidized sulfur per mole of water
        wet_share = min(1.0, max(water_layer_m, 0.0) * 1000.0 / WATER_MOLAR_MASS_KG / sulfate)
    wet_kept = 1.0 if age_s <= 0.0 else min(1.0, lifetimes.wet_s / age_s)
    dry_kept = 1.0 if age_s <= 0.0 else min(1.0, dry / age_s)
    burden_sulfate = sulfate * (wet_share * wet_kept + (1.0 - wet_share) * dry_kept)
    burden_sulfide = sulfide * (wet_kept if wet_full else dry_kept)
    lifetime = lifetimes.wet_s if wet_share >= 0.5 else dry
    excess = total - sulfate - sulfide
    return SulfurBalance(air_mol_m2=excess + burden_sulfate + burden_sulfide, sulfate_mol_m2=sulfate - burden_sulfate,
                         sulfide_mol_m2=sulfide - burden_sulfide, lifetime_s=lifetime)
