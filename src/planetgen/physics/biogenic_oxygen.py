# planetgen/physics/biogenic_oxygen.py

"""
Free oxygen from oxygenic photosynthesis (GEN.855, docs/design/habitability/life-timeline-and-oxygen.md)
=========================================================================================================

A world whose life splits water for its electrons (a life class of `physics/life_classes.py` whose `metabolism` is `"oxygenic"`) fills its air with oxygen; one whose life does
not, does not. The oxygen a world holds is the balance of a source and its sinks, in moles of O2 a year:

* **Source** `S`: the organic carbon its life buries, which leaves the oxygen of the water it was made from behind. It goes with the
  area where light, water and warmth meet (`productivity`): the planet's area, the share the ocean covers, the square root of the
  light of its band it gets (class 1a's photosynthetically active 400 to 700 nm, 1b's blue and 1c's far red, GEN.1744; the production saturates), and a Boltzmann factor for the temperature of the surface, to zero
  outside the liquid range of the sea.
* **Reductant flux** `R`: the volcanic and metamorphic gases (H2, CO, H2S) that take oxygen as they are released, which goes with the
  heat the interior gives off and with how reduced the mantle is (`reductant_flux`).
* **Oxidative weathering** `W`: the rocks at the surface (kerogen, pyrite, ferrous iron) take oxygen at a rate per area of land that goes with the oxygen's
  partial pressure `p` to a power, `W = W0 (p / p0)^m` per area with `m = OXIDATIVE_WEATHERING_EXPONENT` (GEN.1733; linear, the top of the range the sources give: the square root of the COPSE model [R5280] and the fit of [R5281] to a linear form),
  where `W0` and `p0` are Earth's (`EARTH_SURPLUS_MOL_YR`, `EARTH_O2_PA`). The weathering depends on the pressure, not the moles, because the rock meets the air's oxygen
  at its partial pressure: a heavier world with the same oxygen per area has the larger pressure and weathers faster, so its pressure is no larger.
* **A store of reductants** `B`: the ocean's ferrous iron and the crust's sulfides and iron, and the reduced gases of the air itself (H2, CO,
  CH4, NH3, H2S, C2H6: `REDUCED_GASES`), take the first oxygen the surplus makes, so the air holds none until `B / (S - R)` has passed
  (`Profile.delay_gyr`, the time the Great Oxidation took after the first oxygenic photosynthesis). The reduced gases are burned away by
  then: `apply` turns them into their oxides (CO2, N2) and the water they make goes to the ocean.

`dN/dt = S - R - W(p)` once the store is spent, so the air settles to the pressure `p* = p0 ((S - R) A0 / (A (S0 - R0)))^(1/m)` (`balance_pa`: the surplus per area over
Earth's, to the power 1/m) in a few million years, a time short beside the Gyr the schedule runs in, and the world holds `p*` (as the moles that make that partial pressure
among the other gases of its air) when its photosynthesis is older than the delay and nothing before. A world with `S <= R` never holds any. The most the air can hold is the oxygen pressure `mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA`, the
fire limit (GEN.1320: a pressure, so a thin air can reach a larger share than a thick one; never more than `mixture_oxygen.ABIOTIC_O2_MAX` of the molecules).
The whole air, oxygen and all, stays under the pressure the retention regime of the air solve lets the body keep (`surface_air.Air.retention_cap_pa`, the row `d6`; GEN.1340): a body that
holds only 1 or 10 kPa against the escape of its air does not hold 40 kPa because its life made oxygen, and the shares of its gases are kept when the cap takes the excess.

The numbers are those of the one world we know (`EARTH_*`, round figures from the literature of the oxygen cycle: [R0689, R1980], from memory), used as the test point and the scale and not as a preset: a world
with another light, water, heat, redox and gravity gets another oxygen. Every constant is a stand-in in section 6.21 of
`docs/design/research/gap-filling-plan.md` with its tier and tweak level.

A world with no wide sunlit sea but with light and liquid water at its surface (a brine seep, a seasonal melt, a spring) has the small source of `refuge_productivity` (GEN.1723, `physics/surface_refuge.py`).

What it leaves out: the loss of oxygen to space (a few parts in 1e4 of the inventory a Gyr for a habitable-zone world, far under the
weathering), the hydrogen escape that adds oxygen, the phosphorus feedback and, for a world with a wide sea, the land plants that raised Earth's burial (they would
change `S` by a factor of a few, inside the 1.7 to 4 that this module's L3 tweak already stands for; a world with no wide sea has them from the rooted stage on, GEN.1729, in `refuge_productivity`).
"""

import collections
import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.generation import abundance_tier
from planetgen.physics import aerosols, atmosphere, atmosphere_state, constants, envelope_share, gases, habitability_world, mixture, mixture_oxygen, stellar_spectrum, surface_air, surface_refuge
from planetgen.util import log

EARTH_O2_PA = 21.2e3
"""float: The oxygen partial pressure of Earth's air, Pa (20.95 percent of 101.325 kPa), the point the oxidative weathering is scaled from (`p0`). [BO14, T1, L0.]"""

EARTH_BURIAL_MOL_YR = 1.0e13
"""float: The organic carbon Earth buries a year, mol (the source `S`; [R1859], about 1e13). T2, L3."""

EARTH_REDUCTANT_MOL_YR = 3.0e12
"""float: The oxygen Earth's volcanic and metamorphic reductants take a year, mol (`R`; [R1980] few 1e12). T3, L3."""

LAND_SHARE_OF_BURIAL = 0.27
"""float: The share of Earth's organic carbon burial that land plants carry (GEN.1729): they add 25 to 50 percent to the burial ([R5283], C/P 1000 to 2000), so 0.20 to 0.33 of today's total; the middle. Used only for a
world with no wide sea, since `EARTH_BURIAL_MOL_YR` of a sea world already has the plants inside it. One model, so tier 3 (the sources of its flux, [R5284], are images on the page; the land share of burial in [R5285] was not read). [BO15, T3, L0.]"""

LAND_SHARE_OF_BURIAL_RANGE = (0.20, 0.33)
"""tuple: The range of `LAND_SHARE_OF_BURIAL` ([R5283]: 25 percent at C/P 1000 to 50 percent at C/P 2000 over the burial without plants). [BO15, T3, L0.]"""

EARTH_STORE_MOL = 4.0e21
"""float: The oxygen the ocean's iron and the crust's sulfides and iron take before the air holds any, mol (`B`), chosen so that the first
oxygen follows the first photosynthesis by 0.7 Gyr when `S` is Earth's and `R` is Earth's at that age, from the model's own heat history of an
Earth-mass body (GEN.1463: 2.2 times today's heat when the biomarkers appear, 2.7 Ga, falling to 1.9 times at the rise; it was 5e21 while `R` was held
at today's value), the Great Oxidation coming some 0.3 to 1 Gyr after oxygenic photosynthesis began. T3, L1."""

RISE_STEP_GYR = 0.25
"""float: The step, Gyr, of the integral of the surplus `S - R(age)` in `rise_age_gyr`; the heat flow it reads changes by under 10 percent over a step after the first Gyr."""

EARTH_AREA_M2 = 4.0 * math.pi * constants.EARTH_RADIUS_M ** 2
EARTH_OCEAN_FRACTION = 0.71
EARTH_LAND_FRACTION = 1.0 - EARTH_OCEAN_FRACTION
"""float: The share of Earth's surface that is land, 0.29 (the land the plants' share of the burial is over; with the ocean's 0.71 it makes the surface)."""
EARTH_HEAT_W = 4.4e13
"""float: Earth's own heat flow now, W (44 TW)."""

EARTH_MANTLE_DELTA_IW = 3.5
"""float: The redox of Earth's upper mantle now, log10 units above the iron-wustite buffer (the scale `oxygen_fugacity` gives `mantle_delta_iw`)."""

EARTH_SURPLUS_MOL_YR = EARTH_BURIAL_MOL_YR - EARTH_REDUCTANT_MOL_YR
"""float: The surplus of Earth's oxygen budget, mol a year, which Earth's oxidative weathering takes at `EARTH_O2_PA` (`W0`)."""

OXIDATIVE_WEATHERING_EXPONENT = 1.0
"""float: The power of the oxygen partial pressure that oxidative weathering goes with: linear, the top of the sourced range. The baseline of the COPSE model [R5280] (taken from [R5282]) is the square root, which the weathering
model of [R5281] fits to within its own scatter ("consistent with coal oxidation kinetics"); a linear dependency has been used, and the evidence for any dependency has been disputed [R5280]. At the square root the balance pressure
goes with the surplus per area squared and the 16 to 35.5 kPa window of the habitable classes is a surplus of 0.87 to 1.29 of Earth's, too narrow: measured on 400 Abundant systems the oxygen-making sea worlds inside the window were 7 at 0.5
against 24 at 1.0 (11 under the old moles law). At 1.0 the window is 0.75 to 1.67 of Earth's surplus. [BO13, T3, L1: the end of its range, chosen for the width of the window.]"""

OXIDATIVE_WEATHERING_EXPONENT_RANGE = (0.5, 1.0)
"""tuple: The range of `OXIDATIVE_WEATHERING_EXPONENT` of the sources that give a dependence (square root to linear). [BO13, T3, L1.]"""

REDOX_DECADE_PER_UNIT = 0.25
REDOX_FACTOR_RANGE = (0.2, 5.0)
"""A mantle one log unit more reduced than Earth's puts out `10 ** 0.25` times the reductants and holds that times the store, kept inside the range."""

LIGHT_EXPONENT = 0.5
"""float: Production goes as the active light to this power: it saturates (T6, L3)."""

LIGHT_RANGE = (0.0, 2.0)
PAR_NM = (400.0, 700.0)
"""tuple[float]: The photosynthetically active band, nm: class 1a's (chlorophyll a), the band Earth's production is counted in."""

Band = collections.namedtuple("Band", ("nm", "efficiency"))
"""namedtuple: The light an oxygenic class eats: its band `nm` (low, high) and the `efficiency` of using it against class 1a's (`life_classes.json` `energy`: 1b blue 300 to 500 nm at 1.0, 1c far red 700 to 800 nm at 0.8; GEN.1744)."""

PAR_BAND = Band(PAR_NM, 1.0)
"""Band: Class 1a's, and the band of Earth's production: the one a body's production is counted in when no class is named."""

SUN_TEMPERATURE_K = constants.SOLAR_TEMPERATURE_K

TEMPERATURE_ACTIVATION_EV = 0.32
"""float: The activation energy of photosynthesis, eV (about a half of a heterotroph's; T2, L0)."""

REFERENCE_TEMPERATURE_K = 288.0
OCEAN_SPAN_K = (258.0, 275.0, 315.0, 343.0)
"""tuple[float]: The temperatures, K, at which the factor for the sea's liquid range rises from 0 to 1 and falls to 0 again (T6, L3)."""

LIVING_HYDROSPHERES = ("surface ocean", "hycean")
"""tuple: The states with a sunlit sea for the photosynthesis to stand in (a `brine seep` is not one, GEN.908; its life, and a melt's or a spring's, is a refuge of `surface_refuge`, GEN.1723)."""

MINIMUM_PA = 1.0
"""float: Under this partial pressure, Pa, a world is counted as making no oxygen."""

EUKARYOTE_MINIMUM_PA = 200.0
"""float: The oxygen, Pa, the air must hold for the wait for complex cells to start (about 1 percent of Earth's: [R3084]; T3, L0)."""

PLATEAU_PA = 600.0
"""float: The oxygen partial pressure of the plateau, Pa: where the oxygen stands after the first rise (the reductant sink `R` under the source `S`) until the second (GEN.1707). The models of Earth's
plateau put it at 0.06 to 0.25 of today's oxygen, where the proxies say at most 0.001 (`docs/design/research/` oxygen plateau report: 0.001 to 0.25, 0.02 to 5 kPa); 0.03 of Earth's 21 kPa is 0.6 kPa. Tier 6 (educated guess from physics), tweak level 1."""

SECOND_RISE_RATIO = 3.0
"""float: The oxygen rises again to the balance level once the reductant sink `R` is under this share of the source `S`, 1 over the ratio (GEN.1707): the deep ocean and the shelf stop recycling phosphorus
and the oxygen gets a second jump. A single model's number (the reductant input crossing 1.5e13 mol O2-equivalent a year for the first rise and 0.5e13 for the second, [R5252]); the model's
own range for how strongly burial depends on bottom-water oxygen is 10 to 40 percent. Tier 3, tweak level 1."""

REBALANCE_BRACKETS = 40
"""int: The most times `rebalance_co2` doubles the CO2 pressure while it looks for one the sink holds less than (GEN.1589: the bracket of the root; it stops at the retention cap)."""

REBALANCE_BISECTIONS = 60
"""int: The most halvings of that bracket (GEN.1589). Sixty halvings narrow any bracket of an air to a part in 10^15, so a run that uses them all found no root: the sink and the temperature did not
settle (a switch in the air that makes the answer jump), and it is logged."""

REBALANCE_TOLERANCE = 0.01
"""float: The relative width of the bracket (and the relative change of the CO2) under which `rebalance_co2` stops."""

HAZE_PASSES = 12
"""int: The most times `store_air` works the haze of the air and the temperature against each other (GEN.1589), as the air solver does at its first answer; the last answer stands when they do not settle."""

HAZE_TOLERANCE_K = 0.01
"""float: The temperature change under which `store_air` stops working the haze and the temperature against each other. A stand-in numerical tolerance."""

APPLY_RISE = 1.01
"""float: The air takes new oxygen only when the target is this much above what it holds, so that applying the model again changes nothing."""


REDUCED_GASES = {"H2": (0.5, {}), "CO": (0.5, {"CO2": 1.0}), "CH4": (2.0, {"CO2": 1.0}), "NH3": (0.75, {"N2": 0.5}),
                 "H2S": (2.0, {}), "C2H6": (3.5, {"CO2": 2.0})}
"""dict: The gases oxygen burns: formula to (the moles of O2 a mole takes, the other gases it leaves per mole, mol). The water they make is
not in the products: it goes to the ocean, and the sulfur of H2S to the sulfate of the sea (stand-in T5, L0)."""


@dataclass(frozen=True)
class Profile:
    """
    What oxygenic photosynthesis would do on a world.

    Attributes:
        source_mol_yr (float): `S`, the oxygen its life makes a year.
        reductant_mol_yr (float): `R`.
        store_mol (float): `B`, with the reduced gases of the air.
        steady_mol (float): The oxygen, mol, the air holds once the store is spent, within the fire limit; 0 when `S <= R`.
        steady_pa (float): That oxygen's partial pressure, Pa, in the air as it stands once the reduced gases are burned.
        delay_gyr (float or None): How long after the first photosynthesis the air holds that, Gyr; `None` when it never does.
    """

    source_mol_yr: float
    reductant_mol_yr: float
    store_mol: float
    steady_mol: float
    steady_pa: float
    delay_gyr: float | None


NONE = Profile(0.0, 0.0, 0.0, 0.0, 0.0, None)


def _clamp(value, low, high):
    return min(max(value, low), high)


def redox_factor(delta_iw):
    """How many times more reductants than Earth's a mantle of `delta_iw` puts out (and holds in store); 1 for none."""
    if delta_iw is None:
        return 1.0
    low, high = REDOX_FACTOR_RANGE
    return _clamp(10.0 ** (REDOX_DECADE_PER_UNIT * (EARTH_MANTLE_DELTA_IW - delta_iw)), low, high)


def temperature_factor(temperature_k):
    """The share of Earth's production a sea at `temperature_k` allows: a Boltzmann factor inside the liquid range, 0 outside it."""
    if not temperature_k:
        return 0.0
    cold, warm_from, warm_to, hot = OCEAN_SPAN_K
    if temperature_k <= cold or temperature_k >= hot:
        return 0.0
    edge = min((temperature_k - cold) / (warm_from - cold), (hot - temperature_k) / (hot - warm_to), 1.0)
    boltzmann = math.exp(TEMPERATURE_ACTIVATION_EV * constants.ELECTRON_VOLT_J / constants.BOLTZMANN
                         * (1.0 / REFERENCE_TEMPERATURE_K - 1.0 / temperature_k))
    return edge * boltzmann


def active_light_ratio(body, band=PAR_BAND):
    """The light of `band` the body gets, over Earth's photosynthetically active light: the star's flux at the body's distance times the share of a blackbody
    of the star's temperature in the band's `nm`, times its `efficiency` (GEN.1744: class 1b eats 300 to 500 nm, 1c 700 to 800 nm at 0.8, 1a `PAR_NM`). 0 when the star states no light."""
    star = getattr(body, "star", None)
    luminosity = getattr(star, "luminosity", None)
    temperature = getattr(star, "temperature", None)
    distance = getattr(body, "_star_distance_au", None) or getattr(body, "distance", None)
    if not luminosity or not temperature or not distance:
        return 0.0
    here = stellar_spectrum.bolometric_flux_w_m2(luminosity, distance) * float(stellar_spectrum.blackbody_fraction(*band.nm, temperature)) * band.efficiency
    earth = stellar_spectrum.bolometric_flux_w_m2(constants.SOLAR_LUMINOSITY, 1.0) * float(
        stellar_spectrum.blackbody_fraction(*PAR_NM, SUN_TEMPERATURE_K))
    return here / earth


def reductant_flux(body, area_m2):
    """
    `R`, mol/yr: Earth's, times the heat the body's interior and tide give off over Earth's, times `redox_factor`. A body a Class M abundance tier drew with a quiet interior
    (`abundance_tier.quiet_interior`, GEN.1721; never on `science`) gives `tuning.ABUNDANCE_QUIET_SINK_FACTOR` of that.
    """
    flux = (getattr(body, "internal_heat_flux_w_m2", None) or 0.0) + (getattr(body, "tidal_flux_w_m2", None) or 0.0)
    quiet = tuning.ABUNDANCE_QUIET_SINK_FACTOR if abundance_tier.quiet_interior(body) else 1.0
    return EARTH_REDUCTANT_MOL_YR * (flux * area_m2 / EARTH_HEAT_W) * redox_factor(getattr(body, "mantle_delta_iw", None)) * quiet


def reductant_history_factor(body, age_gyr):
    """
    `R(age) / R(now)` of a rocky body (GEN.1463): the heat its interior gave off at `age_gyr` over the heat it gives off at its star's age (`interior_heat.heat_at`, the body's own
    thermal history, the tide stored with it counted as it is), so the reductant that the volcanoes and the crust gave the oxygen to was larger when the body was young, in the
    proportion `reductant_flux` already takes from the heat. 1 for a body with no thermal history stored, none at all, or no heat now.
    """
    if body is None:
        return 1.0
    from planetgen.physics import interior_heat
    now_age = getattr(getattr(body, "star", None), "age", None)
    if not now_age:
        return 1.0
    then, now = interior_heat.heat_at(body, min(max(age_gyr, 0.0), now_age)), interior_heat.heat_at(body, now_age)       # no history beyond the star's age: R is what it is now
    if then is None or now is None:
        return 1.0
    tide = getattr(body, "tidal_flux_w_m2", None) or 0.0
    flux_now = now.own_flux_w_m2 + tide
    return 1.0 if flux_now <= 0.0 else (then.own_flux_w_m2 + tide) / flux_now


def refuge_productivity(body, area_m2, band=PAR_BAND):
    """
    `S`, mol/yr, of the oxygenic life of a body that has no wide sunlit sea but has light and liquid water at its surface, in a brine seep, a seasonal melt or a spring
    (`surface_refuge`, GEN.1723): Earth's burial per area of ocean (no new constant) on the area where those are wet and lit, each place weighed by its production against an ocean's
    (`surface_refuge.Refuge.area`), under the same light and temperature factors as a sea. 0 for a body with none. About 1e-2 of Earth's source for a seep world and 0.14 at the very most (the cover cap over a whole sea that melts all year), so it never
    beats the reductant flux of a world whose interior is Earth's: if a dry world is to hold oxygen the lever is `R`, not this term.

    Once the life has reached the rooted stage (`surface_refuge.ROOTED`) the groundwater-fed places carry land plants (GEN.1729): their wet land (`Refuge.plants`) makes Earth's burial per area of land under plants, `LAND_SHARE_OF_BURIAL`
    of the total over Earth's 0.29 of land, in place of the cover law on those places (`Refuge.plants_cover`). A mostly-land desert world with roots over a quarter of its land makes about 0.2 of Earth's source then, which beats the
    reductants of a quiet interior (under about 0.7 of Earth's) and not an Earth-like one.
    """
    if not getattr(body, "surface_temperature", None) or not area_m2:
        return 0.0
    if getattr(body, "hydrosphere", None) not in ("brine seep", "ice-covered ocean", "ice") and not getattr(body, "aquifer_thickness_km", None):
        return 0.0          # the common case, none of the places a refuge can be in (`surface_refuge.refuge_of`), without reading the body's columns
    refuge = surface_refuge.refuge_of(habitability_world.values_from_body(body))
    if refuge.area <= 0.0 or not refuge.temperature_k:
        return 0.0
    light = _clamp(active_light_ratio(body, band), *LIGHT_RANGE) ** LIGHT_EXPONENT
    sea_like = (refuge.area - refuge.plants_cover) / EARTH_OCEAN_FRACTION
    land_plants = LAND_SHARE_OF_BURIAL * refuge.plants / EARTH_LAND_FRACTION      # GEN.1729: the wet land that roots reach, at Earth's burial per area of land under plants
    return EARTH_BURIAL_MOL_YR * (area_m2 / EARTH_AREA_M2) * (sea_like + land_plants) * light * temperature_factor(refuge.temperature_k)


def productivity(body, area_m2, band=PAR_BAND):
    """`S`, mol/yr, of a body with a sunlit sea (a body with a refuge only: `refuge_productivity`; 0 for one with neither) of the oxygenic class that eats `band` (GEN.1744; class 1a's by default)."""
    if getattr(body, "hydrosphere", None) not in LIVING_HYDROSPHERES:
        return refuge_productivity(body, area_m2, band)
    ocean = getattr(body, "ocean_fraction", None) or 0.0
    light = _clamp(active_light_ratio(body, band), *LIGHT_RANGE) ** LIGHT_EXPONENT
    return (EARTH_BURIAL_MOL_YR * (area_m2 / EARTH_AREA_M2) * min(ocean / EARTH_OCEAN_FRACTION, 1.3) * light
            * temperature_factor(getattr(body, "surface_temperature", None)))


def area_m2(body):
    return 4.0 * math.pi * (body.radius * 1000.0) ** 2


def air_moles(body):
    """The moles of each gas in the air of `body` (formula to mol): the column `P A / (g M)` of the mixture, shared by the mole fractions."""
    partials = atmosphere.gas_partials_pa(body)
    if not partials:
        return {}
    total = sum(partials.values())
    column = total * area_m2(body) / (body.gravity * constants.EARTH_GRAVITY * mixture.mean_molar_mass_kg_mol(partials))
    return {formula: pascal / total * column for formula, pascal in partials.items()}


def burn(moles):
    """
    The air `moles` (formula to mol) once its reduced gases have taken oxygen: `(new moles, mol of O2 they took)`. The water made is
    not kept in the air (it goes to the ocean).
    """
    result = {formula: amount for formula, amount in moles.items() if formula not in REDUCED_GASES}
    used = 0.0
    for formula, (oxygen, products) in REDUCED_GASES.items():
        amount = moles.get(formula, 0.0)
        used += amount * oxygen
        for product, per_mole in products.items():
            result[product] = result.get(product, 0.0) + amount * per_mole
    return result, used


def _fire_limit_mol(body, moles):
    """
    The most oxygen, mol, that the other gases of `moles` leave room for on `body` under the fire limit (GEN.1320): the oxygen's partial pressure may reach
    `mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA` and its share of the molecules `mixture_oxygen.ABIOTIC_O2_MAX`, whichever is lower.

    With `n` mol of oxygen among `rest` mol of other gases of weight `weight` kg, the partial pressure is `n / (n + rest)` of `(weight + n m) g / A`, which the cap
    `c` equals where `m n^2 + (weight - c A / g) n - c A rest / g = 0`.
    """
    rest = sum(amount for formula, amount in moles.items() if formula != "O2")
    if rest <= 0.0:
        return 0.0
    by_share = mixture_oxygen.ABIOTIC_O2_MAX / (1.0 - mixture_oxygen.ABIOTIC_O2_MAX) * rest
    return min(by_share, _oxygen_mol_at(body, moles, mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA))


def _oxygen_mol_at(body, moles, pascal):
    """
    The oxygen, mol, that makes its partial pressure `pascal` on `body` among the other gases of `moles` (formula to mol), 0 when they hold nothing.

    With `n` mol of oxygen among `rest` mol of other gases of weight `weight` kg, the partial pressure is `n / (n + rest)` of `(weight + n m) g / A`, which `pascal`
    equals where `m n^2 + (weight - c A / g) n - c A rest / g = 0` (`c` the pressure).
    """
    rest = sum(amount for formula, amount in moles.items() if formula != "O2")
    if rest <= 0.0:
        return 0.0
    weight = sum(amount * gases.gas(formula).molar_mass for formula, amount in moles.items() if formula != "O2")
    molar = gases.gas("O2").molar_mass
    cap = pascal * area_m2(body) / (body.gravity * constants.EARTH_GRAVITY)       # kg of column weight the pressure equals
    b = weight - cap
    return (-b + math.sqrt(b * b + 4.0 * molar * cap * rest)) / (2.0 * molar)


def stage_target_mol(body, profile):
    """
    The oxygen, mol, that the air of `body` holds at the star's age (GEN.1707), given its `profile` (the balance at the summed source) at the reductant sink `R` the body has now: the balance level `steady_mol` once `R` is
    under `S / SECOND_RISE_RATIO` (the second rise), and until then the plateau at `PLATEAU_PA` (never more than the balance). The sink falls with the interior's heat, so a body that is old enough has crossed both.
    """
    if profile.steady_mol <= 0.0:
        return 0.0
    if profile.reductant_mol_yr <= profile.source_mol_yr / SECOND_RISE_RATIO:
        return profile.steady_mol
    burnt, _ = burn(air_moles(body))
    return min(profile.steady_mol, _oxygen_mol_at(body, burnt, PLATEAU_PA))


def partials_pa(body, moles):
    """The partial pressures, Pa, of `moles` (formula to mol) on `body`: the weight of the column over the area, shared by the fractions."""
    column = sum(moles.values())
    weight = sum(amount * gases.gas(formula).molar_mass for formula, amount in moles.items())
    total = weight * body.gravity * constants.EARTH_GRAVITY / area_m2(body)
    return {formula: amount / column * total for formula, amount in moles.items() if amount > 0.0}


def band_of(life_class):
    """The `Band` of the light the oxygenic class `life_class` eats: its `energy` row's `band_nm` and `efficiency` in `life_classes.json` (GEN.1744)."""
    energy = life_class.energy
    return Band(tuple(energy["band_nm"]), energy["efficiency"])


def source_of(body, area, producers):
    """`S`, mol/yr, of the oxygenic classes `producers`, a list of `(share, band)`: each class carries its share of the body's food production (GEN.1213) in the light of its own band (GEN.1744)."""
    return sum(share * productivity(body, area, band) for share, band in producers)


def balance_pa(surplus_mol_yr, area):
    """
    The oxygen partial pressure, Pa, at which oxidative weathering takes the `surplus_mol_yr` (`S - R`) that a body of `area` m2 makes: Earth's `EARTH_O2_PA` times the surplus
    per area over Earth's, to the power `1 / OXIDATIVE_WEATHERING_EXPONENT` (GEN.1733). 0 when there is no surplus.
    """
    if surplus_mol_yr <= 0.0 or area <= 0.0:
        return 0.0
    per_area = (surplus_mol_yr / area) / (EARTH_SURPLUS_MOL_YR / EARTH_AREA_M2)
    return EARTH_O2_PA * per_area ** (1.0 / OXIDATIVE_WEATHERING_EXPONENT)


def profile_of(body, share=1.0, extra_store_mol=0.0, producers=None):
    """
    The `Profile` of a rocky body (`NONE` for a giant, or a body missing the radius, gravity or air to read). `share` is the part of the body's food
    production that oxygenic life carries (GEN.1213: the classes of a body share its production, so the source `S` goes with their summed share), counted in class 1a's light (`PAR_BAND`);
    `producers`, a list of `(share, band)` of the oxygenic classes (GEN.1744: 1b eats blue, 1c far red), is that sum in each class's own band and replaces `share`;
    `extra_store_mol` is the oxygen, mol, that other reduced things the life made (the methane of the classes before the rise) take first, which adds to the store.
    """
    radius_km, gravity = getattr(body, "radius", None), getattr(body, "gravity", None)
    if getattr(body, "body_type", "t") == "g" or not radius_km or not gravity:
        return NONE
    area = area_m2(body)
    moles = air_moles(body)
    burnt, taken = burn(moles)
    source, reductant = source_of(body, area, producers if producers is not None else [(share, PAR_BAND)]), reductant_flux(body, area)
    store = EARTH_STORE_MOL * (area / EARTH_AREA_M2) * redox_factor(getattr(body, "mantle_delta_iw", None)) + taken + extra_store_mol
    surplus = source - reductant
    if surplus <= 0.0 or not burnt:
        return Profile(source, reductant, store, 0.0, 0.0, None)
    steady = min(_oxygen_mol_at(body, burnt, balance_pa(surplus, area)), _fire_limit_mol(body, burnt))
    pascal = partials_pa(body, {**burnt, "O2": steady}).get("O2", 0.0)
    if pascal < MINIMUM_PA:
        return Profile(source, reductant, store, 0.0, 0.0, None)
    return Profile(source, reductant, store, steady, pascal, store / surplus / 1e9)


def rise_age_gyr(body, starts, extra_store_mol=0.0):
    """
    The age, Gyr, at which the oxygen of several oxygenic classes first rises on `body` (GEN.1213), or `None` when it never does. `starts` is a list of
    `(age_gyr, share)` or `(age_gyr, share, band)`: each class's producers begin at its age and carry its share of the food production, in the light of its own `band` (GEN.1744; class 1a's
    when not given). The surplus `S - R(age)` of the producers
    there are fills the store `B` (the first oxygen goes to the reduced things there), `S` going with the shares of the classes begun, each in its band (`carried`, mol a year), and `R` with the
    heat the body gave off at that age (GEN.1463, `reductant_history_factor`); the oxygen rises when it is full. The surplus is summed in steps of `RISE_STEP_GYR` (the
    last step cut where the store fills) from the first producers' age; past the star's age the reductant is held at its value now.
    """
    starts = sorted((age, share, (rest[0] if rest else PAR_BAND)) for age, share, *rest in starts if age is not None and share > 0.0)
    if not starts:
        return None
    full = profile_of(body, 1.0, extra_store_mol)
    unit = {band: profile_of(body, 1.0, extra_store_mol, producers=[(1.0, band)]).source_mol_yr for band in {band for _, _, band in starts}}       # S of a class that carries it all, in its own band
    if max(unit.values()) <= 0.0:
        return None
    horizon = getattr(getattr(body, "star", None), "age", None) or 0.0       # the thermal history ends at the star's age; past it R is what it is now
    filled, age, carried, index = 0.0, starts[0][0], 0.0, 0
    while True:
        while index < len(starts) and starts[index][0] <= age + 1e-12:
            carried += starts[index][1] * unit[starts[index][2]]
            index += 1
        if age >= horizon and index >= len(starts):
            surplus = carried - full.reductant_mol_yr
            return None if surplus <= 0.0 else age + (full.store_mol - filled) / surplus / 1e9      # the surplus at the horizon runs on
        end = min(age + RISE_STEP_GYR, starts[index][0]) if index < len(starts) else age + RISE_STEP_GYR
        surplus = carried - full.reductant_mol_yr * reductant_history_factor(body, 0.5 * (age + end))
        if surplus > 0.0:
            need = full.store_mol - filled
            if surplus * (end - age) * 1e9 >= need:
                return age + need / surplus / 1e9
            filled += surplus * (end - age) * 1e9
        age = end


def apply(body, target_mol):
    """
    Puts `target_mol` of oxygen into the air of `body`, burning its reduced gases first (`burn`): the stored partial pressures follow from the
    new moles, with the total pressure, molar mass, density and scale height of the new mixture, and the warming and surface temperature
    that mixture gives under the starlight the body absorbed (the closure of `surface_air`, so the stored greenhouse factor still follows
    from the stored air). The oxygen is held within the fire limit and an air that has more already keeps it. Nothing happens, and `False`
    returns, when the air has no reduced gas to burn and holds that oxygen already (within `APPLY_RISE`), or has no air to take it. The doses,
    ozone, colour and score that read the air are the caller's to make again (`generation.life`).
    """
    moles = air_moles(body)
    if not moles or target_mol <= 0.0:
        return False
    burnt, taken = burn(moles)
    have = burnt.get("O2", 0.0)
    wanted = max(have, min(target_mol, _fire_limit_mol(body, burnt)))
    if taken <= 0.0 and wanted <= have * APPLY_RISE:
        return False
    burnt["O2"] = wanted
    air = partials_pa(body, burnt)
    cap = getattr(body, "_retention_cap_pa", None)
    total = sum(air.values())
    if cap is not None and total > cap:        # GEN.1340: the air is held to the pressure its retention regime keeps (as the pressure stage holds it), oxygen and all
        air = {formula: pascal * cap / total for formula, pascal in air.items()}
    if "O2" not in air or air["O2"] < MINIMUM_PA:
        return False
    store_air(body, air)
    return True


def rebalance_co2(body):
    """
    Solves the carbon sink of the air again at the temperature the oxygen rise left (GEN.1309 A). The air was solved with the methane and hydrogen that warmed it; `apply` burned
    them (the methane's carbon is CO2 now, one for one) and the air found another surface temperature (7 to 67 K lower in the worlds this was measured on), but the sink that holds
    the CO2 of a world (continental weathering, seafloor weathering, or the sealed ocean's chemistry: `surface_air.CarbonSink`) was not asked again. It is asked at the new
    temperature: the pressure it holds there against the one it held where the air was solved shifts the CO2 by the same amount, and the greenhouse of the new CO2 gives the temperature
    again. The CO2 the sink holds at the temperature a CO2 gives is lower the more CO2 there is (the greenhouse warms the air and the sink holds less), so the CO2 that equals it is the one
    root of a rising difference, found by bracketing and halving (GEN.1589; the damped passes it replaced stopped after 8 with the answer still falling, 562 K to 432 K in one air, and were
    left where they stood). The root is held to `REBALANCE_TOLERANCE`, and a run that does not settle in `REBALANCE_BISECTIONS` is logged. A body with no sink (no carbon to weather, or an air the oddity
    rows replaced) is left as it is.

    Returns:
        bool: Whether the air changed.
    """
    sink = getattr(body, "_co2_sink", None)
    if sink is None or not sink.solved_bar > 0.0:
        return False
    air = dict(atmosphere.gas_partials_pa(body))
    base = air.get("CO2", 0.0)
    current = base
    cap = getattr(body, "_retention_cap_pa", None)
    room = float("inf") if cap is None else max(cap - sum(pascal for formula, pascal in air.items() if formula != "CO2"), 0.0)   # the retention regime's cap holds the whole air: the carbon dioxide takes what room the other gases (the oxygen among them, set by their own balance) leave

    started = (body.surface_temperature, getattr(body, "greenhouse_factor", None))

    def store(carbon_dioxide):
        """The air with `carbon_dioxide` stored, from the temperature and greenhouse factor the body began with, so every trial works from the same flux and haze start."""
        body.surface_temperature, body.greenhouse_factor = started
        store_air(body, {formula: pascal for formula, pascal in {**air, "CO2": carbon_dioxide}.items() if pascal > 0.0})

    def difference(carbon_dioxide):
        """The CO2 asked of the air less the CO2 the sink holds at the temperature that CO2 gives: below zero while the sink wants more."""
        store(carbon_dioxide)
        return carbon_dioxide - max(base + (sink.kept_at(body.surface_temperature) - sink.solved_bar) * 1.0e5, 0.0)

    low, high = 0.0, min(max(current, MINIMUM_PA), room)
    for _ in range(REBALANCE_BRACKETS):
        if difference(high) >= 0.0 or high >= room:
            break
        low, high = high, min(2.0 * high, room)
    settled = False
    for _ in range(REBALANCE_BISECTIONS):
        middle = 0.5 * (low + high)
        if high - low <= REBALANCE_TOLERANCE * max(middle, MINIMUM_PA):
            settled = True
            break
        if difference(middle) < 0.0:
            low = middle
        else:
            high = middle
    if not settled:
        log.warning(f"The carbon dioxide of a world did not settle against its sink: between {low:.4g} and {high:.4g} Pa after {REBALANCE_BISECTIONS} halvings.")
    answer = min(high, room)
    if abs(answer - current) <= REBALANCE_TOLERANCE * max(current, MINIMUM_PA):
        store(current)          # the trials left the body at the last of them: the air it was given is stored again, solved as the trials were
        return False
    store(answer)
    return True


FREEZE_AGAIN_PASSES = 20
"""int: The most times `store_air` lets the warming and the freeze-out of the new air move each other (the colder air condenses more, the condensed air warms less, and the pair settles in a few passes; a numerical bound, GEN.1564)."""
FREEZE_AGAIN_TOLERANCE = 1.0e-6
"""float: The relative change of a partial pressure under which the freeze-out has settled. A numerical tolerance, not a physical value (GEN.1564)."""


def freeze_out_at(body, air, temperature_k):
    """
    `air` (formula to partial pressure, Pa) once the gases that stand above their saturation pressure at `temperature_k` have condensed (GEN.1564): the same
    `surface_air.freeze_out` the first solve uses, on the moles of the column. A gas with no curve, or above its critical temperature, is left alone. The water
    of the body's ocean (`_ocean_moles`) is the reservoir of the vapor, as in the first solve: it brings the vapor to its humidity from either side, so a sea
    that is not frozen over keeps its air wet at a colder temperature; no other gas has a reservoir here.
    """
    gravity = body.gravity * constants.EARTH_GRAVITY
    area = area_m2(body)
    weight = sum(pascal * gases.gas(formula).molar_mass for formula, pascal in air.items())
    total_pa = sum(air.values())
    if weight <= 0.0 or total_pa <= 0.0:
        return air
    mean_molar_mass = weight / total_pa
    moles = {formula: pascal * area / (gravity * mean_molar_mass) for formula, pascal in air.items()}
    ocean_moles = _ocean_moles(body, area)
    reservoirs = {"H2O": (surface_air.RELATIVE_HUMIDITY, ocean_moles)} if ocean_moles > 0.0 else None
    left, _, _ = surface_air.freeze_out(moles, gravity, area, temperature_k, reservoirs)
    left = {formula: count for formula, count in left.items() if count > 0.0}
    column = sum(left.values())
    if column <= 0.0:
        return {}
    pressure = sum(count * gases.gas(formula).molar_mass for formula, count in left.items()) * gravity / area
    return {formula: count / column * pressure for formula, count in left.items()}


def _ocean_moles(body, area):
    """Moles of water in the body's ocean (its cover times its mean depth, at liquid water's density), 0 for a body with no sea."""
    depth_km = getattr(body, "ocean_depth_km", None) or 0.0
    cover = getattr(body, "ocean_fraction", None) or 0.0
    if getattr(body, "hydrosphere", None) not in ("surface ocean", "hycean", "ice-covered ocean") or depth_km <= 0.0 or cover <= 0.0:
        return 0.0
    return cover * area * depth_km * 1000.0 * constants.LIQUID_WATER_DENSITY_KG_M3 / gases.gas("H2O").molar_mass


def store_air(body, air):
    """
    Gives `body` the air `air` (formula to partial pressure, Pa): the warming and surface temperature that mixture gives under the starlight the
    body absorbed (the closure of `surface_air`, so the stored greenhouse factor still follows from the stored air, haze included: GEN.1589), and the total pressure,
    molar mass, density and scale height of the new mixture. The doses, aerosols, colour and score that read the air are the caller's to make
    again (`generation.life`).

    The new air has a new temperature, and at it a gas may stand above its saturation pressure: the freeze-out of the first solve is run again at that
    temperature (`freeze_out_at`, GEN.1564), and the warming of the air that is left, until the pair settles. A hydrogen-rich air keeps its gases (the
    first solve does not freeze a giant either).
    """
    gravity = body.gravity * constants.EARTH_GRAVITY
    radiated = constants.STEFAN_BOLTZMANN_CONSTANT * body.surface_temperature ** 4 / (1.0 + body.greenhouse_factor)
    warming = _warming(body, air, gravity, radiated)
    for _ in range(FREEZE_AGAIN_PASSES):
        frozen = freeze_out_at(body, air, warming.surface_temperature_k)
        if not frozen:
            break
        if all(abs(frozen.get(formula, 0.0) - pascal) <= FREEZE_AGAIN_TOLERANCE * max(pascal, MINIMUM_PA) for formula, pascal in air.items()):
            break
        air = frozen
        warming = _warming(body, air, gravity, radiated)
    state = atmosphere_state.air_state(air, warming.surface_temperature_k, gravity)
    body.surface_temperature = warming.surface_temperature_k
    body.greenhouse_factor = warming.factor
    body.atmospheric_pressure = state.total_pressure_pa
    body.atm_molar_density = state.mean_molar_mass_kg_mol
    body.atm_density = state.density_kg_m3
    body.scale_height = state.scale_height_m / 1000.0
    atmosphere.store_air(body, air, body.mantle_delta_iw)
    envelope_share.fill(body)       # GEN.1226: the envelope share follows the air the life has changed


def _warming(body, air, gravity, radiated):
    """
    The warming of `air` under the flux `radiated` (W/m2, the flux the stored temperature and greenhouse factor were made under), with the haze of the air as the air solver takes it
    (GEN.1589): the haze takes a share of the absorbed sunlight high up (`aerosols.haze_of`, the body's own draws, kept from its air step as `_air_haze`), the share and the pressure it
    takes it at go to the greenhouse table, and the haze depends on the temperature, so the two are worked against each other (damped, `HAZE_PASSES`). A body that kept no haze
    (a world read back from the database, or an air step with no haze draws) gets the clear-air warming, as before.
    """
    kept = getattr(body, "_air_haze", None)
    if kept is None:
        return atmosphere_state.greenhouse_warming(air, gravity, radiated)
    haze_draws, star_temperature_k, far_uv, internal_flux_w_m2 = kept
    temperature = body.surface_temperature
    warming = None
    for _ in range(HAZE_PASSES):
        haze = aerosols.haze_of(air, gravity, temperature, haze_draws, star_temperature_k, far_uv)
        aloft_share = haze.absorbed_share * max(radiated - internal_flux_w_m2, 0.0) / radiated     # the haze takes a share of the sunlight, not of the heat the body radiates from inside
        warming = atmosphere_state.greenhouse_warming(air, gravity, radiated, aloft_share, haze.pressure_pa)
        found = warming.surface_temperature_k
        if abs(found - temperature) <= HAZE_TOLERANCE_K:
            break
        temperature = 0.5 * (temperature + found)
    return warming
