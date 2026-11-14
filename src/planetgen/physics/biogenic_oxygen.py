# planetgen/physics/biogenic_oxygen.py

"""
Free oxygen from oxygenic photosynthesis (GEN.855, docs/design/habitability/life-timeline-and-oxygen.md)
=========================================================================================================

A world whose life splits water for its electrons (a life class of `physics/life_classes.py` whose `metabolism` is `"oxygenic"`) fills its air with oxygen; one whose life does
not, does not. The oxygen a world holds is the balance of a source and its sinks, in moles of O2 a year:

* **Source** `S`: the organic carbon its life buries, which leaves the oxygen of the water it was made from behind. It goes with the
  area where light, water and warmth meet (`productivity`): the planet's area, the share the ocean covers, the square root of the
  photosynthetically active light it gets (the production saturates), and a Boltzmann factor for the temperature of the surface, to zero
  outside the liquid range of the sea.
* **Reductant flux** `R`: the volcanic and metamorphic gases (H2, CO, H2S) that take oxygen as they are released, which goes with the
  heat the interior gives off and with how reduced the mantle is (`reductant_flux`).
* **Oxidative weathering** `k N`: the rocks at the surface (kerogen, pyrite, ferrous iron) take oxygen in proportion to what the air holds,
  `N` moles, with the rate constant that balances the budget on the world the numbers are taken from (`WEATHERING_PER_YR`).
* **A store of reductants** `B`: the ocean's ferrous iron and the crust's sulfides and iron, and the reduced gases of the air itself (H2, CO,
  CH4, NH3, H2S, C2H6: `REDUCED_GASES`), take the first oxygen the surplus makes, so the air holds none until `B / (S - R)` has passed
  (`Profile.delay_gyr`, the time the Great Oxidation took after the first oxygenic photosynthesis). The reduced gases are burned away by
  then: `apply` turns them into their oxides (CO2, N2) and the water they make goes to the ocean.

`dN/dt = S - R - k N` once the store is spent, so the air settles to `N* = (S - R) / k` in a few million years, a time short beside the
Gyr the schedule runs in, and the world holds `N*` (as the partial pressure `p = N M g / A`) when its photosynthesis is older than the
delay and nothing before. A world with `S <= R` never holds any. The most the air can hold is the oxygen pressure `mixture_oxygen.BIOTIC_O2_PRESSURE_MAX_PA`, the
fire limit (GEN.1320: a pressure, so a thin air can reach a larger share than a thick one; never more than `mixture_oxygen.ABIOTIC_O2_MAX` of the molecules).
The whole air, oxygen and all, stays under the pressure the retention regime of the air solve lets the body keep (`surface_air.Air.retention_cap_pa`, the row `d6`; GEN.1340): a body that
holds only 1 or 10 kPa against the escape of its air does not hold 40 kPa because its life made oxygen, and the shares of its gases are kept when the cap takes the excess.

The numbers are those of the one world we know (`EARTH_*`, round figures from the literature of the oxygen cycle: [R0689, R1980], from memory), used as the test point and the scale and not as a preset: a world
with another light, water, heat, redox and gravity gets another oxygen. Every constant is a stand-in in section 6.21 of
`docs/design/research/gap-filling-plan.md` with its tier and tweak level.

What it leaves out: the loss of oxygen to space (a few parts in 1e4 of the inventory a Gyr for a habitable-zone world, far under the
weathering), the hydrogen escape that adds oxygen, the phosphorus feedback and the land plants that raised Earth's burial (they would
change `S` by a factor of a few, inside the 1.7 to 4 that this module's L3 tweak already stands for).
"""

import math
from dataclasses import dataclass

from planetgen.physics import aerosols, atmosphere, atmosphere_state, constants, envelope_share, gases, mixture, mixture_oxygen, stellar_spectrum
from planetgen.util import log

EARTH_O2_MOL = 3.7e19
"""float: The oxygen of Earth's air, mol (21 percent of 1.77e20 mol of air; 22.8 kPa at its gravity and area)."""

EARTH_BURIAL_MOL_YR = 1.0e13
"""float: The organic carbon Earth buries a year, mol (the source `S`; [R1859], about 1e13). T2, L3."""

EARTH_REDUCTANT_MOL_YR = 3.0e12
"""float: The oxygen Earth's volcanic and metamorphic reductants take a year, mol (`R`; [R1980] few 1e12). T3, L3."""

EARTH_STORE_MOL = 4.0e21
"""float: The oxygen the ocean's iron and the crust's sulfides and iron take before the air holds any, mol (`B`), chosen so that the first
oxygen follows the first photosynthesis by 0.7 Gyr when `S` is Earth's and `R` is Earth's at that age, from the model's own heat history of an
Earth-mass body (GEN.1463: 2.2 times today's heat when the biomarkers appear, 2.7 Ga, falling to 1.9 times at the rise; it was 5e21 while `R` was held
at today's value), the Great Oxidation coming some 0.3 to 1 Gyr after oxygenic photosynthesis began. T3, L1."""

RISE_STEP_GYR = 0.25
"""float: The step, Gyr, of the integral of the surplus `S - R(age)` in `rise_age_gyr`; the heat flow it reads changes by under 10 percent over a step after the first Gyr."""

EARTH_AREA_M2 = 4.0 * math.pi * constants.EARTH_RADIUS_M ** 2
EARTH_OCEAN_FRACTION = 0.71
EARTH_HEAT_W = 4.4e13
"""float: Earth's own heat flow now, W (44 TW)."""

EARTH_MANTLE_DELTA_IW = 3.5
"""float: The redox of Earth's upper mantle now, log10 units above the iron-wustite buffer (the scale `oxygen_fugacity` gives `mantle_delta_iw`)."""

WEATHERING_PER_YR = (EARTH_BURIAL_MOL_YR - EARTH_REDUCTANT_MOL_YR) / EARTH_O2_MOL
"""float: The share of the air's oxygen that weathering takes a year, 1/yr, which balances Earth's budget at `EARTH_O2_MOL` (1.9e-7)."""

REDOX_DECADE_PER_UNIT = 0.25
REDOX_FACTOR_RANGE = (0.2, 5.0)
"""A mantle one log unit more reduced than Earth's puts out `10 ** 0.25` times the reductants and holds that times the store, kept inside the range."""

LIGHT_EXPONENT = 0.5
"""float: Production goes as the active light to this power: it saturates (T6, L3)."""

LIGHT_RANGE = (0.0, 2.0)
PAR_NM = (400.0, 700.0)
"""tuple[float]: The photosynthetically active band, nm."""

SUN_TEMPERATURE_K = constants.SOLAR_TEMPERATURE_K

TEMPERATURE_ACTIVATION_EV = 0.32
"""float: The activation energy of photosynthesis, eV (about a half of a heterotroph's; T2, L0)."""

REFERENCE_TEMPERATURE_K = 288.0
OCEAN_SPAN_K = (258.0, 275.0, 315.0, 343.0)
"""tuple[float]: The temperatures, K, at which the factor for the sea's liquid range rises from 0 to 1 and falls to 0 again (T6, L3)."""

LIVING_HYDROSPHERES = ("surface ocean", "hycean")
"""tuple: The states with a sunlit sea for the photosynthesis to stand in (a `brine seep` is not one; GEN.908)."""

MINIMUM_PA = 1.0
"""float: Under this partial pressure, Pa, a world is counted as making no oxygen."""

EUKARYOTE_MINIMUM_PA = 200.0
"""float: The oxygen, Pa, the air must hold for the wait for complex cells to start (about 1 percent of Earth's: [R3084]; T3, L0)."""

PLATEAU_PA = 600.0
"""float: The oxygen partial pressure of the plateau, Pa: where the oxygen stands after the first rise (the reductant sink `R` under the source `S`) until the second (GEN.1707). The models of Earth's
plateau put it at 0.06 to 0.25 of today's oxygen, where the proxies say at most 0.001 (`docs/design/research/` oxygen plateau report: 0.001 to 0.25, 0.02 to 5 kPa); 0.03 of Earth's 21 kPa is 0.6 kPa. Tier 6 (educated guess from physics), tweak level 1."""

SECOND_RISE_RATIO = 3.0
"""float: The oxygen rises again to the balance level once the reductant sink `R` is under this share of the source `S`, 1 over the ratio (GEN.1707): the deep ocean and the shelf stop recycling phosphorus
and the oxygen gets a second jump. A single model's number (the reductant input crossing 1.5e13 mol O2-equivalent a year for the first rise and 0.5e13 for the second, Alcott, Mills and Poulton 2019 [R5252]); the model's
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


def active_light_ratio(body):
    """The photosynthetically active light the body gets, over Earth's: the star's flux at the body's distance times the share of a blackbody
    of the star's temperature in `PAR_NM`. 0 when the star states no light."""
    star = getattr(body, "star", None)
    luminosity = getattr(star, "luminosity", None)
    temperature = getattr(star, "temperature", None)
    distance = getattr(body, "_star_distance_au", None) or getattr(body, "distance", None)
    if not luminosity or not temperature or not distance:
        return 0.0
    here = stellar_spectrum.bolometric_flux_w_m2(luminosity, distance) * float(stellar_spectrum.blackbody_fraction(*PAR_NM, temperature))
    earth = stellar_spectrum.bolometric_flux_w_m2(constants.SOLAR_LUMINOSITY, 1.0) * float(
        stellar_spectrum.blackbody_fraction(*PAR_NM, SUN_TEMPERATURE_K))
    return here / earth


def reductant_flux(body, area_m2):
    """`R`, mol/yr: Earth's, times the heat the body's interior and tide give off over Earth's, times `redox_factor`."""
    flux = (getattr(body, "internal_heat_flux_w_m2", None) or 0.0) + (getattr(body, "tidal_flux_w_m2", None) or 0.0)
    return EARTH_REDUCTANT_MOL_YR * (flux * area_m2 / EARTH_HEAT_W) * redox_factor(getattr(body, "mantle_delta_iw", None))


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


def productivity(body, area_m2):
    """`S`, mol/yr, of a body with a sunlit sea (0 for one without)."""
    if getattr(body, "hydrosphere", None) not in LIVING_HYDROSPHERES:
        return 0.0
    ocean = getattr(body, "ocean_fraction", None) or 0.0
    light = _clamp(active_light_ratio(body), *LIGHT_RANGE) ** LIGHT_EXPONENT
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


def profile_of(body, share=1.0, extra_store_mol=0.0):
    """
    The `Profile` of a rocky body (`NONE` for a giant, or a body missing the radius, gravity or air to read). `share` is the part of the body's food
    production that oxygenic life carries (GEN.1213: the classes of a body share its production, so the source `S` goes with their summed share);
    `extra_store_mol` is the oxygen, mol, that other reduced things the life made (the methane of the classes before the rise) take first, which adds to the store.
    """
    radius_km, gravity = getattr(body, "radius", None), getattr(body, "gravity", None)
    if getattr(body, "body_type", "t") == "g" or not radius_km or not gravity:
        return NONE
    area = area_m2(body)
    moles = air_moles(body)
    burnt, taken = burn(moles)
    source, reductant = productivity(body, area) * share, reductant_flux(body, area)
    store = EARTH_STORE_MOL * (area / EARTH_AREA_M2) * redox_factor(getattr(body, "mantle_delta_iw", None)) + taken + extra_store_mol
    surplus = source - reductant
    if surplus <= 0.0 or not burnt:
        return Profile(source, reductant, store, 0.0, 0.0, None)
    steady = min(surplus / WEATHERING_PER_YR, _fire_limit_mol(body, burnt))
    pascal = partials_pa(body, {**burnt, "O2": steady}).get("O2", 0.0)
    if pascal < MINIMUM_PA:
        return Profile(source, reductant, store, 0.0, 0.0, None)
    return Profile(source, reductant, store, steady, pascal, store / surplus / 1e9)


def rise_age_gyr(body, starts, extra_store_mol=0.0):
    """
    The age, Gyr, at which the oxygen of several oxygenic classes first rises on `body` (GEN.1213), or `None` when it never does. `starts` is a list of
    `(age_gyr, share)`: each class's producers begin at its age and carry its share of the food production. The surplus `S - R(age)` of the producers
    there are fills the store `B` (the first oxygen goes to the reduced things there), `S` going with the summed share of the classes begun and `R` with the
    heat the body gave off at that age (GEN.1463, `reductant_history_factor`); the oxygen rises when it is full. The surplus is summed in steps of `RISE_STEP_GYR` (the
    last step cut where the store fills) from the first producers' age; past the star's age the reductant is held at its value now.
    """
    starts = sorted((age, share) for age, share in starts if age is not None and share > 0.0)
    if not starts:
        return None
    full = profile_of(body, 1.0, extra_store_mol)
    if full.source_mol_yr <= 0.0:
        return None
    horizon = getattr(getattr(body, "star", None), "age", None) or 0.0       # the thermal history ends at the star's age; past it R is what it is now
    filled, age, carried, index = 0.0, starts[0][0], 0.0, 0
    while True:
        while index < len(starts) and starts[index][0] <= age + 1e-12:
            carried += starts[index][1]
            index += 1
        if age >= horizon and index >= len(starts):
            surplus = full.source_mol_yr * carried - full.reductant_mol_yr
            return None if surplus <= 0.0 else age + (full.store_mol - filled) / surplus / 1e9      # the surplus at the horizon runs on
        end = min(age + RISE_STEP_GYR, starts[index][0]) if index < len(starts) else age + RISE_STEP_GYR
        surplus = full.source_mol_yr * carried - full.reductant_mol_yr * reductant_history_factor(body, 0.5 * (age + end))
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
        log.warning(f"The carbon dioxide of a world did not settle against its sink: between {low:.4g} and {high:.4g} Pa after {REBALANCE_BISECTIONS} halvings (GEN.1589).")
    answer = min(high, room)
    if abs(answer - current) <= REBALANCE_TOLERANCE * max(current, MINIMUM_PA):
        store(current)          # the trials left the body at the last of them: the air it was given is stored again, solved as the trials were
        return False
    store(answer)
    return True


def store_air(body, air):
    """
    Gives `body` the air `air` (formula to partial pressure, Pa): the warming and surface temperature that mixture gives under the starlight the
    body absorbed (the closure of `surface_air`, so the stored greenhouse factor still follows from the stored air, haze included: GEN.1589), and the total pressure,
    molar mass, density and scale height of the new mixture. The doses, aerosols, colour and score that read the air are the caller's to make
    again (`generation.life`).
    """
    gravity = body.gravity * constants.EARTH_GRAVITY
    radiated = constants.STEFAN_BOLTZMANN_CONSTANT * body.surface_temperature ** 4 / (1.0 + body.greenhouse_factor)
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
