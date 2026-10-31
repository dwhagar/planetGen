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
delay and nothing before. A world with `S <= R` never holds any. The most the air can hold is `mixture_oxygen.BIOTIC_O2_MAX` of it (the
fire limit).

The numbers are those of the one world we know (`EARTH_*`, round figures from the literature of the oxygen cycle: Holland 2002; Lenton and
Catling and Claire 2005, from memory), used as the test point and the scale and not as a preset: a world
with another light, water, heat, redox and gravity gets another oxygen. Every constant is a stand-in in section 6.21 of
`docs/design/research/gap-filling-plan.md` with its tier and tweak level.

What it leaves out: the loss of oxygen to space (a few parts in 1e4 of the inventory a Gyr for a habitable-zone world, far under the
weathering), the hydrogen escape that adds oxygen, the phosphorus feedback and the land plants that raised Earth's burial (they would
change `S` by a factor of a few, inside the 1.7 to 4 that this module's L3 tweak already stands for).
"""

import math
from dataclasses import dataclass

from planetgen.physics import atmosphere, atmosphere_state, constants, gases, mixture, mixture_oxygen, stellar_spectrum

EARTH_O2_MOL = 3.7e19
"""float: The oxygen of Earth's air, mol (21 percent of 1.77e20 mol of air; 22.8 kPa at its gravity and area)."""

EARTH_BURIAL_MOL_YR = 1.0e13
"""float: The organic carbon Earth buries a year, mol (the source `S`; Hayes and Waldbauer 2006, about 1e13). T2, L3."""

EARTH_REDUCTANT_MOL_YR = 3.0e12
"""float: The oxygen Earth's volcanic and metamorphic reductants take a year, mol (`R`; Holland 2002, a few 1e12). T3, L3."""

EARTH_STORE_MOL = 5.0e21
"""float: The oxygen the ocean's iron and the crust's sulfides and iron take before the air holds any, mol (`B`), chosen so that the first
oxygen follows the first photosynthesis by 0.7 Gyr when `S` and `R` are Earth's (the Great Oxidation came some 0.3 to 1 Gyr after oxygenic
photosynthesis began). T9, L3."""

EARTH_AREA_M2 = 4.0 * math.pi * (constants.EARTH_RADIUS_KM * 1000.0) ** 2
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

SUN_TEMPERATURE_K = 5772.0

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
"""float: The oxygen, Pa, the air must hold for the wait for complex cells to start (about 1 percent of Earth's: Mills et al. 2014; T3, L0)."""

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


def _fire_limit_mol(moles):
    """The most oxygen, mol, the other gases of `moles` leave room for under `mixture_oxygen.BIOTIC_O2_MAX` of the molecules."""
    rest = sum(amount for formula, amount in moles.items() if formula != "O2")
    share = mixture_oxygen.BIOTIC_O2_MAX
    return share / (1.0 - share) * rest


def partials_pa(body, moles):
    """The partial pressures, Pa, of `moles` (formula to mol) on `body`: the weight of the column over the area, shared by the fractions."""
    column = sum(moles.values())
    weight = sum(amount * gases.gas(formula).molar_mass for formula, amount in moles.items())
    total = weight * body.gravity * constants.EARTH_GRAVITY / area_m2(body)
    return {formula: amount / column * total for formula, amount in moles.items() if amount > 0.0}


def profile_of(body):
    """The `Profile` of a rocky body (`NONE` for a giant, or a body missing the radius, gravity or air to read)."""
    radius_km, gravity = getattr(body, "radius", None), getattr(body, "gravity", None)
    if getattr(body, "body_type", "t") == "g" or not radius_km or not gravity:
        return NONE
    area = area_m2(body)
    moles = air_moles(body)
    burnt, taken = burn(moles)
    source, reductant = productivity(body, area), reductant_flux(body, area)
    store = EARTH_STORE_MOL * (area / EARTH_AREA_M2) * redox_factor(getattr(body, "mantle_delta_iw", None)) + taken
    surplus = source - reductant
    if surplus <= 0.0 or not burnt:
        return Profile(source, reductant, store, 0.0, 0.0, None)
    steady = min(surplus / WEATHERING_PER_YR, _fire_limit_mol(burnt))
    pascal = partials_pa(body, {**burnt, "O2": steady}).get("O2", 0.0)
    if pascal < MINIMUM_PA:
        return Profile(source, reductant, store, 0.0, 0.0, None)
    return Profile(source, reductant, store, steady, pascal, store / surplus / 1e9)


def oxygen_mol(profile, gyr_since_photosynthesis):
    """The oxygen, mol, that photosynthesis `gyr_since_photosynthesis` old has built: the steady amount once the store is spent, 0 before."""
    if profile.delay_gyr is None or gyr_since_photosynthesis < profile.delay_gyr:
        return 0.0
    return profile.steady_mol


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
    wanted = max(have, min(target_mol, _fire_limit_mol(burnt)))
    if taken <= 0.0 and wanted <= have * APPLY_RISE:
        return False
    burnt["O2"] = wanted
    air = partials_pa(body, burnt)
    if "O2" not in air or air["O2"] < MINIMUM_PA:
        return False
    store_air(body, air)
    return True


def store_air(body, air):
    """
    Gives `body` the air `air` (formula to partial pressure, Pa): the warming and surface temperature that mixture gives under the starlight the
    body absorbed (the closure of `surface_air`, so the stored greenhouse factor still follows from the stored air), and the total pressure,
    molar mass, density and scale height of the new mixture. The doses, aerosols, colour and score that read the air are the caller's to make
    again (`generation.life`).
    """
    gravity = body.gravity * constants.EARTH_GRAVITY
    absorbed = constants.STEFAN_BOLTZMANN_CONSTANT * body.surface_temperature ** 4 / (1.0 + body.greenhouse_factor)
    warming = atmosphere_state.greenhouse_warming(air, gravity, absorbed)
    state = atmosphere_state.air_state(air, warming.surface_temperature_k, gravity)
    body.surface_temperature = warming.surface_temperature_k
    body.greenhouse_factor = warming.factor
    body.atmospheric_pressure = state.total_pressure_pa
    body.atm_molar_density = state.mean_molar_mass_kg_mol
    body.atm_density = state.density_kg_m3
    body.scale_height = state.scale_height_m / 1000.0
    atmosphere.store_air(body, air, body.mantle_delta_iw)
