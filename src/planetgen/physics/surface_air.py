# planetgen/physics/surface_air.py

"""
The air of a body: the fixed point of its mixture, pressure, greenhouse and temperature (GEN.411, GEN.422;
docs/design/composition/surface-pressure-supply-law.md section "In the generator")
=========================================================================================================

The surface temperature depends on the air through the greenhouse, the air depends on the temperature through the freeze-out of its
condensable gases and the carbonate-silicate cycle, and the pressure depends on both. This module solves them together, in the order
GEN.422 fixes, repeating until the pressure and the temperature stop moving:

1. *Mixture.* The mixture draw (`mixture_draw.draw_mixture`) at the temperature and pressure of the step, over the element budget the
   interior has degassed (`surface_pressure.degassed_share` of each element's inventory) less what the sinks take: the sulfate and sulfide sinks of the crust take sulfur (`sulfur_balance`, which the ocean then dissolves
   some of, `mixture_sinks`), and the carbonate cycle takes the carbon that weathering locks up at the temperature
   (`surface_pressure.carbonate_pco2_bar`).
2. *Pressure.* The moles of gas are an air mass, `P = g M_air / (4 pi R^2)`, and the gas that stays is what the light-gas escape
   (`mixture_escape`) leaves. The escape stays out until the fixed point has settled without it; it then compares its loss to
   that settled pressure and the fixed point goes on with the loss in it (a loss compared to the pressure of the step before is
   compared to what stays already and takes its share again at every step).
3. *Freeze-out.* What a gas holds above its saturation pressure at the surface temperature is condensate (an ocean, a cap): the
   column loses it, exactly, with the mean molar mass it leaves (`freeze_out`); the retention regime of the recipe (`d`) caps the
   total.
4. *Density and scale height* of what is left (`atmosphere_state.air_state`).
5. *Greenhouse and temperature* of what is left (`atmosphere_state.greenhouse_warming`).

The step is damped, in the logarithm of the pressure and in the temperature, and it stops at `TOLERANCE`. Every draw is made once, at
the start, and the mixture is replayed from the same seed at every step: the regimes of the recipe are drawn again at each step's
temperature (they depend on it), and kept from the `FREEZE_AFTER` th step on so that a regime on the edge of its gate cannot flip the
solution back and forth. A giant, a brown dwarf and an envelope have no surface: their reference level is `REFERENCE_PRESSURE_BAR`,
and the mixture there is drawn from the layer temperature of the body, not from the temperature the step solves, so one step is the
solution. The vapour over a melt or an ice sets its own pressure.

What it does not do: the magma-ocean epoch, impacts, non-thermal escape of the heavy gases (the retention rows stand in), the seafloor
weathering of an ocean world, a tectonic switch of the heat loss, the helium-rain share of a giant's own heat (GEN.810: its temperature at the reference
level radiates the absorbed starlight and the cooling-law heat of `giant_draw.own_flux_w_m2`, GEN.816).
"""

import math
from dataclasses import dataclass, field

from planetgen.physics import (
    aerosols, atmosphere, atmosphere_state, clathrate_brine, constants, gas_equilibrium, gases, giant_draw, greenhouse, hydrosphere, interior_heat, mixture_draw,
    mixture_escape, mixture_rocky,
    mixture_sinks, mixture_subject, radiogenic, saturation, sulfur_balance, surface_pressure, waterworld_carbon,
)
from planetgen import tuning
from planetgen.generation.possibility_table import Subject
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.physics.mixture_base import bracket
from planetgen.util import draw

# --- numerical ----------------------------------------------------------------------------------------------------------------

MAX_STEPS = 60
"""int: The most steps of the fixed point. A body that has not settled by then keeps its last step and says so (`Air.converged`)."""

FREEZE_AFTER = 8
"""int: From this step on the recipe draw of the last step is kept (the regimes can no longer flip)."""

DAMPING = 0.6
"""float: The share of each step taken (in ln P and in T) to start with; a variable whose change has the sign opposite to the last
one's (the solution has been overshot) has its share halved, down to `MIN_DAMPING`."""

MIN_DAMPING = 0.05

MAX_OCEAN_FLIPS = 4
"""int: How many times the ocean may freeze and thaw during the fixed point before the body is taken to cycle between the two. A world
on the edge of its liquid range has no state that is both (an ocean weathers the CO2 away and freezes; the ice leaves the CO2 to
build up and thaws), so it is solved as one that keeps its ocean: the weathering law is continuous in the temperature (Lehmer et al.
2020), and its own temperature dependence takes the weathering away where it is cold."""

MAX_GAIN = 20.0
"""float: The most an extrapolated step may exceed the damped one (`_step`)."""

TOLERANCE_LN_P = 1e-3
"""float: The relative change of the pressure under which the fixed point has settled (the greenhouse solver's own noise, 0.01 K,
moves a carbonate-limited pressure by a few parts in 1e4)."""

TOLERANCE_K = 0.1
"""float: The change of the temperature (K) under which the fixed point has settled. The greenhouse closure's own fixed point stops at
1e-3 K between its steps, and in a wet air near the moist-greenhouse limit its steps contract by 0.98, so its answer is good to about
0.1 K and a tighter outer tolerance only follows its noise."""

AIRLESS_PRESSURE_BAR = 1.0e-12
"""float: The pressure the subject of an airless body's redox is built at (the mixture subject wants a positive one; the redox rows do not
read it)."""

MAX_EQUILIBRIUM_STALLS = 3
"""int: How many times the fixed point may step back after the gas equilibrium fails to converge (see `solve`)."""

VACUUM_PRESSURE_PA = 1.0
"""float: Below this pressure (Pa) the volcanic gas equilibrium does not converge (it is a solver of a dense gas, and a fixed point that
falls this far has no air to equilibrate); a body whose pressure falls here keeps only an exosphere, which is taken as airless. A
planetGen stand-in (tier T5, L0): the exosphere of Mercury or the Moon is 1e-10 Pa and under."""

MAX_SURFACE_TEMPERATURE_K = 2000.0
"""float: The hottest surface the fixed point may reach. A steam or carbon dioxide greenhouse that runs away stops where the surface is
molten rock radiating through its own windows (a radiative limit near 1500 to 1800 K, recalled from the magma-ocean literature, not
read), and the equilibrium of the gas does not converge far above it. A planetGen stand-in (tier T5, L1): the iteration of an
unbounded greenhouse reached 3700 K at 3700 bar for a hot, heavy body."""

MIN_ABSORBED_FLUX_W_M2 = 1.0e-30
"""float: A body in the dark absorbs no light, and the greenhouse closure wants a positive flux; its temperature is the background's."""

FIRST_PRESSURE_PA = 1.0e5
"""float: The pressure the first step's mixture is solved at (bar scale of the equilibrium; the fixed point moves it)."""

STORE_MIN_FRACTION = 1e-9
"""float: A gas below this mole fraction of the air is not kept (the pressure is the sum of the gases kept)."""

SATURATION_TOLERANCE = 1e-9
"""float: How far above its saturation pressure a gas may stand before the freeze-out takes it down."""

FREEZE_PASSES = 100
"""int: The most passes of the freeze-out over the condensable gases."""

# --- the design constants -----------------------------------------------------------------------------------------------------

REFERENCE_PRESSURE_BAR = 1.0
"""float: The pressure of the reference level of a body with no surface: the 1 bar level that the radius and the temperature of a giant
refer to (Jupiter's equatorial radius is the 1 bar one). A convention, not a stand-in."""

GY_S = 3.15576e16
"""float: One billion Julian years, s."""

UNKNOWN_OXIDISED_SHARE = 0.5
"""float: The share of the sulfur taken to be oxidised at the first step, before there is a gas to read it from (`sulfur_balance.oxidised_share_of`
reads the previous step's gas from the second step on). A starting guess of the fixed point, not a stand-in: it changes only the first step."""

LONG_LIVED_FLUX_FRACTION = 0.5
"""float: A body whose interior still carries this share of the heat flow of the Earth twin out of its surface is geologically alive
(the redox of its mantle has moved on to the later upper mantle's: `mixture_rocky.draw_redox`). Stand-in (tier: educated guess from
inferred physics, tweak level L0)."""

RELATIVE_HUMIDITY = 0.77
"""float: The relative humidity of the air over a surface with water (the fixed value of Manabe and Wetherald 1967, 0.77, the Earth's
mean; tier: verified, tweak level L0): the water vapour stands at this share of its saturation pressure while the ocean or ice lasts."""

GIANT_FAMILIES = ("giant", "brown-dwarf", "hot-jupiter", "envelope")
"""tuple: The families with a reference level and no surface."""

VAPOUR_FAMILIES = ("lava", "frost")
"""tuple: The families whose pressure is the vapour pressure of the surface (`MixtureDraw.pressure_bar`)."""


# --- the result -----------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Setting:
    """
    What the fixed point needs of a body.

    Attributes:
        gravity_m_s2, radius_m (float): The surface gravity and the radius.
        mantle_kg (float): The mass of the mantle that holds the delivered volatiles.
        released_j (float): The heat the interior has carried out (`interior_heat.Heat.released_j`); 0 for a body with no interior.
        heat_flux_w_m2 (float): Its heat flow now (the outgassing flux of the carbonate balance).
        ocean_cover (float): The share of the surface a liquid ocean could cover (`hydrosphere`'s hypsometry), 0 for none.
        water_kg (float): The water at the surface (ocean, ice, the part of the mantle's that the surface holds), the reservoir the water
            vapour of the air stands in contact with.
        absorbed_flux_w_m2 (float): The starlight the surface absorbs, averaged over the sphere (`greenhouse.absorbed_flux_w_m2`).
        internal_flux_w_m2 (float): A giant's own heat flux at its reference level (`giant_draw.own_flux_w_m2`, GEN.816), radiated with the
            absorbed starlight; 0 for a body with a surface, whose heat flow is too small to warm the air.
        age_gy (float): The body's age (its star's, else the reference age of a starless body), Gy: the span the degassed sulfur arrived over.
        released_at_seal_j (float): The heat the interior had carried out by `waterworld_carbon.SEAL_TIME_GY` (or the body's age, if younger):
            what has degassed into the closed system of a world whose ocean seals the seafloor.
        mantle_water_kg (float): The water the mantle holds (the bulk water less the surface's); volcanoes bring up at most
            `tuning.MANTLE_OUTGAS_SHARE` of it, and the air draw takes no more hydrogen than that (GEN.837); infinity for a setting that
            does not say.
        bulk_water_fraction (float): The water mass fraction of the whole body (an icy body, from `hydrosphere.ICY_MIN_WATER_FRACTION`,
            holds guest gases in the cages of its ice skin).
        subject_of (callable): `(ts_k, pressure_bar, ocean_present, long_lived_geology) -> Subject` of the mixture draw.
        giant (bool): The body has no surface (a giant).
        haze_draws (HazeDraws or None): The body's haze draws (`aerosols.draw_haze_draws` of `aerosols.haze_stream_of`), `None` for a setting
            whose air makes no haze (the haze takes sunlight aloft, GEN.851).
        star_temperature_k (float or None): The star's effective temperature, whose blackbody the haze's share of the sunlight is weighted by.
    """

    gravity_m_s2: float
    radius_m: float
    mantle_kg: float
    released_j: float
    heat_flux_w_m2: float
    ocean_cover: float
    water_kg: float
    absorbed_flux_w_m2: float
    released_at_seal_j: float
    subject_of: object
    mantle_water_kg: float = math.inf
    age_gy: float = tuning.ROGUE_EARTH_REFERENCE_AGE_GY
    bulk_water_fraction: float = 0.0
    giant: bool = False
    internal_flux_w_m2: float = 0.0
    haze_draws: object = None
    star_temperature_k: object = None

    @property
    def area_m2(self):
        return 4.0 * math.pi * self.radius_m ** 2

    @property
    def radiated_flux_w_m2(self):
        """What the body radiates in balance: the starlight it absorbs and its own heat flux."""
        return self.absorbed_flux_w_m2 + self.internal_flux_w_m2

    @property
    def effective_temperature_k(self):
        return greenhouse.effective_temperature_k(self.radiated_flux_w_m2)


@dataclass(frozen=True)
class Air:
    """
    The solved air.

    Attributes:
        mixture (MixtureDraw): The last step's mixture draw (`family`, `drawn`, `redox_diw`, the recipe).
        partials_pa (dict): Gas-table formula -> partial pressure at the surface, Pa, after the freeze-out (only the gases kept).
        pressure_pa (float): Their sum.
        temperature_k (float): The surface temperature.
        warming (Warming or None): The greenhouse of the last step (`None` for a body with no air).
        haze (Haze): The haze the last step took sunlight aloft with (`aerosols.NO_HAZE` for none).
        state (AirState or None): The density, mean molar mass and scale height.
        condensed_kg (dict): Formula -> the mass the freeze-out took out of the air, kg (the ocean and the caps).
        steps (int): The steps taken.
        converged (bool): The pressure and the temperature settled to `TOLERANCE_LN_P` and `TOLERANCE_K`.
        ocean_present (bool): The last step's ocean flag.
        constants (dict): The design constants drawn (the efficiencies, the cycle's draws, the sulfur lifetimes, the skin of ice).
    """

    mixture: object
    partials_pa: dict
    pressure_pa: float
    temperature_k: float
    warming: object
    state: object
    condensed_kg: dict
    steps: int
    converged: bool
    ocean_present: bool
    constants: dict = field(default_factory=dict)
    haze: object = aerosols.NO_HAZE

    @property
    def greenhouse_factor(self):
        """`G = (T_s / T_e)^4 - 1` of the solved air (0 for none)."""
        return 0.0 if self.warming is None else self.warming.factor


# --- freeze-out -----------------------------------------------------------------------------------------------------------------

def saturation_limit_pa(formula, temperature_k):
    """The most `formula` can stand at over a surface of `temperature_k`, Pa: its saturation pressure; infinite above its critical
    point and for a gas with no curve."""
    if not saturation.has(formula):
        return math.inf
    return saturation.pressure_pa(formula, temperature_k)


def _moles_at_partial(limit_pa, a, mass, rest_moles, rest_weight):
    """The moles of a gas (molar mass `mass`) for which its partial pressure in the column is `limit_pa`, the other gases being fixed:
    the positive root of `a M n^2 + (a W - p) n - p N = 0` (`a = g / A`, `W` and `N` the weight and moles of the rest)."""
    b = a * rest_weight - limit_pa
    discriminant = math.sqrt(b * b + 4.0 * a * mass * limit_pa * rest_moles)
    return (discriminant - b) / (2.0 * a * mass) if b <= 0.0 else 2.0 * limit_pa * rest_moles / (discriminant + b)


def freeze_out(moles, gravity_m_s2, area_m2, temperature_k, reservoirs=None):
    """
    The air that stands over a surface once the gases above their saturation pressure have condensed and the gases in contact with a
    reservoir have come to their vapour pressure.

    A hydrostatic well-mixed column has the pressure `P = g M / A` and a gas the partial pressure `x_i P`. A gas that stands above its
    saturation pressure `p_sat` loses moles until `x_i P = p_sat`, which changes the mean molar mass and so the partial pressures of the
    rest; a gas in contact with a reservoir (the water of an ocean or an ice cap) is brought to `humidity * p_sat` from either side,
    taking moles from the reservoir while it lasts. Each is set in turn by the positive root of the quadratic its condition makes
    (`_moles_at_partial`) until none moves; a gas with no saturation curve, or above its critical temperature, is left alone, and a
    gas in contact with a reservoir there takes the whole reservoir.

    Args:
        moles (dict): Formula -> moles of the gas in the column.
        gravity_m_s2, area_m2 (float): The surface gravity and the area.
        temperature_k (float): The surface temperature.
        reservoirs (dict, optional): Formula -> `(humidity, moles in the reservoir)`.

    Returns:
        tuple: `(moles left, net moles condensed, moles left in the reservoirs)`: the first two by formula (a gas the reservoir fed has a
        negative net), the last a dict by reservoir formula.
    """
    reservoirs = dict(reservoirs or {})
    left = {formula: count for formula, count in moles.items()}
    for formula in reservoirs:
        left.setdefault(formula, 0.0)
    start = dict(left)
    a = gravity_m_s2 / area_m2
    mass = {formula: gases.gas(formula).molar_mass for formula in left}
    limits = {formula: saturation_limit_pa(formula, temperature_k) for formula in left}
    available = {formula: moles_ for formula, (_, moles_) in reservoirs.items()}
    for _ in range(FREEZE_PASSES):
        moved = 0.0
        for formula in sorted(left):
            limit = limits[formula]
            humid = formula in reservoirs
            if limit == math.inf:
                if humid and available[formula] > 0.0:
                    left[formula] += available[formula]
                    available[formula] = 0.0
                    moved = 1.0
                continue
            target = reservoirs[formula][0] * limit if humid else limit
            total = math.fsum(left.values())
            weight = math.fsum(left[f] * mass[f] for f in left)
            if total <= 0.0:
                continue
            partial = left[formula] / total * a * weight
            too_much = partial > target * (1.0 + SATURATION_TOLERANCE)
            too_little = humid and available[formula] > 0.0 and partial < target * (1.0 - SATURATION_TOLERANCE)
            if not (too_much or too_little):
                continue
            new = _moles_at_partial(target, a, mass[formula], total - left[formula], weight - left[formula] * mass[formula])
            if new > left[formula]:
                if not humid:
                    continue
                new = left[formula] + min(new - left[formula], available[formula])
                available[formula] -= new - left[formula]
            else:
                if humid:
                    available[formula] += left[formula] - new
            moved = max(moved, abs(new - left[formula]) / max(left[formula], 1.0e-300))
            left[formula] = new
        if moved < 1.0e-12:
            break
    condensed = {formula: start[formula] - left[formula] for formula in left}
    return left, condensed, available


# --- the solution ----------------------------------------------------------------------------------------------------------------

def liquid_ocean(ocean_cover, temperature_k, pressure_pa):
    """Whether a surface with water that could cover `ocean_cover` of it has a liquid ocean: between the freezing point and the critical
    point of water, and under at least its vapour pressure (above the boiling point at the pressure there is none)."""
    if ocean_cover <= 0.0 or not surface_pressure.LIQUID_WATER_K[0] <= temperature_k < surface_pressure.LIQUID_WATER_K[1]:
        return False
    return saturation.pressure_pa("H2O", temperature_k) <= pressure_pa


def _weathers(setting, water_left, temperature_k):
    """Whether silicate weathering takes carbon at this step: the body has surface water that has not all gone to vapour, and its surface is
    under the critical temperature of water. The cold side is not gated: the weathering law falls off exponentially with the temperature
    (`surface_pressure.carbonate_pco2_bar`), and a hard gate at the freezing point would make the carbon of a frozen world jump to all that
    has been degassed in one step instead of building up until it thaws."""
    return (not setting.giant and setting.ocean_cover > 0.0 and water_left > 0.0
            and temperature_k < surface_pressure.LIQUID_WATER_K[1])


def _dissolve_sulfur(moles, temperature_k, setting):
    """`moles` (formula -> moles in the column) less the sulfur gases an ocean of all the surface water dissolves (Henry's law,
    `mixture_sinks.sulfur_remaining_fraction`), and the moles each lost."""
    column_kg_m2 = setting.water_kg / setting.area_m2
    left, dissolved = dict(moles), {}
    for formula in mixture_sinks.SULFUR_GASES:
        if left.get(formula, 0.0) > 0.0:
            kept = mixture_sinks.sulfur_remaining_fraction(formula, temperature_k, setting.gravity_m_s2, column_kg_m2)
            dissolved[formula] = left[formula] * (1.0 - kept)
            left[formula] *= kept
    return left, dissolved


_STORED_OF = {gases.stored(name).formula: name for name in atmosphere.SPECIES}
"""dict: Gas-table formula -> the stored species name of the ten the planet stores (`atmosphere.SPECIES`)."""


def _hold_in_hydrate(partials_pa, temperature_k, setting, skin_m):
    """
    `partials_pa` less what the water cages of an icy cold body's ice skin hold (`clathrate_brine.hydrate_hold_pa`, GEN.762): a body
    with at least `hydrosphere.ICY_MIN_WATER_FRACTION` of its mass in water whose surface is under the ice point holds the guest gases
    that stand above their dissociation pressure, in a skin of `skin_m` metres of ice.
    """
    if setting.bulk_water_fraction < hydrosphere.ICY_MIN_WATER_FRACTION or temperature_k >= clathrate_brine.ICE_POINT_K:
        return partials_pa
    listed = {_STORED_OF[formula]: pressure for formula, pressure in partials_pa.items() if formula in _STORED_OF}
    held = clathrate_brine.hydrate_hold_pa(listed, temperature_k, setting.gravity_m_s2, skin_m)
    return {formula: pressure - held.get(_STORED_OF.get(formula), 0.0) for formula, pressure in partials_pa.items()}


def _damped(damping, residual, previous, tolerance):
    """The damping of one variable after a step: halved (down to `MIN_DAMPING`) when the change has the sign opposite to the last one's
    (the solution was overshot) by more than the tolerance, so noise under the tolerance does not slow the others."""
    if residual * previous < 0.0 and abs(residual) > tolerance:
        return max(MIN_DAMPING, 0.5 * damping)
    return damping


def _step(damping, residual, history):
    """
    The step of one variable: the damped residual, or, when the last three residuals fall by the same ratio each (a creep toward the
    solution, as when the water vapour and the temperature feed each other), the damped step times `1 / (1 - q)` for that ratio `q`,
    the distance of a geometric series to its limit, at most `MAX_GAIN` times. A residual that does not shrink by a steady ratio (a
    swing, noise, a drift away) is never extrapolated.

    Args:
        damping (float): The damping of the variable.
        residual (float): The residual at the variable (`f(value) - value`).
        history (list): The residuals of the steps before, oldest first (the last two are used).
    """
    plain = damping * residual
    if len(history) < 2 or 0.0 in history[-2:] or residual == 0.0:
        return plain
    first, second = history[-2:]
    q_before, q = second / first, residual / second
    if not (0.0 < q_before < 1.0 and 0.0 < q < 1.0):
        return plain
    if not 0.7 <= (1.0 - q_before) / (1.0 - q) <= 1.4:
        return plain
    return plain * min(MAX_GAIN, 1.0 / (1.0 - q))


def long_lived_geology(heat_flux_w_m2):
    """Whether a body's interior is still carrying `LONG_LIVED_FLUX_FRACTION` of the Earth twin's heat flow out."""
    return heat_flux_w_m2 >= LONG_LIVED_FLUX_FRACTION * surface_pressure.HEAT_FLUX_EARTH_W_M2


@dataclass(frozen=True)
class WaterWorld:
    """
    The carbon sink of a body with an ocean (`waterworld_carbon`): which of the three sinks it has and the draws they take.

    Attributes:
        sealed (bool): The seafloor pressure of the ocean reaches the seal pressure (the mantle no longer feeds the air).
        draws (waterworld_carbon.WaterWorldDraws): The seal pressure and the leached crust thickness.
        law (waterworld_carbon.SeafloorLaw): The seafloor weathering law.
    """

    sealed: bool
    draws: object
    law: object


def water_world_of(setting, design):
    """The `WaterWorld` of a body (its draws recorded in `design`): sealed when the pressure of its ocean's column on the seafloor reaches
    the drawn seal pressure."""
    draws = waterworld_carbon.draw_waterworld(design)
    law = waterworld_carbon.draw_seafloor_law(design)
    column = setting.water_kg / setting.area_m2 if setting.ocean_cover > 0.0 else 0.0
    sealed = waterworld_carbon.is_sealed(waterworld_carbon.seafloor_pressure_pa(column, setting.gravity_m_s2), draws.seal_pa)
    return WaterWorld(sealed=sealed, draws=draws, law=law)


def _shares(setting, efficiencies, cycle, sulfur, ocean, weathers, temperature_k, water_world):
    """The `supply` callable of one step: the share of each delivered element that the surface holds at `temperature_k`. `sulfur` is
    `(lifetimes, age_s, oxidised_share)` for the sulfate and sulfide balance (`sulfur_balance`)."""
    released = setting.released_at_seal_j if water_world.sealed else setting.released_j
    base = surface_pressure.degassed_share(released, setting.mantle_kg)
    moles_per_100_g = setting.mantle_kg / 0.1
    mantle_hydrogen_mol = 2.0 * tuning.MANTLE_OUTGAS_SHARE * setting.mantle_water_kg / gases.gas("H2O").molar_mass
    crust_kg_m2 = sulfur_balance.crust_column_kg_m2(setting.released_j, setting.radius_m)

    def shares(inventory):
        shares_ = {element: base * efficiency for element, efficiency in efficiencies.items()}
        if ocean:   # GEN.524 part 2: hydrogen chloride and fluoride are so soluble that an ocean takes all the volcanoes give (the air keeps none)
            shares_["Cl"] = shares_["F"] = 0.0
        degassed_sulfur = inventory.get("S", 0.0) * shares_.get("S", 0.0) * moles_per_100_g
        if degassed_sulfur > 0.0:
            lifetimes, age_s, oxidised_share = sulfur
            kept = sulfur_balance.balance(degassed_sulfur / setting.area_m2, age_s=age_s, ocean_present=ocean, oxidised_share=oxidised_share,
                                          crust_kg_m2=crust_kg_m2, lifetimes=lifetimes)
            shares_["S"] *= kept.air_mol_m2 * setting.area_m2 / degassed_sulfur
        degassed_hydrogen = inventory.get("H", 0.0) * shares_.get("H", 0.0) * moles_per_100_g
        if degassed_hydrogen > mantle_hydrogen_mol:  # GEN.837: volcanoes bring up no more hydrogen than MANTLE_OUTGAS_SHARE of the mantle's water holds
            shares_["H"] *= mantle_hydrogen_mol / degassed_hydrogen
        carbon = inventory.get("C", 0.0) * shares_.get("C", 0.0) * moles_per_100_g
        if weathers and carbon > 0.0:
            supply_bar = carbon * gases.gas("CO2").molar_mass * setting.gravity_m_s2 / setting.area_m2 / 1.0e5
            kind = waterworld_carbon.sink_kind(water_world.sealed, 1.0 - setting.ocean_cover)
            if kind == waterworld_carbon.CONTINENTAL_WEATHERING:
                kept_bar = surface_pressure.carbonate_pco2_bar(supply_bar, setting.heat_flux_w_m2, 1.0 - setting.ocean_cover,
                                                               lambda p_bar: temperature_k, *cycle)
            elif kind == waterworld_carbon.SEAFLOOR_WEATHERING:
                kept_bar = waterworld_carbon.seafloor_pco2_bar(supply_bar, setting.heat_flux_w_m2, lambda p_bar: temperature_k,
                                                               water_world.law, cycle[2])
            else:
                sodium, calcium = waterworld_carbon.leached_cations_mol_m2(water_world.draws.altered_km)
                kept_bar = waterworld_carbon.ocean_pco2(carbon / setting.area_m2, temperature_k, setting.water_kg / setting.area_m2,
                                                        setting.gravity_m_s2, sodium_mol_m2=sodium, calcium_mol_m2=calcium).pco2_bar
            shares_["C"] *= min(1.0, kept_bar / supply_bar)
        return shares_

    return shares


def _draw_stream(seed):
    return draw.bound(draw.Stream(seed))


def _airless(setting, mixture, steps, design, ocean):
    return Air(mixture=mixture, partials_pa={}, pressure_pa=0.0,
               temperature_k=max(constants.COSMIC_BACKGROUND_TEMPERATURE_K, setting.effective_temperature_k), warming=None,
               state=None, condensed_kg={}, steps=steps, converged=True, ocean_present=ocean, constants=design)


def _column_moles(result, mantle_kg):
    """The moles of each gas in the column of a mixture that has `result.moles` per 100 g of the mantle."""
    total = result.moles * mantle_kg / 0.1
    return {formula: fraction * total for formula, fraction in result.mole_fractions.items()}


def solve(setting, start_temperature_k):
    """
    Solves the air of a body (see the module text).

    Args:
        setting (Setting): The body.
        start_temperature_k (float): The temperature the first step starts from (the airless one).

    Returns:
        Air: The air. Draws from the ambient stream once (the design constants and the seed of the replay).
    """
    design = {}
    efficiencies = surface_pressure.draw_efficiencies(design)
    cycle = surface_pressure.draw_cycle(design)
    lifetimes = sulfur_balance.draw_lifetimes(design)
    oxidised_share = UNKNOWN_OXIDISED_SHARE
    water_world = water_world_of(setting, design)
    skin_m = bracket(design, "hydrate_skin_m", (*clathrate_brine.HYDRATE_SKIN_RANGE_M, "log"))
    seed = draw.getrandbits(64)
    geology = long_lived_geology(setting.heat_flux_w_m2)

    temperature = max(constants.COSMIC_BACKGROUND_TEMPERATURE_K, start_temperature_k)
    ln_pressure = math.log(REFERENCE_PRESSURE_BAR * 1.0e5 if setting.giant else FIRST_PRESSURE_PA)
    recipe_subject = None
    last = None
    water_moles = setting.water_kg / gases.gas("H2O").molar_mass
    reservoir_left = water_moles
    damping = [DAMPING, DAMPING]
    previous_residual = (0.0, 0.0)
    history = [[], []]
    flips, last_ocean, cycling = 0, False, False
    stalls = 0
    escape_bar = None       # the pressure the light-gas escape compares its loss to: set once the air has settled without it
    for step in range(1, MAX_STEPS + 1):
        pressure = math.exp(ln_pressure)
        if cycling:
            ocean = True
        else:
            ocean = (not setting.giant) and reservoir_left > 0.0 and (step == 1 and setting.ocean_cover > 0.0
                                                                       or liquid_ocean(setting.ocean_cover, temperature, pressure))
            flips += ocean != last_ocean and step > 1
            last_ocean = ocean
            cycling = flips >= MAX_OCEAN_FLIPS
            ocean = ocean or cycling
        weathers = _weathers(setting, reservoir_left, temperature)
        subject = setting.subject_of(temperature, pressure / 1.0e5, ocean, geology)
        shares = _shares(setting, efficiencies, cycle, (lifetimes, setting.age_gy * GY_S, oxidised_share), ocean, weathers, temperature,
                         water_world)
        try:
            with _draw_stream(seed):
                mixture = mixture_draw.draw_mixture(subject, recipe_subject=recipe_subject, supply=shares, escape=escape_bar is not None,
                                                    escape_column_bar=escape_bar)
        except gas_equilibrium.EquilibriumError:
            if last is None:
                raise
            if pressure <= VACUUM_PRESSURE_PA:
                return _airless(setting, last[0], step, design, ocean)
            # The equilibrium of the gas can stall at a pressure and temperature the fixed point only passed through: step back to the
            # middle of the way to the last good air and try there, a few times, before calling it a fault.
            stalls += 1
            if stalls > MAX_EQUILIBRIUM_STALLS:
                raise
            ln_pressure = 0.5 * (ln_pressure + math.log(last[2]))
            temperature = 0.5 * (temperature + last[3])
            continue
        if step == FREEZE_AFTER:
            recipe_subject = subject
        family = mixture.family
        if mixture.mole_fractions:
            oxidised_share = sulfur_balance.oxidised_share_of(mixture.mole_fractions)
        if family in ("airless", "coma") or not mixture.mole_fractions:
            return _airless(setting, mixture, step, design, ocean)

        condensed_kg = {}
        if family in GIANT_FAMILIES:
            new_pressure = REFERENCE_PRESSURE_BAR * 1.0e5
            partials = {f: x * new_pressure for f, x in mixture.mole_fractions.items()}
        elif family in VAPOUR_FAMILIES:
            new_pressure = mixture.pressure_bar * 1.0e5
            partials = {f: x * new_pressure for f, x in mixture.mole_fractions.items()}
        else:
            moles = _column_moles(mixture, setting.mantle_kg)
            reservoirs = {"H2O": (RELATIVE_HUMIDITY, water_moles)} if water_moles > 0.0 else None
            left, condensed, remaining = freeze_out(moles, setting.gravity_m_s2, setting.area_m2, temperature, reservoirs)
            reservoir_left = remaining.get("H2O", 0.0)
            if ocean:
                left, dissolved = _dissolve_sulfur(left, temperature, setting)
                for formula, moles_ in dissolved.items():
                    condensed[formula] = condensed.get(formula, 0.0) + moles_
            condensed_kg = {f: n * gases.gas(f).molar_mass for f, n in condensed.items() if n > 0.0}
            weight = math.fsum(n * gases.gas(f).molar_mass for f, n in left.items())
            total = math.fsum(left.values())
            new_pressure = surface_pressure.pressure_pa(weight, setting.gravity_m_s2, setting.radius_m)
            partials = {f: n / total * new_pressure for f, n in left.items()} if total > 0.0 else {}
            cap = _retention_cap_pa(mixture, subject)
            if cap is not None and new_pressure > cap:
                partials = {f: p * cap / new_pressure for f, p in partials.items()}
                new_pressure = cap
            partials = _hold_in_hydrate(partials, temperature, setting, skin_m)
            new_pressure = math.fsum(partials.values())
        floor = new_pressure * STORE_MIN_FRACTION
        partials = {f: p for f, p in partials.items() if p >= floor and p > 0.0}
        new_pressure = math.fsum(partials.values())
        if new_pressure <= 0.0:
            return _airless(setting, mixture, step, design, ocean)

        haze = aerosols.NO_HAZE if setting.haze_draws is None else aerosols.haze_of(partials, setting.gravity_m_s2, temperature, setting.haze_draws,
                                                                                  setting.star_temperature_k)
        # the haze takes a share of the starlight; the flux the air radiates also holds a giant's own heat, so the share of that flux is smaller
        aloft_share = haze.absorbed_share * setting.absorbed_flux_w_m2 / setting.radiated_flux_w_m2
        warming = atmosphere_state.greenhouse_warming(partials, setting.gravity_m_s2, setting.radiated_flux_w_m2, aloft_share, haze.pressure_pa)
        new_temperature = min(max(constants.COSMIC_BACKGROUND_TEMPERATURE_K, warming.surface_temperature_k), MAX_SURFACE_TEMPERATURE_K)
        change_p = abs(math.log(new_pressure) - ln_pressure)
        change_t = abs(new_temperature - temperature)
        last = (mixture, partials, new_pressure, new_temperature, warming, condensed_kg, ocean, haze)
        if setting.giant or (change_p < TOLERANCE_LN_P and change_t < TOLERANCE_K):
            if setting.giant or escape_bar is not None or family not in mixture_escape.FAMILIES:
                break
            # The air has settled without the escape. The loss acts on the air that stands, so it is compared to this pressure (the pressure of a step
            # that has lost its share already would lose it again at every step, and the first steps' reducing mixtures, which the settled air does not
            # keep, would take the whole air); the fixed point goes on with the loss in it.
            escape_bar = new_pressure / 1.0e5
            ln_pressure, temperature = math.log(new_pressure), new_temperature
            damping, previous_residual, history = [DAMPING, DAMPING], (0.0, 0.0), [[], []]
            continue
        residual = (math.log(new_pressure) - ln_pressure, new_temperature - temperature)
        damping = [_damped(d, r, q, tolerance) for d, r, q, tolerance in zip(damping, residual, previous_residual,
                                                                            (TOLERANCE_LN_P, TOLERANCE_K))]
        previous_residual = residual
        step_p = _step(damping[0], residual[0], history[0])
        step_t = _step(damping[1], residual[1], history[1])
        history[0].append(residual[0])
        history[1].append(residual[1])
        if setting.giant or family in VAPOUR_FAMILIES:
            ln_pressure = math.log(new_pressure)
        else:
            ln_pressure += step_p
        temperature = max(constants.COSMIC_BACKGROUND_TEMPERATURE_K, temperature + step_t)
    mixture, partials, pressure, temperature_k, warming, condensed_kg, ocean, haze = last
    state = atmosphere_state.air_state(partials, temperature_k, setting.gravity_m_s2)
    converged = setting.giant or (change_p < TOLERANCE_LN_P and change_t < TOLERANCE_K)
    return Air(mixture=mixture, partials_pa=partials, pressure_pa=pressure, temperature_k=temperature_k, warming=warming, state=state,
               condensed_kg=condensed_kg, steps=step, converged=converged, ocean_present=ocean, constants=design, haze=haze)


def _retention_cap_pa(mixture, subject):
    """The most pressure the retention regime of the recipe (`d6`) lets the air keep, Pa; `None` for no cap."""
    recipe = mixture.recipe
    if recipe is None or "d" not in recipe.groups:
        return None
    limit = RECIPES.derive("d6", Subject(subject.body_class, {**subject.quantities, "drawn.d": recipe.groups["d"].name}))["p_max_bar"]
    return None if limit is None else limit * 1.0e5


# --- a body ---------------------------------------------------------------------------------------------------------------------

def mantle_redox_diw(body, star_distance_au, water_lost):
    """
    The oxygen fugacity of the mantle of a body whose air is not drawn (airless by its class, or by loss), log10 units above
    iron-wuestite: the recipe's redox rows (`mixture_rocky.draw_redox`) at the body's own mass and age. `None` for a giant.
    """
    if body.body_type == "g":
        return None
    subject = mixture_subject.subject_of_body(body, ts_k=body.surface_temperature, pressure_bar=AIRLESS_PRESSURE_BAR,
                                              ocean_present=False, long_lived_geology=False, water_lost=water_lost,
                                              star_distance_au=star_distance_au)
    return mixture_rocky.draw_redox(subject, mixture_draw.draw_recipe(subject), {})


def setting_of(body, *, albedo, starlight_w_m2, star_distance_au, heat, surface_w, water_lost):
    """
    The `Setting` of a planet or moon whose class has an air: from its mass, radius, gravity and bulk composition, the heat its
    interior has released (`interior_heat.Heat`, `None` for a body with no interior), the surface water (`hydrosphere.draw_surface_water`)
    and the starlight it absorbs.
    """
    giant = body.body_type == "g"
    gravity = body.gravity * constants.EARTH_GRAVITY
    rock = radiogenic.rock_of(body)
    flux = 0.0 if heat is None else heat.flux_w_m2
    age_gy = getattr(getattr(body, "star", None), "age", None) or tuning.ROGUE_EARTH_REFERENCE_AGE_GY

    def subject_of(ts_k, pressure_bar, ocean_present, long_lived_geology):
        return mixture_subject.subject_of_body(body, ts_k=ts_k, pressure_bar=pressure_bar, ocean_present=ocean_present,
                                               long_lived_geology=long_lived_geology, water_lost=water_lost,
                                               star_distance_au=star_distance_au)

    return Setting(
        gravity_m_s2=gravity, radius_m=body.radius * 1000.0, mantle_kg=0.0 if giant else rock.silicate * body.mass,
        released_j=0.0 if heat is None else heat.released_j, heat_flux_w_m2=flux,
        ocean_cover=0.0 if giant else hydrosphere.ocean_cover(body, surface_w), water_kg=0.0 if giant else surface_w * body.mass,
        absorbed_flux_w_m2=max(greenhouse.absorbed_flux_w_m2(starlight_w_m2, albedo), MIN_ABSORBED_FLUX_W_M2),
        released_at_seal_j=interior_heat.released_by(body, min(age_gy, waterworld_carbon.SEAL_TIME_GY)),
        mantle_water_kg=0.0 if giant else max(0.0, (body.water_mass_fraction or 0.0) * body.mass - surface_w * body.mass),
        subject_of=subject_of, age_gy=age_gy,
        bulk_water_fraction=0.0 if giant else (body.water_mass_fraction or 0.0), giant=giant,
        internal_flux_w_m2=giant_draw.own_flux_w_m2(body) if giant else 0.0,
        haze_draws=aerosols.draw_haze_draws(aerosols.haze_stream_of(body)),
        star_temperature_k=getattr(getattr(body, "star", None), "temperature", None))
