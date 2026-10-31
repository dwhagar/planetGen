# planetgen/physics/interior_heat.py

"""
The heat flowing out of a rocky body: radiogenic heat plus the heat it was born with (GEN.699)
==============================================================================================

The total heat of the interior is **radiogenic heat plus primordial heat**, and the Urey ratio is what comes out, not a draw. A
body's primordial heat is the energy it kept from its making: the gravitational energy of assembling it (accretion) and of sinking
its iron into a core (differentiation), less what the magma ocean radiated away. It leaves as the mantle cools.

The model is the one-box parameterised convection of Thiriet et al. 2019 as used by Guimond, Rudge and Shorttle 2022
(arXiv:2201.05636, section 2.2; READ) for stagnant-lid planets of 0.1 to 5 Earth masses. The mantle (with the core at its
temperature, `Tc` follows `Tm` in that paper's own runs) has one temperature `Tm`:

    (Mm cm + Mc cc) dTm/dt = -qu Au + Mm h(t)

`h(t)` is the radiogenic heat per kilogram of the bulk silicate at the age `t` (`radiogenic.heat_terms`, the stored concentrations
decayed from formation, 26Al included), so the heat of the rock and the heat of the interior read the same stored content. The
heat flux out of the top is conductive across the upper boundary layer:

    qu = km dT / delta,   delta = d (Ra_crit / Ra)^(1/3),   Ra = alpha rho g dT d^3 / (kappa eta(Tm)),
    dT = a_rh R Tm^2 / Ea,   eta(Tm) = eta0 exp(Ea / R Tm)

with `d` the mantle's thickness, `Ra_crit` 450, `a_rh` 2.44 (Thiriet et al. 2019, from 3D spherical convection runs) and `dT` the
temperature drop across the upper boundary layer, the viscous temperature scale times `a_rh`. Mass enters through the gravity, the
thickness, the area and the heat capacity; composition through the heat sources, the core share and, via the radius, the water
(the dry part's radius is `bulk_structure.radius_fit_earth` or `small_body_radius_earth`, the core's `R CMF^0.5`,
Zeng and Jacobsen 2017 as Guimond et al. use them).

Two things are drawn, the rest follows (`draw_interior`): the viscosity prefactor `eta0` and activation energy `Ea` (the paper's own
uncertainty on the mantle's rheology, uniform over its ranges). The heat the body is born with is computed, not drawn: it starts at
`Ts + E / C` (`E` the gravitational binding energy of the differentiated body, `binding_energy_j`: two uniform shells, the accretion
of the whole and the sinking of the iron; `C` the heat capacity above) and never above the liquidus at the middle of the mantle
(`mantle_ceiling_k`, Manske et al. 2025: a magma ocean radiates away what is above its liquidus, as does the 26Al of the first million
years). The surface is the body's own surface temperature once its climate is known (`update_heat`), 273 K otherwise. A planet starts at the ceiling, a
planetesimal cold (its `E / C` is a kelvin), a moon-sized body in between. The model is nearly memoryless (the paper primes it for 5
Gy and finds the paths parallel), so the start matters for a young body and fades for an old one: a younger body gives more heat.

The lid (GEN.721). Above is a **stagnant lid**: one rigid plate, the convection under it working with the drop `dT` across the
boundary layer above. A **mobile lid** (plate tectonics, `lid_regime`, drawn once after the hydrosphere: `generate_lid`) carries the whole
drop from the mantle to the surface, `dT = Tm - Ts`, through plates that move, with the same `Ra` and a mobile `beta` of 1/3 (boundary-layer
theory, Gurnis 1989 as read in O'Neill 2020; the exponent and the critical Rayleigh number are stand-ins, see `MOBILE_HEAT_FLOW_EXPONENT`).
At the same mantle temperature a mobile lid loses about ten times the heat of a stagnant one, so it cools to a lower temperature (an Earth
twin: about 1,430 K against 1,950 K for the stagnant lid at 4.568 Gy, with the Earth's potential temperature near 1,600 K) and gives more
heat in the first few hundred million years at the middle of the paper's rheologies and less afterwards, when the heat it was born with is spent. Over the paper's rheologies an Earth twin
gives 25 to 32 TW and a Urey ratio of 0.66 to 0.83 with a mobile lid and 27 to 36 TW and 0.59 to 0.76 with a stagnant one, where the
Earth's measured values are 42 to 47 TW and 0.47 to 0.52 (Nimmo et al. 2020): this one box holds the radiogenic heat (about 21 TW) and
the cooling of one temperature and does not reach the measured total with either lid (the heat inside the geoneutrino range, 38.2 +13.6
-12.7 TW). The paper holds the surface at 273 K (a mild effect for its boundary layer); here a body's own surface temperature is used when known. The lid is thin
against the mantle and is not followed, so its radiogenic heat counts with the mantle's. A body whose regime has not been drawn yet (the
hydrosphere is made before it, and reads the heat of a stagnant lid) and a row from before GEN.721 are stagnant.
"""

import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import (atmosphere, bulk_structure, constants, core_heat, giant_heat, icy_tide, magnetism, moon_orbits, radiogenic,
                               spin_torques, tidal_heating, tidal_rheology, unit_registry as units)
from planetgen.util import draw

_kg_to_mearth = units.converter("kg", "Mearth")

BODY_FIELDS = ("mantle_eta0_pa_s", "mantle_activation_kj_mol")
"""tuple: What a body stores of its interior's draws, in `planets` and `moons`."""

OUTPUT_FIELDS = ("urey_ratio", "internal_heat_flux_w_m2", "tidal_flux_w_m2")
"""tuple: What a body stores of the heat that comes out of them (`generate_interior`): the Urey ratio, the heat flow of its own (radiogenic and
primordial) and its tide's (GEN.737; the habitability energy term counts both)."""

REGIME_FIELDS = ("lid_regime",)
"""tuple: What a body stores of its tectonic regime (`generate_lid`, GEN.721): `"mobile"` or `"stagnant"`, `None` for a body with no interior."""

STAGNANT, MOBILE = "stagnant", magnetism.LID_MOBILE
"""str: The two lids. A stagnant lid is one rigid plate; a mobile lid is plate tectonics."""

GAS_CONSTANT = 8.314462618
"""float: J/(mol K)."""
GY_S = 1.0e9 * constants.SECONDS_PER_YEAR
"""float: Seconds in a gigayear."""

MANTLE_DENSITY = 3500.0
"""float: kg/m3 (Thiriet et al. 2019 via Guimond et al. 2022 Table 2; READ)."""
MANTLE_SPECIFIC_HEAT = 1142.0
"""float: J/(kg K), the same table."""
CORE_SPECIFIC_HEAT = 840.0
"""float: J/(kg K), the same table."""
MANTLE_CONDUCTIVITY = 4.0
"""float: W/(m K), the same table."""
MANTLE_EXPANSIVITY = 2.5e-5
"""float: 1/K, the same table."""
MANTLE_DIFFUSIVITY = 1.0e-6
"""float: m2/s, the same table."""
CRITICAL_RAYLEIGH = 450.0
"""float: The upper boundary layer's, the same table (Thiriet et al. 2019)."""
VISCOSITY_SCALE_COEFFICIENT = 2.44
"""float: `a_rh`, the same table: the boundary layer's temperature drop in viscous temperature scales."""
HEAT_FLOW_EXPONENT = 1.0 / 3.0
"""float: `beta`, the same table."""
MOBILE_HEAT_FLOW_EXPONENT = 1.0 / 3.0
"""float: `beta` of a mobile lid: boundary-layer theory, and what Gurnis 1989 found for systems whose lid is viscous with weak zones to
let plates move (read in O'Neill 2020, arXiv:2006.05654 section 1; Moresi and Solomatov 1998 find 0.224 to 0.293 for a yielding lid, and
O'Neill's own runs 0.26). STAND-IN in the gap registry: the paper-range 1/3 is used and the mobile `Ra_crit` below is the stagnant one."""
MOBILE_CRITICAL_RAYLEIGH = CRITICAL_RAYLEIGH
"""float: STAND-IN: the mobile lid's critical Rayleigh number (not read from a source for this model; the stagnant lid's)."""
SURFACE_TEMPERATURE_K = 273.0
"""float: The surface the boundary layer is measured against, the same table."""

MAX_VISCOSITY_EXPONENT = 600.0
"""float: `Ea / R T` beyond which the mantle does not convect (a viscosity of 1e260 times the prefactor): a body with a cold surface, a
small one far from its star, starts at that temperature and warms with its radiogenic heat until convection starts; the conduction that
carries heat out of a mantle that does not convect is not modelled."""

TIDE_NEGLIGIBLE_SHARE = 1.0e-3
"""float: A tide under this share of the body's own heat is added to the output without running the box again (it moves the mantle's
temperature by under a part in a thousand of the heat)."""
TIDE_STEP_LIMIT_K = 10.0
"""float: The most a step may change the mantle's temperature when a tide heats it, K: the tide's `k2/Q` changes by a factor of several over a few tens of kelvin
near the solidus, so the step follows it (GEN.955)."""
TIDE_SAMPLES = 8
"""int: The temperatures between the mantle's and the heat pipe's at which the tide is sampled to see whether it matters at all (GEN.955)."""
STEP_LIMIT_GY = 0.05
"""float: The longest time step, Gy."""
FIRST_STEP_GY = 1.0e-4
CORE_STRIDE = 4          # the box's steps to one of the core's (GEN.838)
"""float: The first, Gy; steps grow by `STEP_GROWTH` until `STEP_LIMIT_GY` or the cooling rate asks for less."""
STEP_GROWTH = 1.15
"""float: The factor steps grow by."""


class InteriorHeatError(ValueError):
    """A body this model cannot follow."""


@dataclass(frozen=True)
class Interior:
    """
    What a body draws of its interior (`BODY_FIELDS`).

    Attributes:
        eta0_pa_s (float): The viscosity prefactor of `eta(T) = eta0 exp(Ea / R T)`, Pa s.
        activation_kj_mol (float): `Ea`, kJ/mol.
    """

    eta0_pa_s: float
    activation_kj_mol: float


@dataclass(frozen=True)
class Heat:
    """
    A body's heat at its age.

    Attributes:
        total_w (float): The heat flowing out of the interior.
        radiogenic_w (float): Of it, from the decay in the bulk silicate now.
        primordial_w (float): The rest, the cooling of the interior: the heat it was born with.
        mantle_temperature_k (float): The mantle's temperature.
        flux_w_m2 (float): `total_w` over the body's surface.
        released_j (float): The heat the interior has carried out since the body formed: the heat flow integrated over its age
            (the melt the mantle has made, and so the gas it has given off, follows it: `degassing`).
        tidal_w (float): Of `total_w`, the tide's (GEN.737, `tidal_heating`): 0 for a body on a circle or with none stored.
        core (core_heat.Core or None): The core's own temperature, boundary heat flow, inner core and dynamo (GEN.838), `None` for a body with no metal.
    """

    total_w: float
    radiogenic_w: float
    primordial_w: float
    mantle_temperature_k: float
    flux_w_m2: float
    released_j: float
    tidal_w: float = 0.0
    core: core_heat.Core | None = None

    @property
    def own_w(self):
        """The heat flowing out that is not the tide's: the radiogenic and the primordial."""
        return max(self.total_w - self.tidal_w, 0.0)

    @property
    def own_flux_w_m2(self):
        """`own_w` over the body's surface."""
        return self.flux_w_m2 * self.own_w / self.total_w if self.total_w > 0.0 else 0.0

    @property
    def urey_ratio(self):
        """The radiogenic share of the heat flowing out that is not the tide's (an output, 0 for a body with no heat)."""
        return self.radiogenic_w / self.own_w if self.own_w > 0.0 else 0.0


def draw_interior(rng=draw):
    """A body's `Interior`: eta0 and Ea uniform over the ranges of Guimond et al. 2022 (`tuning.INTERIOR_VISCOSITY_PREFACTOR_RANGE_PA_S`,
    `INTERIOR_ACTIVATION_RANGE_KJ_MOL`). Always two draws."""
    return Interior(rng.uniform(*tuning.INTERIOR_VISCOSITY_PREFACTOR_RANGE_PA_S), rng.uniform(*tuning.INTERIOR_ACTIVATION_RANGE_KJ_MOL))


def interior_of(body):
    """The `Interior` a body stores, or `None` when it stores none (a giant, or a row from before GEN.699)."""
    values = [getattr(body, name, None) for name in BODY_FIELDS]
    return None if any(value is None for value in values) else Interior(*values)


def binding_energy_j(mass_kg, radius_m, core_mass_kg, core_radius_m):
    """The energy that assembling a body of two uniform shells releases, J: a core of `core_mass_kg` inside `core_radius_m` under a
    mantle to `radius_m` (the integral of G m(r) dm / r, exact for uniform shells). With no core it is 3/5 G M^2 / R, the accretion
    of a uniform body; the difference is the energy of sinking the iron."""
    mantle_mass = mass_kg - core_mass_kg
    core = 0.6 * constants.G * core_mass_kg ** 2 / core_radius_m if core_radius_m > 0.0 else 0.0
    if mantle_mass <= 0.0:
        return core
    volume = 4.0 / 3.0 * math.pi * (radius_m ** 3 - core_radius_m ** 3)
    density = mantle_mass / volume
    a = core_mass_kg - 4.0 / 3.0 * math.pi * density * core_radius_m ** 3
    shell = 4.0 * math.pi * density * (a * (radius_m ** 2 - core_radius_m ** 2) / 2.0
                                       + 4.0 / 3.0 * math.pi * density * (radius_m ** 5 - core_radius_m ** 5) / 5.0)
    return core + constants.G * shell


def dry_radius_km(mass_kg, radius_km, rock):
    """The radius of the rock and iron under the water, ice or envelope: the body's own when it has none, else the dry planet of the
    same dry mass and core share (`bulk_structure.radius_fit_earth`, or `small_body_radius_earth` for a small one); the dry share's
    cube root of the radius where neither fit holds."""
    if rock.wet <= 0.0:
        return radius_km
    dry = rock.silicate + rock.metal
    mass_earth = _kg_to_mearth(mass_kg * dry)
    cmf = rock.metal / dry
    if 0.1 <= mass_earth <= 20.0 and cmf <= 0.8:
        return bulk_structure.radius_fit_earth(mass_earth, cmf) * constants.EARTH_RADIUS_KM
    if 1.5e-4 <= mass_earth <= 0.3:
        return bulk_structure.small_body_radius_earth(mass_earth, cmf) * constants.EARTH_RADIUS_KM
    return radius_km * dry ** (1.0 / 3.0)


def liquidus_k(pressure_gpa):
    """The mantle's liquidus at `pressure_gpa`, K: `a_l (p + p_l)^c_l` (Manske et al. 2025, `tuning.INTERIOR_LIQUIDUS_*`)."""
    return tuning.INTERIOR_LIQUIDUS_COEFFICIENT_K_GPA * (pressure_gpa + tuning.INTERIOR_LIQUIDUS_OFFSET_GPA) ** tuning.INTERIOR_LIQUIDUS_EXPONENT


def mantle_ceiling_k(gravity, thickness_m):
    """The hottest the one-box mantle is, K: the liquidus at the pressure in the middle of the mantle (`MANTLE_DENSITY g d / 2`). 3,100 K for the
    Earth, 1,900 K for a body the size of the Moon."""
    return liquidus_k(MANTLE_DENSITY * gravity * thickness_m / 2.0e9)


@dataclass(frozen=True)
class Mantle:
    """
    The body's rock as the box sees it (`mantle_of`).

    Attributes:
        dry_kg (float): The silicate and metal.
        dry_m (float): Its radius with the water taken off.
        thickness_m (float): The mantle's, `dry_m` less the core.
        gravity (float): At the dry surface, m/s2.
        ceiling_k (float): `mantle_ceiling_k`, the hottest the box lets the mantle be.
        pressure_gpa (float): In the middle of the mantle, where the melt is read.
    """

    dry_kg: float
    dry_m: float
    thickness_m: float
    gravity: float
    ceiling_k: float
    pressure_gpa: float

    @property
    def solidus_k(self):
        return tidal_heating.solidus_k(self.pressure_gpa)

    @property
    def liquidus_k(self):
        """All melt at and above it (`tidal_heating.liquidus_k`)."""
        return tidal_heating.liquidus_k(self.pressure_gpa)

    @property
    def heat_pipe_k(self):
        """The mantle temperature that holds the critical melt fraction (`tidal_heating.heat_pipe_temperature_k`)."""
        return tidal_heating.heat_pipe_temperature_k(self.pressure_gpa)


def mantle_of(mass_kg, radius_km, rock):
    """The `Mantle` of a body of `mass_kg` and `radius_km` made of `rock`; `None` for a body with no silicate."""
    silicate_kg, metal_kg = mass_kg * rock.silicate, mass_kg * rock.metal
    if silicate_kg <= 0.0:
        return None
    dry_m = dry_radius_km(mass_kg, radius_km, rock) * 1000.0
    dry_kg = silicate_kg + metal_kg
    thickness_m = dry_m - dry_m * (metal_kg / dry_kg) ** 0.5
    gravity = constants.G * dry_kg / dry_m ** 2
    pressure_gpa = MANTLE_DENSITY * gravity * thickness_m / 2.0e9
    return Mantle(dry_kg, dry_m, thickness_m, gravity, mantle_ceiling_k(gravity, thickness_m), pressure_gpa)


def _flux_w_m2(temperature_k, interior, gravity, thickness_m, lid=STAGNANT, surface_k=SURFACE_TEMPERATURE_K):
    """`qu`, W/m2, at the mantle temperature. A stagnant lid loses heat across the boundary layer's viscous temperature scale `dT = a_rh R
    Tm^2 / Ea` (the lid is cold and rigid, the convection under it only has that drop to work with); a mobile lid's plates carry the
    whole drop from the mantle to the surface, `dT = Tm - Ts`, with the viscosity of the mantle at `Tm` and a mobile `beta`."""
    activation = interior.activation_kj_mol * 1000.0
    exponent_of_viscosity = activation / (GAS_CONSTANT * temperature_k)
    if exponent_of_viscosity > MAX_VISCOSITY_EXPONENT:
        return 0.0
    viscosity = interior.eta0_pa_s * math.exp(exponent_of_viscosity)
    if lid == MOBILE:
        drop, critical, exponent = temperature_k - surface_k, MOBILE_CRITICAL_RAYLEIGH, MOBILE_HEAT_FLOW_EXPONENT
    else:
        drop, critical, exponent = (VISCOSITY_SCALE_COEFFICIENT * GAS_CONSTANT * temperature_k ** 2 / activation,
                                    CRITICAL_RAYLEIGH, HEAT_FLOW_EXPONENT)
    if drop <= 0.0:
        return 0.0      # no temperature difference to carry heat across: a mobile lid whose mantle is at its surface's temperature
    rayleigh = MANTLE_EXPANSIVITY * MANTLE_DENSITY * gravity * drop * thickness_m ** 3 / (MANTLE_DIFFUSIVITY * viscosity)
    layer = thickness_m * (critical / rayleigh) ** exponent
    return MANTLE_CONDUCTIVITY * drop / layer


def _tide_w(forcing, mantle, kelvin, interior):
    """The tide's power, W, in a body of this `mantle` (its rock, with the water taken off) with the mantle at `kelvin`: the eccentricity tide at the orbital rate
    and, for a planet, the tides of its moons against its spin (`_moons_w`), each with `tidal_rheology`'s `k2/Q` at its frequency, the melt `kelvin` makes
    between the solidus and the liquidus, and the body's drawn viscosity."""
    power = 0.0
    if forcing.eccentricity > 0.0:
        found = tidal_rheology.rock_response(mantle.dry_kg, mantle.dry_m, tidal_heating.mean_motion_rad_s(forcing), kelvin, interior.eta0_pa_s,
                                             interior.activation_kj_mol * 1000.0, mantle.solidus_k, mantle.liquidus_k)
        power += tidal_heating.peale_power_w(found.k2_over_q, forcing, mantle.dry_m)
    return power + _moons_w(forcing, mantle, kelvin, interior)


def _moons_w(forcing, mantle, kelvin, interior):
    """The power the tides of a planet's moons dissipate in it, W: for each moon `-N (spin - n)` with the torque `N` of `spin_torques` at the mantle's `k2/Q` for the
    frequency of that tide (GEN.959; the planet's spin is `forcing.spin_rate_rad_s`, the moons go round at the Kepler rates of `forcing.moons`)."""
    if not forcing.moons:
        return 0.0
    planet = spin_torques.Planet(mantle.dry_kg, mantle.dry_m)
    torques = spin_torques.Torques(planet, (), spin_torques.rheology_dissipation(planet, kelvin, interior.eta0_pa_s, interior.activation_kj_mol * 1000.0,
                                                                                mantle.solidus_k, mantle.liquidus_k),
                                   None, spin_torques.GRAVITATIONAL_TORQUE_FACTOR, moon_orbits.COROTATION_SOFTNESS)
    power = 0.0
    for mass_kg, distance_m, retrograde, rate in forcing.moons:
        power += max(-torques.gravitational(forcing.spin_rate_rad_s, spin_torques.Raiser(mass_kg, distance_m, rate)) * (forcing.spin_rate_rad_s - rate), 0.0)
    return power


def _under_a_magma_ocean(terms, silicate_kg, area_m2, age_gy, surface_k, tide_w):
    """
    The `Heat` of a body whose surface is at or above its mantle's liquidus (`mantle_ceiling_k`): a magma ocean under a hot air. The
    mantle can be no cooler than the surface it lies under and cannot be hotter than its ceiling, so it sits at the surface's temperature
    with no drop to carry heat across; the melt carries out what the body makes, the decay of its silicate and a tide that flexes a
    molten body (`tide_w`, the tide's power at the mantle's temperature, `None` for no tide; a molten mantle dissipates little), and nothing is
    left over from the heat it was born with. The heat released over its age is what it made (GEN.856: the box's flux laws have no
    meaning without a drop, and ran to a zero divisor on a mobile lid and to a mantle cooler than its surface on a stagnant one).
    """
    now_w = silicate_kg * math.fsum(a * math.exp(-rate * age_gy) for a, rate in terms)
    made_w_gy = silicate_kg * math.fsum(a * (age_gy if rate == 0.0 else (1.0 - math.exp(-rate * age_gy)) / rate) for a, rate in terms)
    molten_w = tide_w(surface_k) if tide_w is not None else 0.0
    total_w = now_w + molten_w
    return Heat(total_w, now_w, 0.0, surface_k, total_w / area_m2, made_w_gy * GY_S + molten_w * age_gy * GY_S, molten_w)


def heat(terms, interior, mass_kg, radius_km, age_gy, rock, forcing=None, lid=STAGNANT, surface_k=SURFACE_TEMPERATURE_K):
    """
    The interior's heat at `age_gy`: the one-box model above integrated from the formation of the body.

    Args:
        terms (tuple): `radiogenic.heat_terms(content, rock.aluminium)`: W/kg of bulk silicate and decay constants per Gy.
        interior (Interior): The body's draws.
        mass_kg, radius_km (float): The body's.
        age_gy (float): Its age.
        rock (radiogenic.Rock): What it is made of.
        forcing (tidal_heating.Forcing, optional): The orbit that flexes it (GEN.737). Its tide is a second heat source of the box, the Peale
            power with the `k2/Q` of the mantle at its temperature (`tidal_rheology`, GEN.955), integrated with the cooling. If that tide is at least
            the body's own heat the mantle is never hotter than `Mantle.heat_pipe_k`, where its melt erupts and carries the heat out; if it makes more
            heat there than it conducts the body is melt-limited and holds that temperature (Io).
        lid (str): `STAGNANT` or `MOBILE`, the regime that sets the heat flux law (GEN.721).
        surface_k (float): The surface temperature the mantle's drop is measured against (a mobile lid) and the body starts from (GEN.735).

    Returns:
        Heat: Zero heat for a body with no silicate or no core-and-mantle to convect.
    """
    if age_gy < 0.0:
        raise InteriorHeatError(f"an age of {age_gy} Gy is before formation")
    area_m2 = 4.0 * math.pi * (radius_km * 1000.0) ** 2
    silicate_kg, metal_kg = mass_kg * rock.silicate, mass_kg * rock.metal
    if silicate_kg <= 0.0:
        return Heat(0.0, 0.0, 0.0, SURFACE_TEMPERATURE_K, 0.0, 0.0)
    mantle = mantle_of(mass_kg, radius_km, rock)
    dry_m, dry_kg, thickness_m, gravity, cap_k = mantle.dry_m, mantle.dry_kg, mantle.thickness_m, mantle.gravity, mantle.ceiling_k
    core_m = dry_m * (metal_kg / dry_kg) ** 0.5
    area_dry_m2 = 4.0 * math.pi * dry_m ** 2
    capacity = silicate_kg * MANTLE_SPECIFIC_HEAT + metal_kg * CORE_SPECIFIC_HEAT

    def radiogenic_w(t_gy):
        return silicate_kg * math.fsum(a * math.exp(-rate * t_gy) for a, rate in terms)

    energy_j = binding_energy_j(dry_kg, dry_m, metal_kg, core_m)
    activation = interior.activation_kj_mol * 1000.0

    def tide_w(kelvin):
        """The tide's power with the mantle at `kelvin`: the Peale formula with the rheology's `k2/Q` there (GEN.955)."""
        return _tide_w(forcing, mantle, kelvin, interior)

    if surface_k >= cap_k:
        return _under_a_magma_ocean(terms, silicate_kg, area_m2, age_gy, surface_k, tide_w if forcing is not None else None)

    def conduction_w(kelvin):
        return _flux_w_m2(kelvin, interior, gravity, thickness_m, lid, surface_k) * area_dry_m2

    def cool(tide=None, ceiling_k=cap_k):
        """
        The temperature and released heat at `age_gy`, with `tide` (the tide's power at the mantle's temperature, W) added to the radiogenic. The mantle is
        never hotter than `ceiling_k`: above `cap_k` the melt is radiated away and is not counted, but under it (the heat pipe, `Mantle.heat_pipe_k`) the
        melt erupts and the heat it carries out is released heat. Returns the temperature, the heat released and whether the last step ended on the
        ceiling with more heat made than it could carry (the melt-limited state).
        """
        top_k = min(ceiling_k, cap_k)

        def clipped(kelvin):
            # a strong tide moves the temperature past the cap within one step, and a stage beyond the range of the viscosity law makes no sense
            return min(max(kelvin, surface_k), top_k)

        def rate(radiogenic, kelvin, conduction=None):
            if conduction is None:
                conduction = conduction_w(clipped(kelvin))
            return GY_S * (radiogenic + (tide(clipped(kelvin)) if tide is not None else 0.0) - conduction) / capacity

        start_k = surface_k + energy_j / capacity
        temperature = min(start_k, top_k)
        t, step, released_j, held = 0.0, FIRST_STEP_GY, 0.0, False
        if top_k < cap_k and start_k > top_k:
            released_j += capacity * (start_k - top_k)
        core = core_heat.Stepper(interior, metal_kg, silicate_kg, core_m, dry_m) if metal_kg > 0.0 else None
        if core:
            core.start(temperature)
        core_steps = 0
        # the radiogenic power at the start of a step is that at the end of the one before, and the two middle stages share one time;
        # the flux at the start is the one the last step ended with (PERF.201: the same values, each worked out once)
        radiogenic_start = radiogenic_w(0.0)
        flux = _flux_w_m2(temperature, interior, gravity, thickness_m, lid, surface_k)
        while t < age_gy:
            # the cooling rate's own time scale: the flux changes by Ea / (R T^2) per kelvin (an overestimate by about three)
            stiffness = GY_S * flux * area_dry_m2 / capacity * activation / (GAS_CONSTANT * temperature ** 2)
            h = min(step, age_gy - t, 1.0 / stiffness if stiffness > 0.0 else step)
            k1 = rate(radiogenic_start, temperature, flux * area_dry_m2 if surface_k <= temperature <= top_k else None)
            if tide is not None and k1 != 0.0 and not (temperature >= top_k and k1 > 0.0):   # a mantle held on its ceiling goes no further
                h = min(h, TIDE_STEP_LIMIT_K / abs(k1))
            radiogenic_middle = radiogenic_w(t + 0.5 * h)
            radiogenic_end = radiogenic_w(t + h)
            k2 = rate(radiogenic_middle, temperature + 0.5 * h * k1)
            k3 = rate(radiogenic_middle, temperature + 0.5 * h * k2)
            k4 = rate(radiogenic_end, temperature + h * k3)
            # a step can overshoot the cold end as well as the cap: below the surface the viscosity law runs on a negative or tiny exponent
            # and the flux (and so the stiffness) is astronomical, which collapses the step to nothing (GEN.866)
            before = temperature
            raw_k = temperature + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
            temperature = clipped(raw_k)
            held = raw_k > top_k
            if held and top_k < cap_k:
                released_j += capacity * (raw_k - top_k)
            if core:
                # the core's own time scale is Gy, so it takes one step over a few of the box's (the midpoint rule over their sum)
                if not core_steps:
                    core_from_t, core_from_k, core_h = t, before, 0.0
                core_steps += 1
                core_h += h
                if core_steps == CORE_STRIDE or t + h >= age_gy:
                    core.advance(core_from_t, core_h, core_from_k, temperature)
                    core_steps = 0
            flux_end = _flux_w_m2(temperature, interior, gravity, thickness_m, lid, surface_k)
            released_j += 0.5 * (flux + flux_end) * area_dry_m2 * h * GY_S
            flux = flux_end
            radiogenic_start = radiogenic_end
            t += h
            step = min(step * STEP_GROWTH, STEP_LIMIT_GY)
        return temperature, released_j, held, core.state(temperature) if core else None

    now_w = radiogenic_w(age_gy)

    def result(temperature, released_j, tidal_w, unrun_tidal_w=0.0, core=None):
        """The `Heat` at `temperature`; `unrun_tidal_w` is a tide the box did not run with, flowing out beside the conduction."""
        total_w = conduction_w(temperature) + unrun_tidal_w
        if temperature >= cap_k:
            total_w = max(total_w, now_w + tidal_w)
        return Heat(total_w, now_w, total_w - now_w - tidal_w, temperature, total_w / area_m2, released_j, tidal_w, core)

    temperature, released_j, _, core = cool()
    if forcing is None:
        return result(temperature, released_j, 0.0, core=core)
    own_w = conduction_w(temperature)
    pipe_k = min(mantle.heat_pipe_k, cap_k)
    reach_k = max(pipe_k, temperature)
    biggest_w = max(tide_w(temperature + (reach_k - temperature) * i / TIDE_SAMPLES) for i in range(TIDE_SAMPLES + 1))
    if biggest_w <= TIDE_NEGLIGIBLE_SHARE * max(own_w, now_w):
        unrun_w = tide_w(temperature)
        return result(temperature, released_j, unrun_w, unrun_tidal_w=unrun_w, core=core)
    # a tide that at some temperature up to the heat-pipe temperature is at least the body's own heat can melt the mantle to the critical melt
    # fraction, where it erupts and carries the heat out: the mantle is never hotter than that, and holds it if it makes more heat there than it
    # can conduct (melt-limited; the heat pipe of Io). A weaker tide is a small heat source of the box.
    warm, warm_released_j, held, warm_core = cool(tide_w, pipe_k if biggest_w >= max(own_w, now_w) else cap_k)
    if held:
        tidal_w = tide_w(warm)
        total_w = max(now_w + tidal_w, conduction_w(warm))
        released_j = max(warm_released_j, total_w * age_gy * GY_S)
        return Heat(total_w, now_w, total_w - now_w - tidal_w, warm, total_w / area_m2, released_j, tidal_w, warm_core)
    return result(warm, warm_released_j, tide_w(warm), core=warm_core)


def heat_flow(content, interior, mass_kg, radius_km, age_gy, rock, forcing=None, lid=STAGNANT, surface_k=SURFACE_TEMPERATURE_K):
    """`heat` for a body's stored `radiogenic.Content`: the decay terms of its concentrations and aluminium."""
    return heat(radiogenic.heat_terms(content, rock.aluminium), interior, mass_kg, radius_km, age_gy, rock, forcing, lid, surface_k)


def tide_temperature_k(body):
    """
    The mantle temperature that sets a body's tidal `k2` and `Q`, K (GEN.955): the box run to the star's age with no tide, a stagnant lid and a surface of
    `SURFACE_TEMPERATURE_K`, so it is a function of the body's stored content, mass, radius and age alone and the spin, the circularisation and the heat
    agree whenever each is worked out. The tide warms the real mantle (the box in `heat` runs with it); this is the temperature it starts from. No draws.
    """
    content, interior = radiogenic.content_of(body), interior_of(body)
    age_gy = getattr(body.star, "age", None) or tuning.ROGUE_EARTH_REFERENCE_AGE_GY
    return heat_flow(content, interior, body.mass, body.radius, age_gy, radiogenic.rock_of(body), None, STAGNANT,
                     SURFACE_TEMPERATURE_K).mantle_temperature_k


def tidal_response(body, frequency_rad_s):
    """
    `(k2, Q)` of a body, its tidal Love number and quality factor, at a tide of `frequency_rad_s` (GEN.955): a gas giant's fixed pair
    (`tidal_heating`), an icy body's ice shell's (`icy_tide.k2_over_q` at the heat flux the shell conducts, GEN.960; `inf` for a body with no shell), or for rock `tidal_rheology.rock_response` of its mass and radius at `tide_temperature_k` with the viscosity it drew, the melt that
    temperature makes in its mantle, and the tide's frequency. A cold small body is stiff (Q up to `tidal_rheology.COLD_Q_CEILING`), a mantle near its
    solidus dissipates most, a molten one is weak.
    """
    if getattr(body, "body_type", None) == "g":
        return tidal_heating.GIANT_K2, tidal_heating.GIANT_TIDAL_Q
    if tidal_heating.is_icy(body):
        radius_m = body.radius * 1000.0
        k2 = icy_tide.love_k2(body.mass, radius_m, frequency_rad_s)
        k2_over_q = icy_tide.k2_over_q(body.mass, radius_m, frequency_rad_s, icy_shell_flux_w_m2(body), icy_shell_top_k(body))
        return k2, (k2 / k2_over_q if k2_over_q > 0.0 else math.inf)
    interior = interior_of(body)
    if interior is None or radiogenic.content_of(body) is None:
        raise InteriorHeatError("a rocky body's tidal response needs its stored content and mantle draws")
    mantle = mantle_of(body.mass, body.radius, radiogenic.rock_of(body))
    if mantle is None:
        raise InteriorHeatError("a body with no silicate has no mantle to flex")
    found = tidal_rheology.rock_response(body.mass, body.radius * 1000.0, frequency_rad_s, tide_temperature_k(body), interior.eta0_pa_s,
                                         interior.activation_kj_mol * 1000.0, mantle.solidus_k, mantle.liquidus_k)
    return found.k2, found.q


def generate_interior(body):
    """
    Sets `BODY_FIELDS` on a rocky planet or moon: the two draws of its mantle's rheology (`draw_interior`), made once with the rest of
    its rock (`radiogenic.generate_content` first). A body with no content (a gas giant) has every field `None`, outputs included.
    """
    if radiogenic.content_of(body) is None:
        for name in (*BODY_FIELDS, *OUTPUT_FIELDS, *REGIME_FIELDS):
            setattr(body, name, None)
        return
    interior = draw_interior()
    body.mantle_eta0_pa_s, body.mantle_activation_kj_mol = interior.eta0_pa_s, interior.activation_kj_mol


def released_by(body, age_gy):
    """
    The heat the interior of `body` has carried out by `age_gy` (J): the thermal history of `update_heat` run to that age instead of the
    star's. 0 for a body with no stored content or interior. No draws.
    """
    content, interior = radiogenic.content_of(body), interior_of(body)
    if content is None or interior is None:
        return 0.0
    return heat_flow(content, interior, body.mass, body.radius, age_gy, radiogenic.rock_of(body), None, lid_of(body), surface_of(body)).released_j


def lid_of(body):
    """The lid a body stores (`STAGNANT` or `MOBILE`); `STAGNANT` for one that has none yet (its regime is drawn after its hydrosphere,
    so the heat the hydrosphere reads is a stagnant lid's) or a row from before GEN.721."""
    return MOBILE if getattr(body, "lid_regime", None) == MOBILE else STAGNANT


def icy_shell_top_k(body):
    """The temperature at the top of an icy body's shell, K: its surface temperature once its climate is known, else `icy_tide.TOP_FALLBACK_K`."""
    surface_k = getattr(body, "surface_temperature", None)
    return surface_k if surface_k and surface_k > 0.0 else icy_tide.TOP_FALLBACK_K


def icy_shell_flux_w_m2(body):
    """The heat flux an icy body's shell conducts, W/m2: its own heat (stored once `update_heat` has run, else worked out here with no tide) plus the tide's."""
    own = getattr(body, "internal_heat_flux_w_m2", None)
    if own is None:
        content, interior = radiogenic.content_of(body), interior_of(body)
        own = heat_flow(content, interior, body.mass, body.radius, getattr(body.star, "age", None) or tuning.ROGUE_EARTH_REFERENCE_AGE_GY,
                        radiogenic.rock_of(body), None, lid_of(body), surface_of(body)).own_flux_w_m2 if content is not None and interior is not None else icy_tide.OWN_FLUX_FALLBACK_W_M2
    return max(own + (getattr(body, "tidal_flux_w_m2", None) or 0.0), icy_tide.MIN_FLUX_W_M2)


def surface_of(body):
    """The surface temperature a body's interior is measured against, K: its own once its climate is known, else `SURFACE_TEMPERATURE_K`."""
    surface_k = getattr(body, "surface_temperature", None)
    return surface_k if surface_k and surface_k > 0.0 else SURFACE_TEMPERATURE_K


def generate_lid(body):
    """
    Sets `lid_regime` on a rocky planet or moon: one draw, `MOBILE` when it falls under `magnetism.mobile_lid_chance` (a liquid surface
    ocean to lubricate the plates, and a mass between 0.5 and 5 Earth masses at 0.6, heavier at 0.3: the literature is unsettled), else
    `STAGNANT`. Call after the hydrosphere. A body with no stored content (a gas giant) has `None` and draws nothing. The magnetic field
    and the ground dose read the stored regime, so the three agree (they each drew their own coin before GEN.721).
    """
    if radiogenic.content_of(body) is None:
        body.lid_regime = None
        return
    mass_earth = _kg_to_mearth(body.mass)
    mobile = draw.random() < magnetism.mobile_lid_chance(body, mass_earth)
    body.lid_regime = MOBILE if mobile else STAGNANT


def heat_of(body):
    """
    The `Heat` of a rocky body at its star's age now, with the tide of the orbit stored with it (`tidal_heating.forcing_of`): a rocky body's
    goes through the box (`heat`, which may melt it), an icy body's (a water mass fraction of `hydrosphere.ICY_MIN_WATER_FRACTION` or more)
    heats its ice shell (`icy_tide`, GEN.960) and is not in the box. `None` for a body with no stored content or interior (a gas giant, or a row
    from before GEN.699). No draws.
    """
    content, interior = radiogenic.content_of(body), interior_of(body)
    if content is None or interior is None:
        return None
    age_gy = getattr(body.star, "age", None) or tuning.ROGUE_EARTH_REFERENCE_AGE_GY
    forcing = None if tidal_heating.is_icy(body) else tidal_heating.forcing_of(body)
    return heat_flow(content, interior, body.mass, body.radius, age_gy, radiogenic.rock_of(body), forcing, lid_of(body), surface_of(body))


def core_of(body):
    """The core's own temperature, boundary heat flow, inner core and dynamo (`core_heat.Core`, GEN.838) at its star's age now; `None` for a
    body with no interior or no metal. No draws."""
    found = heat_of(body)
    return None if found is None else found.core


def update_heat(body):
    """
    Runs the thermal history of a rocky body to its star's age now and sets `OUTPUT_FIELDS` (the Urey ratio, the heat flow of its own and the
    tide's, GEN.737); no draws, so it can be repeated when the star's age is adjusted for its planets. The tide is that of the orbit stored
    with the body (`tidal_heating.forcing_of`): a rocky body's goes through the box (`heat`, which may melt it), an icy body's (a water mass
    fraction of `hydrosphere.ICY_MIN_WATER_FRACTION` or more) heats its ice shell (`icy_tide`, GEN.960). Returns the `Heat`, or `None`
    (fields `None`) for a body with no stored content or interior (a row from before GEN.699). A gas giant has no content: it returns `None`,
    and its `internal_heat_flux_w_m2` is its own radiated heat with the helium rain its stored layers say it has (`giant_heat.of`, GEN.810).
    """
    if radiogenic.content_of(body) is None:
        giant = giant_heat.of(body)         # GEN.810: a giant's own heat, with the helium rain its stored layers say it has
        for name in OUTPUT_FIELDS:
            setattr(body, name, None)
        if giant is not None:
            body.internal_heat_flux_w_m2 = giant.flux_w_m2
        return None
    found = heat_of(body)
    if found is None:
        for name in OUTPUT_FIELDS:
            setattr(body, name, None)
        return None
    body.urey_ratio, body.internal_heat_flux_w_m2 = found.urey_ratio, found.own_flux_w_m2
    if tidal_heating.is_icy(body):
        # GEN.960: the shell's tide in equilibrium with the shell it makes (`icy_tide.flux_from_tide`), over its surface
        forcing = tidal_heating.forcing_of(body)
        radius_m = body.radius * 1000.0
        if forcing is not None and forcing.eccentricity > 0.0:
            _, power = icy_tide.flux_from_tide(body.mass, radius_m, forcing.host_mass_kg, forcing.semi_major_axis_m, forcing.eccentricity,
                                               found.own_flux_w_m2, icy_shell_top_k(body))
            body.tidal_flux_w_m2 = power / (4.0 * math.pi * radius_m ** 2)
        else:
            body.tidal_flux_w_m2 = 0.0
    else:
        body.tidal_flux_w_m2 = found.tidal_w / (4.0 * math.pi * (body.radius * 1000.0) ** 2)
    return found
