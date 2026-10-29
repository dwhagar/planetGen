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
(the dry part's radius is `bulk_structure.radius_fit_earth` or `small_body_radius_earth`, the core's `R CMF^0.5`, Zeng and
Jacobsen 2017 as Guimond et al. use them).

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
-12.7 TW, Agostini et al. 2019). The paper holds the surface at 273 K (a mild effect for its boundary layer); here a body's own surface temperature is used when known. The lid is thin
against the mantle and is not followed, so its radiogenic heat counts with the mantle's. A body whose regime has not been drawn yet (the
hydrosphere is made before it, and reads the heat of a stagnant lid) and a row from before GEN.721 are stagnant.
"""

import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import atmosphere, bulk_structure, constants, magnetism, radiogenic, tidal_heating
from planetgen.util import draw

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
"""float: `beta`, the same table (Solomatov 1995)."""
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
STEP_LIMIT_GY = 0.05
"""float: The longest time step, Gy."""
FIRST_STEP_GY = 1.0e-4
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
    """

    total_w: float
    radiogenic_w: float
    primordial_w: float
    mantle_temperature_k: float
    flux_w_m2: float
    released_j: float
    tidal_w: float = 0.0

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
    mass_earth = mass_kg * dry / constants.EARTH_MASS_TO_KG
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
    rayleigh = MANTLE_EXPANSIVITY * MANTLE_DENSITY * gravity * drop * thickness_m ** 3 / (MANTLE_DIFFUSIVITY * viscosity)
    layer = thickness_m * (critical / rayleigh) ** exponent
    return MANTLE_CONDUCTIVITY * drop / layer


def heat(terms, interior, mass_kg, radius_km, age_gy, rock, forcing=None, lid=STAGNANT, surface_k=SURFACE_TEMPERATURE_K):
    """
    The interior's heat at `age_gy`: the one-box model above integrated from the formation of the body.

    Args:
        terms (tuple): `radiogenic.heat_terms(content, rock.aluminium)`: W/kg of bulk silicate and decay constants per Gy.
        interior (Interior): The body's draws.
        mass_kg, radius_km (float): The body's.
        age_gy (float): Its age.
        rock (radiogenic.Rock): What it is made of.
        forcing (tidal_heating.Forcing, optional): The orbit that flexes it (GEN.737). A solid body's tide is a second heat source of the
            box, constant over its age. If that tide is at least the body's own heat and warms the mantle to `tidal_heating.HEAT_PIPE_TEMPERATURE_K`,
            the body is melt-limited: the mantle holds that temperature and the partly molten body is flexed `MOLTEN_K2_OVER_Q` (Io).
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
    dry_m = dry_radius_km(mass_kg, radius_km, rock) * 1000.0
    dry_kg = silicate_kg + metal_kg
    core_m = dry_m * (metal_kg / dry_kg) ** 0.5
    thickness_m = dry_m - core_m
    gravity = constants.G * dry_kg / dry_m ** 2
    area_dry_m2 = 4.0 * math.pi * dry_m ** 2
    capacity = silicate_kg * MANTLE_SPECIFIC_HEAT + metal_kg * CORE_SPECIFIC_HEAT

    def radiogenic_w(t_gy):
        return silicate_kg * math.fsum(a * math.exp(-rate * t_gy) for a, rate in terms)

    energy_j = binding_energy_j(dry_kg, dry_m, metal_kg, core_m)
    activation = interior.activation_kj_mol * 1000.0
    cap_k = mantle_ceiling_k(gravity, thickness_m)

    def conduction_w(kelvin):
        return _flux_w_m2(kelvin, interior, gravity, thickness_m, lid, surface_k) * area_dry_m2

    def cool(source_w):
        """The temperature and released heat at `age_gy` with `source_w` of constant heat added to the radiogenic."""
        temperature = min(surface_k + energy_j / capacity, cap_k)

        def rate(t_gy, kelvin):
            # a strong tide moves the temperature past the cap within one step, and a stage beyond the range of the viscosity law makes no sense
            kelvin = min(max(kelvin, surface_k), cap_k)
            return GY_S * (radiogenic_w(t_gy) + source_w - conduction_w(kelvin)) / capacity

        t, step, released_j = 0.0, FIRST_STEP_GY, 0.0
        while t < age_gy:
            flux = _flux_w_m2(temperature, interior, gravity, thickness_m, lid, surface_k)
            # the cooling rate's own time scale: the flux changes by Ea / (R T^2) per kelvin (an overestimate by about three)
            stiffness = GY_S * flux * area_dry_m2 / capacity * activation / (GAS_CONSTANT * temperature ** 2)
            h = min(step, age_gy - t, 1.0 / stiffness if stiffness > 0.0 else step)
            k1 = rate(t, temperature)
            k2 = rate(t + 0.5 * h, temperature + 0.5 * h * k1)
            k3 = rate(t + 0.5 * h, temperature + 0.5 * h * k2)
            k4 = rate(t + h, temperature + h * k3)
            # a step can overshoot the cold end as well as the cap: below the surface the viscosity law runs on a negative or tiny exponent
            # and the flux (and so the stiffness) is astronomical, which collapses the step to nothing (GEN.866)
            temperature = min(max(temperature + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0, surface_k), cap_k)
            released_j += 0.5 * (flux + _flux_w_m2(temperature, interior, gravity, thickness_m, lid, surface_k)) * area_dry_m2 * h * GY_S
            t += h
            step = min(step * STEP_GROWTH, STEP_LIMIT_GY)
        return temperature, released_j

    now_w = radiogenic_w(age_gy)

    def result(temperature, released_j, tidal_w, unrun_tidal_w=0.0):
        """The `Heat` at `temperature`; `unrun_tidal_w` is a tide the box did not run with, flowing out beside the conduction."""
        total_w = conduction_w(temperature) + unrun_tidal_w
        if temperature >= cap_k:
            total_w = max(total_w, now_w + tidal_w)
        return Heat(total_w, now_w, total_w - now_w - tidal_w, temperature, total_w / area_m2, released_j, tidal_w)

    temperature, released_j = cool(0.0)
    if forcing is None:
        return result(temperature, released_j, 0.0)
    own_w = conduction_w(temperature)
    solid_w = tidal_heating.solid_power_w(forcing, dry_kg, dry_m)
    if solid_w <= TIDE_NEGLIGIBLE_SHARE * max(own_w, now_w):
        return result(temperature, released_j, solid_w, unrun_tidal_w=solid_w)
    warm, warm_released_j = cool(solid_w)
    if solid_w < own_w or warm < tidal_heating.HEAT_PIPE_TEMPERATURE_K:
        return result(warm, warm_released_j, solid_w)
    # melt-limited: the mantle holds the heat-pipe temperature, the melt carries out the tide's heat and the rest
    molten_w = tidal_heating.molten_power_w(forcing, dry_m)
    pipe_k = tidal_heating.HEAT_PIPE_TEMPERATURE_K
    total_w = max(now_w + molten_w, conduction_w(pipe_k))
    released_j = max(warm_released_j, total_w * age_gy * GY_S)
    return Heat(total_w, now_w, total_w - now_w - molten_w, pipe_k, total_w / area_m2, released_j, molten_w)


def heat_flow(content, interior, mass_kg, radius_km, age_gy, rock, forcing=None, lid=STAGNANT, surface_k=SURFACE_TEMPERATURE_K):
    """`heat` for a body's stored `radiogenic.Content`: the decay terms of its concentrations and aluminium."""
    return heat(radiogenic.heat_terms(content, rock.aluminium), interior, mass_kg, radius_km, age_gy, rock, forcing, lid, surface_k)


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
    mass_earth = body.mass / constants.EARTH_MASS_TO_KG
    mobile = draw.random() < magnetism.mobile_lid_chance(body, mass_earth)
    body.lid_regime = MOBILE if mobile else STAGNANT


def update_heat(body):
    """
    Runs the thermal history of a rocky body to its star's age now and sets `OUTPUT_FIELDS` (the Urey ratio, the heat flow of its own and the
    tide's, GEN.737); no draws, so it can be repeated when the star's age is adjusted for its planets. The tide is that of the orbit stored
    with the body (`tidal_heating.forcing_of`): a rocky body's goes through the box (`heat`, which may melt it), an icy body's (a water mass
    fraction of `hydrosphere.ICY_MIN_WATER_FRACTION` or more) heats its ice shell and is `tidal_heating.icy_flux_w_m2`. Returns the `Heat`, or `None`
    (fields `None`) for a body with no stored content or interior (a gas giant, or a row from before GEN.699).
    """
    content, interior = radiogenic.content_of(body), interior_of(body)
    if content is None or interior is None:
        for name in OUTPUT_FIELDS:
            setattr(body, name, None)
        return None
    age_gy = getattr(body.star, "age", None) or tuning.ROGUE_EARTH_REFERENCE_AGE_GY
    forcing = tidal_heating.forcing_of(body)
    icy = tidal_heating.is_icy(body)
    found = heat_flow(content, interior, body.mass, body.radius, age_gy, radiogenic.rock_of(body), None if icy else forcing,
                      lid_of(body), surface_of(body))
    body.urey_ratio, body.internal_heat_flux_w_m2 = found.urey_ratio, found.own_flux_w_m2
    if icy:
        body.tidal_flux_w_m2 = tidal_heating.icy_flux_w_m2(forcing, body.radius) if forcing else None
    else:
        body.tidal_flux_w_m2 = found.tidal_w / (4.0 * math.pi * (body.radius * 1000.0) ** 2) if found.tidal_w > 0.0 else None
    return found
