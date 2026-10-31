# planetgen/physics/hydrosphere.py

"""
Water on planets, moons and rogue planets (GEN.88).

The design is docs/design/habitability/activity-magnetism-radiation-hydrosphere.md
section 5. Both star-system bodies (`generate_hydrosphere`) and rogue
planets (`physics.rogue_surface`) use the functions here.

1. **Inventory.** A rocky body's water mass fraction is drawn with its bulk
   composition (`bulk_draw`, GEN.428) from the delivery law of its mass and where it formed
   (`delivered_water`, GEN.849; Earth's oceans are 2.34e-4). The mantle holds part of it: the surface
   keeps `w / b`, `b` log-uniform 1 to 8, but the mantle never holds more
   than it can store (`mantle_max_water_fraction`: 2 Earth oceans at 0.5 Earth masses to 8 at 5, Guimond et al. 2023, GEN.1147). Gas giants have no
   hydrosphere (every field `None`).
2. **Land.** A two-level hypsometry (GEN.1147, Cowan and Abbot 2014's idea): the surface is ocean crust and continental crust, each a
   Gaussian of height (`HYPSOMETRY`: shares, mean height and spread in km), the relief scaled as `g_E / g`. `D`, the global water layer, fills the
   lowest ground to a level `z`; the share of the surface under `z` is `ocean_fraction`, liquid or frozen, and `land_fraction` the rest. Earth's
   2.7 km layer covers 0.71.
3. **State.** Under water's triple-point pressure and warmer than
   `VACUUM_ICE_MAX_K` the water has sublimated away (`dry`); above the
   boiling curve it is in the air (`vapour`). Below freezing, an ice lid grows to
   `ice_shell_thickness_km` from the body's internal heat; it seals the
   whole column (`ice`) or floats on an ocean (`ice-covered ocean`). The
   lid is ice (917 kg/m3) and is compared with the column by mass; its base
   melts at 273.15 K less 0.0074 K per bar of ice above it. A body whose
   surface water is under `BRINE_MAX_WATER_FRACTION` has its salts
   concentrated: the base melts at the drawn eutectic of a perchlorate brine
   (`clathrate_brine.draw_brine`, 198 to 236 K, GEN.762) and the brine ocean
   keeps that salt's water activity. A liquid surface whose ocean class is the
   chloride brine is a `brine seep`, not a `surface ocean` (GEN.908): so little
   water that its salts fix its chemistry, often under no air at all.
4. **The compressed column and high-pressure ice (GEN.461).** The water over
   the covered area is integrated from the surface down with its density
   rising with the pressure (`water_layer.column`: 746 km deep at 11 GPa for
   a tenth of an Earth mass of water, not the 1047 km of a fixed 1000 kg/m3),
   and the pressure at its bottom is stored (`water_base_pressure_pa`). No
   liquid lies deeper than the liquidus of the bottom temperature (the melting
   curves of ice III, V and VI, recomputed from the Journaux splines); the
   rest of the column is ice VI or VII (`hp_ice_km`), and the ocean is
   ice-sealed: no contact with rock. The depths are the real, compressed
   kilometres, not water-equivalent.
5. **Ocean class** (section 5.4), by rule, first match: ice-sealed;
   chloride brine (surface water under 2e-5 of the mass); acid sulfate
   (oxidized mantle and SO2 over 1 percent of SO2 + CO2); soda (reduced
   or intermediate mantle, CO2 over 10 kPa, ocean over 90 percent of the
   surface); neutral. Each draws its pH and water activity; phosphorus is
   `high` in soda, `starved` when ice-sealed, else `limited`.
6. **Hycean** (section 5.5): a liquid surface ocean under air that is
   mostly hydrogen.

The ranges are design defaults, not fits.
"""

import math

from planetgen import tuning
from planetgen.physics import aquifers, atmosphere, clathrate_brine, constants, gases, runaway_history, saturation, water_layer
from planetgen.util import draw
from planetgen.util.random import log_uniform

HYDROSPHERE_FIELDS = (
    "water_mass_fraction", "hydrosphere", "ocean_fraction", "land_fraction",
    "ocean_depth_km", "ice_shell_km", "hp_ice_km", "water_base_pressure_pa",
    "ocean_class", "ocean_ph", "water_activity", "phosphorus", *aquifers.AQUIFER_FIELDS,
)
"""tuple: What `generate_hydrosphere` sets on a planet or moon, as stored in
`planets` and `moons` (the aquifer layer on a world with no liquid at its surface)."""

BRINE_SEEP = "brine seep"
"""str: The state of a liquid surface whose ocean is a chloride brine (GEN.908): a puddle of concentrated salt water that melts below
273 K under thin or no air, not a sea. It is a liquid surface (the ocean columns apply, and it holds water for a base), but the checks
that ask for a `surface ocean` (a sunlit sea for life, the plate tectonics of a wet world, an induction ocean) do not take it."""
HYDROSPHERE_STATES = ("dry", "vapour", "ice", "ice-covered ocean", "surface ocean", BRINE_SEEP, "hycean")
OCEAN_CLASSES = ("ice-sealed", "chloride brine", "acid sulfate", "soda", "neutral")

WATER_MELTING_K = 273.15
FROZEN_COVER_STATES = ("ice", "ice-covered ocean")
"""tuple: The hydrosphere states whose surface is ice (GEN.833: their cover is the water layer the body holds, `frozen_cover`)."""
LIQUID_WATER_DENSITY_KG_M3 = 1000.0
ICE_DENSITY_KG_M3 = water_layer.ICE_LID_DENSITY_KG_M3
ICE_CONDUCTIVITY_A_W_M = 567.0
"""float: Water ice conducts heat as k = A / T (Turcotte and Schubert)."""
VACUUM_ICE_MAX_K = 220.0
"""float: Under water's triple-point pressure, exposed ice warmer than this
sublimates away over billions of years (the Moon keeps ice only in cold
traps; Mars's caps survive at 150 to 210 K): such a body is `dry`. A
planetGen default."""
MELTING_SLOPE_K_PER_BAR = 0.0074
"""float: Ice Ih's melting point falls this much per bar."""

EARTH_OCEAN_FRACTION = 2.34e-4
MANTLE_CAP_OCEANS = ((0.5, 2.0), (5.0, 8.0))
"""tuple: The most water a mantle stores, in Earth oceans, at two masses in Earth masses (`mantle_max_water_fraction` joins them as a power law and holds
the ends outside): Guimond et al. 2023, whose mantle storage saturates at 2 to 8 Earth oceans for 0.5 to 5 Earth masses (read by way of the GEN.1146
report). Stand-in row HY1."""
MANTLE_BUFFER_RANGE = (1.0, 8.0)
COVER_SNAP = 0.999
"""float: A cover this high is all of the surface: the last tenth of a percent of a Gaussian tail is a seamount, not an island."""
HYPSOMETRY = ((0.6, -4.5, 1.0), (0.4, 0.45, 0.9))
"""tuple: The surface of a rocky body as two Gaussians of height, each `(share of the area, mean height km, spread km)` at Earth's gravity: ocean crust and
continental crust (Cowan and Abbot 2014: two levels, relief as 1 / g; Guimond et al.: continental crust is about 40 percent of Earth's surface).
The four heights and spreads are stand-ins tuned only so that Earth's 2.7 km water layer covers 0.71 of the surface (L2). Stand-in row HY2."""

BRINE_MAX_WATER_FRACTION = 2e-5
ACID_SO2_SHARE = 0.01
SODA_MIN_CO2_KPA = 10.0
SODA_MIN_OCEAN_FRACTION = 0.9
HYCEAN_MIN_H2_SHARE = 0.5

OCEAN_CHEMISTRY = {
    "ice-sealed": ((3.0, 5.5), (0.99, 1.0)),
    "chloride brine": ((3.5, 6.5), (0.30, 0.65)),
    "acid sulfate": ((1.0, 4.5), (0.90, 0.98)),
    "soda": ((9.0, 11.5), (0.92, 0.99)),
    "neutral": ((6.5, 8.5), (0.97, 0.99)),
}
"""dict: Each class's `(pH range, water activity range)`, drawn uniformly
(section 5.4; neutral's water activity is 20 to 50 g/kg of salt)."""

ICY_MIN_WATER_FRACTION = 0.01
"""float: The water mass fraction from which a body's surface is ice (an icy body holds guest gases in the cages of its ice skin, `surface_air`)."""


def water_boiling_k(pressure_pa):
    """Water's boiling point at `pressure_pa`, from the one saturation function (`physics/saturation.py`), capped at its critical
    point: above that there is no liquid, only a supercritical fluid. Under the triple-point pressure it is the sublimation
    temperature, below the melting point."""
    if pressure_pa <= 0:
        return WATER_MELTING_K
    critical = gases.gas("H2O").critical_temperature_k
    if pressure_pa >= saturation.pressure_pa("H2O", critical - 1e-6):
        return critical
    low, high = 1.0, critical
    for _ in range(80):
        middle = 0.5 * (low + high)
        if saturation.pressure_pa("H2O", middle) < pressure_pa:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def ice_shell_thickness_km(flux_w_m2, top_temperature_k, base_temperature_k=WATER_MELTING_K):
    """The conducting lid: with k = A / T, Fourier's law integrates to
    D = (A / F) ln(T_base / T_top). 0 when the top is already at melting. The flux an icy body's shell conducts includes its own tide, which
    in turn depends on this thickness (`icy_tide`, GEN.960): the stored `internal + tidal` flux is the one in equilibrium with it."""
    if top_temperature_k >= base_temperature_k:
        return 0.0
    return ICE_CONDUCTIVITY_A_W_M / flux_w_m2 * math.log(base_temperature_k / top_temperature_k) / 1000.0


def pressure_melted_shell_km(flux_w_m2, top_temperature_k, gravity_m_s2, melting_k=WATER_MELTING_K):
    """`ice_shell_thickness_km` with the base melting point (`melting_k`: pure water's, or a brine's eutectic) lowered by the
    lid's own weight (fix 3 of section 5.3). The fall stops at the ice III triple point (`water_layer.LIQUIDUS_K_MPA`'s first point, 251.2 K): under
    more pressure than that the base meets ice III, not a deeper Ih melt (GEN.1222: with no floor a thick lid of a cold, low-heat body pulled its own
    base under the surface temperature, the shell came out 0 and the body a liquid ocean at 25 K)."""
    base_k = melting_k
    floor_k = min(melting_k, water_layer.LIQUIDUS_K_MPA[0][0])
    shell_km = ice_shell_thickness_km(flux_w_m2, top_temperature_k, base_k)
    for _ in range(3):
        bar = ICE_DENSITY_KG_M3 * gravity_m_s2 * shell_km * 1000.0 / 1e5
        base_k = max(melting_k - MELTING_SLOPE_K_PER_BAR * bar, floor_k)
        shell_km = ice_shell_thickness_km(flux_w_m2, top_temperature_k, base_k)
    return shell_km


def water_layer_depth_km(mass_kg, radius_km, water_mass_fraction):
    """How deep `water_mass_fraction` of the body's mass would lie spread
    over its surface, at liquid water's density (ignores compression): the
    volume the hypsometry of `ocean_fraction` is about. The depth of the
    column itself, compressed, is `water_layer.column`'s."""
    area_m2 = 4 * math.pi * (radius_km * 1000.0) ** 2
    return water_mass_fraction * mass_kg / (LIQUID_WATER_DENSITY_KG_M3 * area_m2) / 1000.0


def _normal_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _normal_pdf(x):
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


def _filled_volume_km(level_km, scale):
    """The water per unit area, km, that fills the surface of `HYPSOMETRY` (relief times `scale`) to `level_km`: the mean of `max(0, level - height)`."""
    total = 0.0
    for share, mean, spread in HYPSOMETRY:
        mean, spread = mean * scale, spread * scale
        a = (level_km - mean) / spread
        total += share * ((level_km - mean) * _normal_cdf(a) + spread * _normal_pdf(a))
    return total


def ocean_fraction(layer_km, gravity_m_s2, icy=False):
    """
    The share of the surface a global water layer of `layer_km` covers (GEN.1147): the level `z` at which the water fills the two-level surface of
    `HYPSOMETRY` (found by bisection, the volume rising with the level), and the share of the surface under `z`. The relief is scaled as 1 / g.
    For `icy` (GEN.833) it is at most Earth's: observed icy relief is a few to 20 km on bodies of any gravity (Enceladus, Callisto, Iapetus), ice
    being strength-limited, so it does not grow as 1 / g below Earth's.
    """
    if layer_km <= 0 or gravity_m_s2 <= 0:
        return 0.0
    scale = constants.EARTH_GRAVITY / gravity_m_s2
    if icy:
        scale = min(scale, 1.0)
    low = min(mean - 6.0 * spread for _share, mean, spread in HYPSOMETRY) * scale
    high = max(mean + 6.0 * spread for _share, mean, spread in HYPSOMETRY) * scale + layer_km
    if _filled_volume_km(high, scale) <= layer_km:
        return 1.0
    for _ in range(60):
        middle = 0.5 * (low + high)
        if _filled_volume_km(middle, scale) < layer_km:
            low = middle
        else:
            high = middle
    level = 0.5 * (low + high)
    covered = sum(share * _normal_cdf((level - mean * scale) / (spread * scale)) for share, mean, spread in HYPSOMETRY)
    return 1.0 if covered >= COVER_SNAP else covered


def frozen_cover(body, gravity_m_s2):
    """
    The share of its surface the water of a frozen body covers as ice (GEN.833): the cover (`ocean_fraction`'s law) of the layer all of the body's
    water makes at ice's density. The surface water (`draw_surface_water`, what the mantle leaves) is the reservoir of a liquid ocean;
    a frozen body holds its water as ice through its outer layers, so the ice that covers it is the water it holds, not the part
    a liquid ocean would have had. A small moon with a thick ice shell is covered (`ocean_fraction`'s `icy` relief).
    """
    area_m2 = 4 * math.pi * (body.radius * 1000.0) ** 2
    layer_km = (body.water_mass_fraction or 0.0) * body.mass / (ICE_DENSITY_KG_M3 * area_m2) / 1000.0
    return ocean_fraction(layer_km, gravity_m_s2, icy=True)


def mantle_max_water_fraction(mass_kg):
    """The most water, as a mass fraction of the body, a mantle stores (GEN.1147): `MANTLE_CAP_OCEANS` joined as a power law of the mass, the ends held outside."""
    (low_mass, low_oceans), (high_mass, high_oceans) = MANTLE_CAP_OCEANS
    earths = mass_kg / constants.EARTH_MASS_TO_KG
    exponent = math.log(high_oceans / low_oceans) / math.log(high_mass / low_mass)
    oceans = low_oceans * (min(max(earths, low_mass), high_mass) / low_mass) ** exponent
    return oceans * EARTH_OCEAN_FRACTION


def surface_water_fraction(water_mass_fraction, mass_kg, rng=draw):
    """Rule 1: what the mantle leaves at the surface of a body of `mass_kg`."""
    buffer = log_uniform(*MANTLE_BUFFER_RANGE, rng=rng)
    return max(water_mass_fraction / buffer, water_mass_fraction - mantle_max_water_fraction(mass_kg))


def ocean_class(body, hp_ice_km, surface_fraction, covered):
    """Rule 5's class for a liquid ocean."""
    if hp_ice_km > 0:
        return "ice-sealed"
    if surface_fraction < BRINE_MAX_WATER_FRACTION:
        return "chloride brine"
    so2 = atmosphere.partial_kpa(body, "so2")
    co2 = atmosphere.partial_kpa(body, "co2")
    redox = getattr(body, "mantle_redox", None)
    if redox == "oxidized" and so2 > ACID_SO2_SHARE * (so2 + co2):
        return "acid sulfate"
    if redox in ("reduced", "intermediate") and co2 > SODA_MIN_CO2_KPA and covered > SODA_MIN_OCEAN_FRACTION:
        return "soda"
    return "neutral"


def _clear(body):
    """Clears what the hydrosphere sets; the water mass fraction is the bulk composition's (`bulk_draw`) and stays."""
    for field in HYDROSPHERE_FIELDS:
        if field != "water_mass_fraction":
            setattr(body, field, None)


def water_loss_kg(body):
    """The water escape has stripped from a body, kg (GEN.639; `runaway_history.loss_kg`; no draws)."""
    return runaway_history.loss_kg(getattr(body, "star", None), getattr(body, "_star_distance_au", None), body.mass, body.radius)


def draw_surface_water(body):
    """
    The water mass fraction at the surface of a rocky body: what its mantle leaves of the bulk water (rule 1, a seeded draw of the
    mantle's buffer), less what escape has stripped (GEN.639), plus the part of the mantle's that volcanism returns. Drawn before the
    air (`planets.calculate_atmospheric_conditions`): the ocean is the reservoir the air's water vapour stands in contact with, and the
    carbonate cycle needs it, so the air and the ocean are solved together.

    Raises:
        ValueError: The body has no water mass fraction (its bulk composition is drawn first).
    """
    w = body.water_mass_fraction
    if w is None:
        raise ValueError("the body has no water mass fraction: its bulk composition is drawn first")
    surface_w = surface_water_fraction(w, body.mass)
    lost_kg = water_loss_kg(body)
    if lost_kg:  # GEN.639: escape strips the air and the oceans, and volcanism later returns part of the mantle's water
        lost = min(surface_w, lost_kg / body.mass)
        mantle = w - surface_w  # volcanism returns part of the mantle's water, never more than was lost
        surface_w = surface_w - lost + min(lost, tuning.MANTLE_OUTGAS_SHARE * mantle)
    return surface_w


def ocean_cover(body, surface_w):
    """The share of the surface a body's surface water (mass fraction `surface_w`) would cover (rule 1 and 2: `ocean_fraction` of the
    layer it makes)."""
    return ocean_fraction(water_layer_depth_km(body.mass, body.radius, surface_w), (body.gravity or 0.0) * constants.EARTH_GRAVITY)


def generate_hydrosphere(body, internal_flux_w_m2, surface_w):
    """Sets `HYDROSPHERE_FIELDS` on a planet or moon from its class, mass,
    radius, surface temperature (K), pressure (Pa) and gases. Any ice lid
    grows to match `internal_flux_w_m2` (plus the body's tide, `tidal_flux_w_m2`, set with the heat model by
    `interior_heat.update_heat`, GEN.737), and a world with no liquid at its surface may hold
    water underground (`aquifers`, GEN.654). `surface_w` is the surface water mass
    fraction (`draw_surface_water`). Call after its atmosphere."""
    _clear(body)
    if body.body_type == "g":
        return
    _generate_surface(body, internal_flux_w_m2, surface_w)
    aquifers.generate(body, surface_w, internal_flux_w_m2)


def _generate_surface(body, internal_flux_w_m2, surface_w):
    """The surface part of `generate_hydrosphere`: the ocean, its ice and its chemistry; returns `surface_w`."""
    w = body.water_mass_fraction
    gravity = (body.gravity or 0.0) * constants.EARTH_GRAVITY
    layer_km = water_layer_depth_km(body.mass, body.radius, surface_w)
    covered = ocean_fraction(layer_km, gravity)
    body.ocean_fraction = covered
    body.land_fraction = 1.0 - covered
    temperature = body.surface_temperature or 0.0
    pressure = body.atmospheric_pressure or 0.0
    if not covered:
        body.hydrosphere = "dry"
        return surface_w

    if pressure < saturation.WATER_TRIPLE_POINT_PA and temperature >= VACUUM_ICE_MAX_K:
        body.hydrosphere = "dry"
        body.ocean_fraction, body.land_fraction = 0.0, 1.0
        return surface_w
    water_kg = surface_w * body.mass
    salt = None
    if temperature >= WATER_MELTING_K:
        if temperature >= water_boiling_k(pressure):
            body.hydrosphere = "vapour"
            body.ocean_fraction, body.land_fraction = 0.0, 1.0
            return surface_w
        column = water_layer.column(body.mass, body.radius, water_kg, covered=covered, bottom_temperature_k=temperature,
                                    pressure_pa=pressure)
        h2 = atmosphere.partial_kpa(body, "h2")
        state = "hycean" if h2 > HYCEAN_MIN_H2_SHARE * pressure / 1000.0 else "surface ocean"
        if not column.liquid_km:    # GEN.1222: under the pressure of a thick air the warm water is high-pressure ice all the way, not an ocean
            state = "ice"
    else:
        flux = internal_flux_w_m2 + (getattr(body, "tidal_flux_w_m2", None) or 0.0)  # GEN.737: the tide stored with the heat model (GEN.700: the heat an ice-covered ocean's life can use)
        melting_k = WATER_MELTING_K
        if surface_w < BRINE_MAX_WATER_FRACTION:
            salt, melting_k = clathrate_brine.draw_brine()  # GEN.762: little water, so the salts are concentrated and the ice melts at their eutectic
        shell = pressure_melted_shell_km(max(flux, 1e-6), max(temperature, 1.0), gravity, melting_k)
        column = water_layer.column(body.mass, body.radius, water_kg, covered=covered, shell_km=shell,
                                    bottom_temperature_k=melting_k, pressure_pa=pressure)
        if not column.liquid_km:
            state = "ice"
        else:  # a brine warmer than its eutectic is liquid at the surface (GEN.762): no lid to cover it
            state = "ice-covered ocean" if column.ice_shell_km else "surface ocean"
    shell_km, liquid_km, hp_km = column.ice_shell_km, column.liquid_km or None, column.hp_ice_km
    if state == "ice" and not shell_km:     # GEN.1222: water that is high-pressure ice all the way (no lid over a liquid): the ice is the shell
        shell_km = hp_km
    body.hydrosphere = state
    if state in FROZEN_COVER_STATES:    # GEN.833: the ice covers what the water the body holds makes, not the toy ocean cover
        body.ocean_fraction = max(covered, frozen_cover(body, gravity))
        body.land_fraction = 1.0 - body.ocean_fraction
    body.ice_shell_km = shell_km or None
    body.hp_ice_km = hp_km or None
    body.ocean_depth_km = liquid_km
    body.water_base_pressure_pa = column.base_pressure_pa
    if not liquid_km:
        return surface_w
    kind = ocean_class(body, hp_km, surface_w, covered)
    if kind == "chloride brine" and state == "surface ocean":
        body.hydrosphere = BRINE_SEEP
    (ph_lo, ph_hi), (aw_lo, aw_hi) = OCEAN_CHEMISTRY[kind]
    body.ocean_class = kind
    body.ocean_ph = draw.uniform(ph_lo, ph_hi)
    body.water_activity = draw.uniform(aw_lo, aw_hi)
    if kind == "chloride brine" and salt is not None and clathrate_brine.water_activity(salt) is not None:
        body.water_activity = clathrate_brine.water_activity(salt)  # GEN.762: the drawn salt's own figure
    body.phosphorus = "high" if kind == "soda" else "starved" if kind == "ice-sealed" else "limited"
    return surface_w
