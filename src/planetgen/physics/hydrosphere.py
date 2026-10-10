# planetgen/physics/hydrosphere.py

"""
Water on planets, moons and rogue planets (GEN.88).

The design is docs/design/activity-magnetism-radiation-hydrosphere.md
section 5. Both star-system bodies (`generate_hydrosphere`) and rogue
planets (`physics.rogue_surface`) use the functions here.

1. **Inventory.** A rocky body's water mass fraction is drawn log-uniformly
   from its class's range (`WATER_MASS_FRACTION`, planetGen's defaults;
   Earth's oceans are 2.34e-4). The mantle holds part of it: the surface
   keeps `w / b`, `b` log-uniform 1 to 8, but the mantle never holds more
   than ten Earth oceans (`MANTLE_MAX_WATER_FRACTION`). Gas giants have no
   hydrosphere (every field `None`).
2. **Land.** A toy hypsometry calibrated to Earth's 0.71 ocean fraction:
   `f_ocean = min(1, sqrt(2 D / L))`, `D` the global water layer and `L =
   10.5 km * (g_E / g)` the relief (Cowan and Abbot 2014 put the threshold
   for any land at about 0.2 percent water). `ocean_fraction` is the share
   the water covers, liquid or frozen; `land_fraction` the rest.
3. **State.** Under water's triple-point pressure and warmer than
   `VACUUM_ICE_MAX_K` the water has sublimated away (`dry`); above the
   boiling curve it is in the air (`vapour`). Below freezing, an ice lid grows to
   `ice_shell_thickness_km` from the body's internal heat; it seals the
   whole column (`ice`) or floats on an ocean (`ice-covered ocean`). The
   lid is ice (917 kg/m3) and the column water-equivalent, so the lid is
   divided by 1.09 when compared; its base melts at 273.15 K less 0.0074 K
   per bar of ice above it.
4. **High-pressure ice.** No liquid lies deeper than the liquidus
   (`max_liquid_depth_km`, Wagner 2011's triple points); the rest of the
   column is ice VI or VII (`hp_ice_km`), and the ocean is ice-sealed: no
   contact with rock.
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

from planetgen.physics import constants
from planetgen.util import draw
from planetgen.util.random import log_uniform

HYDROSPHERE_FIELDS = (
    "water_mass_fraction", "hydrosphere", "ocean_fraction", "land_fraction",
    "ocean_depth_km", "ice_shell_km", "hp_ice_km",
    "ocean_class", "ocean_ph", "water_activity", "phosphorus",
)
"""tuple: What `generate_hydrosphere` sets on a planet or moon, as stored in
`planets` and `moons`."""

HYDROSPHERE_STATES = ("dry", "vapour", "ice", "ice-covered ocean", "surface ocean", "hycean")
OCEAN_CLASSES = ("ice-sealed", "chloride brine", "acid sulfate", "soda", "neutral")

WATER_MELTING_K = 273.15
WATER_CRITICAL_K = 647.1
WATER_BOILING_K = 373.15
WATER_TRIPLE_POINT_PA = 611.657
WATER_VAPORIZATION_J_MOL = 40.7e3
GAS_CONSTANT_J_MOL_K = 8.314
LIQUID_WATER_DENSITY_KG_M3 = 1000.0
DEEP_WATER_DENSITY_KG_M3 = 1050.0
"""float: The compressed column's mean density, for its bottom pressure."""
ICE_DENSITY_KG_M3 = 917.0
HP_ICE_DENSITY_KG_M3 = 1300.0
"""float: Ice VI (1.31 g/cm3); `hp_ice_km` is that thickness."""
ICE_TO_WATER = LIQUID_WATER_DENSITY_KG_M3 / ICE_DENSITY_KG_M3
"""float: 1.09, a km of ice lid in km of water."""
ICE_CONDUCTIVITY_A_W_M = 567.0
"""float: Water ice conducts heat as k = A / T (Turcotte and Schubert)."""
VACUUM_ICE_MAX_K = 220.0
"""float: Under water's triple-point pressure, exposed ice warmer than this
sublimates away over billions of years (the Moon keeps ice only in cold
traps; Mars's caps survive at 150 to 210 K): such a body is `dry`. A
planetGen default."""
MELTING_SLOPE_K_PER_BAR = 0.0074
"""float: Ice Ih's melting point falls this much per bar."""

LIQUIDUS_K_MPA = ((251.0, 207.5), (256.0, 346.0), (273.3, 626.0), (355.0, 2200.0))
"""tuple: Wagner 2011's triple points, `(K, MPa)`: Ih-L-III, III-L-V,
V-L-VI and VI-L-VII. Liquid is stable only below the line through them."""

EARTH_OCEAN_FRACTION = 2.34e-4
MANTLE_MAX_WATER_FRACTION = 10 * EARTH_OCEAN_FRACTION
MANTLE_BUFFER_RANGE = (1.0, 8.0)
EARTH_RELIEF_KM = 10.5

WATER_MASS_FRACTION = {
    "A": (1e-6, 1e-4), "B": (1e-7, 1e-5), "C": (1e-6, 1e-4), "D": (0.05, 0.5),
    "E": (1e-6, 1e-4), "F": (3e-5, 2e-4), "G": (1e-5, 1e-4), "H": (1e-6, 1e-5),
    "K": (3e-5, 1e-3), "L": (3e-5, 5e-4), "M": (1e-4, 5e-4), "N": (1e-7, 1e-5),
    "O": (1e-2, 0.1), "P": (1e-4, 1e-3), "Q": (1e-5, 1e-3), "S": (1e-6, 1e-4),
    "V": (1e-5, 1e-3),
}
"""dict: A rocky class's water mass fraction range (log-uniform; planetGen's
defaults): dry for the molten, barren and hot classes, Earth's for M, under
10 percent ocean for H, an ocean world for O, an icy body (Europa 8
percent, Ganymede about 45) for D."""
DEFAULT_WATER_MASS_FRACTION = (1e-5, 1e-3)

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

ICY_MOON_K2_Q = 3e-4
ICY_MOON_ECCENTRICITY_RANGE = (1e-3, 1e-2)
SMALL_MOON_MASS_EARTH = 1e-3
RESONANCE_BOOST_RANGE = (1.0, 50.0)
ICY_MIN_WATER_FRACTION = 0.01


def water_boiling_k(pressure_pa):
    """Water's boiling point at `pressure_pa` (Clausius-Clapeyron), capped at
    its critical point: above that there is no liquid, only a supercritical
    fluid."""
    if pressure_pa <= 0:
        return WATER_MELTING_K
    inverse = 1 / WATER_BOILING_K - GAS_CONSTANT_J_MOL_K * math.log(pressure_pa / 101325.0) / WATER_VAPORIZATION_J_MOL
    return WATER_CRITICAL_K if inverse <= 0 else min(1 / inverse, WATER_CRITICAL_K)


def ice_shell_thickness_km(flux_w_m2, top_temperature_k, base_temperature_k=WATER_MELTING_K):
    """The conducting lid: with k = A / T, Fourier's law integrates to
    D = (A / F) ln(T_base / T_top). 0 when the top is already at melting."""
    if top_temperature_k >= base_temperature_k:
        return 0.0
    return ICE_CONDUCTIVITY_A_W_M / flux_w_m2 * math.log(base_temperature_k / top_temperature_k) / 1000.0


def pressure_melted_shell_km(flux_w_m2, top_temperature_k, gravity_m_s2):
    """`ice_shell_thickness_km` with the base melting point lowered by the
    lid's own weight (fix 3 of section 5.3)."""
    base_k = WATER_MELTING_K
    shell_km = ice_shell_thickness_km(flux_w_m2, top_temperature_k, base_k)
    for _ in range(3):
        bar = ICE_DENSITY_KG_M3 * gravity_m_s2 * shell_km * 1000.0 / 1e5
        base_k = WATER_MELTING_K - MELTING_SLOPE_K_PER_BAR * bar
        shell_km = ice_shell_thickness_km(flux_w_m2, top_temperature_k, base_k)
    return shell_km


def water_layer_depth_km(mass_kg, radius_km, water_mass_fraction):
    """How deep `water_mass_fraction` of the body's mass would lie spread
    over its surface, at liquid water's density (ignores compression)."""
    area_m2 = 4 * math.pi * (radius_km * 1000.0) ** 2
    return water_mass_fraction * mass_kg / (LIQUID_WATER_DENSITY_KG_M3 * area_m2) / 1000.0


def liquidus_pressure_pa(bottom_temperature_k):
    """The pressure where water at `bottom_temperature_k` freezes to a
    high-pressure ice: linear between Wagner's triple points, extended past
    both ends."""
    points = LIQUIDUS_K_MPA
    if bottom_temperature_k <= points[1][0]:
        (t0, p0), (t1, p1) = points[0], points[1]
    elif bottom_temperature_k >= points[-2][0]:
        (t0, p0), (t1, p1) = points[-2], points[-1]
    else:
        (t0, p0), (t1, p1) = next((a, b) for a, b in zip(points, points[1:]) if bottom_temperature_k <= b[0])
    mpa = p0 + (bottom_temperature_k - t0) * (p1 - p0) / (t1 - t0)
    return max(mpa, points[0][1]) * 1e6


def max_liquid_depth_km(bottom_temperature_k, gravity_m_s2):
    """The deepest liquid column, km of water, before its bottom reaches the
    liquidus (section 5.2: 73 km at Earth gravity and 280 K)."""
    return liquidus_pressure_pa(bottom_temperature_k) / (DEEP_WATER_DENSITY_KG_M3 * gravity_m_s2) / 1000.0


def ocean_fraction(layer_km, gravity_m_s2):
    """Section 5.1's `f_ocean = min(1, sqrt(2 D / L))`."""
    if layer_km <= 0 or gravity_m_s2 <= 0:
        return 0.0
    relief_km = EARTH_RELIEF_KM * constants.EARTH_GRAVITY / gravity_m_s2
    return min(1.0, math.sqrt(2 * layer_km / relief_km))


def split_column(column_km, shell_km, bottom_temperature_k, gravity_m_s2):
    """A water column (km of water) under an ice lid (km of ice): `(ice
    shell km, liquid km, high-pressure ice km)`, liquid `None` when the lid
    or the high-pressure ice takes it all."""
    lid_km = shell_km / ICE_TO_WATER
    if lid_km >= column_km:
        return column_km * ICE_TO_WATER, None, 0.0
    deepest = max_liquid_depth_km(bottom_temperature_k, gravity_m_s2)
    if column_km <= deepest:
        return shell_km, column_km - lid_km, 0.0
    hp_km = (column_km - deepest) * LIQUID_WATER_DENSITY_KG_M3 / HP_ICE_DENSITY_KG_M3
    liquid = deepest - lid_km
    return shell_km, (liquid if liquid > 0 else None), hp_km


def tidal_heating_w(k2_over_q, primary_mass_kg, radius_km, semi_major_axis_m, eccentricity):
    """Peale's eccentricity tide, `E = (21/2) (k2/Q) G M_p^2 n R^5 e^2 /
    a^6` (Europa 1.4e11 W at k2/Q 3e-4)."""
    n = math.sqrt(constants.G * primary_mass_kg / semi_major_axis_m ** 3)
    return (10.5 * k2_over_q * constants.G * primary_mass_kg ** 2 * n * (radius_km * 1000.0) ** 5
            * eccentricity ** 2 / semi_major_axis_m ** 6)


def surface_water_fraction(water_mass_fraction, rng=draw):
    """Rule 1: what the mantle leaves at the surface."""
    buffer = log_uniform(*MANTLE_BUFFER_RANGE, rng=rng)
    return max(water_mass_fraction / buffer, water_mass_fraction - MANTLE_MAX_WATER_FRACTION)


def ocean_class(body, hp_ice_km, surface_fraction, covered):
    """Rule 5's class for a liquid ocean."""
    if hp_ice_km > 0:
        return "ice-sealed"
    if surface_fraction < BRINE_MAX_WATER_FRACTION:
        return "chloride brine"
    so2 = getattr(body, "p_so2_kpa", None) or 0.0
    co2 = getattr(body, "p_co2_kpa", None) or 0.0
    redox = getattr(body, "mantle_redox", None)
    if redox == "oxidized" and so2 > ACID_SO2_SHARE * (so2 + co2):
        return "acid sulfate"
    if redox in ("reduced", "intermediate") and co2 > SODA_MIN_CO2_KPA and covered > SODA_MIN_OCEAN_FRACTION:
        return "soda"
    return "neutral"


def _clear(body):
    for field in HYDROSPHERE_FIELDS:
        setattr(body, field, None)


def generate_hydrosphere(body, internal_flux_w_m2):
    """Sets `HYDROSPHERE_FIELDS` on a planet or moon from its class, mass,
    radius, surface temperature (K), pressure (Pa) and gases. Any ice lid
    grows to match `internal_flux_w_m2` (plus an icy moon's tides,
    `moon_tidal_flux_w_m2`). Call after its atmosphere."""
    _clear(body)
    if body.body_type == "g":
        return
    w = log_uniform(*WATER_MASS_FRACTION.get(body.planet_class, DEFAULT_WATER_MASS_FRACTION))
    surface_w = surface_water_fraction(w)
    gravity = (body.gravity or 0.0) * constants.EARTH_GRAVITY
    layer_km = water_layer_depth_km(body.mass, body.radius, surface_w)
    covered = ocean_fraction(layer_km, gravity)
    body.water_mass_fraction = w
    body.ocean_fraction = covered
    body.land_fraction = 1.0 - covered
    column_km = layer_km / covered if covered else 0.0
    temperature = body.surface_temperature or 0.0
    pressure = body.atmospheric_pressure or 0.0
    if not covered:
        body.hydrosphere = "dry"
        return

    if pressure < WATER_TRIPLE_POINT_PA and temperature >= VACUUM_ICE_MAX_K:
        body.hydrosphere = "dry"
        body.ocean_fraction, body.land_fraction = 0.0, 1.0
        return
    if temperature >= WATER_MELTING_K:
        if temperature >= water_boiling_k(pressure):
            body.hydrosphere = "vapour"
            body.ocean_fraction, body.land_fraction = 0.0, 1.0
            return
        shell_km, liquid_km, hp_km = split_column(column_km, 0.0, temperature, gravity)
        h2 = getattr(body, "p_h2_kpa", None) or 0.0
        state = "hycean" if h2 > HYCEAN_MIN_H2_SHARE * pressure / 1000.0 else "surface ocean"
    else:
        flux = internal_flux_w_m2 + moon_tidal_flux_w_m2(body, getattr(body, "primary_mass_kg", None), w)
        shell = pressure_melted_shell_km(max(flux, 1e-6), max(temperature, 1.0), gravity)
        shell_km, liquid_km, hp_km = split_column(column_km, shell, WATER_MELTING_K, gravity)
        state = "ice-covered ocean" if liquid_km else "ice"
    body.hydrosphere = state
    body.ice_shell_km = shell_km or None
    body.hp_ice_km = hp_km or None
    body.ocean_depth_km = liquid_km
    if not liquid_km:
        return
    kind = ocean_class(body, hp_km, surface_w, covered)
    (ph_lo, ph_hi), (aw_lo, aw_hi) = OCEAN_CHEMISTRY[kind]
    body.ocean_class = kind
    body.ocean_ph = draw.uniform(ph_lo, ph_hi)
    body.water_activity = draw.uniform(aw_lo, aw_hi)
    body.phosphorus = "high" if kind == "soda" else "starved" if kind == "ice-sealed" else "limited"


def moon_tidal_flux_w_m2(moon, primary_mass_kg, water_mass_fraction):
    """An icy moon's eccentricity-tide flux, W/m2 (k2/Q 3e-4, a forced
    eccentricity log-uniform 0.001 to 0.01, small moons boosted 1 to 50x
    by resonance as Enceladus needs); 0 for a rocky moon."""
    if not moon.is_moon or not primary_mass_kg or water_mass_fraction < ICY_MIN_WATER_FRACTION:
        return 0.0
    eccentricity = log_uniform(*ICY_MOON_ECCENTRICITY_RANGE)
    power = tidal_heating_w(ICY_MOON_K2_Q, primary_mass_kg, moon.radius, moon.distance * constants.AU_TO_M,
                            eccentricity)
    if moon.mass / constants.EARTH_MASS_TO_KG < SMALL_MOON_MASS_EARTH:
        power *= log_uniform(*RESONANCE_BOOST_RANGE)
    return power / (4 * math.pi * (moon.radius * 1000.0) ** 2)
