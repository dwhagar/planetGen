# planetgen/physics/aerosols.py

"""
Haze, dust and cloud of an atmosphere (GEN.414; docs/design/atmosphere-colour.md section 10, part 06 sections 3.3 to 3.9).

Three things scatter light in a planet's air besides the gas, each with an optical depth at 550 nm that the body stores
(`tau_cloud`, `tau_haze`, `tau_dust`, `cloud_cover`, `aerosol_draw`, DB.68). `draw_aerosols` makes them from the air column and the
surface temperature; `layers` turns them into the `AerosolLayer`s `physics/atmospheric_radiation.surface_radiation` takes.

**Cloud (10.1).** The air is cooled along a dry adiabat from the ground to the skin temperature (`Ts / 2^(1/4)`), isothermal above. A
condensable gas `i` with mole fraction `x_i` forms a cloud where `x_i p` reaches its saturation pressure (the first crossing going up;
at the ground when the ground is already saturated). The most condensate that can fall out is the vapour above the base less what
the saturated air still holds, `m_c = (M_i / (M g)) integral (x_i - min(x_i, p_sat(T(p)) / p)) dp`, and the in-cloud depth is
`tau = 3 Q_ext m_c / (4 rho_p r_eff)` (Q_ext = 2). The droplet radius is log-uniform on 1 to 20 micron; the depth is log-uniform
from 0.1 to its limit; `cloud_cover` is 0 without a condensate, else uniform on 0.1 to 1. Clouds are grey with a single-scattering
albedo of 1 and an asymmetry of 0.85 (the Monte Carlo table of part 06 3.10).

**Haze (10.2, stand-in).** Above a CH4/CO2 ratio of 0.1 (thick above 0.2; a tenth of the production in hydrogen-rich air) methane
photolysis makes a haze. The supply is diffusion limited, `N K_zz f_CH4 / H_a`, with `N` the air at the level where the methane
column is one Lyman-alpha optical depth, `K_zz` log-uniform on 1e6 to 1e10 cm2/s (a nuisance parameter, decision 12), times a yield
log-uniform on 0.01 to 1. The steady depth is the particles' mass flux over their fall speed (Stokes with the Cunningham slip
correction) summed down the air column, `tau = 3 Q_ext Mdot / (4 rho_p r_eff) integral dz / v_f`, with the particle radius
log-uniform on 0.05 to 0.3 micron and the density of the Titan tholin (1000 kg/m3). The Titan laws (9.6 at 500 nm) are a test
target, never a preset. The extinction goes as `lambda^-a` with `a` between 0.97 and 2.34 (aggregates); the absorption follows the
Khare tholin imaginary index at 450, 550 and 650 nm.

**Dust (10.3, stand-in).** A dry solid surface lofts dust when the friction velocity is above the threshold; the wind has no source,
so `u*/u*_t` is log-uniform on 0.3 to 3 and dust exists above 1, with a depth log-uniform on 0.5 to 1.5, grey, a radius of 1.2 to
1.8 micron, an asymmetry of 0.673 +- 0.081 and the single-scattering albedo of the Wolff 2009 spectra (digitised from its Figures 12
to 14, a labelled stand-in until a table is read).

`aerosol_draw` is one number in [0, 1) that picks every within-bracket optical property (the dust radius and asymmetry, the haze
exponent), so the stored five values repeat the optics with no other column. Draws come from `rng` (an object with `random()`),
always four for the cloud, four for the haze, two for the dust and one for `aerosol_draw`, whether or not the body has the aerosol, so one body's draws never
shift another's.

Units: Pa, K, kg, m, s; wavelengths in nm.
"""

import math
from dataclasses import dataclass

import numpy as np

from planetgen.physics import air_column, atmosphere, atmospheric_radiation, constants, gases
from planetgen.util import draw

SKIN_FACTOR = 2.0 ** -0.25
"""float: The skin (isothermal upper air) temperature as a share of the surface temperature; the grey-atmosphere limit."""

CLOUD_Q_EXT = 2.0
"""float: The extinction efficiency of droplets above a micron (part 06 3.5)."""
CLOUD_RADIUS_UM = (1.0, 20.0)
"""tuple: The cloud droplet effective radius range, micron (open; Earth model droplets 10 to 17.5, Groot 2020)."""
CLOUD_TAU_FLOOR = 0.1
"""float: The lowest drawn cloud depth (atmosphere-colour.md 10.1)."""
CLOUD_COVER_RANGE = (0.1, 1.0)
"""tuple: The cloud cover draw (open)."""
CLOUD_G = 0.85
"""float: The cloud asymmetry factor of the Monte Carlo reflectance table (part 06 3.10)."""
CLOUD_OMEGA = 1.0
"""float: The cloud single-scattering albedo of the same table."""
DRAWS = 11
"""int: What `draw_aerosols` takes from its stream: four for the cloud, four for the haze, two for the dust, one for `aerosol_draw`."""

CONDENSATE_DENSITY_KG_M3 = {
    # part 06 table 2.8: water liquid 997 (recalled), NH3 ice 840 and CH4 ice 490 (virga), CO2 ice 1560 (recalled)
    "H2O": 997.0, "NH3": 840.0, "CH4": 490.0, "CO2": 1560.0,
}
"""dict: The clouds the model forms, by gas-table formula, with the particle density. Other condensables (N2, CO, SO2, H2S, O2) have no
row in part 06 2.8 and form no cloud until GEN.452 adds the condensate registry."""

HAZE_ONSET_CH4_CO2 = 0.1
"""float: The CH4/CO2 ratio above which a methane haze forms (Arney 2016; part 06 2.9)."""
HAZE_THICK_CH4_CO2 = 0.2
HAZE_H2_CUT = 0.1
"""float: The share of the production left in hydrogen-rich air (DeWitt 2009, about a tenth)."""
HAZE_H2_FRACTION = 0.3
"""float: The H2 mole fraction above which the air counts as hydrogen-rich (open)."""
HAZE_YIELD_RANGE = (0.01, 1.0)
HAZE_KZZ_RANGE_M2_S = (1e2, 1e6)
"""tuple: K_zz, 1e6 to 1e10 cm2/s (Gao 2018: 1e6 to 1e8; Adams 2019: 1e8 to 1e10)."""
HAZE_RADIUS_UM = (0.05, 0.3)
HAZE_DENSITY_KG_M3 = 1000.0
"""float: Titan tholin, 1 g/cm3 (assumed by the haze papers; part 06 2.8)."""
HAZE_Q_EXT = 1.0
HAZE_EXPONENT_RANGE = (0.97, 2.34)
HAZE_G = 0.65
"""float: The haze asymmetry factor (a stand-in between the 0.1 of small spheres and the 0.7 of aggregates)."""
HAZE_K_550 = (0.057, 0.0256, 0.012)
"""tuple: The Khare tholin imaginary index at 450, 550 and 650 nm (part 06 2.8)."""
HAZE_OMEGA_550 = 0.87
"""float: The tholin single-scattering albedo at 550 nm for 0.1 micron grains (part 06 2.8: 0.873)."""
LYMAN_ALPHA_CH4_CM2 = 1.8e-17
"""float: The methane cross section at Lyman alpha (recalled), which sets the level the methane is split at."""

DUST_RATIO_RANGE = (0.3, 3.0)
DUST_TAU_RANGE = (0.5, 1.5)
DUST_RADIUS_UM = (1.2, 1.8)
DUST_G = (0.673, 0.081)
DUST_DENSITY_KG_M3 = 2500.0
"""float: The density of a dust grain, kg/m3 (DESIGN: basaltic dust, Mars's is measured at about 2,500 to 3,000)."""
DUST_Q_EXT = 2.5
"""float: The extinction efficiency of a micrometre grain at 550 nm (DESIGN; 2 to 3 for the sizes of `DUST_RADIUS_UM`)."""
"""tuple: The dust asymmetry factor and its spread (Wolff 2009)."""
DUST_OMEGA_NM = (300.0, 380.0, 440.0, 460.0, 500.0, 600.0, 700.0, 800.0, 2480.0)
DUST_OMEGA_R12 = (0.705, 0.705, 0.765, 0.780, 0.825, 0.941, 0.968, 0.965, 0.965)
DUST_OMEGA_R18 = (0.705, 0.705, 0.772, 0.788, 0.831, 0.945, 0.971, 0.972, 0.972)
"""The Wolff 2009 single-scattering albedo for a dust radius of 1.2 and 1.8 micron, digitised from its Figures 12 to 14
(accuracy 0.0005 to 0.001); under 440 nm the UV-band mean of 0.67 to 0.74 (partly verified)."""


@dataclass(frozen=True)
class Aerosols:
    """
    The stored aerosol values of a body.

    Attributes:
        tau_haze (float): Haze optical depth at 550 nm (0 for none).
        tau_dust (float): Dust optical depth at 550 nm.
        tau_cloud (float): In-cloud optical depth at 550 nm (grey).
        cloud_cover (float): The share of the sky under cloud, 0 to 1 (0 with no condensate).
        aerosol_draw (float): The number in [0, 1) that picks the optics inside the brackets.
    """

    tau_haze: float
    tau_dust: float
    tau_cloud: float
    cloud_cover: float
    aerosol_draw: float


NONE = Aerosols(0.0, 0.0, 0.0, 0.0, 0.5)
"""Aerosols: What a body with no air has."""


def _log_uniform(rng, low, high):
    return math.exp(math.log(low) + rng.random() * (math.log(high) - math.log(low)))


def _cp_j_mol_k(column, temperature_k):
    """The mole-fraction-weighted ideal-gas heat capacity of the mixture at `temperature_k`, linear between the table's 100, 300 and 1000 K."""
    total = 0.0
    for formula, pressure in column.partials_pa.items():
        table = gases.gas(formula).cp_j_mol_k
        points = sorted((float(key), value) for key, value in table.items())
        total += pressure / column.pressure_pa * float(np.interp(temperature_k, [p[0] for p in points], [p[1] for p in points]))
    return total


def saturation_pa(formula, temperature_k):
    """
    The vapour pressure over the condensate of `formula` at `temperature_k`, Pa. The two-point Clausius-Clapeyron through
    `atmosphere.CONDENSATION_K` where it has the gas (GEN.452 replaces it with the checked curves), else through the gas table's triple
    point and latent heat (ammonia). Infinite at and above the critical temperature, where nothing condenses.
    """
    gas = gases.gas(formula)
    if gas.critical_temperature_k and temperature_k >= gas.critical_temperature_k:
        return math.inf
    key = formula.lower()
    if key in atmosphere.CONDENSATION_K:
        value = atmosphere.vapour_pressure_kpa(key, temperature_k) * 1000.0
    else:
        heat = gas.latent_heat_sublimation_kj_mol if temperature_k < gas.triple_point_k else gas.latent_heat_vaporisation_kj_mol
        if isinstance(heat, (list, tuple)):
            heat = 0.5 * (heat[0] + heat[1])
        heat = heat or gas.latent_heat_vaporisation_kj_mol
        value = gas.triple_point_pa * math.exp(-heat * 1000.0 / constants.R * (1.0 / temperature_k - 1.0 / gas.triple_point_k))
    return value


def profile(column, temperature_k, levels=240, decades=8.0):
    """
    The pressure and temperature of `levels` points from the ground (index 0) up: the dry adiabat from `temperature_k`, isothermal at the
    skin temperature above.

    Returns:
        tuple: `(pressure_pa, temperature_k)`, arrays.
    """
    skin = temperature_k * SKIN_FACTOR
    cp = _cp_j_mol_k(column, 0.5 * (temperature_k + skin))
    kappa = constants.R / cp
    pressure = column.pressure_pa * 10.0 ** (-np.linspace(0.0, decades, levels))
    temperature = np.maximum(temperature_k * (pressure / column.pressure_pa) ** kappa, skin)
    return pressure, temperature


def condensate_columns_kg_m2(column, temperature_k):
    """
    The most condensate each cloud-forming gas can put in the air, kg per m2 of ground: `{formula: m_c}` for the gases that reach saturation.
    """
    pressure, temperature = profile(column, temperature_k)
    weight = np.gradient(-pressure)
    result = {}
    for formula in CONDENSATE_DENSITY_KG_M3:
        fraction = column.fraction(formula)
        if fraction <= 0.0:
            continue
        saturation = np.array([saturation_pa(formula, t) for t in temperature])
        if not np.any(fraction * pressure >= saturation):
            continue
        vapour = np.minimum(fraction, saturation / pressure)
        grams = gases.gas(formula).molar_mass / column.mean_molar_mass_kg_mol / column.gravity_m_s2
        mass = grams * float(np.sum((fraction - vapour) * weight))
        if mass > 0.0:
            result[formula] = mass
    return result


def cloud_depth_limit(column, temperature_k, radius_um):
    """The in-cloud optical depth at 550 nm if every gas that can condense did, droplets of `radius_um`: the top of the cloud draw."""
    total = 0.0
    for formula, mass in condensate_columns_kg_m2(column, temperature_k).items():
        total += 3.0 * CLOUD_Q_EXT * mass / (4.0 * CONDENSATE_DENSITY_KG_M3[formula] * radius_um * 1e-6)
    return total


def draw_cloud(column, temperature_k, rng):
    """`(tau_cloud, cloud_cover)`: four draws of `rng` whether or not the air condenses."""
    radius = _log_uniform(rng, *CLOUD_RADIUS_UM)
    fraction = rng.random()
    cover = CLOUD_COVER_RANGE[0] + rng.random() * (CLOUD_COVER_RANGE[1] - CLOUD_COVER_RANGE[0])
    rng.random()
    limit = cloud_depth_limit(column, temperature_k, radius)
    if limit <= 0.0:
        return 0.0, 0.0
    low = min(CLOUD_TAU_FLOOR, limit)
    return math.exp(math.log(low) + fraction * (math.log(limit) - math.log(low))), cover


def _viscosity_pa_s(column, temperature_k):
    """The mole-fraction-weighted dilute-gas viscosity of the mixture (the Chapman-Enskog value where the table has constants, else the
    NIST value at 300 K for a polar gas; the gases with neither are left out)."""
    total, weight = 0.0, 0.0
    for formula, pressure in column.partials_pa.items():
        gas = gases.gas(formula)
        try:
            value = gas.viscosity_upa_s(temperature_k)
        except gases.NoData:
            value = gas.viscosity.nist_300_upa_s if gas.viscosity and gas.viscosity.nist_300_upa_s else None
        if value:
            total += pressure * value
            weight += pressure
    return (total / weight if weight else 1.5e1) * 1e-6


def haze_fall_integral_s_m(column, temperature_k, radius_m, density_kg_m3, top_pressure_pa, kzz_m2_s):
    """`integral dz / v_f` from the ground to `top_pressure_pa`, the time in seconds a grain takes to fall through the air
    column: the Stokes speed with the Cunningham slip correction, or the eddy speed `K_zz / H` where that is faster (a grain stays aloft where
    the mixing time `H^2 / K_zz` is shorter than the fall time `H / v`), and hydrostatic `dz = dp / (rho g)`."""
    pressure, temperature = profile(column, temperature_k, levels=400, decades=min(12.0, math.log10(column.pressure_pa / top_pressure_pa)))
    molar = column.mean_molar_mass_kg_mol
    rho = pressure * molar / (constants.R * temperature)
    eta = np.array([_viscosity_pa_s(column, t) for t in temperature])
    path = eta / (0.499 * rho * np.sqrt(8.0 * constants.R * temperature / (math.pi * molar)))
    kn = path / radius_m
    slip = 1.0 + 1.257 * kn + 0.400 * kn * np.exp(-1.10 / kn)
    speed = 2.0 / 9.0 * density_kg_m3 * column.gravity_m_s2 * radius_m ** 2 * slip / eta
    scale_height = constants.R * temperature / (molar * column.gravity_m_s2)
    integrand = 1.0 / (rho * column.gravity_m_s2 * np.maximum(speed, kzz_m2_s / scale_height))
    return float(np.sum(0.5 * (integrand[1:] + integrand[:-1]) * -np.diff(pressure)))


def haze_supply_kg_m2_s(column, temperature_k, kzz_m2_s):
    """The diffusion-limited methane supply as a carbon-bearing mass flux, `N K_zz f_CH4 / H_a`, kg per m2 per s, at the level where the
    methane column is one Lyman-alpha optical depth."""
    fraction = column.fraction("CH4")
    if fraction <= 0.0:
        return 0.0
    skin = temperature_k * SKIN_FACTOR
    molar = column.mean_molar_mass_kg_mol
    number_cm2 = 1.0 / LYMAN_ALPHA_CH4_CM2
    level_pa = number_cm2 * 1e4 / fraction * molar * column.gravity_m_s2 / 6.02214076e23
    level_pa = min(level_pa, column.pressure_pa)
    density = level_pa / (constants.BOLTZMANN * skin)
    scale_height = constants.R * skin / (molar * column.gravity_m_s2)
    flux = density * kzz_m2_s * fraction / scale_height
    return flux * gases.gas("CH4").molar_mass / 6.02214076e23


def haze_gate(column):
    """The share of the full haze production the air allows: 0 below a CH4/CO2 ratio of 0.1, 1 above 0.2 (linear between), a tenth of it in
    hydrogen-rich air."""
    methane = column.fraction("CH4")
    carbon_dioxide = max(column.fraction("CO2"), 1e-12)
    ratio = methane / carbon_dioxide
    share = min(1.0, max(0.0, (ratio - HAZE_ONSET_CH4_CO2) / (HAZE_THICK_CH4_CO2 - HAZE_ONSET_CH4_CO2)))
    if column.fraction("H2") > HAZE_H2_FRACTION:
        share *= HAZE_H2_CUT
    return share


def draw_haze(column, temperature_k, rng):
    """`tau_haze`: the steady haze depth at 550 nm. Four draws of `rng` (K_zz, yield, radius, spare) whether or not the air hazes."""
    kzz = _log_uniform(rng, *HAZE_KZZ_RANGE_M2_S)
    yield_ = _log_uniform(rng, *HAZE_YIELD_RANGE)
    radius = _log_uniform(rng, *HAZE_RADIUS_UM) * 1e-6
    rng.random()
    gate = haze_gate(column)
    if gate <= 0.0:
        return 0.0
    mdot = haze_supply_kg_m2_s(column, temperature_k, kzz) * yield_ * gate
    top = min(column.pressure_pa * 1e-8, 1e-4)
    fall = haze_fall_integral_s_m(column, temperature_k, radius, HAZE_DENSITY_KG_M3, max(top, 1e-12), kzz)
    return 3.0 * HAZE_Q_EXT * mdot * fall / (4.0 * HAZE_DENSITY_KG_M3 * radius)


def dust_loading_mg_m3(tau_dust, scale_height_km):
    """
    The airborne dust at the ground, mg/m3, from the optical depth `tau_dust` and the air's scale height (GEN.658): the
    grains' column mass `4 rho r tau / (3 Q)` (`r` the middle of `DUST_RADIUS_UM`) spread through one scale height. 0 for no
    dust, `None` when the scale height is unknown. Mars at `tau` 0.5 to 1 comes out near the 0.1 mg/m3 of its dust limit.
    """
    if not tau_dust or tau_dust <= 0.0:
        return 0.0
    if not scale_height_km or scale_height_km <= 0.0:
        return None
    radius_m = 0.5 * sum(DUST_RADIUS_UM) * 1e-6
    column_kg_m2 = 4.0 * DUST_DENSITY_KG_M3 * radius_m * tau_dust / (3.0 * DUST_Q_EXT)
    return column_kg_m2 / (scale_height_km * 1000.0) * 1e6


def draw_dust(dry_surface, rng):
    """`tau_dust`: two draws of `rng` (the wind ratio and the depth) whether or not the surface lofts dust."""
    ratio = _log_uniform(rng, *DUST_RATIO_RANGE)
    depth = _log_uniform(rng, *DUST_TAU_RANGE)
    if not dry_surface or ratio < 1.0:
        return 0.0
    return depth


def draw_aerosols(column, temperature_k, rng, *, dry_surface=True):
    """
    The five stored aerosol values of a body.

    Args:
        column (AirColumn or None): The air at the level the light is read at (the ground; a giant's 1 bar level). `None` for no air.
        temperature_k (float): The temperature there.
        rng: An object with `random()`; `DRAWS` draws are taken, the same ten whatever the air.
        dry_surface (bool): Whether the ground is dry solid (no ocean); only then does the wind loft dust.
    """
    if column is None:
        for _ in range(DRAWS):
            rng.random()
        return NONE
    cloud, cover = draw_cloud(column, temperature_k, rng)
    haze = draw_haze(column, temperature_k, rng)
    dust = draw_dust(dry_surface, rng)
    return Aerosols(haze, dust, cloud, cover, rng.random())


def _haze_exponent(draw):
    low, high = HAZE_EXPONENT_RANGE
    return low + draw * (high - low)


def _dust_omega(centres, draw):
    low = np.interp(centres, DUST_OMEGA_NM, DUST_OMEGA_R12)
    high = np.interp(centres, DUST_OMEGA_NM, DUST_OMEGA_R18)
    return (1.0 - draw) * low + draw * high


def _haze_absorption_shape(centres):
    """The tholin's absorption at each wavelength over its value at 550 nm: `k(lambda) / lambda`, log-linear in k between 450 and 650 nm."""
    k = np.exp(np.interp(centres, (450.0, 550.0, 650.0), np.log(HAZE_K_550)))
    return (k / HAZE_K_550[1]) * (550.0 / centres)


def layers(aerosols, centres_nm):
    """
    The `AerosolLayer`s of a body for `atmospheric_radiation.surface_radiation`.

    Returns:
        tuple: `(haze, low_aerosol)`: the haze and dust mixed with the lower air (`haze=`), and the cloud (`low_aerosol=`, with
        `cloud_cover=aerosols.cloud_cover`). Either is `None` where its depth is 0.
    """
    centres = np.asarray(centres_nm, dtype=float)
    parts = []
    if aerosols.tau_haze > 0.0:
        extinction = aerosols.tau_haze * (550.0 / centres) ** _haze_exponent(aerosols.aerosol_draw)
        absorption = aerosols.tau_haze * (1.0 - HAZE_OMEGA_550) * _haze_absorption_shape(centres)
        omega = np.clip(1.0 - absorption / extinction, 0.0, 1.0)
        parts.append((extinction, omega, np.full_like(extinction, HAZE_G)))
    if aerosols.tau_dust > 0.0:
        radius_share = aerosols.aerosol_draw
        g = DUST_G[0] + DUST_G[1] * (2.0 * aerosols.aerosol_draw - 1.0)
        parts.append((np.full_like(centres, aerosols.tau_dust), _dust_omega(centres, radius_share), np.full_like(centres, g)))
    mixed = None
    if parts:
        tau = sum(p[0] for p in parts)
        scattering = sum(p[0] * p[1] for p in parts)
        omega = scattering / tau
        g = sum(p[0] * p[1] * p[2] for p in parts) / scattering
        mixed = atmospheric_radiation.AerosolLayer(tau, omega, g)
    cloud = None
    if aerosols.tau_cloud > 0.0 and aerosols.cloud_cover > 0.0:
        cloud = atmospheric_radiation.grey_layer(aerosols.tau_cloud, CLOUD_OMEGA, CLOUD_G, centres)
    return mixed, cloud


BODY_FIELDS = ("tau_haze", "tau_dust", "tau_cloud", "cloud_cover", "aerosol_draw")
"""tuple: What `generate_aerosols` sets on a planet or moon, as stored in `planets` and `moons` (DB.68)."""


def _clear(body):
    for field in BODY_FIELDS:
        setattr(body, field, None)


def stream_of(body):
    """The body's own stream for these draws, seeded from what is final about it before they are made, so no draw of the run's stream
    moves when they are added, taken out or retuned."""
    return draw.Stream(f"aerosols:{body.is_moon}:{body.planet_class}:{body.mass!r}:{body.radius!r}:{body.distance!r}:{body.surface_temperature!r}")


def generate_aerosols(body):
    """
    Sets `BODY_FIELDS` on a planet or moon from its air, its surface temperature and its ocean (call after `hydrosphere`). A body with no
    air or no temperature has every field `None`.
    """
    _clear(body)
    temperature = getattr(body, "surface_temperature", None)
    column = air_column.column_of(body)
    if column is None or not temperature:
        return
    dry = body.body_type != "g" and not (getattr(body, "ocean_fraction", None) or 0.0)
    drawn = draw_aerosols(column, temperature, stream_of(body), dry_surface=dry)
    for field in BODY_FIELDS:
        setattr(body, field, getattr(drawn, field))


def of_body(body):
    """The `Aerosols` a body stores, or `NONE` for a body that has none (no air, or saved before DB.68)."""
    if any(getattr(body, field, None) is None for field in BODY_FIELDS):
        return NONE
    return Aerosols(*(getattr(body, field) for field in BODY_FIELDS))
