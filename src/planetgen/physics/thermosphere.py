# planetgen/physics/thermosphere.py

"""
The thermosphere and the thermal escape that comes out of it (GEN.738; docs/design/composition/atmospheric-escape-model.md)
==========================================================================================================================

What a small or hot body loses from its air is set by the temperature at the exobase, how thin the air is there, and how much XUV the
star has put out over the body's life. This module takes those three steps and returns, per species, the share of an air that is still
there after the star's XUV history (`kept_shares`), plus the mass-loss rate of a hydrogen envelope (`hydro_based_mass_loss_g_s`).

1. *The exobase temperature* (`exobase_temperature_k`). A thermosphere heated by the star's XUV and cooled by the infrared emission of
   carbon dioxide: the mesopause temperature plus a rise that goes as the XUV flux and falls with the carbon dioxide of the air. Its
   form and constants are fitted to the 1-D hydrodynamic upper-atmosphere model of Johnstone, Guedel, Lammer and Kislyakova 2018
   (arXiv:1806.06897): the Earth's exobase temperature over a hundred-fold range of CO2 (Fig. 12) and the rise at 4.9 times the present
   XUV (section 4.2). Venus, Mars and Titan are not fitted; they are the checks (the tests).
2. *The exobase* (`exobase`): the height where the mean free path equals the scale height, of the air rising from the mesopause on a
   Bates temperature profile; it gives the density of the gas that can leave.
3. *The flux of a species.* The Jeans flux from the exobase (`jeans_flux`), at most the diffusion-limited flux of Hunten 1973 for a
   species lighter than the air (`diffusion_limited_flux`: its binary diffusion coefficients are the hard-sphere kinetic theory of
   Chapman and Enskog, which reproduces the diffusion of H2 in N2 to 16 percent and the 1e28 s-1 H2 flux of Titan to a factor of four).
   The total mass loss is at most the energy-limited rate (`energy_limited_mass_loss_kg_s`, Watson et al. 1981; Erkaev et al. 2007).
4. *The XUV history* (`kept_shares`): the star's `L_XUV` over its life from `activity` (saturated, then decaying), the flux of the body at
   each age, the exobase temperature at each age, and each species lost at a first-order rate set by its flux over its column. The flux
   is held at `XUV_FLUX_MAX_EARTH` before the star decays to it: the model is not validated above it, and the loss of the hotter early
   phase belongs to the retention rows of the pressure stage (the cosmic shoreline, GEN.511), which `mixture_escape` does not repeat.

*Envelopes.* The mass-loss rate of a hydrogen-dominated planet of 1 to 40 Earth masses is the analytic fit of Kubyshkina et al. 2018
(arXiv:1810.06920) to a grid of 7000 hydrodynamic models; it is checked against the nine planets of its Table 2.

Not modelled: non-thermal escape (sputtering, photochemical and ion escape, which dominate Mars's and Venus's loss of hydrogen), the day
and night difference of a slowly rotating body (the temperature here is the dayside of a fast rotator, the Earth model's global mean),
magnetic shielding (GEN.804 wires the stored field), and the cooling by gases other than CO2.
"""

import functools
import math

from planetgen import tuning
from planetgen.physics import activity, constants, gases, ion_escape

# --- the exobase temperature ------------------------------------------------------------------------------------------------

MESOPAUSE_OVER_TEQ = 0.78
"""float: The temperature of the mesopause over the equilibrium temperature (Earth: the 200 K minimum of the thermal structure of
Johnstone et al. 2018 over the equilibrium 255 K). A stand-in for the thermal structure of the lower air, tier 5, tweak level 0."""

EXOBASE_RISE_K = 1095.0
"""float: The rise of the exobase temperature above the mesopause at the present XUV flux of the Earth and no CO2, K (fitted to the
neutral exobase temperatures of Fig. 12 of Johnstone et al. 2018 read from the figure: 1310, 1260, 1155, 985, 540 and 290 K for CO2
mixing ratios at 65 km of 4e-6 to 0.4; tier 3 for the figure's reading, tweak level 3)."""

CO2_COOLING_HALF = 0.0126
"""float: The CO2 mixing ratio at which the rise is halved, with `CO2_COOLING_EXPONENT` (the same fit; mole fraction)."""

CO2_COOLING_EXPONENT = 0.70
"""float: The sharpness of the CO2 cooling, the rise going as `1 / (1 + (x / CO2_COOLING_HALF)^p)` (the same fit)."""

XUV_EXPONENT = 0.96
"""float: The rise goes as the XUV flux to this power: the Earth of Johnstone et al. 2018 at 4.9 times the present XUV (their text: their
model 3 Gyr ago corresponds to the 4.9 times of Tian et al. 2008) has a neutral exobase temperature of 4600 K against 1155 K now, which is
`((4600 - 200) / (1155 - 200))^(1 / ln 4.9) = 0.96`. Tier 3, tweak level 3."""

XUV_FLUX_MAX_EARTH = 4.9
"""float: The highest XUV flux, in the Earth's present, the temperature model is taken to (its published range); the flux is held there for the
younger star (see the module text)."""

# --- the exobase and the flux of a species ---------------------------------------------------------------------------------

MESOPAUSE_PRESSURE_PA = 1.0
"""float: The pressure at the bottom of the thermosphere, 1e-2 mbar, where the thermosphere model of Johnstone et al. 2018 starts its own
exobase estimate (their boundary conditions use 200 K and 1e-2 mbar). Tier 3, tweak level 0."""

BATES_LENGTH_SCALE_HEIGHTS = 2.0
"""float: The length over which the temperature rises from the mesopause to the exobase value, in scale heights of the mesopause (Earth:
12 km). Fitted to the exobase altitudes of Johnstone et al. 2018, Fig. 12 and 13 (475 km now, about 2,000 to 2,500 km at 4,600 K; the profile gives
440 and 2,160 km; read from the figures, tier 3), tweak level 3."""

COLLISION_CROSS_SECTION_M2 = 2e-19
"""float: The cross section of the exobase condition, `n sigma H = 1` (Chamberlain 1963: 2e-15 cm2 is the usual value). Tier 3, level 0."""

MAX_EXOBASE_RADII = 10.0
"""float: The exobase is held inside this many body radii (a numerical bound for a hot air over a small body)."""

MAX_EXOBASE_STEPS = 200_000
"""int: The most steps the exobase search takes (GEN.841). A profile that reaches `MAX_EXOBASE_RADII` takes up to about 15,000 steps over
the whole range of bodies swept (radius 1e6 to 7e7 m, gravity 0.5 to 300 m/s2, 2 to 80 amu, mesopause 5 to 400 K, exobase to 5000 K), so
this bound never acts on a real air; it stops a search that cannot end, and the exobase is then taken at the height reached."""

KINETIC_DIAMETER_PM = {
    "H2": 289.0, "He": 260.0, "Ne": 275.0, "H2O": 265.0, "CH4": 380.0, "CO": 376.0, "N2": 364.0, "O2": 346.0, "Ar": 340.0, "CO2": 330.0,
}
"""dict: Kinetic diameter of a molecule, picometres (Breck 1974, recalled, tier 3): the hard-sphere size of the diffusion coefficient."""

DEFAULT_KINETIC_DIAMETER_PM = 360.0
"""float: The diameter of a gas that is not in `KINETIC_DIAMETER_PM` (about the mean of the heavier molecules; tier 6, level 0)."""

# --- the energy-limited and hydrodynamic rates -----------------------------------------------------------------------------

ENERGY_LIMITED_EFFICIENCY = 0.15
"""float: The share of the absorbed XUV that drives the escape (the grid of Kubyshkina et al. 2018 takes 0.15; Shematovich et al. 2014: 0.10 to
0.20). The water-loss history of GEN.639 keeps its own, lower DESIGN value (`tuning.ENERGY_LIMITED_EFFICIENCY`). Tier 2, tweak level 0."""

HYDRO_BASED_LOW = {"beta": 32.0199, "alpha1": 0.4222, "alpha2": -1.7489, "alpha3": 3.7679, "zeta": -6.8618, "theta": 0.0095}
"""dict: Table 1 of Kubyshkina et al. 2018, for the escape parameter below `exp(Sigma)`: boil-off, driven by the planet's heat and low gravity."""

HYDRO_BASED_HIGH = {"beta": 16.4084, "alpha1": 1.0000, "alpha2": -3.2861, "alpha3": 2.7500, "zeta": -1.2978, "theta": 0.8846}
"""dict: Table 1, at or above `exp(Sigma)`: the XUV-driven, near energy-limited branch."""

HYDRO_BASED_SIGMA = (15.611, -0.578, 1.537, 1.018, 5.564, 0.894)
"""tuple: The constants of Eq. 10, `Sigma = (a + b ln F + c ln d + e ln R) / (f + g ln d)` (F in erg/cm2/s, d in AU, R in Earth radii)."""

HYDRO_BASED_MASS_RANGE_EARTH = (1.0, 40.0)
HYDRO_BASED_TEQ_RANGE_K = (300.0, 2000.0)
"""tuple: The grid the fit covers; outside it the rate is an extrapolation."""

_G = constants.G
_K = constants.BOLTZMANN
_AVOGADRO = 6.02214076e23
_M_HYDROGEN_KG = 1.6735575e-27
_M_EARTH_KG = constants.EARTH_MASS_TO_KG
_R_EARTH_M = constants.EARTH_RADIUS_KM * 1000.0
_AU_M = constants.AU_TO_M
_GY_S = 1e9 * constants.SECONDS_PER_YEAR
_ERG_CM2_S_PER_W_M2 = 1e3


def mesopause_temperature_k(teq_k):
    """The temperature at the bottom of the thermosphere, K."""
    return MESOPAUSE_OVER_TEQ * teq_k


def exobase_temperature_k(xuv_flux_earth, mesopause_k, co2_fraction):
    """
    The neutral exobase temperature, K: the mesopause plus the rise, which goes as the XUV flux and is cooled by CO2.

    Args:
        xuv_flux_earth (float): The XUV flux at the body in the Earth's present (`constants.EARTH_XUV_FLUX_W_M2`), held to the model's range.
        mesopause_k (float): `mesopause_temperature_k`, or the observed value of a body whose lower air is not in radiative equilibrium
            with the star (Titan's haze: 150 K).
        co2_fraction (float): The CO2 mole fraction of the air.
    """
    flux = min(max(xuv_flux_earth, 0.0), XUV_FLUX_MAX_EARTH)
    cooling = 1.0 + (max(co2_fraction, 0.0) / CO2_COOLING_HALF) ** CO2_COOLING_EXPONENT
    return mesopause_k + EXOBASE_RISE_K * flux ** XUV_EXPONENT / cooling


@functools.lru_cache(maxsize=None)
def _molar_mass_kg_mol(formula):
    return gases.gas(formula).molar_mass


def mean_mass_kg(fractions):
    """The mean molecular mass of a mixture (mole fractions, any normalisation), kg."""
    total = sum(fractions.values())
    return sum(f * _molar_mass_kg_mol(formula) for formula, f in fractions.items()) / total / _AVOGADRO


def exobase(radius_m, gravity_ms2, mean_mass, temperature_k, mesopause_k):
    """
    The exobase of the thermosphere: where `n sigma H = 1`, for the air rising from `MESOPAUSE_PRESSURE_PA` at the mesopause temperature
    to the exobase temperature on the profile of Bates (1959), `T(z) = T_exo - (T_exo - T_meso) exp(-z / L)` with `L` `BATES_LENGTH_SCALE_HEIGHTS`
    scale heights of the mesopause, in hydrostatic balance under the gravity that falls off as the inverse square of the radius.

    Returns:
        tuple: `(radius of the exobase, m; number density there, m-3; scale height there, m)`.
    """
    if not all(math.isfinite(value) and value > 0.0 for value in (radius_m, gravity_ms2, mean_mass, temperature_k, mesopause_k)):
        raise ValueError(f"exobase needs finite positive inputs, got radius {radius_m}, gravity {gravity_ms2}, mean mass {mean_mass}, "
                         f"temperature {temperature_k}, mesopause {mesopause_k}")
    length = BATES_LENGTH_SCALE_HEIGHTS * _K * mesopause_k / (mean_mass * gravity_ms2)
    log_density = math.log(MESOPAUSE_PRESSURE_PA / (_K * mesopause_k))
    limit = (MAX_EXOBASE_RADII - 1.0) * radius_m
    height = 0.0
    rise = temperature_k - mesopause_k
    for _ in range(MAX_EXOBASE_STEPS):
        decay = math.exp(-height / length)
        temperature = temperature_k - rise * decay
        scale = _K * temperature / (mean_mass * gravity_ms2 * (radius_m / (radius_m + height)) ** 2)
        if math.exp(log_density) * COLLISION_CROSS_SECTION_M2 * scale <= 1.0 or height >= limit:
            break
        step = 0.25 * min(scale, length)
        gradient = rise * decay / length / temperature
        log_density -= (1.0 / scale + gradient) * step
        height += step
    # past the step bound the air is taken to reach the limit, as for a profile that never closes inside MAX_EXOBASE_RADII
    temperature = temperature_k - (temperature_k - mesopause_k) * math.exp(-height / length)
    scale = _K * temperature / (mean_mass * gravity_ms2 * (radius_m / (radius_m + height)) ** 2)
    return radius_m + height, 1.0 / (COLLISION_CROSS_SECTION_M2 * scale), scale


def jeans_parameter(mass_kg, temperature_k, gm_m3_s2, radius_m):
    """`lambda = G M m / (k T r)`: the escape energy over the thermal energy of a molecule of `mass_kg`."""
    return gm_m3_s2 * mass_kg / (_K * temperature_k * radius_m)


def jeans_flux(number_density_m3, mass_kg, temperature_k, gm_m3_s2, radius_m):
    """The Jeans escape flux of a species through the exobase, m-2 s-1: `n v0 (1 + lambda) exp(-lambda) / (2 sqrt(pi))`."""
    parameter = jeans_parameter(mass_kg, temperature_k, gm_m3_s2, radius_m)
    most_probable = math.sqrt(2.0 * _K * temperature_k / mass_kg)
    return number_density_m3 * most_probable * (1.0 + parameter) * math.exp(-parameter) / (2.0 * math.sqrt(math.pi))


def _diameter_m(formula):
    return KINETIC_DIAMETER_PM.get(formula, DEFAULT_KINETIC_DIAMETER_PM) * 1e-12


@functools.lru_cache(maxsize=4096)
def binary_diffusion_parameter(formula_a, formula_b, temperature_k):
    """`b = n D` of two gases, m-1 s-1: the first approximation of Chapman and Enskog for hard spheres, `3 / (8 sigma^2) sqrt(k T / (2 pi mu))`.
    Cached: the same pairs come back at the same mesopause temperature on every step of a loss history."""
    mass_a = _molar_mass_kg_mol(formula_a) / _AVOGADRO
    mass_b = _molar_mass_kg_mol(formula_b) / _AVOGADRO
    reduced = mass_a * mass_b / (mass_a + mass_b)
    sigma = 0.5 * (_diameter_m(formula_a) + _diameter_m(formula_b))
    return 3.0 / (8.0 * sigma ** 2) * math.sqrt(_K * temperature_k / (2.0 * math.pi * reduced))


def diffusion_limited_flux(formula, fractions, gravity_ms2, temperature_k, mean=None):
    """
    The most of a species lighter than the air that the air can carry up, m-2 s-1 (Hunten 1973): `b f g (m_air - m) / (k T)` with the
    binary diffusion parameter `b` of the species in the rest of the air (Blanc's law). `math.inf` for a species that is not lighter.
    `mean` is `mean_mass_kg(fractions)` when the caller has it already (it asks for every species of one air).
    """
    total = sum(fractions.values())
    if mean is None:
        mean = mean_mass_kg(fractions)
    mass = _molar_mass_kg_mol(formula) / _AVOGADRO
    if mass >= mean:
        return math.inf
    share = fractions.get(formula, 0.0) / total
    others = {name: f for name, f in fractions.items() if name != formula and f > 0.0}
    rest = sum(others.values())
    if rest <= 0.0:
        return math.inf
    inverse = sum(f / rest / binary_diffusion_parameter(formula, name, temperature_k) for name, f in others.items())
    return share * gravity_ms2 * (mean - mass) / (_K * temperature_k * inverse)


def energy_limited_mass_loss_kg_s(xuv_flux_w_m2, absorbing_radius_m, mass_kg, efficiency=ENERGY_LIMITED_EFFICIENCY):
    """The most mass the absorbed XUV can lift out of the planet's well, kg/s: `eps pi F R_xuv^3 / (G M)` (Watson et al. 1981; Erkaev et al. 2007)."""
    return efficiency * math.pi * xuv_flux_w_m2 * absorbing_radius_m ** 3 / (_G * mass_kg)


# --- the XUV flux and its history ------------------------------------------------------------------------------------------

def xuv_flux_earth(subject):
    """
    The XUV flux at the body now, in the Earth's present: the stored `xuv_flux_earth` when the subject has it, else the star's corona at its
    age (`activity.coronal_xuv_at`) over the body's distance (from `insolation` and `luminosity_solar`); 0 for a body with no star.
    """
    stored = subject.get("xuv_flux_earth")
    if stored is not None:
        return stored
    mass_sol, age, luminosity, insolation = (subject.get(k) for k in ("host_mass_solar", "age_gy", "luminosity_solar", "insolation"))
    if not (mass_sol and age and luminosity and insolation):
        return 0.0
    luminosity_w = luminosity * constants.SOLAR_LUMINOSITY
    distance_m = _AU_M * math.sqrt(luminosity / insolation)
    return activity.coronal_xuv_at(mass_sol, age, luminosity_w) / (4.0 * math.pi * distance_m ** 2) / constants.EARTH_XUV_FLUX_W_M2


def xuv_ratio(mass_sol, luminosity_w, age_gy, at_gy, log_lx_lbol=None):
    """The star's `L_XUV` at `at_gy` over its `L_XUV` at `age_gy` (the corona of `activity`); with the star's own stored `log_lx_lbol` the decay is the one
    through it (`activity.effective_saturation_age_gy`, GEN.804), the table's nominal one without."""
    if log_lx_lbol is None:
        now = activity.coronal_xuv_at(mass_sol, age_gy, luminosity_w)
        return activity.coronal_xuv_at(mass_sol, at_gy, luminosity_w) / now if now > 0.0 else 1.0
    t_sat = activity.effective_saturation_age_gy(mass_sol, age_gy, log_lx_lbol)
    now = activity.coronal_xuv_with(mass_sol, age_gy, luminosity_w, t_sat)
    return activity.coronal_xuv_with(mass_sol, at_gy, luminosity_w, t_sat) / now if now > 0.0 else 1.0


# --- the share of an air that is kept --------------------------------------------------------------------------------------

HISTORY_STEPS = 48
"""int: Log-spaced steps of the loss history from `tuning.RUNAWAY_HISTORY_START_YR` to the body's age."""


def kept_shares(fractions, mass_kg, radius_m, teq_k, pressure_pa, age_gy, xuv_flux_now, host_mass_sol=None, luminosity_w=None, mesopause_k=None,
                magnetopause_rp=None, log_lx_lbol=None):
    """
    The share of each species of an air that is left at `age_gy`, after thermal escape through the star's XUV history.

    The air is followed through the history: the columns of its species (moles per area, `f P / (m g)`) fall by the flux each species has at the
    time, and the mean mass, the CO2 that cools the thermosphere and the share each species has in the air that can leave are those of the air
    that is left (a hydrogen-rich air drags the heavier gases only while it still has the hydrogen).

    Args:
        fractions (dict): Formula -> mole fraction of the air as it is supplied.
        mass_kg, radius_m (float): The body's.
        teq_k (float): Its equilibrium temperature (the mesopause is `MESOPAUSE_OVER_TEQ` of it, unless `mesopause_k`).
        pressure_pa (float): The surface pressure (the column that the loss is compared to).
        age_gy (float): The body's age.
        xuv_flux_now (float): `xuv_flux_earth` of the body.
        host_mass_sol, luminosity_w (float, optional): The star's, for the shape of its XUV decay; without them the flux is held at its present value.
        mesopause_k (float, optional): The mesopause temperature where the lower air is not in radiative equilibrium with the star.
        magnetopause_rp (float, optional): The magnetopause in planet radii (GEN.804): the stellar wind carries off ions (`ion_escape`) through the share of the
            surface the field leaves open, 0 for a world with no field (open all over). `None` is a field not yet drawn: the ion loss is not followed.
        log_lx_lbol (float, optional): The star's own `log10 L_X/L_bol` at `age_gy`: the XUV history then follows the decay through it (a star that spun
            slowly or fast kept its saturation for a different time) instead of the table's nominal one.

    Returns:
        tuple: `(kept, present_exobase_k)`, `kept` the formula -> share 0..1 and the temperature that of the air as it is supplied.
    """
    total = sum(fractions.values())
    shares = {formula: f / total for formula, f in fractions.items() if f > 0.0}
    masses = {formula: gases.gas(formula).molar_mass / _AVOGADRO for formula in shares}
    gm = _G * mass_kg
    gravity = gm / radius_m ** 2
    meso = mesopause_k if mesopause_k is not None else mesopause_temperature_k(teq_k)
    mean = mean_mass_kg(shares)
    column0 = {formula: share * pressure_pa / (mean * gravity) for formula, share in shares.items()}
    present = exobase_temperature_k(xuv_flux_now, meso, shares.get("CO2", 0.0))
    column = dict(column0)
    start_gy = min(tuning.RUNAWAY_HISTORY_START_YR / 1e9, 0.5 * age_gy)
    ratio_step = (age_gy / start_gy) ** (1.0 / HISTORY_STEPS)
    t0 = start_gy
    limit_key, limits = None, {}
    for _ in range(HISTORY_STEPS):
        t1 = t0 * ratio_step
        alive = {formula: n for formula, n in column.items() if n > 1e-12 * column0[formula]}
        if not alive:
            break
        count = sum(alive.values())
        now = {formula: n / count for formula, n in alive.items()}
        mean_now = mean_mass_kg(now)
        history = xuv_ratio(host_mass_sol, luminosity_w, age_gy, math.sqrt(t0 * t1), log_lx_lbol) if host_mass_sol and luminosity_w else 1.0
        flux = xuv_flux_now * history
        temperature = exobase_temperature_k(flux, meso, now.get("CO2", 0.0))
        r_exo, density, _ = exobase(radius_m, gravity, mean_now, temperature, meso)
        area_ratio = (r_exo / radius_m) ** 2
        key = (tuple(now.items()), mean_now)
        if key != limit_key:                                                      # a heavy gas that loses nothing leaves the air as it was
            limit_key, limits = key, {}
        rate = {}
        for formula, share in now.items():
            jeans = jeans_flux(share * density, masses[formula], temperature, gm, r_exo) * area_ratio
            if jeans > 0.0:                                                       # `min(0, limit)` is 0 for any limit that is not negative
                limit = limits.get(formula)
                if limit is None:
                    limit = limits[formula] = diffusion_limited_flux(formula, now, gravity, meso, mean_now)
                jeans = min(jeans, limit)
            rate[formula] = jeans
        budget = energy_limited_mass_loss_kg_s(min(flux, XUV_FLUX_MAX_EARTH) * constants.EARTH_XUV_FLUX_W_M2, r_exo, mass_kg)
        demanded = 4.0 * math.pi * radius_m ** 2 * sum(rate[f] * masses[f] for f in rate)
        scale = min(1.0, budget / demanded) if demanded > 0.0 else 1.0
        dt = (t1 - t0) * _GY_S
        ion_kg_s = 0.0 if magnetopause_rp is None else ion_escape.ion_loss_kg_s(
            energy_limited_mass_loss_kg_s(min(flux, XUV_FLUX_MAX_EARTH) * constants.EARTH_XUV_FLUX_W_M2, 1.05 * radius_m, mass_kg, efficiency=1.0),
            magnetopause_rp)
        ion_per_m2 = ion_kg_s / (mean_now * 4.0 * math.pi * radius_m ** 2)     # molecules per m2 s, taken from the species in the proportion they are in
        for formula, value in rate.items():
            column[formula] = alive[formula] * math.exp(-min((scale * value + ion_per_m2 * now[formula]) * dt / alive[formula], 700.0))
        t0 = t1
    return {formula: column[formula] / column0[formula] for formula in shares}, present


# --- the hydrodynamic rate of a hydrogen envelope -------------------------------------------------------------------------

def hydro_based_mass_loss_g_s(mass_earth, radius_earth, teq_k, distance_au, xuv_flux_erg_cm2_s):
    """
    The mass-loss rate of a planet with a hydrogen-dominated envelope, g/s (Kubyshkina et al. 2018, Eq. 9 and Table 1): a power law in
    the XUV flux, the distance and the radius times the escape parameter `Lambda = G M m_H / (k T_eq R)` to a power that depends on the
    distance, with one set of constants where `Lambda < exp(Sigma)` (boil-off) and another above it.
    """
    parameter = _G * mass_earth * _M_EARTH_KG * _M_HYDROGEN_KG / (_K * teq_k * radius_earth * _R_EARTH_M)
    a, b, c, e, f, g = HYDRO_BASED_SIGMA
    ln_flux, ln_distance, ln_radius = math.log(xuv_flux_erg_cm2_s), math.log(distance_au), math.log(radius_earth)
    sigma = (a + b * ln_flux + c * ln_distance + e * ln_radius) / (f + g * ln_distance)
    row = HYDRO_BASED_LOW if math.log(parameter) < sigma else HYDRO_BASED_HIGH  # in logs: a starved flux makes exp(sigma) overflow
    power = row["zeta"] + row["theta"] * ln_distance
    return (math.exp(row["beta"]) * xuv_flux_erg_cm2_s ** row["alpha1"] * distance_au ** row["alpha2"] * radius_earth ** row["alpha3"]
            * parameter ** power)
