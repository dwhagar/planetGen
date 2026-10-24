# planetgen/physics/runaway_history.py

"""
The water a planet or moon lost while its star's light was past the runaway greenhouse limit (GEN.639; design note 4.1 of
docs/design/atmospheres-retention-and-classes.md).

A planet does not lose its oceans when the flux it gets now crosses `orbits.runaway_greenhouse_flux`, it loses them whenever it
did, and a cool star is brightest, and its planets' skies bluest in XUV, when it is young: a 0.3 solar-mass star is 10 to 100
times more luminous on the pre-main-sequence track than it will ever be on the main sequence, and its habitable zone of the
main sequence lies deep inside that early runaway limit. Every body's flux and the star's light are followed from the end of
the class II disk (`tuning.RUNAWAY_HISTORY_START_YR`) to the star's age: the pre-main-sequence track (`pre_main_sequence`),
then a log-linear brightening from the luminosity at the main sequence to the star's present one. Wherever the flux is above
the planet-mass runaway limit the hydrogen of the steam escapes at the energy-limited rate

    Mdot_H = eps pi R^3 F_xuv / (G M),   eps = `tuning.ENERGY_LIMITED_EFFICIENCY`

with the star's XUV (the coronal part from `activity`, the photospheric part above 13.6 eV) at that age and distance, and the
water lost is `tuning.WATER_PER_HYDROGEN_MASS` times the hydrogen. Nothing is drawn: the loss follows from the star, the orbit
and the planet, so a body of the same inputs always loses the same water. The star's own activity scatter enters through the
ratio of its stored present `L_X/L_bol` to the nominal one (`log_lx_lbol`).
"""

import math

from planetgen import tuning
from planetgen.physics import activity, constants, orbits, pre_main_sequence

COOLEST_RUNAWAY_FLUX = min(orbits.runaway_greenhouse_flux(teff, mass)
                           for teff in (orbits.KOPPARAPU_TEFF_RANGE_K[0], 4000.0, 5780.0, orbits.KOPPARAPU_TEFF_RANGE_K[1])
                           for mass in orbits.KOPPARAPU_MASS_RANGE_EARTH)
"""float: The lowest runaway limit of any star and planet mass the three sets span, in the Sun's flux at the Earth: a
body that never got more than this has nothing to follow."""

HISTORY_YERKES = ("V", "VI", "IV", "PMS")
"""tuple: The stars whose past `loss_kg` follows: those still on, or just off, the main sequence. A giant's, white
dwarf's or compact remnant's earlier life is not followed."""

_YEAR_S = constants.SECONDS_PER_YEAR


def _light(mass_sol, age_yr, now):
    """`(luminosity_w, temperature_k)` of the star at `age_yr`: the contraction track until it reaches the main sequence,
    then a log-linear rise to the present values in `now = (age_yr, luminosity_w, temperature_k)`."""
    arrival = pre_main_sequence.contraction_time_yr(mass_sol)
    if age_yr <= arrival:
        on_track = pre_main_sequence._track_state(mass_sol, age_yr)
        return on_track["luminosity_sol"] * constants.SOLAR_LUMINOSITY, on_track["temperature_k"]
    now_age, now_l, now_t = now
    if now_age <= arrival:
        return now_l, now_t
    start = pre_main_sequence._track_state(mass_sol, arrival)
    share = math.log(age_yr / arrival) / math.log(now_age / arrival)
    l0 = start["luminosity_sol"] * constants.SOLAR_LUMINOSITY
    t0 = start["temperature_k"]
    return (math.exp(math.log(l0) + share * (math.log(now_l) - math.log(l0))),
            math.exp(math.log(t0) + share * (math.log(now_t) - math.log(t0))))


def _xuv_w(mass_sol, age_yr, luminosity_w, temperature_k, activity_scale):
    """The star's `L_XUV` at this age: coronal (saturated, then decaying) plus the photosphere above 13.6 eV."""
    t_sat, decay = activity._interp(activity.ACTIVITY_TABLE, mass_sol)
    ratio = activity.lx_lbol_at(age_yr / 1e9, activity.saturated_lx_lbol(mass_sol), t_sat, decay)
    if age_yr / 1e9 > t_sat:
        ratio *= activity_scale
    photosphere = luminosity_w * activity.blackbody_fraction_above(activity.XUV_EDGE_EV, temperature_k)
    return activity.coronal_xuv_w(ratio * luminosity_w) + photosphere


def _activity_scale(star, mass_sol, age_yr):
    """How much of the nominal decaying `L_X/L_bol` this star shows now (its stored `log_lx_lbol` against the nominal
    value at its age), 1 when the star keeps none or is still saturated. Held to a factor of 10 either way."""
    stored = getattr(star, "log_lx_lbol", None)
    if stored is None:
        return 1.0
    t_sat, decay = activity._interp(activity.ACTIVITY_TABLE, mass_sol)
    age_gy = age_yr / 1e9
    if age_gy <= t_sat:
        return 1.0
    nominal = activity.lx_lbol_at(age_gy, activity.saturated_lx_lbol(mass_sol), t_sat, decay)
    return min(10.0, max(0.1, 10 ** stored / nominal))


def loss_kg(star, distance_au, mass_kg, radius_km):
    """
    The water, kg, a body of `mass_kg` and `radius_km` at `distance_au` from `star` lost to the runaway greenhouse since the
    end of its disk. 0 for a star whose past is not followed (`HISTORY_YERKES`), a star younger than the start of the
    history, or a body that never reached the runaway flux.
    """
    if (getattr(star, "yerkes_class", None) not in HISTORY_YERKES or not star.luminosity or not star.temperature
            or not star.mass or not distance_au or distance_au <= 0 or not mass_kg or not radius_km):
        return 0.0
    mass_sol = star.mass / constants.SOLAR_MASS_TO_KG
    if mass_sol >= tuning.PRE_MAIN_SEQUENCE_MAX_MASS_SOL:
        return 0.0
    age_yr = (star.age or 0.0) * 1e9
    start_yr = tuning.RUNAWAY_HISTORY_START_YR
    if age_yr <= start_yr:
        return 0.0
    now = (age_yr, star.luminosity, star.temperature)
    solar_l = constants.SOLAR_LUMINOSITY
    peak = max(_light(mass_sol, start_yr, now)[0], star.luminosity)
    if peak / solar_l / distance_au ** 2 < COOLEST_RUNAWAY_FLUX:
        return 0.0
    planet_mass_earth = mass_kg / constants.EARTH_MASS_TO_KG
    scale = _activity_scale(star, mass_sol, age_yr)
    area_factor = (tuning.ENERGY_LIMITED_EFFICIENCY * math.pi * (radius_km * 1000.0) ** 3
                   / (constants.G * mass_kg))
    distance_m = distance_au * constants.AU_TO_M
    steps = tuning.RUNAWAY_HISTORY_STEPS
    ratio = (age_yr / start_yr) ** (1.0 / steps)
    t0 = start_yr
    hydrogen_kg = 0.0
    for _ in range(steps):
        t1 = t0 * ratio
        middle = math.sqrt(t0 * t1)
        luminosity, temperature = _light(mass_sol, middle, now)
        if luminosity / solar_l / distance_au ** 2 > orbits.runaway_greenhouse_flux(temperature, planet_mass_earth):
            flux = _xuv_w(mass_sol, middle, luminosity, temperature, scale) / (4.0 * math.pi * distance_m ** 2)
            hydrogen_kg += area_factor * flux * (t1 - t0) * _YEAR_S
        t0 = t1
    return hydrogen_kg * tuning.WATER_PER_HYDROGEN_MASS
