# planetgen/physics/pre_main_sequence.py

"""
The pre-main-sequence state of a young star (GEN.242): the track from the birthline to the main sequence, the stage a
star of a given mass and age is in (protostar, pre-main-sequence star, or neither) and the disk it still holds.

The track is worked out, not tabulated. A star shrinks while it radiates, so its energy `E = -k G M^2 / R` (`k` is 3/7
for a convective star, a polytrope of index 3/2, and 3/4 for a radiative one, index 3) changes at the rate it
shines: `dE/dt = -L`. On the Hayashi line a convective star keeps its effective temperature `T`, so `L = 4 pi R^2
sigma T^4` and the radius follows `R^-3 = R0^-3 + 3 k4 t` with `k4 = 4 pi sigma T^4 / (k G M^2)`. That phase ends when
the star's luminosity has fallen to its main-sequence value (the radius it reaches is `R_ms (T_ms / T)^2`); from there
to the main-sequence radius it is radiative, shines at the main-sequence luminosity, and `1/R` grows linearly with
time while the temperature rises. The sum of the two is the Kelvin-Helmholtz time, and it replaces the `31 Myr * m^-2.5` law of the multi-star
note (which gives 10 Gyr at 0.1 solar masses, where models give about 1 Gyr, and had to be capped). Contraction does
not pay for all the light, so `G M^2` is stretched by `tuning.PRE_MAIN_SEQUENCE_ENERGY_SOURCE_FACTOR`.

Three inputs (`tuning`, GEN.529) are fitted to the published BHAC15 tracks [R0254]: the
Hayashi temperature and its mass slope, the birthline radius and its slope, and the energy factor above. Past the
table's 1.4 solar masses they are extended by the same power laws.

Two times (GEN.1174, Boss 5:25 AM EDT 2026-10-31: both times, the class by the radius). `radius_settled_yr` is the age at which the radius
has settled to within 5 percent of its main-sequence value in the BHAC15 tracks, the time a star stops contracting; the model's clock
is stretched onto it, the star is pre-main-sequence until then, and its track ends on the main sequence at that age.
`nuclear_arrival_yr` is MIST's arrival (GEN.947: the heaviest star still in phase -1 of each isochrone age,
`isochrones.pre_main_sequence_arrival_yr`, a nuclear criterion that comes later than the radius for the lightest stars), kept beside it;
where MIST's rows do not reach the mass it is the model's own Kelvin-Helmholtz time.
"""

import functools
import math

from planetgen import tuning
from planetgen.datasets import isochrones
from planetgen.physics import constants, main_sequence, unit_registry as units
from planetgen.util import draw

_msun_to_kg = units.converter("Msun", "kg")

RADIUS_SETTLED_MYR = (
    (0.07, 1010.4), (0.08, 830.5), (0.09, 506.6), (0.10, 402.4), (0.11, 355.5), (0.13, 279.2), (0.15, 240.1), (0.17, 216.1),
    (0.20, 187.5), (0.30, 130.2), (0.40, 105.5), (0.50, 88.2), (0.60, 74.8), (0.70, 62.2), (0.80, 50.3),
)
"""tuple: `(initial mass in solar masses, Myr)`: the age at which the radius of the BHAC15 track ([R0254], solar metallicity,
`BHAC15_tracks+structure`) first falls within 5 percent of its value at 5 Gyr, read for GEN.1174 (the report `gen1169-pre-main-sequence-share-2026-10-31.md`
and its script). Above 0.8 solar masses the radius of the track goes through a minimum and the criterion stops meaning the same thing, so the model's
own Kelvin-Helmholtz time is carried on at the 0.8 ratio (`radius_settled_yr`); registry row LM6 of 6.66."""

PROTOSTAR = "protostar"
PRE_MAIN_SEQUENCE = "PMS"
"""str: The `yerkes_class` of a protostar and of a pre-main-sequence star (Class II and III)."""
CLASS_NAMES = {PROTOSTAR: "Protostar", PRE_MAIN_SEQUENCE: "Pre-Main-Sequence Star"}
"""dict: Display name of each class, in the same place `YERKES_CLASS_NAMES` is used."""
STAGE_CLASSES = (PROTOSTAR, PRE_MAIN_SEQUENCE)
DISK_UNITS = {"disk_kind": "text", "disk_mass_solar": "Msun", "disk_inner_au": "au", "disk_outer_au": "au"}
"""dict: The disk descriptor's columns on `stars`, in `draw_disk`'s key order, with each attribute's native unit."""
DISK_FIELDS = tuple(DISK_UNITS)
"""tuple: The names of `DISK_UNITS`."""
_yr_to_s = units.converter("yr", "s")
_s_to_yr = units.converter("s", "yr")
_m_to_au = units.converter("m", "au")


def hayashi_temperature_k(mass_sol):
    """The effective temperature of a convective star of this mass on its Hayashi line (fitted to BHAC15, `tuning`)."""
    return tuning.HAYASHI_TEMPERATURE_1_SOL_K * mass_sol ** tuning.HAYASHI_TEMPERATURE_MASS_EXPONENT


def birthline_radius_m(mass_sol):
    """The radius a star has when its accretion ends and it starts to contract, metres (fitted to BHAC15, `tuning`)."""
    return (tuning.BIRTHLINE_RADIUS_1_SOL_RSUN * mass_sol ** tuning.BIRTHLINE_RADIUS_MASS_EXPONENT
            * constants.SOLAR_RADIUS_M)


def _binding(mass_sol):
    """`G M^2`, stretched by `PRE_MAIN_SEQUENCE_ENERGY_SOURCE_FACTOR` (contraction pays for only part of the light)."""
    return tuning.PRE_MAIN_SEQUENCE_ENERGY_SOURCE_FACTOR * constants.G * (_msun_to_kg(mass_sol)) ** 2


@functools.lru_cache(maxsize=4096)
def _track_constants(mass_sol, fe_h=0.0):
    """
    The constants of the contraction before any stretch: `r0` (birthline radius), `r_h` (where the Hayashi phase ends), `r_ms`, the
    temperature `t_h` of the Hayashi phase, the main-sequence luminosity `l_ms` (W), the two phase lengths `t1`
    and `t2` (s) and the rate `k4` of the convective phase.
    """
    arrival = main_sequence.zams_state(mass_sol, fe_h)     # GEN.652: the MIST row the star contracts onto
    r_ms = arrival.radius_sol * constants.SOLAR_RADIUS_M
    l_ms = arrival.luminosity_sol * constants.SOLAR_LUMINOSITY
    t_ms = arrival.temperature_k
    t_h = min(hayashi_temperature_k(mass_sol), t_ms)
    r_h = r_ms * (t_ms / t_h) ** 2
    r0 = birthline_radius_m(mass_sol)
    binding = _binding(mass_sol)
    k4 = (4.0 * math.pi * constants.STEFAN_BOLTZMANN_CONSTANT * t_h ** 4
          / (tuning.CONVECTIVE_BINDING_COEFFICIENT * binding))
    t1 = max(0.0, (r_h ** -3 - r0 ** -3) / (3.0 * k4)) if r0 > r_h else 0.0
    start = min(r0, r_h)
    t2 = tuning.RADIATIVE_BINDING_COEFFICIENT * binding * (1.0 / r_ms - 1.0 / start) / l_ms
    t2 = max(t2, 0.0)
    return {"r0": r0, "r_h": r_h, "r_ms": r_ms, "t_h": t_h, "l_ms": l_ms, "t1": t1, "t2": t2, "k4": k4, "start": start}


@functools.lru_cache(maxsize=4096)
def _track(mass_sol, fe_h=0.0):
    """The constants of `_track_constants` and `scale`, the factor that stretches the model's time onto the radius-settled time (GEN.1174)."""
    track = dict(_track_constants(mass_sol, fe_h))
    own_yr = _s_to_yr(track["t1"] + track["t2"])
    track["scale"] = radius_settled_yr(mass_sol, fe_h, own_yr=own_yr) / own_yr if own_yr > 0.0 else 1.0
    return track


def radius_settled_yr(mass_sol, fe_h=0.0, own_yr=None):
    """
    The age at which a star of this initial mass has contracted to the main sequence, years (GEN.1174): the BHAC15 age at which its radius settles
    to within 5 percent (`RADIUS_SETTLED_MYR`, log-log between the masses, flat below 0.07 solar masses); above 0.8 solar masses the model's own
    Kelvin-Helmholtz time (`own_yr`, worked out when not given) scaled by the ratio the two have at 0.8. The tracks are for solar metallicity; the
    time is taken for every [Fe/H] (`fe_h` is only here for the callers' sake).
    """
    masses = [m for m, _ in RADIUS_SETTLED_MYR]
    if mass_sol <= masses[-1]:
        ln_m = math.log(max(mass_sol, masses[0]))
        for (m0, t0), (m1, t1) in zip(RADIUS_SETTLED_MYR, RADIUS_SETTLED_MYR[1:]):
            if mass_sol <= m1:
                share = (ln_m - math.log(m0)) / (math.log(m1) - math.log(m0))
                return 1e6 * math.exp(math.log(t0) + share * (math.log(t1) - math.log(t0)))
    edge = _track_own_yr(masses[-1], fe_h)
    if own_yr is None:
        own_yr = _track_own_yr(mass_sol, fe_h)
    return RADIUS_SETTLED_MYR[-1][1] * 1e6 / edge * own_yr


def _track_own_yr(mass_sol, fe_h):
    """The model's own Kelvin-Helmholtz time, years, before any stretch."""
    track = _track_constants(mass_sol, fe_h)
    return _s_to_yr(track["t1"] + track["t2"])


def nuclear_arrival_yr(mass_sol, fe_h=0.0):
    """
    The age at which a star of this initial mass reaches MIST's main sequence, years (GEN.947): the age at which the isochrone rows hold no contracting star
    of it, where they reach the mass, else the Kelvin-Helmholtz time of the model. A nuclear criterion; it comes later than `radius_settled_yr` for
    the lightest stars (0.73 against 0.19 Gyr at 0.2 solar masses) and the star is called pre-main-sequence by the radius, not by this.
    """
    mist_yr = isochrones.pre_main_sequence_arrival_yr(mass_sol, fe_h)
    return mist_yr if mist_yr is not None else _track_own_yr(mass_sol, fe_h)


def contraction_time_yr(mass_sol, fe_h=0.0):
    """The age at which a star of this mass has contracted to the main sequence, years: when its radius has settled (`radius_settled_yr`, GEN.1174)."""
    track = _track(mass_sol, fe_h)
    return _s_to_yr((track["t1"] + track["t2"]) * track["scale"])


def stage(mass_sol, age_yr, fe_h=0.0):
    """
    `PROTOSTAR`, `PRE_MAIN_SEQUENCE` or None for a star of this initial mass at this age.

    A star of `PRE_MAIN_SEQUENCE_MAX_MASS_SOL` or more is on the main sequence by the time it shows, so it has none.
    Under `PROTOSTAR_MAX_AGE_YR` (Class 0 and I) a star still contracting is a protostar; after that, and until it
    has contracted to the main sequence however long that takes (about a billion years at 0.1 solar masses), a
    pre-main-sequence star.
    """
    if (mass_sol >= tuning.PRE_MAIN_SEQUENCE_MAX_MASS_SOL or age_yr >= tuning.PRE_MAIN_SEQUENCE_LONGEST_YR
            or age_yr >= contraction_time_yr(mass_sol, fe_h)):
        return None
    return PROTOSTAR if age_yr < tuning.PROTOSTAR_MAX_AGE_YR else PRE_MAIN_SEQUENCE


def _track_state(mass_sol, age_yr, fe_h=0.0):
    """The luminosity (solar), temperature and radius (solar) of a star on the track at this age, any age up to its
    contraction time, with no regard to the stage limits of `stage`."""
    track = _track(mass_sol, fe_h)
    t = _yr_to_s(age_yr) / track["scale"]   # the model's time, stretched onto MIST's arrival
    if t <= track["t1"]:
        radius = (track["r0"] ** -3 + 3.0 * track["k4"] * t) ** (-1.0 / 3.0)
        temperature = track["t_h"]
        luminosity = 4.0 * math.pi * radius ** 2 * constants.STEFAN_BOLTZMANN_CONSTANT * temperature ** 4
    else:
        radius = 1.0 / (1.0 / track["start"] + track["l_ms"] * (t - track["t1"])
                        / (tuning.RADIATIVE_BINDING_COEFFICIENT * _binding(mass_sol)))
        luminosity = track["l_ms"]
        temperature = (luminosity / (4.0 * math.pi * radius ** 2 * constants.STEFAN_BOLTZMANN_CONSTANT)) ** 0.25
    return {"luminosity_sol": luminosity / constants.SOLAR_LUMINOSITY, "temperature_k": temperature,
            "radius_sol": radius / constants.SOLAR_RADIUS_M}


def state(mass_sol, age_yr, fe_h=0.0):
    """
    The luminosity, temperature and radius of a contracting star, or None when `stage` is None.

    Returns:
        dict: `yerkes_class` (`PROTOSTAR` or `PRE_MAIN_SEQUENCE`), `luminosity_sol`, `temperature_k`, `radius_km`
            (None: it follows from the first two, like every other star's), `phase_end_gy` (the age it leaves the
            state: the radius-settled time), `nuclear_arrival_gy` (MIST's, GEN.1174), and `disk_class` ("0",
            "I", "II" or "III").
    """
    which = stage(mass_sol, age_yr, fe_h)
    if which is None:
        return None
    on_track = _track_state(mass_sol, age_yr, fe_h)
    return {"yerkes_class": which, "luminosity_sol": on_track["luminosity_sol"],
            "temperature_k": on_track["temperature_k"], "radius_km": None,
            "phase_end_gy": contraction_time_yr(mass_sol, fe_h) / 1.0e9,
            "nuclear_arrival_gy": nuclear_arrival_yr(mass_sol, fe_h) / 1.0e9,
            "disk_class": disk_class(age_yr)}


def disk_class(age_yr):
    """The infrared class of a young star by age: "0", "I", "II", and "III" once its disk is gone."""
    for name, limit in tuning.DISK_CLASS_AGES_YR:
        if age_yr < limit:
            return name
    return "III"


def sublimation_radius_au(luminosity_sol):
    """Where a blackbody grain reaches `DISK_DUST_SUBLIMATION_K` in the star's light, AU: the disk's inner edge."""
    metres = math.sqrt(luminosity_sol * constants.SOLAR_LUMINOSITY
                       / (16.0 * math.pi * constants.STEFAN_BOLTZMANN_CONSTANT * tuning.DISK_DUST_SUBLIMATION_K ** 4))
    return _m_to_au(metres)


def _log_uniform(low, high, rng):
    return math.exp(rng.uniform(math.log(low), math.log(high)))


def draw_disk(mass_sol, age_yr, luminosity_sol, rng=draw):
    """
    The disk of a young star: `disk_kind` ("class-0", "class-i" or "class-ii"), `disk_mass_solar` (a log-uniform
    share of the star's mass, `DISK_MASS_FRACTION_RANGE`), `disk_inner_au` (the dust sublimation radius) and
    `disk_outer_au` (log-uniform over `DISK_OUTER_RADIUS_RANGE_AU`), with no bodies in it; None once the disk is gone
    (Class III), and for a star that is no longer contracting.
    """
    name = disk_class(age_yr)
    if name == "III" or stage(mass_sol, age_yr) is None:
        return None
    return {"disk_kind": f"class-{name.lower()}",
            "disk_mass_solar": mass_sol * _log_uniform(*tuning.DISK_MASS_FRACTION_RANGE, rng),
            "disk_inner_au": sublimation_radius_au(luminosity_sol),
            "disk_outer_au": _log_uniform(*tuning.DISK_OUTER_RADIUS_RANGE_AU, rng)}


def companion_chance(mass_sol, age_yr):
    """
    The share of protostars that have a companion, by their Class (Perseus' [R4379], `tuning`): the
    Class 0 and Class I values; None for a star that is not a protostar, which keeps the field's multiplicity.
    """
    if stage(mass_sol, age_yr) != PROTOSTAR:
        return None
    return tuning.PROTOSTAR_MULTIPLICITY_BY_CLASS.get(disk_class(age_yr))
