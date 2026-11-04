# planetgen/physics/super_puff.py

"""
A super-puff: a low-mass planet with a huge radius, and how it fades with age (GEN.1042, W15 of docs/design/systems/weird-worlds.md, registry 6.76)
=====================================================================================================================================================

Kepler-51 b, c and d (below 0.1 g/cc, 7 to 10 Earth radii, 2 to 6 Earth masses), Kepler-79 d and Kepler-87 c are planets of a few Earth masses whose
radius is that of a Neptune or a Saturn. Three papers say what makes one and what the module takes from each:

* [R2609]: a core of a few Earth masses can hold a thick envelope, a gas-to-core ratio of 7 to 28 percent for
  Kepler-51 b and over 20 percent for the class, when it accretes a dust-free atmosphere in the outer disk and migrates in. This gives the
  *envelope*, `PUFF_ENVELOPE`, and the class (`MAX_MASS_EARTH`, the periods, "an order of magnitude more rare" than smaller planets).
* [R2741] (`bulk_structure.envelope_radius_earth`): the radius of a core under an envelope of that fraction. It is the radius of
  the module everywhere it holds (1 to 30 Earth masses, envelope up to a half): Kepler-51 d comes out at 9.3 Earth radii with a 41 percent
  envelope, Kepler-51 b at 7.1 with 16 percent. Past it (`covers`) the radius is the giant tables' (`giant_radius`), which give WASP-107 b
  (30.5 Earth masses, 0.96 Jupiter radii, 0.134 g/cc, [R3516]) its radius with a 3.7 Earth-mass core.
* [R3079]: obliquity tides in a close planet put heat in the core and inflate the radius by 10 to 100 percent, and
  by up to a factor 2. The inflation ratio is her fit, equation 14, which holds from 20 Earth fluxes up (Kepler-51 c and d are colder and get
  none: "their irradiation fluxes are too small for the parameter domain of our models").

The fade is the sum of two things the repository already has: the age term of the envelope law (`(age / 5 Gyr)^-0.18`, a puffy young planet
contracts) and the loss of the envelope to the star's XUV (`envelope_escape.envelope_left_earth`, which starts at 1 Gyr). The paper gives no
fading law of its own, so there is none here.

Use:

    envelope = envelope_for_radius(mass_earth, radius_earth, flux_earth, age_gy)
    radius = radius_earth(mass_earth, envelope, flux_earth, age_gy)
    later = radius_at_age_earth(mass_earth, envelope, flux_earth, host_mass_solar, luminosity_solar, 5.0)
    is_super_puff(mass_earth, radius_earth_, age_gy, period_d)
"""

import math

from planetgen.physics import bulk_structure, constants, envelope_escape, giant_radius, unit_registry as units

_mearth_to_kg = units.converter("Mearth", "kg")

MAX_MASS_EARTH = 10.0
"""float: A super-puff weighs under this (the design table; Lee and Chiang's class is 2 to 6 Earth masses, row PF1)."""

MAX_AGE_GY = 1.0
"""float: A planet is a super-puff while it is young (the design table: under 1 Gyr, row PF1). Older ones have contracted and are Neptunes."""

PERIOD_RANGE_D = (10.0, 100.0)
"""tuple: The orbital periods a super-puff has, days (the design table, row PF1; Lee and Chiang have it under 50 days)."""

DENSITY_RANGE_G_CM3 = (0.03, 0.1)
"""tuple: The bulk density of a super-puff, g/cc (the design table; Kepler-51 d is 0.038, row PF1). The upper edge is the definition."""

PUFF_ENVELOPE = (0.07, 0.45)
"""tuple: The envelope fractions of a super-puff: Kepler-51 b 7 to 28 percent (Lee and Chiang's quote of Lopez and Fortney), Kepler-79 d 33 to 40,
Kepler-51 d over 30 and the module's solution 41 percent (row PF2)."""

OCCURRENCE_SHARE = 0.005
"""float: The share of sub-Neptunes that are super-puffs (the design table says under 1 percent; Lee and Chiang bound planets over 4 Earth radii
to under a tenth of the smaller ones inside 50 days; row PF3, a stand-in: no survey statistic was read)."""

LOPEZ_MAX_MASS_EARTH = 30.0
"""float: The heaviest planet the law of Lopez and Fortney is used for here (the range of `envelope_escape`)."""

TIDAL_MIN_FLUX_EARTH = 20.0
"""float: The flux from which Millholland's fit holds, Earth fluxes (Sec. 4.1: "F/F_earth >= 20", row PF4)."""

TIDAL_MIN_LUMINOSITY_RATIO = 1.0e-5
"""float: The tidal luminosity over the irradiation over which the radius inflates (Sec. 4.1, row PF4)."""

TIDAL_MAX_RATIO = 2.0
"""float: The largest inflation: "it caps at about 2" (Sec. 4.1, row PF4)."""

_TIDAL_FIT = (0.29, -0.90, 0.91, 3.11, -0.88)
"""tuple: [R3079] equation 14: `R(tides) / R(no tides) = 1 + a (M / 10 Me)^b (f_env / 0.05)^c (F / 100 F_earth)^d (Q / 1e5)^e`, with
`Q = Q' (1 + cos^2 eps) / sin^2 eps`."""

_R_EARTH_M = constants.EARTH_RADIUS_KM * 1.0e3
_INSOLATION_W_M2 = 1361.0


def covers(mass_earth, envelope):
    """Whether the law of Lopez and Fortney holds: 1 to 30 Earth masses and an envelope up to `envelope_escape.MAX_ENVELOPE`."""
    return 1.0 <= mass_earth <= LOPEZ_MAX_MASS_EARTH and 0.0 < envelope <= envelope_escape.MAX_ENVELOPE


def density_g_cm3(mass_earth, radius_earth_):
    """The bulk density, g/cc."""
    volume = 4.0 / 3.0 * math.pi * (radius_earth_ * _R_EARTH_M) ** 3
    return _mearth_to_kg(mass_earth) / volume / 1000.0


def tidal_luminosity_w(host_mass_kg, radius_m, distance_m, q_term):
    """The tidal power in the planet's core with its orbit circular, W ([R3079] eq. 11): `(9 n / 2) (G M*^2 / R) (R / a)^6 / Q`
    with `Q = Q' (1 + cos^2 eps) / sin^2 eps` (`q_term`, about 1e3 to 1e5 for the planets of her tables)."""
    n = math.sqrt(constants.G * host_mass_kg / distance_m ** 3)
    return 4.5 * n / q_term * constants.G * host_mass_kg ** 2 / radius_m * (radius_m / distance_m) ** 6


def tidal_luminosity_ratio(host_mass_kg, radius_m, distance_m, flux_earth, q_term):
    """The tidal power over the starlight the planet intercepts (`pi R^2 F`)."""
    return tidal_luminosity_w(host_mass_kg, radius_m, distance_m, q_term) / (math.pi * radius_m ** 2 * flux_earth * _INSOLATION_W_M2)


def tidal_inflation_ratio(mass_earth, envelope, flux_earth, q_term):
    """
    The radius with the obliquity tides over the radius without them ([R3079] eq. 14), from `TIDAL_MIN_FLUX_EARTH` up, held at
    `TIDAL_MAX_RATIO`; 1 for a colder planet. Her fit is of her own radii, here the ratio multiplies the radius of Lopez and Fortney (row PF4).
    """
    if flux_earth < TIDAL_MIN_FLUX_EARTH:
        return 1.0
    a, b, c, d, e = _TIDAL_FIT
    ratio = 1.0 + a * (mass_earth / 10.0) ** b * (envelope / 0.05) ** c * (flux_earth / 100.0) ** d * (q_term / 1.0e5) ** e
    return min(ratio, TIDAL_MAX_RATIO)


def radius_earth(mass_earth, envelope, flux_earth, age_gy, tidal_ratio=1.0):
    """
    The radius, Earth radii, of a planet of `mass_earth` with an H/He envelope of `envelope` of its mass: the law of Lopez and Fortney times
    `tidal_ratio` where `covers` holds, the giant tables (`giant_radius`, core `mass (1 - envelope)`) past it.
    """
    if covers(mass_earth, envelope):
        return bulk_structure.envelope_radius_earth(mass_earth, envelope, flux_earth, age_gy) * tidal_ratio
    return giant_radius.radius_km(mass_earth, mass_earth * (1.0 - envelope), flux_earth, age_gy) / constants.EARTH_RADIUS_KM


def envelope_for_radius(mass_earth, target_radius_earth, flux_earth, age_gy):
    """The envelope fraction the law of Lopez and Fortney needs for the radius, by bisection; `None` when the radius lies outside what 0.5 of
    the mass can give or the planet is not `covers`."""
    low, high = 0.002, envelope_escape.MAX_ENVELOPE
    if not covers(mass_earth, high):
        return None
    if not (bulk_structure.envelope_radius_earth(mass_earth, low, flux_earth, age_gy) <= target_radius_earth
            <= bulk_structure.envelope_radius_earth(mass_earth, high, flux_earth, age_gy)):
        return None
    for _ in range(60):
        middle = 0.5 * (low + high)
        if bulk_structure.envelope_radius_earth(mass_earth, middle, flux_earth, age_gy) < target_radius_earth:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def radius_at_age_earth(mass_earth, envelope, flux_earth, host_mass_solar, luminosity_solar, age_gy, tidal_ratio=1.0):
    """
    The radius at `age_gy` of a planet that had `envelope` of its `mass_earth` at the start, Earth radii: the envelope the XUV has left it
    (`envelope_escape.envelope_left_earth`, none before 1 Gyr) under the age term of the envelope law. A planet that has lost its envelope
    is its bare core.
    """
    core = mass_earth * (1.0 - envelope)
    left = envelope_escape.envelope_left_earth(core, envelope, flux_earth, host_mass_solar, luminosity_solar, age_gy)
    if left <= 0.0:
        return bulk_structure.core_radius_earth(core)
    mass = core + left
    return radius_earth(mass, left / mass, flux_earth, age_gy, tidal_ratio)


def is_super_puff(mass_earth, radius_earth_, age_gy, period_d):
    """Whether a planet is a super-puff: under `MAX_MASS_EARTH`, younger than `MAX_AGE_GY`, in `PERIOD_RANGE_D` and of a density under the top of
    `DENSITY_RANGE_G_CM3`."""
    return (mass_earth <= MAX_MASS_EARTH and age_gy <= MAX_AGE_GY and PERIOD_RANGE_D[0] <= period_d <= PERIOD_RANGE_D[1]
            and density_g_cm3(mass_earth, radius_earth_) <= DENSITY_RANGE_G_CM3[1])
