# planetgen/physics/super_puff.py

"""
A super-puff: a low-mass planet with a huge radius, and how it fades with age (GEN.1042, W15 of docs/design/systems/weird-worlds.md, registry 6.76)
=====================================================================================================================================================

Kepler-51 b, c and d (below 0.1 g/cc, 7 to 10 Earth radii, 2 to 6 Earth masses), Kepler-79 d and Kepler-87 c are planets of a few Earth masses whose
radius is that of a Neptune or a Saturn. Three papers say what makes one and what the module takes from each:

* [R2609]: a core of a few Earth masses can hold a thick envelope, a gas-to-core ratio of 7 to 28 percent for
  Kepler-51 b and over 20 percent for the class, when it accretes a dust-free atmosphere in the outer disk and migrates in. This gives the
  *envelope*, `PUFF_ENVELOPE`, and the class (`START_MAX_MASS_EARTH`, `START_MAX_PERIOD_D`, "an order of magnitude more rare" than smaller planets).
* [R2741] (`bulk_structure.envelope_radius_earth`): the radius of a core under an envelope of that fraction. It is the radius of
  the module everywhere it holds (1 to 30 Earth masses, envelope up to a half): Kepler-51 d comes out at 9.3 Earth radii with a 41 percent
  envelope, Kepler-51 b at 7.1 with 16 percent. Past it (`covers`) the radius is the giant tables' (`giant_radius`), which give WASP-107 b
  (30.5 Earth masses, 0.96 Jupiter radii, 0.134 g/cc, [R3516]) its radius with a 3.7 Earth-mass core.
* [R3079]: obliquity tides in a close planet put heat in the core and inflate the radius by 10 to 100 percent, and
  by up to a factor 2. The inflation ratio is her fit, equation 14, which holds from 20 Earth fluxes up (Kepler-51 c and d are colder and get
  none: "their irradiation fluxes are too small for the parameter domain of our models").

The fade is the sum of two things the repository already has: the age term of the envelope law (`(age / 5 Gyr)^-0.18`, a puffy young planet
contracts) and the loss of the envelope to the star's XUV (`envelope_escape.envelope_left_earth`, which starts at 1 Gyr). The paper gives no
fading law of its own, so there is none here. What keeps a super-puff inflated past 1 Gyr (Kepler-79 d, 1.9 to 3.4 Gyr; Kepler-27 b and 31 c, 4.7 Gyr) is the obliquity
tide: `envelope_capture.kept` multiplies the radius by `tidal_inflation_ratio` for a planet of the settle pass, whose tilt and tidal quality are known (GEN.1347). A planet with
no tilt (the first pass, before its spin) gets no tide. The tilt is the planet's own: a planet the tides of its star have brought to a synchronous or upright state has a tilt
near zero and no obliquity tide, so the old super-puffs are the ones that kept a tilt.

Use:

    radius = radius_earth(mass_earth, envelope, flux_earth, age_gy)
    is_super_puff(mass_earth, radius_earth_)
"""

import math

from planetgen.physics import bulk_structure, constants, envelope_escape, giant_radius, unit_registry as units

_mearth_to_kg = units.converter("Mearth", "kg")

START_MAX_MASS_EARTH = 10.0
"""float: The heaviest core that can START as a super-puff (`envelope_capture.puff_eligible`, GEN.1343; the design table, Lee and Chiang's class is 2 to 6 Earth masses, row PF1)."""

START_MAX_PERIOD_D = 100.0
"""float: The longest period of a core that can START as a super-puff, days (`envelope_capture.puff_eligible`; the Kepler sub-Neptunes the occurrence share counts, row PF1). It gates the
start only: a planet is a super-puff by its mass and density (`is_super_puff`), at any period (Kepler-51 d is at 130 days, Kepler-87 c at 191)."""

CLASS_MAX_MASS_EARTH = 30.0
"""float: The heaviest planet that is a super-puff (GEN.1470, row PF1): the top of the range of the radius law (`LOPEZ_MAX_MASS_EARTH`), which is also the sub-Saturn mass of the
super-puffs of the survey catalogues ([R5217], abstract read: "sub-Saturn mass planets with extremely low bulk densities"; no number is in the
abstract, so the 30 is the law's). A planet of 10 to 30 Earth masses under the density edge is a puff the start gate never made but the tide or the young age can."""

DENSITY_RANGE_G_CM3 = (0.03, 0.1)
"""tuple: The bulk density of a super-puff, g/cc (the design table; Kepler-51 d is 0.038, row PF1). The upper edge is the definition."""

PUFF_ENVELOPE = (0.07, 0.45)
"""tuple: The envelope fractions of a super-puff: Kepler-51 b 7 to 28 percent (Lee and Chiang's quote of Lopez and Fortney), Kepler-79 d 33 to 40,
Kepler-51 d over 30 and the module's solution 41 percent (row PF2)."""

OCCURRENCE_SHARE = 0.005
"""float: The share of sub-Neptunes that are super-puffs (the design table says under 1 percent; Lee and Chiang bound planets over 4 Earth radii
to under a tenth of the smaller ones inside 50 days; row PF3, a stand-in: no survey statistic was read). The generator's census sits under it (GEN.1343: 0.84 percent of the sub-Neptunes of
2 to 10 Earth masses inside 100 days pass `is_super_puff`)."""

START_SHARE = 0.05
"""float: The share of the sub-Neptune cores (2 to 10 Earth masses, inside 100 days) that start as a super-puff (`envelope_capture.draw_puff_start`, GEN.1343; row PF3, T8, L3): half of the
tenth that Lee and Chiang bound the planets over 4 Earth radii to. No occurrence statistic of cores that accrete a dust-free outer-disk atmosphere was found; the census cannot set it, because the XUV loss and the
age term strip or contract most puffs (started all, 11 percent end under 0.1 g/cc) and the young planets that pass the gate already come from the sub-Neptune start of GEN.1323."""

LOPEZ_MAX_MASS_EARTH = 30.0
"""float: The heaviest planet the law of Lopez and Fortney is used for here (the range of `envelope_escape`)."""

TIDAL_MIN_FLUX_EARTH = 20.0
"""float: The flux from which Millholland's fit holds, Earth fluxes (Sec. 4.1: "F/F_earth >= 20", row PF4)."""

TIDAL_MIN_LUMINOSITY_RATIO = 1.0e-5
"""float: The tidal luminosity over the irradiation over which the radius inflates (Sec. 4.1, row PF4)."""

TIDAL_MAX_RATIO = 2.0
"""float: The largest inflation: "it caps at about 2" (Sec. 4.1, row PF4)."""

TIDAL_Q_PRIME_RANGE = (1.0e3, 1.0e5)
"""tuple: The tidal quality factor `Q'` a sub-Neptune draws, log-uniform (STAND-IN PF5, T6, L1): the range of the planets of [R3079]'s tables (about 1e3 to 1e5; Neptune's own is
1e4 to 1e5). It is drawn per planet by `envelope_capture.tidal_q_prime`; the obliquity is the planet's own (the spin step)."""

MIN_OBLIQUITY_DEG = 0.01
"""float: Under this tilt, degrees, a planet has no obliquity tide (`Q` runs to infinity as the sine falls; the Cassini-state tilts of the generator's synchronous planets are 0.01 to 1 degree)."""

_TIDAL_FIT = (0.29, -0.90, 0.91, 3.11, -0.88)
"""tuple: [R3079] equation 14: `R(tides) / R(no tides) = 1 + a (M / 10 Me)^b (f_env / 0.05)^c (F / 100 F_earth)^d (Q / 1e5)^e`, with
`Q = Q' (1 + cos^2 eps) / sin^2 eps`."""

_R_EARTH_M = constants.EARTH_RADIUS_M
_INSOLATION_W_M2 = constants.SOLAR_CONSTANT_W_M2


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


def is_super_puff(mass_earth, radius_earth_):
    """Whether a planet is a super-puff: up to `CLASS_MAX_MASS_EARTH` and of a density under the top of `DENSITY_RANGE_G_CM3`, at any age and at any period (GEN.1470: the density is the
    definition; GEN.1347: the obliquity tide keeps one inflated, and Kepler-79 d is 1.9 to 3.4 Gyr old; Kepler-51 d is at 130 days)."""
    return mass_earth <= CLASS_MAX_MASS_EARTH and density_g_cm3(mass_earth, radius_earth_) <= DENSITY_RANGE_G_CM3[1]


def q_term(q_prime, obliquity_deg):
    """The `Q` of [R3079]: `Q' (1 + cos^2 eps) / sin^2 eps`, or `math.inf` for a planet whose obliquity is under `MIN_OBLIQUITY_DEG` (no obliquity tide)."""
    epsilon = math.radians(obliquity_deg)
    sine = math.sin(epsilon)
    if abs(sine) < math.sin(math.radians(MIN_OBLIQUITY_DEG)):
        return math.inf
    return q_prime * (1.0 + math.cos(epsilon) ** 2) / sine ** 2
