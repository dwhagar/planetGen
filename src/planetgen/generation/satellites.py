# planetgen/generation/satellites.py

"""
The satellite galaxies of a galaxy, drawn from its seed and its halo (GEN.369, GEN.624; docs/design/galaxy/intergalactic-and-halo-contents.md
5.5 and 5.13, row 10; docs/design/galaxy/galaxy-laws-research.md Q3, Q4, Q6 and Q7). Nothing here is a real galaxy's satellite: the
Milky Way's own rows are reference data for the plausibility tests (`tests/reference_satellites.py`).

The chain, each link from the one before it and each scatter drawn through `unit_seed(galaxy_seed, "satellites", <link>)`:

1. **Subhaloes** of the host's halo mass `Mh` (200 times critical, `galaxy.mass_chain`): the number above a mass fraction
   `mu = Msub / Mh` is `N(> mu) = (mu / MU1)^A exp[-(mu / MU_CUT)^B]` (Boylan-Kolchin et al. 2010, arXiv:0911.4484, Eqs. 7
   and 8; about 640 above `mu = 1e-5`), the count drawn as a Poisson number and each mass from that law. The mass is taken
   as the peak mass, the one a galaxy forms by.
2. **Occupation**: a subhalo hosts a galaxy with the probability `0.5 [1 + erf((log Mpeak - log M50) / (sqrt(2) sigma))]`
   (Nadler et al. 2020, arXiv:1912.03303, Eq. 3: `M50` at its 95 percent upper limit, 8.5e7 solar masses, `sigma` 0.03 dex).
3. **Stellar mass** from the peak mass by the stellar-to-halo relation of the mass chain (Moster, Naab and White 2013,
   `mass_chain.median_stellar_mass`, extrapolated down to the dwarf end, which the research note marks as an underestimate
   of the faintest galaxies), scattered by 0.19 dex (Nadler's 95 percent upper limit for satellites). A galaxy under
   `MIN_STELLAR_MASS_SOLAR` is not made.
4. **Size** from the stellar mass: Carlsten et al. 2021 (ELVES I, arXiv:2105.03435, Eq. 5) `log r_e = 1.071 + 0.247 log M*`
   with 0.181 dex below `log M* = 8.5`, and above it the late-type size-mass law of the galaxy shape (Shen et al. 2003,
   `galaxy.shape_draw.size_median_kpc`), which meets it there to 5 percent.
5. **Metallicity**: `<[Fe/H]> = -1.69 + 0.30 log(M* / 1e6)` (Kirby et al. 2013, ApJ 779, 102) with 0.17 dex about it (0.09
   for the gas-rich ones).
6. **Profile and age**: exponential for a gas-rich galaxy (`M*` over `EXPONENTIAL_MIN_STELLAR_MASS_SOLAR`), Plummer for a
   dwarf spheroidal; the mean age falls with mass along `AGE_AT_1E5_GY` and `AGE_SLOPE_GY_PER_DEX`, STAND-INS anchored on
   the quenching times of Weisz et al. 2015 (arXiv:1503.05195: the faintest quenched 10 Gyr ago, a few-1e8 galaxy still
   forming).
7. **Gas by distance (GEN.1193)**: dwarfs of the Local Group within about 270 kpc of the Milky Way or Andromeda mostly hold no
   HI (under 1e5 solar masses) and those beyond mostly do, 1e5 to 1e8 (Grcevich and Putman 2009, arXiv:0901.4975, abstract and
   section 3.4 read: the cutoff lies between 250 and 280 kpc, 14 of the 16 dwarfs beyond 280 kpc are detected, and the ram
   pressure of the host's halo gas is the likely cause). The cutoff scales with the host's virial radius
   (`STRIPPING_RADIUS_MILKY_WAY_KPC`). Inside it a satellite lighter than `STRIPPING_RETAINED_MIN_STELLAR_MASS_SOLAR` (the
   Magellanic Clouds, left out of that survey for their mass, are the lightest that kept theirs) has lost its gas: its
   mean age is raised to the oldest at which it forms no stars. Beyond it, a share `GAS_DETECTED_BEYOND_STRIPPING` of the
   satellites keeps its gas: its mean age is lowered to the oldest at which it still forms stars, if it forms at any.
   The age is drawn first from the mass law (6), so this moves only those near the star formation threshold.
8. **Place**: the distance from the host's centre from an Einasto profile far less concentrated than the dark matter
   (`EINASTO_ALPHA`, `EINASTO_R_MINUS2_OVER_R200`: Nadler 2020, a calibration on one halo), between `INNER_OVER_R200` of
   the virial radius (STAND-IN: the host's disc destroys what comes closer) and the virial radius, in an isotropic
   direction. Nothing lies beyond the virial radius, and the radius is the halo's, not the disc's.

The velocities are zero: the sector gate reads only the position, and the satellites' orbits belong with the orbital
update (GEN.339). A galaxy with no halo (`mass_chain.has_halo`) has no satellites.
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.galaxy import components, mass_chain, satellite_history, seed as galaxySeed, shape_draw
from planetgen.physics import component_profiles
from planetgen.util import draw
from planetgen.util.poisson import poisson_count

SUBHALO_MU1 = 0.01
SUBHALO_MU_CUT = 0.096
SUBHALO_SLOPE = -0.935
SUBHALO_CUTOFF_EXPONENT = 1.29
"""float: The constants of the subhalo count `N(> mu) = (mu / MU1)^a exp[-(mu / MU_CUT)^b]` (Boylan-Kolchin 2010 Eqs. 7, 8)."""

SUBHALO_MU_MAX = 0.2
"""float: STAND-IN: the heaviest subhalo drawn, as a fraction of the host's halo: a satellite is the lighter of a pair, and a
subhalo near the host's own mass is a merger, not a satellite."""

SATELLITE_MAX_STELLAR_FRACTION = 0.25
"""float: STAND-IN: the most stellar mass a satellite holds, as a fraction of its host's (the Large Magellanic Cloud is about a
tenth of the Milky Way's): a host drawn light for its halo would otherwise out-weigh itself with a satellite."""

OCCUPATION_M50_SOLAR = 8.5e7
OCCUPATION_SIGMA_DEX = 0.03
"""float: The centre and width of the peak-mass occupation curve (Nadler 2020)."""
OCCUPATION_FLOOR_SOLAR = 1.0e6
"""float: The lightest peak mass a detectable galaxy needs (Nadler 2020, section 7.5)."""

STELLAR_SCATTER_DEX = 0.19
MIN_STELLAR_MASS_SOLAR = 1.0e3
"""float: The lightest satellite made, solar masses of stars: about the faintest known satellites (the ultra-faint dwarfs
hold a few hundred to a few thousand). STAND-IN: below it a population is a handful of stars in a cell."""

CARLSTEN_INTERCEPT = 1.071
CARLSTEN_SLOPE = 0.247
CARLSTEN_SCATTER_DEX = 0.181
SIZE_JOIN_LOG_MASS = 8.5
"""float: The `log10 M*` where the dwarf size law hands over to the disc size law."""

KIRBY_FE_H_AT_1E6 = -1.69
KIRBY_FE_H_SLOPE = 0.30
KIRBY_SCATTER_DEX = 0.17
KIRBY_IRREGULAR_SCATTER_DEX = 0.09

EXPONENTIAL_MIN_STELLAR_MASS_SOLAR = 1.0e8
"""float: STAND-IN: the stellar mass from which a satellite is a gas-rich, still star-forming exponential disc (Weisz 2015:
about 40 percent of the galaxies near 1e8 are quenched today, all of the fainter ones)."""
AGE_AT_1E5_GY = 11.5
AGE_SLOPE_GY_PER_DEX = -1.2
AGE_SCATTER_GY = 1.5
AGE_RANGE_GY = (2.0, 13.0)
"""float: STAND-INS for the mean stellar age of a satellite against its stellar mass (see the module docstring)."""

MILKY_WAY_STELLAR_MASS_SOLAR = 5.0e10
"""float: The Milky Way's stellar mass (as `galaxy.structure` takes it), whose median halo sets the stripping radius's scale."""
STRIPPING_RADIUS_MILKY_WAY_KPC = 270.0
"""float: The distance from a host inside which its dwarfs have lost their gas, for the Milky Way's halo: Grcevich and Putman
2009 (arXiv:0901.4975, abstract and section 3.4 read) put the cutoff between 250 and 280 kpc. Scaled with the host's virial
radius [I]: the ram pressure that strips them is the halo gas's."""
STRIPPING_RETAINED_MIN_STELLAR_MASS_SOLAR = 4.0e8
"""float: STAND-IN: the lightest satellite that keeps its gas inside the stripping radius. The Small Magellanic Cloud (4.6e8
solar masses of stars, 64 kpc) holds its HI, and the survey leaves both Clouds out for their large mass; nothing lighter is
known to. Between 2e8 and 4.6e8."""
STRIPPING_AGE_MARGIN_GY = 0.01
"""float: How far past the oldest age that forms stars a stripped satellite's age is put, Gy."""
GAS_DETECTED_BEYOND_STRIPPING = 0.875
"""float: The share of satellites beyond the stripping radius that hold gas: 14 of the 16 dwarfs beyond 280 kpc are detected
(Tucana and Cetus are not; Grcevich and Putman 2009, section 3.4 read)."""

EINASTO_ALPHA = 0.678
EINASTO_R_MINUS2_OVER_R200 = 0.81
INNER_OVER_R200 = 0.05
RADIAL_STEPS = 4000


def subhalo_count_above(mu):
    """The mean number of subhaloes above the mass fraction `mu` of the host's halo."""
    return (mu / SUBHALO_MU1) ** SUBHALO_SLOPE * math.exp(-((mu / SUBHALO_MU_CUT) ** SUBHALO_CUTOFF_EXPONENT))


def occupation_probability(peak_mass_solar):
    """The probability that a subhalo of this peak mass hosts a galaxy (Nadler 2020 Eq. 3)."""
    if peak_mass_solar < OCCUPATION_FLOOR_SOLAR:
        return 0.0
    along = (math.log10(peak_mass_solar) - math.log10(OCCUPATION_M50_SOLAR)) / (math.sqrt(2.0) * OCCUPATION_SIGMA_DEX)
    return 0.5 * (1.0 + math.erf(along))


def lightest_peak_mass_solar():
    """The lightest peak mass worth drawing: the one whose median stellar mass, three scatters (`STELLAR_SCATTER_DEX`)
    up, reaches `MIN_STELLAR_MASS_SOLAR` (a lighter halo cannot make a galaxy of the floor's mass), and at least the
    occupation floor."""
    target = math.log10(MIN_STELLAR_MASS_SOLAR) - 3.0 * STELLAR_SCATTER_DEX
    low, high = 6.0, 14.0
    for _ in range(60):
        middle = 0.5 * (low + high)
        if math.log10(mass_chain.median_stellar_mass(10.0 ** middle)) < target:
            low = middle
        else:
            high = middle
    return max(10.0 ** (0.5 * (low + high)), OCCUPATION_FLOOR_SOLAR)


def _stream(galaxy_seed, link):
    return draw.Stream(galaxySeed.unit_seed(galaxy_seed or mass_chain.FALLBACK_SEED, "satellites", link))


def _scatter(stream, sigma, limit=3.0):
    for _ in range(50):
        value = stream.gauss(0.0, sigma)
        if abs(value) <= limit * sigma:
            return value
    return 0.0


def _subhalo_fraction(stream, mu_min):
    """A mass fraction from `mu_min` to `SUBHALO_MU_MAX` from the count law: the inverse of `N(> mu)` over that range by
    bisection in log mu."""
    top = subhalo_count_above(SUBHALO_MU_MAX)
    target = top + stream.random() * (subhalo_count_above(mu_min) - top)
    low, high = math.log10(mu_min), math.log10(SUBHALO_MU_MAX)
    for _ in range(60):
        middle = 0.5 * (low + high)
        if subhalo_count_above(10.0 ** middle) > target:
            low = middle
        else:
            high = middle
    return 10.0 ** (0.5 * (low + high))


def _radius_table(r200_pc):
    """`(radii, cumulative share)` of the Einasto distribution of subhaloes between the inner limit and the virial radius."""
    r_minus2 = EINASTO_R_MINUS2_OVER_R200 * r200_pc
    radii = np.linspace(INNER_OVER_R200 * r200_pc, r200_pc, RADIAL_STEPS)
    rho = np.exp(-(2.0 / EINASTO_ALPHA) * ((radii / r_minus2) ** EINASTO_ALPHA - 1.0))
    shell = radii ** 2 * rho
    cumulative = np.concatenate(([0.0], np.cumsum(0.5 * (shell[1:] + shell[:-1]) * np.diff(radii))))
    return radii, cumulative / cumulative[-1]


def _direction(stream):
    z = stream.uniform(-1.0, 1.0)
    azimuth = stream.uniform(0.0, 2.0 * math.pi)
    xy = math.sqrt(max(1.0 - z * z, 0.0))
    return xy * math.cos(azimuth), xy * math.sin(azimuth), z


def half_light_radius_pc(stellar_mass_solar, deviation):
    """A satellite's half-light radius, parsecs, at `deviation` standard deviations about the size law of its mass (the
    dwarf law below `SIZE_JOIN_LOG_MASS`, the disc law above)."""
    log_mass = math.log10(stellar_mass_solar)
    if log_mass < SIZE_JOIN_LOG_MASS:
        return 10.0 ** (CARLSTEN_INTERCEPT + CARLSTEN_SLOPE * log_mass + deviation * CARLSTEN_SCATTER_DEX)
    sigma_dex = shape_draw.size_sigma_ln(stellar_mass_solar) / math.log(10.0)
    return 1000.0 * shape_draw.size_median_kpc(stellar_mass_solar) * 10.0 ** (deviation * sigma_dex)


def mean_fe_h(stellar_mass_solar):
    """The mean [Fe/H] of a satellite of this stellar mass (Kirby et al. 2013)."""
    return KIRBY_FE_H_AT_1E6 + KIRBY_FE_H_SLOPE * math.log10(stellar_mass_solar / 1.0e6)


def mean_age_gy(stellar_mass_solar):
    """The median mean stellar age, Gyr, of a satellite of this stellar mass (`AGE_AT_1E5_GY`)."""
    return AGE_AT_1E5_GY + AGE_SLOPE_GY_PER_DEX * (math.log10(stellar_mass_solar) - 5.0)


def stripping_radius_kpc(virial_radius_kpc):
    """The distance from the host inside which its dwarfs have lost their gas: the Milky Way's, scaled with the virial radius."""
    return (STRIPPING_RADIUS_MILKY_WAY_KPC * virial_radius_kpc
            / mass_chain.virial_radius_kpc(mass_chain.median_halo_mass(MILKY_WAY_STELLAR_MASS_SOLAR)))


def age_with_gas_by_distance(mass, half_mass_radius_pc, age_gy, distance_pc, stripping_radius_pc, keeps_gas_draw):
    """
    A satellite's mean age with its gas decided by its distance from the host (see the module docstring, link 7): inside
    `stripping_radius_pc` a satellite lighter than `STRIPPING_RETAINED_MIN_STELLAR_MASS_SOLAR` that would form stars has
    its age raised to just past the oldest at which it forms them; beyond it a satellite whose `keeps_gas_draw` (uniform
    in [0, 1)) is under `GAS_DETECTED_BEYOND_STRIPPING` that forms none has its age lowered to the oldest at which it forms
    them, if there is one.
    """
    youngest, oldest = AGE_RANGE_GY
    limit = satellite_history.oldest_forming_age_gy(mass, half_mass_radius_pc, youngest, oldest)
    if distance_pc < stripping_radius_pc:
        if mass < STRIPPING_RETAINED_MIN_STELLAR_MASS_SOLAR and limit is not None and age_gy <= limit:
            return min(limit + STRIPPING_AGE_MARGIN_GY, oldest)
        return age_gy
    if keeps_gas_draw < GAS_DETECTED_BEYOND_STRIPPING and limit is not None and age_gy > limit:
        return limit
    return age_gy


def _seed(galaxy_seed, name):
    text = galaxy_seed.hex() if isinstance(galaxy_seed, (bytes, bytearray)) else str(galaxy_seed)
    return draw.Stream(f"satellite:{text}:{name}").getrandbits(63)


def draw_satellites(galaxy_seed, chain, epoch_unix=0.0):
    """
    The galaxy's satellite galaxies: the same seed and halo give the same population.

    Args:
        galaxy_seed (bytes or int): The galaxy's seed.
        chain (galaxy.mass_chain.MassChain): The galaxy's stellar mass and halo (`halo_mass_solar` and `virial_radius_kpc`
            are `None` for a galaxy with no halo, which has no satellites).
        epoch_unix (float): When the positions hold.

    Returns:
        list[Component]: One `satellite` component per galaxy, the heaviest first, named "Generated satellite 001" and on.
    """
    if chain.halo_mass_solar is None:
        return []
    halo_mass = chain.halo_mass_solar
    r200_pc = chain.virial_radius_kpc * 1000.0
    mu_min = lightest_peak_mass_solar() / halo_mass
    if mu_min >= SUBHALO_MU_MAX:
        return []
    count = poisson_count(subhalo_count_above(mu_min) - subhalo_count_above(SUBHALO_MU_MAX), _stream(galaxy_seed, "count"))
    mass_stream, occupation = _stream(galaxy_seed, "mass"), _stream(galaxy_seed, "occupation")
    stellar, size = _stream(galaxy_seed, "stellar-mass"), _stream(galaxy_seed, "size")
    metals, ages = _stream(galaxy_seed, "metallicity"), _stream(galaxy_seed, "age")
    retention = _stream(galaxy_seed, "gas-retention")
    places, directions = _stream(galaxy_seed, "radius"), _stream(galaxy_seed, "direction")
    radii, shares = _radius_table(r200_pc)
    drawn = []
    for _ in range(count):
        peak = _subhalo_fraction(mass_stream, mu_min) * halo_mass
        hosts_galaxy = occupation.random() < occupation_probability(peak)
        mass = min(mass_chain.median_stellar_mass(peak) * 10.0 ** _scatter(stellar, STELLAR_SCATTER_DEX),
                   SATELLITE_MAX_STELLAR_FRACTION * chain.stellar_mass_solar)
        distance = float(np.interp(places.random(), shares, radii))
        direction = _direction(directions)
        deviation = _scatter(size, 1.0)
        fe_deviation, age_deviation = _scatter(metals, 1.0), _scatter(ages, 1.0)
        keeps_gas_draw = retention.random()
        if hosts_galaxy and mass >= MIN_STELLAR_MASS_SOLAR:
            drawn.append((mass, distance, direction, deviation, fe_deviation, age_deviation, keeps_gas_draw))
    stripping_radius_pc = stripping_radius_kpc(chain.virial_radius_kpc) * 1000.0
    found = []
    for index, (mass, distance, direction, deviation, fe_deviation, age_deviation, keeps_gas_draw) in enumerate(
            sorted(drawn, key=lambda row: -row[0])):
        name = f"Generated satellite {index + 1:03d}"
        gas_rich = mass >= EXPONENTIAL_MIN_STELLAR_MASS_SOLAR
        profile = "exponential" if gas_rich else "plummer"
        radius = half_light_radius_pc(mass, deviation)
        low, high = tuning.METALLICITY_RANGE_FE_H
        scatter = KIRBY_IRREGULAR_SCATTER_DEX if gas_rich else KIRBY_SCATTER_DEX
        fe_h = min(max(mean_fe_h(mass) + fe_deviation * scatter, low), high)
        age = min(max(mean_age_gy(mass) + age_deviation * AGE_SCATTER_GY, AGE_RANGE_GY[0]), AGE_RANGE_GY[1])
        age = age_with_gas_by_distance(mass, radius, age, distance, stripping_radius_pc, keeps_gas_draw)
        found.append(components.Component(
            kind="satellite", name=name, source="generated",
            center_pc=tuple(distance * c for c in direction), velocity_kms=(0.0, 0.0, 0.0), epoch_unix=epoch_unix,
            mass_solar=mass, half_mass_radius_pc=radius, profile=profile, age_gy=age, metallicity_fe_h=fe_h,
            seed=_seed(galaxy_seed, name), tidal_radius_pc=component_profiles.cut_radius_pc(profile, radius)))
    return found
