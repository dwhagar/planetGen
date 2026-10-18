# planetgen/galaxy/structure_draw.py

"""
The structure of a spiral galaxy drawn from its seed (GEN.532, GEN.533; docs/design/galaxy-laws-research.md Q1, Q3,
Q5 and docs/design/gap-filling-plan.md 4.3, 4.4, 4.9, 4.11), the part of the shape that is not a stored column
(`galaxy.structure`). The same seed gives the same structure; each link draws from `unit_seed(galaxy_seed,
"structure", <link>)` through the one random wrapper (`planetgen.util.draw.Stream`), so the draws neither depend on
each other nor on the order they run in.

1. **Circular speed** `Vc` from the baryonic Tully-Fisher relation `M_b = 47 Vf^4` (McGaugh 2012, arXiv:1107.2934)
   with `M_b` taken as the stellar mass of the mass chain (`galaxy.mass_chain`), scattered by 0.0375 dex (the 0.15
   dex intrinsic upper limit in mass, over the slope 4). STAND-IN until the potential of GEN.338 gives the circular
   speed from the enclosed mass; the gas is not added.
2. **Thick disk** (fills GR3 and the thick-disk rows of the gap plan; Boss's cards of 12:02 PM EDT and 1:12 PM EDT):
   - mass over the thin disk's: the centre is drawn log-uniformly between the two fits (both verified, "widen and keep
     both"): `0.53 (Vc / 100)^-2.3` (Yoachim and Dalcanton 2006, astro-ph/0508460, sections 4.4.5 and 4.4.6) and
     `10^2.113 Vc^-1.086 = 0.87 (Vc / 100)^-1.086` (Comeron, Salo and Knapen 2018, A&A 610, A5, VizieR J/A+A/610/A5,
     124 galaxies, scatter 0.268 dex), and scattered by Comeron's 0.268 dex, the wider of the two;
   - scale length over the thin disk's: Comeron 2018's fit for the outermost exponential sections,
     `10^(-1.638 + 0.817 log Vc)` with 0.316 dex (104 galaxies; newest of two verified sources, the older, Yoachim and
     Dalcanton's `1.15 (Vc / 100)^-0.2`, stays in the record);
   - scale height over the thin disk's: drawn log-uniformly from 2.5 (Yoachim and Dalcanton, section 4.4) to 4.2
     (Comeron 2018's median, 4.24; "widen and keep both") with Comeron's 0.13 dex.
   The hot component (thick disk plus classical bulge) bounds the mean mass ratio from above: `M_h / (M* - M_h)` with
   `log M_h = 0.489 + 0.9161 log M*` (Comeron et al. 2021, arXiv:2101.04478 Eq. 1); galaxy by galaxy 39 percent of
   Comeron's 124 exceed it, so the bound holds for the mean relation or with a margin of 0.13 dex (a test).
3. **Stellar halo** (fill GR1): its mass over the stellar mass is `f = 0.05 (M* / 5e10)^0.79` (placeholder,
   Harmsen et al. 2017, arXiv:1611.05448 Table 1 and D'Souza et al. 2014, arXiv:1404.2123 Fig. 16) drawn with a
   scatter of 1.01 dex (Dragonfly: Merritt et al. 2016, arXiv:1606.08847; the widest of it and Harmsen's 0.33 dex,
   Boss's rule of 1:12 PM EDT), cut so the halo is 0.1 percent to 30 percent of the galaxy's stars (the upper end is
   D'Souza 2014's disc-dominated galaxies; below it the Dragonfly galaxies with no detected halo).
4. **Motions** (fills GR4, GR5): the birth-gas dispersion 6.5 +/- 1 km/s (a molecular stand-in, GEN.536; the atomic value is 10 +/- 2, Tamburro 2009), the disc's central vertical
   dispersion `0.248 V` with 15 percent (Martinsson 2013), the age exponent 0.441 +/- 0.03 (Sharma and Hayden 2020),
   and the bulge's scatter 0.097 dex about `Vc / sigma_0 = 2 (1 - B/T)` (Courteau 2007).
5. **Bar** (GEN.533; sources and conflicts in `docs/design/gen533-bar-sources.md`): whether the galaxy has a bar, the
   bar's (or an unbarred bulge's) axis ratios and the bar's angle; see `draw_bar`.

All the scatters are log-normal (normal in the log of the quantity), cut at 3 sigma so no seed draws a number the
sources do not span.
"""

import math

from planetgen import tuning
from planetgen.galaxy import mass_chain, seed as galaxySeed
from planetgen.galaxy.structure import Structure
from planetgen.util import draw

BTFR_MASS_OVER_V4 = 47.0
"""float: McGaugh 2012: `M_b = 47 Vf^4`, solar masses and km/s."""
BTFR_V_SCATTER_DEX = 0.15 / 4.0
"""float: The intrinsic scatter's upper limit, 0.15 dex in mass, as dex in the circular speed."""

YOACHIM_MASS_RATIO_AT_100 = 0.53
YOACHIM_MASS_RATIO_SLOPE = -2.3
COMERON_MASS_RATIO_INTERCEPT = 2.113
COMERON_MASS_RATIO_SLOPE = -1.086
COMERON_MASS_RATIO_SCATTER_DEX = 0.268
COMERON_LENGTH_RATIO_INTERCEPT = -1.638
COMERON_LENGTH_RATIO_SLOPE = 0.817
COMERON_LENGTH_RATIO_SCATTER_DEX = 0.316
THICK_HEIGHT_RATIO_RANGE = (2.5, 4.2)
THICK_HEIGHT_RATIO_SCATTER_DEX = 0.13

HOT_COMPONENT_INTERCEPT = 0.489
HOT_COMPONENT_SLOPE = 0.9161
HOT_COMPONENT_SCATTER_DEX = 0.13

HALO_RATIO_AT_5E10 = 0.05
HALO_RATIO_SLOPE = 0.79
HALO_RATIO_SCATTER_DEX = 1.01
HALO_RATIO_RANGE = (0.001 / 0.999, 0.30 / 0.70)

GAS_DISPERSION_KMS = (6.5, 1.0)
"""tuple: `(centre, scatter)` km/s of the cold gas the stars are born from. A STAND-IN until GEN.536 reads the measured
CO (molecular) dispersions of nearby spirals: stars form in molecular clouds, whose dispersion is lower than the atomic
gas's 10 +/- 2 km/s (Tamburro et al. 2009, HI at r25 in 11 THINGS discs) that this value replaced. At 6.5 km/s the
young disc is 0.22 of the old disc's height (65 pc for 300 pc; observed about 50 pc), where 10 km/s gave 0.40 and put
twice the molecular clouds in the Galaxy that the catalogues hold. Not read from a paper; tunable."""
DISC_CENTRAL_DISPERSION_OVER_V = (0.248, 0.15)
"""tuple: `(centre, relative scatter)`: Martinsson et al. 2013, `0.248 +/- 0.038`, 15 percent."""
DISC_AGE_EXPONENT = (0.441, 0.03)
BULGE_DISPERSION_SCATTER_DEX = 0.097
SCATTER_LIMIT = 3.0
"""float: How many sigma a scatter may reach."""


def _stream(galaxy_seed, link):
    return draw.Stream(galaxySeed.unit_seed(galaxy_seed or mass_chain.FALLBACK_SEED, "structure", link))


def _scatter(stream, sigma):
    """A normal draw of width `sigma`, redrawn while it lies beyond `SCATTER_LIMIT` sigma (a few draws; the last is
    clipped)."""
    for _ in range(50):
        value = stream.gauss(0.0, sigma)
        if abs(value) <= SCATTER_LIMIT * sigma:
            return value
    return 0.0


def circular_speed_kms(stellar_mass_solar):
    """The median circular speed, km/s, of a galaxy of `stellar_mass_solar` (the baryonic Tully-Fisher relation)."""
    return (stellar_mass_solar / BTFR_MASS_OVER_V4) ** 0.25


def thick_mass_ratio_bounds(vc_kms):
    """`(Yoachim and Dalcanton, Comeron)`: the two fits' thick-over-thin mass ratio at `vc_kms`."""
    yoachim = YOACHIM_MASS_RATIO_AT_100 * (vc_kms / 100.0) ** YOACHIM_MASS_RATIO_SLOPE
    comeron = 10.0 ** (COMERON_MASS_RATIO_INTERCEPT + COMERON_MASS_RATIO_SLOPE * math.log10(vc_kms))
    return yoachim, comeron


def thick_length_ratio_median(vc_kms):
    """Comeron et al. 2018's median thick-over-thin scale length ratio at `vc_kms`."""
    return 10.0 ** (COMERON_LENGTH_RATIO_INTERCEPT + COMERON_LENGTH_RATIO_SLOPE * math.log10(vc_kms))


def hot_component_bound(stellar_mass_solar):
    """The most the thick disk can weigh over the thin disk's, if the hot component (thick disk and classical bulge)
    is all thick disk: `M_h / (M* - M_h)` (Comeron et al. 2021 Eq. 1)."""
    hot = 10.0 ** (HOT_COMPONENT_INTERCEPT + HOT_COMPONENT_SLOPE * math.log10(stellar_mass_solar)) / stellar_mass_solar
    return hot / (1.0 - hot)


def halo_ratio_median(stellar_mass_solar):
    """The placeholder median of the stellar halo's mass over the galaxy's stellar mass."""
    return HALO_RATIO_AT_5E10 * (stellar_mass_solar / 5.0e10) ** HALO_RATIO_SLOPE


def draw_structure(galaxy_seed):
    """
    The structure the galaxy of `galaxy_seed` has: the same seed gives the same structure.

    Args:
        galaxy_seed (bytes or None): The galaxy's 16-byte seed; `None` uses the mass chain's fallback seed.

    Returns:
        structure.Structure
    """
    mass = mass_chain.draw_mass_chain(galaxy_seed).stellar_mass_solar
    vc = circular_speed_kms(mass) * 10.0 ** _scatter(_stream(galaxy_seed, "speed"), BTFR_V_SCATTER_DEX)

    low, high = thick_mass_ratio_bounds(vc)
    centre = _stream(galaxy_seed, "thick-mass-centre").uniform(math.log10(min(low, high)), math.log10(max(low, high)))
    mass_ratio = 10.0 ** (centre + _scatter(_stream(galaxy_seed, "thick-mass"), COMERON_MASS_RATIO_SCATTER_DEX))
    length_ratio = thick_length_ratio_median(vc) * 10.0 ** _scatter(
        _stream(galaxy_seed, "thick-length"), COMERON_LENGTH_RATIO_SCATTER_DEX)
    height_centre = _stream(galaxy_seed, "thick-height-centre").uniform(
        math.log10(THICK_HEIGHT_RATIO_RANGE[0]), math.log10(THICK_HEIGHT_RATIO_RANGE[1]))
    height_ratio = 10.0 ** (height_centre + _scatter(_stream(galaxy_seed, "thick-height"), THICK_HEIGHT_RATIO_SCATTER_DEX))

    halo = halo_ratio_median(mass) * 10.0 ** _scatter(_stream(galaxy_seed, "halo"), HALO_RATIO_SCATTER_DEX)
    halo = min(max(halo, HALO_RATIO_RANGE[0]), HALO_RATIO_RANGE[1])

    gas = GAS_DISPERSION_KMS[0] + _scatter(_stream(galaxy_seed, "gas-dispersion"), GAS_DISPERSION_KMS[1])
    centre_over_v, relative = DISC_CENTRAL_DISPERSION_OVER_V
    central = vc * centre_over_v * math.exp(_scatter(_stream(galaxy_seed, "disc-dispersion"), relative))
    exponent = DISC_AGE_EXPONENT[0] + _scatter(_stream(galaxy_seed, "disc-heating"), DISC_AGE_EXPONENT[1])
    bulge_scatter = _scatter(_stream(galaxy_seed, "bulge-dispersion"), BULGE_DISPERSION_SCATTER_DEX)
    bar = draw_bar(galaxy_seed, mass)
    return Structure(
        vc_kms=vc, thick_mass_ratio=mass_ratio, thick_length_ratio=length_ratio, thick_height_ratio=height_ratio,
        bar_present=bar[0], bar_axis_ratio_y=bar[1], bar_axis_ratio_z=bar[2], bar_angle_deg=bar[3],
        halo_stellar_ratio=halo,
        disc_gas_dispersion_kms=gas, disc_central_dispersion_kms=max(central, gas * 1.05), disc_age_exponent=exponent,
        bulge_dispersion_scatter_dex=bulge_scatter)


BAR_FRACTION_COEFFICIENTS = (-82.2, 17.1, -0.88)
"""tuple: Erwin 2018 (arXiv:1711.04867) Eq. 6 and Section 3: the fraction of S4G spirals with a bar (strong or weak) is
`1 / (1 + exp(-(a + b1 m + b2 m^2)))`, `m = log10(M*/Msun)`, fitted over m = 8 to 12."""
BAR_ELLIPTICITY = (0.53, 0.15)
BAR_ELLIPTICITY_RANGE = (0.15, 0.90)
"""tuple: Diaz-Garcia et al. 2016 (arXiv:1509.06743) Table 3 and the VizieR table J/A+A/587/A160 (504 bars): the
maximum isophotal ellipticity of a bar, mean 0.534 and standard deviation 0.148, independent of the host's mass and
Hubble type; its range 0.10 to 0.92 cut at Aguerri et al. 2009's smallest detected bar, 0.16. The model's `y0/x0`
is `1 - e` (the isodensity contour in the plane has that axis ratio): a STAND-IN for the bar body, which Gadotti 2011
(median 0.63, SD 0.10) finds more elongated than its isophotes."""
BOXY_PEANUT_COEFFICIENTS = (-42.79, 4.13)
"""tuple: Erwin and Debattista 2017 (arXiv:1703.01602) Sect. 4.1 Table 2: the fraction of barred galaxies whose bar has a
boxy or peanut bulge, `1 / (1 + exp(-(a + b m)))` (half at m = 10.4)."""
BOXY_PEANUT_HEIGHT = (0.23, 0.06)
BOXY_PEANUT_HEIGHT_RANGE = (0.12, 0.38)
"""tuple: The boxy/peanut bar's `z0/x0`: STAND-IN, the peanut's half-length over the bar's (0.42 +/- 0.09, Erwin and
Debattista 2017) times its height over its length (0.5 to 0.6, Ciambur and Graham 2016, arXiv:1603.00019): 0.23, 25
percent scatter; no external-galaxy measurement of the ratio was found."""
THIN_BAR_HEIGHT = (0.071, 0.02)
THIN_BAR_HEIGHT_RANGE = (0.05, 0.15)
"""tuple: A bar without a peanut: thickness over length 1 / (14 +/- 4) in edge-on near-infrared images (Luetticke et al.
2000, astro-ph/0008478); the upper end is the simulated 0.1 to 0.2 (Martinez-Valpuesta et al. 2006)."""
UNBARRED_AXIS_RATIO_Y = (0.83, 0.21)
UNBARRED_AXIS_RATIO_Y_RANGE = (0.30, 1.0)
UNBARRED_AXIS_RATIO_Z = (0.54, 0.19)
UNBARRED_AXIS_RATIO_Z_RANGE = (0.19, 1.0)
"""tuple: The bulge of an unbarred spiral: B/A and C/A of the 41 unbarred CALIFA bulges of Costantin et al. 2018
(arXiv:1710.05222) Table C.1 (means and standard deviations computed from the table; the truncation to the table's range
and C/A not above B/A are a STAND-IN). Mendez-Abreu et al. 2010 give a thicker bulge (C/A near 0.8); the newest source
is taken and the older kept in the record."""


def _clamped_normal(stream, centre, sigma, limits):
    return min(max(stream.gauss(centre, sigma), limits[0]), limits[1])


def bar_fraction(stellar_mass_solar):
    """The share of spirals of this stellar mass that have a bar (Erwin 2018 Eq. 6), the same for early and late
    stages (the S4G split, 0.63 against 0.66, is 1.1 sigma: Diaz-Garcia et al. 2016 Table 1)."""
    m = math.log10(stellar_mass_solar)
    a, b1, b2 = BAR_FRACTION_COEFFICIENTS
    return 1.0 / (1.0 + math.exp(-(a + b1 * m + b2 * m * m)))


def boxy_peanut_fraction(stellar_mass_solar):
    """The share of barred galaxies of this stellar mass whose bar has a boxy or peanut bulge (Erwin and Debattista
    2017)."""
    a, b = BOXY_PEANUT_COEFFICIENTS
    return 1.0 / (1.0 + math.exp(-(a + b * math.log10(stellar_mass_solar))))


def draw_bar(galaxy_seed, stellar_mass_solar):
    """
    Whether the galaxy has a bar and the axes of its bulge: the same seed gives the same bar.

    A barred galaxy's bulge is the bar: `y0/x0` from the bar's ellipticity and `z0/x0` from a boxy/peanut or a thin
    bar. An unbarred galaxy's is a spheroid with the flattenings of the CALIFA bulges. The bar's angle from the
    Sun-centre line is uniform on 0 to 180 degrees (a bar has no direction; by symmetry nothing prefers one: a
    STAND-IN, no measured test was found).

    Returns:
        tuple: `(present, y0/x0, z0/x0, angle in degrees)`.
    """
    present = _stream(galaxy_seed, "bar").random() < bar_fraction(stellar_mass_solar)
    angle = _stream(galaxy_seed, "bar-angle").uniform(0.0, 180.0)
    if present:
        ellipticity = _clamped_normal(_stream(galaxy_seed, "bar-width"), *BAR_ELLIPTICITY, BAR_ELLIPTICITY_RANGE)
        axis_y = 1.0 - ellipticity
        peanut = _stream(galaxy_seed, "bar-peanut").random() < boxy_peanut_fraction(stellar_mass_solar)
        centre, sigma = BOXY_PEANUT_HEIGHT if peanut else THIN_BAR_HEIGHT
        limits = BOXY_PEANUT_HEIGHT_RANGE if peanut else THIN_BAR_HEIGHT_RANGE
        axis_z = _clamped_normal(_stream(galaxy_seed, "bar-height"), centre, sigma, limits)
    else:
        axis_y = _clamped_normal(_stream(galaxy_seed, "bulge-width"), *UNBARRED_AXIS_RATIO_Y, UNBARRED_AXIS_RATIO_Y_RANGE)
        axis_z = _clamped_normal(_stream(galaxy_seed, "bulge-height"), *UNBARRED_AXIS_RATIO_Z, UNBARRED_AXIS_RATIO_Z_RANGE)
    return present, axis_y, min(axis_z, axis_y), angle
