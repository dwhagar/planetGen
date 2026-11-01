# planetgen/galaxy/structure_draw.py

"""
The structure of a spiral galaxy drawn from its seed (GEN.532, GEN.533; docs/design/galaxy/galaxy-laws-research.md Q1, Q3,
Q5 and docs/design/research/gap-filling-plan.md 4.3, 4.4, 4.9, 4.11), the part of the shape that is not a stored column
(`galaxy.structure`). The same seed gives the same structure; each link draws from `unit_seed(galaxy_seed,
"structure", <link>)` through the one random wrapper (`planetgen.util.draw.Stream`), so the draws neither depend on
each other nor on the order they run in.

1. **Circular speed** `Vc` from the baryonic Tully-Fisher relation `M_b = 47 Vf^4` [R2992]
   with `M_b` taken as the stellar mass of the mass chain (`galaxy.mass_chain`), scattered by 0.0375 dex (the 0.15
   dex intrinsic upper limit in mass, over the slope 4). STAND-IN until the potential of GEN.338 gives the circular
   speed from the enclosed mass; the gas is not added.
2. **Thick disk** (fills GR3 and the thick-disk rows of the gap plan; Boss's cards of 12:02 PM EDT and 1:12 PM EDT):
   - mass over the thin disk's: the centre is drawn log-uniformly between the two fits (both verified, "widen and keep
     both"): `0.53 (Vc / 100)^-2.3` ([R4824], astro-ph/0508460, sections 4.4.5 and 4.4.6) and
     `10^2.113 Vc^-1.086 = 0.87 (Vc / 100)^-1.086` ([R0865], VizieR J/A+A/610/A5,
     124 galaxies, scatter 0.268 dex), and scattered by Comeron's 0.268 dex, the wider of the two;
   - scale length over the thin disk's: [R0865]'s fit for the outermost exponential sections,
     `10^(-1.638 + 0.817 log Vc)` with 0.316 dex (104 galaxies; newest of two verified sources, the older, Yoachim and
     Dalcanton's `1.15 (Vc / 100)^-0.2`, stays in the record);
   - scale height over the thin disk's: drawn log-uniformly from 2.5 (Yoachim and Dalcanton, section 4.4) to 4.2
     ([R0865]'s median, 4.24; "widen and keep both") with Comeron's 0.13 dex.
   The hot component (thick disk plus classical bulge) bounds the mean mass ratio from above: `M_h / (M* - M_h)` with
   `log M_h = 0.489 + 0.9161 log M*` ([R0860: Eq. 1]); galaxy by galaxy 39 percent of
   Comeron's 124 exceed it, so the bound holds for the mean relation or with a margin of 0.13 dex (a test).
3. **Stellar halo** (fill GR1): its mass over the stellar mass is `f = 0.05 (M* / 5e10)^0.79` (placeholder,
   [R1816: Table 1] and [R1169: Fig. 16]) drawn with a
   scatter of 1.01 dex (Dragonfly: [R3041]; the widest of it and Harmsen's 0.33 dex,
   Boss's rule of 1:12 PM EDT), cut so the halo is 0.1 percent to 30 percent of the galaxy's stars (the upper end is
   [R1169]'s disc-dominated galaxies; below it the Dragonfly galaxies with no detected halo). A galaxy under
   `dials.FULL_STRUCTURE_MIN_STELLAR_MASS_SOLAR` (3e8, Boss 7:58 PM EDT 10-21) runs none of this: it has no stellar halo, no thick
   disc, no bar and no nuclear structures.
4. **Motions** (fills GR4, GR5): the birth-gas dispersion 6.5 +/- 1 km/s (the molecular CO widths of [R3114], Sun 2020, [R2830], GEN.536; the atomic value is 10 +/- 2, [R4298]), the disc's central vertical
   dispersion `0.248 V` with 15 percent [R2917], the age exponent 0.441 +/- 0.03 [R4043],
   and the bulge's scatter 0.097 dex about `Vc / sigma_0 = 2 (1 - B/T)` [R0908].
5. **Black holes** (GEN.563; `draw_black_hole_shares`): the shares of the massive stars that fail to explode at solar and at
   zero metallicity, which `galaxy.remnant_density` turns into a rate per star at a place.
6. **Bar** (GEN.533; sources and conflicts in `docs/design/galaxy/gen533-bar-sources.md`): whether the galaxy has a bar, the
   bar's (or an unbarred bulge's) axis ratios and the bar's angle; see `draw_bar`.
7. **Nuclear disc** (GEN.628): the barred galaxy's disc mass from its stellar mass, the centre log-uniform between the
   two fits of [R1487] ([R1487], Eqs. 1 and 2), scattered by a placeholder 0.3 dex;
   none for an unbarred galaxy or one under `dials.FULL_STRUCTURE_MIN_STELLAR_MASS_SOLAR`. It no longer follows the
   nuclear cluster's mass.

All the scatters are log-normal (normal in the log of the quantity), cut at 3 sigma so no seed draws a number the
sources do not span.
"""

import math

from planetgen.galaxy import dials, mass_chain, natal_kick, seed as galaxySeed
from planetgen.galaxy.structure import REFERENCE_STRUCTURE, Structure
from planetgen.physics import black_hole_formation
from planetgen.util import draw

BTFR_MASS_OVER_V4 = 47.0
"""float: [R2992]: `M_b = 47 Vf^4`, solar masses and km/s."""
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
"""tuple: `(centre, scatter)` km/s of the cold gas the stars are born from. Stars form in molecular clouds, whose line-of-sight
CO width in nearby spirals is 7.3 +/- 1.7 km/s per pixel ([R3114], 9 THINGS and HERACLES galaxies), 6.3 to 7.5
km/s mass-weighted and 4.4 to 4.9 km/s area-weighted in galaxy discs at 90 and 150 pc (Sun et al. 2020, PHANGS-ALMA, 70
galaxies), 7.1 to 7.8 km/s near the Sun [R2830], and at most 12.0 +/- 3.9 km/s from spectra stacked over inclined
discs [R0625], against the atomic gas's 10 +/- 2 km/s ([R4298], HI at r25 in 11 THINGS
discs). The centre 6.5 km/s sits inside 7.0 +/- 1.7 and gives a young disc 0.22 of the old disc's height (66 pc for 300 pc),
inside the 35 to 70 pc that Milky Way OB stars, HII regions, Cepheids and young open clusters show (GEN.1182; 7.0 gives
73 pc, 10 km/s 119 pc and twice the molecular clouds the catalogues hold). Confident unverified (the CO widths are read, the
step to the stars' birth dispersion is physics), tweak L1 (the weighting); reports/gen536-co-velocity-dispersion-2026-10-31.md."""
DISC_CENTRAL_DISPERSION_OVER_V = (0.248, 0.15)
"""tuple: `(centre, relative scatter)`: [R2917], `0.248 +/- 0.038`, 15 percent."""
DISC_AGE_EXPONENT = (0.441, 0.03)
BULGE_DISPERSION_SCATTER_DEX = 0.097
NUCLEAR_CLUSTER_MASS_AT_1E9 = 6.51
"""float: `log10` of a nuclear star cluster's mass, solar masses, in a galaxy of 1e9 solar masses of stars ([R3273: Eq. 1]: `log M_NSC = 0.48 log(M* / 1e9) + 6.51`)."""
NUCLEAR_CLUSTER_SLOPE_RANGE = (0.48, 1.36)
"""tuple: The slope of that relation: uniform between Neumayer's 0.48 and [R1538]'s 1.36 ([R1538]:
1.0 late types, 1.36 early types), the two verified papers that disagree; both stay (Boss, gap plan decision 2, GX3)."""
NUCLEAR_CLUSTER_SCATTER_DEX = 0.6
"""float: Its scatter, dex ([R3273], mostly intrinsic)."""
NUCLEAR_DISC_MASS_AT_1E9 = ((7.62, 1.32), (7.99, 0.86))
"""tuple: `(log10 of the nuclear disc's mass, solar masses, in a galaxy of 1e9 solar masses of stars, slope)` of the two
fits of [R1487] ([R1487], to appear in MNRAS Letters, Eqs. 1 and 2, 16 TIMER barred
galaxies, bootstrapped errors 0.44 and 0.25 on the intercept and slope): `log M_ND = 1.32 log(M* / 1e9) + 7.62` from the
photometric decompositions and `0.86 log(M* / 1e9) + 7.99` from the spectroscopic star formation histories, which sit
about 0.5 dex lower (the photometric masses may be biased high by the other components). Two verified fits that disagree:
the centre is drawn log-uniform between them (the widest range, as for the thick disc), both kept. The Milky Way's own
nuclear disc is not in either fit; they put it at 7e9 and 3e9 against the 1e9 to 1.5e9 measured."""
NUCLEAR_DISC_SCATTER_DEX = 0.3
"""PLACEHOLDER (gap plan 6): the scatter of a nuclear disc's mass about the fit, dex. The paper prints the parameters'
bootstrapped errors but not the scatter of the 16 points, which a figure holds; replaced when it is read."""
NUCLEATION_PEAK = 0.85
NUCLEATION_PEAK_MASS_SOLAR = 1.0e9
NUCLEATION_FALL_MASS_SOLAR = 1.0e11
NUCLEATION_FALL_FRACTION = 0.2
"""PLACEHOLDER (gap plan 6): the share of galaxies with a nuclear star cluster. [R3273] section 4.1 gives 80 percent
or more near 1e9 solar masses (and in over 300 late-type spirals) and a decline above it, without numbers: the share is
0.85 up to 1e9, falls in log mass to 0.2 at 1e11 and stays there. Replaced when the occupation fractions are read."""
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


def draw_circular_speed_kms(galaxy_seed, stellar_mass_solar):
    """The galaxy's circular speed, km/s: the baryonic Tully-Fisher relation's median with its scatter."""
    return circular_speed_kms(stellar_mass_solar) * 10.0 ** _scatter(_stream(galaxy_seed, "speed"), BTFR_V_SCATTER_DEX)


def thick_mass_ratio_bounds(vc_kms):
    """`(Yoachim and Dalcanton, Comeron)`: the two fits' thick-over-thin mass ratio at `vc_kms`."""
    yoachim = YOACHIM_MASS_RATIO_AT_100 * (vc_kms / 100.0) ** YOACHIM_MASS_RATIO_SLOPE
    comeron = 10.0 ** (COMERON_MASS_RATIO_INTERCEPT + COMERON_MASS_RATIO_SLOPE * math.log10(vc_kms))
    return yoachim, comeron


def thick_length_ratio_median(vc_kms):
    """[R0865]'s median thick-over-thin scale length ratio at `vc_kms`."""
    return 10.0 ** (COMERON_LENGTH_RATIO_INTERCEPT + COMERON_LENGTH_RATIO_SLOPE * math.log10(vc_kms))


def hot_component_bound(stellar_mass_solar):
    """The most the thick disk can weigh over the thin disk's, if the hot component (thick disk and classical bulge)
    is all thick disk: `M_h / (M* - M_h)` ([R0860: Eq. 1])."""
    hot = 10.0 ** (HOT_COMPONENT_INTERCEPT + HOT_COMPONENT_SLOPE * math.log10(stellar_mass_solar)) / stellar_mass_solar
    return hot / (1.0 - hot)


def halo_ratio_median(stellar_mass_solar):
    """The placeholder median of the stellar halo's mass over the galaxy's stellar mass."""
    return HALO_RATIO_AT_5E10 * (stellar_mass_solar / 5.0e10) ** HALO_RATIO_SLOPE


def nucleation_fraction(stellar_mass_solar):
    """The share of galaxies of `stellar_mass_solar` that hold a nuclear star cluster (a placeholder, see above)."""
    span = math.log10(NUCLEATION_FALL_MASS_SOLAR / NUCLEATION_PEAK_MASS_SOLAR)
    along = (math.log10(stellar_mass_solar / NUCLEATION_PEAK_MASS_SOLAR)) / span
    along = min(max(along, 0.0), 1.0)
    return NUCLEATION_PEAK + (NUCLEATION_FALL_FRACTION - NUCLEATION_PEAK) * along


def nuclear_cluster_mass_solar(stellar_mass_solar, slope, deviation):
    """A nuclear star cluster's mass, solar masses, in a galaxy of `stellar_mass_solar`: the power law of Neumayer and
    Georgiev (a slope between theirs), `deviation` standard deviations of 0.6 dex about it."""
    return 10.0 ** (NUCLEAR_CLUSTER_MASS_AT_1E9 + slope * math.log10(stellar_mass_solar / 1.0e9)
                    + deviation * NUCLEAR_CLUSTER_SCATTER_DEX)


def nuclear_disc_mass_bounds_solar(stellar_mass_solar):
    """`(photometric, spectroscopic)`: the two fits' nuclear disc mass, solar masses, in a galaxy of `stellar_mass_solar`."""
    return tuple(10.0 ** (intercept + slope * math.log10(stellar_mass_solar / 1.0e9))
                 for intercept, slope in NUCLEAR_DISC_MASS_AT_1E9)


def draw_nuclear_disc_mass_solar(galaxy_seed, stellar_mass_solar):
    """A nuclear stellar disc's mass, solar masses, in a galaxy of `stellar_mass_solar`: the centre log-uniform between the
    two fits (`NUCLEAR_DISC_MASS_AT_1E9`), scattered by `NUCLEAR_DISC_SCATTER_DEX`."""
    low, high = sorted(nuclear_disc_mass_bounds_solar(stellar_mass_solar))
    centre = _stream(galaxy_seed, "nuclear-disc-centre").uniform(math.log10(low), math.log10(high))
    return 10.0 ** (centre + _scatter(_stream(galaxy_seed, "nuclear-disc"), NUCLEAR_DISC_SCATTER_DEX))


def draw_structure(galaxy_seed, stellar_mass_solar=None, stellar_halo_share=None):
    """
    The structure the galaxy of `galaxy_seed` has: the same seed gives the same structure.

    Args:
        galaxy_seed (bytes or None): The galaxy's 16-byte seed; `None` uses the mass chain's fallback seed.
        stellar_mass_solar (float, optional): A stellar mass the plan pinned (GEN.603); the seed's own residuals are
            kept, so the same seed at a new mass has the same draws about the laws.
        stellar_halo_share (float, optional): A stellar halo share of the galaxy's stars the plan pinned (the Halo dial,
            GEN.606); it replaces the halo draw only. Refused for a galaxy with no halo.

    Returns:
        structure.Structure
    """
    mass = mass_chain.draw_mass_chain(galaxy_seed, stellar_mass_solar).stellar_mass_solar
    vc = draw_circular_speed_kms(galaxy_seed, mass)

    low, high = thick_mass_ratio_bounds(vc)
    centre = _stream(galaxy_seed, "thick-mass-centre").uniform(math.log10(min(low, high)), math.log10(max(low, high)))
    mass_ratio = 10.0 ** (centre + _scatter(_stream(galaxy_seed, "thick-mass"), COMERON_MASS_RATIO_SCATTER_DEX))
    length_ratio = thick_length_ratio_median(vc) * 10.0 ** _scatter(
        _stream(galaxy_seed, "thick-length"), COMERON_LENGTH_RATIO_SCATTER_DEX)
    height_centre = _stream(galaxy_seed, "thick-height-centre").uniform(
        math.log10(THICK_HEIGHT_RATIO_RANGE[0]), math.log10(THICK_HEIGHT_RATIO_RANGE[1]))
    height_ratio = 10.0 ** (height_centre + _scatter(_stream(galaxy_seed, "thick-height"), THICK_HEIGHT_RATIO_SCATTER_DEX))

    full = mass_chain.has_full_structure(mass)
    if not full:      # Boss, 7:58 PM EDT 10-21: one disc plus a bulge under the mass the thick disc and halo laws were fitted to
        mass_ratio = 0.0
    if stellar_halo_share is not None:
        if not full:
            raise ValueError(f"a galaxy of {mass:.3g} solar masses of stars has no stellar halo to pin "
                             f"(none under {dials.FULL_STRUCTURE_MIN_STELLAR_MASS_SOLAR:.3g})")
        halo = stellar_halo_share
    elif full:
        halo = halo_ratio_median(mass) * 10.0 ** _scatter(_stream(galaxy_seed, "halo"), HALO_RATIO_SCATTER_DEX)
        halo = min(max(halo, HALO_RATIO_RANGE[0]), HALO_RATIO_RANGE[1])
    else:
        halo = 0.0

    gas = GAS_DISPERSION_KMS[0] + _scatter(_stream(galaxy_seed, "gas-dispersion"), GAS_DISPERSION_KMS[1])
    centre_over_v, relative = DISC_CENTRAL_DISPERSION_OVER_V
    central = vc * centre_over_v * math.exp(_scatter(_stream(galaxy_seed, "disc-dispersion"), relative))
    exponent = DISC_AGE_EXPONENT[0] + _scatter(_stream(galaxy_seed, "disc-heating"), DISC_AGE_EXPONENT[1])
    bulge_scatter = _scatter(_stream(galaxy_seed, "bulge-dispersion"), BULGE_DISPERSION_SCATTER_DEX)
    bar = draw_bar(galaxy_seed, mass)
    if not full:
        bar = (False, *bar[1:])
    from planetgen.galaxy import shape_draw   # the stage is drawn there (a late import: shape_draw imports this module)
    classical = shape_draw.stage_of(galaxy_seed, stellar_mass_solar) == "early"
    black_hole = _scatter(_stream(galaxy_seed, "black-hole"), 1.0)
    if _stream(galaxy_seed, "nucleation").random() < nucleation_fraction(mass):
        slope = _stream(galaxy_seed, "nuclear-slope").uniform(*NUCLEAR_CLUSTER_SLOPE_RANGE)
        cluster = nuclear_cluster_mass_solar(mass, slope, _scatter(_stream(galaxy_seed, "nuclear-cluster"), 1.0))
    else:
        cluster = 0.0
    if not full:
        cluster = 0.0
    # Nuclear discs are the bar's (every TIMER galaxy, the sample the law is fitted to, is barred; GEN.282's ring needs
    # one too): a galaxy with no bar has none.
    disc = draw_nuclear_disc_mass_solar(galaxy_seed, mass) if bar[0] else 0.0
    black_holes = draw_black_hole_shares(galaxy_seed)
    kick = _kick_inputs(galaxy_seed, stellar_mass_solar)
    return Structure(
        vc_kms=vc, thick_mass_ratio=mass_ratio, thick_length_ratio=length_ratio, thick_height_ratio=height_ratio,
        bar_present=bar[0], bar_axis_ratio_y=bar[1], bar_axis_ratio_z=bar[2], bar_angle_deg=bar[3],
        halo_stellar_ratio=halo,
        disc_gas_dispersion_kms=gas, disc_central_dispersion_kms=max(central, gas * 1.05), disc_age_exponent=exponent,
        bulge_dispersion_scatter_dex=bulge_scatter,
        stellar_mass_solar=mass, classical_bulge=classical, black_hole_deviation=black_hole,
        nuclear_cluster_mass_solar=cluster, nuclear_disc_mass_solar=disc,
        black_hole_solar_share=black_holes[0], black_hole_metal_free_share=black_holes[1],
        neutron_star_unbound_share=natal_kick.neutron_star_unbound_share(*kick),
        neutron_star_ejecta_speeds_kms=natal_kick.ejecta_speeds_kms(*kick))


def _kick_inputs(galaxy_seed, stellar_mass_solar=None):
    """`(mass chain, circular speed, disc scale length in kpc)` the seed's natal kick is measured against (a pinned stellar
    mass replaces the seed's)."""
    from planetgen.galaxy import shape_draw   # shape_draw imports this module
    seed = galaxy_seed or mass_chain.FALLBACK_SEED
    chain = mass_chain.draw_mass_chain(seed, stellar_mass_solar)
    length_kpc = shape_draw.drawn_disk_scale_length_pc(seed, chain.stellar_mass_solar) / 1000.0
    return chain, draw_circular_speed_kms(seed, chain.stellar_mass_solar), length_kpc


def draw_neutron_star_unbound_share(galaxy_seed, stellar_mass_solar=None):
    """The share of the neutron stars the galaxy forms that its kicks eject (GEN.554, `galaxy.natal_kick`): the natal kick
    law against the escape speed of the seed's halo (none under the halo threshold), stars and disc, the disc's scale
    length the one the seed draws (a shape that fixes its own scale length keeps this share). A pinned stellar mass
    replaces the seed's."""
    return natal_kick.neutron_star_unbound_share(*_kick_inputs(galaxy_seed, stellar_mass_solar))


def draw_black_hole_shares(galaxy_seed):
    """
    The galaxy's `(solar share, metal-free share)` of the massive stars that fail to explode and leave a black hole
    (GEN.563, `physics.black_hole_formation.draw_failed_shares`): two ordered uniform draws, the solar one first, each over
    its widest verified range and cut at the neutron star budget of a galaxy of the Milky Way's masses
    (`black_hole_formation.share_cap` at the reference structure's unbound share, the same for every galaxy: a galaxy
    that ejects most of its neutron stars would otherwise have no room for black holes).
    """
    cap = black_hole_formation.share_cap(REFERENCE_STRUCTURE.neutron_star_unbound_share)
    return black_hole_formation.draw_failed_shares(
        _stream(galaxy_seed, "black-hole-solar").random(), _stream(galaxy_seed, "black-hole-metal-free").random(), cap)


BAR_FRACTION_COEFFICIENTS = (-82.2, 17.1, -0.88)
"""tuple: [R1288: Eq. 6] and Section 3: the fraction of S4G spirals with a bar (strong or weak) is
`1 / (1 + exp(-(a + b1 m + b2 m^2)))`, `m = log10(M*/Msun)`, fitted over m = 8 to 12."""
BAR_ELLIPTICITY = (0.53, 0.15)
BAR_ELLIPTICITY_RANGE = (0.15, 0.90)
"""tuple: [R1102: Table 3] and the VizieR table J/A+A/587/A160 (504 bars): the
maximum isophotal ellipticity of a bar, mean 0.534 and standard deviation 0.148, independent of the host's mass and
Hubble type; its range 0.10 to 0.92 cut at [R0045]'s smallest detected bar, 0.16. The model's `y0/x0`
is `1 - e` (the isodensity contour in the plane has that axis ratio): a STAND-IN for the bar body, which [R1486]
(median 0.63, SD 0.10) finds more elongated than its isophotes."""
BOXY_PEANUT_COEFFICIENTS = (-42.79, 4.13)
"""tuple: [R1290] Sect. 4.1 Table 2: the fraction of barred galaxies whose bar has a
boxy or peanut bulge, `1 / (1 + exp(-(a + b m)))` (half at m = 10.4)."""
BOXY_PEANUT_HEIGHT = (0.23, 0.06)
BOXY_PEANUT_HEIGHT_RANGE = (0.12, 0.38)
"""tuple: The boxy/peanut bar's `z0/x0`: STAND-IN, the peanut's half-length over the bar's (0.42 +/- 0.09, [R1030]) times its height over its length (0.5 to 0.6, [R0799]): 0.23, 25
percent scatter; no external-galaxy measurement of the ratio was found."""
THIN_BAR_HEIGHT = (0.071, 0.02)
THIN_BAR_HEIGHT_RANGE = (0.05, 0.15)
"""tuple: A bar without a peanut: thickness over length 1 / (14 +/- 4) in edge-on near-infrared images (Luetticke et al.
2000, astro-ph/0008478); the upper end is the simulated 0.1 to 0.2 [R2909]."""
UNBARRED_AXIS_RATIO_Y = (0.83, 0.21)
UNBARRED_AXIS_RATIO_Y_RANGE = (0.30, 1.0)
UNBARRED_AXIS_RATIO_Z = (0.54, 0.19)
UNBARRED_AXIS_RATIO_Z_RANGE = (0.19, 1.0)
"""tuple: The bulge of an unbarred spiral: B/A and C/A of the 41 unbarred CALIFA bulges of [R0902] Table C.1 (means and standard deviations computed from the table; the truncation to the table's range
and C/A not above B/A are a STAND-IN). [R3034] give a thicker bulge (C/A near 0.8); the newest source
is taken and the older kept in the record."""


def _clamped_normal(stream, centre, sigma, limits):
    return min(max(stream.gauss(centre, sigma), limits[0]), limits[1])


def bar_fraction(stellar_mass_solar):
    """The share of spirals of this stellar mass that have a bar ([R1288: Eq. 6]), the same for early and late
    stages (the S4G split, 0.63 against 0.66, is 1.1 sigma: [R1102: Table 1])."""
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
