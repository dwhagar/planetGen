# planetgen/galaxy/structure.py

"""
The structure of a galaxy that its shape row does not store (GEN.532, GEN.533): the thick disk, the bar, the stellar
halo and the motions of the populations, drawn from the galaxy seed by `galaxy.structure_draw` and rebuilt from it
every time the shape is loaded, so no column holds them (DB.66 stores the drawn values later).

`Structure` is one more field of `density.GalaxyShape`, so every reader of the density model sees the galaxy's own
numbers; `REFERENCE_STRUCTURE` is the structure of a shape built with no seed (a test, the math check): the
published structure of the Milky Way, a check target and not a default the generator draws from.

The populations' scale heights and dispersions follow ONE law, here, from the dispersion of the stars born with the
gas and the way a disc's dispersion grows with age:

- A thin-disc star of age `t` has a vertical dispersion `sigma_z(t)^2 = sigma_gas^2 + sigma_age^2 (t / tau)^(2 beta)`:
  the dispersion of the gas it was born from (molecular, 6.5 +/- 1 km/s inside the CO widths of 7.3 +/- 1.7 km/s, Mogotsi et al.
  2016, and 4.4 to 7.5 km/s, Sun et al. 2020, GEN.536; the atomic gas's is 10 +/- 2 km/s: [R4298]) plus a heating that grows as `t^beta` (`beta` = 0.441 +/- 0.007:
  [R4043: Table 2]; 0.5 +/- 0.1: [R2795]).
- `sigma_age` is fixed by the disc's central vertical dispersion, which is a fraction of the circular speed
  (`0.248 +/- 0.038 V`, [R2917], 30 spirals): the root mean square of `sigma_z` over
  a star formation flat in age from 0 to `tau` is that central value.
- The dispersion falls as `exp(-R / 2 h_R)` [R2917]; the model's profiles are at the reference radius
  (`tuning.REFERENCE_RADIUS_SCALE_LENGTHS` scale lengths), where every population is taken.
- A population's scale height goes as its dispersion squared (an isothermal sheet: `h = sigma^2 / (pi G Sigma)`, one
  surface density for all), so its ratio to the old disc's is `sigma_z(population)^2 / sigma_z(old)^2`.
- The profile's dispersion is the three-axis total: the gas part is isotropic (`3 sigma_gas^2`), the heated part has
  the radial dispersion `sigma_z / 0.54` ([R4043]: `sigma_z / sigma_R` = 0.54 at 10 Gyr, inside the
  0.5 to 0.6 of [R0201]) and the azimuthal `sigma_R / sqrt(2)` (a flat
  rotation curve). STAND-IN: one anisotropy for every age (Aumer and Binney give 0.41 at 1 Gyr to 0.57 at 10 Gyr).

The bulge's dispersion is the central one, from the circular speed and the bulge share: `Vc / sigma_0 = 2 (1 - B/T)`
for spirals ([R0908: Eq. 5] and the section after it), scattered
by 0.097 dex (the scatter of the relation about its fit read from the paper's Fig. 2, Gemini check, gap plan 7.2a);
the profile uses its three-axis total, `sqrt(3) sigma_0`. The stellar halo's is the Jeans equation for a tracer of
density slope `alpha` in a flat rotation curve with isotropic orbits, `sigma_r = Vc / sqrt(alpha)` ([R0908: Eq. 2]; [R4801]), three-axis total `sqrt(3) sigma_r`. Clusters, satellites and the nuclear
structures keep the one-axis `tuning.STELLAR_POPULATION_DISPERSION_KMS` (GEN.403, GEN.369, GEN.280), which the profile
states as a three-axis total like the others (GEN.284).

Every fill here is a labelled stand-in or a derivation (gap plan, section 6): see `structure_draw`.
"""

import math
from collections import namedtuple

from planetgen import tuning
from planetgen.galaxy import mass_chain, natal_kick
from planetgen.physics import black_hole_formation

Structure = namedtuple(
    "Structure",
    [
        "vc_kms",
        "thick_mass_ratio", "thick_length_ratio", "thick_height_ratio",
        "bar_present", "bar_axis_ratio_y", "bar_axis_ratio_z", "bar_angle_deg",
        "halo_stellar_ratio",
        "disc_gas_dispersion_kms", "disc_central_dispersion_kms", "disc_age_exponent",
        "bulge_dispersion_scatter_dex",
        "stellar_mass_solar", "classical_bulge", "black_hole_deviation", "nuclear_cluster_mass_solar",
        "nuclear_disc_mass_solar",
        "black_hole_solar_share", "black_hole_metal_free_share", "neutron_star_unbound_share",
        "neutron_star_ejecta_speeds_kms",
    ],
)
"""A galaxy's thick disk, bar, stellar halo and population motions (not stored; see the module docstring).

`vc_kms` is the circular speed; `thick_mass_ratio`, `thick_length_ratio` and `thick_height_ratio` are the thick disk's
mass, scale length and scale height over the thin disk's; `bar_present` says whether the bulge is a bar (a boxy
component along `bar_angle_deg` from the Sun-centre line with axes `bar_axis_ratio_y` and `bar_axis_ratio_z` of its
length) or a spheroid with those axis ratios; `halo_stellar_ratio` is the stellar halo's mass over the galaxy's stellar
mass; the three `disc_*` fields are the gas dispersion, the central vertical dispersion of the disc and the age
exponent of its heating; `bulge_dispersion_scatter_dex` is the bulge's dispersion draw about its law.

`stellar_mass_solar` is the galaxy's stellar mass; `classical_bulge` says whether its bulge is a classical one (the
black hole follows the bulge's mass) or not (it follows the stellar mass); `black_hole_deviation` is the central black
hole's draw about its law in standard deviations; `nuclear_cluster_mass_solar` is the nuclear star cluster's mass,
zero for a galaxy that has none (GEN.402); `nuclear_disc_mass_solar` is the nuclear stellar disc's mass, drawn from the galaxy's stellar mass
(GEN.628) and zero for a galaxy without a bar; the two `black_hole_*_share` are the shares of the massive stars that fail to explode,
at solar and at zero metallicity (GEN.563, `physics.black_hole_formation`); `neutron_star_unbound_share` is the share of the
neutron stars formed that the galaxy ejects (GEN.554, `galaxy.natal_kick`: the kick law against the escape speed), which the black hole
shares' cap and the neutron stars' bound count both read; `neutron_star_ejecta_speeds_kms` are the equal-share quantiles of the speed
those ejected neutron stars have far out (`natal_kick.ejecta_speeds_kms`), empty for a galaxy that ejects none."""

REFERENCE_RADIUS_OVER_TWO_LENGTHS = tuning.REFERENCE_RADIUS_SCALE_LENGTHS / 2.0
"""float: `R / 2 h_R` at the reference radius, the exponent of the fall of the dispersion."""
ANISOTROPY_SIGMA_Z_OVER_SIGMA_R = 0.54
"""float: The heated part's vertical over radial dispersion ([R4043]; a stand-in for every age)."""
AZIMUTHAL_OVER_RADIAL = math.sqrt(0.5)
"""float: `sigma_phi / sigma_R` for a flat rotation curve (the epicycle ratio `kappa / 2 Omega`)."""
HEATED_TOTAL_OVER_VERTICAL = math.sqrt(
    1.0 + (1.0 + AZIMUTHAL_OVER_RADIAL ** 2) / ANISOTROPY_SIGMA_Z_OVER_SIGMA_R ** 2)
"""float: The three-axis total of the heated part over its vertical dispersion (2.48)."""
STELLAR_HALO_DENSITY_SLOPE = (tuning.HALO_INNER_SLOPE + tuning.HALO_OUTER_SLOPE) / 2.0
"""float: The tracer slope the Jeans equation takes (3.5: the middle of the model's inner and outer slopes)."""

_MILKY_WAY_SOLAR_RADIUS_OVER_LENGTH = tuning.REFERENCE_RADIUS_SCALE_LENGTHS
_MILKY_WAY_THICK_LENGTH_RATIO = 2000.0 / 2600.0
_MILKY_WAY_THICK_HEIGHT_RATIO = 900.0 / 300.0
_MILKY_WAY_THICK_LOCAL_DENSITY_RATIO = 0.04


def _milky_way_thick_mass_ratio():
    """The thick disk's mass over the thin disk's of the Milky Way structure of BHG16 (a scale height of 900 pc and
    a scale length of 2.0 kpc against 300 pc and 2.6 kpc, 4% of the thin disk's density in the plane at the Sun): 18%
    (BHG16 give 6e9 against 3.5e10 solar masses)."""
    ratio = _MILKY_WAY_THICK_LENGTH_RATIO
    sun = _MILKY_WAY_SOLAR_RADIUS_OVER_LENGTH
    amplitude = _MILKY_WAY_THICK_LOCAL_DENSITY_RATIO * math.exp(sun / ratio - sun)
    return amplitude * ratio ** 2 * _MILKY_WAY_THICK_HEIGHT_RATIO


_MILKY_WAY_DISC_SCALE_LENGTH_KPC = 2.6
_MILKY_WAY_HALO_MASS_SOLAR = mass_chain.median_halo_mass(5.0e10)
_MILKY_WAY_CHAIN = mass_chain.MassChain(
    5.0e10, _MILKY_WAY_HALO_MASS_SOLAR, mass_chain.virial_radius_kpc(_MILKY_WAY_HALO_MASS_SOLAR),
    mass_chain.median_concentration(_MILKY_WAY_HALO_MASS_SOLAR))
_MILKY_WAY_NEUTRON_STAR_UNBOUND_SHARE = natal_kick.neutron_star_unbound_share(
    _MILKY_WAY_CHAIN, 220.0, _MILKY_WAY_DISC_SCALE_LENGTH_KPC)
_REFERENCE_BLACK_HOLE_SHARES = black_hole_formation.reference_failed_shares(_MILKY_WAY_NEUTRON_STAR_UNBOUND_SHARE)

REFERENCE_STRUCTURE = Structure(
    vc_kms=220.0,
    thick_mass_ratio=_milky_way_thick_mass_ratio(),
    thick_length_ratio=_MILKY_WAY_THICK_LENGTH_RATIO,
    thick_height_ratio=_MILKY_WAY_THICK_HEIGHT_RATIO,
    bar_present=True, bar_axis_ratio_y=620.0 / 1580.0, bar_axis_ratio_z=430.0 / 1580.0, bar_angle_deg=27.0,
    halo_stellar_ratio=1.4e9 / 5.0e10,
    disc_gas_dispersion_kms=6.5, disc_central_dispersion_kms=0.248 * 220.0, disc_age_exponent=0.441,
    bulge_dispersion_scatter_dex=0.0,
    stellar_mass_solar=5.0e10, classical_bulge=False, black_hole_deviation=0.0, nuclear_cluster_mass_solar=2.5e7,
    nuclear_disc_mass_solar=1.5e9,
    black_hole_solar_share=_REFERENCE_BLACK_HOLE_SHARES[0], black_hole_metal_free_share=_REFERENCE_BLACK_HOLE_SHARES[1],
    neutron_star_unbound_share=_MILKY_WAY_NEUTRON_STAR_UNBOUND_SHARE,
    neutron_star_ejecta_speeds_kms=natal_kick.ejecta_speeds_kms(_MILKY_WAY_CHAIN, 220.0, _MILKY_WAY_DISC_SCALE_LENGTH_KPC),
)
"""Structure: The Milky Way's published structure, the check target of the math check and the structure of a shape
built with no seed. The thick disk and the bar are the [R0422] (ARA&A 54:529) and Dwek et al.
1995 (ApJ 445:716, fit G2) numbers, the bar's angle [R4649] (MNRAS 435:1874) 27 degrees from the
Sun-centre line (BHG16: 28 to 33), the halo 1.4e9 solar masses of the 5e10 of stars [R1026], the disc's motions the laws above at `Vc` = 220 km/s, the neutron stars' unbound share the kick law's at
the median halo of 5e10 solar masses of stars and a 2.6 kpc disc, and the black hole shares the middle of the draws (a
stand-in: no measurement gives a galaxy's)."""


def stellar_halo_system_fraction(structure):
    """The stellar halo's systems over the visible components' systems: its mass over the rest of the galaxy's stars
    (the halo is `halo_stellar_ratio` of them all, GEN.404)."""
    return structure.halo_stellar_ratio / (1.0 - structure.halo_stellar_ratio)


def _heated_square(structure, name):
    """The squared vertical dispersion the heating gives population `name` at the solar-like radius, km^2 s^-2 (the
    root mean square over the population's ages of `sigma_age^2 (t / tau)^(2 beta)` times the radial fall squared)."""
    low, high = tuning.THIN_DISC_OLD_WINDOW_GY if name == "old" else tuning.STELLAR_POPULATION_AGE_RANGES_GY[name]    # GEN.1025: the thin disc's
    tau = tuning.STAR_FORMATION_AGE_RANGE_GY[1]
    power = 1.0 + 2.0 * structure.disc_age_exponent
    at_tau = (structure.disc_central_dispersion_kms ** 2 - structure.disc_gas_dispersion_kms ** 2) * power
    mean = (high ** power - low ** power) / (power * (high - low)) / tau ** (2.0 * structure.disc_age_exponent)
    return at_tau * mean * math.exp(-2.0 * REFERENCE_RADIUS_OVER_TWO_LENGTHS)


def disc_vertical_dispersion_kms(structure, name):
    """The vertical dispersion of the thin-disc population `name` ("young", "intermediate" or "old") at the
    solar-like radius, km/s."""
    return math.sqrt(structure.disc_gas_dispersion_kms ** 2 + _heated_square(structure, name))


def population_height_ratios(structure):
    """
    The scale height of each thin-disc population over the old population's (the old one is 1): the squared vertical
    dispersions' ratio.

    Returns:
        dict: "young", "intermediate", "old" -> ratio.
    """
    sigma = {name: disc_vertical_dispersion_kms(structure, name) ** 2 for name in ("young", "intermediate", "old")}
    return {name: value / sigma["old"] for name, value in sigma.items()}


def disc_dispersion_kms(structure, name):
    """The three-axis total dispersion of the thin-disc population `name` at the solar-like radius, km/s."""
    return math.sqrt(3.0 * structure.disc_gas_dispersion_kms ** 2
                     + HEATED_TOTAL_OVER_VERTICAL ** 2 * _heated_square(structure, name))


def halo_dispersion_kms(structure):
    """The stellar halo's three-axis total dispersion: `sqrt(3) Vc / sqrt(alpha)` (the Jeans equation, flat curve,
    isotropic orbits)."""
    return math.sqrt(3.0 / STELLAR_HALO_DENSITY_SLOPE) * structure.vc_kms


def bulge_dispersion_kms(structure, bulge_to_total):
    """The bulge's three-axis total dispersion: `sqrt(3) sigma_0`, `sigma_0 = Vc / (2 (1 - B/T))` times the draw's
    scatter [R0908]."""
    sigma_0 = structure.vc_kms / (2.0 * (1.0 - bulge_to_total)) * 10.0 ** structure.bulge_dispersion_scatter_dex
    return math.sqrt(3.0) * sigma_0


BLACK_HOLE_BULGE_MASS_AT_1E11 = 0.49e9
"""float: The black hole's mass, solar masses, in a classical bulge of 1e11 solar masses (Kormendy and Ho 2013,
[R2456: Eq. 10]: `M_BH / 1e9 = 0.49 (M_bulge / 1e11)^1.16`)."""
BLACK_HOLE_BULGE_SLOPE = 1.16
BLACK_HOLE_BULGE_SCATTER_DEX = 0.29
"""float: The scatter of that relation, dex (Kormendy and Ho 2013, classical bulges and ellipticals)."""
BLACK_HOLE_STELLAR_INTERCEPT = 7.45
BLACK_HOLE_STELLAR_SLOPE = 1.05
BLACK_HOLE_STELLAR_SCATTER_DEX = 0.24
"""float: The relation of a galaxy without a classical bulge, `log M_BH = 7.45 + 1.05 log(M* / 1e11)`, intrinsic
scatter 0.24 dex ([R3723: Eq. 5], broad-line active nuclei: virial masses, a different population
from the bulge relation's dynamical ones)."""
M_SIGMA_INTERCEPT = -0.510
M_SIGMA_SLOPE = 4.377
"""float: `log(M_BH / 1e9) = -0.510 + 4.377 log(sigma / 200 km/s)` (Kormendy and Ho 2013, Eq. 3)."""


def black_hole_mass_solar(structure, bulge_to_total):
    """
    The galaxy's central black hole, solar masses: the mass in a classical bulge's `bulge_to_total` share of the stars
    (Kormendy and Ho 2013), or the stellar mass without one [R3723], times the galaxy's draw about
    the law.
    """
    if structure.classical_bulge:
        bulge = max(bulge_to_total, 1.0e-6) * structure.stellar_mass_solar
        mass = BLACK_HOLE_BULGE_MASS_AT_1E11 * (bulge / 1.0e11) ** BLACK_HOLE_BULGE_SLOPE
        scatter = BLACK_HOLE_BULGE_SCATTER_DEX
    else:
        mass = 10.0 ** (BLACK_HOLE_STELLAR_INTERCEPT
                        + BLACK_HOLE_STELLAR_SLOPE * math.log10(structure.stellar_mass_solar / 1.0e11))
        scatter = BLACK_HOLE_STELLAR_SCATTER_DEX
    return mass * 10.0 ** (structure.black_hole_deviation * scatter)


def black_hole_dispersion_kms(mass_solar):
    """The bulge velocity dispersion, km/s, of a galaxy whose central black hole weighs `mass_solar`: the M-sigma
    relation inverted (Kormendy and Ho 2013, Eq. 3; its scatter is the hole's own draw)."""
    return 200.0 * 10.0 ** ((math.log10(mass_solar / 1.0e9) - M_SIGMA_INTERCEPT) / M_SIGMA_SLOPE)
