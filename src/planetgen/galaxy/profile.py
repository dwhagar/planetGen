# planetgen/galaxy/profile.py

"""
The Sector Profile (GEN.235)
============================

The program has the same rules everywhere in the sphere, and each sector carries the parameters that tune those
rules to where it is (docs/design/galaxy/intergalactic-and-halo-contents.md section 5). This module computes those
parameters: `sector_profile(address, skeleton)` for a grid cell, `profile_at(position_pc, ...)` for a point.

It is a pure function of the galaxy's shape, its sector edge, the address and (for the superbubble a place lies in) the
galaxy's seed. It reads no database, draws no random number of the run's and is closed form, so the same code answers for an uncharted preview, a scatter and a fill, and two calls
give the same profile. Nothing here generates anything: a profile of a cell nobody asked for is computed, shown and
dropped.

A profile holds:

* `region`: the term of the density that dominates (`nuclear`, `bulge`, `disk`, `halo`, `cluster` or `satellite`), a label for pages and
  not a rule;
* `density` and `expected_systems`: today's `relative_density` and the expected systems per cell, `E x density`,
  which is the mean of the sector's Poisson draw;
* `mix`: the density by population (`young`, `intermediate`, `old`, `bulge`, `halo`, `cluster`, `satellite`,
  `nuclear_cluster`, `nuclear_disc`), summing to `density`;
* `populations`: for each, its age window, [Fe/H] distribution, velocity dispersion and what its net motion follows;
* `gas`: the cold-gas fraction against the solar circle and the ambient medium with its density, and the superbubble
  (GEN.254) the place lies in, if any;
* `dark_matter`: the NFW halo's density at the cell and the mass in it.

* `kind_factors` (GEN.236): the factor by which the place scales the rate per star of each rare kind
  (`kind_factors_at`).
* `ejecta` (GEN.238): the neutron stars and black holes a cell holds that were ejected from the galaxy, on average
  (`ejecta_field.expected_in_cell`), which the fill adds to the cell's draw.

Later items add to it (the components, GEN.237; the scatter-limit tier, GEN.279): each adds a field and its source, and
none changes a rule.
"""

import functools
import math
from dataclasses import dataclass, field

import numpy as np

from planetgen import tuning
from planetgen.galaxy import remnant_density as remnantDensity
from planetgen.galaxy import component_field as componentField
from planetgen.galaxy import density as galaxyDensity
from planetgen.galaxy import ejecta_field, nebula_field, spherical, superbubble_field
from planetgen.galaxy import structure as galaxyStructure
from planetgen.physics import constants, metallicity
from planetgen.physics.units import pc_to_ly

POPULATIONS = ("young", "intermediate", "old", "bulge", "halo", "cluster", "satellite", "nuclear_cluster",
               "nuclear_disc", "satellite_young", "satellite_intermediate", *tuning.OPEN_CLUSTER_POPULATIONS)
"""tuple: The populations of a profile's mix, in the order `galaxy.density.population_densities` has them."""

_DISK = ("young", "intermediate", "old", *tuning.OPEN_CLUSTER_POPULATIONS)
"""tuple: The populations that make up the disk's share of a region: a young open cluster lies in the disk."""

_YOUNG_OPEN_CLUSTER_POPULATIONS = tuple(
    name for name in tuning.OPEN_CLUSTER_POPULATIONS
    if tuning.STELLAR_POPULATION_AGE_RANGES_GY[name][1] <= tuning.STELLAR_POPULATION_AGE_RANGES_GY["young"][1])
"""tuple: The open clusters' populations whose stars are all younger than the disk's `young` window ends: they count as
young stars where the young share matters (the supernova remnants, the molecular clouds)."""


@dataclass(frozen=True)
class PopulationProfile:
    """
    One stellar population at a place.

    Attributes:
        name (str): One of `POPULATIONS`.
        density (float): Its share of the cell's `relative_density`.
        share (float): `density` over the cell's total (0 when the cell holds none).
        age_gy (tuple): `(youngest, oldest)` age of its stars, Gyr.
        fe_h (tuple): `((weight, mean, sigma), ...)`, the [Fe/H] distribution of its stars at this place, dex
            (one peak; two for a disk or halo place that two components reach, one per component for a cluster or
            satellite place).
        dispersion_kms (float): Velocity dispersion, km/s: the three-axis total for every population (the disc, bulge and
            halo from the galaxy's own structure; the others the one-axis values of `tuning` times
            `constants.ISOTROPIC_TOTAL_OVER_ONE_AXIS`). `place_terms.from_profile` reads the one-axis value back.
        rotation (str): What its net motion follows: `disk`, `bar`, `orbit` or `none`.
    """
    name: str
    density: float
    share: float
    age_gy: tuple
    fe_h: tuple
    dispersion_kms: float
    rotation: str


@dataclass(frozen=True)
class GasProfile:
    """
    The gas in a cell.

    Attributes:
        cold_fraction (float): Molecular gas here against the solar circle's average, capped at 1: 0 in the
            halo, 1 where the disk holds at least as much as the Sun's neighbourhood.
        ambient (str): The densest ambient medium: `disk`, `superbubble`, `hot-halo` or `intergalactic`.
        ambient_density_cm3 (float): Its density, H atoms per cm3.
        mass_solar (float): The ambient gas in the cell, solar masses.
        superbubble (Superbubble or None): The cavity the place lies in (GEN.254; `superbubble_field`): there the
            disk's gas is its thin, ionised gas, `cold_fraction` is 0 and no molecular cloud is kept. `None` outside
            one, or when no seed was given.
    """
    cold_fraction: float
    ambient: str
    ambient_density_cm3: float
    mass_solar: float
    superbubble: object = None


@dataclass(frozen=True)
class DarkMatter:
    """The dark-matter halo at a cell: `density_solar_pc3` at its centre and `mass_solar` in the cell. Gravity only."""
    density_solar_pc3: float
    mass_solar: float


@dataclass(frozen=True)
class SectorProfile:
    """
    The parameters of one place. See the module docstring.

    Attributes:
        address (tuple or None): `(shell, band, wedge)`; `None` for a point.
        position_pc (tuple): The centre, galaxy-frame parsecs.
        edge_pc (float): The sector edge the profile was made at.
        volume_pc3 (float): The cell's volume.
        shape (GalaxyShape): The galaxy's shape (the fill's molecular cloud field reads its gas from it).
        region (str): `nuclear`, `bulge`, `disk`, `halo`, `cluster` or `satellite`.
        density (float): `relative_density` here, never negative.
        expected_systems (float): The mean of the cell's Poisson draw, `E x density` (`E` the systems a cell
            holds at density 1).
        mix (dict): Population -> density (each `>= 0`, adding up to `density` except where a shape's negative
            term was counted as none).
        populations (dict): Population -> `PopulationProfile`.
        gas (GasProfile): The cell's gas.
        dark_matter (DarkMatter): The cell's dark matter.
        kind_factors (dict): Kind -> the factor the place applies to that kind's rate per star (`kind_factors_at`);
            a kind not in it is not scaled by place.
        type_ia_share (float): The share of the place's supernova remnants that are Type Ia (`type_ia_share_at`): 1
            where no young population is, `tuning.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE` against the solar circle's
            mix.
        ejecta (dict): Kind -> the mean number of ejected objects of that kind in the cell (`ejecta_field`, GEN.238):
            neutron stars and black holes the galaxy lost to a kick, which stream through every cell at a density
            that falls as 1 / r^2. Not part of the stars' rate: the fill adds it to the cell's own Poisson draw.
    """
    address: tuple | None
    position_pc: tuple
    edge_pc: float
    volume_pc3: float
    shape: object
    region: str
    density: float
    expected_systems: float
    mix: dict
    populations: dict
    gas: GasProfile
    dark_matter: DarkMatter
    kind_factors: dict = field(default_factory=dict)
    type_ia_share: float = tuning.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE
    ejecta: dict = field(default_factory=dict)

    @property
    def center_pc(self):
        """The centre, as the fill reads it."""
        return self.position_pc


def _region(mix):
    """The label of the term that dominates `mix`: the disk's three populations count together."""
    terms = {"disk": sum(mix[name] for name in _DISK), "bulge": mix["bulge"], "halo": mix["halo"],
             "cluster": mix["cluster"],
             "satellite": mix["satellite"] + mix["satellite_young"] + mix["satellite_intermediate"],
             "nuclear": mix["nuclear_cluster"] + mix["nuclear_disc"]}
    best = max(terms.values())
    if best <= 0.0:
        return "halo"
    for name in ("nuclear", "cluster", "satellite", "bulge", "disk", "halo"):    # a tie goes to the more specific place
        if terms[name] >= best:
            return name
    return "halo"


def _dispersion_kms(name, shape):
    """
    The velocity dispersion of population `name` for this galaxy, km/s (three-axis total): the thin-disc populations',
    the bulge's and the halo's follow from the galaxy's own structure (`galaxy.structure`, GEN.532); the others are
    the one-axis `tuning.STELLAR_POPULATION_DISPERSION_KMS` times `constants.ISOTROPIC_TOTAL_OVER_ONE_AXIS` (GEN.284).
    """
    structure = shape.structure
    if name in ("young", "intermediate", "old") or name in tuning.OPEN_CLUSTER_POPULATIONS:
        return galaxyStructure.disc_dispersion_kms(structure, name)       # a young cluster's bin is heated as its age
    if name == "bulge":
        return galaxyStructure.bulge_dispersion_kms(structure, galaxyDensity.bulge_to_total(shape))
    if name == "halo":
        return galaxyStructure.halo_dispersion_kms(structure)
    return tuning.STELLAR_POPULATION_DISPERSION_KMS[name] * constants.ISOTROPIC_TOTAL_OVER_ONE_AXIS


def _fe_h(name, distance_ly, shape=None, position_pc=None):
    """`((weight, mean, sigma), ...)` of a population's [Fe/H] `distance_ly` from the centre. A cluster or satellite
    place takes the [Fe/H] of the components that reach it (`ComponentField.fe_h_peaks`) and, where it has none, the
    generic peaks of `tuning`."""
    if name in componentField.POPULATIONS_IN_THE_GATE:
        found = () if shape is None or shape.components is None else shape.components.fe_h_peaks(position_pc, name)
        if found:
            return found
        if name in tuning.OPEN_CLUSTER_POPULATIONS:      # a place no cluster reaches: the disk's stars of that age
            name = metallicity.open_cluster_disk_population(name)
        else:
            return metallicity.generic_peaks(name)
    mean, sigma = metallicity.mean_fe_h(name, distance_ly)
    return ((1.0, mean, sigma),)


@functools.lru_cache(maxsize=8)
def _reference_young_share(shape):
    """The young stars' share of the stars on the reference ring on the plane (the solar circle, where the rates of the
    gas-born kinds are quoted), averaged around it, arms and gaps alike."""
    radius = nebula_field.REFERENCE_RADIUS_SCALE_LENGTHS * shape.disk_scale_length_pc
    shares = []
    for step in range(36):
        angle = 2.0 * math.pi * step / 36
        mix = galaxyDensity.population_densities((radius * math.cos(angle), radius * math.sin(angle), 0.0), shape)
        total = sum(max(value, 0.0) for value in mix.values())
        shares.append(max(mix["young"], 0.0) / total if total > 0.0 else 0.0)
    return sum(shares) / len(shares)


def _young_and_old_against_the_solar_circle(populations, shape):
    """`(young, old)`, the young stars' share of the place's stars and the share that is not young, each divided by
    the same share around the solar circle (`_reference_young_share`): 1 and 1 for the solar circle's mix."""
    reference = _reference_young_share(shape)
    young_share = (populations["young"].share + populations["satellite_young"].share
                   + sum(populations[name].share for name in _YOUNG_OPEN_CLUSTER_POPULATIONS))
    young = 0.0 if reference <= 0.0 else young_share / reference
    old = (1.0 - young_share) / (1.0 - reference) if reference < 1.0 else 0.0
    return young, old


def planetary_nebula_factor_array(x, y, z, shape):
    """
    `kind_factors_at`'s `planetary-nebula` factor at arrays of points (GEN.528): the share of the stars that are not
    young against the solar circle's (`_reference_young_share`), capped at 1, from the continuous populations' array
    form. The phenomenon scatters weigh a longitude bin by it, so a place where the young are most of the stars (an
    arm's crest, the inner disk's thin layer) makes fewer planetary nebulae than an old one, as the sector fill does.
    """
    mix = galaxyDensity.population_densities_array(x, y, z, shape)
    total = sum(np.maximum(value, 0.0) for value in mix.values())
    young = np.maximum(mix["young"], 0.0) + np.maximum(mix["satellite_young"], 0.0)
    share = np.where(total > 0.0, young / np.where(total > 0.0, total, 1.0), 0.0)
    reference = _reference_young_share(shape)
    if reference >= 1.0:
        return np.zeros_like(share)
    return np.minimum(1.0, (1.0 - share) / (1.0 - reference))


def supernova_progenitor_factor(young_share, shape):
    """
    The factor on the supernova remnants of stars of which `young_share` are young, against the rate per star at the solar circle (GEN.1191,
    GEN.1198): a core-collapse remnant comes from a young star and a Type Ia from one that is not, so it is the solar circle's odds of each
    (`tuning.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE`, 3 to 1) times the young share and the old share against the solar circle's. 1 for the solar
    circle's own mix, above 1 where stars form faster per star than the disk (a Cloud), 0 for stars with none young: they form none now, so
    they hold none of the gas a remnant needs to be seen in. Not capped.
    """
    if young_share <= 0.0:
        return 0.0
    reference = _reference_young_share(shape)
    chance = tuning.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE
    young_ratio = 0.0 if reference <= 0.0 else young_share / reference
    old_ratio = (1.0 - young_share) / (1.0 - reference) if reference < 1.0 else 0.0
    return (1.0 - chance) * young_ratio + chance * old_ratio


def type_ia_share_at(populations, shape):
    """
    The share of a place's supernova remnants that are Type Ia (GEN.319). Core-collapse remnants come from the young
    stars and Type Ia from the stars that are not young, so the odds of a Type Ia remnant are the solar circle's odds
    (`tuning.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE`, three core collapses to a Type Ia) times the place's old share over
    its young share, each against the solar circle's. A place with no young stars has only Type Ia (a halo, a cluster,
    an old satellite: no neutron star or black hole in the remnant), one with only young stars has none.

    Args:
        populations (dict): Population -> `PopulationProfile`.
        shape (GalaxyShape): The galaxy's shape.

    Returns:
        float: In 0 to 1.
    """
    young, old = _young_and_old_against_the_solar_circle(populations, shape)
    chance = tuning.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE
    type_ia, core_collapse = chance * old, (1.0 - chance) * young
    total = type_ia + core_collapse
    if total <= 0.0:
        return 1.0
    return type_ia / total


_SATELLITE_POPULATIONS = ("satellite_young", "satellite_intermediate", "satellite")


def _supernova_remnant_factor(populations, gas_fraction, shape):
    """The `supernova-remnant` factor of a place: the disk's gas fraction (0 in clusters) for the share of its stars that are not a
    satellite's, and the satellite's own progenitor factor for the rest (GEN.1198)."""
    satellite = sum(populations[name].share for name in _SATELLITE_POPULATIONS)
    others = gas_fraction * (1.0 - populations["cluster"].share)
    if satellite <= 0.0:
        return others
    own = supernova_progenitor_factor(populations["satellite_young"].share / satellite, shape)
    return (1.0 - satellite) * others + satellite * own


def kind_factors_at(populations, gas_fraction, cold_fraction, shape):
    """
    The factors by which a place scales the rate per star of each kind that is not a plain per-star rate (GEN.236;
    docs/design/galaxy/intergalactic-and-halo-contents.md section 5). The rates themselves stay per star
    (`tuning.phenomenon_rate_per_star`); a kind absent here is not scaled.

    * **1 everywhere**: `brown-dwarf` (every population makes it).
    * **`black-hole`** (GEN.563): the share of the place's massive stars that failed to explode and have made their
      remnant (`remnant_density.place_factor_from_peaks`): the populations' shares, ages and [Fe/H] peaks, so a
      metal-poor place makes more. Its rate per star is `remnant_density.base_rate_per_star()`. The 9,000 L_sun
      floor under the bright stars is not a rate and stays. A runaway flag is a speed against the neighbours, so it
      needs no factor either.
    * **`neutron-star`** (GEN.554): what the black holes leave of the massive stars that have made their remnant, less
      the galaxy's unbound share (`remnant_density.place_factor_from_peaks`, `kind="neutron-star"`): a metal-rich place
      makes more. A bound remnant is born where massive stars were and kicked into the disk's thickness, so it follows
      the stars; the kick's thickening belongs to the disk populations' motion, not to the rate.
    * **`supernova-remnant`**: the gas fraction (the remnant of a core collapse shines only where the shell meets gas),
      0 in the halo and in clusters. Not zeroed in a superbubble: the hot cavity is where the explosions went off. The
      stars of a satellite are not the disk's (GEN.1198): their share of the place takes the satellite's factor
      (`supernova_progenitor_factor` of their own young share, uncapped, as the component pass has it, GEN.1191)
      in place of the gas fraction, so a Cloud's cell holds remnants at the rate its own star formation gives.
    * **`planetary-nebula`** (GEN.319): the share of the place's stars that are not young against the solar circle's,
      capped at 1: the young are too young to have made one, every older population makes them, so 1 in the disk
      (within 3%) and in every older place.
    * **`molecular-cloud`**: the cold gas fraction times the young share against the solar circle's, so the disk's
      young, gas-rich places carry them and the old and halo ones none.
    * **`nebula`** (the emission and reflection nebulae grown around a sector's own hot stars): the cold gas fraction.
      The hot host star's own age is the young share here, so it is not counted twice.
    * **`rogue-planet` and `comet`**: the planet occurrence at the place's [Fe/H], against the disk's
      (`metallicity.planet_occurrence_factor`), averaged over the stars' [Fe/H] distributions.

    Args:
        populations (dict): Population -> `PopulationProfile`.
        gas_fraction (float): The molecular gas against the solar circle's, capped at 1, a superbubble ignored.
        cold_fraction (float): The same, 0 inside a superbubble (`GasProfile.cold_fraction`).
        shape (GalaxyShape): The galaxy's shape.

    Returns:
        dict: Kind -> factor, each in 0 to 1 but a satellite's supernova remnants, which scale with its star formation.
    """
    total = sum(pop.share for pop in populations.values())
    occurrence = 1.0
    if total > 0.0:
        occurrence = sum(pop.share * metallicity.expected_factor(metallicity.planet_occurrence_factor, pop.fe_h)
                         for pop in populations.values()) / total
    young, old = _young_and_old_against_the_solar_circle(populations, shape)
    young = min(1.0, young)
    return {
        "neutron-star": remnantDensity.place_factor_from_peaks(populations, shape.structure, "neutron-star"),
        "black-hole": remnantDensity.place_factor_from_peaks(populations, shape.structure),
        "planetary-nebula": min(1.0, old), "brown-dwarf": 1.0,
        "supernova-remnant": _supernova_remnant_factor(populations, gas_fraction, shape),
        "molecular-cloud": cold_fraction * young,
        "nebula": cold_fraction,
        "rogue-planet": occurrence, "comet": occurrence,
    }


def _scale(shape):
    """How much larger than the Milky Way's the galaxy is, by its disk scale length."""
    return shape.disk_scale_length_pc / tuning.HALO_SCALE_LENGTH_PC


def dark_matter_density(radius_pc, shape):
    """
    The NFW halo's density, solar masses per pc3, at `radius_pc` from the centre: `M_h / (4 pi r_h^3) / (u (1 + u)^2)`,
    `u = r / r_h`, with the Milky Way's `M_h` and `r_h` (`tuning`) scaled by the galaxy's disk scale length.
    """
    scale = _scale(shape)
    r_h = tuning.DARK_MATTER_SCALE_RADIUS_PC * scale
    u = max(radius_pc, 1.0) / r_h
    return tuning.DARK_MATTER_SCALE_MASS_SOLAR * scale / (4.0 * math.pi * r_h ** 3) / (u * (1.0 + u) ** 2)


def _ambient(position_pc, shape, cold_fraction, superbubble=None):
    """`(medium, density per cm3)` of the densest of the disk's ISM (a superbubble's gas inside one), the hot halo and
    the intergalactic medium."""
    scale = _scale(shape)
    radius = math.sqrt(sum(c * c for c in position_pc))
    core = tuning.HOT_HALO_CORE_PC * scale
    hot = tuning.HOT_HALO_DENSITY_CM3 * (max(radius, core) / (tuning.HOT_HALO_REFERENCE_PC * scale)) ** -tuning.HOT_HALO_SLOPE
    media = {
        "disk": tuning.DISK_ISM_DENSITY_CM3 * cold_fraction,
        "hot-halo": hot,
        "intergalactic": tuning.COSMIC_MEAN_DENSITY_CM3 * tuning.LOCAL_GROUP_GAS_OVERDENSITY,
    }
    if superbubble is not None:
        media = {"superbubble": superbubble.ambient_density_cm3, "hot-halo": hot, "intergalactic": media["intergalactic"]}
    name = max(media, key=media.get)
    return name, media[name]


def profile_at(position_pc, shape, edge_pc, expected_at_density_1, address=None, volume_pc3=None, galaxy_seed=None):
    """
    The profile of the place at `position_pc`.

    Args:
        position_pc (tuple): `(x, y, z)`, galaxy-frame parsecs.
        shape (GalaxyShape): The galaxy's shape.
        edge_pc (float): The sector edge, parsecs.
        expected_at_density_1 (float): The systems a sector holds at `relative_density` 1
            (`galaxy.skeleton.expected_system_count_at_density_1`; stored as the skeleton's
            `expected_system_count_at_density_1`).
        address (tuple, optional): The cell's `(shell, band, wedge)`, for the record.
        volume_pc3 (float, optional): The cell's volume; a cube of `edge_pc` when not given.
        galaxy_seed (bytes, optional): The galaxy's seed, for the superbubble field (GEN.254); without it the place
            is taken to lie in none.

    Returns:
        SectorProfile
    """
    position_pc = tuple(float(c) for c in position_pc)
    volume = edge_pc ** 3 if volume_pc3 is None else volume_pc3
    mix = {name: max(value, 0.0) for name, value in galaxyDensity.population_densities(position_pc, shape).items()}
    density = max(galaxyDensity.relative_density(position_pc, shape), 0.0)
    total = sum(mix.values())
    distance_ly = pc_to_ly(math.sqrt(sum(c * c for c in position_pc)))
    ages = tuning.STELLAR_POPULATION_AGE_RANGES_GY
    populations = {
        name: PopulationProfile(
            name=name, density=mix[name], share=mix[name] / total if total > 0.0 else 0.0, age_gy=tuple(ages[name]),
            fe_h=_fe_h(name, distance_ly, shape, position_pc), dispersion_kms=_dispersion_kms(name, shape),
            rotation=tuning.STELLAR_POPULATION_ROTATION[name])
        for name in POPULATIONS}
    bubble = None if galaxy_seed is None else superbubble_field.containing(galaxy_seed, shape, position_pc)
    gas_fraction = min(1.0, max(nebula_field.gas_factor(position_pc, shape), 0.0))
    cold = 0.0 if bubble is not None else gas_fraction
    medium, n_cm3 = _ambient(position_pc, shape, cold, bubble)
    gas = GasProfile(cold_fraction=cold, ambient=medium, ambient_density_cm3=n_cm3,
                     mass_solar=n_cm3 * tuning.GAS_SOLAR_MASS_PER_PC3_PER_CM3 * volume, superbubble=bubble)
    rho = dark_matter_density(math.sqrt(sum(c * c for c in position_pc)), shape)
    return SectorProfile(
        address=None if address is None else tuple(address), position_pc=position_pc, edge_pc=edge_pc,
        volume_pc3=volume, shape=shape, region=_region(mix), density=density,
        expected_systems=expected_at_density_1 * density, mix=mix, populations=populations, gas=gas,
        dark_matter=DarkMatter(density_solar_pc3=rho, mass_solar=rho * volume),
        kind_factors=kind_factors_at(populations, gas_fraction, cold, shape),
        type_ia_share=type_ia_share_at(populations, shape),
        ejecta=ejecta_field.expected_in_cell(position_pc, volume, shape, expected_at_density_1, edge_pc))


def sector_profile(address, skeleton):
    """
    The profile of the grid cell `address`.

    Args:
        address (tuple): `(shell, band, wedge)`.
        skeleton (store.GalaxySkeletonInfo): The galaxy's stored skeleton (its `shape`, `edge_pc` and
            `expected_system_count_at_density_1`).

    Raises:
        ValueError: If the address is not a cell of the grid.
    """
    shell, band, wedge = address
    position = spherical.sector_position_pc(shell, band, wedge, skeleton.edge_pc)
    return profile_at(position, skeleton.shape, skeleton.edge_pc, skeleton.expected_system_count_at_density_1,
                      address=address, volume_pc3=spherical.cell_volume(shell, band, wedge, skeleton.edge_pc),
                      galaxy_seed=skeleton.galaxy_seed)

