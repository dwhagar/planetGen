# planetgen/generation/nebula_contents.py

"""
The stars a nebula needs, made with it (GEN.320, GEN.241)
=========================================================

A nebula of a class that has stars gets them as real systems, inside its bounds, from seeds that depend on the nebula
or the star alone, never on the order a sector is filled in:

* **Hosted nebulae (C, D, E, F, G).** Whether a hot star has a nebula, and which class, is rolled on a stream of the
  star's own (`hosted_seed`), so it is the same for a star the scatter pre-placed and for one made during a fill, and a
  scatter row can say it before any fill (`hosted_nebula_of`).
* **Planetary nebulae.** The central white dwarf is made with the nebula; `PLANETARY_NEBULA_COMPANION_CHANCE` of them
  gain a close companion (`planetary_central_config`).
* **Bok globules (P) and star-forming cores (Q).** Their protostars come from the cloud's own seed (`cloud_members`),
  positions in the galaxy frame, so a cloud across several sectors needs no coordination: each sector builds only the
  members whose points fall inside it (`run_sector.add_cloud_protostars`).
* **Emission regions D and E (GEN.321).** The cluster behind the seed star is its ionizing members: the stars of 5.5
  solar masses and up that are still alive, which the nebula's `hosts` text counts by spectral type. `emission_members`
  makes them real from the seed star's own point (the nebula's seed): a mass in each type's range from the IMF, the
  nebula's age, and a place in a Plummer sphere of half-mass radius 1 pc (D) or 2 pc (E), the profile and sampler of
  the satellites (`physics.component_profiles`), clipped to the nebula's radius. The members are a function of the
  galaxy, the nebula and the host's point alone, so every sector that reaches the region gives the same ones and builds
  those whose points fall in its cube (`run_sector.add_emission_members`); a sector finds the regions of its neighbours'
  pre-placed hot stars from their rows (`row_host`, `emission_reach_pc`), whether or not that neighbour was filled.
  Only the ionizing members are real: the rest of the cluster stays the sector's ordinary background density.
* Supernova remnants keep their neutron star or black hole through `SupernovaRemnant` itself.
"""

import math
import re
import types

import numpy as np

from planetgen import tuning
from planetgen.galaxy import seed as galaxySeed
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula, choose_weighted_class
from planetgen.physics import component_profiles, hii, stellar_evolution
from planetgen.physics.units import ly_to_pc
from planetgen.util import draw

FALLBACK_SEED = bytes(16)
"""bytes: The seed of a galaxy planned before it had one, as `nebula_field` takes it."""

POSITION_ROUND_MPC = 1
"""int: A cloud's centre is rounded to this many milliparsecs before it names the cloud's seed."""


def hosted_seed(star):
    """The seed of the stream a star's hosted-nebula roll and nebula are drawn from: the star's own mass, age and
    temperature, so a pre-placed star and the same star made in a fill agree. The numbers go in as floats: a star
    made in a fill holds a whole-kelvin temperature as an integer, the same star read back from its `bright_stars` row
    (`row_host`) as a float, and `repr` tells `36900` from `36900.0`."""
    return f"hosted-nebula:{float(star.mass)!r}:{float(star.age)!r}:{float(star.temperature)!r}"


def hosted_nebula_of(star, chance_factor=1.0, markdown=False):
    """
    The nebula a main-sequence or Wolf-Rayet star lights (`tuning.NEBULA_HOST_RULES`), or None: the star rolls its row's chance
    (scaled by the sector profile's `"nebula"` factor `chance_factor`) on its own stream, and on a hit the class is
    drawn among those the star is young enough for and the nebula built from the same stream.

    Args:
        star: A `Star` (or anything with `type`, `age` in Gy, `mass` and `temperature`).
        chance_factor (float): The sector profile's `"nebula"` factor (GEN.236).
        markdown (bool): The nebula's `SystemConfig.MARKDOWN`.

    Returns:
        Nebula or None.
    """
    match = re.fullmatch(r"(?:([OBAFGKM])([0-9])V|(W[NCO])([0-9]))", (getattr(star, "type", "") or "").split(" ", 1)[0])
    if match is None:
        return None
    if match.group(3):
        # GEN.243: a Wolf-Rayet star (WN6) is matched on the letter W, and the nebula keeps its subtype.
        letter, subclass, host = "W", int(match.group(4)), (match.group(3), int(match.group(4)))
    else:
        letter, subclass = match.group(1), int(match.group(2))
        host = (letter, subclass)
    age_myr = star.age * 1000.0
    for rule_letter, low, high, classes, chance in tuning.NEBULA_HOST_RULES:
        if letter != rule_letter or not (low <= subclass <= high):
            continue
        young = tuple(c for c in classes if age_myr <= tuning.NEBULA_CLASS_MAX_AGE_MYR[c])
        if not young:
            return None
        with draw.bound(hosted_seed(star)):
            if draw.random() >= chance * chance_factor:
                return None
            config = SystemConfig()
            config.MARKDOWN = markdown
            nebula_class = choose_weighted_class(young, "Star-hosted nebula class")
            return Nebula(config, nebula_class=nebula_class, host=host, age_myr=age_myr)
    return None


def planetary_central_config(config):
    """Gives the `SystemConfig` of a planetary nebula's central white dwarf its companion roll: a close companion in
    `PLANETARY_NEBULA_COMPANION_CHANCE` of nebulae (drawn from the stream now bound), none in the rest."""
    config.BINARY_SYSTEM = draw.random() < tuning.PLANETARY_NEBULA_COMPANION_CHANCE
    config.WIDE_BINARY = False
    return config


def _point_seed(galaxy_seed, kind, point_pc):
    """A seed from the galaxy's, a name and a point in whole milliparsecs."""
    key = tuple(int(round(c * 1000.0 / POSITION_ROUND_MPC)) for c in point_pc)
    return galaxySeed.short_seed(galaxy_seed or FALLBACK_SEED, kind, key)


def cloud_seed(galaxy_seed, center_pc):
    """The seed of one field cloud's contents: the galaxy's, and the cloud's centre in whole milliparsecs."""
    return _point_seed(galaxy_seed, "nebula-members", center_pc)


def emission_seed(galaxy_seed, host_position_pc):
    """The seed of an emission region's members: the galaxy's, and the seed star's point in whole milliparsecs."""
    return _point_seed(galaxy_seed, "emission-members", host_position_pc)


class CloudMember:
    """One protostar of a field cloud: its galaxy-frame point `position_pc`, its initial mass `mass_sol`, its age
    `age_yr` and the `seed` its system is built from."""

    __slots__ = ("position_pc", "mass_sol", "age_yr", "seed")

    def __init__(self, position_pc, mass_sol, age_yr, seed):
        self.position_pc = position_pc
        self.mass_sol = mass_sol
        self.age_yr = age_yr
        self.seed = seed


def _protostar_count(rule, rng):
    if "chance" in rule:
        return 1 if rng.random() < rule["chance"] else 0
    limit = math.exp(-rule["mean"])
    while True:  # zero-truncated Poisson: a core that forms stars has at least one
        count, product = 0, rng.random()
        while product > limit:
            count += 1
            product *= rng.random()
        if count >= 1:
            return min(count, tuning.NEBULA_CLOUD_PROTOSTAR_MAX)


def cloud_members(galaxy_seed, nebula, center_pc):
    """
    The protostars of a field cloud, from its seed alone: none for a class without any, one in `chance` of Bok globules
    (P) at the centre, a clump of one to four for a star-forming core (Q), Gaussian about the centre and kept inside
    the cloud's radius. A member's mass is drawn from the IMF below `PRE_MAIN_SEQUENCE_MAX_MASS_SOL` and its age
    uniformly under `PROTOSTAR_MAX_AGE_YR` (Class 0 and I).

    Returns:
        list: `CloudMember`s, the same every time for the same galaxy and cloud.
    """
    rule = tuning.NEBULA_CLOUD_PROTOSTARS.get(nebula.nebula_class)
    if rule is None:
        return []
    base = cloud_seed(galaxy_seed, center_pc)
    rng = draw.Stream(f"cloud-members:{base}")
    radius_pc = ly_to_pc(nebula.radius_ly)
    sigma = radius_pc * tuning.NEBULA_CLOUD_CLUMP_FRACTION
    members = []
    for n in range(_protostar_count(rule, rng)):
        while True:
            offset = tuple(rng.gauss(0.0, sigma) for _ in range(3)) if "mean" in rule else (0.0, 0.0, 0.0)
            if math.sqrt(sum(c * c for c in offset)) <= radius_pc:
                break
        mass = stellar_evolution.sample_imf_mass_sol(None, tuning.PRE_MAIN_SEQUENCE_MAX_MASS_SOL, rng)
        age = rng.uniform(0.0, tuning.PROTOSTAR_MAX_AGE_YR)
        members.append(CloudMember(tuple(center_pc[i] + offset[i] for i in range(3)), mass, age, f"{base}:{n}"))
    return members


EMISSION_CLASSES = tuple(tuning.NEBULA_CLUSTER_HALF_MASS_PC)
"""tuple: The classes whose ionizing members are real systems (D and E)."""

_HOST_COUNT = re.compile(r"(\w+) x(\d+)")
_DEAD_MEMBER_TRIES = 200


def _plummer_table():
    return component_profiles.profile_table("plummer")


def emission_reach_pc():
    """How far from its seed star a member of any emission region can lie, parsecs: the Plummer table's cut radius for
    the widest class."""
    table = _plummer_table()
    return max(tuning.NEBULA_CLUSTER_HALF_MASS_PC.values()) / table.half_mass_x * table.tidal_x


def host_name_of(star):
    """The spectral code of a star (`"O5V"`), the first word of its type, as the nebula's `hosts` text names the seed."""
    return (getattr(star, "type", "") or "").split(" ", 1)[0]


def member_counts(nebula, host_name):
    """
    The ionizing members of an emission region by spectral type, hottest first, `[(type, count), ...]`: the `hosts` text
    less the seed star itself (one of `host_name`).
    """
    if nebula.nebula_class not in EMISSION_CLASSES or not nebula.hosts:
        return []
    counts = [[name, int(count)] for name, count in _HOST_COUNT.findall(nebula.hosts)]
    for entry in counts:
        if entry[0] == host_name:
            entry[1] -= 1
            break
    return [(name, count) for name, count in counts if count > 0]


def _member_mass(name, age_myr, rng):
    """A mass for a member of type `name`: the IMF within the type's range, redrawn while the star would be dead by
    `age_myr`. The `hosts` text counted only survivors, so a type's range always holds one that lives."""
    low, high = hii.mass_range_of_type(name, tuning.NEBULA_SEED_MASS_RANGE_SOL["E"][1])
    age_gy = age_myr / 1000.0
    for _ in range(_DEAD_MEMBER_TRIES):
        mass = stellar_evolution.sample_imf_mass_sol(low, high, rng)
        if stellar_evolution.main_sequence_lifetime_gy(mass) > age_gy:
            return mass
    return low


def _plummer_offset(radius_pc, half_mass_pc, rng):
    """An offset from the Plummer sphere's centre, parsecs: a radius by inverse transform of the profile's mass table,
    redrawn while it lies outside the nebula, a direction uniform on the sphere."""
    table = _plummer_table()
    scale = half_mass_pc / table.half_mass_x
    for _ in range(_DEAD_MEMBER_TRIES):
        r = scale * float(np.interp(rng.random(), table.mass_fraction, table.x))
        if r <= radius_pc:
            break
    else:
        r = radius_pc * rng.random() ** (1.0 / 3.0)
    cos_colat = 2.0 * rng.random() - 1.0
    in_plane = r * math.sqrt(max(0.0, 1.0 - cos_colat * cos_colat))
    lon = 2.0 * math.pi * rng.random()
    return (in_plane * math.cos(lon), in_plane * math.sin(lon), r * cos_colat)


def emission_members(galaxy_seed, nebula, host_name, host_position_pc):
    """
    The ionizing members of an emission region (class D or E) other than its seed star, from the galaxy, the nebula and
    the seed star's point alone: one per count of the nebula's `hosts` text less the seed, each a mass in its type's
    range, the nebula's age, and a point in a Plummer sphere about the seed star (half-mass radius
    `tuning.NEBULA_CLUSTER_HALF_MASS_PC`) inside the nebula's radius.

    Args:
        galaxy_seed (bytes): The galaxy's seed.
        nebula (Nebula): The region.
        host_name (str): The seed star's spectral code.
        host_position_pc (tuple): The seed star's galaxy-frame point, parsecs.

    Returns:
        list: `CloudMember`s, the same every time for the same galaxy, nebula and point. Empty for another class.
    """
    counts = member_counts(nebula, host_name)
    if not counts:
        return []
    base = emission_seed(galaxy_seed, host_position_pc)
    rng = draw.Stream(f"emission-members:{base}")
    radius_pc = ly_to_pc(nebula.radius_ly)
    half_mass_pc = tuning.NEBULA_CLUSTER_HALF_MASS_PC[nebula.nebula_class]
    members = []
    for name, count in counts:
        for _ in range(count):
            mass = _member_mass(name, nebula.age_myr, rng)
            offset = _plummer_offset(radius_pc, half_mass_pc, rng)
            members.append(CloudMember(tuple(host_position_pc[i] + offset[i] for i in range(3)), mass,
                                       nebula.age_myr * 1.0e6, f"{base}:{len(members)}"))
    return members


def row_host(row, chance_factor=1.0):
    """
    The emission region a pre-placed hot star would light, from its `bright_stars` row alone: `(nebula, host_name,
    host_position_pc)`, or None when it lights none or one of another class. The star's roll is on its own stream
    (`hosted_nebula_of`), so this is the nebula its sector builds, and the row's point is its place.
    """
    star = types.SimpleNamespace(type=row["star_type"], age=row["age_gy"], mass=row["mass_kg"],
                                 temperature=row["temperature_k"])
    nebula = hosted_nebula_of(star, chance_factor)
    if nebula is None or nebula.nebula_class not in EMISSION_CLASSES:
        return None
    point = tuple(row[f"position_{axis}_mpc"] / 1000.0 for axis in "xyz")
    return nebula, host_name_of(star), point
