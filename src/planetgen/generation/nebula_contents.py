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
* Supernova remnants keep their neutron star or black hole through `SupernovaRemnant` itself.
"""

import math
import re

from planetgen import tuning
from planetgen.galaxy import seed as galaxySeed
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula, choose_weighted_class
from planetgen.physics import stellar_evolution
from planetgen.physics.units import ly_to_pc
from planetgen.util import draw

FALLBACK_SEED = bytes(16)
"""bytes: The seed of a galaxy planned before it had one, as `nebula_field` takes it."""

POSITION_ROUND_MPC = 1
"""int: A cloud's centre is rounded to this many milliparsecs before it names the cloud's seed."""


def hosted_seed(star):
    """The seed of the stream a star's hosted-nebula roll and nebula are drawn from: the star's own mass, age and
    temperature, so a pre-placed star and the same star made in a fill agree."""
    return f"hosted-nebula:{star.mass!r}:{star.age!r}:{star.temperature!r}"


def hosted_nebula_of(star, chance_factor=1.0, markdown=False):
    """
    The nebula a main-sequence star lights (`tuning.NEBULA_HOST_RULES`), or None: the star rolls its row's chance
    (scaled by the sector profile's `"nebula"` factor `chance_factor`) on its own stream, and on a hit the class is
    drawn among those the star is young enough for and the nebula built from the same stream.

    Args:
        star: A `Star` (or anything with `type`, `age` in Gy, `mass` and `temperature`).
        chance_factor (float): The sector profile's `"nebula"` factor (GEN.236).
        markdown (bool): The nebula's `SystemConfig.MARKDOWN`.

    Returns:
        Nebula or None.
    """
    match = re.fullmatch(r"([OBAFGKM])([0-9])V", (getattr(star, "type", "") or "").split(" ", 1)[0])
    if match is None:
        return None
    letter, subclass = match.group(1), int(match.group(2))
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
            return Nebula(config, nebula_class=nebula_class, host=(letter, subclass), age_myr=age_myr)
    return None


def planetary_central_config(config):
    """Gives the `SystemConfig` of a planetary nebula's central white dwarf its companion roll: a close companion in
    `PLANETARY_NEBULA_COMPANION_CHANCE` of nebulae (drawn from the stream now bound), none in the rest."""
    config.BINARY_SYSTEM = draw.random() < tuning.PLANETARY_NEBULA_COMPANION_CHANCE
    config.WIDE_BINARY = False
    return config


def cloud_seed(galaxy_seed, center_pc):
    """The seed of one field cloud's contents: the galaxy's, and the cloud's centre in whole milliparsecs."""
    key = tuple(int(round(c * 1000.0 / POSITION_ROUND_MPC)) for c in center_pc)
    return galaxySeed.short_seed(galaxy_seed or FALLBACK_SEED, "nebula-members", key)


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
