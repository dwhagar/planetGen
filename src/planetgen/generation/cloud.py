# planetgen/generation/cloud.py

"""
The cloud a system sits in while it is generated (GEN.357)
==========================================================

Until now a system learned its nebula or supernova remnant only after it was saved (`store.refresh_containment` sets
`inside_nebula_id`/`inside_remnant_id`, and `store.load_star_system` hangs the cloud on the loaded system). Generation
knows it sooner: the galaxy's field clouds are drawn from the seed alone (`galaxy.nebula_field`), so a sector can list the
ones that reach it before it places a system, and `SpaceSector.add_system` hands each system the cloud holding its position
(`apply_in_sector`, then `StarSystem.apply_cloud`).

The cloud is the same dict the store gives a loaded system (`cloud_dict`, built by `store.surrounding_cloud` from a row and
by `cloud_of_nebula` from a drawn `Nebula`), so a system's cloud has one shape wherever it came from:

====================  =====================================================================================
`type`                `"nebula"` or `"supernova_remnant"`
`id`, `name`          the stored row's; `None` for a field cloud no sector has saved yet (its name is drawn when it is
                      saved, and drawing it here would move the draws of everything after)
`class`, `descriptor` the class letter, and the nebula's family or the remnant's morphology
`density_cm3`, `temperature_k`, `extinction_av`, `radius_ly`
`age_myr`             the age of the hosts (nebula) or of the remnant, Myr; `None` for a cloud without an age
`center`              the cloud's centre in the galaxy frame, parsecs; `None` for a cloud not placed
`host_stars`          the stars that light it: dicts of `type`, `count`, `luminosity_sol`, `teff_k`, `position_pc`
====================  =====================================================================================

Only the field clouds (classes M, N, P and Q, a pure function of the galaxy seed) are found while a sector is generated.
A nebula or remnant another sector stored is not: reading it would make a sector depend on which neighbour was saved
first (GEN.39 makes a sector's draws independent of the run's order), and the post-hoc containment pass already covers it.
GEN.99's regenerate-in-place calls `StarSystem.apply_cloud` with the cloud the store gives it.

A system outside every cloud draws nothing and is changed in nothing. A system inside one draws from a stream of its own
(`fill_stages.CLOUDS`, one per placement), so neither it nor the systems after it move any other draw.

Out of scope here (GEN.358): the rules that change planets (photoevaporation cut, planetary-nebula discs, the swallowed
giant flag, orbit expansion). `apply_cloud` stamps `StarSystem.cloud_effects` and regenerates nothing yet. The nebula's
composition flags of the design sketch are not made: composition comes through the composition chain (GEN.428, GEN.440),
and `enrichment_k` is a factor for it to read, not a second set of composition fields.
"""

import dataclasses
import math
import re

from planetgen import tuning
from planetgen.galaxy import containment, frame
from planetgen.physics import hii, main_sequence
from planetgen.physics.units import ly_to_pc
from planetgen.util import detmath, draw
from planetgen.generation import fill_stages

LISTED_KEYS = ("type", "id", "name", "class", "descriptor", "density_cm3", "temperature_k")
"""tuple: The keys of the cloud the API and the pages show for the cloud a row is inside (`query.containing_cloud`); the
rest of the dict is for generation."""

_HOST_COUNT = re.compile(r"(\w+) x(\d+)")

SOLAR_LUMINOSITY_W = 3.828e26
"""float: Watts in one solar luminosity (IAU 2015 nominal)."""

_PC_CM = 3.0857e18


@dataclasses.dataclass(frozen=True)
class CloudEffects:
    """
    What a cloud does to a system inside it (`StarSystem.cloud_effects`), derived at generation and not stored.

    Attributes:
        cloud_id (int or None): The cloud's stored id, `None` for a field cloud not saved yet.
        cloud_class (str): The cloud's class letter.
        g0 (float): The far-ultraviolet flux of the cloud's hosts at the system, Habing units (0 without hosts).
        radiation_tier (int): The tier the habitability index adds to the system's dose tier, 0 to 3: the larger of
            the class's own tier (`tuning.CLOUD_CLASS_RADIATION_TIER`) and the G0 tier (`g0_tier`).
        enrichment_k (float or None): 26Al/27Al in units of `tuning.CLOUD_ENRICHMENT_REFERENCE_RATIO`, drawn for
            the classes in `tuning.CLOUD_ENRICHMENT_CLASSES`; `None` for the rest.
    """

    cloud_id: object
    cloud_class: str
    g0: float
    radiation_tier: int
    enrichment_k: object


def cloud_dict(kind, ident, name, cloud_class, descriptor, density_cm3, temperature_k, extinction_av, radius_ly, age_myr,
               center_pc, hosts=None):
    """
    A cloud as the dict a system's `surrounding_cloud` holds (see the module docstring); `hosts` is a nebula's
    `hosts` text (`"O5V x1, O9V x4"`) or `None`.
    """
    center = None if center_pc is None else tuple(center_pc)
    return {"type": kind, "id": ident, "name": name, "class": cloud_class, "descriptor": descriptor,
            "density_cm3": density_cm3, "temperature_k": temperature_k, "extinction_av": extinction_av,
            "radius_ly": radius_ly, "age_myr": age_myr, "center": center,
            "host_stars": host_stars(hosts, center)}


def cloud_of_nebula(nebula, center_pc, ident=None, name=None):
    """The cloud dict of a drawn `Nebula` at `center_pc` (galaxy frame, parsecs); `ident` and `name` are the stored
    row's, when there is one."""
    return cloud_dict("nebula", ident, name, nebula.nebula_class, nebula.nebula_type, nebula.density_cm3,
                      nebula.temperature_k, nebula.extinction_av, nebula.radius_ly, nebula.age_myr, center_pc,
                      getattr(nebula, "hosts", None))


def host_stars(hosts, position_pc):
    """
    The stars that light a nebula, from its `hosts` text, as dicts of `type`, `count`, `luminosity_sol`, `teff_k` and
    `position_pc`. Every host sits at `position_pc` (the cloud's centre): a seeded cloud has no stored positions for
    its stars, so a system's distance to them is its distance to the centre, within the cluster's half-mass radius
    (`tuning.NEBULA_CLUSTER_HALF_MASS_PC`) of the truth. A type off the O3V to B5V ladder (a Wolf-Rayet star) has no
    luminosity here and is left out; so is a cloud without a position.
    """
    if not hosts or position_pc is None:
        return []
    mass_of = {name: mass for name, mass, _log_q in hii.TYPE_ANCHORS}
    stars = []
    for kind, count in _HOST_COUNT.findall(hosts):
        if kind not in mass_of:
            continue
        mass = mass_of[kind]
        arrival = main_sequence.zams_state(mass)
        stars.append({"type": kind, "count": int(count), "luminosity_sol": arrival.luminosity_sol,
                      "teff_k": arrival.temperature_k,
                      "position_pc": tuple(position_pc)})
    return stars


def fuv_fraction(teff_k):
    """The share of a star's light in the far ultraviolet at `teff_k` (`tuning.CLOUD_FUV_FRACTION_BY_TEFF_K`)."""
    rows = sorted(tuning.CLOUD_FUV_FRACTION_BY_TEFF_K.items())
    if teff_k <= rows[0][0]:
        return rows[0][1]
    for (low_t, low_f), (high_t, high_f) in zip(rows, rows[1:]):
        if teff_k <= high_t:
            return low_f + (teff_k - low_t) / (high_t - low_t) * (high_f - low_f)
    return rows[-1][1]


def g0_at(cloud, point_pc):
    """The far-ultraviolet flux, Habing units, of a cloud's host stars at `point_pc` (galaxy frame): the sum of
    `f_FUV L / (4 pi d^2)` over every host, each no nearer than `tuning.CLOUD_HOST_MIN_DISTANCE_PC`. Unattenuated, so an
    upper bound (the design note ignores the gas's extinction in this first version)."""
    total = 0.0
    for star in cloud["host_stars"]:
        distance_pc = max(detmath.dist(point_pc, star["position_pc"]), tuning.CLOUD_HOST_MIN_DISTANCE_PC)
        flux = (fuv_fraction(star["teff_k"]) * star["luminosity_sol"] * SOLAR_LUMINOSITY_W * 1.0e7
                / (4.0 * math.pi * (distance_pc * _PC_CM) ** 2))
        total += star["count"] * flux / tuning.HABING_FLUX_ERG_CM2_S
    return total


def g0_tier(g0):
    """How many of `tuning.CLOUD_G0_TIER_EDGES` the flux `g0` reaches, 0 to 3."""
    return sum(1 for edge in tuning.CLOUD_G0_TIER_EDGES if g0 >= edge)


def draw_enrichment(cloud_class, rng):
    """The 26Al enrichment `k` of a system in a cloud of the class (design note, "Rule table": log-normal about the field
    median, or about the rich median with the class's chance), or `None` for a class that does not enrich."""
    if cloud_class not in tuning.CLOUD_ENRICHMENT_CLASSES:
        return None
    rich = rng.random() < tuning.CLOUD_ENRICHMENT_RICH_CHANCE.get(cloud_class, 0.0)
    median = tuning.CLOUD_ENRICHMENT_RICH_MEDIAN if rich else tuning.CLOUD_ENRICHMENT_FIELD_MEDIAN
    return median * 10.0 ** rng.gauss(0.0, tuning.CLOUD_ENRICHMENT_SIGMA_DEX)


def effects_of(cloud, point_pc, rng):
    """The `CloudEffects` of `cloud` on a system at `point_pc` (galaxy frame, parsecs), drawing from `rng`."""
    g0 = g0_at(cloud, point_pc)
    tier = max(tuning.CLOUD_CLASS_RADIATION_TIER.get(cloud["class"], 0), g0_tier(g0))
    return CloudEffects(cloud_id=cloud["id"], cloud_class=cloud["class"], g0=g0, radiation_tier=tier,
                        enrichment_k=draw_enrichment(cloud["class"], rng))


def field_containers(field_nebulae):
    """The containers (`galaxy.containment`) of the `(Nebula, center_pc)` pairs of a sector's field clouds, each with its
    cloud dict under `cloud`. They are not stored, so their `id` is `None`."""
    return [{"column": "inside_nebula_id", "id": None, "center": tuple(center), "radius_pc": ly_to_pc(nebula.radius_ly),
             "shape": nebula.get_shape, "cloud": cloud_of_nebula(nebula, center)}
            for nebula, center in field_nebulae]


def begin(sector, cloud_field, field_nebulae):
    """
    Keeps on `sector` the field clouds that reach it (`nebula_field.clouds_reaching(*cloud_field)`, from the galaxy seed
    alone), found before it places a system: `field_nebulae` (stored later by the save), `cloud_field` and
    `cloud_containers`, which `apply_in_sector` tests every placement against.
    """
    sector.cloud_field = cloud_field
    sector.field_nebulae = field_nebulae
    sector.cloud_containers = field_containers(field_nebulae)


def apply_in_sector(sector, star_system, position_ly):
    """
    Gives `star_system`, just placed at the sector-local `position_ly`, the cloud holding that point (the innermost,
    `containment.innermost_container`), if any, through `StarSystem.apply_cloud`. Draws from the clouds stream of this
    placement (`fill_stages.CLOUDS`), never from the one bound now. Returns the cloud dict or `None`.
    """
    center_pc = sector.cloud_field[2]
    point = frame.local_to_galaxy_pc(center_pc, tuple(ly_to_pc(c) for c in position_ly))
    found = containment.innermost_container(point, sector.cloud_containers)
    if found is None:
        return None
    seed = fill_stages.stream_seed(sector, fill_stages.CLOUDS, len(sector.entries))
    star_system.apply_cloud(found["cloud"], point, rng=draw if seed is None else draw.Stream(seed))
    fill_stages.count(sector, fill_stages.CLOUDS, "star-system")
    return found["cloud"]
