# planetgen/generation/pulsar_planets.py

"""
Planets of a millisecond pulsar (GEN.317; docs/design/systems/exotic-environments-planets-and-compact-binaries.md section 2.1)
=====================================================================================================================

Whether a pulsar holds planets at all is the remnant satellite table's (`phenomena.hosts`: 0.7% of millisecond pulsars,
none for a young pulsar or any other remnant). The ones that do have one of three origins, drawn here from the four
known systems (`tuning.PULSAR_PLANET_ORIGIN_WEIGHTS`, low confidence):

* **disk** (25%): the carbon-rich fallback disk of a disrupted white dwarf, as PSR B1257+12's. One to three planets
  between 0.1 and 0.6 AU, masses log-uniform from 0.01 to 10 Earth masses, each period near a low-order ratio to the
  one inside it. The mass range is drawn whole (GEN.593: the heaviest known, PSR B1257+12's 4.3 Earth masses, is inside it; the radius follows from the mass by the
  terrestrial law of `bulk_draw`, so no radius ceiling cuts the top). They are carbon-rich (the disk's C/O, `bulk_draw.REMNANT_DISK_RATIOS`) and have no moons. The pulsar
  wind (`magnetism.pulsar_wind_pressure_pa`) crushes their magnetosphere to the surface, and the radiation rating
  (`habitability_world.compact_host_name`) always puts them in the lowest tier.
* **captured-giant** (25%): a giant taken in a globular cluster (PSR B1620-26 b), one planet of 1 to 4 Jupiter masses at 20 to
  50 AU (`tuning.PULSAR_GIANT_MASS_JUPITER_RANGE`, `PULSAR_GIANT_DISTANCE_AU_RANGE`), built through the ordinary planet classes (a giant that formed round another star, so it is
  not a disk body and keeps what the planet classes give a giant of that mass at that distance). Built only for a pulsar of a cluster (GEN.163: `star.cluster_id`, the cluster the scatter drew it in,
  `phenomenon_scatter.cluster_id`): a field pulsar that draws it stays bare, as the design rule says ("restricted to millisecond pulsars whose position lies inside a cluster").
* **ablated-companion** (50%): the stripped core of a companion of 0.001 to 0.05 solar masses in an orbit of hours to
  a day (PSR J1719-1438 b, J2322-2650 b). Not built: the planet classes hold no sub-stellar carbon body, and the
  companion is a pair member the pair builder makes (GEN.555, `tuning.PULSAR_PLANET_PREREQUISITES`).

The stand-in of an origin that is not built (GEN.593): the pulsar stays a bare pulsar (the system pass counts it `:empty`), so a planet-holding pulsar holds planets
only for the origins built and the rest of the 0.7% show none, the honest shortfall until the prerequisite is built; the log line of the draw names the origin and what it waits for.
`tuning.PULSAR_PLANET_BUILT_ORIGINS` is the one list to extend as each is built; the captured giant (GEN.163) is built, and the ablated companion (GEN.555) is the one left.
"""

import math

from planetgen import tuning
from planetgen.generation.planet import Planet
from planetgen.util import draw
from planetgen.util import log
from planetgen.physics import unit_registry as units

_mearth_to_kg = units.converter("Mearth", "kg")
_mjupiter_to_kg = units.converter("Mjup", "kg")


def applies(star):
    """Whether `star` is a millisecond pulsar: the only star whose planets this module draws."""
    return getattr(star, "pulsar_type", None) == "millisecond"


def draw_origin(system_config):
    """The origin of a pulsar's planets (`tuning.PULSAR_PLANET_ORIGIN_WEIGHTS`); a run that asks for planets
    (`PLANETS` true) gets the disk, the one origin that is built."""
    if system_config.PLANETS is True:
        return "disk"
    origins = list(tuning.PULSAR_PLANET_ORIGIN_WEIGHTS)
    weights = [tuning.PULSAR_PLANET_ORIGIN_WEIGHTS[name] for name in origins]
    origin = draw.choices(origins, weights=weights)[0]
    log.choice("Pulsar planet origin", origin, f"weights {dict(zip(origins, weights))}")
    return origin


def _disk_distances_au(count):
    """`count` distances, AU, in order: the inner planet's log-uniform in a span that keeps the outermost inside the
    disk's outer edge, each next one at the period ratio drawn for it (a low-order ratio scattered by
    `PULSAR_DISK_PERIOD_RATIO_SCATTER`) with Kepler's third law (the planets weigh nothing against the pulsar)."""
    ratios = [draw.choice(tuning.PULSAR_DISK_PERIOD_RATIOS) * (1.0 + draw.gauss(0.0, tuning.PULSAR_DISK_PERIOD_RATIO_SCATTER))
              for _ in range(count - 1)]
    steps = [ratio ** (2.0 / 3.0) for ratio in ratios]
    span = math.prod(steps)
    inner, outer = tuning.PULSAR_DISK_DISTANCE_AU_RANGE
    high = max(inner, outer / span)
    first = math.exp(draw.uniform(math.log(inner), math.log(high)))
    distances = [first]
    for step in steps:
        distances.append(distances[-1] * step)
    return distances


def _captured_giant(system_config, star, habitable_zone, orbit_ceiling_au):
    """
    The captured giant of a pulsar of a globular cluster (PSR B1620-26 b): one planet of `PULSAR_GIANT_MASS_JUPITER_RANGE` Jupiter masses at `PULSAR_GIANT_DISTANCE_AU_RANGE`, both log-uniform. A pulsar
    of the field (`cluster_id` `None`) gets none, and nor does one whose draw lands past `orbit_ceiling_au`.
    """
    if getattr(star, "cluster_id", None) is None:
        log.choice("Pulsar planets", "none", "origin captured-giant needs a pulsar of a globular cluster and this one is of the field; the pulsar stays bare")
        return []
    low, high = tuning.PULSAR_GIANT_MASS_JUPITER_RANGE
    mass_jupiter = math.exp(draw.uniform(math.log(low), math.log(high)))
    near, far = tuning.PULSAR_GIANT_DISTANCE_AU_RANGE
    distance = math.exp(draw.uniform(math.log(near), math.log(far)))
    if distance > orbit_ceiling_au:
        log.choice("Pulsar planets", "none", f"origin captured-giant drew {distance:.3g} AU, past the orbit ceiling of {orbit_ceiling_au:.3g} AU")
        return []
    planet = Planet(system_config, star, habitable_zone, distance, mass=_mjupiter_to_kg(mass_jupiter))
    log.choice("Pulsar planets", "1 captured giant", f"origin captured-giant, {mass_jupiter:.3g} Jupiter masses at {distance:.3g} AU")
    return [planet]


def generate(system_config, star, habitable_zone, orbit_ceiling_au):
    """
    The planets of the millisecond pulsar `star`, in orbital order: none for an origin that is not built, or when
    `PLANETS` is false. `NUM_ORBITS` sets the disk's count when given.

    Args:
        system_config (SystemConfig): The system's config.
        star (NeutronStar): The pulsar.
        habitable_zone (tuple): The pulsar's habitable zone, AU.
        orbit_ceiling_au (float): Nothing is placed past it.

    Returns:
        list: `Planet` instances.
    """
    if system_config.PLANETS is False:
        return []
    origin = draw_origin(system_config)
    if origin not in tuning.PULSAR_PLANET_BUILT_ORIGINS:
        log.choice("Pulsar planets", "none", f"origin {origin} is not built yet (waits for {tuning.PULSAR_PLANET_PREREQUISITES[origin]}); the pulsar stays bare")
        return []
    if origin == "captured-giant":
        return _captured_giant(system_config, star, habitable_zone, orbit_ceiling_au)
    count = system_config.NUM_ORBITS
    if count is None:
        count = draw.randint(*tuning.PULSAR_DISK_PLANET_COUNT)
    low, high = tuning.PULSAR_DISK_MASS_EARTH_RANGE
    planets = []
    for distance in _disk_distances_au(count) if count > 0 else []:
        if distance > orbit_ceiling_au:
            break
        mass_earth = math.exp(draw.uniform(math.log(low), math.log(high)))
        planets.append(Planet(system_config, star, habitable_zone, distance,
                              mass=_mearth_to_kg(mass_earth), moon_count=0,
                              formation_distance_au=0.0, late_core=True))      # the disk is the debris of a disrupted star: no gas disk, so no hydrogen envelope (GEN.1323 gave it one)
    log.choice("Pulsar planets", f"{len(planets)} disk planets",
               f"origin disk, {count} drawn, masses {low:.3g} to {high:.3g} Earth masses")
    return planets
