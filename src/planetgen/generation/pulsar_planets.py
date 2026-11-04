# planetgen/generation/pulsar_planets.py

"""
Planets of a millisecond pulsar (GEN.317; docs/design/systems/exotic-environments-planets-and-compact-binaries.md section 2.1)
=====================================================================================================================

Whether a pulsar holds planets at all is the remnant satellite table's (`phenomena.hosts`: 0.7% of millisecond pulsars,
none for a young pulsar or any other remnant). The ones that do have one of three origins, drawn here from the four
known systems (`tuning.PULSAR_PLANET_ORIGIN_WEIGHTS`, low confidence):

* **disc** (25%): the carbon-rich fallback disc of a disrupted white dwarf, as PSR B1257+12's. One to three planets
  between 0.1 and 0.6 AU, masses log-uniform from 0.01 to 10 Earth masses, each period near a low-order ratio to the
  one inside it. They are carbon-rich (the disc's C/O, `bulk_draw.REMNANT_DISC_RATIOS`) and have no moons. The pulsar
  wind (`magnetism.pulsar_wind_pressure_pa`) crushes their magnetosphere to the surface, and the radiation rating
  (`habitability_world.compact_host_name`) always puts them in the lowest tier.
* **captured-giant** (25%): a giant taken in a globular cluster (PSR B1620-26 b), 1 to 4 Jupiter masses at 20 to
  50 AU. Not built: it needs the cluster giant switch (GEN.163), which stays off, so this origin gives the pulsar no
  planets.
* **ablated-companion** (50%): the stripped core of a companion of 0.001 to 0.05 solar masses in an orbit of hours to
  a day (PSR J1719-1438 b, J2322-2650 b). Not built: the planet classes hold no sub-stellar carbon body, and the
  companion is a pair member the pair builder makes (GEN.555), so this origin gives the pulsar no planets until then.

A pulsar whose origin is not built stays a bare pulsar (the system pass counts it `:empty`).
"""

import math

from planetgen import tuning
from planetgen.generation.planet import Planet
from planetgen.physics import constants, relations
from planetgen.util import draw
from planetgen.util import log
from planetgen.physics import unit_registry as units

_kg_to_mearth = units.converter("kg", "Mearth")
_mearth_to_kg = units.converter("Mearth", "kg")

_km_to_m = units.converter("km", "m")
_m_to_km = units.converter("m", "km")


def applies(star):
    """Whether `star` is a millisecond pulsar: the only star whose planets this module draws."""
    return getattr(star, "pulsar_type", None) == "millisecond"


def draw_origin(system_config):
    """The origin of a pulsar's planets (`tuning.PULSAR_PLANET_ORIGIN_WEIGHTS`); a run that asks for planets
    (`PLANETS` true) gets the disc, the one origin that is built."""
    if system_config.PLANETS is True:
        return "disc"
    origins = list(tuning.PULSAR_PLANET_ORIGIN_WEIGHTS)
    weights = [tuning.PULSAR_PLANET_ORIGIN_WEIGHTS[name] for name in origins]
    origin = draw.choices(origins, weights=weights)[0]
    log.choice("Pulsar planet origin", origin, f"weights {dict(zip(origins, weights))}")
    return origin


def disc_mass_range_earth():
    """`(low, high)` Earth masses of a disc planet: the design range, cut at the largest disc body (a sphere of `tuning.PULSAR_DISC_MAX_RADIUS_KM`
    at the mean of its density range)."""
    low, high = tuning.PULSAR_DISC_MASS_EARTH_RANGE
    low_density, high_density = constants.PLANET_DENSITY["t"]
    max_radius_km = tuning.PULSAR_DISC_MAX_RADIUS_KM
    ceiling = _mass_earth(max_radius_km, (low_density + high_density) / 2.0)
    return low, max(low, min(high, ceiling))


def _mass_earth(radius_km, density_g_cm3):
    """The mass, Earth masses, of a sphere of `radius_km` and `density_g_cm3`."""
    volume_m3 = relations.sphere_volume(_km_to_m(radius_km))
    return _kg_to_mearth(volume_m3 * density_g_cm3 * 1000.0)


def _radius_km(mass_earth, density_g_cm3):
    """The radius, km, of a sphere of `mass_earth` and `density_g_cm3`."""
    volume_m3 = _mearth_to_kg(mass_earth) / (density_g_cm3 * 1000.0)
    return _m_to_km((3.0 * volume_m3 / (4.0 * math.pi)) ** (1.0 / 3.0))


def _disc_distances_au(count):
    """`count` distances, AU, in order: the inner planet's log-uniform in a span that keeps the outermost inside the
    disc's outer edge, each next one at the period ratio drawn for it (a low-order ratio scattered by
    `PULSAR_DISC_PERIOD_RATIO_SCATTER`) with Kepler's third law (the planets weigh nothing against the pulsar)."""
    ratios = [draw.choice(tuning.PULSAR_DISC_PERIOD_RATIOS) * (1.0 + draw.gauss(0.0, tuning.PULSAR_DISC_PERIOD_RATIO_SCATTER))
              for _ in range(count - 1)]
    steps = [ratio ** (2.0 / 3.0) for ratio in ratios]
    span = math.prod(steps)
    inner, outer = tuning.PULSAR_DISC_DISTANCE_AU_RANGE
    high = max(inner, outer / span)
    first = math.exp(draw.uniform(math.log(inner), math.log(high)))
    distances = [first]
    for step in steps:
        distances.append(distances[-1] * step)
    return distances


def generate(system_config, star, habitable_zone, orbit_ceiling_au):
    """
    The planets of the millisecond pulsar `star`, in orbital order: none for an origin that is not built, or when
    `PLANETS` is false. `NUM_ORBITS` sets the disc's count when given.

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
        return []
    count = system_config.NUM_ORBITS
    if count is None:
        count = draw.randint(*tuning.PULSAR_DISC_PLANET_COUNT)
    low, high = disc_mass_range_earth()
    low_density, high_density = constants.PLANET_DENSITY["t"]
    max_radius_km = tuning.PULSAR_DISC_MAX_RADIUS_KM
    planets = []
    for distance in _disc_distances_au(count) if count > 0 else []:
        if distance > orbit_ceiling_au:
            break
        mass_earth = math.exp(draw.uniform(math.log(low), math.log(high)))
        radius_km = min(max_radius_km, _radius_km(mass_earth, (low_density + high_density) / 2.0))
        planets.append(Planet(system_config, star, habitable_zone, distance, radius=radius_km,
                              mass=_mearth_to_kg(mass_earth), moon_count=0,
                              formation_distance_au=0.0))
    log.choice("Pulsar planets", f"{len(planets)} disc planets",
               f"origin disc, {count} drawn, masses {low:.3g} to {high:.3g} Earth masses")
    return planets
