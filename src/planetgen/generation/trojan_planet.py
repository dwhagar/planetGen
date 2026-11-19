# planetgen/generation/trojan_planet.py

"""
Whole planets at L4 and L5 of a giant (GEN.1151 slice A, docs/design/systems/lagrange-trojans.md)
=================================================================================================

A giant that grew where it orbits may hold a body of planet mass 60 degrees ahead of it (L4) or behind it (L5), grown from the solids of its 1:1 zone while the gas disk lasted. `build` asks the formation
gate of `physics/trojan_bodies.py` for each giant: the solids the disk held in the two tadpole zones (the nebula's solid surface density where the giant formed, the architecture's disk factor and the star's metals, the
same disk the planets and the belts were made from), the share of them that grows (TB3) over the share of the growth time the star's gas disk allowed (TB4), held to the room the mass rule leaves. A ceiling under the lightest
planet (`tuning.PLANET_MIN_MASS_EARTH`) gives no body; a higher one draws a mass log-uniform between the two. No rate is set, and none is read from the searches: they only bound it (a planet-mass trojan among 2244 Kepler
candidates: none). The body then has to stay: the giant's neighbours and the companion stars leave a share of such a zone stable (`trojans.surviving_share`, the share GEN.1150 gives the swarm), taken here as the chance it stays. A giant that moved from where it formed
loses the zone and holds none.

The planet is made by the planet generator from a stream of its own with the giant's orbit (`Planet(trojan_of=(giant, point))`: the plane, node, periapsis and eccentricity of the giant, 60 degrees ahead or behind in mean
anomaly), as a belt member promoted to a planet is (`belt_planets`), and put in the star's list after its giant. It is made once the system is final, before the irregular moons, so the spacing pass never sees two planets at one
orbit. It is stored as a row of `planets` with `lagrange_point` and the link to its giant (`lagrange_host_id`).
"""

from planetgen import tuning
from planetgen.generation import belt_birth, remnant_orbits
from planetgen.generation.planet import Planet
from planetgen.physics import bulk_draw, disk_truncation, planet_mass, trojan_bodies, trojans
from planetgen.physics import unit_registry as units
from planetgen.util import draw, log

ORIGIN_SENTENCE = "Grown from the solids of the {point} zone of {host}, which it shares an orbit with."

def _disk_scatter(system, star):
    """The factor the star's planets' disk was scattered by (`planet_architecture.Architecture.disk_scatter`, GEN.1261); 1.0 for a star the architecture did not plan."""
    return next((architecture.disk_scatter for planned, architecture in getattr(system, "architectures", ()) if planned is star), 1.0)


def zone_ceiling(system, star, giant, giants):
    """
    `(ceiling in Earth masses, solids in the zone in Earth masses)` of what the zone of `giant` grew, or `(0.0, 0.0)` for a giant that moved from where it formed. `giants` are the star's giants, nearest first.
    """
    formed_au = giant.formation_distance_au or giant.distance
    if not trojan_bodies.zone_kept(giant.distance, formed_au):
        return 0.0, 0.0
    sigma = planet_mass.solids_surface_density_gcm2(formed_au, remnant_orbits.formation_mass_kg(star), star.snow_line_au, bulk_draw.metal_ratio_to_solar(star), _disk_scatter(system, star))
    star_earth, giant_earth = units.convert(star.mass, "kg", "Mearth"), units.convert(giant.mass, "kg", "Mearth")
    solids = trojan_bodies.zone_solids_earth(sigma, formed_au, giant_earth, star_earth)
    life_myr = disk_truncation.disk_time_yr(star) / 1.0e6
    return trojan_bodies.ceiling_earth(solids, life_myr, star_earth, giant_earth), solids


def survival_of(star, giant, giants):
    """The share of such a zone the giant's neighbours and the companion stars leave stable (`trojans.surviving_share`); the chance the body stays."""
    neighbours = [(other.mass, other.distance) for other in giants if other is not giant]
    return trojans.surviving_share(star.mass, giant.mass, giant.distance, neighbours, getattr(star, "a_crit_au", None))


def build(system, star, bodies):
    """
    The planets at L4 and L5 of the giants among `bodies` (one star's list, final): at most one a giant, each from the giant's own stream (`trojan-planet:` its distance, mass and the star's), nearest the star first.
    Not placed in `bodies`: the caller puts them in the list and settles them (`StarSystem._add_trojan_planets`).
    """
    giants = sorted((body for body in bodies if getattr(body, "body_type", None) == 'g' and getattr(body, "lagrange_point", None) is None), key=lambda body: body.distance)
    made = []
    for giant in giants:
        stream = draw.Stream(f"trojan-planet:{giant.distance!r}:{giant.mass!r}:{star.mass!r}")
        ceiling, solids = zone_ceiling(system, star, giant, giants)
        mass_earth = trojan_bodies.draw_mass_earth(ceiling, tuning.PLANET_MIN_MASS_EARTH, stream.random())
        leading, stays = trojans.leading_share(stream.random()), stream.random()      # the same draws whether or not a body forms, so the stream keeps its place
        if mass_earth is None or stays >= survival_of(star, giant, giants):
            continue
        point = "L4" if stream.random() < leading else "L5"
        with draw.bound(draw.Stream(f"trojan-planet-make:{giant.distance!r}:{giant.mass!r}:{star.mass!r}:{point}")):
            planet = Planet(system.system_config, star, star.habitable_zone, giant.distance, mass=units.convert(mass_earth, "Mearth", "kg"), trojan_of=(giant, point))
        if planet.planet_class is None:
            continue
        log.choice("Trojan planet", f"{point} of the giant at {giant.distance:.4g} AU", f"{mass_earth:.3g} Earth masses from {solids:.3g} in the zone (ceiling {ceiling:.3g}), class {planet.planet_class}")
        made.append(planet)
    return made
