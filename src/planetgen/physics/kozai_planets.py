# planetgen/physics/kozai_planets.py

"""
Kozai-Lidov pumping of a planet by the distant bodies of its system (GEN.1447, the planet and sibling-planet cases; slice 0, pure and unwired)
===============================================================================================================================================

`kozai_bodies` screens a moon as the inner body of the triple (its planet, the moon, the planet's star). This module sets up the same screen for a **planet**: the inner pair is the planet and its
host star, and the third body is any body farther out whose orbit is inclined to the planet's by 39.2 to 140.8 degrees:

* a **companion star** or a dark member of the tree, at its axis and eccentricity about the host (`radiation_field.Companion`), on the plane of the orbit node that joins it to the host;
* an **outer sibling planet**, a giant that outweighs the planet and lies past three times its axis (`kozai_bodies.HIERARCHY_RATIO`), on its own plane.

Nothing is drawn from the run's stream and nothing is stored or changed: `screen_planets` takes plain numbers (`PlanetInput`, `Perturber`) and returns, for each planet, the `kozai_bodies.Kozai` of the
perturber that pulls it fastest, the list of the ones that cycle, and the ledger rate of every perturber (1 / Kozai time, effect 5). The wiring (a later slice, a lock item) applies the outcome
to the planet's eccentricity and inclination at the star's age and drops a planet whose maximum puts its periapsis inside the star or its Roche limit.

* **The other precessions** a cycle has to beat are the siblings' and relativity's: `secular_chain.precession_rates` of the star's planets, the rate that function was built to give this screen (GEN.1450);
  an outer sibling that is itself a perturber is in that rate too, so it quenches a cycle of its own pull exactly as the planet-planet theory says.
* **The starting argument of periapsis** is drawn from a stream of the planet's own numbers (`kozai-planet:` plus its mass, axis and the host's mass), so a planet made again is the same and no other draw moves.
* **Which perturber:** the screen is run for each, and the planet takes the one that cycles with the shortest Kozai time (`cycling_perturbers` lists the rest); with none cycling the planet is quiet.

Registered stand-ins: the reuse of KZ1 to KZ6 (section 6.262) for the triple, and KZP1 to KZP2 (section 6.279) for the choice of the strongest perturber and the outer-sibling mass rule.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, kozai_bodies, rings, secular_chain, stability
from planetgen.util import draw

SIBLING_MASS_RATIO = 1.0
"""float: An outer sibling counts as a perturber only when it is at least this many times the planet's mass (a heavier outer body; a lighter one is a test particle's equal and quenches rather than pumps). [Stand-in KZP2, tier "educated guess from physics", tweak L1.]"""


@dataclass(frozen=True)
class PlanetInput:
    """A planet as the screen reads it: a `key` (any hashable), its mass (kg) and radius (km), its axis (AU) and eccentricity, and its plane (inclination and node, radians, in the system's frame)."""

    key: object
    mass_kg: float
    radius_km: float
    axis_au: float
    eccentricity: float
    inclination_rad: float
    node_rad: float


@dataclass(frozen=True)
class Perturber:
    """A body farther out than the planets it pulls: a `key`, its `kind` (`companion` or `sibling`), mass (solar masses), axis (AU) about the host, eccentricity and plane (radians)."""

    key: object
    kind: str
    mass_solar: float
    axis_au: float
    eccentricity: float
    inclination_rad: float
    node_rad: float


@dataclass(frozen=True)
class PlanetKozai:
    """
    What the screen made of one planet.

    Attributes:
        key: The planet's key.
        outcome (kozai_bodies.Kozai): The pull that matters most: the shortest-Kozai-time cycle, else the quiet screen with the shortest Kozai time (None when no perturber reached the screen).
        perturber: The key of that perturber (None when there is none).
        cycling_perturbers (tuple): The keys of every perturber whose pull runs a cycle, shortest Kozai time first.
        ledger_rates (dict): `{perturber key: 1 / Kozai time}` per year, for every perturber that passed the hierarchy test (the quenched ones too, per design rule 6).
    """

    key: object
    outcome: object
    perturber: object
    cycling_perturbers: tuple
    ledger_rates: dict


def _mass_solar(kg):
    return kg / constants.SOLAR_MASS_TO_KG


def sibling_perturbers(planets, index):
    """The planets of `planets` that pull `planets[index]` as a Kozai perturber: farther than the hierarchy ratio and at least `SIBLING_MASS_RATIO` times its mass."""
    inner = planets[index]
    return [Perturber(other.key, "sibling", _mass_solar(other.mass_kg), other.axis_au, other.eccentricity, other.inclination_rad, other.node_rad)
            for other in planets
            if other.axis_au >= inner.axis_au * kozai_bodies.HIERARCHY_RATIO and other.mass_kg >= SIBLING_MASS_RATIO * inner.mass_kg]


def screen_planets(host_mass_kg, host_radius_km, age_yr, planets, companions=()):
    """
    The Kozai-Lidov outcome of every planet of one star.

    Args:
        host_mass_kg, host_radius_km (float): The star the planets orbit.
        age_yr (float): The system's age, years.
        planets (sequence of PlanetInput): The star's planets (belts and moons left out), in any order.
        companions (sequence of Perturber): The companions the star's planets see, with the plane of each one's orbit (the sibling perturbers are found here).

    Returns:
        list: One `PlanetKozai` per planet, in the order given.
    """
    if not host_mass_kg > 0.0 or not host_radius_km > 0.0 or not age_yr > 0.0:
        raise ValueError("the host mass, the host radius and the age must be positive")
    planets = list(planets)
    if not planets:
        return []
    host_solar = _mass_solar(host_mass_kg)
    apsidal = secular_chain.precession_rates([p.axis_au for p in planets], [p.mass_kg for p in planets], host_mass_kg)[0]
    host_density = host_mass_kg / (4.0 / 3.0 * math.pi * (host_radius_km * 1000.0) ** 3)
    results = []
    for index, planet in enumerate(planets):
        density = planet.mass_kg / (4.0 / 3.0 * math.pi * (planet.radius_km * 1000.0) ** 3)
        extent_km = max(host_radius_km + planet.radius_km, rings.roche_limit_km(host_radius_km, host_density, density))
        extent_au = extent_km * 1000.0 / constants.AU_TO_M
        w0 = draw.Stream(f"kozai-planet:{float(planet.mass_kg)!r}:{float(planet.axis_au)!r}:{float(host_mass_kg)!r}").uniform(0.0, 2.0 * math.pi)
        quiet = []
        cycling = []
        rates = {}
        for perturber in (*companions, *sibling_perturbers(planets, index)):
            if perturber.axis_au / planet.axis_au < kozai_bodies.HIERARCHY_RATIO:
                continue
            mutual = stability.mutual_inclination(planet.inclination_rad, planet.node_rad, perturber.inclination_rad, perturber.node_rad)
            outcome = kozai_bodies.screen(host_solar, _mass_solar(planet.mass_kg), perturber.mass_solar, planet.axis_au, planet.eccentricity, perturber.axis_au,
                                          perturber.eccentricity, mutual, age_yr, extent_au, w0_rad=w0, other_rate_rad_yr=apsidal[index])
            if outcome.timescale_yr > 0.0:
                rates[perturber.key] = 1.0 / outcome.timescale_yr
            (cycling if outcome.cycles else quiet).append((outcome.timescale_yr, perturber.key, outcome))
        cycling.sort(key=lambda item: item[0])
        quiet.sort(key=lambda item: item[0])
        chosen = (cycling or quiet or [(None, None, None)])[0]
        results.append(PlanetKozai(planet.key, chosen[2], chosen[1], tuple(key for _t, key, _o in cycling), rates))
    return results
