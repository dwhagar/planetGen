# tests/audit/ledger_storage_census.py

"""
GEN.1465: how many effect-ledger rows a system would store at each floor, measured on generated systems.

The effect code does not fill the ledger yet (each lane's item calls `system.effect_ledger.add` as it lands), so this audit stands in a *proxy* contribution
for each of the seven effects from the bodies' stored masses, sizes and orbits, with the hierarchical separations of the generator (a moon about its planet,
a planet about its host, a companion star at the binary separation). The proxy has the right dependence on distance and mass for each effect (light and UV
over d^2, tide over d^3, heating over the square of the tide, secular, precession and eccentricity from the mass ratio and the semi-major axis ratio), which is
what sets how many sources lie within a factor of 1e2, 1e3 or 1e6 of the largest; it does not give the value of any effect. Run:

    python -m tests.audit.ledger_storage_census 300 [first_seed]

and it prints the mean rows a system, bytes a system and gigabytes at the 200 pc default for each floor.
"""

import math
import statistics
import sys

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

AU_KM = 1.495978707e8
FLOORS = (("largest only", 1.0), ("1e-2", 1e-2), ("1e-3", 1e-3), ("1e-6", 1e-6), ("full", 1e-37))
ROW_BYTES = 56
SYSTEMS_200_PC = 7.45e6
EFFECTS = ("light", "uv", "tide", "heating", "secular", "precession", "eccentricity")


class Body:
    """A body with its place in the hierarchy: `host` the star index it orbits (None for a star), `parent` the planet a moon orbits, `r` its distance from what it orbits (AU)."""

    def __init__(self, obj, kind, host, r, mass, radius_km, parent=None, luminosity=0.0, temperature=5772.0, index=None):
        self.obj, self.kind, self.host, self.r, self.mass, self.radius_km = obj, kind, host, r, mass, radius_km
        self.index = index
        self.parent, self.luminosity, self.temperature = parent, luminosity, temperature


def _b_half(alpha):
    """The Laplace coefficient b_{3/2}^{(1)}(alpha) by 64-point quadrature (alpha under 1)."""
    alpha = min(alpha, 0.99)
    total = 0.0
    for k in range(64):
        psi = (k + 0.5) * math.pi / 64
        total += math.cos(psi) / (1 - 2 * alpha * math.cos(psi) + alpha * alpha) ** 1.5
    return 2 * total / 64


def _g2(ra, rb):
    """Mean 1/d^2 between two bodies on circles of radii ra and rb about the same center: 1/(R^2 - d^2), held at 5 percent of R^2."""
    big, small = max(ra, rb), min(ra, rb)
    return 1.0 / max(big * big - small * small, 0.05 * big * big)


def bodies_of(system):
    stars = list(system.stars)
    out = []
    sep = 0.1
    if system.binary_type == "wide":
        sep = system.wide_binary.separation_au
    elif system.binary_type is None and len(stars) > 1:
        sep = 100.0
    for index, star in enumerate(stars):
        out.append(Body(star, "star", None, 0.0, star.mass, star.radius, luminosity=star.luminosity, temperature=getattr(star, "temperature", 5772.0) or 5772.0, index=index))
    hosts = [(system.planets, 0)]
    if system.secondary_planets:
        hosts.append((system.secondary_planets, 1 if len(stars) > 1 else 0))
    for i, group in enumerate(system.member_planets):
        hosts.append((group, min(2 + i, len(stars) - 1)))
    for group, host in hosts:
        for planet in group:
            kind = "belt" if planet.body_type == "a" else "planet"
            mass = getattr(planet, "mass_kg", None) or getattr(planet, "mass", 0.0)
            radius = 0.0 if kind == "belt" else planet.radius
            parent = Body(planet, kind, host, planet.distance, mass, radius)
            out.append(parent)
            for moon in getattr(planet, "moons", None) or []:
                out.append(Body(moon, "moon", host, moon.distance, moon.mass, moon.radius, parent=parent))
    for comet in list(system.comets) + list(system.secondary_comets):
        out.append(Body(comet, "comet", 0, comet.distance_au, comet.mass_kg, 1.0))
    return out, sep


def _d2g(a, b, sep):
    """(mean 1/d^2, mean 1/d^3) between bodies a and b, by the generator's hierarchy."""
    if a.kind == "star" or b.kind == "star":
        star, other = (a, b) if a.kind == "star" else (b, a)
        if other.kind == "star":
            r = sep
        elif other.kind == "moon":
            r = other.parent.r
        else:
            r = other.r
        if other.kind != "star" and star is not _host_star_of(other, star):
            r = max(sep, r)
        g2 = 1.0 / (r * r)
        return g2, g2 ** 1.5
    if a.parent is b or b.parent is a:
        r = a.r if a.parent is b else b.r
        g2 = 1.0 / (r * r)
        return g2, g2 ** 1.5
    if a.parent is not None and a.parent is b.parent:
        g2 = _g2(a.r, b.r)
        return g2, g2 ** 1.5
    if a.host != b.host:
        g2 = 1.0 / (sep * sep)
        return g2, g2 ** 1.5
    ra = a.parent.r if a.parent is not None else a.r
    rb = b.parent.r if b.parent is not None else b.r
    g2 = _g2(ra, rb)
    return g2, g2 ** 1.5


def _host_star_of(body, star):
    """`star` when it is the star `body` orbits (by index), else None; stars are compared by their place in the system's list."""
    return star if getattr(star, "index", None) == body.host else None


def _secular(a, b, bodies, sep):
    """The proxy of how strongly b forces a's orbit (mass ratio times a function of the semi-major axis ratio), 0 when they do not share a central body."""
    if a.kind == "moon":
        center_mass = a.parent.mass
        if b.parent is a.parent and b.kind == "moon":
            lo, hi = sorted((a.r, b.r))
            alpha = lo / hi
            return (b.mass / center_mass) * alpha * alpha * _b_half(alpha) / 4.0
        if b is a.parent:
            return 0.0
        outer = a.parent.r if (b.kind == "star" and getattr(b, "index", None) == a.host) else max(b.r, a.parent.r, sep if b.kind == "star" else 0.0)
        return (b.mass / center_mass) * (a.r / outer) ** 3 if outer > 0 else 0.0
    host_star = next((x for x in bodies if x.kind == "star" and x.index == a.host), bodies[0])
    center_mass = host_star.mass
    if b.kind == "star":
        if b is host_star:
            return 0.0
        return (b.mass / center_mass) * (a.r / sep) ** 3
    if b.kind == "moon":
        return 0.0
    if b.host != a.host:
        return (b.mass / center_mass) * (a.r / sep) ** 3
    lo, hi = sorted((a.r, b.r))
    if hi <= 0 or lo <= 0:
        return 0.0
    alpha = lo / hi
    return (b.mass / center_mass) * alpha * alpha * _b_half(alpha) / 4.0


def contributions(system):
    """`{effect: {receiver index: [value, ...]}}` of the proxy contributions."""
    bodies, sep = bodies_of(system)
    out = {name: {} for name in EFFECTS}
    host_mass = {}
    for i, body in enumerate(bodies):
        if body.kind == "star":
            host_mass[i] = body.mass
    for i, a in enumerate(bodies):
        light, uv, tide, heating, secular, precession, eccentricity = ([] for _ in EFFECTS)
        for j, b in enumerate(bodies):
            if i == j or a.mass <= 0:
                continue
            g2, g3 = _d2g(a, b, sep)
            if b.kind == "star":
                light.append(b.luminosity * g2)
                uv.append(b.luminosity * (b.temperature / 5772.0) ** 2 * g2)
            if b.mass > 0:
                value = b.mass * max(a.radius_km, 1.0) * g3
                tide.append(value)
                heating.append(value * value)
            if a.kind != "star":
                value = _secular(a, b, bodies, sep)
                if value > 0:
                    secular.append(value)
                    precession.append(value)
                    eccentricity.append(value * min(a.r, b.r if b.r > 0 else a.r) / max(a.r, b.r, 1e-12))
        for name, values in zip(EFFECTS, (light, uv, tide, heating, secular, precession, eccentricity)):
            values = [v for v in values if v > 0]
            if values:
                out[name][i] = values
    return out


def rows_at(contribs, floor, only=None):
    """Rows stored at `floor`: each source over the floor (the largest always) and one summary row for each (body, effect); `only` names one effect."""
    total = 0
    for name, by_body in contribs.items():
        if only is not None and name != only:
            continue
        for values in by_body.values():
            top = max(values)
            total += 1 + sum(1 for v in values if v != top and v / top >= floor) + 1
    return total


def main(count=300, first=1000):
    rows = {label: [] for label, _ in FLOORS}
    per_effect = {name: {label: [] for label, _ in FLOORS} for name in EFFECTS}
    bodies = []
    for seed in range(first, first + count):
        with draw.bound(seed):
            system = StarSystem(system_config=SystemConfig())
        contribs = contributions(system)
        bodies.append(len(bodies_of(system)[0]))
        for label, floor in FLOORS:
            rows[label].append(rows_at(contribs, floor))
            for name in EFFECTS:
                per_effect[name][label].append(rows_at(contribs, floor, only=name))
    print(f"{count} default systems, seeds {first}.., mean bodies {statistics.mean(bodies):.2f}")
    print("floor           mean rows   p90 rows   max rows   KB/system   GB at 200 pc (7.45 million systems)")
    for label, _ in FLOORS:
        data = rows[label]
        mean = statistics.mean(data)
        p90 = sorted(data)[int(0.9 * len(data))]
        print(f"{label:<14} {mean:>10.1f} {p90:>10} {max(data):>10} {mean * ROW_BYTES / 1024:>11.2f} {mean * ROW_BYTES * SYSTEMS_200_PC / 1e9:>12.1f}")
    print("mean rows a system by effect at each floor:")
    print("effect         " + " ".join(f"{label:>13}" for label, _ in FLOORS))
    for name in EFFECTS:
        print(f"{name:<14} " + " ".join(f"{statistics.mean(per_effect[name][label]):>13.1f}" for label, _ in FLOORS))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 300, int(sys.argv[2]) if len(sys.argv) > 2 else 1000)
