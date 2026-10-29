# tests/audit/sweep.py

"""
The body sweep (part F of the prevalence suite, TEST.209, TEST.210, TEST.211-style rule sets): one shared set of generated planets and
moons that the consistency sweeps read, so the four rule sets (atmosphere and temperature, bulk composition, possibility tables, classes
from conditions) run over one build and not four.

`entries()` builds, for every zone (h, e, c) and each of the sweep's host stars, `PLANETS_PER_ZONE` planets at a distance in the zone,
and keeps the moons each planet drew; every body is an `Entry`. The class of a planet is derived from what it comes out as (GEN.646),
so the sweep draws bodies and lets the classes fall where they do. The draws are seeded from the host type, so a repeat gives the same
bodies.

A rule is a `Rule`: `check(entry)` returns the list of violations of the entry (empty when it holds, and when the rule does not apply),
and `breaker(entry)` returns a copy of the entry's body with exactly that rule's quantity broken, or `None` for an entry the rule does
not apply to. `run_rules` applies every rule to every entry; `prove_rules_can_fail` finds, for every rule, an entry it applies to and
checks that the broken copy is reported. A sweep with a rule that cannot fail is not a sweep.
"""

import copy
import functools
import zlib
from dataclasses import dataclass, field

HOST_TYPES = ("M2V", "G2V", "F5V")
"""tuple: The host stars of the sweep: a red dwarf, a Sun-like star and an F star."""

ZONES = "hec"

PLANETS_PER_ZONE = 14
"""int: The planets drawn in each zone round each host star."""


@dataclass
class Entry:
    """One body of the sweep: a planet (`parent` is `None`) or a moon (`parent` is its planet)."""

    body: object
    parent: object
    star: object
    label: str

    @property
    def is_moon(self):
        return self.parent is not None


@dataclass(frozen=True)
class Rule:
    """A named rule: `check(entry)` -> violations, `breaker(entry)` -> a copy of the body with that rule broken (or `None`)."""

    name: str
    check: callable
    breaker: callable
    note: str = ""


def seed_for(*parts):
    return zlib.crc32(":".join(str(part) for part in parts).encode()) or 1


@functools.lru_cache(maxsize=None)
def entries(host_types=HOST_TYPES):
    """Every body of the sweep, planets first and each planet's moons after it."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.planet import Planet
    from planetgen.generation.star import Star
    from planetgen.util import draw

    out = []
    for host in host_types:
        draw.set_run_seed(seed_for("sweep", host))
        config = SystemConfig()
        config.STAR_TYPE = host
        star = Star(config)
        inner, outer = star.habitable_zone
        for zone in ZONES:
            distance = {"h": inner * 0.5, "c": outer * 2.0}.get(zone, (inner + outer) / 2.0)
            for index in range(PLANETS_PER_ZONE):
                with draw.bound(f"sweep:{host}:{zone}:{index}"):
                    planet = Planet(SystemConfig(), star, star.habitable_zone, distance)
                label = f"{host}/{planet.planet_class}{zone}{index}"
                out.append(Entry(planet, None, star, label))
                for moon_index, moon in enumerate(planet.moons):
                    out.append(Entry(moon, planet, star, f"{label}/moon{moon_index}"))
    return tuple(out)


SYSTEM_SEEDS = 24
"""int: Whole generated systems added to the sweep: their bodies have been through the history the class-by-class builds skip."""


@functools.lru_cache(maxsize=None)
def system_entries(count=SYSTEM_SEEDS):
    """The planets and moons of `count` generated systems (seeded), as entries; belts are left out."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    out = []
    for index in range(count):
        draw.set_run_seed(seed_for("sweep-system", index))
        system = StarSystem(system_config=SystemConfig())
        for planet in system.all_planets:
            if planet.body_type == "a":
                continue
            out.append(Entry(planet, None, planet.star, f"system{index}/{planet.name}"))
            for moon in planet.moons:
                out.append(Entry(moon, planet, planet.star, f"system{index}/{moon.name}"))
    return tuple(out)


def everything():
    """The class-by-class builds and the generated systems together."""
    return entries() + system_entries()


def broken_copy(body, **changes):
    """A shallow copy of `body` with `changes` set (a dict attribute is copied, so the original is never touched)."""
    clone = copy.copy(body)
    for key, value in changes.items():
        setattr(clone, key, value)
    return clone


def run_rules(rules, sweep=None):
    """`{rule name: [(entry label, violation), ...]}` over every entry; a rule that holds everywhere has an empty list."""
    sweep = everything() if sweep is None else sweep
    found = {rule.name: [] for rule in rules}
    for entry in sweep:
        for rule in rules:
            for violation in rule.check(entry):
                found[rule.name].append((entry.label, violation))
    return found


def prove_rules_can_fail(rules, sweep=None):
    """`{rule name: True/False}`: for each rule, the first entry its breaker applies to, with the body broken, is reported by the
    rule's check. `None` for a rule no entry applies to (a sweep that never exercises it)."""
    sweep = everything() if sweep is None else sweep
    result = {}
    for rule in rules:
        result[rule.name] = None
        for entry in sweep:
            broken = rule.breaker(entry)
            if broken is None:
                continue
            clone = Entry(broken, entry.parent, entry.star, entry.label)
            result[rule.name] = bool(rule.check(clone))
            break
    return result
