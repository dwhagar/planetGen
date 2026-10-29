# planetgen/generation/class_sample.py

"""
A sample that holds every class (GEN.791)
=========================================

Random generation hardly ever draws some classes (under 30 planets of O, P, T, V, E, K, L, N, F and G in 400 systems), so a
check that needs a number of bodies of every class (the property sheet fill audit, `describe/sheet.py`, UX.208) cannot
use a plain sample. `fill` tops a sample up with systems generated under `SystemConfig.FORCE_CLASS`, which gives each one at
least one planet of the asked class, until every class has enough planets or the tries run out. A test and admin helper:
default generation never sets the switch.

Use:

    systems, counts, short = class_sample.fill(per_class=30, seed=7)
"""

import collections

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

TRIES_PER_BODY = 12
"""int: How many forced systems `fill` may generate for each body it still needs of a class before it gives up on it."""


def planet_counts(systems):
    """`Counter` of the planet classes over the planets of `systems` (belts and comets have none)."""
    counts = collections.Counter()
    for system in systems:
        counts.update(planet.planet_class for planet in system.planets if hasattr(planet, "planet_class"))
    return counts


def fill(per_class=30, seed=0, letters=None, systems=(), tries_per_body=TRIES_PER_BODY):
    """
    Generates systems until `systems` plus the new ones hold `per_class` planets of every class in `letters`.

    Args:
        per_class (int): The number of planets of each class to hold.
        seed (int): The seed of the stream the new systems are drawn from (the same arguments give the same systems).
        letters (iterable or None): The classes to fill; `None` means every planet class.
        systems (iterable): Systems already sampled; they count, and are returned first.
        tries_per_body (int): The tries (forced systems that kept no planet of the class are not counted as bodies)
            allowed per planet a class still needs.

    Returns:
        tuple: `(systems, counts, short)`: the list of systems, the `Counter` of planet classes over them and the
            letters that did not reach `per_class` within the tries.
    """
    letters = sorted(letters if letters is not None else tuning.PLANET_CLASSES)
    systems = list(systems)
    counts = planet_counts(systems)
    short = []
    with draw.bound(draw.Stream(seed)):
        for letter in letters:
            budget = max(1, per_class - counts[letter]) * tries_per_body
            while counts[letter] < per_class and budget > 0:
                budget -= 1
                config = SystemConfig()
                config.FORCE_CLASS = letter
                system = StarSystem(system_config=config)
                if system.unmet_requirements:
                    continue
                systems.append(system)
                counts.update(planet.planet_class for planet in system.planets if hasattr(planet, "planet_class"))
            if counts[letter] < per_class:
                short.append(letter)
    return systems, counts, short
