# planetgen/physics/delivered_water.py

"""
The water a rocky body is born with, from where it formed and how heavy it is (GEN.849; docs/design/composition/volatile-delivery-formation-distance.md)
==========================================================================================================================================================

The Bern population synthesis (`data/ngpps_formation_water.json`, GEN.802) says how much of a planet's mass is water against the distance
it formed at, in units of its star's snow line (`x = a_form / a_snow`), and its mass: the share of the planets that are `dry` (under 0.025
percent water), `low` (to 0.16 percent), `high` (to 10 percent) or `rich` (over 10 percent, the hycean limit). The atmosphere recipes b1 to b4
already read that law for the water an air is delivered. The body's own water used to be a second, unrelated draw (a log-uniform range inside
the snow line and 3 percent of worlds with an ocean's worth); now the body takes its class from the same law and a water share inside the
class, so the air and the body agree.

1. `class_shares(mass_earth, snow_ratio)` is the law: the recipes' rows b1 to b4 as `{class: share}` (the mass held to the population's 0.5 to
   30 Earth masses, the distance to its 0.03 to 4 snow lines, where the law was fitted).
2. `sample_fraction(mass_earth, snow_ratio)` draws the class from those shares, then the water share log-uniform inside the class: from
   `tuning.DELIVERED_WATER_DRY_MIN` to the first class edge for `dry`, between the edges for `low` and `high`, from the hycean limit to
   `tuning.DELIVERED_WATER_RICH_MAX` for `rich`. Two draws, always.

The class edges are the model's (`classes.edges_f`: the nearest mapping of Ortenzi's 0, 0.05 and 0.5 percent cases); the brackets inside a
class and the two ends are stand-ins (registry rows DW1 to DW3, `docs/design/research/gap-filling-plan.md` section 6.41).
"""

import json
import os

from planetgen import tuning
from planetgen.generation.possibility_table import Subject
from planetgen.physics import elements
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.util import draw
from planetgen.util.random import log_uniform

with open(os.path.join(elements.DATA_DIRECTORY, "ngpps_formation_water.json"), encoding="utf-8") as _handle:
    _MODEL = json.load(_handle)

CLASSES = tuple(_MODEL["classes"]["names"])
"""tuple: The water classes, driest first (`dry`, `low`, `high`, `rich`)."""

EDGES = tuple(_MODEL["classes"]["edges_f"])
"""tuple: The water mass fractions between the classes (the model's own)."""

ROWS = dict(zip(("b1", "b2", "b3", "b4"), CLASSES))
"""dict: The recipe row that is each class."""

MASS_RANGE_EARTH = tuple(_MODEL["model"]["mass_range_earth"])
X_RANGE = tuple(_MODEL["model"]["x_range"])
"""tuple: The population's masses (Earth masses) and distances (snow lines), where the law was fitted."""


def _clamp(value, low, high):
    return min(max(value, low), high)


def class_shares(mass_earth, snow_ratio):
    """`{class: share}` of the bodies of `mass_earth` Earth masses formed at `snow_ratio` snow lines, adding up to 1 (the law of the recipes b1 to b4)."""
    mass = _clamp(mass_earth, *MASS_RANGE_EARTH)
    subject = Subject("terrestrial", {"mass_earth": mass, "snow_ratio": _clamp(snow_ratio, *X_RANGE), "const:b:c_b4": 1.0,
                                      "core_radius_earth": mass ** 0.27})
    table = dict(RECIPES.probabilities("b", subject))
    return {ROWS[row]: table.get(row, 0.0) for row in ROWS}


def bracket(water_class):
    """`(low, high)` of the water mass fraction a class draws from."""
    index = CLASSES.index(water_class)
    low = tuning.DELIVERED_WATER_DRY_MIN if index == 0 else EDGES[index - 1]
    high = EDGES[index] if index < len(EDGES) else tuning.DELIVERED_WATER_RICH_MAX
    return low, high


def sample_fraction(mass_earth, snow_ratio, rng=draw):
    """
    The water mass fraction of a rocky body of `mass_earth` Earth masses that formed at `snow_ratio` snow lines: the class from
    `class_shares`, then a log-uniform share inside it. Two draws (a uniform for the class, one for the share).
    """
    shares = class_shares(mass_earth, snow_ratio)
    pick = rng.random()
    running = 0.0
    chosen = CLASSES[-1]
    for name in CLASSES:
        running += shares[name]
        if pick < running:
            chosen = name
            break
    return log_uniform(*bracket(chosen), rng=rng)
