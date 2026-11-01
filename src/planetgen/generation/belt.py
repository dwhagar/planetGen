# planetgen/generation/belt.py

"""
Asteroid Belt Generation
========================

This module contains the `AsteroidBelt` class, which provides a representation
for asteroid belts within a star system. It also includes utility functions
related to planet mass ranges, which are used by both planets and asteroid belts.

`format_composition_summary` is a shared, module-level extraction of `AsteroidBelt`'s own
formatting logic, reused as-is by `asteroidFieldData.AsteroidField` -- a
standalone asteroid field encountered in open space is physically the same
thing (density + mineral composition), just without a host star/orbit --
so that field-specific class delegates to these rather than duplicating
this module's logic.
"""


import functools
import weakref

from planetgen.generation.config import SystemConfig
from planetgen.physics import belt_ledger, belt_mass, belt_origin, bulk_draw, radiogenic, small_body_composition
from planetgen.util import draw
from planetgen.util import log
from planetgen.util.serialization import fields_from_dict, fields_to_dict


def format_composition_summary(composition):
    """
    Builds the human-readable composition summary sentence fragment (e.g.
    "high concentrations of iron, moderate concentrations of nickel, and
    trace amounts of platinum") for a (component, concentration) list in
    the shape `small_body_composition.describe` produces.

    Args:
        composition (list): A list of (component, concentration) tuples.

    Returns:
        str: The composition summary phrase, or "various unknown
            materials" if `composition` is empty.
    """
    composition_phrases = []
    for component, concentration in composition:
        if concentration == "trace":
            composition_phrases.append(f"trace amounts of {component}")
        else:
            composition_phrases.append(f"{concentration} concentrations of {component}")

    if not composition_phrases:
        return "various unknown materials"  # Fallback if no components are selected
    if len(composition_phrases) == 1:
        return composition_phrases[0]
    if len(composition_phrases) == 2:
        return f"{composition_phrases[0]} and {composition_phrases[1]}"
    return ", ".join(composition_phrases[:-1]) + f", and {composition_phrases[-1]}"


def _star_al26_ratio(system_config, star):
    """The aluminium-26 the star was born with (`radiogenic.star_al26_ratio`); a star with no [Fe/H] or age (a binary's proxy) takes its galaxy's mean."""
    shape = getattr(system_config, "GALAXY_SHAPE", None)
    ratio = radiogenic.star_al26_ratio(star, shape)
    return radiogenic.al26_mean(shape, getattr(star, "metallicity_fe_h", None)) if ratio is None else ratio


def present_state(star, distance, lower_limit, upper_limit, primordial_earth=None):
    """
    `(largest body km, mass kg, density word)` of a belt between `lower_limit` and `upper_limit` AU around `star` at its age: the largest body drawn from the size law (`belt_ledger.ICY_LAW` for a
    belt out past the snow line, `ROCKY_LAW` inside it) from a stream of the belt's own and settled with the mass the disc's solids leave after `star.age` of collisions (`belt_ledger.belt_state`,
    GEN.1000), from `primordial_earth` Earth masses when the giants have already acted on it (`belt_origin`, GEN.1068) or from the whole disc otherwise. `belt_birth` reads it to keep only the belts
    that are left to count.
    """
    return _present_state(star.mass, star.snow_line_au, star.age, bulk_draw.metal_ratio_to_solar(star), distance, lower_limit, upper_limit, primordial_earth)


@functools.lru_cache(maxsize=4096)
def _present_state(star_mass_kg, snow_line_au, age_gy, metal_ratio, distance, lower_limit, upper_limit, primordial_earth):
    """`present_state` on the numbers it reads (so the same belt asked for again, by the builder and by the belt itself, is settled once)."""
    law = belt_ledger.ICY_LAW if distance >= snow_line_au else belt_ledger.ROCKY_LAW
    stream = draw.Stream(f"belt-largest:{distance!r}:{lower_limit!r}:{star_mass_kg!r}")   # not the outer limit: a belt cut back keeps its own draw
    mass_kg, density, settled = belt_ledger.belt_state(lower_limit, upper_limit, star_mass_kg, snow_line_au, metal_ratio, age_gy, law, stream, primordial_earth)
    return (0.0 if settled is None else settled.largest_km), mass_kg, density


_STARS = weakref.WeakKeyDictionary()
"""The star each belt of this run was made around (`AsteroidBelt.shorten_to` reads it). Kept apart from the belt so it is never one of its stored attributes."""


class AsteroidBelt:
    """
    A basic class to store information for an asteroid belt.

    This class serves as a simple container for the properties of an asteroid
    belt, primarily its orbital distance and boundaries.

    `body_type` is always `'a'` here -- see `Planet`'s class docstring for how
    this is a shared discriminator with `Planet.body_type` (values `'t'`/`'g'`),
    used by `StarSystem` to tell the two classes apart within its single
    mixed `planets` list.

    Attributes:
        distance (float): The average distance of the asteroid belt from the star in AU.
        lower_limit (float): The inner boundary of the asteroid belt in AU.
        upper_limit (float): The outer boundary of the asteroid belt in AU.
        body_type (str): Always 'a' for an asteroid belt.
        density (str): The density of the asteroid belt ('dense', 'sparse', 'typical'): its mass per area against the median belt of its age (GEN.477).
        mass_kg (float): What its disc's solids leave after the star's age of collisions (`belt_ledger.belt_state`, [R4785]).
        chondrite_group (str): The nearest meteorite group of the mean composition (`small_body_composition.chondrite_group`, GEN.475): `CI`, `CM`, `H`, `L` ..., or `none` far from every group.
        composition (list): A list of (component, concentration) tuples -- contrast
                            `Planet.composition`, a descriptive string (same name,
                            unrelated shape, on the two classes). The words of the
                            makeup below (`small_body_composition.describe`).
        makeup (SmallBodyMakeup or None): The drawn composition the words come from (GEN.442):
                            taxon shares, grain mixture, porosity, bulk density and albedo.
                            None for a belt read back from a store that predates it.
    """

    SERIALIZABLE_FIELDS = ["distance", "lower_limit", "upper_limit", "body_type", "density", "mass_kg", "chondrite_group", "origin", "co_causes",
                           "formation_snow_line_au", "primordial_mass_earth", "cleared_mass_earth", "scattered_mass_earth", "stirred", "origin_history",
                           "origin_sentence", "scatter_fates"]
    """
    Every attribute set by `__init__`, excluding `system_config` (a shared
    back-reference), `composition` (handled separately in `to_dict`/
    `from_dict`, since it's a list of tuples -- JSON has no tuple type) and `makeup`.
    """

    origin = None
    co_causes = ()
    formation_snow_line_au = None
    primordial_mass_earth = None
    cleared_mass_earth = None
    scattered_mass_earth = None
    stirred = None
    origin_history = ()
    origin_sentence = None
    scatter_fates = None
    """Where the belt came from (GEN.1068, `physics/belt_origin.py`): the formation cause (`gap`, `tail`, `snow-line-bump`, `outer-disc`, `destroyed-planet`, `collisional-debris`), the causes that
    coincide with it, the snow line when the system formed, the mass the belt formed with, what the giants cleared from it and what an instability scattered (Earth masses), whether a giant stirs it,
    the sentences of what shaped it, its sentence for the page and, for an outer disc, the shares the scatter step gave each population (`physics/scatter.py`). All empty for a belt a `SLOTS` entry asked for."""

    def __init__(self, system_config: SystemConfig, star, distance, lower_limit, upper_limit, origin=None, scatter_fates=None): # Added system_config
        """
        Initializes an AsteroidBelt object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object for the system.
            star: The star (or binary proxy) the belt orbits; its luminosity and the aluminium-26 it was born with decide the taxa (GEN.473).
            distance (float): The average distance from the star in AU.
            lower_limit (float): The inner boundary of the belt in AU.
            upper_limit (float): The outer boundary of the belt in AU.
            origin (belt_origin.BeltOrigin, optional): Where the belt came from; its mass after the giants is the mass the belt begins to grind from (GEN.1068).
            scatter_fates (dict, optional): The shares of an outer disc in each population after the scatter step (`scatter.scatter_shares`).
        """
        self.system_config = system_config # Store system_config
        self.distance = distance
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        self.body_type = 'a'
        largest_km, self.mass_kg, self.density = present_state(star, distance, lower_limit, upper_limit, None if origin is None else origin.mass_earth)
        _STARS[self] = star       # what `shorten_to` needs; not stored, a belt read back from a store is not shortened
        if origin is not None:
            self.origin, self.co_causes, self.formation_snow_line_au = origin.origin, list(origin.co_causes), origin.formation_snow_line_au
            self.primordial_mass_earth, self.cleared_mass_earth, self.scattered_mass_earth = origin.primordial_mass_earth, origin.cleared_earth, origin.scattered_earth
            self.stirred, self.origin_history, self.origin_sentence = origin.stirred, [event.text for event in origin.history], origin.sentence
            self.scatter_fates = scatter_fates
        log.choice("Asteroid belt density", self.density, f"{self.mass_kg:.4g} kg left of the disc's solids after {star.age:.4g} Gy, its largest body {largest_km:.4g} km (stand-in)")
        self.makeup = small_body_composition.draw_belt(distance, star.luminosity, _star_al26_ratio(system_config, star))
        self.composition = small_body_composition.describe(self.makeup)
        self.chondrite_group = self.makeup.detail["group"]

    def shorten_to(self, upper_limit):
        """
        Cuts the belt back to end at `upper_limit` AU (the orbits the star holds, once its age is final): the mass it formed with, what the giants cleared and what was scattered shrink with
        the disc's solids between the limits, and the mass left after the star's age is worked out again. Returns False, and changes nothing, for a belt that was not made in this run or
        whose stretch would be narrower than `belt_origin.MIN_BELT_RATIO`.
        """
        star = _STARS.get(self)
        if star is None or self.origin is None or upper_limit >= self.upper_limit:
            return self.upper_limit <= upper_limit
        if upper_limit < self.lower_limit * belt_origin.MIN_BELT_RATIO:
            return False
        metal = bulk_draw.metal_ratio_to_solar(star)
        solids = lambda upper: belt_mass.primordial_mass_earth(self.lower_limit, upper, star.mass, star.snow_line_au, metal)
        share = solids(upper_limit) / solids(self.upper_limit)
        self.primordial_mass_earth *= share
        self.cleared_mass_earth *= share
        self.scattered_mass_earth *= share
        self.upper_limit = upper_limit
        left_earth = self.primordial_mass_earth - self.cleared_mass_earth - self.scattered_mass_earth
        _, self.mass_kg, self.density = present_state(star, self.distance, self.lower_limit, upper_limit, left_earth)
        return True

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this belt's properties, per
        `SERIALIZABLE_FIELDS`, with `composition` stored as a list of
        2-element lists (JSON has no tuple type).

        Returns:
            dict: One entry per field in `SERIALIZABLE_FIELDS`, plus
                 `composition`.
        """
        data = fields_to_dict(self, self.SERIALIZABLE_FIELDS)
        data["composition"] = [list(pair) for pair in self.composition]
        data["makeup"] = self.makeup.to_dict() if self.makeup is not None else None
        return data

    @classmethod
    def from_dict(cls, data, system_config):
        """
        Reconstructs an `AsteroidBelt` from a dict in the shape `to_dict()`
        produces, without re-running generation (`__init__` is bypassed via
        `object.__new__`).

        Args:
            data (dict): A dict in the shape `to_dict()` produces.
            system_config (SystemConfig): The system's shared config.

        Returns:
            AsteroidBelt: The reconstructed belt.
        """
        belt = object.__new__(cls)
        belt.system_config = system_config
        fields_from_dict(belt, data, cls.SERIALIZABLE_FIELDS)
        belt.composition = [tuple(pair) for pair in data["composition"]]
        belt.makeup = small_body_composition.SmallBodyMakeup.from_dict(data["makeup"]) if data.get("makeup") else None
        return belt

    def get_composition_summary(self):
        """
        Builds the human-readable composition summary sentence fragment
        (e.g. "high concentrations of iron, moderate concentrations of
        nickel, and trace amounts of platinum") as its own method so the
        database persistence layer can store the same as-published text
        (`asteroid_belts.composition_summary`) without duplicating this
        formatting logic (see `planetgen/db/store.py`). Belts have no
        properties-dict data table the way stars/planets do, so this is
        their equivalent "searchable, as-published" text.

        Returns:
            str: The composition summary phrase, or
                "various unknown materials" if no components were selected.
        """
        return format_composition_summary(self.composition)

from planetgen.generation import object_timing  # noqa: E402 -- registers the generators below (PERF.92)
object_timing.track(AsteroidBelt, "belt")
