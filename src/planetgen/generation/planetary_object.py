# planetgen/generation/planetary_object.py

"""
Planetary Object
================

The family class of the bodies that are planets by size and nature (GEN.565; the `PlanetaryObject` of
docs/design/architecture/mass-object-class.md, section 10): a planet bound to a star, a moon bound to a planet, and a
rogue planet bound to nothing. What they share is what a planet *is*: a class letter, a bulk-composition description,
whether it is a gas giant, and what it orbits. A subclass keeps its own draw, fields and stored names.

A moon is not a class of its own: it is a planetary object whose parent is a planet. A rogue planet is one with no parent.

Nothing here draws a random number.
"""

from planetgen.generation.mass_object import MassObject


class PlanetaryObject(MassObject):
    """A planet, a moon or a rogue planet."""

    GIANT_LETTER = "g"
    """str: The `body_type` / `planet_type` letter of a gas giant (the other kind is `'t'`)."""

    planet_class = None
    """str or None: Its `PLANET_CLASSES` letter."""
    composition = None
    """str or None: A descriptive bulk-composition string."""
    is_moon = False
    """bool: Whether its parent is a planet."""

    @property
    def giant_kind(self):
        """The `'g'` or `'t'` letter the family stores under its own name (`body_type`, `planet_type`)."""
        raise NotImplementedError

    @property
    def is_giant(self):
        """Whether it is a gas giant (a brown dwarf drawn as a rogue planet counts, as its `planet_type` does)."""
        return self.giant_kind == self.GIANT_LETTER

    @property
    def parent_kind(self):
        """What it orbits: `"star"` for a planet, `"planet"` for a moon, `None` for a rogue planet."""
        raise NotImplementedError
