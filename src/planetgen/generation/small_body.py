# planetgen/generation/small_body.py

"""
Small Body
==========

The family class of the small bodies that have no atmosphere to speak of and are described by a size and a
composition list (GEN.565e; the `SmallBody` of docs/design/architecture/mass-object-class.md, section 10): a comet bound to a star,
a comet passing through on an unbound path, and a field of asteroids. What the three share is the preamble of their
constructors (the shared configuration, a name drawn on first read, whether the name was given) and the save and load of
a composition kept as a list of strings. A subclass keeps its own draw and its own fields.

Nothing here draws a random number.
"""

from planetgen.generation.mass_object import MassObject
from planetgen.names.wordsalad import LazySaladName
from planetgen.physics.small_body_composition import SmallBodyMakeup


class SmallBody(LazySaladName, MassObject):
    """A comet, an interstellar comet or an asteroid field: a small body with a size and a composition."""

    composition = None
    """list or None: What it is made of, as a list of component strings (an asteroid field keeps a list of
    `(component, share)` pairs, and overrides `_encode` and `_decode`)."""

    def __init__(self, system_config, name=None):
        """
        Args:
            system_config (SystemConfig): The shared SystemConfig object (its settings).
            name (str, optional): An explicit name. Random if omitted, drawn on first read (PERF.43); a given name is
                kept over an object ID (GEN.64).
        """
        self.system_config = system_config
        self._init_salad_name(name)  # drawn on first read (PERF.43)
        self.name_given = bool(name)  # a given name is kept over an object ID (GEN.64)

    makeup = None
    """SmallBodyMakeup or None: What it is made of as the draw made it (GEN.442): the grain mixture, the porosity, the bulk
    density. `composition` holds the words taken from it. None for a body read back from a store that has not got it."""

    def _encode(self, data):
        data["composition"] = list(self.composition)
        data["makeup"] = self.makeup.to_dict() if self.makeup is not None else None

    def _decode(self, data):
        self.composition = list(data["composition"])
        self.makeup = SmallBodyMakeup.from_dict(data["makeup"]) if data.get("makeup") else None
