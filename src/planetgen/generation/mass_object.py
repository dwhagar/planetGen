# planetgen/generation/mass_object.py

"""
Mass Object
===========

The core every body with a mass shares (GEN.565, the class the unified
star-like objects design calls `MassObject`): the galactic motion of a body that is not bound to a star's frame, and the
one `to_dict`/`from_dict` pair that each family used to write out by hand.

A family class (`Star`, `Planet`, `Comet`, `RoguePlanet`, ...) lists its
own `SERIALIZABLE_FIELDS` and, where a value is not plain JSON, overrides
`_encode` (object to dict) and `_decode` (dict back onto the object). Nothing
here draws a random number except `set_galactic_motion`, which draws what
the older copies drew, in the same place in each constructor.
"""

from planetgen.galaxy.galactic_orbit import generate_galactic_orbit_fields
from planetgen.util.serialization import fields_from_dict, fields_to_dict

GALACTIC_MOTION_FIELDS = [
    "galactic_orbital_speed_kms", "galactic_orbital_period_gy",
    "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
]
"""The four fields `set_galactic_motion` sets: the circular speed and period about the galactic centre, the phase
along that orbit, and the shortest interval between orbital updates."""


class MassObject:
    """The core of a body with a mass: galactic motion and serialization from a field list. A family whose name is
    a salad drawn on first read lists `LazySaladName` before it."""

    SERIALIZABLE_FIELDS = []
    """The attributes a family persists; each family sets its own."""

    def set_galactic_motion(self, galactic_center_dist_ly=None, galactic_orbital_phase_deg=None):
        """Sets the four `GALACTIC_MOTION_FIELDS` from the body's distance from the galactic centre (none: the
        generator's default distance) and a given phase (none: drawn)."""
        (self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy,
         self.galactic_orbital_phase_deg, self.galactic_min_update_interval_years) = \
            generate_galactic_orbit_fields(galactic_center_dist_ly, galactic_orbital_phase_deg)

    def to_dict(self):
        """A JSON-serializable dict of the body's `SERIALIZABLE_FIELDS`, with `_encode` applied."""
        data = fields_to_dict(self, self.SERIALIZABLE_FIELDS)
        self._encode(data)
        return data

    def _encode(self, data):
        """Turns the values of `data` that are not plain JSON into ones that are; the family overrides."""

    @classmethod
    def _rebuild(cls, data, system_config, **context):
        """The body for `data`, without running `__init__` (that would draw it again): its shared config, any
        `context` attributes (what a family threads in rather than stores, such as a planet's star), the
        `SERIALIZABLE_FIELDS` set, and `_decode` applied."""
        body = object.__new__(cls)
        body.system_config = system_config
        body.__dict__.update(context)
        fields_from_dict(body, data, cls.SERIALIZABLE_FIELDS)
        body._decode(data)
        return body

    def _decode(self, data):
        """Turns the values `_encode` changed back; the family overrides."""

    @classmethod
    def from_dict(cls, data, system_config):
        """Reconstructs the body from a dict in the shape `to_dict` produces."""
        return cls._rebuild(data, system_config)
