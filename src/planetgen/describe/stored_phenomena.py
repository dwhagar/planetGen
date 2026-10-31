# planetgen/describe/stored_phenomena.py

"""
The Document of a stored phenomenon (UX.185). `describe_phenomenon` (`describe/phenomena.py`) takes a freshly generated
object; a phenomenon page only has the row `db.query.phenomenon_detail` answers. `stored_phenomenon_document` wraps that
row so the same builders read it, filling in what the row does not hold but its columns give (the speed round the galaxy
from the velocity, a rogue planet's kind from its type and mass bin, a quasar's brightness against a Milky Way). The
words come from the builders, so a stored page and a generated one cannot say different things.
"""

import math

from planetgen import tuning
from planetgen.describe import build, model, phenomena
from planetgen.physics import brown_dwarf as brownDwarf, unit_registry

_msun_to_kg = unit_registry.converter("Msun", "kg")


TYPES = ("nebula", "quasar", "rogue_planet", "interstellar_comet", "asteroid_field", "supernova_remnant", "black_hole",
         "neutron_star")
"""tuple: The phenomenon types with a paragraph."""

_BUILDERS = {
    "nebula": phenomena._nebula, "quasar": phenomena._quasar, "rogue_planet": phenomena._rogue_planet,
    "interstellar_comet": phenomena._interstellar_comet, "asteroid_field": phenomena._asteroid_field,
    "supernova_remnant": phenomena._supernova_remnant,
    "black_hole": lambda row: model.section("black-hole", None, row.name, build._black_hole_blocks(row)),
    "neutron_star": lambda row: model.section("neutron-star", None, row.name, build._neutron_star_blocks(row)),
}
"""dict: The section builder of each type (the generated object's own, `describe/phenomena.py`)."""

_DEFAULTS = {"reflex_offset_x": 0.0, "reflex_offset_y": 0.0, "reflex_offset_z": 0.0, "has_moons": False,
             "has_liquid_water": False, "ice_shell_thickness_km": 0.0, "ocean_depth_km": 0.0, "hosts": None,
             "q_total": None}
"""dict: What a builder reads that a table may leave out."""


class _Row:
    """A stored row with the attribute names the generated object has; a name the row lacks is `None`."""

    def __init__(self, values):
        self.__dict__.update(values)

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        return None


def _speed(detail):
    parts = [detail.get(f"velocity_{axis}_kms") for axis in "xyz"]
    return None if any(part is None for part in parts) else math.sqrt(sum(part * part for part in parts))


def _rogue_kind(detail):
    if detail.get("mass_bin") == "brown-dwarf":
        return "Brown Dwarf"
    return "Gas Giant" if detail.get("planet_type") == "g" else "Terrestrial"


def _fields(phenomenon_type, detail):
    """The values to add to `detail` for the builder of `phenomenon_type`."""
    extra = {"galactic_orbital_speed_kms": _speed(detail)}
    if phenomenon_type == "quasar":
        extra["galaxy_luminosity_multiple"] = (detail.get("luminosity_w") or 0.0) / tuning.MILKY_WAY_STELLAR_LUMINOSITY_W
    elif phenomenon_type == "rogue_planet":
        extra["kind_label"] = _rogue_kind(detail)
        extra["spectral_type"] = (brownDwarf.spectral_type(detail["effective_temperature_k"])
                                  if detail.get("mass_bin") == "brown-dwarf" and detail.get("effective_temperature_k")
                                  else None)
    elif phenomenon_type == "interstellar_comet":
        summary = detail.get("composition_summary") or ""
        extra["get_composition_summary"] = lambda: summary
    elif phenomenon_type == "asteroid_field":
        extra["composition"] = [(part["component"], part["concentration"]) for part in detail.get("composition") or []]
    elif phenomenon_type == "black_hole":
        extra.update(type=f"{(detail.get('mass_class') or 'stellar').replace('-', ' ').title()} Black Hole",
                     mass=_msun_to_kg(detail.get("mass_solar") or 0.0), radius=detail.get("event_horizon_radius_km"),
                     age=detail.get("age_gy"), temperature=detail.get("temperature_k"), name=detail["name"])
    elif phenomenon_type == "neutron_star":
        magnetar = detail.get("pulsar_type") == "magnetar"
        extra.update(type="Neutron Star (Magnetar)" if magnetar else "Neutron Star",
                     mass=_msun_to_kg(detail.get("mass_solar") or 0.0), radius=detail.get("radius_km"),
                     age=detail.get("age_gy"), name=detail["name"])
    elif phenomenon_type == "supernova_remnant":
        extra["compact_remnant"] = None    # the remnant's centre is its own page; its name is `compact_remnant_name`
    return extra


def stored_phenomenon_document(phenomenon_type, detail):
    """
    The Document of the stored phenomenon `detail` (`db.query.phenomenon_detail`) of `phenomenon_type`, one section with
    its facts table and sentences; `None` for a type with no paragraph, or a row too bare to describe (a column the words
    need is empty, as in a row stored before that column existed).
    """
    builder = _BUILDERS.get(phenomenon_type)
    if builder is None:
        return None
    row = _Row({**_DEFAULTS, **detail, **_fields(phenomenon_type, detail)})
    try:
        return model.document("phenomenon", None, detail["name"], [builder(row)])
    except (KeyError, TypeError, ValueError):
        return None
