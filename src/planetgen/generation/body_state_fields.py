# planetgen/generation/body_state_fields.py

"""
Which stored fields of a star, planet or moon are its birth and which are its present state (GEN.582, light: the classification and its tests)
==============================================================================================================================================

The unified orbited-body plan (GEN.550, docs/design/architecture/unified-system-barycentre.md, Boss's card S10) splits what a body is from
what it has become: the **birth** values are the draws the body was made with (birth mass, age, [Fe/H], composition, formation distance and the
orbit-geometry draws); the **present** state is everything worked out from them about the primary the body now orbits, at the age it now has. A
remake (a body moved to a new primary or orbit, a star made again at another age) keeps the birth values and works the present state out again.

There is no Field table yet (GEN.565 moved to 8.2), so this module is a plain list that GEN.565 can absorb. It classifies; it reads and changes nothing, and no generated
value moves. A field is **present** unless it is named in `STAR_BIRTH` or `PLANET_BIRTH`, so a field added later is present until someone names it birth, and the tests
below fail on a name that no longer exists.

Use:

    is_birth(star, "age")                  # True
    present_fields(Planet)                 # the stored fields a remake must work out again
"""

from planetgen.generation.stellar_object import StellarObject
from planetgen.physics import bulk_composition, core_composition, giant_interior, interior_heat, moon_composition, radiogenic

STAR_BIRTH = (
    "name", "mass", "age", "metallicity_fe_h", "lifespan", "nuclear_arrival_age_gy",
    "spin_axis_x", "spin_axis_y", "spin_axis_z", "axial_tilt_deg",
    "galactic_orbital_speed_kms", "galactic_orbital_period_gy", "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
)
"""tuple: The stored Star fields that are birth values. `mass` is the mass the star is stored at (the initial mass, GEN.1611 slice 1 keeps it beside the present one when wind or mass loss moves it);
`lifespan` and `nuclear_arrival_age_gy` follow from the birth mass and [Fe/H]; the spin axis and tilt are draws of the star's making. Everything else of a star follows its age (the radius,
temperature, luminosity, class and type, the habitable zone and snow line that follow the luminosity, the rotation period that slows with age, the XUV and flare numbers, the disk), or its
companions (`a_crit_au`, the perimeter, the reflex offset), and is present."""

PLANET_BIRTH = (
    "name", "body_type", "is_moon", "mass", "composition", "flavor_key", "flavor_text_count",
    "born_beyond_snow_line", "collision_count",
    "distance", "eccentricity", "orbital_inclination_deg", "orbital_ascending_node_deg", "host_eccentricity",
    *radiogenic.BODY_FIELDS, *core_composition.BODY_FIELDS, *giant_interior.BODY_FIELDS, *moon_composition.BODY_FIELDS, *interior_heat.BODY_FIELDS,
    *(field for field in bulk_composition.NEW_FIELDS if not field.startswith("envelope_")),
)
"""tuple: The stored Planet (and Moon: a moon is a `Planet` with `is_moon`) fields that are birth values: the identity, the mass, the bulk composition the body was drawn with (rock oxides, core
light elements, ices, porosity; not the envelope shares, which the escape wears down), the radiogenic inventory, the giant's interior and core, the moon's disk and origin, the mantle
rheology, the draw of the formation (`born_beyond_snow_line`, `collision_count`) and the orbit-geometry draws (S10): semi-major axis, eccentricity, inclination and node, as first drawn. The
secular exchange and the tide move the stored orbit from them later; the draw record that keeps the first values is GEN.582's other half. The radius, density, gravity and volume are
present: an envelope world's radius follows its age. Everything not named here (the air and gases, the climate, the hydrosphere, the surface, the heat and tide outputs, the magnetism,
the light the body gets, the spin it has settled to, the rings, the class and the life) is worked out from the birth values about the primary and is present."""


def birth_fields(kind):
    """The birth fields of `kind`, `Star` or `Planet` (a moon is a Planet), as a tuple."""
    return STAR_BIRTH if kind is StellarObject or getattr(kind, "__name__", "") == "Star" else PLANET_BIRTH


def stored_fields(kind):
    """Every stored field of `kind`: its `SERIALIZABLE_FIELDS`."""
    return tuple(kind.SERIALIZABLE_FIELDS)


def present_fields(kind):
    """The stored fields of `kind` that are present state: what a remake must work out again."""
    birth = set(birth_fields(kind))
    return tuple(field for field in stored_fields(kind) if field not in birth)


def is_birth(kind, field):
    """Whether `field` of `kind` is a birth value."""
    return field in birth_fields(kind)
