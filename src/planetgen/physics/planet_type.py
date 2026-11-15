# planetgen/physics/planet_type.py

"""
The type of a planet or moon: terrestrial, transitional or gas (GEN.1219; design: reports/class-type-and-domain-scoring-design-2026-10-31.md)
=========================================================================================================================================

The first qualification of a body (Boss, 10:31 AM EDT 2026-10-31), before any class: how its mass is shared between a core and an atmosphere.
The one number that cuts is the envelope's share of the mass, `envelope_mass_fraction` (the hydrogen-helium envelope of `bulk_structure`, [R2741]):

* **terrestrial**: under `TERRESTRIAL_MAX_ENVELOPE` (0.5 percent, about the 1.6 Earth radii of the rocky limit): mostly rock; any air is secondary and may be
  absent (Venus's 92 bar is 1e-4 of its mass);
* **transitional**: from there up to `GAS_MIN_ENVELOPE`: a real core under a real envelope (sub-Neptunes, gas dwarfs and, by Boss's choice of 10:36 AM EDT,
  the ice giants);
* **super-puff** (GEN.1343, W15): a transitional planet that passes the gate of `super_puff.is_super_puff` (up to 30 Earth masses and under 0.1 g/cc, at any period and any age), the thick-envelope
  tail of the sub-Neptunes (Kepler-51 b, c and d);
* **gas**: `GAS_MIN_ENVELOPE` (50 percent, the line `bulk_structure.GIANT_ENVELOPE_FRACTION` already holds) or more: the atmosphere outweighs the core.

A body whose envelope share is not stored is typed by its `body_type` (g gas, t terrestrial), and the type then says `known=False` in `type_of`.
The type is stored on every planet and moon as `planet_type` (`fill`, called with the class, `planets.apply_class`, GEN.1219), after the envelope share is stored for
every body (`envelope_share`, GEN.1226), so a hydrogen-air dwarf is typed correctly. `tests/planet_type_cli.py` prints how a generated sample splits.
"""

from planetgen.physics import constants, super_puff

TERRESTRIAL = "terrestrial"
TRANSITIONAL = "transitional"
GAS = "gas"
SUPER_PUFF = "super-puff"
TYPES = (TERRESTRIAL, TRANSITIONAL, SUPER_PUFF, GAS)

TERRESTRIAL_MAX_ENVELOPE = 0.005
"""float: STAND-IN PT1 (Boss's choice, L0): an envelope under this share of the mass is a terrestrial body (about 1.6 Earth radii; [R1477, R3800])."""
GAS_MIN_ENVELOPE = 0.5
"""float: STAND-IN PT2 (Boss's choice, L0): an envelope of this share of the mass or more is a gas body (the line of `bulk_structure.GIANT_ENVELOPE_FRACTION`)."""


def type_of(envelope_fraction, body_type=None):
    """`(type, known)`: the type of a body from its envelope share of the mass, or from `body_type` ('t' or 'g') when the share is not stored."""
    if envelope_fraction is None:
        return (GAS if body_type == "g" else TERRESTRIAL), False
    if envelope_fraction < TERRESTRIAL_MAX_ENVELOPE:
        return TERRESTRIAL, True
    if envelope_fraction < GAS_MIN_ENVELOPE:
        return TRANSITIONAL, True
    return GAS, True


def planet_type(body):
    """The type of a planet or moon (an object, or a column-named mapping). A transitional planet of a star that passes the gate of `super_puff.is_super_puff`
    (GEN.1343, W15) is a super-puff: a planet of up to 30 Earth masses under 0.1 g/cc, at any period and any age (GEN.1470, GEN.1347). A column-named mapping has no star to judge by, so a transitional one is a super-puff when its stored `planet_type` says so."""
    stored = hasattr(body, "get")
    get = body.get if stored else lambda name: getattr(body, name, None)
    kind = type_of(get("envelope_mass_fraction"), get("body_type"))[0]
    if kind == TRANSITIONAL and stored and get("planet_type") == SUPER_PUFF:
        return SUPER_PUFF       # a stored row has no star or radius law to judge by, but it keeps the type its build gave it (GEN.1601)
    if kind == TRANSITIONAL and not stored and not getattr(body, "is_moon", False) and getattr(body, "star", None) is not None:
        if super_puff.is_super_puff(body.mass / constants.EARTH_MASS_TO_KG, body.radius / constants.EARTH_RADIUS_KM):
            return SUPER_PUFF
    return kind


def has_ground(body):
    """
    Whether a planet or moon has a surface to stand on (GEN.1334): a gas body has none, and neither has a rocky one under an envelope of `TERRESTRIAL_MAX_ENVELOPE` or more of
    its mass (a sub-Neptune, a gas dwarf), whose pressure rises through the hydrogen and helium down to the core; a body whose envelope share is not stored has its
    `body_type`'s.
    """
    get = body.get if hasattr(body, "get") else lambda name: getattr(body, name, None)
    return type_of(get("envelope_mass_fraction"), get("body_type"))[0] == TERRESTRIAL


SURFACE_PRESSURE_LABEL = "Surface pressure"
"""str: What a page calls the stored pressure of a body with a ground."""

REFERENCE_PRESSURE_LABEL = "Pressure at the 1 bar reference level"
"""str: What a page calls the stored pressure of a body with no ground (UX.297, GEN.1532): a giant, or a rocky planet under a deep envelope (a sub-Neptune, a gas dwarf), has no
surface to stand on, so its air is solved at the 1 bar level (`surface_air.REFERENCE_PRESSURE_BAR`, GEN.1334) and the stored pressure is that convention, not a surface."""


def has_ground_type(type_text):
    """Whether a stored `planet_type` text is a body with a ground (GEN.1334): terrestrial, or no type stored (`None`, an older row), which is read as having one."""
    return type_text in (None, TERRESTRIAL)


def pressure_label(type_text):
    """The label of a body's stored pressure from its stored `planet_type` text: the surface pressure for a ground, the reference-level pressure for none (UX.297)."""
    return SURFACE_PRESSURE_LABEL if has_ground_type(type_text) else REFERENCE_PRESSURE_LABEL


FIELDS = ("planet_type",)
"""tuple: What a planet or a moon stores of this module, as `planets` and `moons` columns (GEN.1219)."""


def fill(body):
    """Sets `planet_type` on a planet or a moon from its stored envelope share (`type_of`)."""
    body.planet_type = planet_type(body)
