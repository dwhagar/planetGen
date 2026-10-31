# planetgen/physics/planet_type.py

"""
The type of a planet or moon: terrestrial, transitional or gas (GEN.1219; design: reports/class-type-and-domain-scoring-design-2026-10-31.md)
=========================================================================================================================================

The first qualification of a body (Boss, 10:31 AM EDT 2026-10-31), before any class: how its mass is shared between a core and an atmosphere.
The one number that cuts is the envelope's share of the mass, `envelope_mass_fraction` (the hydrogen-helium envelope of `bulk_structure`, Lopez
and Fortney 2014):

* **terrestrial**: under `TERRESTRIAL_MAX_ENVELOPE` (0.5 percent, about the 1.6 Earth radii of the rocky limit): mostly rock; any air is secondary and may be
  absent (Venus's 92 bar is 1e-4 of its mass);
* **transitional**: from there up to `GAS_MIN_ENVELOPE`: a real core under a real envelope (sub-Neptunes, gas dwarfs and, by Boss's choice of 10:36 AM EDT,
  the ice giants);
* **gas**: `GAS_MIN_ENVELOPE` (50 percent, the line `bulk_structure.GIANT_ENVELOPE_FRACTION` already holds) or more: the atmosphere outweighs the core.

A body whose envelope share is not stored is typed by its `body_type` (g gas, t terrestrial), and the type then says `known=False` in `type_of`.
Nothing stores or reads the type yet (GEN.1219 stores it with GEN.1207's epoch); `tests/planet_type_cli.py` prints how a generated sample splits.
"""

TERRESTRIAL = "terrestrial"
TRANSITIONAL = "transitional"
GAS = "gas"
TYPES = (TERRESTRIAL, TRANSITIONAL, GAS)

TERRESTRIAL_MAX_ENVELOPE = 0.005
"""float: STAND-IN PT1 (Boss's choice, L0): an envelope under this share of the mass is a terrestrial body (about 1.6 Earth radii; Rogers 2015, Fulton and Petigura 2018)."""
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
    """The type of a planet or moon (an object, or a column-named mapping)."""
    get = body.get if hasattr(body, "get") else lambda name: getattr(body, name, None)
    return type_of(get("envelope_mass_fraction"), get("body_type"))[0]
