# planetgen/physics/bulk_composition.py

"""
The bulk composition of a planet or moon (GEN.428, docs/design/composition/bulk-composition.md, the contracts of
docs/design/composition/composition-spine-contracts.md)
=============================================================================================================

What a body is made of in bulk, as mass fractions of the whole body, read through one accessor so no
lane reads the columns that hold it. The draw is `physics/bulk_draw.py` (GEN.428)
and the derived radius, density and gravity are `physics/bulk_structure.py` (GEN.429). The columns (DB.69) are on `planets`
and `moons` from schema 144; a body saved before has none of them and `bulk_fractions` returns `None` for it.

The five fractions add up to 1 (to `SUM_TOLERANCE`):

- `metal`: the iron-rich core, light elements included.
- `rock`: the silicate mantle and crust (silicon carbide and graphite when `carbon_rock`).
- `water`: H2O in every phase (the body's `water_mass_fraction`).
- `ice`: every other ice (CO2, NH3, CH4, CO, N2, H2S, CH3OH, organics) as one share.
- `envelope`: the hydrogen and helium gas, with `envelope_z` of heavy elements inside it.

Use:

    bulk = bulk_composition.bulk_fractions(planet)       # BulkFractions, or None for a body saved before schema 144
    if bulk is not None:
        dry_core_fraction = bulk.metal / (bulk.metal + bulk.rock)
"""

from dataclasses import dataclass

SUM_TOLERANCE = 1e-9
"""float: How far the five fractions may be from 1."""

ORIGINS = ("drawn", "inverse")
"""tuple: `drawn` is route 2 (composition first), `inverse` route 1 (worked back from a given mass and radius)."""

FRACTION_ATTRIBUTES = {
    "metal": "metal_mass_fraction",
    "rock": "rock_mass_fraction",
    "water": "water_mass_fraction",
    "ice": "ice_mass_fraction",
    "envelope": "envelope_mass_fraction",
}
"""dict: Each whole-body fraction and the body attribute (and column) that holds it."""

DETAIL_ATTRIBUTES = {
    "rock_feo": "rock_feo_fraction",
    "rock_al2o3": "rock_al2o3_fraction",
    "rock_mg_si": "rock_mg_si_ratio",
    "envelope_helium": "envelope_helium_fraction",
    "envelope_z": "envelope_z_fraction",
    "porosity": "porosity",
}
"""dict: The optional descriptors and their attributes. `None` where the body has none (a rock fraction of an
envelope-only giant, a porosity of a planet)."""

FLAG_ATTRIBUTES = {"carbon_rock": "carbon_rock", "origin": "composition_origin"}

BODY_FIELDS = (*FRACTION_ATTRIBUTES.values(), *DETAIL_ATTRIBUTES.values(), *FLAG_ATTRIBUTES.values())
"""tuple: Every attribute (and column) a body holds for its bulk composition."""


NEW_FIELDS = tuple(name for name in BODY_FIELDS if name != "water_mass_fraction")
"""tuple: The attributes DB.69 added to `planets` and `moons`; `water_mass_fraction` was already a column (GEN.88) and `hydrosphere`
keeps it."""

COLUMN_FIELDS = tuple(name for name in NEW_FIELDS if name != "carbon_rock")
"""tuple: `NEW_FIELDS` that store as they are; `carbon_rock` is a flag stored as a tristate."""


class CompositionError(ValueError):
    """Fractions that do not describe a body."""


@dataclass(frozen=True)
class BulkFractions:
    """
    A body's bulk composition. `metal + rock + water + ice + envelope` is 1.

    Attributes:
        metal, rock, water, ice, envelope (float): Shares of the whole body's mass (see the module text).
        rock_feo (float or None): FeO, a mass fraction of the rock.
        rock_al2o3 (float or None): Al2O3, a mass fraction of the rock (the Earth's bulk silicate has 0.045).
        rock_mg_si (float or None): Molar Mg/Si of the rock.
        envelope_helium (float or None): Y, a mass fraction of the envelope.
        envelope_z (float or None): Heavy elements, a mass fraction of the envelope.
        porosity (float or None): Void share of a small body's volume.
        carbon_rock (bool): The rock is silicon carbide and graphite, not silicate.
        origin (str): One of `ORIGINS`.
    """

    metal: float
    rock: float
    water: float
    ice: float
    envelope: float
    rock_feo: float | None = None
    rock_al2o3: float | None = None
    rock_mg_si: float | None = None
    envelope_helium: float | None = None
    envelope_z: float | None = None
    porosity: float | None = None
    carbon_rock: bool = False
    origin: str = "drawn"

    def __post_init__(self):
        shares = self.shares()
        if any(value < 0.0 for value in shares.values()):
            raise CompositionError(f"a bulk fraction is negative: {shares}")
        if abs(sum(shares.values()) - 1.0) > SUM_TOLERANCE:
            raise CompositionError(f"the bulk fractions add up to {sum(shares.values())!r}, not 1: {shares}")
        if self.origin not in ORIGINS:
            raise CompositionError(f"origin {self.origin!r} is not one of {ORIGINS}")

    def shares(self):
        """The five whole-body fractions as a dict, in the order of `FRACTION_ATTRIBUTES`."""
        return {name: getattr(self, name) for name in FRACTION_ATTRIBUTES}

    @property
    def dry_core_fraction(self):
        """The metal's share of the dry part (metal plus rock): the CMF of the rocky-planet literature. `None` for
        a body with neither."""
        dry = self.metal + self.rock
        return self.metal / dry if dry > 0.0 else None

    def to_attributes(self):
        """The body attributes (and column values) that hold this composition, keyed by `BODY_FIELDS`."""
        values = {}
        for name, attribute in (*FRACTION_ATTRIBUTES.items(), *DETAIL_ATTRIBUTES.items(), *FLAG_ATTRIBUTES.items()):
            values[attribute] = getattr(self, name)
        return values


def bulk_fractions(body):
    """
    `body`'s bulk composition.

    Returns:
        BulkFractions or None: `None` for a body with no stored composition (saved before GEN.428, or an object
        that has none). A caller then uses its own stand-in and says so.

    Raises:
        CompositionError: The stored fractions do not add up to 1 or one is negative.
    """
    values = {}
    for name, attribute in FRACTION_ATTRIBUTES.items():
        value = getattr(body, attribute, None)
        if value is None:
            return None
        values[name] = float(value)
    for name, attribute in DETAIL_ATTRIBUTES.items():
        value = getattr(body, attribute, None)
        values[name] = None if value is None else float(value)
    values["carbon_rock"] = bool(getattr(body, "carbon_rock", False))
    values["origin"] = getattr(body, "composition_origin", None) or "drawn"
    return BulkFractions(**values)
