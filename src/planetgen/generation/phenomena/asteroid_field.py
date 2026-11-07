# planetgen/generation/phenomena/asteroid_field.py

"""
Asteroid Field Generation
=========================

This module contains the `AsteroidField` class, a standalone descriptive
phenomenon representing a field of asteroid debris drifting in open
interstellar space, encountered independently of any particular star
system -- unlike `asteroidData.AsteroidBelt`, which always orbits a star
within a `StarSystem`. Physically these are the same kind of object
(density + mineral composition), just without a host star/orbit, so this
class draws its composition from the same taxa (`small_body_composition.draw_field`, GEN.477) and
writes it with `asteroidData.format_composition_summary`.

Like `nebulaData.Nebula` and the other standalone exotic phenomena, an
`AsteroidField` is generated on demand by `phenomenonGen.py`'s separate,
rarer exotic-phenomenon mode, never by `StarSystem._generate_planets`'s
normal per-slot rolls.
"""

import math

from planetgen.generation.belt import format_composition_summary
from planetgen.generation.config import SystemConfig
from planetgen import tuning
from planetgen.physics import constants, field_dispersal, small_body_composition
from planetgen.util import draw
from planetgen.util import log
from planetgen.generation.mass_object import GALACTIC_MOTION_UNITS
from planetgen.generation.small_body import SmallBody
from planetgen.physics import unit_registry as units

_ly_to_au = units.converter("ly", "au")


AU_PER_LY = _ly_to_au(1.0)
"""float: Astronomical units in one light-year, exact from the two definitions in `physics/constants.py` (GEN.572), for
the size digit."""


RING_RADIUS_AU = 50.0
"""float: The radius, AU, of the ring an `annulus` field was made as (the belt of the star that exploded, the design's 50 AU belt); a `ball` starts from a point. [Stand-in HF1, tier T6, tweak level L1.]"""


def radius_ly_at(age_yr, rms_speed_kms, shape="ball"):
    """The rms radius of a field `age_yr` after it was made, light-years: the field spreads from a point (a `ball`) or from a ring of `RING_RADIUS_AU` (an `annulus`) at the rms speed of its pieces
    (`field_dispersal.rms_radius_m`), so a field's radius is a result of its age (GEN.1103) and bring-to-epoch advances it by the same formula."""
    ring_m = RING_RADIUS_AU * constants.AU_M if shape == "annulus" else 0.0
    return field_dispersal.rms_radius_m(age_yr, rms_speed_kms * 1000.0, ring_m) / constants.LIGHTYEAR_M


def asteroid_field_size_digit(radius_ly):
    """The size digit of a field's class: floor(log10(radius in AU)),
    clamped to `ASTEROID_FIELD_SIZE_DIGIT_RANGE` (1-4)."""
    low, high = tuning.ASTEROID_FIELD_SIZE_DIGIT_RANGE
    return max(low, min(high, math.floor(math.log10(radius_ly * AU_PER_LY))))


def asteroid_field_class(composition_family, density, radius_ly):
    """
    A field's full class, letter plus size digit (e.g. `"C3"`).

    Args:
        composition_family (str): An `ASTEROID_FIELD_COMPOSITIONS` key.
        density (str): `"sparse"`, `"typical"` or `"dense"`.
        radius_ly (float): The field's radius, light-years.
    """
    letter = tuning.ASTEROID_FIELD_COMPOSITIONS[composition_family]["letters"][density]
    return f"{letter}{asteroid_field_size_digit(radius_ly)}"


def composition_family_for_letter(letter):
    """The composition family a class letter belongs to."""
    for family, data in tuning.ASTEROID_FIELD_COMPOSITIONS.items():
        if letter in data["letters"].values():
            return family
    raise ValueError(f"unknown asteroid field class letter {letter!r}")


def asteroid_field_designation(field_class, sector_code, index):
    """
    An asteroid field's designation (v40): `AF <class>-<sector>-<nn>`,
    e.g. `AF C3-FE81000A2B-01`, where `<sector>` is the sector it was
    found in (its grid designation, or its name) and `<nn>` counts that
    sector's fields from 01. `AF <class>-<nn>` for one found in no
    sector.
    """
    if sector_code:
        return f"AF {field_class}-{sector_code}-{index:02d}"
    return f"AF {field_class}-{index:02d}"


class AsteroidField(SmallBody):
    """
    A basic class to store information for a standalone field of asteroid
    debris drifting in open space, bound to no star.

    Attributes:
        name (str): A generated or explicitly given name for the field.
        field_class (str): The class, a letter from composition and
            density plus a size digit (e.g. `"C3"`;
            `tuning.ASTEROID_FIELD_COMPOSITIONS`).
        composition_family (str): The family the letter comes from
            (`"carbonaceous"`, `"stony"`, `"metallic"`, ...).
        density (str): The density of the field ('dense', 'sparse', 'typical'
            -- the same three levels `AsteroidBelt.density` uses).
        composition (list): A list of (component, concentration) tuples,
            generated the same way `AsteroidBelt.composition` is (the words of `makeup`).
        chondrite_group (str): The nearest meteorite group of the mean composition (`small_body_composition.chondrite_group`, GEN.475), `none` far from every group.
        makeup (SmallBodyMakeup or None): The drawn composition (GEN.477); None for a field read from a store without it.
        radius_ly (float): The field's rms radius, in light-years: a result of `birth_age_yr` and `rms_speed_kms` (`radius_ly_at`).
        rms_speed_kms (float): The rms speed of its pieces, km/s (how fast it spreads).
        birth_age_yr (float): Its age at the sector's epoch, years.
        shape (str): `ball` (or `annulus`), what the cloud was made as.
        origin (str or None): How it came to be (`unbound remnant` for a top-up field), None for one that claims none.
        history (str or None): The sentence of where it came from, None for a field with no history.
        largest_body_kg (float or None): The mass of its largest piece, None where it is not known.
        galactic_orbital_speed_kms (float): Circular orbital speed around
            the galactic center, km/s -- an asteroid field is still
            gravitationally part of the galaxy even though it isn't bound
            to any star (see `galactic_orbit.generate_galactic_orbit_fields`, and
            `nebulaData.Nebula`'s identical fields' docstring).
        galactic_orbital_period_gy (float): Orbital period, billions of years.
        galactic_orbital_phase_deg (float): Current angular position around
            that orbit -- advanced over time by `store.advance_orbital_phases`.
        galactic_min_update_interval_years (float): Floating-point update
            guard for `galactic_orbital_phase_deg`.
    """

    FIELD_UNITS = {
        "name": "text", "field_class": "text", "composition_family": "text", "density": "text", "radius_ly": "ly",
        "chondrite_group": "text", "rms_speed_kms": "km_s", "birth_age_yr": "yr", "shape": "text", "origin": "text", "history": "text", "largest_body_kg": "kg",
        **GALACTIC_MOTION_UNITS,
    }
    SERIALIZABLE_FIELDS = list(FIELD_UNITS)
    """Every attribute set by `__init__`, with its native unit (GEN.570), excluding `system_config` (a
    shared back-reference) and `composition` (handled separately in
    `to_dict`/`from_dict`, since it's a list of tuples -- JSON has no
    tuple type -- mirroring `AsteroidBelt.SERIALIZABLE_FIELDS`'s identical
    exclusion)."""

    def __init__(self, system_config: SystemConfig, name=None):
        """
        Initializes an AsteroidField object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object
                (its settings).
            name (str, optional): An explicit name. Random if omitted.
        """
        # The name is a placeholder: `store.insert_asteroid_field` replaces it with the
        # field's designation (`asteroid_field_designation`, v40).
        super().__init__(system_config, name=name)

        self.density = draw.choice(["dense", "sparse", "typical"])
        log.choice("Asteroid field density", self.density, "uniform draw among dense/sparse/typical")
        self.composition_family = small_body_composition.draw_field_family()
        log.choice("Asteroid field composition", self.composition_family,
                   "gated draw of the asteroid-field-family table (weights are the stand-in frequencies)")
        self.makeup = small_body_composition.draw_field(self.composition_family)
        self.composition = small_body_composition.describe(self.makeup)
        self.chondrite_group = self.makeup.detail["group"]
        self.rms_speed_kms = math.exp(draw.uniform(*(math.log(speed) for speed in tuning.ASTEROID_FIELD_SPEED_KMS_RANGE)))
        low_age = field_dispersal.age_yr_at_radius(tuning.ASTEROID_FIELD_RADIUS_RANGE_LY[0] * constants.LIGHTYEAR_M, self.rms_speed_kms * 1000.0)
        high_age = field_dispersal.retirement_age_yr(self.rms_speed_kms * 1000.0)
        self.birth_age_yr = draw.uniform(low_age, high_age)         # births are uniform in time: the age is uniform up to the retirement age, and the radius follows from it
        self.radius_ly = radius_ly_at(self.birth_age_yr, self.rms_speed_kms)
        log.choice("Asteroid field age", self.birth_age_yr, f"uniform between the age that gives the smallest field and the age at which it is retired, at {self.rms_speed_kms:.3g} km/s")
        self.shape, self.origin, self.history, self.largest_body_kg = "ball", None, None, None
        self.field_class = asteroid_field_class(self.composition_family, self.density, self.radius_ly)

        self.set_galactic_motion()

    def _encode(self, data):
        """`composition` as a list of 2-element lists (JSON has no tuple type), as `AsteroidBelt.to_dict` does, and the makeup."""
        data["composition"] = [list(pair) for pair in self.composition]
        data["makeup"] = self.makeup.to_dict() if self.makeup is not None else None

    def _decode(self, data):
        self.composition = [tuple(pair) for pair in data["composition"]]
        self.makeup = small_body_composition.SmallBodyMakeup.from_dict(data["makeup"]) if data.get("makeup") else None

    @classmethod
    def from_dict(cls, data, system_config):
        """Reconstructs an `AsteroidField` from a dict in the shape `to_dict()` produces, without re-running
        generation."""
        data = dict(data)
        if data.get("field_class") is None:
            # Saved before classes existed (schema v38): a random mineral
            # mix is honestly "mixed".
            data["composition_family"] = "mixed"
            data["field_class"] = asteroid_field_class("mixed", data["density"], data["radius_ly"])
        return cls._rebuild(data, system_config)

    def get_composition_summary(self):
        """
        Builds the human-readable composition summary phrase, the same
        as-published-text role `AsteroidBelt.get_composition_summary`
        plays for the database persistence layer.

        Returns:
            str: The composition summary phrase.
        """
        return format_composition_summary(self.composition)

from planetgen.generation import object_timing  # noqa: E402 -- registers the generators below (PERF.92)
object_timing.track(AsteroidField, "asteroid-field")
