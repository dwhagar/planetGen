# stellarObjects/cometData.py

"""
Star-Bound Comet Generation
===========================

This module contains the `Comet` class: a comet gravitationally bound to
a star, propagated via real two-body Kepler/Barker orbital mechanics (see
`keplerMotion.py`) -- contrast `roguePlanetData.InterstellarComet`, an
unbound object passing through on a fixed hyperbolic trajectory,
encountered independently of any star system.

See docs/design/comet-orbital-realism.md for the full research/design
writeup this implements. In short, real comets split by orbital
eccentricity into three physically distinct categories; this module adds
the two that were missing (`InterstellarComet` already covers the third,
`e > 1`, hyperbolic case):

- "elliptical" (`0 <= eccentricity < 1`): a periodic, bound orbit that
  returns every `orbital_period_years` -- further tagged with a
  `period_class` (`program_constants.COMET_PERIOD_CLASSES`) purely for
  descriptive/plausibility flavor (Jupiter-family/Halley-type/long-period
  -- propagation itself is identical for every subtype).
- "parabolic" (`eccentricity` just under 1, always propagated as an exact
  parabola): a marginally-bound, single-apparition orbit that escapes
  after this one perihelion passage and never returns, without ever
  having been an interstellar visitor -- it originated in this system's
  own outer reaches (its Oort-Cloud analog).

This is a straightforward two-body model: `Comet` is generated given its
primary's mass (in solar masses) and rolls its own orbital elements
independently -- it does not (yet) model gravitational capture or a
planetary slingshot, both of which require a genuine three-body close
encounter (see docs/design/comet-orbital-realism.md's Context section for
why that's out of scope for this pass).
"""

import math
import random

from .config import SystemConfig
from . import keplerMotion, physical_constants, program_constants
from .names import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from .planetPhysics import calculate_orbital_period_years
from .roguePlanetData import format_comet_composition_summary
from .serialization import fields_from_dict, fields_to_dict
from .utils import generate_phoneme_salad_name, minimum_update_interval_years, reseed_rng, years_to_time_string

PERIOD_CLASS_LABELS = {
    "jupiter_family": "Jupiter-family",
    "halley_type": "Halley-type",
    "long_period": "long-period",
}
"""dict: Human-readable label per `program_constants.COMET_PERIOD_CLASSES`
key, for `Comet.to_paragraph_list`'s descriptive text -- presentation
only, not consulted by generation or propagation."""


def _activity_chance(perihelion_distance_au):
    """
    Linearly interpolates the chance a comet currently shows a coma/tail,
    from `program_constants.COMET_ACTIVITY_MAX_CHANCE` at
    `perihelion_distance_au == 0` down to `COMET_ACTIVITY_MIN_CHANCE` at
    or beyond `COMET_ACTIVITY_PERIHELION_THRESHOLD_AU` -- real cometary
    activity is driven by solar heating at the comet's *closest* approach,
    not by which orbit_type/period_class it happens to be (see those
    constants' own docstrings and docs/design/comet-orbital-realism.md).

    Args:
        perihelion_distance_au (float): Perihelion distance `q`, in AU.

    Returns:
        float: Activity chance, in `[COMET_ACTIVITY_MIN_CHANCE,
              COMET_ACTIVITY_MAX_CHANCE]`.
    """
    threshold = program_constants.COMET_ACTIVITY_PERIHELION_THRESHOLD_AU
    max_chance = program_constants.COMET_ACTIVITY_MAX_CHANCE
    min_chance = program_constants.COMET_ACTIVITY_MIN_CHANCE
    fraction = min(1.0, max(0.0, perihelion_distance_au / threshold))
    return max_chance - fraction * (max_chance - min_chance)


class Comet:
    """
    A comet gravitationally bound to a star, on either a periodic
    elliptical orbit or a single-apparition parabolic one (see this
    module's own docstring).

    Attributes:
        name (str): A generated or explicitly given name.
        orbit_type (str): `"elliptical"` or `"parabolic"`.
        period_class (str or None): One of
            `program_constants.COMET_PERIOD_CLASSES`'s keys for an
            elliptical comet (`"jupiter_family"`/`"halley_type"`/
            `"long_period"`); `None` for a parabolic comet (no period to
            classify).
        nucleus_diameter_km (float): Nucleus diameter, in kilometers.
        composition (list): A list of composition component strings,
            sampled from `program_constants.COMET_COMPOSITION`.
        perihelion_distance_au (float): Perihelion distance `q`, in AU.
        eccentricity (float): Orbital eccentricity.
        inclination_deg (float): Orbital plane tilt, in degrees.
        arg_periapsis_deg (float): Argument of periapsis, in degrees.
        ascending_node_deg (float): Longitude of the ascending node, in
            degrees.
        orbital_period_years (float or None): Orbital period, in years --
            `None` for a parabolic comet (infinite/undefined).
        mean_anomaly_deg (float or None): Current mean anomaly, in
            degrees -- the eccentric-orbit analog of
            `Planet.orbital_phase_deg`, advancing linearly with time
            (`0-360`, wraps). `None` for a parabolic comet.
        parabolic_mean_anomaly (float or None): Current parabolic mean
            anomaly (see `keplerMotion.parabolic_mean_anomaly`) -- `None`
            for an elliptical comet. Unlike `mean_anomaly_deg`, this does
            NOT wrap: a parabolic pass is a one-shot event, not periodic,
            and may be negative (still approaching perihelion).
        min_update_interval_years (float or None): The floating-point
            update guard `mean_anomaly_deg` needs (see
            `utils.minimum_update_interval_years`) -- `None` for a
            parabolic comet, which has no periodic angle to guard the
            same way (see `parabolic_mean_anomaly`'s own docstring).
        primary_mass_solar (float): The host star's mass, in solar
            masses -- stored directly (rather than a `Star` back-
            reference) so this class stays decoupled from the full `Star`
            object graph, the same way `AsteroidBelt` needs no star
            reference of its own.
        is_active (bool): Whether it currently shows a coma/tail from
            sublimating ices (see `_activity_chance`).
        distance_au (float): Current distance from the star, in AU.
        position_x_au / position_y_au / position_z_au (float): Current
            3D position relative to the star, in AU.
        orbital_speed_kms (float): Current orbital speed, in km/s.
    """

    SERIALIZABLE_FIELDS = [
        "name", "orbit_type", "period_class", "nucleus_diameter_km",
        "perihelion_distance_au", "eccentricity", "inclination_deg",
        "arg_periapsis_deg", "ascending_node_deg", "orbital_period_years",
        "mean_anomaly_deg", "parabolic_mean_anomaly", "min_update_interval_years",
        "primary_mass_solar", "is_active",
        "distance_au", "position_x_au", "position_y_au", "position_z_au", "orbital_speed_kms",
    ]
    """Every attribute set by `__init__`, excluding `system_config` (a
    shared back-reference) and `composition` (handled separately in
    `to_dict`/`from_dict`, mirroring `InterstellarComet.SERIALIZABLE_FIELDS`'s
    identical exclusion)."""

    def __init__(self, system_config: SystemConfig, primary_mass_solar, name=None, orbit_type=None):
        """
        Initializes a Comet object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object
                (used only for its `MARKDOWN` flag, via `to_paragraph_list`).
            primary_mass_solar (float): The host star's mass, in solar
                masses.
            name (str, optional): An explicit name. Random if omitted.
            orbit_type (str, optional): Force `"elliptical"` or
                `"parabolic"` rather than rolling
                `program_constants.COMET_PARABOLIC_CHANCE`.
        """
        reseed_rng()
        self.system_config = system_config
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)
        self.primary_mass_solar = primary_mass_solar

        self.orbit_type = orbit_type if orbit_type else (
            "parabolic" if random.random() < program_constants.COMET_PARABOLIC_CHANCE else "elliptical"
        )

        self.nucleus_diameter_km = random.uniform(*program_constants.BOUND_COMET_NUCLEUS_DIAMETER_RANGE_KM)
        num_components = min(3, len(program_constants.COMET_COMPOSITION))
        self.composition = random.sample(program_constants.COMET_COMPOSITION, k=num_components)

        self.perihelion_distance_au = random.uniform(*program_constants.COMET_PERIHELION_DISTANCE_RANGE_AU)
        self.arg_periapsis_deg = random.uniform(0, 360)
        self.ascending_node_deg = random.uniform(0, 360)

        if self.orbit_type == "elliptical":
            class_names = list(program_constants.COMET_PERIOD_CLASSES.keys())
            weights = [program_constants.COMET_PERIOD_CLASSES[c]["weight"] for c in class_names]
            self.period_class = random.choices(class_names, weights=weights, k=1)[0]
            class_data = program_constants.COMET_PERIOD_CLASSES[self.period_class]

            self.eccentricity = random.uniform(*class_data["eccentricity_range"])
            self.inclination_deg = random.uniform(0, class_data["inclination_max_deg"])

            semi_major_axis_au = self.perihelion_distance_au / (1 - self.eccentricity)
            primary_mass_kg = primary_mass_solar * physical_constants.SOLAR_MASS_TO_KG
            self.orbital_period_years = calculate_orbital_period_years(semi_major_axis_au, primary_mass_kg)

            self.mean_anomaly_deg = random.uniform(0, 360)
            self.parabolic_mean_anomaly = None
            self.min_update_interval_years = minimum_update_interval_years(self.orbital_period_years)
        else:
            self.period_class = None
            self.eccentricity = random.uniform(*program_constants.PARABOLIC_COMET_ECCENTRICITY_RANGE)
            self.inclination_deg = random.uniform(0, program_constants.PARABOLIC_COMET_INCLINATION_MAX_DEG)

            self.orbital_period_years = None
            self.mean_anomaly_deg = None
            # A parabolic comet is "encountered" somewhere around its one
            # perihelion passage, not necessarily exactly at it -- drawn
            # symmetrically around zero (negative: still approaching;
            # positive: receding) so generated instances vary realistically
            # in current distance/activity rather than all starting frozen
            # at perihelion itself.
            self.parabolic_mean_anomaly = random.uniform(-3.0, 3.0)
            self.min_update_interval_years = None

        self.update_orbital_state()
        self.is_active = random.random() < _activity_chance(self.perihelion_distance_au)

    def update_orbital_state(self):
        """
        Recomputes `distance_au`, `position_x/y/z_au`, and
        `orbital_speed_kms` from this comet's current orbital elements and
        anomaly (`mean_anomaly_deg` or `parabolic_mean_anomaly`, whichever
        applies to `orbit_type`), via `keplerMotion.comet_orbital_state`.

        Called at generation time (right after the anomaly is rolled) and
        meant to be called again by any future periodic updater once
        `mean_anomaly_deg`/`parabolic_mean_anomaly` has been advanced by
        elapsed time -- the same "recompute anything derived from current
        orbital position" role `planetPhysics.update_orbital_position`
        plays for a planet/moon after its own `orbital_phase_deg` changes.
        """
        mean_anomaly_rad = math.radians(self.mean_anomaly_deg) if self.mean_anomaly_deg is not None else None
        state = keplerMotion.comet_orbital_state(
            self.orbit_type, self.perihelion_distance_au, self.eccentricity,
            self.inclination_deg, self.arg_periapsis_deg, self.ascending_node_deg,
            self.primary_mass_solar,
            mean_anomaly_rad=mean_anomaly_rad,
            parabolic_mean_anomaly_value=self.parabolic_mean_anomaly,
            orbital_period_years=self.orbital_period_years,
        )
        self.distance_au = state["distance_au"]
        self.position_x_au = state["position_x_au"]
        self.position_y_au = state["position_y_au"]
        self.position_z_au = state["position_z_au"]
        self.orbital_speed_kms = state["orbital_speed_kms"]

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this comet's properties, per
        `SERIALIZABLE_FIELDS`, plus `composition`.

        Returns:
            dict: One entry per field in `SERIALIZABLE_FIELDS`, plus
                 `composition`.
        """
        data = fields_to_dict(self, self.SERIALIZABLE_FIELDS)
        data["composition"] = list(self.composition)
        return data

    @classmethod
    def from_dict(cls, data, system_config):
        """
        Reconstructs a `Comet` from a dict in the shape `to_dict()`
        produces, without re-running generation (`__init__` is bypassed
        via `object.__new__`).

        Args:
            data (dict): A dict in the shape `to_dict()` produces.
            system_config (SystemConfig): The shared config.

        Returns:
            Comet: The reconstructed comet.
        """
        comet = object.__new__(cls)
        comet.system_config = system_config
        fields_from_dict(comet, data, cls.SERIALIZABLE_FIELDS)
        comet.composition = list(data["composition"])
        return comet

    def get_composition_summary(self):
        """
        Builds the human-readable composition summary phrase, the same
        as-published-text role `AsteroidBelt.get_composition_summary`/
        `InterstellarComet.get_composition_summary` play elsewhere.

        Returns:
            str: The composition summary phrase.
        """
        return format_comet_composition_summary(self.composition)

    def to_paragraph_list(self):
        """
        Generates a list of descriptive paragraphs for the comet: a header
        naming it, and a description of its size, composition, orbit, and
        current position/activity.

        Returns:
            list: A list of strings, where each string is a paragraph
                  describing the comet.
        """
        header_level = '##' if self.system_config.MARKDOWN else '=='
        kind_label = "Parabolic Comet" if self.orbit_type == "parabolic" else "Comet"
        header = f"{header_level} {self.name} ({kind_label}) {header_level if not self.system_config.MARKDOWN else ''}".rstrip()

        activity = (
            "It currently shows an active coma and tail as its surface ices sublimate."
            if self.is_active else
            "It shows no current activity, its surface ices dormant for now."
        )

        if self.orbit_type == "elliptical":
            period_label = PERIOD_CLASS_LABELS.get(self.period_class, "periodic")
            orbit_sentence = (
                f"It follows a {period_label} elliptical orbit around the star, with a perihelion of "
                f"{self.perihelion_distance_au:.3f} AU, an eccentricity of {self.eccentricity:.3f}, and a "
                f"period of {years_to_time_string(self.orbital_period_years)}, returning to the inner system "
                f"every orbit."
            )
        else:
            orbit_sentence = (
                f"It is on a marginally unbound, near-parabolic orbit with a perihelion of "
                f"{self.perihelion_distance_au:.3f} AU and an eccentricity of {self.eccentricity:.4f} -- having "
                f"originated in this system's own outer reaches, it will not return after this passage."
            )

        description = (
            f"{self.name} is a comet bound to this system, with a nucleus roughly "
            f"{self.nucleus_diameter_km:.2f} km across, composed of {self.get_composition_summary()}. "
            f"{orbit_sentence} It is currently {self.distance_au:.3f} AU from the star, moving at "
            f"{self.orbital_speed_kms:.1f} km/s. {activity}"
        )

        return [header, description]

    def __str__(self):
        """
        Returns a string representation of the comet, with paragraphs
        separated by double newlines.
        """
        return "\n\n".join(self.to_paragraph_list())
