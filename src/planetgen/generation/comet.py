# planetgen/generation/comet.py

"""
Star-Bound Comet Generation
===========================

This module contains the `Comet` class: a comet gravitationally bound to
a star, propagated via real two-body Kepler/Barker orbital mechanics (see
`kepler.py`) -- contrast `roguePlanetData.InterstellarComet`, an
unbound object passing through on a fixed hyperbolic trajectory,
encountered independently of any star system.

See docs/design/comet-orbital-realism.md for the full research/design
writeup this implements. In short, real comets split by orbital
eccentricity into three physically distinct categories; this module adds
the two that were missing (`InterstellarComet` already covers the third,
`e > 1`, hyperbolic case):

- "elliptical" (`0 <= eccentricity < 1`): a periodic, bound orbit that
  returns every `orbital_period_years` -- further tagged with a
  `period_class` (`tuning.COMET_PERIOD_CLASSES`) purely for
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

from planetgen.generation.config import SystemConfig
from planetgen.physics import constants as physical_constants, kepler
from planetgen.physics.position import HoldsOrbitPosition, axis_property, velocity_axis_property
from planetgen import tuning
from planetgen.util import draw
from planetgen.util import log
from planetgen.names.wordlists import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from planetgen.physics.planets import calculate_orbital_period_years
from planetgen.generation.phenomena.rogue import format_comet_composition_summary
from planetgen.util.serialization import fields_from_dict, fields_to_dict
from planetgen.names.wordsalad import generate_phoneme_salad_name
from planetgen.physics.orbits import minimum_update_interval_years
from planetgen.util.format import format_distance_au, format_period_years, format_speed_kms

PERIOD_CLASS_LABELS = {
    "jupiter_family": "Jupiter-family",
    "halley_type": "Halley-type",
    "long_period": "long-period",
}
"""dict: Human-readable label per `tuning.COMET_PERIOD_CLASSES`
key, for `Comet.to_paragraph_list`'s descriptive text -- presentation
only, not consulted by generation or propagation."""


PERIODIC_COMET_MAX_PERIOD_YEARS = 200
"""float: A bound comet with an orbital period under this many years is
periodic (`P/`); every other star-bound comet is `C/` -- the IAU's own
split (GEN.13, v40)."""


def comet_designation(host_name, index, comet):
    """
    A star-bound comet's designation (v40): `P/<host>-<n>` for a periodic
    comet (elliptical, period under `PERIODIC_COMET_MAX_PERIOD_YEARS`), or
    `C/<host>-<n>` otherwise, after the real IAU prefixes. `<host>` is the
    star it orbits (the system's name for a single star or close pair)
    and `<n>` counts that star's comets from 1.
    """
    periodic = (comet.orbit_type == "elliptical" and comet.orbital_period_years is not None
                and comet.orbital_period_years < PERIODIC_COMET_MAX_PERIOD_YEARS)
    return f"{'P' if periodic else 'C'}/{host_name}-{index}"


def rename_comet_designation(name, old_host, new_host):
    """
    `name` with its host swapped when it's a `comet_designation` whose
    host is `old_host` or starts with `old_host` followed by a space
    (`bodyNames.rename_prefix`'s rule), else `None`.
    """
    if name is None or len(name) < 3 or name[1] != "/" or name[0] not in "PC" or "-" not in name:
        return None
    host, _, number = name[2:].rpartition("-")
    if not number.isdigit():
        return None
    if host == old_host:
        renamed = new_host
    elif host.startswith(old_host + " "):
        renamed = new_host + host[len(old_host):]
    else:
        return None
    return f"{name[:2]}{renamed}-{number}"


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
    threshold = tuning.COMET_ACTIVITY_PERIHELION_THRESHOLD_AU
    max_chance = tuning.COMET_ACTIVITY_MAX_CHANCE
    min_chance = tuning.COMET_ACTIVITY_MIN_CHANCE
    fraction = min(1.0, max(0.0, perihelion_distance_au / threshold))
    return max_chance - fraction * (max_chance - min_chance)


class Comet(HoldsOrbitPosition):
    """
    A comet gravitationally bound to a star, on either a periodic
    elliptical orbit or a single-apparition parabolic one (see this
    module's own docstring).

    Attributes:
        name (str): A generated or explicitly given name.
        orbit_type (str): `"elliptical"` or `"parabolic"`.
        period_class (str or None): One of
            `tuning.COMET_PERIOD_CLASSES`'s keys for an
            elliptical comet (`"jupiter_family"`/`"halley_type"`/
            `"long_period"`); `None` for a parabolic comet (no period to
            classify).
        nucleus_diameter_km (float): Nucleus diameter, in kilometers.
        composition (list): A list of composition component strings,
            sampled from `tuning.COMET_COMPOSITION`.
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
            anomaly (see `kepler.parabolic_mean_anomaly`) -- `None`
            for an elliptical comet. Unlike `mean_anomaly_deg`, this does
            NOT wrap: a parabolic pass is a one-shot event, not periodic,
            and may be negative (still approaching perihelion).
        min_update_interval_years (float or None): The floating-point
            update guard `mean_anomaly_deg` needs (see
            `orbits.minimum_update_interval_years`) -- `None` for a
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
        velocity_x_kms / velocity_y_kms / velocity_z_kms (float): Current
            velocity relative to the star, km/s (`orbital_speed_kms` long).
    """

    SERIALIZABLE_FIELDS = [
        "name", "orbit_type", "period_class", "nucleus_diameter_km",
        "perihelion_distance_au", "eccentricity", "inclination_deg",
        "arg_periapsis_deg", "ascending_node_deg", "orbital_period_years",
        "mean_anomaly_deg", "parabolic_mean_anomaly", "min_update_interval_years",
        "primary_mass_solar", "is_active",
        "distance_au", "position_x_au", "position_y_au", "position_z_au", "orbital_speed_kms",
        "velocity_x_kms", "velocity_y_kms", "velocity_z_kms",
    ]
    """Every attribute set by `__init__`, excluding `system_config` (a
    shared back-reference) and `composition` (handled separately in
    `to_dict`/`from_dict`, mirroring `InterstellarComet.SERIALIZABLE_FIELDS`'s
    identical exclusion)."""

    position_x_au = axis_property(0)
    position_y_au = axis_property(1)
    position_z_au = axis_property(2)
    velocity_x_kms = velocity_axis_property(0)
    velocity_y_kms = velocity_axis_property(1)
    velocity_z_kms = velocity_axis_property(2)
    """float: This comet's offset from its star, AU: the "system" frame of
    `spatial`, its `SpatialPosition3D` (GEN.74), which `set_position_au` moves."""

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
                `tuning.COMET_PARABOLIC_CHANCE`.
        """
        self.system_config = system_config
        # A placeholder: `_db.insert_star_system` replaces it with the
        # comet's designation (`comet_designation`, v40).
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)
        self.name_given = bool(name)  # a given name is kept over an object ID (GEN.64)
        self.primary_mass_solar = primary_mass_solar

        if orbit_type:
            self.orbit_type = orbit_type
            log.debug(f"Comet orbit type: {orbit_type!r} (forced)")
        else:
            self.orbit_type = (
                "parabolic" if draw.random() < tuning.COMET_PARABOLIC_CHANCE else "elliptical"
            )
            log.choice("Comet orbit type", self.orbit_type,
                       f"roll against COMET_PARABOLIC_CHANCE ({tuning.COMET_PARABOLIC_CHANCE})")

        self.nucleus_diameter_km = draw.uniform(*tuning.BOUND_COMET_NUCLEUS_DIAMETER_RANGE_KM)
        num_components = min(3, len(tuning.COMET_COMPOSITION))
        self.composition = draw.sample(tuning.COMET_COMPOSITION, k=num_components)

        self.perihelion_distance_au = draw.uniform(*tuning.COMET_PERIHELION_DISTANCE_RANGE_AU)
        self.arg_periapsis_deg = draw.uniform(0, 360)
        self.ascending_node_deg = draw.uniform(0, 360)

        if self.orbit_type == "elliptical":
            class_names = list(tuning.COMET_PERIOD_CLASSES.keys())
            weights = [tuning.COMET_PERIOD_CLASSES[c]["weight"] for c in class_names]
            self.period_class = draw.choices(class_names, weights=weights, k=1)[0]
            log.choice("Comet period class", self.period_class,
                       f"weighted draw among {class_names} (weights {weights})")
            class_data = tuning.COMET_PERIOD_CLASSES[self.period_class]

            self.eccentricity = draw.uniform(*class_data["eccentricity_range"])
            self.inclination_deg = draw.uniform(0, class_data["inclination_max_deg"])

            semi_major_axis_au = self.perihelion_distance_au / (1 - self.eccentricity)
            primary_mass_kg = primary_mass_solar * physical_constants.SOLAR_MASS_TO_KG
            self.orbital_period_years = calculate_orbital_period_years(semi_major_axis_au, primary_mass_kg)

            self.mean_anomaly_deg = draw.uniform(0, 360)
            self.parabolic_mean_anomaly = None
            self.min_update_interval_years = minimum_update_interval_years(self.orbital_period_years)
        else:
            self.period_class = None
            self.eccentricity = draw.uniform(*tuning.PARABOLIC_COMET_ECCENTRICITY_RANGE)
            self.inclination_deg = draw.uniform(0, tuning.PARABOLIC_COMET_INCLINATION_MAX_DEG)

            self.orbital_period_years = None
            self.mean_anomaly_deg = None
            # A parabolic comet is "encountered" somewhere around its one
            # perihelion passage, not necessarily exactly at it -- drawn
            # symmetrically around zero (negative: still approaching;
            # positive: receding) so generated instances vary realistically
            # in current distance/activity rather than all starting frozen
            # at perihelion itself.
            self.parabolic_mean_anomaly = draw.uniform(-3.0, 3.0)
            self.min_update_interval_years = None

        self.update_orbital_state()
        activity_chance = _activity_chance(self.perihelion_distance_au)
        self.is_active = draw.random() < activity_chance
        log.choice("Comet activity", self.is_active,
                   f"roll against activity chance {activity_chance:.4g} at perihelion "
                   f"{self.perihelion_distance_au:.4g} AU")

    def orbit_mu_au3_per_year2(self):
        """AU^3/yr^2 of the comet's host star: 4 pi^2 per solar mass."""
        return 4.0 * math.pi ** 2 * self.primary_mass_solar

    def update_orbital_state(self):
        """
        Recomputes `distance_au`, `position_x/y/z_au`, and
        `orbital_speed_kms` from this comet's current orbital elements and
        anomaly (`mean_anomaly_deg` or `parabolic_mean_anomaly`, whichever
        applies to `orbit_type`), via `kepler.comet_orbital_state`.

        Called at generation time (right after the anomaly is rolled) and
        meant to be called again by any future periodic updater once
        `mean_anomaly_deg`/`parabolic_mean_anomaly` has been advanced by
        elapsed time -- the same "recompute anything derived from current
        orbital position" role `planetPhysics.update_orbital_position`
        plays for a planet/moon after its own `orbital_phase_deg` changes.
        """
        mean_anomaly_rad = math.radians(self.mean_anomaly_deg) if self.mean_anomaly_deg is not None else None
        state = kepler.comet_orbital_state(
            self.orbit_type, self.perihelion_distance_au, self.eccentricity,
            self.inclination_deg, self.arg_periapsis_deg, self.ascending_node_deg,
            self.primary_mass_solar,
            mean_anomaly_rad=mean_anomaly_rad,
            parabolic_mean_anomaly_value=self.parabolic_mean_anomaly,
            orbital_period_years=self.orbital_period_years,
        )
        self.distance_au = state["distance_au"]
        self.set_position_au(state["position_x_au"], state["position_y_au"], state["position_z_au"])
        self.orbital_speed_kms = state["orbital_speed_kms"]
        to_kms = physical_constants.AU_TO_KM / physical_constants.SECONDS_PER_YEAR
        self.set_velocity_kms(state["velocity_x_au_per_year"] * to_kms, state["velocity_y_au_per_year"] * to_kms,
                              state["velocity_z_au_per_year"] * to_kms)

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
                f"{format_distance_au(self.perihelion_distance_au)}, an eccentricity of {self.eccentricity:.3f}, and a "
                f"period of {format_period_years(self.orbital_period_years)}, returning to the inner system "
                f"every orbit."
            )
        else:
            orbit_sentence = (
                f"It is on a marginally unbound, near-parabolic orbit with a perihelion of "
                f"{format_distance_au(self.perihelion_distance_au)} and an eccentricity of {self.eccentricity:.4f} -- having "
                f"originated in this system's own outer reaches, it will not return after this passage."
            )

        description = (
            f"{self.name} is a comet bound to this system, with a nucleus roughly "
            f"{self.nucleus_diameter_km:.2f} km across, composed of {self.get_composition_summary()}. "
            f"{orbit_sentence} It is currently {format_distance_au(self.distance_au)} from the star, moving at "
            f"{format_speed_kms(self.orbital_speed_kms)}. {activity}"
        )

        return [header, description]

    def __str__(self):
        """
        Returns a string representation of the comet, with paragraphs
        separated by double newlines.
        """
        return "\n\n".join(self.to_paragraph_list())
