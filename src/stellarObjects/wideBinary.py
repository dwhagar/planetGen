# stellarObjects/wideBinary.py

"""
Wide (S-type) Binary Pair
==========================

This module defines `WideBinaryPair`, which represents the orbital
relationship between the two stars of an **S-type (wide) binary** --
two stars far enough apart that each keeps its own separate identity (its
own mass, luminosity, habitable zone) and can independently host its own
planets, with each star's maximum stable planetary orbit limited by the
other star's gravitational perturbation.

This is deliberately NOT a variant of `doubleStar.BinaryStarProxy` and does
NOT subclass `Star`. `BinaryStarProxy` exists specifically to merge two
stars into one *effective* star (combined mass/luminosity/habitable zone)
for P-type (circumbinary) planet placement -- the opposite model from an
S-type pair, where the two stars never merge. Putting a non-merging,
non-`Star`-shaped orbital-elements class inside `doubleStar.py` would
misrepresent the physics that module's own class hierarchy already commits
to, so this lives in its own sibling module instead, the same way
`starData.py`/`planetData.py`/`asteroidData.py` each own one concept.

`WideBinaryPair` itself only carries the pair's own two-body orbital
elements (separation, eccentricity, inclination, ascending node, phase,
period, speed, periapsis/apoapsis, and a circular-approximation live
position -- see `__init__`'s docstring for why the live position stays
circular even though the orbit's eccentricity is real and used elsewhere).
It does not carry each star's own habitable zone, planets, or any other
per-star state -- `primary`/`secondary` are plain `Star` instances, stored
here only as a convenience back-reference (e.g. for rendering), never
merged or wrapped.
"""

import math
import random

from . import physical_constants
from .planetPhysics import calculate_orbital_period_years
from .serialization import fields_from_dict, fields_to_dict
from .starData import Star
from .utils import (circular_orbital_speed_kms, holman_wiegert_critical_semimajor_axis,
                    minimum_update_interval_years, orbital_position_au,
                    properties_to_string, sample_wide_binary_eccentricity,
                    sample_wide_binary_separation_au, to_scientific_notation,
                    years_to_time_string)


class WideBinaryPair:
    """
    The orbital relationship between the two stars of an S-type (wide)
    binary.

    Attributes:
        primary (Star): The more massive of the pair.
        secondary (Star): The less massive of the pair.
        separation_au (float): The pair's own orbital semi-major axis, in AU.
        eccentricity (float): The pair's own orbital eccentricity.
        period_years (float): The pair's own mutual orbital period, in years.
        speed_kms (float): The pair's own mutual orbital speed, in km/s
            (see `__init__`'s docstring -- a circular-orbit speed, since
            live position tracking stays circular).
        periapsis_au (float): Closest approach between the two stars,
            `separation_au * (1 - eccentricity)`.
        apoapsis_au (float): Farthest separation between the two stars,
            `separation_au * (1 + eccentricity)`.
        inclination_deg (float): The mutual orbit's plane tilt, in degrees.
        ascending_node_deg (float): The mutual orbit's longitude of the
            ascending node, in degrees.
        phase_deg (float): The secondary's current position within its
            orbit around the primary, in degrees.
        min_update_interval_years (float): See `utils.minimum_update_interval_years`.
        position_x_au/position_y_au/position_z_au (float): The secondary's
            current position relative to the primary, in AU.
    """

    SERIALIZABLE_FIELDS = [
        "separation_au", "eccentricity", "period_years", "speed_kms",
        "periapsis_au", "apoapsis_au", "inclination_deg",
        "ascending_node_deg", "phase_deg", "min_update_interval_years",
        "position_x_au", "position_y_au", "position_z_au",
        "primary_position_x_au", "primary_position_y_au", "primary_position_z_au",
        "secondary_position_x_au", "secondary_position_y_au", "secondary_position_z_au",
        "secondary_mass_fraction",
    ]
    """
    Orbital-elements state only -- deliberately excludes `primary`/
    `secondary` (the two `Star` instances serialize on their own, as
    `StarSystem.to_dict`'s own `star`/`secondary_star` keys, so nesting
    them again here would duplicate that data) and `a_crit_au` for either
    star (stored directly on each `Star` instance -- see
    `Star.a_crit_au`'s own docstring -- rather than duplicated here).

    `primary_position_*_au`/`secondary_position_*_au` (schema v18) are each
    star's own offset from the pair's barycenter -- the same proper
    two-body treatment `doubleStar.BinaryStarProxy` gets, alongside the
    pre-existing `position_x/y/z_au` (still the secondary's position
    relative to the primary, unchanged). `secondary_mass_fraction` is the
    constant `secondary.mass / (primary.mass + secondary.mass)`, stored so
    `_db.advance_orbital_phases` never needs to join back to `stars`.
    """

    def __init__(self, system_config, primary: Star, secondary: Star):
        """
        Generates a fresh wide-binary orbital relationship between two
        already-generated stars, including each star's own Holman &
        Wiegert (1999) planetary stability limit.

        Args:
            system_config (SystemConfig): The shared SystemConfig object
                for the system -- threaded through for display formatting
                only (`get_table_properties`/`to_paragraph_list`); nothing
                about generation itself depends on it.
            primary (Star): One of the pair's two stars.
            secondary (Star): The other. Swapped with `primary` if it
                turns out to be the more massive of the two, mirroring
                `BinaryStarProxy.__init__`'s identical convention -- the
                caller need not know which is heavier ahead of time.
        """
        self.system_config = system_config

        if secondary.mass > primary.mass:
            primary, secondary = secondary, primary
        self.primary = primary
        self.secondary = secondary

        # Real wide-binary separations are observed roughly log-uniform
        # over several decades, and their eccentricities -- unlike the
        # close/P-type pair, which tidally circularizes -- are broadly
        # consistent with a "thermal" distribution; see
        # program_constants.WIDE_BINARY_SEPARATION_MIN_AU/MAX_AU and
        # WIDE_BINARY_ECCENTRICITY_MAX for the full justification.
        self.separation_au = sample_wide_binary_separation_au()
        self.eccentricity = sample_wide_binary_eccentricity()

        self.periapsis_au = self.separation_au * (1 - self.eccentricity)
        self.apoapsis_au = self.separation_au * (1 + self.eccentricity)

        total_mass_kg = primary.mass + secondary.mass
        self.period_years = calculate_orbital_period_years(self.separation_au, total_mass_kg)
        # Real orbital speed varies with eccentric anomaly (vis-viva), not
        # constant the way a circular orbit's is -- but this generator's
        # engine only ever tracks circular-orbit positions/speeds for every
        # other body (planets, moons, and the existing P-type pair's own
        # mutual orbit), including in the database's raw-SQL phase-advance
        # mirror (`_db.advance_orbital_phases`), which can't express a
        # Kepler-equation solve as a plain SQL UPDATE. Eccentricity here is
        # real and used for `periapsis_au`/`apoapsis_au` above and for the
        # Holman & Wiegert stability limit below -- it just isn't threaded
        # into live position/speed tracking, a deliberate, documented
        # simplification rather than an oversight.
        self.speed_kms = circular_orbital_speed_kms(self.separation_au, self.period_years)
        self.min_update_interval_years = minimum_update_interval_years(self.period_years)

        # No protoplanetary-disk reason for the pair's own orbital plane to
        # align with anything -- same full-sphere convention already used
        # for the close/P-type pair's own mutual orbit
        # (doubleStar.BinaryStarProxy.__init__).
        self.inclination_deg = random.uniform(0, 180)
        self.ascending_node_deg = random.uniform(0, 360)
        self.phase_deg = random.uniform(0, 360)
        self.position_x_au, self.position_y_au, self.position_z_au = orbital_position_au(
            self.separation_au, self.inclination_deg, self.ascending_node_deg, self.phase_deg
        )

        # Each star's own offset from the pair's barycenter (AU) -- see
        # doubleStar.BinaryStarProxy.__init__'s identical comment for why
        # both stars get their own offset rather than one sitting fixed.
        self.secondary_mass_fraction = secondary.mass / total_mass_kg
        primary_mass_fraction = 1.0 - self.secondary_mass_fraction
        self.primary_position_x_au = -self.secondary_mass_fraction * self.position_x_au
        self.primary_position_y_au = -self.secondary_mass_fraction * self.position_y_au
        self.primary_position_z_au = -self.secondary_mass_fraction * self.position_z_au
        self.secondary_position_x_au = primary_mass_fraction * self.position_x_au
        self.secondary_position_y_au = primary_mass_fraction * self.position_y_au
        self.secondary_position_z_au = primary_mass_fraction * self.position_z_au

        # Each star's own critical semi-major axis is evaluated from THAT
        # star's perspective -- `mu` is always the *other* star's mass
        # fraction -- so the primary and secondary generally get different
        # limits unless the two masses are equal. Stored directly on each
        # `Star` (see `Star.a_crit_au`), not duplicated here.
        primary.a_crit_au = holman_wiegert_critical_semimajor_axis(
            self.separation_au, secondary.mass / total_mass_kg, self.eccentricity)
        secondary.a_crit_au = holman_wiegert_critical_semimajor_axis(
            self.separation_au, primary.mass / total_mass_kg, self.eccentricity)

    @property
    def primary_a_crit_au(self):
        """float: Convenience alias for `self.primary.a_crit_au`."""
        return self.primary.a_crit_au

    @property
    def secondary_a_crit_au(self):
        """float: Convenience alias for `self.secondary.a_crit_au`."""
        return self.secondary.a_crit_au

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this pair's orbital-elements
        state, per `SERIALIZABLE_FIELDS`. Does NOT include `primary`/
        `secondary` -- see `SERIALIZABLE_FIELDS`'s own docstring.

        Returns:
            dict: One entry per field in `SERIALIZABLE_FIELDS`.
        """
        return fields_to_dict(self, self.SERIALIZABLE_FIELDS)

    @classmethod
    def from_dict(cls, data, system_config, primary: Star, secondary: Star):
        """
        Reconstructs a `WideBinaryPair` from a dict in the shape
        `to_dict()` produces, without re-running generation (`__init__` is
        bypassed via `object.__new__`).

        Args:
            data (dict): A dict in the shape `to_dict()` produces.
            system_config (SystemConfig): The system's shared config.
            primary (Star): The pair's already-reconstructed primary star
                (via `Star.from_dict`) -- `StarSystem.from_dict` resolves
                this once and threads it in here, rather than this method
                reconstructing its own copy.
            secondary (Star): Ditto, for the secondary.

        Returns:
            WideBinaryPair: The reconstructed pair.
        """
        pair = object.__new__(cls)
        pair.system_config = system_config
        pair.primary = primary
        pair.secondary = secondary
        fields_from_dict(pair, data, cls.SERIALIZABLE_FIELDS)
        return pair

    def get_table_properties(self):
        """
        Builds the "Wide Binary System Data" property dict -- the exact
        key/value pairs `to_paragraph_list` renders into the pair's own
        relationship data table, distinct from either constituent star's
        own `Star.get_table_properties()` table.

        Returns:
            dict: Keys `separation`, `eccentricity`, `mutual_orbit`,
                 `primary_limit`, `secondary_limit`, each an
                 already-formatted display string.
        """
        separation_km = self.separation_au * physical_constants.AU_TO_KM
        separation_km_scientific = to_scientific_notation(self.system_config, separation_km)
        periapsis_km_scientific = to_scientific_notation(self.system_config, self.periapsis_au * physical_constants.AU_TO_KM)
        apoapsis_km_scientific = to_scientific_notation(self.system_config, self.apoapsis_au * physical_constants.AU_TO_KM)

        separation_string = (
            f"{separation_km_scientific} km ({self.separation_au:,.1f} AU), "
            f"ranging from {periapsis_km_scientific} km ({self.periapsis_au:,.1f} AU) at periapsis "
            f"to {apoapsis_km_scientific} km ({self.apoapsis_au:,.1f} AU) at apoapsis"
        )
        mutual_orbit_string = (
            f"{self.speed_kms:,.2f} km/s "
            f"({years_to_time_string(self.period_years)} per orbit)"
        )
        primary_offset_km = math.sqrt(
            self.primary_position_x_au ** 2 + self.primary_position_y_au ** 2 + self.primary_position_z_au ** 2
        ) * physical_constants.AU_TO_KM
        secondary_offset_km = math.sqrt(
            self.secondary_position_x_au ** 2 + self.secondary_position_y_au ** 2 + self.secondary_position_z_au ** 2
        ) * physical_constants.AU_TO_KM
        wobble_string = (
            f"{self.primary.name}: {to_scientific_notation(self.system_config, primary_offset_km)} km, "
            f"{self.secondary.name}: {to_scientific_notation(self.system_config, secondary_offset_km)} km "
            f"from the barycenter"
        )

        return {
            "separation": separation_string,
            "eccentricity": f"{self.eccentricity:.3f}",
            "mutual_orbit": mutual_orbit_string,
            "wobble": wobble_string,
            "primary_limit": f"{self.primary.name}'s planetary system is stable out to {self.primary_a_crit_au:.2f} AU due to {self.secondary.name}'s gravity.",
            "secondary_limit": f"{self.secondary.name}'s planetary system is stable out to {self.secondary_a_crit_au:.2f} AU due to {self.primary.name}'s gravity.",
        }

    def to_paragraph_list(self):
        """
        Generates a list of descriptive paragraphs for the wide binary
        pair's own orbital relationship. Does NOT include either
        constituent star's own individual details (those are each star's
        own `Star.to_paragraph_list()`, rendered separately -- see
        `systemData.StarSystem.__str__`'s wide-binary branch).

        Returns:
            list: A single-element list: [the pair's own data block].
        """
        pair_properties = self.get_table_properties()
        markdown_key_map = {
            "separation": "Stellar Separation", "eccentricity": "Orbital Eccentricity",
            "mutual_orbit": "Mutual Orbit", "wobble": "Barycenter Offset",
            "primary_limit": "Primary's Planetary Limit",
            "secondary_limit": "Secondary's Planetary Limit",
        }
        return [properties_to_string(self.system_config, pair_properties, "Wide Binary System Data", markdown_key_map=markdown_key_map)]
