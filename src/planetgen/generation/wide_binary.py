# planetgen/generation/wide_binary.py

"""
Wide (S-type) Binary Pair
==========================

This module defines `WideBinaryPair`, which represents the orbital
relationship between the two stars of an **S-type (wide) binary** --
two stars far enough apart that each keeps its own separate identity (its
own mass, luminosity, habitable zone) and can independently host its own
planets, with each star's maximum stable planetary orbit limited by the
other star's gravitational perturbation.

This is deliberately NOT a variant of `binary.BinaryStarProxy` and does
NOT subclass `Star`. `BinaryStarProxy` exists specifically to merge two
stars into one *effective* star (combined mass/luminosity/habitable zone)
for P-type (circumbinary) planet placement -- the opposite model from an
S-type pair, where the two stars never merge. Putting a non-merging,
non-`Star`-shaped orbital-elements class inside `generation/binary.py` would
misrepresent the physics that module's own class hierarchy already commits
to, so this lives in its own sibling module instead, the same way
`generation/star.py`/`generation/planet.py`/`asteroidData.py` each own one concept.

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

from planetgen import tuning
from planetgen.generation import companion_classes, place_terms
from planetgen.physics.planets import calculate_orbital_period_years
from planetgen.physics.system_tree import node_relative_au
from planetgen.util import draw
from planetgen.util.serialization import fields_from_dict, fields_to_dict
from planetgen.generation.star import Star
from planetgen.physics.orbits import (
    circular_orbital_speed_kms, hill_radius_in_pair_au, holman_wiegert_critical_semimajor_axis,
    minimum_update_interval_years, orbital_position_au,
)
from planetgen.physics import unit_registry as units

_km_to_au = units.converter("km", "au")


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
        min_update_interval_years (float): See `orbits.minimum_update_interval_years`.
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

    `primary_position_*_au`/`secondary_position_*_au` (schema v20) are each
    star's own offset from the pair's barycenter -- the same proper
    two-body treatment `binary.BinaryStarProxy` gets, alongside the
    pre-existing `position_x/y/z_au` (still the secondary's position
    relative to the primary, unchanged). `secondary_mass_fraction` is the
    constant `secondary.mass / (primary.mass + secondary.mass)`, stored so
    `store.advance_orbital_phases` never needs to join back to `stars`.
    """

    def __init__(self, system_config, primary: Star, secondary: Star, inner_boost=1.0, orbit=None):
        """
        Generates a fresh wide-binary orbital relationship between two
        already-generated stars, including each star's own
        Holman & Wiegert (1999) planetary stability limit.

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
            inner_boost (float): The weight of the 0.25 to 10 AU range of separations after the metallicity trend
                (`companion_classes.CompanionOdds.inner_boost`); 1 with none.
            orbit (OrbitNode, optional): The pair's stored orbit (a fill that builds a system from the tree its scatter
                row stored, GEN.315): its semi-major axis, eccentricity, plane, phase and period are the pair's and
                nothing is drawn for them. `None` draws them.
        """
        self.system_config = system_config

        if secondary.mass > primary.mass:
            primary, secondary = secondary, primary
        self.primary = primary
        self.secondary = secondary

        # The separation follows Moe & Di Stefano's companions per decade of
        # period for the primary's mass (GEN.448), and the eccentricities -- unlike the
        # close/P-type pair, which tidally circularizes -- are broadly
        # consistent with a "thermal" distribution; see
        # tuning.WIDE_BINARY_SEPARATION_MIN_AU and
        # WIDE_BINARY_ECCENTRICITY_MAX for the full justification. The widest
        # axis is the place's survival limit (GEN.103), not a constant.
        self.max_separation_au = wide_pair_limit_au(system_config, primary, secondary)
        """float: The widest axis the pair's place allowed when it was drawn (not saved)."""
        if orbit is None:
            self.separation_au = sample_wide_binary_separation_au(
                self.max_separation_au, place_terms.born_mass_solar(primary),
                place_terms.born_mass_solar(primary) + place_terms.born_mass_solar(secondary), inner_boost)
            self.eccentricity = sample_wide_binary_eccentricity()
        else:
            self.separation_au = _km_to_au(orbit.separation_km)
            self.eccentricity = orbit.eccentricity

        self.periapsis_au = self.separation_au * (1 - self.eccentricity)
        self.apoapsis_au = self.separation_au * (1 + self.eccentricity)

        total_mass_kg = primary.mass + secondary.mass
        self.period_years = (calculate_orbital_period_years(self.separation_au, total_mass_kg) if orbit is None
                             else orbit.period_years)
        # Real orbital speed varies with eccentric anomaly (vis-viva), not
        # constant the way a circular orbit's is -- but this generator's
        # engine only ever tracks circular-orbit positions/speeds for every
        # other body (planets, moons, and the existing P-type pair's own
        # mutual orbit), including in the database's raw-SQL phase-advance
        # mirror (`store.advance_orbital_phases`), which can't express a
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
        # (binary.BinaryStarProxy.__init__).
        if orbit is None:
            self.inclination_deg = draw.uniform(0, 180)
            self.ascending_node_deg = draw.uniform(0, 360)
            self.phase_deg = draw.uniform(0, 360)
        else:
            self.inclination_deg = orbit.inclination_deg
            self.ascending_node_deg = orbit.ascending_node_deg
            self.phase_deg = orbit.mean_anomaly_deg
        if orbit is None:
            self.position_x_au, self.position_y_au, self.position_z_au = orbital_position_au(
                self.separation_au, self.inclination_deg, self.ascending_node_deg, self.phase_deg
            )
        else:
            # A pair built from a tree sits where its node puts the secondary (the real ellipse, GEN.709), not on the
            # circle the phase is tracked by.
            self.position_x_au, self.position_y_au, self.position_z_au = node_relative_au(orbit, orbit.epoch_unix)

        # Each star's own offset from the pair's barycenter (AU) -- see
        # binary.BinaryStarProxy.__init__'s identical comment for why
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
        # limits unless the two masses are equal. It is the tighter of the
        # Holman & Wiegert fit and the star's own Hill sphere in the pair at
        # periapsis (GEN.228): the fit is the usual limit, the Hill sphere
        # the hard one nothing can stay bound beyond. Stored directly on
        # each `Star` (see `Star.a_crit_au`), not duplicated here.
        primary.a_crit_au = min(
            holman_wiegert_critical_semimajor_axis(
                self.separation_au, secondary.mass / total_mass_kg, self.eccentricity),
            hill_radius_in_pair_au(self.separation_au, primary.mass / total_mass_kg, self.eccentricity))
        secondary.a_crit_au = min(
            holman_wiegert_critical_semimajor_axis(
                self.separation_au, primary.mass / total_mass_kg, self.eccentricity),
            hill_radius_in_pair_au(self.separation_au, secondary.mass / total_mass_kg, self.eccentricity))

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
        # Same "heavier star is the pair's primary" convention `__init__`
        # applies -- a `StarSystem`'s own `primary_star`/`secondary_star`
        # (and so the `stars` rows' `role`) keep generation order, which
        # `__init__` may have swapped, so the pair's orientation has to be
        # re-derived from the masses rather than taken from the caller.
        if secondary.mass > primary.mass:
            primary, secondary = secondary, primary
        pair = object.__new__(cls)
        pair.system_config = system_config
        pair.primary = primary
        pair.secondary = secondary
        fields_from_dict(pair, data, cls.SERIALIZABLE_FIELDS)
        return pair

def wide_pair_limit_au(system_config, primary, secondary):
    """
    The widest axis a wide pair of these two stars keeps in its place (GEN.103,
    `place_terms.wide_axis_limit_au`): the place is `system_config.PLACE_TERMS`, or the Sun's neighbourhood for a
    system with none; the masses are the stars' at birth and the age the primary's.
    """
    terms = getattr(system_config, "PLACE_TERMS", None) or place_terms.SOLAR_NEIGHBOURHOOD
    total_mass_solar = place_terms.born_mass_solar(primary) + place_terms.born_mass_solar(secondary)
    limit = place_terms.wide_axis_limit_au(terms, total_mass_solar, primary.age)
    return min(limit, tuning.WIDE_BINARY_SEPARATION_CEILING_AU)


def sample_wide_binary_separation_au(max_au, born_mass_solar=1.0, total_mass_solar=2.0, inner_boost=1.0):
    """
    Draws an S-type (wide) binary's separation (semi-major axis) from Moe and Di Stefano's companions per decade of
    period for a primary of this mass (`generation.companion_classes`, GEN.448), between
    `tuning.WIDE_BINARY_SEPARATION_MIN_AU` and `max_au` (the place's survival limit, `wide_pair_limit_au`; a limit under
    the lower bound gives the lower bound, which only a pair forced wide reaches). One uniform draw.

    Args:
        max_au (float): The widest axis the pair keeps, AU.
        born_mass_solar (float): The primary's birth mass, solar masses (it picks the distribution).
        total_mass_solar (float): The pair's birth mass in total (it sets the period of an axis).
        inner_boost (float): The weight of the 0.25 to 10 AU range after the metallicity trend (1 with none).

    Returns:
        float: A separation, in AU.
    """
    return companion_classes.sample_separation_au(born_mass_solar, total_mass_solar, max_au, inner_boost)


def sample_wide_binary_eccentricity():
    """
    Draws an S-type (wide) binary's orbital eccentricity from a "thermal"
    distribution, `f(e) = 2e`, capped at
    `tuning.WIDE_BINARY_ECCENTRICITY_MAX` -- see that constant's
    own docstring for why wide pairs (unlike the close/P-type pair) keep a
    realistic, generally non-zero eccentricity.

    Derivation: the thermal PDF `f(e) = 2e` restricted to `[0, e_max]` and
    renormalized is still exactly proportional to `e` (just rescaled), so
    its CDF is `F(e) = (e / e_max)^2` and the closed-form inverse-CDF
    sample is `e = e_max * sqrt(u)`, `u ~ Uniform(0, 1)` -- no rejection
    sampling needed.

    Returns:
        float: An eccentricity, in [0, `WIDE_BINARY_ECCENTRICITY_MAX`).
    """
    return tuning.WIDE_BINARY_ECCENTRICITY_MAX * math.sqrt(draw.random())
