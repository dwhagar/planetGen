# planetgen/generation/phenomena/quasar.py

"""
Quasar Generation
=================

This module contains the `Quasar` class: a galaxy's own central
supermassive black hole in its active, quasar phase -- gas spiralling in
through an accretion disk so hot it outshines every star in the galaxy
combined.

Unlike every other exotic phenomenon, a quasar is not something found
scattered through space: a galaxy has exactly one nucleus, at its
dynamical center. The generator therefore only ever places one there
(`generate.add_galactic_nucleus`, rolled against
`tuning.QUASAR_ACTIVE_NUCLEUS_CHANCE`), and a quasar has no
galactic orbit of its own -- it *is* the point everything else orbits.

Every derived quantity follows from two draws, the black hole's mass and
its Eddington ratio:
  - luminosity = Eddington ratio * the Eddington limit for that mass,
  - accretion rate = luminosity / (radiative efficiency * c^2),
  - broad-line-region radius from the reverberation-mapped
    radius-luminosity relation (see `program_constants`' `QUASAR_*`
    constants for sources).
"""


from planetgen.galaxy import shape_draw as shapeDraw
from planetgen.generation.config import SystemConfig
from planetgen.physics import constants as physical_constants, relations, unit_registry as units
from planetgen import tuning
from planetgen.util import draw
from planetgen.util import log
from planetgen.util.serialization import fields_from_dict, fields_to_dict
from planetgen.names.wordsalad import LazySaladName
from planetgen.util.random import log_uniform



class Quasar(LazySaladName):
    """
    A galaxy's active nucleus.

    Attributes:
        name (str): A generated or explicitly given name.
        black_hole_mass_solar (float): The central black hole's mass, in
            solar masses.
        event_horizon_radius_km (float): Its Schwarzschild radius, in km.
        eddington_ratio (float): Luminosity as a fraction of the
            Eddington limit.
        luminosity_w (float): Bolometric luminosity, in watts.
        accretion_rate_solar_per_year (float): Mass swallowed per year, in
            solar masses.
        broad_line_region_light_days (float): Radius of the fast-moving
            gas cloud around the disk that gives quasars their broad
            emission lines, in light-days.
        is_radio_loud (bool): Whether it launches relativistic jets.
        jet_length_ly (float or None): The jets' extent, in light-years
            (`None` when radio-quiet).
        active_age_years (float): How long this episode of activity has
            been running.
    """

    SERIALIZABLE_FIELDS = [
        "name", "black_hole_mass_solar", "event_horizon_radius_km", "eddington_ratio",
        "luminosity_w", "accretion_rate_solar_per_year", "broad_line_region_light_days",
        "is_radio_loud", "jet_length_ly", "active_age_years",
    ]
    """Every attribute set by `__init__`, excluding `system_config` (a
    shared back-reference)."""

    def __init__(self, system_config: SystemConfig, name=None, black_hole_mass_solar=None):
        """
        Initializes a Quasar object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object
                (its settings).
            name (str, optional): An explicit name. Random if omitted.
            black_hole_mass_solar (float, optional): The galaxy's central black hole (`galaxy.structure`, GEN.402);
                a quasar made on its own is the nucleus of a galaxy drawn for it (`shape_draw.random_nucleus_mass_solar`).
        """
        self.system_config = system_config
        self._init_salad_name(name)  # drawn on first read (PERF.43)
        self.name_given = bool(name)  # a given name is kept over an object ID (GEN.64)

        self.black_hole_mass_solar = black_hole_mass_solar or shapeDraw.random_nucleus_mass_solar()
        mass_kg = _msun_to_kg(self.black_hole_mass_solar)
        self.event_horizon_radius_km = relations.schwarzschild_radius(mass_kg) / 1000

        self.eddington_ratio = log_uniform(*tuning.QUASAR_EDDINGTON_RATIO_RANGE)
        self.luminosity_w = (
            self.eddington_ratio * tuning.EDDINGTON_LUMINOSITY_W_PER_SOLAR_MASS
            * self.black_hole_mass_solar
        )

        accretion_kg_s = self.luminosity_w / (
            tuning.QUASAR_RADIATIVE_EFFICIENCY * physical_constants.SPEED_OF_LIGHT_M_S ** 2
        )
        self.accretion_rate_solar_per_year = (
            _kg_to_msun(accretion_kg_s * physical_constants.SECONDS_PER_YEAR)
        )

        # lambda L_5100 in erg/s (1 W = 1e7 erg/s) for the Bentz et al.
        # radius-luminosity relation.
        l5100_erg_s = self.luminosity_w * 1e7 / tuning.QUASAR_BOLOMETRIC_CORRECTION_5100
        self.broad_line_region_light_days = (
            tuning.QUASAR_BLR_RADIUS_LIGHT_DAYS_AT_1E44
            * (l5100_erg_s / 1e44) ** tuning.QUASAR_BLR_RADIUS_LUMINOSITY_SLOPE
        )

        self.is_radio_loud = draw.random() < tuning.QUASAR_RADIO_LOUD_CHANCE
        log.choice("Quasar radio loudness", "radio-loud" if self.is_radio_loud else "radio-quiet",
                   f"roll against QUASAR_RADIO_LOUD_CHANCE ({tuning.QUASAR_RADIO_LOUD_CHANCE})")
        self.jet_length_ly = (
            log_uniform(*tuning.QUASAR_JET_LENGTH_RANGE_LY) if self.is_radio_loud else None
        )

        self.active_age_years = log_uniform(*tuning.QUASAR_ACTIVE_AGE_RANGE_YEARS)

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this quasar's properties, per
        `SERIALIZABLE_FIELDS`.

        Returns:
            dict: One entry per field in `SERIALIZABLE_FIELDS`.
        """
        return fields_to_dict(self, self.SERIALIZABLE_FIELDS)

    @classmethod
    def from_dict(cls, data, system_config):
        """
        Reconstructs a `Quasar` from a dict in the shape `to_dict()`
        produces, without re-running generation (`__init__` is bypassed
        via `object.__new__`).

        Args:
            data (dict): A dict in the shape `to_dict()` produces.
            system_config (SystemConfig): The shared config.

        Returns:
            Quasar: The reconstructed quasar.
        """
        quasar = object.__new__(cls)
        quasar.system_config = system_config
        fields_from_dict(quasar, data, cls.SERIALIZABLE_FIELDS)
        return quasar

    @property
    def galaxy_luminosity_multiple(self):
        """float: How many times brighter than a Milky-Way-sized galaxy's
        combined starlight this quasar shines."""
        return self.luminosity_w / tuning.MILKY_WAY_STELLAR_LUMINOSITY_W

from planetgen.generation import object_timing  # noqa: E402 -- registers the generators below (PERF.92)

_kg_to_msun = units.converter("kg", "Msun")
_msun_to_kg = units.converter("Msun", "kg")
object_timing.track(Quasar, "quasar")
