# stellarObjects/quasarData.py

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
`program_constants.QUASAR_ACTIVE_NUCLEUS_CHANCE`), and a quasar has no
galactic orbit of its own -- it *is* the point everything else orbits.

Every derived quantity follows from two draws, the black hole's mass and
its Eddington ratio:
  - luminosity = Eddington ratio * the Eddington limit for that mass,
  - accretion rate = luminosity / (radiative efficiency * c^2),
  - broad-line-region radius from the reverberation-mapped
    radius-luminosity relation (see `program_constants`' `QUASAR_*`
    constants for sources).
"""

import math
import random

from .config import SystemConfig
from .names import STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES
from . import log, physical_constants, program_constants
from .serialization import fields_from_dict, fields_to_dict
from .utils import generate_phoneme_salad_name, reseed_rng


def _log_uniform(low, high):
    """A value drawn uniformly in log space between `low` and `high`."""
    return math.exp(random.uniform(math.log(low), math.log(high)))


class Quasar:
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

    def __init__(self, system_config: SystemConfig, name=None):
        """
        Initializes a Quasar object.

        Args:
            system_config (SystemConfig): The shared SystemConfig object
                (used only for its `MARKDOWN` flag, via `to_paragraph_list`).
            name (str, optional): An explicit name. Random if omitted.
        """
        reseed_rng()
        self.system_config = system_config
        self.name = name if name else generate_phoneme_salad_name(STAR_NAMES, STAR_PREFIXES, STAR_SUFFIXES)

        self.black_hole_mass_solar = _log_uniform(*program_constants.QUASAR_BLACK_HOLE_MASS_RANGE_SOLAR)
        mass_kg = self.black_hole_mass_solar * physical_constants.SOLAR_MASS_TO_KG
        self.event_horizon_radius_km = (
            2 * physical_constants.G * mass_kg / physical_constants.SPEED_OF_LIGHT_M_S ** 2 / 1000
        )

        self.eddington_ratio = _log_uniform(*program_constants.QUASAR_EDDINGTON_RATIO_RANGE)
        self.luminosity_w = (
            self.eddington_ratio * program_constants.EDDINGTON_LUMINOSITY_W_PER_SOLAR_MASS
            * self.black_hole_mass_solar
        )

        accretion_kg_s = self.luminosity_w / (
            program_constants.QUASAR_RADIATIVE_EFFICIENCY * physical_constants.SPEED_OF_LIGHT_M_S ** 2
        )
        self.accretion_rate_solar_per_year = (
            accretion_kg_s * physical_constants.SECONDS_PER_YEAR / physical_constants.SOLAR_MASS_TO_KG
        )

        # lambda L_5100 in erg/s (1 W = 1e7 erg/s) for the Bentz et al.
        # radius-luminosity relation.
        l5100_erg_s = self.luminosity_w * 1e7 / program_constants.QUASAR_BOLOMETRIC_CORRECTION_5100
        self.broad_line_region_light_days = (
            program_constants.QUASAR_BLR_RADIUS_LIGHT_DAYS_AT_1E44
            * (l5100_erg_s / 1e44) ** program_constants.QUASAR_BLR_RADIUS_LUMINOSITY_SLOPE
        )

        self.is_radio_loud = random.random() < program_constants.QUASAR_RADIO_LOUD_CHANCE
        log.choice("Quasar radio loudness", "radio-loud" if self.is_radio_loud else "radio-quiet",
                   f"roll against QUASAR_RADIO_LOUD_CHANCE ({program_constants.QUASAR_RADIO_LOUD_CHANCE})")
        self.jet_length_ly = (
            _log_uniform(*program_constants.QUASAR_JET_LENGTH_RANGE_LY) if self.is_radio_loud else None
        )

        self.active_age_years = _log_uniform(*program_constants.QUASAR_ACTIVE_AGE_RANGE_YEARS)

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
        return self.luminosity_w / program_constants.MILKY_WAY_STELLAR_LUMINOSITY_W

    def to_paragraph_list(self):
        """
        Generates a list of descriptive paragraphs for the quasar: a
        header naming it, its engine (black hole, accretion, output), and
        its jets if it has any.

        Returns:
            list: A list of strings, where each string is a paragraph
                  describing the quasar.
        """
        header_level = '##' if self.system_config.MARKDOWN else '=='
        header = f"{header_level} {self.name} (Quasar) {header_level if not self.system_config.MARKDOWN else ''}".rstrip()

        description = (
            f"{self.name} is the galaxy's active nucleus: a supermassive black hole of "
            f"{self.black_hole_mass_solar:.2e} solar masses, its event horizon "
            f"{self.event_horizon_radius_km / physical_constants.AU_TO_KM:,.1f} AU in radius, swallowing "
            f"about {self.accretion_rate_solar_per_year:,.1f} solar masses of gas every year. Its accretion disk "
            f"radiates {self.luminosity_w:.2e} W ({self.eddington_ratio:.0%} of its Eddington limit), "
            f"roughly {self.galaxy_luminosity_multiple:,.0f} times the combined light of every star in a "
            f"galaxy like the Milky Way. Gas in its broad-line region, about "
            f"{self.broad_line_region_light_days:,.0f} light-days out, orbits fast enough to smear its "
            f"emission lines thousands of kilometers per second wide. This phase of activity began roughly "
            f"{self.active_age_years:,.0f} years ago."
        )

        paragraphs = [header, description]

        if self.is_radio_loud:
            paragraphs.append(
                f"It is radio-loud: twin relativistic jets punch out of the galaxy entirely, feeding radio "
                f"lobes some {self.jet_length_ly:,.0f} light-years from end to end."
            )
        else:
            paragraphs.append("It is radio-quiet, with no large-scale jets.")

        return paragraphs

    def __str__(self):
        """
        Returns a string representation of the quasar, with paragraphs
        separated by double newlines.
        """
        return "\n\n".join(self.to_paragraph_list())
