# stellarObjects/compactRemnant.py

"""
Compact Stellar Remnants
========================

This module defines `BlackHole` and `NeutronStar`, the two real end states
of a massive star's core collapse (supernova) that can, per
`docs/design/exotic-phenomena.md`, optionally anchor a full `StarSystem` in
place of an ordinary `Star` -- e.g. a pulsar with a fallback-disk planet
(real examples exist, PSR B1257+12), or a lone stellar-mass black hole with
nothing left orbiting it at all.

Both classes subclass `Star` the same way `doubleStar.BinaryStarProxy`
does: `_skip_property_init=True` skips `Star.__init__`'s own
`generate_star()`-driven property block, and this module fills in the same
attribute surface (`mass`, `radius`, `temperature`, `luminosity`, `age`,
`lifespan`, `habitable_zone`, `system_perimeter`, `heliosphere_radius`,
`galactic_orbital_*`) with values appropriate to a compact remnant instead.
This lets a `CompactRemnant` drop into every place `StarSystem` expects a
`Star` -- orbit placement, serialization, `adjust_age_for_planets` --
without any of that code needing to know compact remnants exist.

A compact remnant's `luminosity` is zero (or a small accretion-disk/thermal
value), which makes `utils.calculate_habitable_zone` naturally collapse its
`habitable_zone` to (0, 0) AU -- there is no thermal habitable zone around a
dark object -- without any special-casing in `StarSystem._generate_planets`.
Likewise, `StarSystem._estimate_max_objects_from_disk_physics`'s snow-line
-driven ceiling naturally comes out to zero planets for a zero-luminosity
remnant, matching the real rarity of confirmed planets around black
holes/neutron stars; `SystemConfig.PLANETS`/`NUM_ORBITS` can still force
orbiting bodies explicitly (see `phenomenonGen.py --num-orbits`).

Neither class models a remnant anchoring a binary system (a still-living
stellar companion, e.g. Cygnus X-1) -- `StarSystem.__init__`'s binary-
generation step is skipped entirely whenever a `compact_remnant` is passed
in. That is flagged as a natural follow-up, not required for this feature's
first version.
"""

import math
import random

from .config import SystemConfig
from . import physical_constants, program_constants
from .serialization import fields_from_dict, fields_to_dict
from .starData import Star
from .utils import (calculate_habitable_zone, format_age_string, format_galactic_orbit,
                    format_length_km, format_relative_to_sol, generate_galactic_orbit_fields,
                    properties_to_string, reseed_rng, to_scientific_notation)


class CompactRemnant(Star):
    """
    Shared base for `BlackHole`/`NeutronStar` -- both are `Star` subclasses
    (see module docstring) that generate an entirely different set of
    physical properties than `Star.generate_star()` would, then reuse this
    class's `_finish_init`/`_generate_remnant_age_and_lifespan`/`to_dict`/
    `from_dict` for the derived attributes and (de)serialization every
    `Star`-like object needs.

    Not instantiated directly -- always `BlackHole` or `NeutronStar`.
    """

    def __init__(self, system_config: SystemConfig, name=None, galactic_center_dist_ly=None):
        """
        Args:
            system_config (SystemConfig): The shared SystemConfig object.
            name (str, optional): Name for the remnant. Random if omitted.
            galactic_center_dist_ly (float, optional): See
                `Star.calculate_system_perimeter`'s docstring for the same
                parameter -- threaded through to this remnant's own Hill-
                sphere/galactic-orbit calculations.
        """
        super().__init__(system_config, name=name, _skip_property_init=True)
        self.galactic_center_dist_ly = galactic_center_dist_ly

    def _generate_remnant_age_and_lifespan(self):
        """
        Draws an age since the core-collapse supernova that formed this
        remnant, in billions of years
        (`program_constants.COMPACT_REMNANT_AGE_RANGE_GY`). Lifespan is
        `float('inf')` -- the same convention `Star` uses for white dwarfs
        -- since a black hole/neutron star doesn't have a further
        evolutionary endpoint this generator models; it simply persists.

        Returns:
            tuple: `(age_gy, lifespan_gy)`.
        """
        reseed_rng()
        age = random.uniform(*program_constants.COMPACT_REMNANT_AGE_RANGE_GY)
        return age, float('inf')

    def _finish_init(self, galactic_orbital_phase_deg=None):
        """
        Computes the derived attributes every `Star`-like object needs,
        once a subclass's `__init__` has already set `mass`, `radius`,
        `temperature`, `luminosity`, `type`, and `yerkes_class`. Mirrors
        the tail end of `Star.__init__`, with `heliosphere_radius` fixed
        at 0.0 rather than computed: `Star.calculate_heliosphere`'s
        stellar-wind mass-loss models (Nieuwenhuijzen & de Jager,
        radiation-driven O/B winds, coronal cool-dwarf winds) don't apply
        to a compact remnant's completely different physics (relativistic
        pulsar winds, accretion-disk outflows), and this generator has no
        model for those instead.

        Args:
            galactic_orbital_phase_deg (float, optional): See
                `Star.__init__`'s docstring for the same parameter.
        """
        self.habitable_zone = calculate_habitable_zone(self.luminosity)
        self.system_perimeter = self.calculate_system_perimeter(self.galactic_center_dist_ly)
        self.heliosphere_radius = 0.0
        (self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy,
         self.galactic_orbital_phase_deg, self.galactic_min_update_interval_years) = \
            generate_galactic_orbit_fields(self.galactic_center_dist_ly, galactic_orbital_phase_deg)

    def to_dict(self):
        """
        Returns a JSON-serializable dict of this remnant's properties, per
        `self.SERIALIZABLE_FIELDS` (the subclass's own list). Mirrors
        `Star.to_dict`'s `habitable_zone`/`lifespan` handling exactly.

        Returns:
            dict: One entry per field in `SERIALIZABLE_FIELDS`.
        """
        data = fields_to_dict(self, self.SERIALIZABLE_FIELDS)
        data["habitable_zone"] = list(self.habitable_zone)
        data["lifespan"] = None if self.lifespan == float('inf') else self.lifespan
        return data

    @classmethod
    def from_dict(cls, data, system_config):
        """
        Reconstructs a `BlackHole`/`NeutronStar` from a dict in the shape
        `to_dict()` produces, without re-running generation (`__init__` is
        bypassed via `object.__new__`, mirroring `Star.from_dict`).

        Args:
            data (dict): A dict in the shape `to_dict()` produces.
            system_config (SystemConfig): The system's shared config.

        Returns:
            CompactRemnant: The reconstructed remnant (concrete subclass
                matching `cls`).
        """
        remnant = object.__new__(cls)
        remnant.system_config = system_config
        fields_from_dict(remnant, data, cls.SERIALIZABLE_FIELDS)
        remnant.habitable_zone = tuple(data["habitable_zone"])
        remnant.lifespan = float('inf') if data["lifespan"] is None else data["lifespan"]
        return remnant

    def __str__(self):
        """Joins `to_paragraph_list()` with blank lines, like every other
        renderable object in this package."""
        return "\n\n".join(self.to_paragraph_list())


class BlackHole(CompactRemnant):
    """
    A stellar-mass (or, rarely, intermediate-mass) black hole: the
    collapsed core of a massive star that exceeded the maximum mass a
    neutron star can support. See module docstring for how it stands in
    for `Star` within a `StarSystem`.

    Attributes:
        mass_solar (float): Mass in solar masses.
        event_horizon_radius_km (float): Schwarzschild radius (`2GM/c^2`),
            in km -- also `self.radius`, the same attribute name `Star`
            uses for its own physical radius.
        spin (float): Dimensionless spin parameter a* in [0, 1).
        has_accretion_disk (bool): Whether this black hole retains a
            visible, actively-accreting disk.
    """

    SERIALIZABLE_FIELDS = Star.SERIALIZABLE_FIELDS + [
        "mass_solar", "event_horizon_radius_km", "spin", "has_accretion_disk",
    ]
    """Extends `Star.SERIALIZABLE_FIELDS` (minus `yerkes_class`'s normal
    meaning, repurposed here as the literal marker `"BH"` so `stars.
    yerkes_class` -- `NOT NULL` -- still holds a value) with this class's
    own fields. `mass`/`radius` are still serialized via the inherited
    names (they ARE `mass_solar`-in-kg / `event_horizon_radius_km`, kept
    as plain `Star`-compatible attributes so `StarSystem`/`_db.insert_star`
    need no special-casing)."""

    def __init__(self, system_config: SystemConfig, name=None, galactic_center_dist_ly=None,
                 galactic_orbital_phase_deg=None):
        super().__init__(system_config, name=name, galactic_center_dist_ly=galactic_center_dist_ly)
        reseed_rng()

        if random.random() < program_constants.BLACK_HOLE_INTERMEDIATE_MASS_CHANCE:
            self.mass_solar = random.uniform(*program_constants.BLACK_HOLE_INTERMEDIATE_MASS_RANGE_SOLAR)
        else:
            self.mass_solar = random.uniform(*program_constants.BLACK_HOLE_MASS_RANGE_SOLAR)
        self.mass = self.mass_solar * physical_constants.SOLAR_MASS_TO_KG

        # Schwarzschild radius: r_s = 2GM/c^2 (the non-rotating event
        # horizon radius -- spin technically shrinks this for a Kerr black
        # hole, but the Schwarzschild value is used here as the simple,
        # standard reference radius).
        schwarzschild_radius_m = (
            2 * physical_constants.G * self.mass / physical_constants.SPEED_OF_LIGHT_M_S ** 2
        )
        self.event_horizon_radius_km = schwarzschild_radius_m / 1000
        self.radius = self.event_horizon_radius_km

        self.spin = random.uniform(*program_constants.BLACK_HOLE_SPIN_RANGE)
        self.has_accretion_disk = random.random() < program_constants.BLACK_HOLE_ACCRETION_DISK_CHANCE

        if self.has_accretion_disk:
            # Eddington luminosity, L_edd = 1.26e31 * (M/Msun) W (standard
            # formula for the maximum luminosity a spherically-accreting
            # mass can sustain) -- an active disk is modeled as some
            # sub-Eddington fraction of this, and the inner-disk
            # temperature as a representative soft-X-ray value, both
            # order-of-magnitude flavor rather than a full accretion-disk
            # model (e.g. Shakura & Sunyaev 1973).
            eddington_luminosity_w = 1.26e31 * self.mass_solar
            self.luminosity = random.uniform(0.0001, 0.05) * eddington_luminosity_w
            self.temperature = random.uniform(1e5, 1e7)
        else:
            self.luminosity = 0.0
            self.temperature = 0.0

        mass_class = "Intermediate-Mass" if self.mass_solar >= program_constants.BLACK_HOLE_INTERMEDIATE_MASS_RANGE_SOLAR[0] else "Stellar-Mass"
        self.type = f"{mass_class} Black Hole"
        self.yerkes_class = "BH"

        self.age, self.lifespan = self._generate_remnant_age_and_lifespan()
        self._finish_init(galactic_orbital_phase_deg)

    def get_table_properties(self):
        """
        Builds the "Black Hole Data" property dict, mirroring
        `Star.get_table_properties`'s role for the database persistence
        layer and `to_paragraph_list`.

        Returns:
            dict: Keys `type`, `mass`, `event_horizon`, `spin`, `disk`,
                 `orbit`, `loc`, each an already-formatted display string.
        """
        mass_string = format_relative_to_sol(self.system_config, self.mass, physical_constants.SOLAR_MASS_TO_KG, "kg")
        radius_string = format_length_km(
            self.system_config, self.radius,
            program_constants.RADIUS_KM_SCIENTIFIC_NOTATION_THRESHOLD,
            program_constants.ROUND_RADIUS_KM, program_constants.SCIENTIFIC_NOTATION_DECIMAL_PLACES,
        )
        orbit_string = format_galactic_orbit(self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy)
        properties = {
            "type": self.type,
            "mass": mass_string,
            "event_horizon": radius_string,
            "spin": f"{self.spin:.3f} (dimensionless)",
            "disk": "Active accretion disk detected" if self.has_accretion_disk else "None detected",
            "orbit": orbit_string,
            "loc": self.name,
        }
        if self.reflex_offset_x or self.reflex_offset_y or self.reflex_offset_z:
            offset_km = math.sqrt(
                self.reflex_offset_x ** 2 + self.reflex_offset_y ** 2 + self.reflex_offset_z ** 2
            ) * physical_constants.AU_TO_KM
            properties["wobble"] = f"{to_scientific_notation(self.system_config, offset_km)} km from its nominal position, pulled by its own planets"
        return properties

    def to_paragraph_list(self):
        """
        Generates descriptive paragraphs for this black hole, analogous to
        `Star.to_paragraph_list`: a data block, then a narrative sentence
        covering its formation and current accretion state.

        Returns:
            list: A list of paragraph strings.
        """
        paragraphs = []
        properties = self.get_table_properties()
        markdown_key_map = {
            "type": "Type", "mass": "Mass", "event_horizon": "Event Horizon Radius",
            "spin": "Spin Parameter", "disk": "Accretion Disk", "orbit": "Galactic Orbit",
            "wobble": "Planetary Wobble", "loc": "Location",
        }
        paragraphs.append(properties_to_string(self.system_config, properties, "Black Hole Data", markdown_key_map=markdown_key_map))

        age_str = format_age_string(self.age)
        sentence = (
            f"{self.name} is the collapsed core of a massive star, formed in a supernova approximately "
            f"{age_str} ago. Nothing, not even light, can escape from within its event horizon."
        )
        if self.has_accretion_disk:
            sentence += (
                f" A faint accretion disk of infalling matter still surrounds it, its inner edge glowing at "
                f"roughly {self.temperature:,.0f} K."
            )
        else:
            sentence += " No companion star or debris disk remains; only its gravity betrays its presence."
        paragraphs.append(sentence)

        return paragraphs


class NeutronStar(CompactRemnant):
    """
    A neutron star: the collapsed, degenerate-neutron core left behind by a
    massive star's supernova when the remnant mass falls short of the
    threshold that would instead produce a black hole. See module
    docstring for how it stands in for `Star` within a `StarSystem`.

    Attributes:
        mass_solar (float): Mass in solar masses.
        spin_period_ms (float): Rotational period in milliseconds.
        magnetic_field_gauss (float): Surface magnetic field strength.
        pulsar_type (str): `"young"`, `"millisecond"`, or `"non-pulsing"`.
        surface_temperature_k (float): Surface temperature in Kelvin --
            also `self.temperature`.
    """

    SERIALIZABLE_FIELDS = Star.SERIALIZABLE_FIELDS + [
        "mass_solar", "spin_period_ms", "magnetic_field_gauss", "pulsar_type", "surface_temperature_k",
    ]
    """See `BlackHole.SERIALIZABLE_FIELDS`'s docstring for the same
    `Star.SERIALIZABLE_FIELDS`-extension convention; `yerkes_class` here is
    the literal marker `"NS"`."""

    def __init__(self, system_config: SystemConfig, name=None, galactic_center_dist_ly=None,
                 galactic_orbital_phase_deg=None):
        super().__init__(system_config, name=name, galactic_center_dist_ly=galactic_center_dist_ly)
        reseed_rng()

        self.mass_solar = random.uniform(*program_constants.NEUTRON_STAR_MASS_RANGE_SOLAR)
        self.mass = self.mass_solar * physical_constants.SOLAR_MASS_TO_KG
        self.radius = random.uniform(*program_constants.NEUTRON_STAR_RADIUS_RANGE_KM)
        self.surface_temperature_k = random.uniform(*program_constants.NEUTRON_STAR_SURFACE_TEMPERATURE_RANGE_K)
        self.temperature = self.surface_temperature_k

        is_pulsar = random.random() < program_constants.NEUTRON_STAR_PULSAR_CHANCE
        if is_pulsar and random.random() < program_constants.PULSAR_MILLISECOND_CHANCE:
            self.pulsar_type = "millisecond"
            self.spin_period_ms = random.uniform(*program_constants.PULSAR_SPIN_PERIOD_MS_RANGE_MILLISECOND)
            self.magnetic_field_gauss = random.uniform(*program_constants.PULSAR_MAGNETIC_FIELD_GAUSS_RANGE_MILLISECOND)
        elif is_pulsar:
            self.pulsar_type = "young"
            self.spin_period_ms = random.uniform(*program_constants.PULSAR_SPIN_PERIOD_MS_RANGE_YOUNG)
            self.magnetic_field_gauss = random.uniform(*program_constants.PULSAR_MAGNETIC_FIELD_GAUSS_RANGE_YOUNG)
        else:
            self.pulsar_type = "non-pulsing"
            self.spin_period_ms = random.uniform(*program_constants.PULSAR_SPIN_PERIOD_MS_RANGE_YOUNG)
            self.magnetic_field_gauss = random.uniform(*program_constants.PULSAR_MAGNETIC_FIELD_GAUSS_RANGE_YOUNG)

        # Blackbody thermal luminosity from the surface (Stefan-Boltzmann
        # law), the same physical relationship Star.generate_star's white-
        # dwarf branch and calculate_heliosphere both already lean on
        # elsewhere in this package -- a neutron star's small radius keeps
        # this a tiny fraction of Sol's luminosity even at its hottest.
        radius_m = self.radius * physical_constants.KM_TO_M_FACTOR
        self.luminosity = (
            physical_constants.STEFAN_BOLTZMANN_CONSTANT * 4 * math.pi * radius_m ** 2
            * self.surface_temperature_k ** 4
        )

        self.type = "Neutron Star" if self.pulsar_type == "non-pulsing" else f"Neutron Star ({self.pulsar_type.capitalize()} Pulsar)"
        self.yerkes_class = "NS"

        self.age, self.lifespan = self._generate_remnant_age_and_lifespan()
        self._finish_init(galactic_orbital_phase_deg)

    def get_table_properties(self):
        """
        Builds the "Neutron Star Data" property dict, mirroring
        `Star.get_table_properties`'s role for the database persistence
        layer and `to_paragraph_list`.

        Returns:
            dict: Keys `type`, `mass`, `radius`, `spin_period`,
                 `magnetic_field`, `surface_temp`, `orbit`, `loc`, each an
                 already-formatted display string.
        """
        mass_string = format_relative_to_sol(self.system_config, self.mass, physical_constants.SOLAR_MASS_TO_KG, "kg")
        radius_string = format_length_km(
            self.system_config, self.radius,
            program_constants.RADIUS_KM_SCIENTIFIC_NOTATION_THRESHOLD,
            program_constants.ROUND_RADIUS_KM, program_constants.SCIENTIFIC_NOTATION_DECIMAL_PLACES,
        )
        orbit_string = format_galactic_orbit(self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy)
        properties = {
            "type": self.type,
            "mass": mass_string,
            "radius": radius_string,
            "spin_period": f"{self.spin_period_ms:.2f} ms",
            "magnetic_field": f"{self.magnetic_field_gauss:.2e} G",
            "surface_temp": f"{self.surface_temperature_k:,.0f} K",
            "orbit": orbit_string,
            "loc": self.name,
        }
        if self.reflex_offset_x or self.reflex_offset_y or self.reflex_offset_z:
            offset_km = math.sqrt(
                self.reflex_offset_x ** 2 + self.reflex_offset_y ** 2 + self.reflex_offset_z ** 2
            ) * physical_constants.AU_TO_KM
            properties["wobble"] = f"{to_scientific_notation(self.system_config, offset_km)} km from its nominal position, pulled by its own planets"
        return properties

    def to_paragraph_list(self):
        """
        Generates descriptive paragraphs for this neutron star, analogous
        to `Star.to_paragraph_list`: a data block, then a narrative
        sentence covering its formation and pulsar behavior (if any).

        Returns:
            list: A list of paragraph strings.
        """
        paragraphs = []
        properties = self.get_table_properties()
        markdown_key_map = {
            "type": "Type", "mass": "Mass", "radius": "Radius",
            "spin_period": "Spin Period", "magnetic_field": "Magnetic Field",
            "surface_temp": "Surface Temperature", "orbit": "Galactic Orbit",
            "wobble": "Planetary Wobble", "loc": "Location",
        }
        paragraphs.append(properties_to_string(self.system_config, properties, "Neutron Star Data", markdown_key_map=markdown_key_map))

        age_str = format_age_string(self.age)
        sentence = (
            f"{self.name} is a neutron star, the collapsed core of a massive star left behind by a supernova "
            f"approximately {age_str} ago, packing more mass than the Sun into a sphere only "
            f"{self.radius:,.1f} km across."
        )
        if self.pulsar_type != "non-pulsing":
            sentence += (
                f" It is a {self.pulsar_type} pulsar, sweeping a beam of radiation past any observer once every "
                f"{self.spin_period_ms:.2f} milliseconds."
            )
        else:
            sentence += " Whatever rotational beam it once had no longer sweeps across this vantage point."
        paragraphs.append(sentence)

        return paragraphs
