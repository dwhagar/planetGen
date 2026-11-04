# planetgen/generation/phenomena/compact_remnant.py

"""
Compact Stellar Remnants
========================

This module defines `BlackHole` and `NeutronStar`, the two real end states
of a massive star's core collapse (supernova) that can, per
`docs/design/systems/anomalies.md` and `docs/design/systems/multistar-and-compact-systems.md`
(section 7), optionally anchor a full `StarSystem` in
place of an ordinary `Star` -- e.g. a pulsar with a fallback-disk planet
(real examples exist, PSR B1257+12), or a lone stellar-mass black hole with
nothing left orbiting it at all.

Both classes subclass `StellarObject` (`generation/stellar_object.py`), the
family class `Star` and `binary.BinaryStarProxy` share, and this module fills
in the same attribute surface (`mass`, `radius`, `temperature`, `luminosity`, `age`,
`lifespan`, `habitable_zone`, `system_perimeter`, `heliosphere_radius`,
`galactic_orbital_*`) with values appropriate to a compact remnant instead.
This lets a `CompactRemnant` drop into every place `StarSystem` expects a
`Star` -- orbit placement, serialization, `adjust_age_for_planets` --
without any of that code needing to know compact remnants exist.

A compact remnant's `luminosity` is tiny (a black hole's Hawking
radiation, a neutron star's thermal glow) unless it has an accretion disk,
which makes `orbits.calculate_habitable_zone` naturally collapse its
`habitable_zone` to (0, 0) AU -- there is no thermal habitable zone around a
dark object -- without any special-casing in `StarSystem._generate_planets`.
Likewise, the planet architecture's disc (`StarSystem._plan_architecture`,
out to a multiple of the snow line) naturally comes out to zero planets for
such a dim remnant, matching the real rarity of confirmed planets around black
holes/neutron stars; `SystemConfig.PLANETS`/`NUM_ORBITS` can still force
orbiting bodies explicitly (see `phenomenonGen.py --num-orbits`).

Neither class models a remnant anchoring a binary system (a still-living
stellar companion, e.g. Cygnus X-1) -- `StarSystem.__init__`'s binary-
generation step is skipped entirely whenever a `compact_remnant` is passed
in. That is flagged as a natural follow-up, not required for this feature's
first version.
"""

import math

from planetgen.generation.config import SystemConfig
from planetgen.physics import activity, constants, formation, neutron_star, relations, remnant_mass, spin
from planetgen import tuning
from planetgen.util import draw
from planetgen.util import log
from planetgen.generation.mass_object import MassObject
from planetgen.generation.stellar_object import StellarObject
from planetgen.galaxy import shape_draw as shapeDraw, structure as galaxyStructure
from planetgen.galaxy.remnant_distribution import pulsar_radial_factor
from planetgen.physics.orbits import calculate_habitable_zone
from planetgen.physics.remnant_progenitors import collapse_age_gy
from planetgen.util.random import log_uniform
from planetgen.physics import unit_registry as units

_kms_to_ms = units.converter("km_s", "m_s")
_m_to_au = units.converter("m", "au")
_au_to_pc = units.converter("au", "pc")
_ly_to_au = units.converter("ly", "au")
_km_to_m = units.converter("km", "m")



def hawking_temperature_k(mass_solar):
    """
    A black hole's Hawking temperature, T = hbar c^3 / (8 pi G M k_B):
    about 6.2e-8 K for one solar mass, falling as 1/M (GEN.82).
    """
    mass_kg = _msun_to_kg(mass_solar)
    return (constants.REDUCED_PLANCK * constants.SPEED_OF_LIGHT_M_S ** 3
            / (8 * math.pi * constants.G * mass_kg * constants.BOLTZMANN))


def hawking_luminosity_w(mass_solar):
    """
    A black hole's Hawking luminosity, L = hbar c^6 / (15360 pi G^2 M^2)
    (photons only): about 9e-29 W for one solar mass, falling as 1/M^2
    (GEN.82).
    """
    mass_kg = _msun_to_kg(mass_solar)
    return (constants.REDUCED_PLANCK * constants.SPEED_OF_LIGHT_M_S ** 6
            / (15360 * math.pi * constants.G ** 2 * mass_kg ** 2))


def infer_black_hole_mass_class(mass_solar):
    """The mass class a black hole of `mass_solar` falls in, for rows
    saved before `mass_class` was stored (schema v36)."""
    if mass_solar >= tuning.BLACK_HOLE_SUPERMASSIVE_MIN_MASS_SOLAR:
        return "supermassive"
    if mass_solar >= tuning.BLACK_HOLE_INTERMEDIATE_MASS_RANGE_SOLAR[0]:
        return "intermediate"
    return "stellar"


class CompactRemnant(StellarObject):
    """
    Shared base for `BlackHole`/`NeutronStar` -- both are `StellarObject` subclasses
    (see module docstring) that generate an entirely different set of
    physical properties than `Star.generate_star()` would, then reuse this
    class's `_finish_init`/`_generate_remnant_age_and_lifespan`/`to_dict`/
    `from_dict` for the derived attributes and (de)serialization every
    `Star`-like object needs.

    Not instantiated directly -- always `BlackHole` or `NeutronStar`.
    """

    # A standalone black hole's or neutron star's name is a draft:
    # `store.insert_black_hole`/`insert_neutron_star` reserve it through the
    # system-name registry (v40), as quasars and rogue planets do.
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
        super().__init__(system_config, name=name, galactic_center_dist_ly=galactic_center_dist_ly)

    ejecta_velocity_kms = None
    """tuple or None: The galactic-axes velocity of a remnant the galaxy ejected (`eject`, GEN.238), km/s; `None` for a
    remnant that orbits the galaxy."""

    def eject(self, velocity_kms):
        """
        Makes this remnant one the galaxy ejected (GEN.238): it moves on a straight line at `velocity_kms`, galactic
        axes, and has no galactic orbit, so its orbit fields are `None`, as the galaxy's central black hole's are.
        """
        self.ejecta_velocity_kms = tuple(velocity_kms)
        self.galactic_orbital_speed_kms = self.galactic_orbital_period_gy = None
        self.galactic_orbital_phase_deg = self.galactic_min_update_interval_years = None

    def _generate_remnant_age_and_lifespan(self, kind=None):
        """
        Draws an age since the core-collapse supernova that formed this
        remnant, in billions of years. A neutron star or stellar black hole
        (`kind`, GEN.553) takes the age of a remnant in its population: the
        progenitor's birth from the population's star formation history less
        the time it took to collapse (`remnant_progenitors.collapse_age_gy`); a
        black hole of another class, which no star collapsed to, draws one
        uniform over `tuning.COMPACT_REMNANT_AGE_RANGE_GY`. Lifespan is
        `float('inf')` -- the same convention `Star` uses for white dwarfs
        -- since a black hole/neutron star doesn't have a further
        evolutionary endpoint this generator models; it simply persists.

        Returns:
            tuple: `(age_gy, lifespan_gy)`.
        """
        config = self.system_config
        if kind is None:
            age = draw.uniform(*tuning.COMPACT_REMNANT_AGE_RANGE_GY)
        else:
            age = collapse_age_gy(kind, config.AGE, draw, getattr(config, "POPULATION", None),
                                  getattr(config, "AGE_LAW", None))
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
        self.habitable_zone = calculate_habitable_zone(self.luminosity, self.temperature)
        self.system_perimeter = self.calculate_system_perimeter(self.galactic_center_dist_ly)
        self.heliosphere_radius = 0.0
        self.set_galactic_motion(self.galactic_center_dist_ly, galactic_orbital_phase_deg)
        # GEN.458: this remnant's own snow line, from its own stream like a star's.
        self.snow_line_au = formation.draw_snow_line_au(
            self.luminosity, draw.Stream(f"snow_line:{self.mass!r}:{self.age!r}:{self.galactic_center_dist_ly!r}"))
        # GEN.104: the spin axis, a Rayleigh tilt from the galactic pole.
        self.rotation_period_hours = self._rotation_period_hours()
        spin.set_spin(self, spin.GALACTIC_POLE, spin.rayleigh_tilt_deg())
        # GEN.86: no corona or flares; what the remnant gives off as XUV.
        activity.generate_remnant_activity(self, self._xuv_w())

    def _rotation_period_hours(self):
        """This remnant's rotation period, hours (GEN.104)."""
        raise NotImplementedError

    def _xuv_w(self):
        """This remnant's XUV luminosity, W (GEN.86): a black hole's disk
        (or Hawking) luminosity; a neutron star overrides it."""
        return self.luminosity or 0.0


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

    FIELD_UNITS = {
        **StellarObject.FIELD_UNITS, "mass_solar": "Msun", "event_horizon_radius_km": "km", "spin": "ratio",
        "has_accretion_disk": "flag", "mass_class": "text",
    }
    SERIALIZABLE_FIELDS = list(FIELD_UNITS)
    """Extends `StellarObject.SERIALIZABLE_FIELDS` (minus `yerkes_class`'s normal
    meaning, repurposed here as the literal marker `"BH"` so `stars.
    yerkes_class` -- `NOT NULL` -- still holds a value) with this class's
    own fields. `mass`/`radius` are still serialized via the inherited
    names (they ARE `mass_solar`-in-kg / `event_horizon_radius_km`, kept
    as plain `Star`-compatible attributes so `StarSystem`/`store.insert_star`
    need no special-casing)."""

    def __init__(self, system_config: SystemConfig, name=None, galactic_center_dist_ly=None,
                 galactic_orbital_phase_deg=None, mass_class=None, mass_range=None, mass_solar=None):
        """
        Args:
            mass_class (str, optional): `"supermassive"` for a galaxy's
                quiescent central black hole (`generate.add_galactic_nucleus`
                -- it sits at the center, so it has no galactic orbit).
                `"stellar"` or `"intermediate"` takes that class without a
                roll (the phenomenon scatter draws them as separate kinds,
                GEN.166). `None` (the default) rolls stellar-mass, or
                intermediate-mass with `BLACK_HOLE_INTERMEDIATE_MASS_CHANCE`.
            mass_range (tuple, optional): `(low, high)` solar masses inside
                the class's own range: the mass is drawn from the class's
                law truncated to it (the scatter's mass cut, GEN.166).
            mass_solar (float, optional): The mass itself, when the phenomenon
                scatter already drew it (MAP.165); no mass is drawn then.
        """
        super().__init__(system_config, name=name, galactic_center_dist_ly=galactic_center_dist_ly)
        self.name_given = bool(name)  # a given name is kept over an object ID (GEN.64)

        if mass_class not in (None, "supermassive", "stellar", "intermediate"):
            raise ValueError(f"mass_class must be None, 'supermassive', 'stellar' or 'intermediate', "
                             f"got {mass_class!r}")
        if mass_class == "supermassive":
            self.mass_class = "supermassive"
            self.mass_solar = mass_solar or shapeDraw.random_nucleus_mass_solar()
        else:
            if mass_class is None:
                intermediate = draw.random() < tuning.BLACK_HOLE_INTERMEDIATE_MASS_CHANCE
                log.choice("Black hole mass regime", "intermediate-mass" if intermediate else "stellar-mass",
                           f"roll against BLACK_HOLE_INTERMEDIATE_MASS_CHANCE "
                           f"({tuning.BLACK_HOLE_INTERMEDIATE_MASS_CHANCE})")
                mass_class = "intermediate" if intermediate else "stellar"
            self.mass_class = mass_class
            if mass_solar:
                self.mass_solar = mass_solar
            elif mass_class == "intermediate":
                self.mass_solar = log_uniform(*remnant_mass.within(tuning.BLACK_HOLE_INTERMEDIATE_MASS_RANGE_SOLAR, mass_range))
            else:
                self.mass_solar = remnant_mass.draw_mass_sol("black-hole", mass_range)
        self.mass = _msun_to_kg(self.mass_solar)

        # Schwarzschild radius: r_s = 2GM/c^2 (the non-rotating event
        # horizon radius -- spin technically shrinks this for a Kerr black
        # hole, but the Schwarzschild value is used here as the simple,
        # standard reference radius).
        schwarzschild_radius_m = relations.schwarzschild_radius(self.mass)
        self.event_horizon_radius_km = schwarzschild_radius_m / 1000
        self.radius = self.event_horizon_radius_km

        self.spin = spin.black_hole_spin()
        if self.mass_class == "supermassive":
            # Always some accretion flow, far below Eddington (Sgr A*).
            self.has_accretion_disk = True
            eddington_ratio = log_uniform(*tuning.BLACK_HOLE_SUPERMASSIVE_EDDINGTON_RATIO_RANGE)
            self.luminosity = (
                eddington_ratio * tuning.EDDINGTON_LUMINOSITY_W_PER_SOLAR_MASS * self.mass_solar
            )
            self.temperature = draw.uniform(1e5, 1e7)
        else:
            self.has_accretion_disk = draw.random() < tuning.BLACK_HOLE_ACCRETION_DISK_CHANCE
            log.choice("Accretion disk", self.has_accretion_disk,
                       f"roll against BLACK_HOLE_ACCRETION_DISK_CHANCE "
                       f"({tuning.BLACK_HOLE_ACCRETION_DISK_CHANCE})")

            if self.has_accretion_disk:
                # Eddington luminosity, L_edd = 1.26e31 * (M/Msun) W (standard
                # formula for the maximum luminosity a spherically-accreting
                # mass can sustain) -- an active disk is modeled as some
                # sub-Eddington fraction of this, and the inner-disk
                # temperature as a representative soft-X-ray value, both
                # order-of-magnitude flavor rather than a full accretion-disk
                # model (e.g. [R4037]).
                eddington_luminosity_w = 1.26e31 * self.mass_solar
                self.luminosity = draw.uniform(0.0001, 0.05) * eddington_luminosity_w
                self.temperature = draw.uniform(1e5, 1e7)
            else:
                # No disk: all it gives off is Hawking radiation, real but
                # tiny (it once stored a flat 0, GEN.82).
                self.luminosity = hawking_luminosity_w(self.mass_solar)
                self.temperature = hawking_temperature_k(self.mass_solar)

        self.type = f"{self.mass_class_label} Black Hole"
        self.yerkes_class = "BH"

        self.age, self.lifespan = self._generate_remnant_age_and_lifespan(
            "black-hole" if self.mass_class == "stellar" else None)
        self._finish_init(galactic_orbital_phase_deg)
        if self.mass_class == "supermassive":
            # The galaxy's center: nothing to orbit, and its reach is its
            # sphere of influence (G*M / sigma^2), not a galactic Hill sphere.
            (self.galactic_orbital_speed_kms, self.galactic_orbital_period_gy,
             self.galactic_orbital_phase_deg, self.galactic_min_update_interval_years) = None, None, None, None
            sigma_m_s = _kms_to_ms(galaxyStructure.black_hole_dispersion_kms(self.mass_solar))
            self.system_perimeter = _m_to_au(constants.G * self.mass / sigma_m_s ** 2)

    @property
    def mass_class_label(self):
        """`"Stellar-Mass"`, `"Intermediate-Mass"` or `"Supermassive"`."""
        return {"stellar": "Stellar-Mass", "intermediate": "Intermediate-Mass",
                "supermassive": "Supermassive"}[self.mass_class]

    @classmethod
    def from_dict(cls, data, system_config):
        """`CompactRemnant.from_dict`, with `mass_class` inferred from the
        mass for data saved before it existed (schema v36)."""
        black_hole = super().from_dict(data, system_config)
        if getattr(black_hole, "mass_class", None) is None:
            black_hole.mass_class = infer_black_hole_mass_class(black_hole.mass_solar)
        return black_hole

    @property
    def spin_rate_rad_s(self):
        """The horizon's rotation rate, rad/s; zero for a* = 0 (no rotation)."""
        return MassObject.spin_rate_rad_s.fget(self) or 0.0

    def _rotation_period_hours(self):
        """The horizon's rotation period, hours; `None` for a* = 0 (no
        rotation)."""
        period_s = spin.black_hole_horizon_period_s(self.mass, self.spin)
        return period_s / 3600.0 if math.isfinite(period_s) else None


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
        pulsar_type (str): `"young"`, `"millisecond"`, `"magnetar"` or
            `"non-pulsing"` (GEN.153).
        surface_temperature_k (float): Surface temperature in Kelvin --
            also `self.temperature`.
    """

    FIELD_UNITS = {
        **StellarObject.FIELD_UNITS, "mass_solar": "Msun", "spin_period_ms": "ms", "magnetic_field_gauss": "G",
        "pulsar_type": "text", "surface_temperature_k": "K",
    }
    SERIALIZABLE_FIELDS = list(FIELD_UNITS)
    """See `BlackHole.SERIALIZABLE_FIELDS`'s docstring for the same
    `StellarObject.SERIALIZABLE_FIELDS`-extension convention; `yerkes_class` here is
    the literal marker `"NS"`."""

    def _pulsar_chance(self):
        """The chance this neutron star is an active young pulsar: it falls
        off with age as `exp(-age / PULSAR_ACTIVE_LIFETIME_YR)`, scaled by
        `PULSAR_RATE_SCALE` and the Lorimer radial profile at its
        galactocentric distance (GEN.132, unchanged when the distance is
        unknown)."""
        chance = math.exp(-self.age * 1e9 / tuning.PULSAR_ACTIVE_LIFETIME_YR) * tuning.PULSAR_RATE_SCALE
        if not self.galactic_center_dist_ly:
            return min(1.0, chance)
        radius_pc = _au_to_pc(_ly_to_au(self.galactic_center_dist_ly))
        return min(1.0, chance * pulsar_radial_factor(radius_pc))

    def __init__(self, system_config: SystemConfig, name=None, galactic_center_dist_ly=None,
                 galactic_orbital_phase_deg=None, mass_range=None, mass_solar=None, age_gy=None, magnetar=False):
        """`mass_range` (`(low, high)` solar masses, optional) truncates the
        mass draw to that part of `NEUTRON_STAR_MASS_RANGE_SOLAR` (the
        scatter's mass cut, GEN.166). `mass_solar` is a mass the phenomenon
        scatter already drew (MAP.165): no mass is drawn then. `age_gy` is
        the age when it is known (the core of a supernova remnant is as old
        as the remnant); drawn when omitted. `magnetar` makes it one (a
        scatter row of the magnetar class, GEN.278), with an age within a
        magnetar's active life unless `age_gy` is given."""
        super().__init__(system_config, name=name, galactic_center_dist_ly=galactic_center_dist_ly)
        self.name_given = bool(name)  # a given name is kept over an object ID (GEN.64)

        self.mass_solar = mass_solar or remnant_mass.draw_mass_sol("neutron-star", mass_range)
        self.mass = _msun_to_kg(self.mass_solar)
        self.radius = neutron_star.radius_km(self.mass_solar, draw.gauss())
        self.surface_temperature_k = draw.uniform(*tuning.NEUTRON_STAR_SURFACE_TEMPERATURE_RANGE_K)
        self.temperature = self.surface_temperature_k
        if age_gy is None and magnetar:
            self.age, self.lifespan = draw.uniform(*tuning.MAGNETAR_AGE_RANGE_GY), float('inf')
        elif age_gy is None:
            self.age, self.lifespan = self._generate_remnant_age_and_lifespan("neutron-star")
        else:
            self.age, self.lifespan = age_gy, float('inf')

        # A star younger than a magnetar's active life may have been born one (a
        # scattered magnetar is one by its row); failing that it pulses by the age-dependent chance, and an old one
        # only as a recycled millisecond pulsar.
        is_magnetar = magnetar or (self.age * 1e9 < tuning.MAGNETAR_ACTIVE_LIFETIME_YR
                                   and draw.random() < tuning.MAGNETAR_BIRTH_FRACTION)
        is_young_pulsar = not is_magnetar and draw.random() < self._pulsar_chance()
        is_recycled = (not is_magnetar and not is_young_pulsar
                       and draw.random() < tuning.PULSAR_RECYCLED_CHANCE * tuning.PULSAR_RATE_SCALE)
        if is_magnetar:
            self.pulsar_type = "magnetar"
            self.spin_period_ms = draw.uniform(*tuning.MAGNETAR_SPIN_PERIOD_MS_RANGE)
            low, high = tuning.MAGNETAR_MAGNETIC_FIELD_GAUSS_RANGE
            self.magnetic_field_gauss = math.exp(draw.uniform(math.log(low), math.log(high)))
        elif is_recycled:
            self.pulsar_type = "millisecond"
            self.spin_period_ms = draw.uniform(*tuning.PULSAR_SPIN_PERIOD_MS_RANGE_MILLISECOND)
            self.magnetic_field_gauss = draw.uniform(*tuning.PULSAR_MAGNETIC_FIELD_GAUSS_RANGE_MILLISECOND)
        elif is_young_pulsar:
            self.pulsar_type = "young"
            self.spin_period_ms = draw.uniform(*tuning.PULSAR_SPIN_PERIOD_MS_RANGE_YOUNG)
            self.magnetic_field_gauss = draw.uniform(*tuning.PULSAR_MAGNETIC_FIELD_GAUSS_RANGE_YOUNG)
        else:
            self.pulsar_type = "non-pulsing"
            self.spin_period_ms = draw.uniform(*tuning.PULSAR_SPIN_PERIOD_MS_RANGE_YOUNG)
            self.magnetic_field_gauss = draw.uniform(*tuning.PULSAR_MAGNETIC_FIELD_GAUSS_RANGE_YOUNG)
        log.choice("Pulsar type", self.pulsar_type,
                   f"age {self.age:.3g} Gy: magnetar={is_magnetar} (MAGNETAR_BIRTH_FRACTION="
                   f"{tuning.MAGNETAR_BIRTH_FRACTION}), young pulsar={is_young_pulsar} "
                   f"(PULSAR_ACTIVE_LIFETIME_YR={tuning.PULSAR_ACTIVE_LIFETIME_YR:g}), "
                   f"recycled={is_recycled} (PULSAR_RECYCLED_CHANCE={tuning.PULSAR_RECYCLED_CHANCE})")

        # Blackbody thermal luminosity from the surface (Stefan-Boltzmann
        # law), the same physical relationship Star.generate_star's white-
        # dwarf branch and calculate_heliosphere both already lean on
        # elsewhere in this package -- a neutron star's small radius keeps
        # this a tiny fraction of Sol's luminosity even at its hottest.
        radius_m = _km_to_m(self.radius)
        self.luminosity = (
            constants.STEFAN_BOLTZMANN_CONSTANT * 4 * math.pi * radius_m ** 2
            * self.surface_temperature_k ** 4
        )

        if self.pulsar_type == "non-pulsing":
            self.type = "Neutron Star"
        elif self.pulsar_type == "magnetar":
            self.type = "Neutron Star (Magnetar)"
        else:
            self.type = f"Neutron Star ({self.pulsar_type.capitalize()} Pulsar)"
        self.yerkes_class = "NS"

        self._finish_init(galactic_orbital_phase_deg)

    def _rotation_period_hours(self):
        """`spin_period_ms` in hours."""
        return self.spin_period_ms / 3.6e6

    @property
    def spin_down_power_w(self):
        """The magnetic-dipole spin-down power, W (`activity.neutron_star_spin_down_w`), the pulsar wind's power
        and the source of its X-rays (GEN.317)."""
        return activity.neutron_star_spin_down_w(self)

    def _xuv_w(self):
        """Thermal X-rays or spin-down X-rays (`activity.neutron_star_xuv_w`)."""
        return activity.neutron_star_xuv_w(self)

from planetgen.generation import object_timing  # noqa: E402 -- registers the generators below (PERF.92)

_msun_to_kg = units.converter("Msun", "kg")
object_timing.track(BlackHole, "black-hole")
object_timing.track(NeutronStar, "neutron-star")
