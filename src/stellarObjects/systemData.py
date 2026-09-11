# stellarObjects/systemData.py

"""
Star System Generation
======================

This module contains the `StarSystem` class, which is used to generate and
represent a full star system, including a central star and a list of planets.
The `StarSystem` class orchestrates the creation of the star and its planets,
applying a set of procedural generation rules to create a diverse and
scientifically grounded star system.

The generation process can be customized through a variety of parameters,
allowing for the creation of specific types of systems, such as those with
habitable worlds, asteroid belts, or a particular star type. The class also
includes methods for validating the system's orbital mechanics, counting the
number of celestial objects, and generating a detailed string representation of
the system.
"""

import copy
import math
import random

from .asteroidData import AsteroidBelt
from .config import SystemConfig
from .doubleStar import BinaryStarProxy
from . import physical_constants, planetLife, planetPhysics, program_constants
from .planetData import Planet
from .starData import Star
from .utils import (
    disk_surface_density_scale,
    mmsn_surface_density_gcm2,
    isolation_mass_kg,
    snow_line_au,
    to_paragraph,
)

# Tracks the shape of `StarSystem.to_dict()`'s output (the serialized
# object-graph -- see TODO.md's Phase 1), independent of
# `stellarObjects._db.SCHEMA_VERSION`, which tracks the *database's* DDL
# structure instead. A JSON export carries this number so `from_dict` can
# raise a clear error if a file's shape is newer than the code understands,
# rather than a confusing `KeyError`/`AttributeError` partway through
# reconstruction.
SERIALIZATION_SCHEMA_VERSION = 1


class StarSystem:
    """
    A class representing a star system, containing a central star and a list of planets.

    The `StarSystem` class is the main entry point for generating a complete star
    system. It initializes a central star and then procedurally generates a list
    of planets and other celestial objects orbiting it. The generation process
    can be influenced by a variety of flags, allowing for fine-tuned control over
    the final system's characteristics.

    The class handles the complex logic of placing planets in stable orbits,
    ensuring that the system is physically plausible. It also includes methods
    for counting the number of different types of objects in the system and
    generating a detailed, human-readable summary of the system's properties.

    Attributes:
        star (Star): The central star of the system.
        planets (list): A list of `Planet` and `AsteroidBelt` objects orbiting the star.
        planet_count (int): The total number of planets in the system.
        belt_count (int): The total number of asteroid belts in the system.
        moon_count (int): The total number of moons in the system.
        hab_count (int): The total number of potentially habitable worlds.
        m_count (int): The total number of Class M worlds.
        system_flavor_count (int): Tracks the total number of flavor texts added across the system.
        system_flavor_text (str or None): The system-level flavor text sentence,
            decided once at generation time (or None if the roll didn't select one).
    """

    def __init__(self, system_config: SystemConfig, galactic_center_dist_ly=None):
        """
        Initializes a StarSystem object, generating a star and its planets.

        This constructor orchestrates the entire star system generation process.
        It begins by creating a `Star` instance, which can be customized with
        flags to control its size and type. It then estimates the number of
        celestial objects the system can support and proceeds to generate them,
        placing them in orbit around the star.

        The generation logic can be influenced by several parameters, allowing for
        fine-tuned control over the final system's characteristics. For example,
        `config.HABITABLE_WORLD = True` ensures that at least one habitable planet
        is generated, while `config.PLANETS = False` can be used to create a star
        with no orbiting bodies. `config.NUM_ORBITS` and `config.SLOTS` allow the
        number of orbits and the exact contents of specific orbital slots to be
        specified explicitly, overriding the random placement logic for those
        slots. The placement of planets is done sequentially, with each new
        planet's orbit being determined based on the position of the previous one
        to ensure a degree of realism in orbital spacing.

        If `config.NAME` is provided, it will be used as the name for the star
        system, overriding the default random name generation.

        Every `Planet` (and moon) is generated without life data — see
        `planetPhysics`'s module docstring. Once all planets/moons exist,
        orbits are validated, and the star's age has been finalized via
        `Star.adjust_age_for_planets`, this constructor makes one pass over
        every planet and moon and applies `planetLife.apply_life_data` and
        `planetLife.decide_flavor_text` to each, so evolutionary timelines
        reflect the star's final age and flavor text is a fixed fact by the
        time anything is ever rendered. The system-level flavor text
        (`self.system_flavor_text`) is decided the same way, once, right
        before this pass.

        Args:
            galactic_center_dist_ly (float, optional): This system's actual
                distance from the galactic center, in light-years -- passed
                straight through to every `Star`/`BinaryStarProxy` this
                system creates (see `Star.calculate_system_perimeter`'s
                docstring). `None` (the default) leaves every constituent
                star's Hill-sphere calculation on the fixed
                `physical_constants.GALACTIC_CENTER_DISTANCE_LY` constant --
                the correct behavior for a system with no galaxy placement
                (e.g. `sectorGen.py`'s own standalone CLI).
        """
        self.system_config = system_config # Assign the passed SystemConfig instance
        # Rolled once here (not left for each Star/BinaryStarProxy to roll
        # its own) and passed identically to every constituent star below --
        # a binary pair's AU-scale separation is negligible next to its
        # light-year-scale galactic orbit radius, so both stars (and the
        # proxy standing in for the pair) move around the galaxy together,
        # sharing one phase, not three independent ones. See
        # `Star.__init__`'s `galactic_orbital_phase_deg` docstring.
        galactic_orbital_phase_deg = random.uniform(0, 360)
        self.star = Star(self.system_config, name=self.system_config.NAME,
                          galactic_center_dist_ly=galactic_center_dist_ly,
                          galactic_orbital_phase_deg=galactic_orbital_phase_deg) # Pass system_config and use its NAME
        self.primary_star = self.star # For single star systems, the primary is the star
        self.planets = []
        self.stars = [self.primary_star] # Keep track of individual stars

        if self.system_config.BINARY_SYSTEM:
            # Create a copy of the system_config for the secondary star
            secondary_star_config = copy.deepcopy(self.system_config)
            # Ensure LARGE_STAR is not forced for the secondary star
            secondary_star_config.LARGE_STAR = False

            # Logic to generate a secondary star (e.g., random mass relative to primary)
            # This is new generation logic, but contained.
            # For simplicity, secondary mass is a fraction of primary mass.
            secondary_mass_factor = random.uniform(0.1, 0.8)
            secondary_mass = self.primary_star.mass * secondary_mass_factor
            # Create secondary star, potentially with a different name or type if desired
            self.secondary_star = Star(secondary_star_config, name=f"{self.primary_star.name} B",
                                        mass_override=secondary_mass,
                                        galactic_center_dist_ly=galactic_center_dist_ly,
                                        galactic_orbital_phase_deg=galactic_orbital_phase_deg)
            self.stars.append(self.secondary_star)
            self.star = BinaryStarProxy(self.system_config, self.primary_star, self.secondary_star,
                                         galactic_center_dist_ly=galactic_center_dist_ly,
                                         galactic_orbital_phase_deg=galactic_orbital_phase_deg) # self.star now points to the proxy

        # Removed: self.system_flavor_count = 0 # Initialize system flavor count
        # `validate_system`'s orbital-overlap correction can, in rare cases,
        # push a deliberately-placed guaranteed body (the forced habitable
        # world, or an asteroid belt's own guaranteed fallback slot) into a
        # zone/spacing its own placement no longer supports -- retry the
        # whole placement (same star, fresh positions) rather than accept a
        # system that silently fails a requirement `system_config` explicitly
        # asked for. See `_generate_planets`/`MAX_SYSTEM_GENERATION_ATTEMPTS`.
        for _attempt in range(program_constants.MAX_SYSTEM_GENERATION_ATTEMPTS):
            self.planets = []
            self._generate_planets()
            self.validate_system()

            habitable_satisfied = self.system_config.HABITABLE_WORLD is not True or self.count_habitable()[0] > 0
            belt_satisfied = self.system_config.ASTEROID_BELT is not True or self.count_objects()[1] > 0
            if habitable_satisfied and belt_satisfied:
                break

        self.star.adjust_age_for_planets(self.planets)

        # Decided once, here, at generation time -- __str__ (which may be called
        # more than once) reads this rather than re-rolling and double-counting
        # system_flavor_count on every render.
        self.system_flavor_text = None
        if random.random() < program_constants.FLAVOR_CHANCE_SYSTEM and self.system_config.system_flavor_count < program_constants.MAX_FLAVOR_TOTAL:
            self.system_flavor_text = random.choice(program_constants.SYSTEM_FLAVOR)
            self.system_config.system_flavor_count += 1

        # All planets and moons are generated without life data (see planetPhysics's
        # module docstring). Apply it now, in one pass over the finished system, so
        # evolutionary timelines are computed against the star's final, planet-adjusted
        # age rather than its provisional pre-adjustment one. Flavor text is decided
        # in the same pass (after life data, since it reads evolutionary_data) so it
        # too is a fixed, pre-rendered fact rather than something __str__ rolls.
        for obj in self.planets:
            if obj.body_type == 'a': # Skip asteroid belts; they carry no life data.
                continue
            planetLife.apply_life_data(obj)
            planetLife.decide_flavor_text(obj)
            for moon in obj.moons:
                planetLife.apply_life_data(moon)
                planetLife.decide_flavor_text(moon)

        self.planet_count, self.belt_count, self.moon_count = self.count_objects()
        self.hab_count, self.m_count = self.count_habitable()

    def _generate_planets(self):
        """
        Populates `self.planets` with a fresh, sequentially-placed set of
        planets/asteroid belts around `self.star` -- extracted out of
        `__init__` so it can be retried wholesale (same star, fresh
        positions) when the result doesn't satisfy an explicitly requested
        `HABITABLE_WORLD`/`ASTEROID_BELT` guarantee (see `__init__`'s own
        retry loop). Appends directly to `self.planets`, which the caller
        is responsible for resetting to `[]` first.
        """
        system_objects = self.estimate_num_objects()
        star_factor = self.star.mass / physical_constants.SOLAR_MASS_TO_KG

        required_objects = 0
        if self.system_config.HABITABLE_WORLD is True:
            required_objects += 1
        if self.system_config.ASTEROID_BELT is True:
            required_objects += 1

        if system_objects < required_objects:
            system_objects = required_objects

        slots = self.system_config.SLOTS or []

        if system_objects > 0:
            belt_index = random.randint(0, system_objects - 1) if self.system_config.ASTEROID_BELT is True else -1
            found_hab = False
            found_belt = False
            i = -1

            while i < system_objects - 1:
                i += 1
                last_asteroid = False
                slot_spec = slots[i] if i < len(slots) else None
                prev_slot_explicit = i > 0 and (i - 1) < len(slots) and slots[i - 1] is not None

                if i > 0:
                    last_planet = self.planets[i - 1]
                    random_buffer = random.uniform(0, star_factor)
                    if last_planet.body_type == 'a':
                        estimated_distance = last_planet.upper_limit + random_buffer * 2
                        last_asteroid = True
                    else:
                        estimated_distance = (last_planet.distance + last_planet.min_orbit_distance) + random_buffer

                    # Each successive slot's minimum spacing scales with the
                    # previous object's own Hill radius, which itself scales
                    # with its distance -- for a system with many objects
                    # (especially several giant planets around a massive
                    # star), that compounds geometrically. `system_perimeter`
                    # (the star's own Hill sphere *relative to the galaxy*,
                    # see `Star.calculate_system_perimeter`) is the real
                    # physical boundary past which nothing is gravitationally
                    # bound to this star at all -- stop adding slots once
                    # the sequential spacing would place one beyond it,
                    # rather than letting that compounding run unbounded.
                    # A system that runs out of stable room this way simply
                    # ends up with fewer objects than `system_objects`
                    # estimated, the same physically-honest outcome a real
                    # protoplanetary disk of finite extent would produce.
                    if estimated_distance > self.star.system_perimeter:
                        break
                else:
                    estimated_distance = program_constants.INITIAL_PLANET_DISTANCE_FACTOR * star_factor

                hz = self.star.habitable_zone[0] < estimated_distance < self.star.habitable_zone[1]

                # An explicit per-slot specification takes priority over all of the
                # normal random/forced generation logic below.
                if slot_spec is not None:
                    obj = self.generate_slot_object(slot_spec, estimated_distance)
                    if getattr(obj, 'planet_class', None) in program_constants.HABITABLE_PLANET_CLASSES:
                        found_hab = True
                    if getattr(obj, 'body_type', None) == 'a':
                        found_belt = True
                    self.planets.append(obj)
                    continue

                if self.system_config.HABITABLE_WORLD is True and not found_hab:
                    if not hz and i == 0:
                        if (estimated_distance > self.star.habitable_zone[1] or
                                0 < self.star.habitable_zone[0] - estimated_distance < 0.2 or system_objects == 1):
                            estimated_distance = self._forced_habitable_distance()
                            hz = True
                    elif not hz and i > 0:
                        last_planet = self.planets[i - 1]
                        beyond_hz = last_planet.body_type == 'a' and last_planet.upper_limit > self.star.habitable_zone[1] or \
                                    last_planet.body_type != 'a' and (last_planet.distance + last_planet.min_orbit_distance >
                                                                 self.star.habitable_zone[1])

                        # Never retroactively overwrite a belt if ASTEROID_BELT is
                        # required: this could be the belt satisfying that guarantee
                        # (found_belt is set the moment it's placed, not re-checked
                        # here), and this mechanism has no way to know whether another
                        # one exists elsewhere to fall back on.
                        belt_is_protected = last_planet.body_type == 'a' and self.system_config.ASTEROID_BELT is True

                        if beyond_hz and not prev_slot_explicit and not belt_is_protected:
                            estimated_distance = self._forced_habitable_distance()
                            planet = Planet(self.system_config, self.star, self.star.habitable_zone, estimated_distance, # Pass system_config
                                            planet_class="M")
                            self.planets[i - 1] = planet
                            i -= 1
                            found_hab = True
                            continue
                        elif i == system_objects - 1:
                            estimated_distance = self._forced_habitable_distance()
                            hz = True

                    if hz:
                        planet = Planet(self.system_config, self.star, self.star.habitable_zone, estimated_distance, # Pass system_config
                                        planet_class="M")
                        found_hab = True
                        self.planets.append(planet)
                        continue

                # Guaranteed last-resort fallback, mirroring HABITABLE_WORLD's own
                # end-of-loop guarantee: if a belt is still required and hasn't
                # happened yet by the trailing slot(s), force it there regardless of
                # the soft hz/last_asteroid avoidance below (validate_system corrects
                # any resulting spacing/overlap afterward). When HABITABLE_WORLD is
                # ALSO still pending, that claims the very last slot for itself (its
                # own unconditional `continue` above would otherwise steal this one
                # out from under the belt), so the belt's fallback claims the slot
                # just before it instead, guaranteeing each requirement its own slot.
                hab_still_pending = self.system_config.HABITABLE_WORLD is True and not found_hab
                belt_fallback_index = (system_objects - 2) if hab_still_pending else (system_objects - 1)
                force_belt = self.system_config.ASTEROID_BELT is True and not found_belt and i >= belt_fallback_index

                if self.system_config.ASTEROID_BELT is not False and (
                    force_belt or (not last_asteroid and not hz and (random.random() < program_constants.ASTEROID_BELT_PROBABILITY or i == belt_index))
                ):
                    min_distance = estimated_distance
                    max_distance = estimated_distance * random.uniform(program_constants.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MIN, program_constants.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MAX)
                    self.planets.append(AsteroidBelt(self.system_config, estimated_distance, min_distance, max_distance)) # Pass system_config
                    found_belt = True
                else:
                    planet = Planet(self.system_config, self.star, self.star.habitable_zone, estimated_distance) # Pass system_config
                    if planet.planet_class == "M":
                        found_hab = True
                    self.planets.append(planet)

    def to_dict(self):
        """
        Returns a JSON-serializable dict of the entire generated object
        graph: this system's config, star (polymorphic -- includes nested
        primary/secondary if binary), every planet/belt (in orbital order,
        each recursively including its own moons), and the system-level
        flavor text (see TODO.md's Phase 0 fix, which moved this to
        generation time so it's a plain, idempotent read here).

        Deliberately omits `planet_count`/`belt_count`/`moon_count`/
        `hab_count`/`m_count` (bookkeeping recomputed on load via
        `count_objects`/`count_habitable`, never trusted from disk) and
        `stars`/`primary_star`/`secondary_star` (resolvable from `star`
        alone -- see `from_dict`).

        Returns:
            dict: `schema_version`, `system_config`, `star`, `is_binary`,
                 `planets`, `system_flavor_text`.
        """
        return {
            "schema_version": SERIALIZATION_SCHEMA_VERSION,
            "system_config": self.system_config.to_dict(),
            "star": self.star.to_dict(),
            "is_binary": isinstance(self.star, BinaryStarProxy),
            "planets": [obj.to_dict() for obj in self.planets],
            "system_flavor_text": self.system_flavor_text,
        }

    @classmethod
    def from_dict(cls, data):
        """
        Reconstructs a `StarSystem` from a dict in the shape `to_dict()`
        produces, without re-running any generation (`__init__` is bypassed
        via `object.__new__`) -- this is a pure, faithful replay of
        already-decided data, not a new roll. There is deliberately no
        seed-based replay anywhere in this design: generation mixes the
        unseedable `secrets` module with the seedable `random` module (see
        `spaceSector.py`'s own module docstring), so a seed alone could
        never reproduce a system -- the actual decided values are stored
        and read back directly instead.

        This is the single place that resolves both shared back-references
        once and re-attaches the same instances everywhere: `system_config`
        is built once and threaded into every child; `star` is built once
        (via a `is_binary` discriminator choosing `Star.from_dict` vs.
        `BinaryStarProxy.from_dict`) and threaded into every top-level
        `Planet`/`AsteroidBelt`. For a binary system, the secondary star is
        deliberately reattached to the *same* shared `system_config` as
        everything else, collapsing the generation-time-only asymmetry
        where `StarSystem.__init__` gives the secondary its own deep-copied
        config (`LARGE_STAR` forced `False`) -- that flag is only ever
        consulted during `generate_star()`, so it has no meaning once a
        star already exists to be reloaded.

        Each item in `planets` is dispatched to `AsteroidBelt.from_dict` or
        `Planet.from_dict` based on its own `body_type` (`'a'` vs.
        `'t'`/`'g'`), the same discriminator `StarSystem` itself already
        uses to tell the two apart at runtime.

        Args:
            data (dict): A dict in the shape `to_dict()` produces.

        Returns:
            StarSystem: The reconstructed system.

        Raises:
            ValueError: If `data["schema_version"]` is newer than this code
                       understands.
        """
        schema_version = data.get("schema_version", SERIALIZATION_SCHEMA_VERSION)
        if schema_version > SERIALIZATION_SCHEMA_VERSION:
            raise ValueError(
                f"StarSystem.from_dict: schema_version {schema_version} is newer "
                f"than this code understands (max {SERIALIZATION_SCHEMA_VERSION})."
            )

        system_config = SystemConfig.from_dict(data["system_config"])

        if data["is_binary"]:
            star = BinaryStarProxy.from_dict(data["star"], system_config)
        else:
            star = Star.from_dict(data["star"], system_config)

        system = object.__new__(cls)
        system.system_config = system_config
        system.star = star

        if isinstance(star, BinaryStarProxy):
            system.primary_star = star._primary
            system.secondary_star = star._secondary
            system.stars = [star._primary, star._secondary]
        else:
            system.primary_star = star
            system.stars = [star]

        system.planets = []
        for obj_data in data["planets"]:
            if obj_data.get("body_type") == "a":
                obj = AsteroidBelt.from_dict(obj_data, system_config)
            else:
                obj = Planet.from_dict(obj_data, star, system_config)
            system.planets.append(obj)

        system.system_flavor_text = data.get("system_flavor_text")

        system.planet_count, system.belt_count, system.moon_count = system.count_objects()
        system.hab_count, system.m_count = system.count_habitable()

        return system

    def generate_slot_object(self, slot_spec, estimated_distance):
        """
        Builds the celestial object explicitly requested for one orbital slot
        by `system_config.SLOTS`.

        Args:
            slot_spec (dict): The slot specification, with a required "type"
                              key ("planet" or "asteroid_belt") and optional
                              "planet_class"/"moons" keys (see `SystemConfig.SLOTS`).
            estimated_distance (float): The orbital distance (in AU) computed
                                        for this slot by the generation loop.

        Returns:
            Planet or AsteroidBelt: The generated object for this slot.

        Raises:
            ValueError: If `slot_spec["type"]` is not "planet" or "asteroid_belt".
        """
        slot_type = slot_spec.get("type", "planet")

        if slot_type == "asteroid_belt":
            min_distance = estimated_distance
            max_distance = estimated_distance * random.uniform(program_constants.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MIN, program_constants.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MAX)
            return AsteroidBelt(self.system_config, estimated_distance, min_distance, max_distance)

        if slot_type == "planet":
            planet_class = slot_spec.get("planet_class")
            distance = self.calculate_distance_for_class(planet_class, estimated_distance)
            return Planet(self.system_config, self.star, self.star.habitable_zone, distance,
                         planet_class=planet_class, moon_count=slot_spec.get("moons"))

        raise ValueError(f"Invalid slot type '{slot_type}'; expected 'planet' or 'asteroid_belt'.")

    def calculate_distance_for_class(self, planet_class, estimated_distance):
        """
        Adjusts an orbital distance so that it falls in a zone (hot, cold, or
        ecosphere) that actually supports a requested `planet_class`.

        The generation loop picks `estimated_distance` before knowing what
        the slot will contain, so a user-requested class (e.g. "M", which
        only exists in the ecosphere) may not be valid at that distance. This
        nudges the distance into a supporting zone, the same way the existing
        forced-habitable-world logic snaps a planet's distance into the
        habitable zone, so an explicit `planet_class` request doesn't fail
        with an "Invalid planet class for this zone" error just because of
        where its slot happened to land in the orbit sequence.

        Args:
            planet_class (str or None): The requested planet class, or None
                                        if the slot doesn't specify one.
            estimated_distance (float): The orbital distance (in AU) computed
                                        for this slot by the generation loop.

        Returns:
            float: `estimated_distance`, or an adjusted distance (in AU) that
                  falls within a zone supporting `planet_class`.
        """
        class_data = program_constants.PLANET_CLASSES.get(planet_class)
        if class_data is None:
            return estimated_distance

        inner, outer = self.star.habitable_zone
        if estimated_distance < inner:
            zone = 'h'
        elif estimated_distance > outer:
            zone = 'c'
        else:
            zone = 'e'

        if class_data.get(zone):
            return estimated_distance

        if class_data.get('e'):
            return self._distance_within_zone_with_margin(inner, outer)
        if class_data.get('h'):
            return random.uniform(inner * 0.05, inner * 0.95)
        if class_data.get('c'):
            return outer * random.uniform(1.05, 3.0)

        # No zone supports this class; leave the distance as-is and let
        # planetPhysics raise its usual, clearer validation error.
        return estimated_distance

    def _distance_within_zone_with_margin(self, inner, outer):
        """
        Draws a distance uniformly within `[inner, outer]`, leaving a
        safety margin at both edges against `validate_system`'s own
        `MIN_ASTEROID_BELT_SEPARATION` nudge.

        A planet deliberately placed right at a zone boundary (by the
        forced-habitable-world logic, via `_forced_habitable_distance`, or
        an explicit slot request via `calculate_distance_for_class`) could
        otherwise get pushed across that boundary later, when
        `validate_system` enforces minimum separation from a neighboring
        asteroid belt -- which `planetPhysics.reconcile_zone_and_class`
        would then have to reclassify away from the very class this draw
        was placing it as. Reserving this margin up front avoids that in
        the common case; margin is capped at a quarter of the zone's own
        width so a pathologically narrow zone still gets a valid
        (non-empty) range to draw from.

        Args:
            inner (float): Zone's inner bound, in AU.
            outer (float): Zone's outer bound, in AU.

        Returns:
            float: A distance in AU, safely inside `[inner, outer]`.
        """
        margin = min(program_constants.MIN_ASTEROID_BELT_SEPARATION, (outer - inner) / 4)
        return random.uniform(inner + margin, outer - margin)

    def _forced_habitable_distance(self):
        """
        Draws a distance for a planet the generation loop is forcing into
        the habitable zone (`HABITABLE_WORLD=True`), via
        `_distance_within_zone_with_margin`.

        Returns:
            float: A distance in AU, safely inside the star's
                  `habitable_zone`.
        """
        inner, outer = self.star.habitable_zone
        return self._distance_within_zone_with_margin(inner, outer)

    def count_objects(self):
        """
        Counts the number of planets, asteroid belts, and moons in the system.

        This method iterates through the list of celestial objects in the system
        and tallies the number of each type. It distinguishes between planets and
        asteroid belts, and also counts the total number of moons orbiting the
        planets. This information is used to provide a summary of the system's
        composition in the `__str__` method.

        Returns:
            tuple: A tuple containing the total number of planets, asteroid belts,
                   and moons in the system, in that order.
        """
        planet_counter, moon_counter, belt_counter = 0, 0, 0
        for planet in self.planets:
            if planet.body_type == 'a':
                belt_counter += 1
            else:
                planet_counter += 1
                moon_counter += len(planet.moons)
        return planet_counter, belt_counter, moon_counter

    def count_habitable(self):
        """
        Counts the number of potentially habitable worlds in the system.

        This method iterates through all planets and their moons, checking their
        classification against `program_constants.HABITABLE_PLANET_CLASSES` to
        determine if they are potentially habitable, and also keeps a separate
        count of Class M worlds, which are considered the most Earth-like. This
        is used for the system summary output.

        Returns:
            tuple: A tuple containing the total number of potentially habitable
                   worlds and the total number of Class M worlds, in that order.
        """
        habitable_classes = program_constants.HABITABLE_PLANET_CLASSES
        hab_count, m_count = 0, 0
        for planet in self.planets:
            if planet.body_type != 'a':
                if planet.planet_class in habitable_classes:
                    hab_count += 1
                if planet.planet_class == "M":
                    m_count += 1
                for moon in planet.moons:
                    if moon.planet_class in habitable_classes:
                        hab_count += 1
                    if moon.planet_class == "M":
                        m_count += 1
        return hab_count, m_count

    def estimate_num_objects(self):
        """
        Estimates the number of objects in a star system from real
        protoplanetary-disk physics, rather than an arbitrary curve fit to
        the star's mass alone.

        `_estimate_max_objects_from_disk_physics` derives a ceiling from
        the star's own disk-mass budget (scaled from the Sun's own
        Minimum Mass Solar Nebula -- see that method's docstring) rather
        than a logarithmic fudge factor. The final number is a random
        value between the minimum and that ceiling, or the ceiling/minimum
        itself if `MAX_PLANETS` is True/False -- this part of the
        contract is unchanged from before.

        `NUM_ORBITS`, when set, bypasses this estimate entirely and is
        returned as-is. `PLANETS` set to False forces zero objects; set to
        True, it raises the minimum considered to 1.

        Returns:
            int: The estimated number of objects to be generated in the system.
        """
        if self.system_config.PLANETS is False:
            return 0

        if self.system_config.NUM_ORBITS is not None:
            return self.system_config.NUM_ORBITS

        max_objects = self._estimate_max_objects_from_disk_physics()

        min_objects = 1 if self.system_config.PLANETS is True else 0
        max_objects = max(max_objects, min_objects)

        if self.system_config.MAX_PLANETS is True:
            return max_objects
        if self.system_config.MAX_PLANETS is False:
            return min_objects
        return random.randint(min_objects, max_objects)

    def _estimate_max_objects_from_disk_physics(self):
        """
        Derives a physically-grounded ceiling on planet/belt count from
        this system's own protoplanetary disk, instead of an arbitrary
        curve fit to stellar mass.

        Real planet formation gives a natural way to turn "how much solid
        material did this star's disk have" into "how many planets could
        that have made": embryos grow by clearing their own feeding zone
        until they reach their oligarchic-growth *isolation mass*
        (`utils.isolation_mass_kg` -- Lissauer 1993; Kokubo & Ida 2000,
        2002), then space apart from their neighbors by the same *mutual*
        Hill radius `_mutual_min_distance_au` already enforces during
        placement (`program_constants.MUTUAL_HILL_RADII_SEPARATION`) --
        so this walk and that spacing rule are provably consistent with
        each other, unlike the old log-mass formula, which had no
        relationship to it at all.

        The walk starts at the same inner-edge distance
        `_generate_planets` itself seeds its first slot at
        (`program_constants.INITIAL_PLANET_DISTANCE_FACTOR * star_factor`)
        and steps outward -- at each step, computing the local isolation
        mass from the disk's surface density there (scaled for this star
        via `utils.disk_surface_density_scale`, boosted beyond the snow
        line via `utils.snow_line_au`), then advancing by that embryo's
        own mutual-Hill-radius feeding zone -- until reaching the disk's
        outer edge, at
        `program_constants.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER` times
        the snow line (real disks are truncated far short of
        `self.star.system_perimeter`'s galactic-tidal scale; see that
        constant's docstring), or `ABSOLUTE_MAX_SYSTEM_OBJECTS` isolation-
        mass slots, whichever comes first.

        Not every one of those oligarchs survives as a final planet --
        real N-body integrations of the subsequent giant-impact phase
        (Chambers 2001) show most merge or get ejected -- so the raw slot
        count is scaled down by
        `program_constants.GIANT_IMPACT_SURVIVAL_FRACTION` before being
        returned.

        Returns:
            int: The physically-derived ceiling on this system's object
                count (`max_objects`, as `estimate_num_objects` already
                names it).
        """
        star_mass_kg = self.star.mass
        star_factor = star_mass_kg / physical_constants.SOLAR_MASS_TO_KG
        snow_line = snow_line_au(self.star.luminosity)
        density_scale = disk_surface_density_scale(star_mass_kg)
        outer_edge_au = min(
            program_constants.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER * snow_line,
            self.star.system_perimeter,
        )

        distance_au = program_constants.INITIAL_PLANET_DISTANCE_FACTOR * star_factor
        oligarch_count = 0
        while (
            distance_au < outer_edge_au
            and oligarch_count < program_constants.ABSOLUTE_MAX_SYSTEM_OBJECTS
        ):
            surface_density = mmsn_surface_density_gcm2(distance_au, snow_line, density_scale)
            embryo_mass_kg = isolation_mass_kg(distance_au, surface_density, star_mass_kg)

            # Two neighboring oligarchs of roughly the same isolation mass
            # (2 * embryo_mass_kg is the pair's combined mass) -- the same
            # kappa/clamp `_mutual_min_distance_au` uses, so this step and
            # that later spacing check agree on how far apart is "enough."
            kappa = program_constants.MUTUAL_HILL_RADII_SEPARATION * (
                (2 * embryo_mass_kg) / (3 * star_mass_kg)
            ) ** (1 / 3)
            kappa = min(kappa, 1.8)

            distance_au *= (1 + kappa / 2) / (1 - kappa / 2)
            oligarch_count += 1

        max_objects = math.ceil(oligarch_count * program_constants.GIANT_IMPACT_SURVIVAL_FRACTION)
        return min(max_objects, program_constants.ABSOLUTE_MAX_SYSTEM_OBJECTS)

    def validate_system(self):
        """
        Validates and adjusts the distances of stellar objects to prevent orbital overlap.

        This method iterates through the generated planets and other celestial
        objects, checking for any orbital overlaps. If two objects are too close
        to each other, it adjusts their orbits to ensure a safe distance is
        maintained. This is crucial for creating a physically plausible and stable
        star system.

        The method accounts for the different types of objects, such as planets
        and asteroid belts, and applies appropriate corrections to ensure that
        their orbits are not just non-overlapping, but also realistically spaced.
        Two adjacent planets' minimum separation is enforced via their
        *mutual* Hill radius (`_mutual_min_separation_au`), not either
        one's own individual Hill radius alone -- see
        `program_constants.MUTUAL_HILL_RADII_SEPARATION`'s docstring for
        why. A belt has no mass/Hill-radius concept of its own, so any
        correction involving one still falls back to the fixed
        `MIN_ASTEROID_BELT_SEPARATION` or the single real planet's own
        `min_orbit_distance`, whichever case applies.
        If an adjustment is made to a planet (as opposed to an asteroid belt,
        which carries no class/climate of its own), `_reconcile_moved_planet`
        re-derives its zone from the corrected distance and, if its
        already-rolled class is no longer valid there, regenerates the class
        and everything derived from it -- composition, radius, mass,
        density, atmosphere, period, gravity, and orbital motion -- along
        with the same treatment for each of its moons. Before this, a
        pushed-out planet could keep reporting a class the corrected
        distance no longer physically supports at all (e.g. an "Earth-like"
        Class M at a Class-M-invalid distance, sometimes thousands of AU
        past the actual habitable zone -- see docs/TODO.md's now-resolved
        "validate_system can strand a planet outside its own zone" entry).
        """
        if len(self.planets) < 2:
            return

        for i in range(1, len(self.planets)):
            planet = self.planets[i]
            last_planet = self.planets[i - 1]

            if last_planet.body_type == 'a':
                distance_to_last = planet.distance - last_planet.upper_limit
            else:
                distance_to_last = planet.distance - last_planet.distance

            if distance_to_last < 0:
                additional_correction = abs(distance_to_last) + last_planet.distance
            else:
                additional_correction = 0

            if planet.body_type == 'a':
                if last_planet.body_type == 'a':
                    if distance_to_last < program_constants.MIN_ASTEROID_BELT_SEPARATION:
                        planet.distance += program_constants.MIN_ASTEROID_BELT_SEPARATION + additional_correction
                        planet.upper_limit += program_constants.MIN_ASTEROID_BELT_SEPARATION + additional_correction
                        planet.lower_limit += program_constants.MIN_ASTEROID_BELT_SEPARATION + additional_correction
                elif distance_to_last < last_planet.min_orbit_distance:
                    planet.distance += last_planet.min_orbit_distance + additional_correction
                    planet.upper_limit += last_planet.min_orbit_distance + additional_correction
                    planet.lower_limit += last_planet.min_orbit_distance + additional_correction
            else:
                if last_planet.body_type == 'a':
                    if distance_to_last < program_constants.MIN_ASTEROID_BELT_SEPARATION:
                        planet.distance += program_constants.MIN_ASTEROID_BELT_SEPARATION + additional_correction
                        self._reconcile_moved_planet(planet)
                else:
                    # Bounded, not a single shot: reclassifying `planet`
                    # (inside `_reconcile_moved_planet`) can change its own
                    # `mass`, which the mutual Hill radius depends on --
                    # re-deriving the requirement against the *new* mass
                    # can call for a further push, so retry until a push
                    # doesn't trigger another reclassification (in
                    # practice at most 2 iterations; the cap is just a
                    # guard against a pathological cycle).
                    for _ in range(3):
                        min_distance = self._mutual_min_distance_au(planet, last_planet)
                        if planet.distance >= min_distance:
                            break
                        planet.distance = min_distance
                        if not self._reconcile_moved_planet(planet):
                            break

    def _mutual_min_distance_au(self, planet, last_planet):
        """
        The closest `planet` (the later/farther of the pair) can stably
        sit to `last_planet` (fixed -- it's already been placed and
        won't move again this pass), via their *mutual* Hill radius
        (`utils.mutual_hill_radius_m`) rather than either one's own
        individual Hill radius alone -- see
        `program_constants.MUTUAL_HILL_RADII_SEPARATION`'s docstring for
        the stability-literature basis.

        Solved in closed form for `planet`'s own distance rather than
        evaluated once at the current (pre-correction) positions and
        added on top: the mutual Hill radius depends on the *average* of
        the two distances, so a naive "compute now, then push `planet`
        out by that much" approximation understates the requirement --
        moving `planet` out raises the average, which raises the
        requirement further, and `MUTUAL_HILL_RADII_SEPARATION`'s margin
        (10x) is large enough that this understatement is measurable, not
        just a rounding-level slip (confirmed by a real
        `assert_no_orbital_overlap` failure before this closed form
        replaced the naive version).

        Derivation: let d = `last_planet.distance` (fixed this pass), x =
        `planet`'s target distance, and kappa =
        `MUTUAL_HILL_RADII_SEPARATION * ((m_planet + m_last) / (3 *
        M_star)) ** (1/3)`. The stability requirement `x - d >= kappa *
        (x + d) / 2` rearranges to `x >= d * (1 + kappa/2) / (1 -
        kappa/2)`.

        Args:
            planet (Planet): The planet whose distance may need raising.
            last_planet (Planet): The fixed reference planet, closer to
                                  the star.

        Returns:
            float: The minimum distance (AU) `planet` can sit at,
                  measured from the star -- not a gap.
        """
        kappa = program_constants.MUTUAL_HILL_RADII_SEPARATION * (
            (planet.mass + last_planet.mass) / (3 * self.star.mass)
        ) ** (1 / 3)
        # A pair whose combined mass is a large enough fraction of the
        # star's own that kappa approaches/exceeds 2 has no finite stable
        # separation under this linear model at all -- clamp well short
        # of that so the formula below always returns a large-but-finite,
        # rather than negative or infinite, distance.
        kappa = min(kappa, 1.8)
        return last_planet.distance * (1 + kappa / 2) / (1 - kappa / 2)

    def _reconcile_moved_planet(self, planet):
        """
        Call right after `validate_system` moves a top-level planet's
        `distance` to resolve an orbital overlap.

        A class is only a valid description of a body at the distance it
        actually ends up at -- pushing `distance` out to fix spacing can
        carry a planet into a zone ('h'/'e'/'c') its already-rolled class
        was never valid in (e.g. an "Earth-like" Class M shoved out past
        the habitable zone into the cold zone), which would otherwise
        report a class the corrected position no longer supports at all.
        `planetPhysics.reconcile_zone_and_class` re-derives the zone from
        the new `distance` and, only if the existing class no longer fits
        there, regenerates the class and everything derived from it
        (radius/mass/density/atmosphere/period/gravity/orbital motion) --
        physically, this is the same thing real orbital migration does:
        change the body's actual final conditions, not just its assumed
        ones.

        Every moon is reconciled the same way: a moon's own zone is always
        its parent's (see `planetPhysics.generate_moons`' `zone_override`),
        and its climate depends on the parent's distance (its
        `distance_override`), so a moon whose parent just got reclassified
        needs the identical treatment even though the moon's own orbit
        around its parent never changed. A reclassified parent can also
        come out with a different `mass` than before -- and Kepler's third
        law makes a moon's own `period` (and the position/orbital speed/
        rotation-period-if-tidally-locked derived from it) depend on the
        *parent's* mass, not just the moon's own (unchanged) distance
        around it, so every moon's period needs refreshing whenever the
        parent was reclassified, even a moon that didn't need
        reclassifying itself -- via `generate_orbital_motion_properties`
        rather than only `update_orbital_position`, since a *moon's*
        `rotation_period_hours` can itself be period-derived (a tidally
        locked moon's day equals its orbit) in a way an ordinary planet's
        never is, and a period change can flip whether that still holds.

        Returns:
            bool: True if `planet` itself was reclassified (and so may
                 have come out with a different `mass` than before --
                 relevant to `validate_system`'s mutual-Hill-radius
                 spacing check, which depends on both planets' masses);
                 False if its existing class remained valid.
        """
        reclassified = planetPhysics.reconcile_zone_and_class(planet, planet.star.mass)
        if not reclassified:
            planetPhysics.calculate_atmospheric_conditions(planet)
            planet.period = planetPhysics.calculate_orbital_period_years(planet.distance, planet.star.mass)
            planetPhysics.update_orbital_position(planet)

        for moon in planet.moons:
            moon_reclassified = planetPhysics.reconcile_zone_and_class(
                moon, planet.mass, distance_override=planet.distance
            )
            if not moon_reclassified:
                planetPhysics.calculate_atmospheric_conditions(moon, planet.distance)
                if reclassified:
                    moon.period = planetPhysics.calculate_orbital_period_years(moon.distance, planet.mass)
                    planetPhysics.generate_orbital_motion_properties(moon, planet.mass)

        return reclassified

    def __str__(self):
        """
        Generates a string output for the system data, including a summary and
        details for each celestial body.

        This method compiles a comprehensive summary of the entire star system,
        including details about the central star and each of its orbiting
        objects. The output is formatted as a human-readable string, suitable for
        display in a console or for writing to a file.

        The method provides a high-level overview of the system, including the
        total number of planets, asteroid belts, and moons, as well as the number
        of potentially habitable worlds. It then lists each celestial body in
        order, providing a detailed description of its properties. The final
        output may also include a category tag for wiki-based systems.

        Returns:
            str: A formatted string representing the entire star system.
        """
        all_output_parts = []

        # Add level 1 header for the system name
        if self.system_config.MARKDOWN:
            all_output_parts.append(f"# {self.star.name}\n\n")
        else:
            all_output_parts.append(f"= {self.star.name} =\n\n")

        # Generate system summary sentences (needed before star details for binary systems)
        system_summary_sentences = []
        segments = []
        if 0 < self.planet_count != self.m_count and self.hab_count != self.planet_count:
            segments.append(f"{self.planet_count} planet{'s' if self.planet_count > 1 else ''}")
        if self.belt_count > 0:
            segments.append(f"{self.belt_count} asteroid belt{'s' if self.belt_count > 1 else ''}")
        if self.moon_count > 0:
            segments.append(f"{self.moon_count} moon{'s' if self.moon_count > 1 else ''}")

        if segments:
            system_string = "This system contains " + ", ".join(segments)
            if len(segments) == 2:
                system_string = system_string.replace(f", {segments[-1]}", f" and {segments[-1]}")
            elif len(segments) > 2:
                system_string = system_string.replace(f", {segments[-1]}", f", and {segments[-1]}")
            system_string += "."
            system_summary_sentences.append(system_string)
        elif self.planet_count == 0 and self.belt_count == 0 and self.moon_count == 0:
            # No segment was built and the system is genuinely empty (as opposed to having
            # its planet count omitted above to avoid redundancy with the habitability
            # sentence that follows).
            system_summary_sentences.append("There are no stellar objects in this system.")

        if self.m_count == 1:
            m_string = "1 of which is class M"
        elif self.m_count > 1:
            m_string = f"{self.m_count} of which are class M"
        else:
            m_string = "none of which are class M"

        if self.hab_count == 1 and self.m_count < 1:
            system_summary_sentences.append("There is 1 potentially habitable world in the system.")
        elif self.hab_count == 1 and self.m_count == 1:
            system_summary_sentences.append("There is 1 class M world in the system.")
        elif self.hab_count > 1 and self.hab_count == self.m_count:
            system_summary_sentences.append(f"There are {self.hab_count} class M worlds in the system.")
        elif self.hab_count > 1:
            system_summary_sentences.append(f"There are {self.hab_count} potentially habitable worlds ({m_string}).")
        else:
            system_summary_sentences.append("There are no potentially habitable worlds in this system.")

        perimeter_ly = self.star.system_perimeter * physical_constants.AU_TO_LY
        heliosphere_ly = self.star.heliosphere_radius * physical_constants.AU_TO_LY
        if heliosphere_ly < 0.1:
            heliosphere_text = f"{self.star.heliosphere_radius:.4f} AU"
        else:
            heliosphere_text = f"{heliosphere_ly:.4f} light-years"

        system_summary_sentences.append(
            f"The star's stellar wind creates a bubble, known as the heliosphere, which extends out to approximately {heliosphere_text}.")
        system_summary_sentences.append(
            f"Beyond this, the star's gravitational influence extends out to a distance of {perimeter_ly:.2f} light-years, marking the ultimate edge of the system.")

        combined_system_summary_paragraph = to_paragraph(system_summary_sentences)


        if isinstance(self.star, BinaryStarProxy):
            # Get combined binary system data and age from the proxy
            # BinaryStarProxy.to_paragraph_list() returns [data_block, age_sentence]
            binary_proxy_paragraphs = self.star.to_paragraph_list()

            # 1. Append the combined binary system data block
            all_output_parts.append(binary_proxy_paragraphs[0]) # Binary System Data
            all_output_parts.append('\n\n')

            # 2. Append the system summary
            all_output_parts.append(combined_system_summary_paragraph)
            all_output_parts.append('\n\n')

            # 3. Append the combined age sentence for the binary system
            all_output_parts.append(binary_proxy_paragraphs[1]) # Binary System Age sentence

            # Flavor text was already decided at generation time; this is a pure
            # read so repeated renders don't re-roll or double-count it.
            if self.system_flavor_text:
                all_output_parts.append(f"\n\nSensors show {self.system_flavor_text}")

            # 4. Append individual star details for primary and secondary stars
            for star_obj in self.stars: # self.stars contains primary and secondary Star objects
                all_output_parts.append('\n\n') # Blank line after age sentence
                header_level = '===' if not self.system_config.MARKDOWN else '###'
                all_output_parts.append(f"{header_level} {star_obj.name} {header_level if not self.system_config.MARKDOWN else ''}".rstrip())
                all_output_parts.append('\n') # Add a newline after the header

                # Each individual star's to_paragraph_list() returns [data_block, age_sentence]
                individual_star_details = star_obj.to_paragraph_list()
                all_output_parts.append(individual_star_details[0]) # Individual Star Data block
                all_output_parts.append('\n\n') # Blank line after data block
                all_output_parts.append(individual_star_details[1]) # Individual Star Age sentence

        else:
            # For single star:
            # Get single star data and age from the star object
            # Star.to_paragraph_list() returns [data_block, age_sentence]
            single_star_paragraphs = self.star.to_paragraph_list()

            # 1. Append the single star data block
            all_output_parts.append(single_star_paragraphs[0])
            all_output_parts.append('\n\n')
            # 2. Append the single star age sentence
            all_output_parts.append(single_star_paragraphs[1])
            all_output_parts.append('\n\n')
            # 3. Append the system summary
            all_output_parts.append(combined_system_summary_paragraph)

            # Flavor text was already decided at generation time; this is a pure
            # read so repeated renders don't re-roll or double-count it.
            if self.system_flavor_text:
                all_output_parts.append(f"\n\nSensors show {self.system_flavor_text}")

        # Add planet/belt paragraphs, each separated by a double newline from the previous.
        if self.planets:
            for planet in self.planets:
                all_output_parts.append('\n\n' + '\n\n'.join(planet.to_paragraph_list()))

        # Add category tag if not markdown
        if not self.system_config.MARKDOWN:
            all_output_parts.append('\n\n' + '[[Category:Star Systems]]')

        return ''.join(all_output_parts)