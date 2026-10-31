# planetgen/generation/system.py

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

import dataclasses
import math

from planetgen.generation.belt import AsteroidBelt
from planetgen.names.bodies import CLOSE_PAIR_LETTERS, generate_star_word, name_bodies, wide_pair_first_word
from planetgen.generation.comet import Comet
from planetgen.generation.config import SystemConfig
from planetgen.generation.binary import BinaryStarProxy
from planetgen.generation import (
    cloud as cloudContext, life as planetLife, member_record, multistar, place_terms, prevalence, stored_trees, validation,
)
from planetgen.physics import stellar_evolution as evolution
from planetgen.physics import habitability, insolation, interior_heat, magnetism, mass_loss, ocean_tides, pre_main_sequence, radiation, rings, wolf_rayet
from planetgen.physics import planets as planetPhysics
from planetgen import tuning
from planetgen.util import draw
from planetgen.util import log
from planetgen.generation import habitable_request, pulsar_planets, remnant_orbits
from planetgen.generation.planet import DEFER_MOONS, Planet
from planetgen.generation.star import Star, adjust_pair_age_for_planets
from planetgen.physics import bulk_draw, clearance_ledger, planet_architecture, planet_zones
from planetgen.physics.system_tree import OrbitNode, TreeMember
from planetgen.physics.orbits import calculate_reflex_offset
from planetgen.generation.wide_binary import WideBinaryPair
from planetgen.physics import unit_registry as units

_km_to_au = units.converter("km", "au")

HELIOSPHERE_COMPRESSION_NOTED = 0.95
"""float: The system text mentions the surrounding cloud's squeeze only when
it leaves the heliosphere below this share of its open-space size."""

# Tracks the shape of `StarSystem.to_dict()`'s output (the serialized
# object-graph -- see TODO.md's Phase 1), independent of
# `planetgen.db.store.SCHEMA_VERSION`, which tracks the *database's* DDL
# structure instead. A JSON export carries this number so `from_dict` can
# raise a clear error if a file's shape is newer than the code understands,
# rather than a confusing `KeyError`/`AttributeError` partway through
# reconstruction.
#
# Bumped 1 -> 2 to add S-type (wide) binary support: `binary_type`
# (None | "close" | "wide", now the authoritative discriminator --
# `is_binary` is kept for backward compatibility and now means "this system
# has two stars", true for either configuration), `secondary_star` and
# `wide_binary` (both null unless `binary_type == "wide"`), and
# `secondary_planets` (the wide binary secondary's own independently
# generated planet list; empty otherwise). A pre-v2 save has no
# `binary_type` key at all -- `from_dict` treats that as "close" whenever
# its `is_binary` is true, since P-type was the only binary configuration
# that existed before this version.
#
# Bumped 2 -> 3 to add native star-bound comets (`cometData.Comet` -- see
# `_generate_comets`): `comets` and `secondary_comets`, the same
# "primary's own list" / "wide binary secondary's own independently
# generated list" split `planets`/`secondary_planets` already uses, kept
# as their own top-level keys rather than folded into `planets` (a
# `Comet`'s `distance_au` is a continuously-varying current position, not
# a fixed orbital-slot `distance` the way a `Planet`'s/`AsteroidBelt`'s
# is -- see `_generate_comets`'s own docstring). A pre-v3 save has neither
# key -- `from_dict` treats a missing key as `[]`, the same fallback
# `secondary_planets` already gets for a pre-wide-binary save.
#
# Bumped 3 -> 4 to add `orbit_nodes` (GEN.311): the stored tree of orbits the stars hang in, a list of dicts, one per node
# (`physics.system_tree.OrbitNode`); `[]` for a single star or a system that was never saved. A pre-v4 save has no key and
# reads as no nodes.
#
# Bumped 4 -> 5 to add `runaway_origin` (GEN.329): how a runaway star came to move so fast, `"supernova"` or `"dynamical"`
# (`None` for a system that is not a runaway). A pre-v5 save has no key and reads as `None`.
#
# Bumped 5 -> 6 to add `member_stars` and `star_slots` (GEN.630): the stars beyond the pair of a system of more than two stars
# and where each star hangs in the tree (`None`/`[]` for a single star and a pair). A pre-v6 save has neither key and reads
# as no members.
#
# Bumped 6 -> 7 to add `member_planets` (GEN.332): the planets and belts of each of `member_stars`, a list of lists in the
# members' order. A pre-v7 save has no key and reads as no planets of the members.
#
# Bumped 7 -> 8 to add `dark_members` (GEN.1013): the brown dwarfs, neutron stars and black holes that hang in the tree but
# are not stars of `stars`, a list of dicts (`stored_trees.DarkMember`). A pre-v8 save has no key and reads as none.
SERIALIZATION_SCHEMA_VERSION = 8


def _exceeds_orbit_ceiling(obj, orbit_ceiling_au):
    """
    True if `obj`'s own outer edge -- `upper_limit` for an `AsteroidBelt`,
    or `distance + min_orbit_distance` (5x its own Hill radius) for a
    `Planet` -- extends past `orbit_ceiling_au`.

    Used by `StarSystem._generate_planets` to enforce that a newly placed
    body's own sphere of influence, not just its central distance, stays
    inside a star's stability ceiling (see `StarSystem._orbit_ceiling_au`)
    -- mirroring the same `distance + min_orbit_distance` idiom already
    used to space successive planets around one star (see
    `physics/planets.py`'s `min_orbit_distance`, `5 * hill_radius`). Without
    this, a candidate's *center* could sit just inside a tight,
    companion-driven S-type limit while its own Hill sphere still reached
    past it.

    Args:
        obj (Planet or AsteroidBelt): The most recently placed body.
        orbit_ceiling_au (float): The boundary to check against, in AU.

    Returns:
        bool: True if `obj` extends past `orbit_ceiling_au`.
    """
    return validation.exceeds_orbit_ceiling(obj, orbit_ceiling_au)


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
        comets (list): A list of native `cometData.Comet` objects bound to
            this star -- kept separate from `planets` since a comet's
            `distance_au` is a continuously-varying current position, not
            a fixed orbital-slot `distance` (see `_generate_comets`).
        planet_count (int): The total number of planets in the system.
        belt_count (int): The total number of asteroid belts in the system.
        comet_count (int): The total number of comets in the system.
        moon_count (int): The total number of moons in the system.
        hab_count (int): The total number of habitable worlds (overall PHI-6 Blue or Green, `habitability.is_habitable`).
        habitable_class_count (int): The planets and moons of a class that can carry life (`count_habitable_class_worlds`).
        habitable_request (str or None): How the system answers `system_config.HABITABLE_WORLD` (`habitable_request_status`).
        m_count (int): The total number of Class M worlds.
        system_flavor_count (int): Tracks the total number of flavor texts added across the system.
        system_flavor_text (str or None): The system-level flavor text sentence,
            decided once at generation time (or None if the roll didn't select one).
    """

    orbit_nodes = ()
    """tuple of `physics.system_tree.OrbitNode`: The tree of orbits the stars hang in (GEN.311). Empty until the system is
    saved (a generated binary's one node is made from its pair then, `db.orbit_trees.binary_nodes`) or for a single star;
    a loaded system has the tree that was stored. A system built from a tree (drawn, GEN.631, or a stored row's,
    `member_tree`) has its nodes from the start, a drawn tree's at `stored_trees.GENERATED_EPOCH_UNIX` until the save."""

    _ledger = None
    """Ledger or None: The clearance ledger of the stars as they are now (`ledger`); reset whenever the companions change."""

    placed_star = None
    """Star: The star the system was built around, the one a pre-placed row or a sector's draw made first (the tree's member 0).
    It is `stars[0]` unless a companion outweighs it, when the heavier takes the primary's place (GEN.631) and the placed star
    is a secondary or a member; what the place says of the star (its natal nebula, GEN.530) is said of this one (GEN.710).
    Set when the system is built, None for one loaded from the database."""

    member_stars = ()
    """tuple of Star: The stars beyond the pair (`stars[2:]`) of a system of more than two stars (GEN.630); empty for a
    single star and a pair. Each is stored as a `member` row of
    `stars`; its planets are `member_planets` (GEN.332)."""

    member_planets = ()
    """tuple of list: The planets and belts of each of `member_stars`, in that order (GEN.332): bodies in the S-type window
    round the member (`physics.planet_zones`), drawn like a wide pair's secondary's and never under a `SLOTS` list or a
    forced habitable world. Empty for a system with no members. Members have no comets yet."""

    star_slots = None
    """tuple or None: Where each star of `stars` hangs in `orbit_nodes`, as `(node_no, slot)` in `stars` order (slot 0 or
    1: the child of that node); `None` for a single star and for a pair on the one root node (the primary in slot 0, the
    secondary in slot 1). A system with `member_stars` or `dark_members` says it."""

    dark_members = ()
    """tuple of `stored_trees.DarkMember`: The members of the tree that are not stars of `stars` (GEN.1013): a brown dwarf, a
    neutron star, a black hole, hung in `orbit_nodes` with the mass and radius the tree drew, so the orbits, the node masses
    and the gravity are the drawn tree's. They have no planets and give no light here. Their `member_no` is the number in the
    plan, which `star_slots` does not use; the tree places them by `(parent_node_no, slot)`."""

    tree_built = False
    """bool: The system was built from the tree its scatter row stored and holds every member of it (GEN.315): the
    save deletes the row's stored members and nodes. False for a system generated whole and for a stored tree the
    build could not hold (a compact or brown-dwarf member)."""

    runaway_class = None
    """str or None: `"runaway"` or `"hypervelocity"` when the system moves
    unusually fast through the galaxy (`generate.flag_fast_stars`, schema
    v37); `None` for an ordinary one."""

    runaway_speed_kms = None
    """float or None: That speed, km/s; `None` unless `runaway_class`."""

    runaway_decided = False
    """bool: A pre-placed O or B star whose seed already decided whether it is a runaway (GEN.329): the fill's field roll
    (`run_sector.flag_fast_stars`) leaves it alone."""

    runaway_origin = None
    """str or None: How a `"runaway"` star came to move so fast (GEN.329): `"supernova"` (its companion exploded) or
    `"dynamical"` (an encounter in a young group), recorded in the flags of its `stars` row. `None` for any other system;
    a hypervelocity star has none (the Hills mechanism)."""

    runaway_direction = None
    """tuple or None: The unit vector, galactic axes, it moves along beyond
    the galaxy's rotation (`galaxy.system_position.peculiar_velocity_ms`);
    `None` unless `runaway_class`."""

    unmet_requirements = ()
    """list of str: What `system_config` forced (a habitable world, an
    asteroid belt) that every placement attempt failed to give this
    system, e.g. a habitable world around a hot O star; empty when nothing
    was missed. A single-system run refuses to save such a system
    (`generate.run_system`, GEN.49)."""

    def __init__(self, system_config: SystemConfig, galactic_center_dist_ly=None, compact_remnant=None,
                 primary_star_params=None, member_tree=None):
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
        `physics.planets`'s module docstring. Once all planets/moons exist,
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
                star's Hill-sphere calculation on the radius that holds half of the
                galaxy's disc stars (`galaxy.scales`) --
                the correct behavior for a system with no galaxy placement
                (e.g. `generation/run_sector.py`'s own standalone CLI).
            compact_remnant (BlackHole or NeutronStar, optional): A
                pre-built compact stellar remnant (see
                `generation/phenomena/compact_remnant.py`'s module docstring) to anchor this
                system in place of generating a normal `Star`. `None` (the
                default) generates an ordinary `Star` as usual. When given,
                the binary-generation step below is skipped entirely --
                this generator doesn't model a compact remnant anchoring a
                binary (e.g. a still-living stellar companion) -- so
                `self.binary_type` stays `None` and `self.star` is exactly
                `compact_remnant`, unwrapped. Used by
                `phenomenonGen.py --anchor-system`.
            primary_star_params (dict, optional): A
                `stellarEvolution.star_params` dict to build the primary
                from with `Star.from_params`, rolling nothing about it: a
                pre-placed bright star, or a star drawn for its position's
                stellar population (`star_population`). Everything else
                (companion, planets, moons, belts, comets) is generated as
                usual around it.
            member_tree (stored_trees.StoredPlan, optional): The tree a scattered row stored (GEN.315), with
                `primary_star_params`: the system takes its other stars and the orbits between them from the tree and
                draws none of its own (`BINARY_SYSTEM` and `WIDE_BINARY` are not read). Each star is built from its
                stored seed, birth mass and the primary's age; a member that is not a star (a brown dwarf, a neutron star, a black
                hole) hangs in the tree as a `dark_members` entry (GEN.1013). A plan that is not `complete` (a tree with a
                member of a kind nothing builds) builds the primary alone. `None`: the system draws its own tree of stars
                (`_drawn_plan`, GEN.631).
        """
        self.system_config = system_config # Assign the passed SystemConfig instance
        # A run's prevalences (GEN.52) decide the measured features' flags
        # for this system before anything is drawn.
        prevalence.resolve(self.system_config)
        # Rolled once here (not left for each Star/BinaryStarProxy to roll
        # its own) and passed identically to every constituent star below --
        # a binary pair's AU-scale separation is negligible next to its
        # light-year-scale galactic orbit radius, so both stars (and the
        # proxy standing in for the pair) move around the galaxy together,
        # sharing one phase, not three independent ones. See
        # `Star.__init__`'s `galactic_orbital_phase_deg` docstring.
        galactic_orbital_phase_deg = draw.uniform(0, 360)
        if compact_remnant is not None:
            self.star = compact_remnant
        elif primary_star_params is not None:
            self.star = Star.from_params(primary_star_params, self.system_config, name=self.system_config.NAME,
                                         galactic_center_dist_ly=galactic_center_dist_ly,
                                         galactic_orbital_phase_deg=galactic_orbital_phase_deg)
        else:
            with log.timed_phase("primary star generation"):
                self.star = Star(self.system_config, name=self.system_config.NAME,
                                  galactic_center_dist_ly=galactic_center_dist_ly,
                                  galactic_orbital_phase_deg=galactic_orbital_phase_deg) # Pass system_config and use its NAME
        self.primary_star = self.star # For single star systems, the primary is the star
        self.placed_star = self.star
        self.planets = []
        self.secondary_planets = [] # Only populated for an S-type (wide) binary's secondary star
        self.binary_type = None # None | "close" (P-type/circumbinary) | "wide" (S-type)
        self.wide_binary = None # A WideBinaryPair, only when binary_type == "wide"
        # The nebula or supernova remnant around the system (the dict of `generation.cloud`): set by
        # `store.load_star_system`, and by `apply_cloud` when generation placed the system inside one (GEN.357).
        self.surrounding_cloud = None
        self.cloud_effects = None  # the `cloud.CloudEffects` of that cloud, set by `apply_cloud`; derived, not stored
        self.stars = [self.primary_star] # Keep track of individual stars

        self.tree_built = member_tree is not None and member_tree.complete
        placed_star = self.primary_star
        drawn = member_tree is None and compact_remnant is None
        plan = stored_trees.SINGLE if member_tree is None else member_tree
        wanted_kind = self._wanted_pair_kind() if drawn else None
        attempts = tuning.COMPANION_DRAW_ATTEMPTS if drawn else 1
        for attempt in range(attempts):
            if drawn:
                plan = self._drawn_plan()
            self._add_companions(plan, galactic_center_dist_ly, galactic_orbital_phase_deg)
            if not drawn or self._companions_fit(wanted_kind):
                break
            if attempt + 1 < attempts:
                self._remove_companions(placed_star)
        else:
            # No draw fit what the directives and the place asked: a forced pair is kept as the last draw made it (of the
            # kind `WIDE_BINARY` asks for, whatever its axis), any other system is the star alone.
            if (self._forced_star_count() or 1) < 2:
                self._remove_companions(placed_star)
            elif self.system_config.WIDE_BINARY is not None and self.binary_type != wanted_kind:
                self._remove_companions(placed_star)
                self._add_companions(plan, galactic_center_dist_ly, galactic_orbital_phase_deg, kind=wanted_kind)

        # Removed: self.system_flavor_count = 0 # Initialize system flavor count
        # `validate_system`'s orbital-overlap correction can, in rare cases,
        # push a deliberately-placed guaranteed body (the forced habitable
        # world, or an asteroid belt's own guaranteed fallback slot) into a
        # zone/spacing its own placement no longer supports -- retry the
        # whole placement (same star, fresh positions) rather than accept a
        # system that silently fails a requirement `system_config` explicitly
        # asked for. See `_generate_planets`/`MAX_SYSTEM_GENERATION_ATTEMPTS`.
        #
        # For an S-type binary, the primary and secondary are each generated
        # independently via their own `_generate_planets` call -- the
        # secondary always gets `apply_guarantees=False` (see that method's
        # docstring), so a HABITABLE_WORLD/ASTEROID_BELT guarantee or an
        # explicit SLOTS spec unambiguously describes the primary's own
        # list only. `_validate_cross_star_clearance` then prunes either
        # list's outermost planet if the two stars' own disks would
        # otherwise gravitationally encroach on each other -- see that
        # method's docstring.
        # TEST.313: the age settled below only rises, and the survival limit falls with age (GEN.325), so the final
        # trim can take the forced belt or habitable world that the attempts above kept. The attempts then run again
        # once, placing under the settled age, before the request is given up as unmet.
        self.build_ledger()
        self._request_unreachable = False
        for _round in range(2):
            for _attempt in range(tuning.MAX_SYSTEM_GENERATION_ATTEMPTS):
                with log.timed_phase(f"planet generation attempt {_attempt + 1}"):
                    primary_ceiling_au = self._orbit_ceiling_au(self.star)
                    self.planets = self._generate_planets(self.star, self.star.habitable_zone, primary_ceiling_au)
                    self._clear_circumbinary_floor(self.planets)
                    self.validate_system(self.planets)
                    # validate_system only ever pushes a body's distance further
                    # OUTWARD to resolve a spacing collision with its inner
                    # neighbor, never closer -- for the galactic-Hill-sphere
                    # ceiling this correction is negligible next to a light-year-
                    # scale boundary, but a much tighter S-type a_crit ceiling
                    # makes it a real possibility, so re-trim afterward using each
                    # body's own final position.
                    self._trim_to_orbit_ceiling(self.planets, primary_ceiling_au)
                    self._apply_star_history(self.planets, self.star)
                    # GEN.318: a white dwarf's survivors move out; the Hill spacing is over the lighter star now
                    if remnant_orbits.widen_for_mass_loss(self.planets, self.star, self._explicit_request()):
                        self.validate_system(self.planets)
                        self._trim_to_orbit_ceiling(self.planets, primary_ceiling_au)

                    self.secondary_planets = []
                    if self.binary_type == "wide":
                        secondary_ceiling_au = self._orbit_ceiling_au(self.secondary_star)
                        self.secondary_planets = self._generate_planets(
                            self.secondary_star, self.secondary_star.habitable_zone,
                            secondary_ceiling_au, apply_guarantees=False,
                        )
                        self.validate_system(self.secondary_planets)
                        self._trim_to_orbit_ceiling(self.secondary_planets, secondary_ceiling_au)
                        self._apply_star_history(self.secondary_planets, self.secondary_star)
                        if remnant_orbits.widen_for_mass_loss(self.secondary_planets, self.secondary_star,
                                                              self._explicit_request()):
                            self.validate_system(self.secondary_planets)
                            self._trim_to_orbit_ceiling(self.secondary_planets, secondary_ceiling_au)
                        # Runs last, once both lists reflect their own truly final
                        # (post-validate_system, post-trim) positions -- see this
                        # method's own docstring.
                        self._validate_cross_star_clearance()

                    self.member_planets = [self._generate_member_planets(member) for member in self.member_stars]

                    habitable_satisfied = self.system_config.HABITABLE_WORLD is not True or habitable_request.count(self.planets) > 0
                    belt_satisfied = self.system_config.ASTEROID_BELT is not True or self.count_objects(self.planets)[1] > 0
                    class_satisfied = self._forced_class_present()
                if habitable_satisfied and belt_satisfied and class_satisfied:
                    log.debug(f"Planet generation attempt {_attempt + 1}: kept ({len(self.planets)} primary + "
                              f"{len(self.secondary_planets)} secondary bodies; HABITABLE_WORLD="
                              f"{self.system_config.HABITABLE_WORLD}, ASTEROID_BELT={self.system_config.ASTEROID_BELT})")
                    break
                draw.redrawn("system.planets_regenerated")
                log.debug(f"Planet generation attempt {_attempt + 1}/{tuning.MAX_SYSTEM_GENERATION_ATTEMPTS}: "
                          f"retrying (habitable world required and found: {habitable_satisfied}, asteroid belt "
                          f"required and found: {belt_satisfied})")
                if self._request_unreachable:
                    break

            if remnant_orbits.is_post_agb_pair(self.star):
                remnant_orbits.post_agb_disk(self.primary_star, self._orbit_floor_au(self.star))

            # A pair's two stars share one age (GEN.53): a close pair's proxy
            # settles it for both, a wide pair's from both stars' own planets.
            if self.binary_type == "wide":
                adjust_pair_age_for_planets(self.primary_star, self.secondary_star, self.planets, self.secondary_planets,
                                            self.member_hosts())
            elif self.binary_type == "close":
                self.star.adjust_age_for_planets(self.planets, self.member_hosts())
            else:
                self.star.adjust_age_for_planets(self.planets)

            # The age above only ever rises, and the place's survival limit falls with age (GEN.325), so a body placed
            # under the provisional age can sit past the ceiling of the final one: trim to it, and judge the forced
            # habitable world and belt on what is left.
            self._trim_to_orbit_ceiling(self.planets, self._orbit_ceiling_au(self.star))
            if self.binary_type == "wide":
                self._trim_to_orbit_ceiling(self.secondary_planets, self._orbit_ceiling_au(self.secondary_star))
            for member, bodies in self.member_hosts():
                self._trim_to_orbit_ceiling(bodies, self._orbit_ceiling_au(member))
            self._settle_to_final_age()
            habitable_satisfied = self.system_config.HABITABLE_WORLD is not True or habitable_request.count(self.planets) > 0
            belt_satisfied = self.system_config.ASTEROID_BELT is not True or self.count_objects(self.planets)[1] > 0
            class_satisfied = self._forced_class_present()
            if (habitable_satisfied and belt_satisfied and class_satisfied) or _round or self._request_unreachable:
                break
            draw.redrawn("system.planets_regenerated_after_age")
        self.unmet_requirements = []
        if not habitable_satisfied:
            self.unmet_requirements.append("a habitable world")
        if not belt_satisfied:
            self.unmet_requirements.append("an asteroid belt")
        if not class_satisfied:
            self.unmet_requirements.append(f"a planet of class {self.system_config.FORCE_CLASS}")
        forced = self._forced_star_count()
        if forced is not None and len(self.stars) != forced and forced > 1:
            # A star under the hydrogen-burning limit has no lighter star to be its companion (GEN.631); a brown dwarf is no
            # star of the count (GEN.899), so the refusal says so (GEN.706).
            wanted = "a companion star" if forced == 2 else f"{forced} stars"
            if place_terms.born_mass_solar(self.primary_star) < multistar.companion_floor_solar(
                    getattr(self.primary_star, "metallicity_fe_h", None)):
                wanted += " (a primary under the hydrogen-burning limit has no lighter star, and a brown dwarf is no star of the count)"
            self.unmet_requirements.append(wanted)
        if self.unmet_requirements:
            log.debug(f"Planet generation gave up after {tuning.MAX_SYSTEM_GENERATION_ATTEMPTS} attempts "
                      f"without {' or '.join(self.unmet_requirements)}")

        # Comets don't participate in the planet-slot retry loop above --
        # see _generate_comets' own docstring for why -- so they're rolled
        # once here, after that loop has settled on a final planet/belt
        # layout, rather than inside it.
        with log.timed_phase("comet generation"):
            wide = self.binary_type == "wide"
            self.comets = self._generate_comets(
                self.star, validation.orbit_ceiling_au(self.primary_star) if wide else None)
            self.secondary_comets = self._generate_comets(
                self.secondary_star, validation.orbit_ceiling_au(self.secondary_star)) if wide else []
            # A wide pair a few AU apart leaves each star a stable limit under any comet's orbit (GEN.448): a forced
            # comet then has no room, like a forced habitable world or belt.
            if self.system_config.COMETS is True and not (self.comets or self.secondary_comets):
                self.unmet_requirements = [*self.unmet_requirements, "a comet"]

        # Each hosting star's own reflex-offset "wobble" (schema v20) --
        # see `Star.reflex_offset_x`'s own docstring and
        # `orbits.calculate_reflex_offset`. Positions are now final (planet
        # placement/validation above is done), so this is computed once,
        # here, same as everything else in this pass.
        self.binary_planetary_wobble_x = self.binary_planetary_wobble_y = self.binary_planetary_wobble_z = 0.0
        # A planet's own wobble from its moons was made when the moons were; their positions have moved since (the orbits were
        # settled and the bodies made again at the star's final age), so it is made again from the final ones.
        for planet in self.all_planets:
            if planet.body_type != "a":
                validation.refresh_planet_reflex(planet)
        if self.binary_type == "close":
            # Circumbinary planets orbit the merged proxy, not either
            # individual star (star_id IS NULL -- see schema.sql), so their
            # combined pull is modeled as one shared wobble on the whole
            # pair rather than split between primary/secondary -- splitting
            # it would require solving a real 3+-body problem, out of scope
            # (see the "Investigated and explicitly out of scope" plan note
            # on a related, larger case).
            self.binary_planetary_wobble_x, self.binary_planetary_wobble_y, self.binary_planetary_wobble_z = \
                calculate_reflex_offset(self.star.mass, self._reflex_children(self.planets))
        else:
            self.primary_star.reflex_offset_x, self.primary_star.reflex_offset_y, self.primary_star.reflex_offset_z = \
                calculate_reflex_offset(self.primary_star.mass, self._reflex_children(self.planets))
            if self.binary_type == "wide":
                self.secondary_star.reflex_offset_x, self.secondary_star.reflex_offset_y, self.secondary_star.reflex_offset_z = \
                    calculate_reflex_offset(self.secondary_star.mass, self._reflex_children(self.secondary_planets))
        for member, bodies in self.member_hosts():
            member.reflex_offset_x, member.reflex_offset_y, member.reflex_offset_z = \
                calculate_reflex_offset(member.mass, self._reflex_children(bodies))

        # Decided once, here, at generation time -- __str__ (which may be called
        # more than once) reads this rather than re-rolling and double-counting
        # system_flavor_count on every render.
        self.system_flavor_text = None
        if draw.random() < tuning.FLAVOR_CHANCE_SYSTEM and self.system_config.system_flavor_count < tuning.MAX_FLAVOR_TOTAL:
            self.system_flavor_text = draw.choice(tuning.SYSTEM_FLAVOR)
            self.system_config.system_flavor_count += 1
            log.choice("System flavor text", self.system_flavor_text,
                       f"roll passed FLAVOR_CHANCE_SYSTEM ({tuning.FLAVOR_CHANCE_SYSTEM}) and "
                       f"system_flavor_count ({self.system_config.system_flavor_count - 1}) was under "
                       f"MAX_FLAVOR_TOTAL ({tuning.MAX_FLAVOR_TOTAL})")
        else:
            log.debug("System flavor text: none (roll failed FLAVOR_CHANCE_SYSTEM or MAX_FLAVOR_TOTAL reached)")

        # All planets and moons are generated without life data (see physics.planets's
        # module docstring). Apply it now, in one pass over the finished system, so
        # evolutionary timelines are computed against the star's final, planet-adjusted
        # age rather than its provisional pre-adjustment one. Flavor text is decided
        # in the same pass (after life data, since it reads life_stage) so it
        # too is a fixed, pre-rendered fact rather than something __str__ rolls.
        with log.timed_phase("snow line pass"):
            for obj in self.all_planets:
                if obj.body_type == 'a':  # A belt is not a body; GEN.442 gives it the stored line to read.
                    continue
                obj.born_beyond_snow_line = bool(obj.distance >= obj.star.snow_line_au)
                for moon in obj.moons:
                    moon.born_beyond_snow_line = obj.born_beyond_snow_line

        with log.timed_phase("ring pass"):
            for obj in self.all_planets:
                rings.apply(obj)  # GEN.732, GEN.1057: after the moons and the snow line are final
                for moon in getattr(obj, "moons", ()):  # a belt has none
                    rings.apply(moon, obj)  # GEN.836

        with log.timed_phase("life data pass"):
            insolation.refresh_system(self)  # GEN.148: the light each planet sums over the stars, from their final orbits
            for obj in self.all_planets:
                if obj.body_type == 'a': # Skip asteroid belts; they carry no life data.
                    continue
                interior_heat.update_heat(obj)  # GEN.763: the thermal history to the star's final age (no draws)
                planetPhysics.refresh_at_final_age(obj)  # GEN.839: the ice shell and the volcanism gate follow it (the body's own stream)
                planetPhysics.apply_class(obj)  # GEN.689: the light each body gets is final only now
                magnetism.generate_stellar_induction(obj)  # GEN.981: the star's field at the orbit, from the final water and core
                for moon in obj.moons:
                    planetPhysics.apply_class(moon)
                planetLife.apply_life_data(obj)
                planetLife.decide_flavor_text(obj)
                for moon in obj.moons:
                    interior_heat.update_heat(moon)
                    planetPhysics.refresh_at_final_age(moon)
                    magnetism.generate_induced_field(moon, obj)  # GEN.178: from the final water and orbit
                    radiation.apply_parent_belts(moon, obj)  # GEN.983: the parent's belts at the moon's orbit; the life pass rebuilds the dose from it
                    planetLife.apply_life_data(moon)
                    planetLife.decide_flavor_text(moon)
            for obj in self.all_planets:  # GEN.1034: the final pass made the water again (a brine seep can end as an ice-covered ocean), and the tides follow it
                if obj.body_type == 'a':
                    continue
                ocean_tides.generate(obj)
                for moon in obj.moons:
                    ocean_tides.generate(moon, obj)
            # The final pass can change a class, so a forced class is judged again on the final classes: one the final pass took away is
            # unmet, one it gave is met.
            wanted = f"a planet of class {self.system_config.FORCE_CLASS}"
            present = self._forced_class_present()
            if not present and wanted not in self.unmet_requirements:
                self.unmet_requirements.append(wanted)
            elif present and wanted in self.unmet_requirements:
                self.unmet_requirements.remove(wanted)
            # The same for a habitable world (GEN.916): life changes the air of a body after the retry judged it, and with it the
            # temperature and the water, so the request is judged again on the finished bodies.
            if self.system_config.HABITABLE_WORLD is True:
                wanted = "a habitable world"
                present = habitable_request.count(self.planets) > 0
                if not present and wanted not in self.unmet_requirements:
                    self.unmet_requirements.append(wanted)
                elif present and wanted in self.unmet_requirements:
                    self.unmet_requirements.remove(wanted)

        # The system takes the name the star was generated (or configured)
        # with; every star, planet and moon inside it is named from that --
        # see names/bodies.py.
        self._name = self.star.name
        self.star_words = None
        if self.binary_type == "wide":
            self.star_words = self._draw_star_words(wide_pair_first_word(self._name))
        self.assign_names()

        self.planet_count, self.belt_count, self.moon_count = self.count_objects()
        self.hab_count, self.m_count = self.count_habitable()
        self.habitable_class_count = self.count_habitable_class_worlds()
        self.habitable_request = self.habitable_request_status()
        self.comet_count = self.count_comets()
        log.debug(f"System {self.name!r} finished: binary={self.binary_type or 'no'}, "
                  f"{self.planet_count} planets, {self.belt_count} belts, {self.moon_count} moons, "
                  f"{self.comet_count} comets, {self.hab_count} habitable")

    @property
    def name(self):
        """The system's own name. A single star (or a close pair's proxy)
        carries it, so it's read from `star`; a wide pair's `star` is its
        primary, whose name adds its own word (see `names/bodies.py`), so
        that system keeps its name separately."""
        if self.binary_type == "wide" and getattr(self, "_name", None):
            return self._name
        return self.star.name

    @name.setter
    def name(self, value):
        self._name = value
        if self.binary_type != "wide":
            self.star.name = value

    @property
    def all_planets(self):
        """list: Every planet and belt of the system: the primary's (or a close pair's), the wide secondary's, then each
        member star's."""
        return [*self.planets, *self.secondary_planets, *(body for bodies in self.member_planets for body in bodies)]

    def apply_cloud(self, cloud, point_pc, rng=draw):
        """
        Puts the system inside `cloud` (GEN.357): the cloud becomes `surrounding_cloud` and its effects at `point_pc`, the
        system's position in the galaxy frame (parsecs), become `cloud_effects` (`cloud.effects_of`: the G0 of the cloud's
        hosts, the radiation tier and the 26Al enrichment). Called once the position is chosen, by `cloud.apply_in_sector`
        while a sector is generated and by GEN.99's regenerate in place. It leaves the planets as drawn: the rules that
        change them are GEN.358's, and nothing here is stored.

        Args:
            cloud (dict): As `generation.cloud` describes it.
            point_pc (tuple): The system's `(x, y, z)` in the galaxy frame, parsecs.
            rng: Where the enrichment is drawn from (`cloud.apply_in_sector` gives it a stream of its own).

        Returns:
            cloud.CloudEffects: The effects, also on `cloud_effects`.
        """
        self.surrounding_cloud = cloud
        self.cloud_effects = cloudContext.effects_of(cloud, point_pc, rng)
        return self.cloud_effects

    def member_hosts(self):
        """`[(star, planets)]` for each of `member_stars`."""
        return list(zip(self.member_stars, self.member_planets or [[] for _ in self.member_stars]))

    def _tree_members(self):
        """The `TreeMember`s of the system's tree: each star of `stars` (the member number is its index there), then each of
        `dark_members` (numbered on from the stars), which weigh in the tree and give no light."""
        if not self.orbit_nodes:
            return [TreeMember(0, None, 0, _kg_to_msun(self.star.mass))]
        slots = self.star_slots or ((0, 0), (0, 1))
        members = [TreeMember(number, node_no, slot, _kg_to_msun(star.mass))
                   for number, (star, (node_no, slot)) in enumerate(zip(self.stars, slots))]
        members += [TreeMember(len(self.stars) + index, dark.parent_node_no, dark.slot, dark.mass_solar)
                    for index, dark in enumerate(self.dark_members)]
        return members

    @property
    def ledger(self):
        """
        The clearance ledger of the stars as they are now (`clearance_ledger.Ledger`, GEN.1014): the window of every host (each star
        of a wide pair or a further member, and a close pair's merged proxy), from the engulfment floor and the Holman and Wiegert
        limits of the tree to the tightest S-type limit above it. It is the one place the floor and the ceiling of a star's planets
        come from; derived from the tree and the stars, never stored. Built when first read after the companions change.
        """
        if self._ledger is None:
            explicit = self._explicit_request()
            self._ledger = clearance_ledger.Ledger(
                self.orbit_nodes, self._tree_members(),
                {number: self._engulfment_radius_au(star, explicit) for number, star in enumerate(self.stars)})
        return self._ledger

    def _ledger_host(self, star):
        """The `(kind, number)` of the ledger host that is `star`: a close pair's proxy is its pair node, any other star its number."""
        if self.binary_type == "close" and star is self.star:
            return planet_zones.P_TYPE, self.pair_node_no()
        return planet_zones.S_TYPE, next((number for number, other in enumerate(self.stars) if other is star), 0)

    def _planet_zones(self):
        """The windows of the ledger: the S-type window of every star in `stars` order and the P-type window of every node, each a
        `planet_zones.Zone`. A dark member has a window too (`Zone.host_no` from `len(stars)` on), which holds no planets."""
        return self.ledger.zones()

    def sibling_stars(self, number):
        """`(node, star numbers)`: the node star `number` of `stars` hangs from and the stars of the subtree beside it (a
        dark member beside it gives no light and is not listed)."""
        node, beside = planet_zones.sibling_of(self.orbit_nodes, self._tree_members(), number)
        return node, sorted(index for index in beside if index < len(self.stars))

    def build_ledger(self):
        """
        Builds the clearance ledger again and gives each star that has a companion above it, and a close pair's proxy, the outer limit
        of its window (GEN.332, GEN.1014): a member's S-type window against what hangs beside it, and a close pair's P-type window up
        to the S-type limit of the pair against the stars beside it, each held by every level above. Both are `a_crit_au`, which
        `validation.orbit_ceiling_au` reads for a wide pair's stars as well; the ledger is the only place that works it out. A lone
        star has none.
        """
        self._ledger = None
        if not self.orbit_nodes or (len(self.stars) < 2 and not self.dark_members):
            return
        ledger = self.ledger
        hosts = [(planet_zones.S_TYPE, number, star) for number, star in enumerate(self.stars)
                 if self.binary_type != "close" or number >= 2]
        if self.binary_type == "close":
            hosts.append((planet_zones.P_TYPE, self.pair_node_no(), self.star))
        for kind, number, star in hosts:
            limit = ledger.stability_limit_au((kind, number))
            if limit is not None:
                star.a_crit_au = limit

    def pair_node_no(self):
        """The node of the pair (`primary_star` and `secondary_star`): the root of a plain pair, else the lowest node above both
        (the node they hang from, unless a dark member or a sub-pair is between them), `None` for a single star."""
        if self.binary_type is None:
            return None
        if self.star_slots is None:
            return self.orbit_nodes[0].node_no
        parent = {node.node_no: node.parent_node_no for node in self.orbit_nodes}
        above = set()
        node_no = self.star_slots[0][0]
        while node_no is not None:
            above.add(node_no)
            node_no = parent[node_no]
        node_no = self.star_slots[1][0]
        while node_no not in above:
            node_no = parent[node_no]
        return node_no

    @staticmethod
    def _draw_star_words(first_word):
        """A wide pair's two star words: the primary's from the star name
        lists, the secondary's from the small/little/child ones (GEN.62),
        each distinct from the other and from the shared first word."""
        primary_word = generate_star_word(exclude=(first_word,))
        return primary_word, generate_star_word(exclude=(first_word, primary_word), companion=True)

    def assign_names(self, name=None):
        """
        Names every star, planet and moon in this system from the system
        name (see `names/bodies.py`'s module docstring for the scheme).

        Args:
            name (str, optional): A new system name to adopt first --
                `store.insert_star_system` passes the name uniqueness
                reservation settled on, which may carry a decoration the
                generated name didn't.
        """
        if name is not None:
            self.name = name
        system_name = self.name
        if self.binary_type is None:
            self.star.name = system_name
            hosts = [(system_name, self.planets)]
        else:
            if self.binary_type == "close":
                # One system, so its stars are its A and B (GEN.62).
                primary_letter, secondary_letter = CLOSE_PAIR_LETTERS
                self.primary_star.name = f"{system_name} {primary_letter}"
                self.secondary_star.name = f"{system_name} {secondary_letter}"
                self.star.name = system_name
                hosts = [(system_name, self.planets)]
            else:
                first_word = wide_pair_first_word(system_name)
                if getattr(self, "star_words", None) is None:
                    # A system reloaded from an export written before star words existed.
                    self.star_words = self._draw_star_words(first_word)
                primary_word, secondary_word = self.star_words
                # Two words each, the first shared; the primary's planets
                # take that first word, the secondary's its own (GEN.62).
                self.primary_star.name = f"{first_word} {primary_word}"
                self.secondary_star.name = f"{first_word} {secondary_word}"
                hosts = [(first_word, self.planets), (secondary_word, self.secondary_planets)]
        for index, member in enumerate(self.member_stars):    # C, D, ... after the pair's A and B (GEN.631)
            member.name = f"{system_name} {chr(ord(CLOSE_PAIR_LETTERS[1]) + 1 + index)}"
            hosts.append((member.name, self.member_planets[index] if self.member_planets else []))    # GEN.332
        for prefix, bodies in hosts:
            name_bodies(prefix, bodies)

    def _generate_comets(self, star, ceiling_au=None):
        """
        Rolls this star's own native comet population (see
        `cometData.Comet` -- contrast the always-standalone, unbound
        `roguePlanetData.InterstellarComet`, which `StarSystem` never
        generates).

        Comets deliberately don't participate in the planet-slot
        generation/retry loop above (`_generate_planets`/
        `validate_system`/`_orbit_ceiling_au`): a `Comet`'s `distance_au`
        is its current, continuously-varying position along a highly
        eccentric/arbitrarily inclined orbit, not a fixed slot in the
        orderly, Hill-sphere-spaced sequence a `Planet`'s/`AsteroidBelt`'s
        own `distance` occupies, so there's no orbital-overlap invariant
        for a comet to satisfy there. Instead this rolls independently,
        once per star, honoring `SystemConfig.COMETS`'s tri-state contract
        the same way `ASTEROID_BELT`/`HABITABLE_WORLD` do elsewhere:
        `True` forces at least `tuning.SYSTEM_COMET_COUNT_RANGE[0]`
        comets, `False` forces none, `None` rolls
        `tuning.SYSTEM_COMET_CHANCE`.

        Args:
            star: The `Star` or `BinaryStarProxy` this comet population is
                 bound to -- its `.mass` (the pair's own combined mass for
                 a close/P-type binary proxy, via `BinaryStarProxy.mass`)
                 drives each comet's own orbital period, the same
                 `star.mass` a `Planet`'s own period calculation reads.
            ceiling_au (float, optional): For a star of a wide binary, its
                 stable limit (`validation.orbit_ceiling_au`): every comet
                 is an elliptical one whose whole orbit stays inside it
                 (GEN.228), re-rolled up to
                 `tuning.COMET_CEILING_ATTEMPTS` times and dropped when
                 none fits, so no comet's path crosses the companion. Without
                 it (a single star or a close pair) an elliptical comet's
                 orbit stays inside the star's own Hill sphere
                 (`validation.orbit_ceiling_au`, GEN.591) and a parabolic
                 comet leaves it by design.

        Returns:
            list: Newly generated `Comet` objects, zero or more.
        """
        if self.system_config.COMETS is False:
            return []

        has_comets = self.system_config.COMETS is True or draw.random() < prevalence.scaled_chance(
            self.system_config, "comets", tuning.SYSTEM_COMET_CHANCE)
        if not has_comets:
            return []

        count = draw.randint(*tuning.SYSTEM_COMET_COUNT_RANGE)
        primary_mass_solar = _kg_to_msun(star.mass)
        if ceiling_au is None:
            # A closed orbit stays inside the star's own Hill sphere (GEN.591); a parabolic comet leaves it by design.
            hill_au = validation.orbit_ceiling_au(star)
            return [Comet(self.system_config, primary_mass_solar, ceiling_au=hill_au, star=star) for _ in range(count)]
        comets = []
        for _ in range(count):
            for _attempt in range(tuning.COMET_CEILING_ATTEMPTS):
                comet = Comet(self.system_config, primary_mass_solar, orbit_type="elliptical", ceiling_au=ceiling_au, star=star)
                if validation.comet_fits_ceiling(comet, ceiling_au):
                    comets.append(comet)
                    break
                draw.redrawn("comet.ceiling")
        return comets

    @staticmethod
    def _reflex_children(planets):
        """
        Builds the `(mass_kg, x_au, y_au, z_au)` list `utils.
        calculate_reflex_offset` expects, out of a `planets`-shaped list
        (`self.planets`/`self.secondary_planets`) -- excludes
        `AsteroidBelt` entries (`body_type == 'a'`), which carry no mass
        anywhere in this codebase (see `AsteroidBelt`'s own class
        docstring), the same exclusion `store.advance_orbital_phases`'s own
        `planets`-table `UPDATE`s get for free (asteroid belts are stored
        in their own separate table).

        Args:
            planets (list): `Planet`/`AsteroidBelt` objects.

        Returns:
            list: `(mass_kg, x_au, y_au, z_au)` tuples, one per real planet.
        """
        return [(p.mass, p.position_x, p.position_y, p.position_z) for p in planets if p.body_type != 'a']

    def _special_companion_chance(self):
        """
        The chance of a companion the primary's own kind sets in place of the star-count table's (`None`: the table's):
        a protostar has the multiplicity of Perseus' Class 0 and Class I protostars (GEN.242), a Wolf-Rayet star a companion
        in 52% (WN) or 90% (WC and WO) of cases (GEN.243).
        """
        primary = self.primary_star
        initial_mass = getattr(primary, "initial_mass_sol", None)
        protostar_chance = (None if initial_mass is None or primary.age is None else
                            pre_main_sequence.companion_chance(initial_mass, primary.age * 1.0e9))
        if protostar_chance is not None:
            return protostar_chance
        if getattr(primary, "yerkes_class", None) == wolf_rayet.CLASS:
            return wolf_rayet.companion_chance(primary.type)
        return None

    def _wanted_pair_kind(self):
        """
        The kind of pair (`close` or `wide`) this system's directives ask for, `None` for whatever the tree draws:
        `WIDE_BINARY` forces it, and a `wide_binary` prevalence (ADM.45) moves the share of wide pairs among the pairs.
        """
        config = self.system_config
        if config.WIDE_BINARY is not None:
            return "wide" if config.WIDE_BINARY else "close"
        if prevalence.percent(config, "wide_binary"):
            chance = prevalence.scaled_chance(config, "wide_binary", prevalence.USUAL_SHARES["wide_binary"])
            return "wide" if draw.random() < chance else "close"
        return None

    def _drawn_plan(self):
        """
        The tree of stars around the primary this system draws for itself (`multistar.draw_member_tree`, at the odds of the
        `star-count` table), as a `stored_trees.StoredPlan`, every member in it: the ones a build cannot make from a birth mass
        and an age (a brown dwarf, a neutron star, a black hole) stay as dark members (GEN.1013). `BINARY_SYSTEM` forces no companion
        (`False`) or a pair (`True`); a protostar's or Wolf-Rayet star's own chance of a companion (`_special_companion_chance`)
        and a `binary_system` prevalence move how many systems have one.
        """
        config = self.system_config
        primary = self.primary_star
        if config.BINARY_SYSTEM is False or primary.age is None:
            return stored_trees.SINGLE
        count = self._forced_star_count()
        at_least_two = False
        special = self._special_companion_chance()
        percent = prevalence.percent(config, "binary_system")
        placed = multistar.Placed("star", _kg_to_msun(primary.mass), place_terms.born_mass_solar(primary),
                                  fe_h=getattr(primary, "metallicity_fe_h", None), record=member_record.of_star(primary))
        if count is not None:
            pass
        elif special is not None:
            if draw.random() < special:
                at_least_two = True
            else:
                count = 1
        elif percent < 0.0:
            if draw.random() < -percent / 100.0:
                count = 1
        elif percent > 0.0:
            weights = multistar.star_count_weights(placed.birth_mass_solar)
            share = 1.0 - weights[0] / sum(weights)     # this star's usual chance of a companion
            at_least_two = draw.random() < min(1.0, (percent / 100.0) * share / max(1.0 - share, 1.0e-9))
        age_yr = primary.age * 1.0e9

        # The widest axis the place keeps (GEN.103), for a pair of the primary and a companion of the mean mass ratio.
        terms = getattr(config, "PLACE_TERMS", None) or place_terms.SOLAR_NEIGHBOURHOOD
        limit_au = place_terms.wide_axis_limit_au(
            terms, placed.birth_mass_solar * (1.0 + sum(tuning.BINARY_MASS_RATIO_RANGE) / 2.0), primary.age)

        # A star of a type the directives named can be older than the model's lifetime for its mass: the tree is drawn at
        # the end of that lifetime, the star itself keeps its age.
        draw_age_yr = min(age_yr, 0.95 * evolution.main_sequence_lifetime_gy(placed.birth_mass_solar) * 1.0e9)

        def tree(n_stars):
            for tree_age_yr in (age_yr, draw_age_yr):
                try:
                    return multistar.draw_member_tree(
                        placed, tree_age_yr, max_stars=config.MAX_STARS, n_stars=n_stars,
                        prevalence=prevalence.multi_star_scale(config), max_axis_au=limit_au, pair_placed=True)
                except ValueError:    # the placed star has collapsed by this age
                    draw.redrawn("tree.drawn_at_earlier_age")
                    continue
            return multistar.draw_member_tree(placed, 0.0, n_stars=1)

        drawn = stored_trees.drawn(tree(count))
        if at_least_two and len(drawn.lit_members) < 2:
            draw.redrawn("tree.second_star_forced")
            drawn = stored_trees.drawn(tree(2))
        return drawn

    def _forced_star_count(self):
        """The number of stars the directives force (`STARS`, else two for `BINARY_SYSTEM`), `None` to leave it to the table."""
        config = self.system_config
        if config.STARS is not None:
            return config.STARS
        return 2 if config.BINARY_SYSTEM is True else None

    def _companions_fit(self, wanted_kind):
        """Whether the companions just added are what the directives and the place need (see `__init__`'s loop)."""
        config = self.system_config
        forced = self._forced_star_count()
        if forced is not None and len(self.stars) != forced:
            draw.redrawn("companions.star_count")
            return False
        if self.binary_type is None and not self.dark_members:
            return True
        if wanted_kind is not None and self.binary_type != wanted_kind:
            draw.redrawn("companions.pair_kind")
            return False
        if not self._built_stars_clear():
            draw.redrawn("companions.overlap")
            return False
        # A close pair whose whole habitable zone lies inside its circumbinary stability limit can't host a habitable world.
        if (self.binary_type == "close" and config.HABITABLE_WORLD is True and config.WIDE_BINARY is None
                and self._orbit_floor_au(self.star) >= self.star.habitable_zone[1]):
            draw.redrawn("companions.circumbinary_floor")
            return False
        return True

    def _built_stars_clear(self):
        """
        Whether the stars as built clear the orbits of the tree they hang in. A drawn tree is built again from each member's
        birth mass and the system's age, so a member the tree evolved on (a white dwarf, a stripped core) can come out a giant
        wider than its orbit: the pair, or a third star about it, then lies inside a star (`multistar.built_tree_clears`).
        """
        if not self.orbit_nodes:
            return True
        slots = self.star_slots or tuple((self.orbit_nodes[0].node_no, slot) for slot in (0, 1))
        stars = [(parent, _km_to_au(star.radius), _kg_to_msun(star.mass))
                 for star, (parent, _slot) in zip(self.stars, slots)]
        stars += [(dark.parent_node_no, dark.radius_au or 0.0, dark.mass_solar) for dark in self.dark_members]
        return multistar.built_tree_clears(self.orbit_nodes, stars)

    def _remove_companions(self, placed_star):
        """Takes the companions back out: the system is `placed_star` alone again."""
        self.primary_star = self.star = placed_star
        self.stars = [placed_star]
        self.__dict__.pop("secondary_star", None)
        self.binary_type = None
        self.wide_binary = None
        self._ledger = None
        self.member_stars = ()
        self.dark_members = ()
        self.star_slots = None
        self.orbit_nodes = ()

    def _add_companions(self, plan, galactic_center_dist_ly, galactic_orbital_phase_deg, kind=None):
        """
        Builds the other stars of `plan` around the primary (member 0 of the tree, already built) and hangs them all in the
        plan's nodes (GEN.315, GEN.631): each from its own seed, birth mass and the primary's age and [Fe/H]. The pair is
        the two built stars of the plan's pair node (`StoredPlan.pair_members`), the heavier the primary; every other built
        star is one of `member_stars`. A member that is not a star (a brown dwarf, a neutron star, a black hole) hangs in
        the nodes as one of `dark_members` (GEN.1013), so the orbits and the mass of every node are the tree's. Sets
        `orbit_nodes` and, for a tree of more than two members or fewer than two stars, `star_slots`.

        Each star is made from its member record (`member_record.make`, GEN.1009) without evolving it again: the tree drew the
        orbits from the same record, so the star is the member the orbits were drawn for. A plan whose built star is not the
        member (its mass or radius differs: the tree records a change a birth mass and an age do not give) leaves the primary
        alone and the plan unspent (`tree_built` false). A plan with no members builds nothing. `kind` forces the kind of the
        pair (`close` or `wide`) whatever its axis.
        """
        if not plan.members:
            return
        config = self.system_config
        placed = self.primary_star
        built = {0: placed}
        for member in plan.lit_members[1:]:
            record = member_record.make(member.seed, member.initial_mass_solar,
                                        member.age_yr or placed.age * 1.0e9, placed.metallicity_fe_h)
            with log.timed_phase("companion star generation"), draw.bound(member.seed):
                built[member.member_no] = Star.from_record(
                    record, config, name=placed.name, galactic_center_dist_ly=galactic_center_dist_ly,
                    galactic_orbital_phase_deg=galactic_orbital_phase_deg)
            if not stored_trees.matches(built[member.member_no], member):
                log.choice("Stored companion", "left out", f"a star of its birth mass {member.initial_mass_solar:.4g}"
                           f" Msun at {placed.age:.4g} Gy is not the stored {member.mass_solar:.4g} Msun")
                self.tree_built = False
                return
        nodes = list(plan.nodes)
        slots = {member.member_no: (member.parent_node_no, member.slot) for member in plan.members}
        self.dark_members = plan.dark_members
        pair = plan.pair_members()
        if pair is None:
            # Fewer than two built stars: the star hangs in the tree with the dark members about it, and is a single star.
            self.orbit_nodes = tuple(nodes)
            self.star_slots = (slots[0],)
            log.choice("Binary type", "none", f"a single star with {len(plan.dark_members)} dark members in its tree")
            return
        pair_no = plan.pair_node_no()
        low, high = pair
        first, second = low, high
        if built[second].mass > built[first].mass:
            # The primary is by definition the heavier star, so the pair's roles swap: the heavier takes slot 0, and the
            # orbit of the slot-1 star about it is the old one turned half way round.
            log.choice("Binary primary", "swapped",
                       f"secondary {_kg_to_msun(built[second].mass):.3g} Msun "
                       f"outweighs primary {_kg_to_msun(built[first].mass):.3g} Msun")
            first, second = second, first
            if slots[low][0] == slots[high][0] == pair_no:
                slots[low], slots[high] = slots[high], slots[low]
                nodes[pair_no] = dataclasses.replace(
                    nodes[pair_no], arg_periapsis_deg=(nodes[pair_no].arg_periapsis_deg + 180.0) % 360.0)
        # The pair is a P-type pair (`close`) where its axis is under the least a wide pair has, else an S-type one (`wide`);
        # its node says which, the other nodes of the tree keep the kind the tree drew them with.
        # A P-type pair is tracked on a circle (`BinaryStarProxy`).
        pair_node = nodes[pair_no]
        if kind is None:
            kind = "close" if _km_to_au(pair_node.separation_km) < tuning.WIDE_BINARY_SEPARATION_MIN_AU else "wide"
        if kind == "close":
            # (a circle has no periapsis: the angle from the ascending node is the argument of periapsis plus the anomaly,
            # which is the phase the pair is tracked by, `BinaryStarProxy`)
            pair_node = dataclasses.replace(
                pair_node, kind="close", eccentricity=0.0, arg_periapsis_deg=0.0,
                mean_anomaly_deg=(pair_node.mean_anomaly_deg + pair_node.arg_periapsis_deg) % 360.0)
        else:
            # (the wide pair holds eccentricities under `WIDE_BINARY_ECCENTRICITY_MAX`, where its stability limit is good)
            pair_node = dataclasses.replace(pair_node, kind="wide", eccentricity=min(
                pair_node.eccentricity, tuning.WIDE_BINARY_ECCENTRICITY_MAX * (1.0 - 1.0e-9)))
        nodes[pair_no] = pair_node
        others = [member for member in plan.lit_members if member.member_no not in (first, second)]
        direct = slots[first][0] == slots[second][0] == pair_no
        primary, secondary = built[first], built[second]
        self.primary_star, self.secondary_star = primary, secondary
        self.star = primary
        self.member_stars = tuple(built[member.member_no] for member in others)
        self.stars = [primary, secondary, *self.member_stars]
        self.orbit_nodes = tuple(nodes)
        self._ledger = None
        if others or self.dark_members or not direct:
            self.star_slots = tuple(slots[member_no] for member_no in (first, second, *(m.member_no for m in others)))
        # Which of the two real binary configurations the pair is was drawn with its node.
        self.binary_type = pair_node.kind
        log.choice("Binary type", self.binary_type, f"the {pair_node.kind} node of a tree of {len(plan.members)} stars")
        # (the pair classes take the orbit of the second about the first: a pair that is not the two children of its node, with
        # the first under the node's slot 1, has it turned half way round)
        if not direct and plan.slot_under(next(m for m in plan.members if m.member_no == first), pair_no) == 1:
            pair_node = dataclasses.replace(pair_node, arg_periapsis_deg=(pair_node.arg_periapsis_deg + 180.0) % 360.0)
        if self.binary_type == "close":
            self.star = BinaryStarProxy(config, primary, secondary, galactic_center_dist_ly=galactic_center_dist_ly,
                                        galactic_orbital_phase_deg=galactic_orbital_phase_deg, orbit=pair_node)
        else:
            # An S-type pair never merges into one effective star -- self.star stays the primary Star itself, and
            # WideBinaryPair.__init__ sets primary_star.a_crit_au / secondary_star.a_crit_au directly on the two Star
            # instances (see Star.a_crit_au's own docstring).
            self.wide_binary = WideBinaryPair(config, primary, secondary, orbit=pair_node)

    def _disk_outer_edge_au(self, star):
        """
        How far out `_generate_planets` places random slots around `star`:
        its protoplanetary disk's outer edge, the same
        `DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER` times the snow line
        `_plan_architecture` builds its oligarch chain out to (GEN.37).
        Without this edge the slots would run planets far past the disk (a red dwarf's
        out to hundreds of times its habitable zone). An evolved star or
        remnant whose present light puts that edge inside its orbit floor
        (`_orbit_floor_au`) no longer describes its disk, so it gets no
        edge.
        """
        edge_au = tuning.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER * star.snow_line_au
        return edge_au if edge_au > self._orbit_floor_au(star) else math.inf

    def _orbit_floor_au(self, star):
        """
        The innermost stable orbit around `star`, in AU, from the ledger's window of it: for a P-type (close) binary's merged proxy,
        Holman & Wiegert's circumbinary critical semi-major axis (about 2-2.4 times the pair's separation), since anything closer is
        inside, or torn apart by, the two stars' own orbit; at least `_engulfment_radius_au` for any star, since a giant (or a white
        dwarf's progenitor) has removed everything inside that.
        """
        return self.ledger.zone(self._ledger_host(star)).inner_au

    def _clear_circumbinary_floor(self, planets):
        """
        Pushes the innermost body of `planets` out past `_orbit_floor_au`
        when an explicit slot or a forced habitable-zone placement put it
        inside a close binary's stable limit; `validate_system`, which runs
        next, then re-spaces everything beyond it. A planet moved this way
        is reconciled to its new zone like any `validate_system` move.
        """
        if not planets:
            return
        floor_au = self._orbit_floor_au(self.star)
        body = planets[0]
        if body.body_type == 'a':
            shift = floor_au - body.lower_limit
            if shift > 0:
                body.distance += shift
                body.lower_limit += shift
                body.upper_limit += shift
        elif body.distance < floor_au:
            body.distance = floor_au
            self._reconcile_moved_planet(body)
        else:
            return
        draw.redrawn("planet.circumbinary_floor_push")
        log.debug(f"Innermost body moved out to the circumbinary stability limit ({floor_au:.4g} AU)")

    def _generate_planets(self, star, habitable_zone, orbit_ceiling_au, apply_guarantees=True):
        """
        Builds a fresh, sequentially-placed list of planets/asteroid belts
        around `star` -- extracted out of `__init__` so it can be retried
        wholesale (same star, fresh positions) when the result doesn't
        satisfy an explicitly requested `HABITABLE_WORLD`/`ASTEROID_BELT`
        guarantee (see `__init__`'s own retry loop), and so it can be
        called independently for each star of an S-type (wide) binary.

        Args:
            star (Star): The star to generate this list around -- `self.star`
                for a single star or a P-type (close) binary's merged
                proxy, or `self.primary_star`/`self.secondary_star`
                individually for an S-type (wide) binary.
            habitable_zone (tuple): `(inner_au, outer_au)` -- `star`'s own
                habitable zone, passed explicitly rather than re-read from
                `star.habitable_zone` inline so every call site is
                self-documenting about which zone this list is placed
                against.
            orbit_ceiling_au (float): The AU distance beyond which no
                further slot is placed -- see `_orbit_ceiling_au`, which
                every caller uses to compute this (a star's own
                `system_perimeter`, or for an S-type binary star, the
                tighter of that and its companion-driven
                `a_crit_au`). This method itself doesn't need to know
                binaries exist at all.
            apply_guarantees (bool): Whether `system_config.HABITABLE_WORLD`/
                `ASTEROID_BELT`'s forced placement and `system_config.SLOTS`'s
                explicit per-orbit spec apply to this call. True (the
                default) for a single star, a P-type proxy, and an S-type
                binary's primary. False only for an S-type binary's
                SECONDARY: its list is always fully, independently random,
                never a `SLOTS` override or a forced habitable-world/
                asteroid-belt placement -- otherwise a hand-authored
                `SLOTS` scenario or `HABITABLE_WORLD=True` guarantee would
                be ambiguous between "the primary's Nth slot" and "the
                secondary's Nth slot".

        Returns:
            list: Newly generated `Planet`/`AsteroidBelt` objects, in
                 orbital order. Does not touch `self.planets` or
                 `self.secondary_planets` -- the caller assigns the result
                 to whichever one applies.
        """
        if pulsar_planets.applies(star):
            return pulsar_planets.generate(self.system_config, star, habitable_zone, orbit_ceiling_au)
        if remnant_orbits.is_post_agb_pair(star):
            return []  # a planetary nebula's close pair holds a circumbinary disk and no planets (GEN.318)

        planets = []
        plan, system_objects = self._plan_architecture(star, orbit_ceiling_au)

        required_objects = 0
        if apply_guarantees and self.system_config.HABITABLE_WORLD is True:
            required_objects += 1
        if apply_guarantees and self.system_config.ASTEROID_BELT is True:
            required_objects += 1

        if system_objects < required_objects:
            system_objects = required_objects

        forced_class = self.system_config.FORCE_CLASS if apply_guarantees else None
        if forced_class and self.system_config.SLOTS:
            raise ValueError("FORCE_CLASS cannot be combined with SLOTS: a slot already names its class.")
        if forced_class and forced_class not in tuning.PLANET_CLASSES:
            raise ValueError(f"FORCE_CLASS {forced_class!r} is not a planet class.")
        if forced_class:
            system_objects = max(system_objects, 1)
        slots = (self.system_config.SLOTS or []) if apply_guarantees else []
        if forced_class and not slots:
            slots = [{"type": "planet", "planet_class": forced_class}]
        disk_edge_au = self._disk_outer_edge_au(star)

        if system_objects > 0:
            belt_index = draw.randint(0, system_objects - 1) if self.system_config.ASTEROID_BELT is True else -1
            if belt_index != -1:
                log.choice("Guaranteed asteroid belt slot", belt_index,
                           f"random slot in [0, {system_objects - 1}] (ASTEROID_BELT is forced True)")
            found_hab = False
            found_belt = False
            i = -1

            while i < system_objects - 1:
                i += 1
                last_asteroid = False
                slot_spec = slots[i] if i < len(slots) else None
                prev_slot_explicit = i > 0 and (i - 1) < len(slots) and slots[i - 1] is not None

                seed = plan.seed(i)
                if i > 0:
                    last_planet = planets[i - 1]
                    # The planet architecture (GEN.728, GEN.729, GEN.730) gives each slot its distance and mass: oligarchs
                    # spaced in mutual Hill radii and thinned by instability, so the gaps are what the disc built
                    # (`validate_system` keeps the stable spacing of the whole). A belt takes the place of a light planet
                    # in its slot (`BELT_REPLACES_BELOW_EARTH`).
                    if last_planet.body_type == 'a':
                        estimated_distance = max(seed.distance_au,
                                                 last_planet.upper_limit * tuning.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MIN)
                        last_asteroid = True
                    else:
                        estimated_distance = max(last_planet.distance, seed.distance_au)

                    # Each successive slot's minimum spacing scales with the
                    # previous object's own Hill radius, which itself scales
                    # with its distance -- for a system with many objects
                    # (especially several giant planets around a massive
                    # star), that compounds geometrically. `orbit_ceiling_au`
                    # (the star's own Hill sphere relative to the galaxy, or
                    # for an S-type binary star, the tighter of that and its
                    # companion-driven stability limit -- see
                    # `_orbit_ceiling_au`) is the real physical boundary past
                    # which nothing is gravitationally bound to this star at
                    # all -- stop adding slots once the sequential spacing
                    # would place one beyond it, rather than letting that
                    # compounding run unbounded. A system that runs out of
                    # stable room this way simply ends up with fewer objects
                    # than `system_objects` estimated, the same
                    # physically-honest outcome a real protoplanetary disk of
                    # finite extent would produce.
                    if estimated_distance > orbit_ceiling_au:
                        break
                    # Nothing forms past the disk's outer edge (GEN.37), but a
                    # slot a guarantee or SLOTS still needs is placed anyway.
                    pending = (slot_spec is not None or (apply_guarantees and (
                        (self.system_config.HABITABLE_WORLD is True and not found_hab)
                        or (self.system_config.ASTEROID_BELT is True and not found_belt))))
                    if estimated_distance > disk_edge_au and not pending:
                        break
                else:
                    # A close binary's first slot starts no closer than its
                    # circumbinary stability limit (0 for any other star).
                    estimated_distance = max(seed.distance_au, self._orbit_floor_au(star))

                hz = habitable_zone[0] < estimated_distance < habitable_zone[1]

                # An explicit per-slot specification takes priority over all of the
                # normal random/forced generation logic below.
                if slot_spec is not None:
                    obj = self.generate_slot_object(star, slot_spec, estimated_distance, planets=planets)
                    if getattr(obj, 'body_type', None) != 'a' and habitable_request.meets_with_moons(obj):
                        found_hab = True
                    if getattr(obj, 'body_type', None) == 'a':
                        found_belt = True
                    planets.append(obj)
                    if _exceeds_orbit_ceiling(planets[-1], orbit_ceiling_au):
                        planets.pop()
                        break
                    continue

                if apply_guarantees and self.system_config.HABITABLE_WORLD is True and not found_hab:
                    if not hz and i == 0:
                        if (estimated_distance > habitable_zone[1] or
                                0 < habitable_zone[0] - estimated_distance < 0.2 or system_objects == 1):
                            estimated_distance = self._forced_habitable_distance(planets)
                            hz = True
                    elif not hz and i > 0:
                        last_planet = planets[i - 1]
                        beyond_hz = last_planet.body_type == 'a' and last_planet.upper_limit > habitable_zone[1] or \
                                    last_planet.body_type != 'a' and (last_planet.distance + last_planet.min_orbit_distance >
                                                                 habitable_zone[1])

                        # Never retroactively overwrite a belt if ASTEROID_BELT is
                        # required: this could be the belt satisfying that guarantee
                        # (found_belt is set the moment it's placed, not re-checked
                        # here), and this mechanism has no way to know whether another
                        # one exists elsewhere to fall back on.
                        belt_is_protected = last_planet.body_type == 'a' and self.system_config.ASTEROID_BELT is True

                        if beyond_hz and not prev_slot_explicit and not belt_is_protected:
                            # planets[i - 1] itself is about to be discarded and
                            # replaced below -- exclude it from belt-avoidance so
                            # a belt being overwritten here doesn't spuriously
                            # constrain where its own replacement can land.
                            estimated_distance = self._forced_habitable_distance(planets[:i - 1])
                            planet = self.build_planet(star, habitable_zone, estimated_distance, life=True, planned_objects=system_objects)
                            planets[i - 1] = planet
                            i -= 1
                            found_hab = True
                            if _exceeds_orbit_ceiling(planets[-1], orbit_ceiling_au):
                                planets.pop()
                                break
                            continue
                        elif i == system_objects - 1:
                            estimated_distance = self._forced_habitable_distance(planets)
                            hz = True

                    if hz:
                        planet = self.build_planet(star, habitable_zone, estimated_distance, life=True, planned_objects=system_objects)
                        found_hab = True
                        planets.append(planet)
                        if _exceeds_orbit_ceiling(planets[-1], orbit_ceiling_au):
                            planets.pop()
                            break
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
                hab_still_pending = apply_guarantees and self.system_config.HABITABLE_WORLD is True and not found_hab
                belt_fallback_index = (system_objects - 2) if hab_still_pending else (system_objects - 1)
                force_belt = apply_guarantees and self.system_config.ASTEROID_BELT is True and not found_belt and i >= belt_fallback_index

                belt_fits_slot = not seed.runaway and seed.mass_kg < _mearth_to_kg(tuning.BELT_REPLACES_BELOW_EARTH)
                if self.system_config.ASTEROID_BELT is not False and (
                    force_belt or i == belt_index
                    or (belt_fits_slot and not last_asteroid and not hz and draw.random() < tuning.ASTEROID_BELT_PROBABILITY)
                ):
                    if force_belt:
                        belt_reason = "forced (last-resort fallback: ASTEROID_BELT required but not yet placed)"
                    elif i == belt_index:
                        belt_reason = "matched the pre-selected guaranteed belt slot"
                    else:
                        belt_reason = f"random roll passed ASTEROID_BELT_PROBABILITY ({tuning.ASTEROID_BELT_PROBABILITY})"
                    log.choice("Orbital slot object", "asteroid belt", f"slot {i}: {belt_reason}")
                    min_distance = estimated_distance
                    max_distance = estimated_distance * draw.uniform(tuning.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MIN, tuning.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MAX)
                    planets.append(AsteroidBelt(self.system_config, star, estimated_distance, min_distance, max_distance)) # Pass system_config
                    found_belt = True
                else:
                    planet = Planet(self.system_config, star, habitable_zone, estimated_distance, # Pass system_config
                                    planned_objects=system_objects, mass=seed.mass_kg, formation_distance_au=seed.formation_au)
                    if habitable_request.meets_with_moons(planet):
                        found_hab = True
                    planets.append(planet)

                if _exceeds_orbit_ceiling(planets[-1], orbit_ceiling_au):
                    planets.pop()
                    break

        return planets

    def _generate_member_planets(self, star):
        """
        The planets and belts of a member star (GEN.332): the same draw as a wide pair's secondary's, inside the S-type
        window of the tree (`build_ledger`), never under a `SLOTS` list or a forced habitable world.
        """
        ceiling_au = self._orbit_ceiling_au(star)
        bodies = self._generate_planets(star, star.habitable_zone, ceiling_au, apply_guarantees=False)
        self.validate_system(bodies)
        self._trim_to_orbit_ceiling(bodies, ceiling_au)
        self._apply_star_history(bodies, star)
        if remnant_orbits.widen_for_mass_loss(bodies, star, self._explicit_request()):
            self.validate_system(bodies)
            self._trim_to_orbit_ceiling(bodies, ceiling_au)
        return bodies

    def _orbit_ceiling_au(self, star):
        """
        The AU distance beyond which `_generate_planets` stops placing
        slots around `star` -- that star's own Hill sphere relative to the
        galaxy (`system_perimeter`), or, for a star that is one constituent
        of an S-type (wide) binary, the tighter of that and its
        companion-driven Holman & Wiegert stability limit (`a_crit_au`,
        see `Star.a_crit_au`'s own docstring). `a_crit_au` is `None` for a
        single star or either star of a P-type (close) binary (whose
        merged `BinaryStarProxy` never sets it), so this is a no-op
        fallback to `system_perimeter` alone in those cases.

        The place's passing systems cut it further (GEN.325): the largest orbit whose expected number of ionising passes
        over the star's age is one (`place_terms.planet_limit_au`), a few hundred AU in the Sun's neighbourhood and a few
        AU in a nuclear disc. Comets keep their own ceiling (`validation.orbit_ceiling_au`).

        Args:
            star (Star): The star to compute this boundary for.

        Returns:
            float: The orbit ceiling, in AU.
        """
        terms = getattr(self.system_config, "PLACE_TERMS", None) or place_terms.SOLAR_NEIGHBOURHOOD
        return min(validation.orbit_ceiling_au(star), place_terms.planet_limit_au(terms, getattr(star, "age", None)))

    @staticmethod
    def _engulfment_radius_au(star, explicit_request=False):
        """
        The distance (AU) inside which `star`'s history has removed every
        planet: `GIANT_ENGULFMENT_RADIUS_FACTOR` times the present radius of
        a subgiant, giant or supergiant, and the progenitor's engulfment radius
        (`mass_loss.engulfment_radius_au`, GEN.318) around a white dwarf (whose
        progenitor's giant phases reached that far). A close binary's proxy takes
        the larger of its two stars. 0 for anything else. A white dwarf of a forced
        type (no model progenitor) leaves nothing out when a habitable world or a
        `SLOTS` list is asked for (`explicit_request`): the request names a world
        around it, and its giant phase is no part of what was specified.
        """
        radius_au = 0.0
        for member in getattr(star, "stars", None) or [star]:
            yerkes = getattr(member, "yerkes_class", None)
            if yerkes in ("VII", "D"):
                if explicit_request and member.initial_mass_sol is None:
                    continue
                progenitor_sol = mass_loss.progenitor_mass_sol(_kg_to_msun(member.mass),
                                                               member.initial_mass_sol)
                radius_au = max(radius_au, mass_loss.engulfment_radius_au(progenitor_sol))
            elif yerkes not in (None, "V", "VI", "PMS", "protostar"):
                radius_au = max(radius_au, _km_to_au(tuning.GIANT_ENGULFMENT_RADIUS_FACTOR
                                          * member.radius))
        return radius_au

    def _explicit_request(self):
        """Whether the config asks for a habitable world or a class, or lists its slots."""
        return (self.system_config.HABITABLE_WORLD is True or bool(self.system_config.SLOTS)
                or bool(self.system_config.FORCE_CLASS))

    def _forced_class_present(self):
        """Whether the system holds a planet of `system_config.FORCE_CLASS` (True when none is forced)."""
        letter = self.system_config.FORCE_CLASS
        return not letter or any(getattr(planet, "planet_class", None) == letter for planet in self.planets)

    def _apply_star_history(self, planets, star):
        """
        Removes, in place, the bodies `star`'s age and history rule out: a
        star younger than `PLANET_MIN_STAR_AGE_GY` keeps its belts (a debris
        disk) but no planets yet, and every body whose inner edge lies
        inside `_engulfment_radius_au` is gone. An explicitly required
        habitable world or `SLOTS` list overrides the age rule: a random
        star is never that young when a habitable world is required
        (`stellarEvolution.sample_living_star`), so this only lets a
        specified young type (an O star) keep what was asked for.
        """
        explicit_request = self._explicit_request()
        engulfed_au = self._engulfment_radius_au(star, explicit_request)
        too_young = (star.age is not None and star.age < tuning.PLANET_MIN_STAR_AGE_GY
                     and not explicit_request)
        kept = []
        for body in planets:
            inner_edge = body.lower_limit if body.body_type == 'a' else body.distance
            if inner_edge < engulfed_au or (too_young and body.body_type != 'a'):
                continue
            kept.append(body)
        if len(kept) != len(planets):
            log.debug(f"Removed {len(planets) - len(kept)} bodies ruled out by the star's age or history "
                      f"(age {star.age:.3g} Gyr, engulfed inside {engulfed_au:.3g} AU)")
        planets[:] = kept

    def _trim_to_orbit_ceiling(self, planets, orbit_ceiling_au):
        """
        Removes trailing planets/belts from `planets` that `validate_system`'s
        own overlap correction pushed past `orbit_ceiling_au`.

        `_generate_planets` already keeps every body it places within
        `orbit_ceiling_au` (via `_exceeds_orbit_ceiling`), but
        `validate_system` (called afterward, to resolve spacing collisions
        between neighbors) only ever pushes a body's `distance` further
        OUTWARD, never closer -- for the galactic-Hill-sphere ceiling that
        correction was never practically reachable (negligible next to a
        light-year-scale boundary), but a much tighter S-type `a_crit_au`
        ceiling makes it a real possibility. Trims from the end (outermost
        first) since an outward-only correction cannot newly push an
        *interior* body past the ceiling without its outer neighbor
        already having been pushed past it first.

        Args:
            planets (list): The list to trim in place (`self.planets` or
                `self.secondary_planets`).
            orbit_ceiling_au (float): The boundary to check against, in AU
                (see `_orbit_ceiling_au`).
        """
        validation.trim_to_orbit_ceiling(planets, orbit_ceiling_au)

    def _validate_cross_star_clearance(self):
        """
        For an S-type (wide) binary, prunes the outermost planet from
        whichever star's list gravitationally encroaches on the other
        star's own outermost planet, using a Gladman (1993)-derived mutual
        Hill radius between the two systems' outermost bodies.

        Each star's own Holman & Wiegert (1999) `a_crit_au` (see
        `_orbit_ceiling_au`) only guarantees a *test particle* is long-term
        stable there -- a real generated planet carries its own mass and
        Hill sphere, so the two stars' outermost planets (the ones nearest
        the gap between the stars) can still encroach on each other even
        though each individually respects its own star's `a_crit_au`. This
        is a physically-motivated EXTENSION of Gladman's mutual-Hill-radius
        criterion, not a literal application of it: Gladman's own
        derivation assumes both bodies orbit the SAME central star, but
        here they orbit two different ones. `M_relevant` (Gladman's
        "central mass") is taken as the pair's COMBINED mass
        (`primary_star.mass + secondary_star.mass`) -- the same
        "combined mass governs the pair's local dynamics" reasoning
        `WideBinaryPair` already applies via `mu` in the Holman & Wiegert
        formula, and a reasonable proxy specifically near the gap between
        the two disks, where both stars' gravity is comparably significant
        (that's the entire reason `a_crit_au` exists in the first place).

        Computes the worst-case center-to-center separation between the two
        stars' outermost planets as `the pair's periapsis -
        outer_primary.distance - outer_secondary.distance` -- the conservative case where the two
        planets' independent, uncorrelated `orbital_phase_deg` happen to
        point directly at each other along the binary axis. This is
        deliberately the worst case, not an average: nothing in this
        generator correlates a primary-planet's instantaneous phase against
        a secondary-planet's (each is rolled independently, and
        `planetgen.cli.orbits` advances them independently over time too), so
        it's the only separation value that stays valid across the whole
        orbital cycle, not just at generation time.

        If that worst-case gap is smaller than
        `constants.GLADMAN_MUTUAL_HILL_STABILITY_FACTOR` mutual
        Hill radii, the offending outermost planet is removed -- starting
        with whichever star's outermost planet has the least margin to its
        own `a_crit_au` -- and the check repeats until clearance holds or
        one star's planet list is exhausted of real planets. Framed
        physically as the companion star's gravity "clearing"/"capturing"
        that orbital slot, the same spirit as `_generate_planets`'s own
        documented "runs out of stable room" behavior.

        A star's "outermost" object can be an `AsteroidBelt` as well as a
        `Planet` -- found as whichever object has the greatest outer edge
        (`upper_limit` for a belt, `distance` for a planet), not just the
        list's last element, since the belt-aware placement fix in
        `_distance_avoiding_belts` no longer guarantees a strictly
        distance-sorted list. A belt has no discrete mass/Hill sphere, so
        the Gladman mutual-Hill-radius criterion below only applies when
        BOTH stars' outermost objects are real planets; if either one is a
        belt, this falls back to a plain geometric non-overlap requirement
        (worst-case gap must clear `tuning.
        MIN_ASTEROID_BELT_SEPARATION`, the same fixed buffer
        `validate_system` already uses between same-star neighbors),
        since a belt has no mass to derive a dynamical-stability margin
        from in the first place.

        As the equal-mass, circular-orbit sanity check in
        `orbits.holman_wiegert_critical_semimajor_axis`'s own docstring
        shows (`a_crit_au` sitting at roughly 0.27-0.30 * a_bin for that
        case, comfortably under half of `a_bin`), this should rarely
        trigger for ordinary wide, low-eccentricity pairs -- it exists as a
        safety net for tight/eccentric edge cases and `MAX_PLANETS=True`-
        forced systems, where both stars' lists are pushed all the way out
        to their own `a_crit_au`.

        Called from `__init__`'s retry loop only after both lists reflect
        their own truly final positions -- each star's own `validate_system`
        pass and `_trim_to_orbit_ceiling` trim have already run -- so the
        worst-case gap below is computed from real, final distances rather
        than positions `validate_system` might still push outward.

        Mutates `self.planets`/`self.secondary_planets` in place. No-op if
        `binary_type != "wide"` or either list is empty.
        """
        validation.cross_star_clearance(self)

    def to_dict(self):
        """
        Returns a JSON-serializable dict of the entire generated object
        graph: this system's config, star (polymorphic -- includes nested
        primary/secondary if a P-type/close binary), every planet/belt (in
        orbital order, each recursively including its own moons), and the
        system-level flavor text (see TODO.md's Phase 0 fix, which moved
        this to generation time so it's a plain, idempotent read here).

        `binary_type` (`None`/`"close"`/`"wide"`) is the authoritative
        discriminator going forward; `is_binary` is kept alongside it for
        backward compatibility and now means "this system has two stars",
        true for either configuration (a strict superset of its old,
        P-type-only meaning, so nothing that already reads it as "does
        this system have a companion star" breaks). `secondary_star` and
        `wide_binary` are only non-null when `binary_type == "wide"` (a
        P-type binary's secondary is already nested inside `star`, the
        `BinaryStarProxy` dict). `secondary_planets` is the wide binary
        secondary's own independently-generated planet list -- empty for
        every other configuration.

        Deliberately omits `planet_count`/`belt_count`/`moon_count`/
        `hab_count`/`m_count` (bookkeeping recomputed on load via
        `count_objects`/`count_habitable`, never trusted from disk) and
        `stars`/`primary_star` (resolvable from `star`/`secondary_star`
        alone -- see `from_dict`).

        Note: a system anchored by a `phenomena.compact_remnant.BlackHole`/
        `NeutronStar` (see `StarSystem.__init__`'s `compact_remnant`
        parameter) serializes its `star` dict with that subclass's own
        extra fields intact, but `from_dict` has no `star_kind`
        discriminator to reconstruct it as anything other than a plain
        `Star` -- round-tripping such a system through JSON is not
        currently supported (`phenomenonGen.py`, the only caller that
        creates one, never serializes to JSON -- only `store.py`'s direct,
        non-`to_dict`-based read/write path, via `insert_star_system`/
        `load_star_system`, persists such a system today).

        Returns:
            dict: `schema_version`, `name`, `star_words`, `system_config`, `star`, `is_binary`,
                 `binary_type`, `secondary_star`, `wide_binary`, `planets`,
                 `secondary_planets`, `comets`, `secondary_comets`,
                 `system_flavor_text`.
        """
        is_wide = self.binary_type == "wide"
        return {
            "schema_version": SERIALIZATION_SCHEMA_VERSION,
            "name": self.name,
            "star_words": list(self.star_words) if getattr(self, "star_words", None) else None,
            "system_config": self.system_config.to_dict(),
            "star": self.star.to_dict(),
            "is_binary": self.binary_type is not None,
            "binary_type": self.binary_type,
            "secondary_star": self.secondary_star.to_dict() if is_wide else None,
            "wide_binary": self.wide_binary.to_dict() if is_wide else None,
            "planets": [obj.to_dict() for obj in self.planets],
            "secondary_planets": [obj.to_dict() for obj in self.secondary_planets],
            "comets": [comet.to_dict() for comet in self.comets],
            "secondary_comets": [comet.to_dict() for comet in self.secondary_comets],
            "system_flavor_text": self.system_flavor_text,
            "binary_planetary_wobble_x": self.binary_planetary_wobble_x,
            "binary_planetary_wobble_y": self.binary_planetary_wobble_y,
            "binary_planetary_wobble_z": self.binary_planetary_wobble_z,
            "runaway_class": self.runaway_class,
            "runaway_speed_kms": self.runaway_speed_kms,
            "runaway_origin": self.runaway_origin,
            "runaway_direction": list(self.runaway_direction) if self.runaway_direction is not None else None,
            "orbit_nodes": [dataclasses.asdict(node) for node in self.orbit_nodes],
            "member_stars": [star.to_dict() for star in self.member_stars],
            "member_planets": [[obj.to_dict() for obj in bodies] for bodies in self.member_planets],
            "star_slots": [list(place) for place in self.star_slots] if self.star_slots is not None else None,
            "dark_members": [dataclasses.asdict(member) for member in self.dark_members],
        }

    @classmethod
    def from_dict(cls, data):
        """
        Reconstructs a `StarSystem` from a dict in the shape `to_dict()`
        produces, without re-running any generation (`__init__` is bypassed
        via `object.__new__`) -- this is a pure, faithful replay of
        already-decided data, not a new roll. A seed reproduces a system
        only as part of its galaxy and on the release that made it
        (GEN.39), so the actual decided values are stored and read back
        directly instead.

        This is the single place that resolves both shared back-references
        once and re-attaches the same instances everywhere: `system_config`
        is built once and threaded into every child; `star` is built once
        (via a `binary_type` discriminator choosing `Star.from_dict` vs.
        `BinaryStarProxy.from_dict`) and threaded into every top-level
        `Planet`/`AsteroidBelt` of `self.planets`. A pre-v2 save has no
        `binary_type` key at all -- treated as `"close"` whenever its
        (older, P-type-only) `is_binary` is true, since P-type was the only
        binary configuration that existed before S-type support. For a
        P-type (close) binary, the secondary star is deliberately
        reattached to the *same* shared `system_config` as everything else,
        collapsing the generation-time-only asymmetry where
        `StarSystem.__init__` gives the secondary its own deep-copied
        config (`LARGE_STAR` forced `False`) -- that flag is only ever
        consulted during `generate_star()`, so it has no meaning once a
        star already exists to be reloaded. For an S-type (wide) binary,
        `star` IS the primary (never wrapped), and the secondary is
        reconstructed separately from `data["secondary_star"]`.

        Each item in `planets`/`secondary_planets` is dispatched to
        `AsteroidBelt.from_dict` or `Planet.from_dict` based on its own
        `body_type` (`'a'` vs. `'t'`/`'g'`), the same discriminator
        `StarSystem` itself already uses to tell the two apart at runtime.

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

        binary_type = data.get("binary_type")
        if binary_type is None and data.get("is_binary"):
            binary_type = "close"  # Pre-v2 save: P-type was the only binary configuration that existed.

        if binary_type == "close":
            star = BinaryStarProxy.from_dict(data["star"], system_config)
        else:
            star = Star.from_dict(data["star"], system_config)

        system = object.__new__(cls)
        system.system_config = system_config
        system.star = star
        system.binary_type = binary_type

        if binary_type == "close":
            system.primary_star = star._primary
            system.secondary_star = star._secondary
            system.stars = [star._primary, star._secondary]
            system.wide_binary = None
        elif binary_type == "wide":
            system.primary_star = star
            system.secondary_star = Star.from_dict(data["secondary_star"], system_config)
            system.stars = [system.primary_star, system.secondary_star]
            system.wide_binary = WideBinaryPair.from_dict(
                data["wide_binary"], system_config, system.primary_star, system.secondary_star
            )
        else:
            system.primary_star = star
            system.stars = [star]
            system.wide_binary = None

        def _load_objects(obj_data_list, obj_star):
            objects = []
            for obj_data in obj_data_list:
                if obj_data.get("body_type") == "a":
                    obj = AsteroidBelt.from_dict(obj_data, system_config)
                else:
                    obj = Planet.from_dict(obj_data, obj_star, system_config)
                objects.append(obj)
            return objects

        system.planets = _load_objects(data["planets"], star)
        system.secondary_planets = (
            _load_objects(data.get("secondary_planets", []), system.secondary_star)
            if binary_type == "wide" else []
        )

        system.comets = [Comet.from_dict(c, system_config) for c in data.get("comets", [])]
        system.secondary_comets = (
            [Comet.from_dict(c, system_config) for c in data.get("secondary_comets", [])]
            if binary_type == "wide" else []
        )

        system.orbit_nodes = tuple(OrbitNode(**node) for node in data.get("orbit_nodes", ()))
        system.member_stars = tuple(Star.from_dict(member, system_config) for member in data.get("member_stars", ()))
        system.stars = system.stars + list(system.member_stars)
        slots = data.get("star_slots")
        system.star_slots = tuple(tuple(place) for place in slots) if slots is not None else None
        system.dark_members = tuple(stored_trees.DarkMember(**member) for member in data.get("dark_members", ()))
        system.member_planets = [_load_objects(bodies, member) for member, bodies in zip(
            system.member_stars, data.get("member_planets") or [[] for _ in system.member_stars])]
        system.build_ledger()
        system.system_flavor_text = data.get("system_flavor_text")
        system.runaway_class = data.get("runaway_class")
        system.runaway_speed_kms = data.get("runaway_speed_kms")
        system.runaway_origin = data.get("runaway_origin")
        direction = data.get("runaway_direction")
        system.runaway_direction = tuple(direction) if direction is not None else None
        # Absent from an export written before names/bodies.py -- the star's
        # own name was the system name then.
        system._name = data.get("name") or star.name
        star_words = data.get("star_words")
        system.star_words = tuple(star_words) if star_words else None

        # Schema v20: absent on a pre-v20 save -- 0.0 (no wobble) is the
        # correct default there, same "was implicitly zero before this
        # field existed" reasoning as Star.reflex_offset_x's own docstring.
        system.binary_planetary_wobble_x = data.get("binary_planetary_wobble_x", 0.0)
        system.binary_planetary_wobble_y = data.get("binary_planetary_wobble_y", 0.0)
        system.binary_planetary_wobble_z = data.get("binary_planetary_wobble_z", 0.0)

        system.planet_count, system.belt_count, system.moon_count = system.count_objects()
        system.hab_count, system.m_count = system.count_habitable()
        system.habitable_class_count = system.count_habitable_class_worlds()
        system.habitable_request = system.habitable_request_status()
        system.comet_count = system.count_comets()

        return system

    def generate_slot_object(self, star, slot_spec, estimated_distance, planets=None):
        """
        Builds the celestial object explicitly requested for one orbital slot
        by `system_config.SLOTS`.

        Args:
            star: The star the slot orbits (a belt reads its stored snow line, GEN.442).
            slot_spec (dict): The slot specification, with a required "type"
                              key ("planet" or "asteroid_belt") and optional
                              "planet_class"/"moons" keys (see `SystemConfig.SLOTS`).
            estimated_distance (float): The orbital distance (in AU) computed
                                        for this slot by the generation loop.
            planets (list, optional): Objects already placed earlier in this
                same loop, passed through to `calculate_distance_for_class`
                so an explicit class request doesn't land inside an
                already-placed belt (see `_distance_avoiding_belts`).

        Returns:
            Planet or AsteroidBelt: The generated object for this slot.

        Raises:
            ValueError: If `slot_spec["type"]` is not "planet" or "asteroid_belt", or a planet slot names an unknown class or a class that
                can carry life while `HABITABLE_WORLD` is False.
        """
        slot_type = slot_spec.get("type", "planet")

        if slot_type == "asteroid_belt":
            min_distance = estimated_distance
            max_distance = estimated_distance * draw.uniform(tuning.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MIN, tuning.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MAX)
            return AsteroidBelt(self.system_config, star, estimated_distance, min_distance, max_distance)

        if slot_type == "planet":
            planet_class = slot_spec.get("planet_class")
            distance = self.calculate_distance_for_class(planet_class, estimated_distance, planets=planets)
            return self.build_planet(self.star, self.star.habitable_zone, distance, planet_class=planet_class,
                                     moon_count=slot_spec.get("moons"))

        raise ValueError(f"Invalid slot type '{slot_type}'; expected 'planet' or 'asteroid_belt'.")

    def calculate_distance_for_class(self, planet_class, estimated_distance, planets=None):
        """
        Adjusts an orbital distance for a requested `planet_class` (GEN.688). A class is derived from a finished body's conditions,
        so no zone is named by it; the one request that has a place is a class that can carry life (`HABITABLE_PLANET_CLASSES`),
        which only a world in the ecosphere can be, so it is snapped into the ecosphere (the same way the forced habitable world
        is). Any other class keeps `estimated_distance`, and the body is drawn there (`build_planet`).

        Args:
            planet_class (str or None): The requested planet class, or None if the slot doesn't specify one.
            estimated_distance (float): The orbital distance (in AU) computed for this slot by the generation loop.
            planets (list, optional): Objects already placed earlier in this same loop -- passed through to
                `_distance_within_zone_with_margin` so the adjusted distance doesn't land inside an already-placed belt.

        Returns:
            float: `estimated_distance`, or a distance (in AU) in the ecosphere for a class that can carry life.
        """
        if planet_class not in tuning.HABITABLE_PLANET_CLASSES:
            return estimated_distance
        inner, outer = self.star.habitable_zone
        if inner < estimated_distance < outer:
            return estimated_distance
        return self._distance_within_zone_with_margin(inner, outer, planets=planets)

    def build_planet(self, star, habitable_zone, distance, planet_class=None, life=False, moon_count=None, planned_objects=None):
        """
        A planet at `distance`, drawn again until it is the one a request names (GEN.688): its derived class is `planet_class`, or,
        for `life`, its solved conditions are the ones a habitable request asks for (GEN.916, `habitable_request.meets`: temperate,
        at least 10 kPa of air, a surface ocean), never a class letter. After `HABITABLE_REQUEST_AIM_AFTER_ATTEMPTS` misses each next
        draw is aimed by the temperature the last one solved, inside the zone (GEN.911, `habitable_request.aimed_distance`). Each draw is a
        whole planet (a new mass, make and air from the stream). With no request this is one `Planet`. The last draw is kept when
        `REQUEST_MAX_PLANET_ATTEMPTS` draws miss, and the system's own retry (`MAX_SYSTEM_GENERATION_ATTEMPTS`) then judges the request.

        Args:
            star (Star): The star the planet orbits.
            habitable_zone (tuple): `(inner_au, outer_au)`.
            distance (float): Its orbital distance in AU.
            planet_class (str, optional): The class the finished planet must have.
            life (bool): The finished planet must have a class that can carry life.
            moon_count (int, optional): As `Planet`'s; the moons are added to the planet that is kept, not to every draw.

        Returns:
            Planet: The planet.
        """
        if planet_class is not None and planet_class not in tuning.PLANET_CLASSES:
            raise ValueError(f"Unknown planet class {planet_class!r}; expected one of {''.join(sorted(tuning.PLANET_CLASSES))}.")
        if planet_class in tuning.HABITABLE_PLANET_CLASSES and self.system_config.HABITABLE_WORLD is False:
            raise ValueError(f"A slot asks for a class {planet_class} planet, which can carry life, but HABITABLE_WORLD is False.")
        if planet_class is None and not life:
            return Planet(self.system_config, star, habitable_zone, distance, moon_count=moon_count, planned_objects=planned_objects)
        inner, outer = habitable_zone
        floor_au = self._orbit_floor_au(star)
        if floor_au < outer:
            inner = max(inner, floor_au)
        for attempt in range(tuning.REQUEST_MAX_PLANET_ATTEMPTS):
            # A request for a habitable world is first left to chance; where chance does not give one (a dim star's ecosphere, a heavy
            # one's) the draws that follow are of an Earth-sized planet, the mass at which a world holds its air and its water. A request
            # for a class is left to chance too, and then draws a mass log-uniformly over all the masses a planet has (the architecture's
            # disc makes few ice giants at a cold slot, so chance alone can miss that class for the whole run of draws).
            given = {}
            if life and attempt >= tuning.HABITABLE_REQUEST_FREE_ATTEMPTS:
                given["mass"] = _mearth_to_kg(math.exp(draw.uniform(*(math.log(m) for m in tuning.HABITABLE_REQUEST_MASS_EARTH))))
            elif planet_class is not None and attempt >= tuning.CLASS_REQUEST_FREE_ATTEMPTS:
                given["mass"] = _mearth_to_kg(math.exp(draw.uniform(*(math.log(m) for m in tuning.CLASS_REQUEST_MASS_EARTH))))
            planet = Planet(self.system_config, star, habitable_zone, distance, moon_count=DEFER_MOONS, planned_objects=planned_objects, **given)
            if (planet_class is not None and planet.planet_class == planet_class) or (life and habitable_request.meets_with_moons(planet)):
                break
            draw.redrawn("planet.request_missed")
            if life and attempt + 1 >= tuning.HABITABLE_REQUEST_AIM_AFTER_ATTEMPTS:
                # GEN.911: the next draw goes where the light would give the temperature this one solved the target
                distance = habitable_request.aimed_distance(distance, getattr(planet, "surface_temperature", None), inner, outer)
        else:
            # No draw at this place made what was asked (a white dwarf's ecosphere holds no air): drawing the whole system again
            # at the same star would miss the same way, so the request is left unmet.
            self._request_unreachable = True
        planet.add_moons(moon_count)
        return planet

    def _distance_within_zone_with_margin(self, inner, outer, planets=None):
        """
        Draws a distance uniformly within `[inner, outer]`, leaving a
        safety margin at both edges against `validate_system`'s own
        `MIN_ASTEROID_BELT_SEPARATION` nudge.

        A planet deliberately placed right at a zone boundary (by the
        forced-habitable-world logic, via `_forced_habitable_distance`, or
        an explicit slot request via `calculate_distance_for_class`) could
        otherwise get pushed across that boundary later, when
        `validate_system` enforces minimum separation from a neighboring
        asteroid belt -- which `physics.planets.reconcile_zone`
        would then have to redraw away from the very class this draw
        was placing it as. Reserving this margin up front avoids that in
        the common case; margin is capped at a quarter of the zone's own
        width so a pathologically narrow zone still gets a valid
        (non-empty) range to draw from.

        Args:
            inner (float): Zone's inner bound, in AU.
            outer (float): Zone's outer bound, in AU.
            planets (list, optional): Objects already placed earlier in
                this same orbit-placement loop (`AsteroidBelt`/`Planet`
                mix). When given, any already-placed belt's own
                `[lower_limit, upper_limit]` span is excluded from the
                draw (see `_distance_avoiding_belts`) -- this is what
                keeps a forced-habitable-world/explicit-class placement
                from landing *inside* an asteroid belt that happened to
                be placed earlier at a smaller loop index but a
                numerically closer distance (the loop's `estimated_distance`
                sequence isn't monotonic once a zone-forced placement
                overrides it, so "placed earlier" doesn't imply "closer
                in").  `None` (the default) keeps the old belt-blind
                behavior, e.g. for callers with no `planets` list yet.

        Returns:
            float: A distance in AU, safely inside `[inner, outer]`, and
                  outside every already-placed belt's span when `planets`
                  is given.
        """
        margin = min(tuning.MIN_ASTEROID_BELT_SEPARATION, (outer - inner) / 4)
        lo, hi = inner + margin, outer - margin
        if planets:
            return self._distance_avoiding_belts(lo, hi, planets)
        return draw.uniform(lo, hi)

    def _distance_avoiding_belts(self, lo, hi, planets):
        """
        Draws a distance uniformly within `[lo, hi]`, excluding any
        already-placed `AsteroidBelt`'s own `[lower_limit, upper_limit]`
        span -- the belt-awareness `_distance_within_zone_with_margin` and
        `calculate_distance_for_class` need so a zone-forced planet
        placement can't land inside a belt that already occupies part of
        that zone (see `docs/TODO.md`'s "generated system had two
        overlapping asteroid belts, and a planet orbit inside a belt"
        entry -- this is the second half of that fix; the first half is
        `_validate_cross_star_clearance` no longer skipping belts).

        Falls back to a plain uniform draw across `[lo, hi]` if belts
        already cover the whole range (a rare, pathological case) --
        `validate_system`'s own adjacent-pair spacing correction remains
        the backstop for whatever residual overlap that leaves.

        Args:
            lo (float): Range lower bound, in AU.
            hi (float): Range upper bound, in AU.
            planets (list): Already-placed `AsteroidBelt`/`Planet` objects
                to check against (only belts, `body_type == 'a'`, actually
                constrain anything here).

        Returns:
            float: A distance in AU, inside `[lo, hi]` and outside every
                  belt span in `planets` when that leaves any room at all.
        """
        if hi <= lo:
            return lo

        belt_spans = sorted(
            (max(lo, obj.lower_limit), min(hi, obj.upper_limit))
            for obj in planets
            if obj.body_type == 'a' and obj.lower_limit < hi and obj.upper_limit > lo
        )

        free_spans = []
        cursor = lo
        for span_lo, span_hi in belt_spans:
            if span_lo > cursor:
                free_spans.append((cursor, span_lo))
            cursor = max(cursor, span_hi)
        if cursor < hi:
            free_spans.append((cursor, hi))

        if not free_spans:
            return draw.uniform(lo, hi)

        weights = [span_hi - span_lo for span_lo, span_hi in free_spans]
        pick = draw.uniform(0, sum(weights))
        for (span_lo, span_hi), weight in zip(free_spans, weights):
            if pick <= weight:
                return draw.uniform(span_lo, span_hi)
            pick -= weight
        return draw.uniform(*free_spans[-1])

    def _forced_habitable_distance(self, planets=None):
        """
        Draws a distance for a planet the generation loop is forcing into
        the habitable zone (`HABITABLE_WORLD=True`), via
        `_distance_within_zone_with_margin`.

        Args:
            planets (list, optional): Already-placed objects in this same
                loop -- passed through so the draw avoids landing inside
                an already-placed belt (see
                `_distance_within_zone_with_margin`).

        Returns:
            float: A distance in AU, safely inside the star's
                  `habitable_zone`.
        """
        inner, outer = self.star.habitable_zone
        # A close binary's habitable zone can start inside its circumbinary
        # stability limit; only the stable part of the zone is usable.
        floor_au = self._orbit_floor_au(self.star)
        if floor_au < outer:
            inner = max(inner, floor_au)
        return self._distance_within_zone_with_margin(inner, outer, planets=planets)

    def count_objects(self, planets=None):
        """
        Counts the number of planets, asteroid belts, and moons in the system.

        This method iterates through the list of celestial objects in the system
        and tallies the number of each type. It distinguishes between planets and
        asteroid belts, and also counts the total number of moons orbiting the
        planets. This information is used to provide a summary of the system's
        composition in the `__str__` method.

        Args:
            planets (list, optional): Defaults to `self.planets +
                self.secondary_planets` -- the whole system's objects
                across both stars (for a single star or P-type binary,
                `self.secondary_planets` is always `[]`, so this is
                identical to `self.planets` alone). Pass `self.planets` or
                `self.secondary_planets` explicitly to count just one
                star's own list, e.g. to check a HABITABLE_WORLD/
                ASTEROID_BELT guarantee against the primary alone (see
                `__init__`'s retry loop).

        Returns:
            tuple: A tuple containing the total number of planets, asteroid belts,
                   and moons in the system, in that order.
        """
        if planets is None:
            planets = self.all_planets
        planet_counter, moon_counter, belt_counter = 0, 0, 0
        for planet in planets:
            if planet.body_type == 'a':
                belt_counter += 1
            else:
                planet_counter += 1
                moon_counter += len(planet.moons)
        return planet_counter, belt_counter, moon_counter

    def count_habitable(self, planets=None):
        """
        Counts the habitable worlds in the system: the planets and moons whose overall PHI-6 score is Blue or Green
        (`habitability.is_habitable`; Boss, 2026-10-23), and separately the Class M worlds, which are considered the
        most Earth-like. This is used for the system summary output.

        Args:
            planets (list, optional): See `count_objects`'s identical
                parameter -- defaults to both stars' combined lists.

        Returns:
            tuple: A tuple containing the total number of habitable worlds and
                   the total number of Class M worlds, in that order.
        """
        if planets is None:
            planets = self.all_planets
        hab_count, m_count = 0, 0
        for planet in planets:
            if planet.body_type != 'a':
                for body in (planet, *planet.moons):
                    if habitability.is_habitable(getattr(body, "phi4", None)):
                        hab_count += 1
                    if body.planet_class == "M":
                        m_count += 1
        return hab_count, m_count

    def habitable_request_status(self, planets=None):
        """
        How the system answers a request for a habitable world (`system_config.HABITABLE_WORLD`), `None` when none was made:
        `"aim"` when a planet or moon is habitable (overall PHI-6 Blue or Green), `"fallback"` when the best it has is overall
        Yellow with Suitability at least Yellow (Boss, 2026-10-23: accepted when the request must fall short), else
        `"unmet"`. The generation retry still steers by life class (`count_habitable_class_worlds`), so `"unmet"` is possible
        until the world is regenerated under restricted draws.
        """
        if getattr(getattr(self, "system_config", None), "HABITABLE_WORLD", None) is not True:
            return None
        if planets is None:
            planets = self.all_planets
        best = "unmet"
        for planet in planets:
            if planet.body_type == 'a':
                continue
            for body in (planet, *planet.moons):
                result = habitability.habitable_request_result(
                    getattr(body, "phi4", None), getattr(body, "phi4_suitability", None))
                if result == "aim":
                    return "aim"
                if result == "fallback":
                    best = "fallback"
        return best

    def _settle_to_final_age(self):
        """
        Makes what the star's age touches again, now that it is final (GEN.961, GEN.985, GEN.992; the age only ever rises, once, after the
        planets are in, `Star.adjust_age_for_planets`). For every planet, in order: a giant sized from the radius tables is sized again at
        that age and the mean light of its orbit and its interior and air are drawn again (`planets.settle_giant`); the tides of its orbit
        and its captured moons' are taken again (`planets.settle_orbit`); a rocky planet and the moons of any planet made at a younger age
        have their air, class, spin and what follows them made again (`planets.refinish_rocky`). Every draw is from a stream of the body's
        own, so the run's stream does not move, and a body whose star's age did not move is left as drawn. It runs before the forced
        habitable world, belt and class are judged, so they are judged on the finished classes.
        """
        with log.timed_phase("settle pass"):
            for obj in self.all_planets:
                if obj.body_type == 'a':  # a belt is not a body
                    continue
                if obj.body_type == 'g' and not obj.is_moon and planetPhysics.settle_giant(obj):
                    log.debug(f"{obj.name}: settled to {obj.radius:.0f} km at the star's age of {obj.star.age:.3g} Gyr")
                for body in (obj, *obj.moons):
                    if planetPhysics.settle_orbit(body):
                        log.debug(f"{body.name}: eccentricity settled to {body.eccentricity:.4f} at the star's age of {obj.star.age:.3g} Gyr")
                made = planetPhysics.refinish_rocky(obj)
                if made:
                    log.debug(f"{obj.name}: {made} bod{'y' if made == 1 else 'ies'} made again at the star's age of {obj.star.age:.3g} Gyr")

    def count_habitable_class_worlds(self, planets=None):
        """
        Counts the planets and moons whose class can carry life (`tuning.HABITABLE_PLANET_CLASSES`). A report column only: a request
        for a habitable world steers by conditions (GEN.916, `habitable_request.count`), and `count_habitable` says how a body scored.

        Args:
            planets (list, optional): See `count_objects`'s identical parameter.

        Returns:
            int: The number of planets and moons of a life class.
        """
        if planets is None:
            planets = self.all_planets
        return sum(1 for planet in planets if planet.body_type != 'a'
                   for body in (planet, *planet.moons) if body.planet_class in tuning.HABITABLE_PLANET_CLASSES)

    def count_comets(self, comets=None):
        """
        Counts native comets in the system (see `_generate_comets`).

        Kept as its own method, separate from `count_objects`, since
        `self.comets`/`self.secondary_comets` is its own list, not part of
        `self.planets`/`self.secondary_planets` -- see `_generate_comets`'s
        own docstring for why a comet doesn't share that list.

        Args:
            comets (list, optional): Defaults to `self.comets +
                self.secondary_comets` -- both stars' combined comet
                population (for a single star or P-type binary,
                `self.secondary_comets` is always `[]`).

        Returns:
            int: The total number of comets.
        """
        if comets is None:
            comets = self.comets + self.secondary_comets
        return len(comets)

    def _plan_architecture(self, star, orbit_ceiling_au):
        """
        The planets `star` is to hold and how many bodies the system is to have (GEN.728, GEN.729, GEN.730, GEN.828, GEN.829, GEN.920,
        `physics/planet_architecture.py`): an inner system drawn at the observed share, a chain of oligarchs thinned by dynamical instability,
        gas giants and migration, all from the star and its disc. The planets of the chain are the system's slots, in order, each with its
        distance and mass; a system asked to hold more has the next oligarchs out added (`planet_architecture.Plan.seed`).

        `PLANETS` False gives no object, `NUM_ORBITS` exactly that many (the first of the chain, then more oligarchs). `MAX_PLANETS` True
        gives the star an inner system and a starting count two widths of its spread above the mean (the most its disc holds), False the
        fewest (none, or one when `PLANETS` is True); a `max_planets` prevalence (GEN.52) scales the chance of the most. Planets inside
        the star's orbit floor (`_orbit_floor_au`) are left out.

        Args:
            star (Star): The star -- `self.star` for a single star or P-type binary, or individually `self.primary_star` or
                `self.secondary_star` for an S-type binary (see `_generate_planets`, this method's only caller).
            orbit_ceiling_au (float): The orbit ceiling of the star (`_orbit_ceiling_au`).

        Returns:
            tuple: The `planet_architecture.Plan` and the number of bodies the system holds.
        """
        config = self.system_config
        fe_h = getattr(star, "metallicity_fe_h", None)
        fe_h = 0.0 if fe_h is None else fe_h
        metal_ratio = bulk_draw.metal_ratio_to_solar(star)
        disk_edge_au = self._disk_outer_edge_au(star)
        if math.isinf(disk_edge_au):
            disk_edge_au = tuning.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER * star.snow_line_au
        outer_au = min(orbit_ceiling_au, disk_edge_au)

        sublimation_au = planet_architecture.sublimation_au(star.luminosity or 0.0)

        def plan_of(architecture, limit=None):
            floor_au = self._orbit_floor_au(star)
            if floor_au > 0.0:
                architecture = architecture._replace(planets=tuple(s for s in architecture.planets if s.distance_au >= floor_au))
            return planet_architecture.Plan(architecture, star.mass, star.snow_line_au, metal_ratio, fe_h, limit, sublimation_au)

        def built(maximum=False):
            return planet_architecture.build(star.mass, star.snow_line_au, metal_ratio, fe_h, outer_au, maximum=maximum, floor_au=sublimation_au)

        if config.PLANETS is False:
            return plan_of(planet_architecture.Architecture((), 1.0, tuning.OLIGARCH_SEPARATION_MEDIAN, 0, planet_architecture.Lost(0.0, 0.0, 0, 0.0, 0, (), (), (), 0.0, 0.0)), 0), 0
        minimum = 1 if config.PLANETS is True else 0
        if config.NUM_ORBITS is not None:
            return plan_of(built(), config.NUM_ORBITS), config.NUM_ORBITS
        if config.MAX_PLANETS is False:
            return plan_of(built(), minimum), minimum
        architecture = built(maximum=True) if config.MAX_PLANETS is True else built()
        if config.MAX_PLANETS is not True and prevalence.percent(config, "max_planets"):
            if draw.random() < prevalence.scaled_chance(config, "max_planets", 1.0 / (len(architecture.planets) + 1)):
                architecture = built(maximum=True)
        plan = plan_of(architecture)
        return plan, max(len(plan), minimum)

    def validate_system(self, planets=None):
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
        *mutual* Hill radius (`validation.min_distance_after_au`), not either
        one's own individual Hill radius alone -- see
        `tuning.MUTUAL_HILL_RADII_SEPARATION`'s docstring for
        why. A belt has no mass/Hill-radius concept of its own, so any
        correction involving one uses the single real planet's own
        `min_orbit_distance` (5 Hill radii) on whichever side of the belt
        the planet is (`validation.min_distance_after_au` for a planet after a
        belt), or the fixed `MIN_ASTEROID_BELT_SEPARATION` between two
        belts.
        If an adjustment is made to a planet (as opposed to an asteroid belt,
        which carries no class/climate of its own), `_reconcile_moved_planet`
        re-derives its zone from the corrected distance and makes its climate and class again there, along with the same
        treatment for each of its moons: the class is derived from the conditions of the corrected position.

        Args:
            planets (list, optional): Defaults to `self.planets`. Pass
                `self.secondary_planets` to separately validate an S-type
                binary's secondary list -- the two lists orbit different
                stars and are spaced/validated independently; a *cross*-star
                clearance check (which does matter) is handled separately
                by `_validate_cross_star_clearance`, called before either
                list reaches this method (see `__init__`'s retry loop).
        """
        if planets is None:
            planets = self.planets
        validation.space_orbits(planets)

    def _reconcile_moved_planet(self, planet):
        """
        Call right after `validate_system` moves a top-level planet's `distance` to resolve an orbital overlap: its zone, climate
        and class (derived from its conditions) follow the corrected position, and so do its moons' (`validation.reconcile_moved_planet`).
        """
        validation.reconcile_moved_planet(planet)

from planetgen.generation import object_timing  # noqa: E402 -- registers the generators below (PERF.92)

_kg_to_msun = units.converter("kg", "Msun")
_mearth_to_kg = units.converter("Mearth", "kg")
object_timing.track(StarSystem, "system")
