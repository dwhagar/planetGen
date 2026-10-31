# planetgen/generation/stellar_object.py

"""
Stellar Object
==============

The family class of everything that shines, or once shone, as a star or was made by one's collapse (GEN.565d; the
`StellarObject` of docs/design/architecture/mass-object-class.md, section 10): what a main-sequence star, a giant, a white dwarf,
a binary's merged proxy, a neutron star and a black hole share. It holds the persisted field list, the class-level
defaults a reloaded object needs, the save and load of its two non-JSON fields, the age rules that keep a star's
class true to its planets, and the system perimeter and galactic orbit every one of them has. `Star` adds the
procedural draw of a main-sequence, giant or dwarf star on top of it; `CompactRemnant` and `BinaryStarProxy` fill in
their own properties and used to do it by telling `Star.__init__` to skip its own (`_skip_property_init`), a
workaround that went with this class.

Nothing here draws a random number except `MassObject.set_galactic_motion`, where each family's constructor calls
it.
"""

from planetgen import tuning
from planetgen.generation.mass_object import GALACTIC_MOTION_UNITS, MassObject
from planetgen.galaxy import enclosed_mass, scales
from planetgen.galaxy.galactic_orbit import calculate_galactic_orbit
from planetgen.names.wordsalad import LazySaladName
from planetgen.physics import activity, age_needs, constants, pre_main_sequence, spin
from planetgen.physics.orbits import calculate_hill_sphere
from planetgen.util import draw
from planetgen.physics import unit_registry as units

_ly_to_m = units.converter("ly", "m")
_m_to_pc = units.converter("m", "pc")
_m_to_au = units.converter("m", "au")


class StellarObject(LazySaladName, MassObject):
    """
    A star-like object: its name, its mass and radius, its luminosity, temperature and type, its age and lifespan,
    its habitable zone and system perimeter, its galactic orbit, its spin and activity and its [Fe/H]. A subclass
    sets them: `Star` by drawing a star, `CompactRemnant` by its collapse, `BinaryStarProxy` from the pair it stands
    for.
    """

    def __init__(self, system_config, name=None, galactic_center_dist_ly=None):
        """
        Args:
            system_config (SystemConfig): The shared SystemConfig object for the system.
            name (str, optional): The name to give; one is drawn on first read when omitted (PERF.43).
            galactic_center_dist_ly (float, optional): This object's distance from the galactic centre, in
                light-years, threaded into `calculate_system_perimeter`; `None` uses the radius that holds half of
                the galaxy's disc stars (`galaxy.scales`).
        """
        self.system_config = system_config
        self._init_salad_name(name)  # drawn on first read (PERF.43)
        self.galactic_center_dist_ly = galactic_center_dist_ly

    FIELD_UNITS = {
        "name": "text", "type": "text", "yerkes_class": "text", "mass": "kg", "radius": "km", "temperature": "K",
        "luminosity": "W", "age": "Gyr", "lifespan": "Gyr", "habitable_zone": "au", "system_perimeter": "au",
        "heliosphere_radius": "au", **GALACTIC_MOTION_UNITS,
        "a_crit_au": "au", "reflex_offset_x": "au", "reflex_offset_y": "au", "reflex_offset_z": "au",
        "rotation_period_hours": "h", **spin.SPIN_UNITS, **activity.STAR_ACTIVITY_UNITS, "metallicity_fe_h": "dex",
        "snow_line_au": "au", **pre_main_sequence.DISK_UNITS,
    }
    """
    Every attribute set by `__init__`/`generate_star`, excluding `system_config` (a shared back-reference threaded
    into `from_dict` rather than serialized redundantly on every star), with its native unit (GEN.570). The habitable
    zone is a pair (its inner and outer edge) of AU.
    """

    SERIALIZABLE_FIELDS = list(FIELD_UNITS)
    """The names of `FIELD_UNITS`, in the order they persist."""

    metallicity_fe_h = None
    """float or None: [Fe/H] in dex (GEN.158); class-level default for the same reason as `a_crit_au`."""

    snow_line_au = None
    """float: This star's own snow line in AU (GEN.458), drawn once around the blackbody position for its luminosity
    (`formation.draw_snow_line_au`) and stored; the water share, the small-body zones and the moons read it. Set by
    the end of every star's and compact remnant's `__init__`; None only before that."""

    disk_kind = None
    disk_mass_solar = None
    disk_inner_au = None
    disk_outer_au = None
    """The disk descriptor of a protostar or a young pre-main-sequence star (GEN.242, `pre_main_sequence.draw_disk`):
    `disk_kind` ("class-0", "class-i" or "class-ii"), its mass in solar masses and its inner and outer edge in AU;
    None for every other star. A descriptor only: no body is placed in it. Class-level defaults for the same reason as
    `a_crit_au`."""

    a_crit_au = None
    """
    float or None: Class-level default (rather than only an `__init__`
    instance assignment) so a `Star` reconstructed via `from_dict`'s
    `object.__new__`+`fields_from_dict` path -- which leaves a field
    untouched when its key is absent from the loaded data (see
    `serialization.fields_from_dict`) -- still reads back `None` for any
    system saved before this field existed, rather than raising
    `AttributeError`. Populated only for a star that is one constituent of
    an S-type (wide) binary: this star's own Holman & Wiegert (1999)
    critical semi-major axis (see
    `orbits.holman_wiegert_critical_semimajor_axis`), the maximum orbit
    distance that stays long-term stable against its companion's
    perturbation. `None` for a single star or either star of a P-type
    (close) binary, where no such companion-driven limit applies (a P-type
    pair's planets orbit the merged `BinaryStarProxy`, which never sets
    this at all). Set directly by `StarSystem.__init__` after both stars
    of a wide pair exist -- a `Star` never computes this on its own, since
    it depends on the companion star's mass/separation/eccentricity, none
    of which a lone `Star` instance knows about.
    """

    initial_mass_sol = None
    """float or None: The initial (zero-age) mass in Msun of a star drawn
    from the population model (`_generate_from_population_model`); None for
    a specified-type or reloaded star. A binary secondary's mass is a
    fraction of the primary's initial mass, since the primary may have lost
    mass since (a white dwarf)."""

    phase_end_age_gy = None
    """float or None: The age (Gy) at which a population-model star leaves
    its present evolutionary phase; `adjust_age_for_planets` keeps its age
    below it so the star's class stays true. None otherwise."""

    _model_age_lifespan = None
    _model_fe_h = None
    _model_rotation_hours = None

    reflex_offset_x = reflex_offset_y = reflex_offset_z = 0.0
    """
    float: Class-level defaults (schema v20), same "still readable after
    `from_dict`'s `object.__new__`+`fields_from_dict` path" reasoning as
    `a_crit_au` above. This star's own displacement (AU) from its nominal
    fixed point, caused by the combined gravitational pull of every planet
    orbiting it directly -- the "wobble" half of a proper two-body
    treatment, alongside the planets' own unchanged `position_x/y/z` (see
    `orbits.calculate_reflex_offset`'s docstring). Zero for a star with no
    planets. Set directly by `StarSystem.__init__` once `self.planets`
    exists -- a `Star` never computes this on its own construction, the
    same reasoning `a_crit_au` documents above. Applies identically to a
    `BlackHole`/`NeutronStar` (both subclass `Star`), consistent with
    `store.advance_orbital_phases`' class-blind bulk `UPDATE`s.
    """

    def _encode(self, data):
        """`habitable_zone` is stored as a list (JSON has no tuple type). `lifespan` is stored as `None` in place of
        `float('inf')` (white dwarfs) -- the same convention `planetgen/db/store.py` uses via `_lifespan_gy` --
        since `float('inf')` round-trips through the standard `json` module as a non-standard `Infinity` token that
        not every JSON consumer accepts."""
        data["habitable_zone"] = list(self.habitable_zone)
        data["lifespan"] = None if self.lifespan == float('inf') else self.lifespan

    def _decode(self, data):
        self.habitable_zone = tuple(data["habitable_zone"])
        self.lifespan = float('inf') if data["lifespan"] is None else data["lifespan"]


    def planet_age_requirement_gy(self, planets):
        """
        The youngest age (Gy) this star can have and still host `planets`:
        the highest minimum age any planet's stored conditions need on one of
        this star's `supported_evolutionary_scales` (`physics.age_needs`, GEN.686).
        0.0 when none needs any. See `adjust_age_for_planets`.
        """
        supported_scales = supported_evolutionary_scales(self)
        required = 0.0
        for planet in planets:
            if hasattr(planet, 'planet_class'):
                required = max(required, age_needs.minimum_age_gy(planet, supported_scales))
        return required

    def age_ceiling_gy(self):
        """
        The oldest (Gy) this star can be made by `adjust_age_for_planets`
        without contradicting what was generated for it, never below its
        present age: a population-model star stays inside its present phase
        (a white dwarf's age sets its cooling, so it can't move at all); a
        specified-type star stays inside its lifespan. Both are capped at
        the age of the universe. A binary's companion passes its ceiling to
        the star that sets the pair's shared age (see
        `adjust_pair_age_for_planets`).
        """
        if self.phase_end_age_gy == float('inf'):
            return self.age
        if self.phase_end_age_gy is not None:
            limit = self.phase_end_age_gy
        elif self.lifespan != float('inf'):
            limit = self.lifespan
        else:
            return float('inf')
        return max(self.age, min(limit * tuning.MAX_PLANET_AGE_ADJUSTMENT_FACTOR,
                                 tuning.UNIVERSE_AGE_GY))

    def adjust_age_for_planets(self, planets, min_required_age_gy=0.0, age_ceiling_gy=None):
        """
        Adjusts the star's age to ensure it is at least as old as the minimum
        age required for the most "evolved" planet in its system. This method
        should be called by the `StarSystem` class after planets have been generated.

        The star's age is a crucial factor in planetary evolution. This function
        iterates through all generated planets, identifies their minimum required
        ages based on their `planet_class` and the star's `supported_evolutionary_scales`,
        and then ensures the star's age meets or exceeds the highest of these
        minimum requirements. This prevents scenarios where a young star might
        host an "old" planet (e.g., a terraformed world requiring billions of years
        of development).

        If the star's lifespan is too short to accommodate the required
        planetary age even at its very end, the star's age is set to be near
        the end of its lifespan (the closest it can get), reflecting a system
        that has evolved as rapidly as physically possible for its type. If
        the requirement *is* reachable within the lifespan, `config.AGE`
        ("young"/"old") still biases where in the achievable window the age
        falls, the same way it biases the star's initial age (the population
        model's `AGE` window) — a forced minimum age
        should narrow the age range, not silently discard a young/old
        request entirely.

        Args:
            planets (list): A list of `Planet` objects in the system.
            min_required_age_gy (float): A further minimum age, from planets
                judged by another star (a binary companion's own).
            age_ceiling_gy (float, optional): A further maximum age (a
                binary companion's `age_ceiling_gy`); the age then also
                never drops below its present value, since the companion was
                born at it.
        """
        original_age = self.age
        if self.phase_end_age_gy == float('inf'):
            # A population-model white dwarf's age sets its cooling
            # luminosity and temperature; it stays as drawn.
            return

        min_required_age_for_system = max(self.planet_age_requirement_gy(planets), min_required_age_gy)

        # Ensure the star's age is at least the minimum required by its planets
        if self.age < min_required_age_for_system:
            if self.phase_end_age_gy is not None and self.phase_end_age_gy != float('inf'):
                # A population-model star can't age past its present phase
                # without its class, luminosity and radius becoming wrong.
                max_reachable_age = min(
                    self.phase_end_age_gy * tuning.MAX_PLANET_AGE_ADJUSTMENT_FACTOR,
                    tuning.UNIVERSE_AGE_GY,
                )
            elif self.lifespan != float('inf'):
                # Capped at the age of the universe, since a long-lived class's lifespan-derived ceiling can be far
                # larger than any star could actually have reached yet.
                max_reachable_age = min(
                    self.lifespan * tuning.MAX_PLANET_AGE_ADJUSTMENT_FACTOR,
                    tuning.UNIVERSE_AGE_GY,
                )
            else:
                max_reachable_age = min_required_age_for_system + tuning.WHITE_DWARF_AGE_ADDITION_GY # Add 5 GY for WD if no upper bound
            if age_ceiling_gy is not None:
                max_reachable_age = max(min(max_reachable_age, age_ceiling_gy), original_age)

            if max_reachable_age < min_required_age_for_system:
                # The requirement is physically unreachable within this star's lifespan
                # (e.g. a habitable world's minimum age exceeds a short-lived O/B/A star's
                # entire life). The best it can do is sit near the very end of its life --
                # there's no real "young" vs "old" choice left to honor here.
                # `min(..., max_reachable_age)` keeps this "near the end of its
                # reachable life" for a long-lived class too: lifespan * 0.85
                # can itself be far larger than max_reachable_age once that's
                # capped at UNIVERSE_AGE_GY, which would otherwise widen this
                # into a near-full-lifespan random range instead of "near the end."
                self.age = draw.uniform(
                    min(self.lifespan * tuning.UNREACHABLE_PLANET_AGE_MIN_LIFESPAN_RATIO, max_reachable_age),
                    max_reachable_age,
                )
                if self.phase_end_age_gy is not None or age_ceiling_gy is not None:
                    # Never younger than it was: that could move a
                    # population-model star (or a binary companion born at
                    # this age) back into an earlier phase.
                    self.age = max(self.age, original_age)
            else:
                # The requirement is reachable within the lifespan: bias where in the
                # achievable [min_required_age_for_system, max_reachable_age] window the
                # age falls, mirroring the young/old bias used for the initial age roll.
                low, high = min_required_age_for_system, max_reachable_age
                if self.system_config.AGE == "old":
                    low = low + (high - low) * tuning.OLD_STAR_AGE_LIFESPAN_RATIO
                elif self.system_config.AGE == "young":
                    high = high - (high - low) * (1 - tuning.YOUNG_STAR_AGE_LIFESPAN_RATIO)
                self.age = draw.uniform(low, high)

        # Ensure age doesn't exceed lifespan (unless lifespan is infinite)
        if self.lifespan != float('inf') and self.age >= self.lifespan:
            self.age = self.lifespan * tuning.MAX_PLANET_AGE_ADJUSTMENT_FACTOR # Star is near end of life
        if age_ceiling_gy is not None:
            self.age = max(min(self.age, age_ceiling_gy), original_age)

    def calculate_system_perimeter(self, galactic_center_dist_ly=None):
        """
        Calculates the Hill sphere for the star system relative to the galaxy.

        The Hill sphere defines the gravitational sphere of influence of the star
        within the larger galactic environment. Any object orbiting the star
        beyond this radius would be more strongly perturbed by the galaxy's
        gravity than by the star's, making it unlikely to be gravitationally
        bound to the star system.

        This method uses the star's mass, the mass of its galaxy inside the star's own distance from the
        galactic center (`galaxy.enclosed_mass`, GEN.370) and that distance to compute this boundary.
        The result is converted from meters to Astronomical Units (AU) for
        convenience and consistency with other system-scale measurements.

        Args:
            galactic_center_dist_ly (float, optional): This star's actual
                distance from the galactic center, in light-years -- derived
                from its owning sector's `galactic_radius_pc` once sectors
                are placed in galaxy-space (see
                `docs/design/galaxy/galaxy-coordinate-system.md` section 5). `None`
                (the default) falls back to the radius that holds half of the galaxy's disc stars (`galaxy.scales`) --
                the correct behavior for a system generated outside any
                galaxy context (e.g. `generation/run_sector.py`'s own standalone CLI,
                whose sectors are never placed in a galaxy).

        Returns:
            float: The radius of the Hill sphere in Astronomical Units (AU).
        """
        if galactic_center_dist_ly is None:
            galactic_center_dist_ly = scales.default_galactocentric_ly()
        galactic_center_dist_m = _ly_to_m(max(galactic_center_dist_ly, constants.MIN_TIDAL_DISTANCE_LY))
        hill_radius_m = calculate_hill_sphere(galactic_center_dist_m, self.mass, enclosed_mass.current_enclosed_mass_kg(_m_to_pc(galactic_center_dist_m)))
        return _m_to_au(hill_radius_m)

    def calculate_galactic_orbit(self, galactic_center_dist_ly=None):
        """
        Calculates this star system's circular orbital speed and period
        around the galactic center, via `galactic_orbit.calculate_galactic_orbit`.

        Unlike `calculate_system_perimeter`, the result doesn't depend on
        this star's own mass -- see that function's docstring -- so this
        method exists purely to thread the star's `galactic_center_dist_ly`
        fallback convention through, the same way `calculate_system_perimeter`
        does.

        Args:
            galactic_center_dist_ly (float, optional): This system's actual
                distance from the galactic center, in light-years -- see
                `calculate_system_perimeter`'s docstring for the same
                parameter. `None` (the default) uses the radius that holds half of the galaxy's disc stars (`galaxy.scales`).

        Returns:
            tuple: `(orbital_speed_kms, orbital_period_gy)` -- see
                  `galactic_orbit.calculate_galactic_orbit`.
        """
        if galactic_center_dist_ly is None:
            galactic_center_dist_ly = scales.default_galactocentric_ly()
        return calculate_galactic_orbit(galactic_center_dist_ly)


SPECTRAL_LETTERS = "OBAFGKM"
"""str: The spectral letters of a star whose life the paces are judged for; a remnant or a substellar object (any other letter) has none."""


def supported_evolutionary_scales(star):
    """
    The evolutionary paces (`tuning.EVOLUTIONARY_TIMELINES`: `fast`, `normal`, `slow`) a planet of `star` can follow to a technological
    civilization, from the time the star gives (GEN.942): the paces whose last milestone comes at or before the end of the star's own
    lifespan, the star model's, not a table of lifespans by spectral letter. A white dwarf's lifespan is infinite (it just cools forever),
    so the time it has had so far, its age, stands in. A star that gives less time than the fastest pace needs still has `fast`, as
    `get_evolutionary_state` does. A remnant or substellar object (a letter outside `SPECTRAL_LETTERS`) has none.

    Args:
        star: A `Star` or `BinaryStarProxy` instance.

    Returns:
        list[str]: The paces, fastest first; empty for a star with no spectral letter of the main sequence.
    """
    reference_star = star._primary if hasattr(star, '_primary') else star
    letter = reference_star.type[0].upper()
    if letter not in SPECTRAL_LETTERS:
        return []
    time_budget = reference_star.age if reference_star.lifespan == float('inf') else reference_star.lifespan
    reachable = [scale for scale in ["fast", "normal", "slow"]
                 if tuning.EVOLUTIONARY_TIMELINES[scale]['technological_civilization'] <= time_budget]
    return reachable or ["fast"]
