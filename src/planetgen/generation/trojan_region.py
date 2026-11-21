# planetgen/generation/trojan_region.py

"""
The trojan swarm of a giant, as the regions at its L4 and L5 (GEN.1150, docs/design/systems/lagrange-trojans.md)
=================================================================================================================

A giant that swept an outer disk (`belt_birth.hand_on`) holds the part of what it swept that the swept split gives the 1:1 zone (`swept_hand_on.trojan_part`, the scatter's `trojans` row, less what its
neighbours leave it, GEN.1455). `build` makes that mass the **leading** (L4) and **trailing** (L5) `TrojanRegion` of the giant (`trojan_swarm.draw_regions`: the split between the two, the libration of
each region's bodies, its mean eccentricity and inclination). A region is a belt-like object: the mass the star's age of collisions leaves of its share, its largest body, its composition (the
belt's, drawn at the giant's distance), and a stretch of orbit radii about the giant's own (the region's mean eccentricity either side). Nothing here sets how often a giant has a swarm; a giant
that swept no outer disk has none, and a swarm under `belt_birth.MIN_PRESENT_MASS_EARTH` is not kept (as a belt of that little is not).

A region is made once the system is final (after the spacing pass, the irregular moons and the comets), so the pass never sees a body sharing a giant's orbit, and it is placed in the list of its host.
It has no members in this slice. It is stored as a row of `asteroid_belts` whose parent is its host planet (`parent_type` `planet`), with the columns `lagrange_point`, `libration_amplitude_deg`,
`libration_phase_deg`, `mean_eccentricity` and `mean_inclination_deg`.
"""

from planetgen.generation import belt_birth
from planetgen.generation.belt import AsteroidBelt, _STARS, lit_luminosity_w, present_state, star_al26_ratio
from planetgen.physics import bulk_draw, small_body_composition, trojan_swarm
from planetgen.physics import unit_registry as units
from planetgen.util import draw, log

ORIGIN = "trojan-swarm"
"""str: The `origin` of a trojan region (`asteroid_belts.origin_cause`)."""


class TrojanRegion(AsteroidBelt):
    """
    One region of a giant's trojan swarm: a belt-like object at the giant's L4 or L5.

    Attributes beyond a belt's:
        host (Planet or None): The giant it librates about. Not stored: the row's parent link names it, and a region read back is given it by the loader.
        lagrange_point (str): `"L4"` (leading) or `"L5"` (trailing).
        libration_amplitude_deg, libration_phase_deg (float): The libration of the region's bodies about the point (`trojan_swarm.longitude_offset_deg`).
        mean_eccentricity, mean_inclination_deg (float): The region's mean eccentricity and its mean inclination to the giant's plane.
    """

    SERIALIZABLE_FIELDS = [*AsteroidBelt.SERIALIZABLE_FIELDS, "lagrange_point", "libration_amplitude_deg", "libration_phase_deg", "mean_eccentricity", "mean_inclination_deg"]
    """Every attribute a belt keeps and the five of the region."""

    host = None
    host_db_id = None
    """int or None: The `planets.id` of the host a region read back from a store names as its parent, until the loader gives it the planet. Not stored here."""

    def __init__(self, system_config, star, host, region):
        """
        Args:
            system_config (SystemConfig): The system's shared config.
            star: The star the host orbits.
            host (Planet): The giant.
            region (trojan_swarm.Region): The region: its point, its mass (kg), its libration, eccentricity and inclination.
        """
        self.system_config = system_config
        self.host = host
        self.body_type = 'a'
        self.distance = host.distance
        self.lower_limit = host.distance * (1.0 - region.eccentricity)
        self.upper_limit = host.distance * (1.0 + region.eccentricity)
        self.origin = ORIGIN
        self.lagrange_point = region.point
        self.libration_amplitude_deg, self.libration_phase_deg = region.libration_amplitude_deg, region.libration_phase_deg
        self.mean_eccentricity, self.mean_inclination_deg = region.eccentricity, region.inclination_deg
        self.primordial_mass_earth = units.convert(region.mass, "kg", "Mearth")
        self.cleared_mass_earth = self.scattered_mass_earth = 0.0
        self.co_causes = []
        self.origin_history = []
        self.origin_sentence = (f"The {'leading' if region.point == 'L4' else 'trailing'} trojan region of the giant at {host.distance:.3g} AU: what its sweep of the outer disk left at its {region.point} point.")
        _STARS[self] = star
        self.largest_body_km, self.mass_kg, self.density = present_state(star, self.distance, self.lower_limit, self.upper_limit, self.primordial_mass_earth)
        log.choice("Trojan region", region.point, f"{self.mass_kg:.4g} kg left of {region.mass:.4g} kg after {star.age:.4g} Gy, libration {region.libration_amplitude_deg:.3g} degrees")
        with draw.bound(draw.Stream(f"trojan-region-makeup:{host.distance!r}:{host.mass!r}:{star.mass!r}:{region.point}")):        # the region's own stream: the draws of what follows do not move with it
            self.makeup = small_body_composition.draw_belt(self.distance, lit_luminosity_w(star, self.lower_limit, self.upper_limit), star_al26_ratio(system_config, star), host=bulk_draw.host_ratios(star))
        self.composition = small_body_composition.describe(self.makeup)
        self.chondrite_group = self.makeup.detail["group"]

    def _left_earth(self):
        """The mass the region grinds from, Earth masses: its share of the swarm."""
        return self.primordial_mass_earth

    def _draw_members(self):
        """A region has no members in this slice (GEN.1150): it holds a mass, a makeup and a largest body."""
        return []


def build(bodies, star, system_config):
    """
    The `TrojanRegion`s of the giants among `bodies` (one star's list, final): for each outer disk the outermost giant swept (`belt_birth.hand_on` left a `swarm_kg`), the two regions of the swarm that giant holds, nearest the
    star first, leading before trailing. A region under `belt_birth.MIN_PRESENT_MASS_EARTH` once the star's age has ground it is not kept. The draw is the swarm's own stream, keyed by the giant (its distance,
    its mass and the star's), so no other draw moves.
    """
    host_mass = star.mass
    found = []
    giants = [body for body in bodies if getattr(body, "body_type", None) == 'g']
    for disk in (body for body in bodies if getattr(body, "body_type", None) == 'a' and getattr(body, "swarm_kg", 0.0) > 0.0):
        if not giants:
            continue
        host = max(giants, key=lambda giant: giant.distance)         # the scatter step served the outermost giant (`belt_birth._fates`)
        stream = draw.Stream(f"trojan-swarm:{host.distance!r}:{host.mass!r}:{host_mass!r}")
        for region in trojan_swarm.draw_regions(host_mass, host.mass, disk.swarm_kg, stream):
            made = TrojanRegion(system_config, star, host, region)
            if units.convert(made.mass_kg, "kg", "Mearth") >= belt_birth.MIN_PRESENT_MASS_EARTH:
                found.append(made)
    return sorted(found, key=lambda made: (made.distance, made.lagrange_point))
