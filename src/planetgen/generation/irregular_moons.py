# planetgen/generation/irregular_moons.py

"""
The irregular moons of a planet (GEN.1083, docs/design/systems/captured-asteroid-moons.md sections 4 and 6)
===========================================================================================================

A second tier of moons above the regular ones. A regular moon is drawn by `physics.planets.generate_moons` from the planet's own formation (its mass is at least `tuning.MOON_MIN_MASS_KG`);
an irregular moon is a small body the planet caught later, in a window, from the outer disk of its system or from its own feeding zone (`physics.capture_windows`,
`physics.irregular_population`). It needs what the regular moons are built before and the belts are built after (the outer disk, the instabilities of the planet architecture), so
`StarSystem` builds the irregulars of a star's planets once its belts stand (`for_planet`), and a planet keeps them in `planet.irregulars`.

For one planet:

1. The zone: Hill radius, stable limits and the inner edge from the outermost regular moon (`irregular_zone`). No room, no irregulars.
2. The pool: what the encounter window and the gas window caught of the outer disk's primordial mass (`irregular_population.make_pool`), ground down for the age of the star
   (`irregular_population.draw_irregulars`), and the one massive body binary exchange may have brought (`draw_exchange`).
3. For each stored body, from a stream of its own (so a rank is the same body whatever the others are): where it came from and what it is made of (`captured_composition`), and its orbit
   (`irregular_orbits.draw_orbit`) drawn again until it is one that stays (`irregular_zone.verdict`), up to `ORBIT_TRIES`.

Nothing about a planet is drawn from the global stream: every number comes from `draw.Stream` keys made of the star's mass and the planet's orbit and mass, so the same planet gets the same
irregulars. Stand-ins are IW1 to IW4 of `docs/design/research/gap-filling-plan.md` section 6.191.
"""

import dataclasses
import math
from dataclasses import dataclass, field

from planetgen.physics import body_status, captured_composition, compaction, irregular_orbits, irregular_population, irregular_zone, moon_composition
from planetgen.util import draw

RESERVOIR_BINARY_FRACTION = 0.3
"""float: The share of the planetesimals of the reservoir that are binaries, which decides the chance of a capture by exchange (`irregular_population.exchange_chance`). About a third of the cold classical Kuiper Belt objects are
binaries (recalled, not read). [Stand-in IW2, tier T3, tweak level L1.]"""
ENCOUNTER_REACH_AU = 3.0
"""float: How far from a planet's orbit a belt's nearest edge may lie and still give the planet bodies, AU, at a planet with Jupiter's Hill radius (the 3 AU encounter sphere of NDD14's `f_3AU`); it scales with the
planet's Hill radius (`irregular_population.REFERENCE_HILL_AU`). [Stand-in CD2, tier T2, tweak level L1.]"""
ORBIT_TRIES = 400
"""int: How many orbits are drawn for one stored body before it is given up as lost (the pool counts bodies that stayed, so a body is drawn again until its orbit stays). [Stand-in IW3, tier T6, tweak level L0.]"""
COMPOSITION_MASS_CAP_KG = compaction.TAPER_END_KG * 0.999
"""float: The heaviest mass the composition rows reach (`captured_composition.compose` stops at `compaction.TAPER_END_KG`); a heavier body (an exchange capture of a Triton) has the composition of this mass. [Stand-in IW4, tier T6, tweak level L2.]"""


MOON_KEYS = ("rank", "window", "mass_kg", "radius_km", "density_kg_m3", "taxon", "source_distance_au", "ice_loss", "rock", "metal", "water", "ice",
             "semi_major_axis_au", "eccentricity", "inclination_deg", "ascending_node_deg", "arg_periapsis_deg", "mean_anomaly_deg", "belt_order", "belt_rank")
"""tuple: The keys of an `IrregularMoon.to_dict`, which are also the keys of the stored row (`db.store`) bar the distances' units."""


@dataclass(frozen=True)
class IrregularMoon:
    """
    One stored irregular moon.

    Attributes:
        rank (int): 0 for the largest of the planet's stored bodies when it was built; the rank names the body for its planet and stays when a smaller body is lost.
        window (str): `gas`, `encounter` or `exchange` (`captured_composition.WINDOWS`).
        mass_kg, radius_km, density_kg_m3 (float): Its mass and, from the composition it was given, its radius and bulk density.
        taxon (str): The class letter of the body at its source (`physics/belt_member.py`).
        source_distance_au (float): Where in the disk it came from.
        ice_loss (float): The share of its water ice that capture removed.
        rock, metal, water, ice (float): The shares of its mass, after the loss.
        orbit (irregular_orbits.Orbit): Its orbit about the planet at the epoch of the star's age.
        belt_order (int or None): For a body caught from an inner belt, the belt's number among the belts of the planet's star, nearest first (`belts_of`); None for one from the outer disk or the planet's own zone.
        belt_rank (int or None): For a body that was a stored member of that belt, its rank there: it keeps the designation `<belt ID>.<rank hex>` and the belt keeps the row, marked as a moon. None for a body of
            the belt below its stored members.
    """
    rank: int
    window: str
    mass_kg: float
    radius_km: float
    density_kg_m3: float
    taxon: str
    source_distance_au: float
    ice_loss: float
    rock: float
    metal: float
    water: float
    ice: float
    orbit: object
    belt_order: object = None
    belt_rank: object = None

    @property
    def retrograde(self):
        return self.orbit.retrograde

    def to_dict(self):
        """The moon as a dict of `MOON_KEYS` (the orbit's six elements flattened in)."""
        o = self.orbit
        return {"rank": self.rank, "window": self.window, "mass_kg": self.mass_kg, "radius_km": self.radius_km, "density_kg_m3": self.density_kg_m3, "taxon": self.taxon,
                "source_distance_au": self.source_distance_au, "ice_loss": self.ice_loss, "rock": self.rock, "metal": self.metal, "water": self.water, "ice": self.ice,
                "semi_major_axis_au": o.semi_major_axis_au, "eccentricity": o.eccentricity, "inclination_deg": o.inclination_deg, "ascending_node_deg": o.ascending_node_deg,
                "arg_periapsis_deg": o.arg_periapsis_deg, "mean_anomaly_deg": o.mean_anomaly_deg, "belt_order": self.belt_order, "belt_rank": self.belt_rank}

    @classmethod
    def from_dict(cls, data):
        orbit = irregular_orbits.Orbit(data["semi_major_axis_au"], data["eccentricity"], data["inclination_deg"], data["ascending_node_deg"], data["arg_periapsis_deg"],
                                       data["mean_anomaly_deg"])
        return cls(data["rank"], data["window"], data["mass_kg"], data["radius_km"], data["density_kg_m3"], data["taxon"], data["source_distance_au"], data["ice_loss"],
                   data["rock"], data["metal"], data["water"], data["ice"], orbit, data.get("belt_order"), data.get("belt_rank"))


@dataclass(frozen=True)
class IrregularSet:
    """
    The irregular moons of a planet and what was not stored.

    Attributes:
        moons (tuple): The stored `IrregularMoon`s, largest first.
        remainder_kg (float): The expected mass of the bodies not stored (`irregular_population.Irregulars.remainder_kg`).
        above_floor (int): How many bodies are above the storage floor (`irregular_population.STORE_FLOOR_KG`).
        lost (int): Stored bodies lost: given up for want of an orbit that stays (`ORBIT_TRIES`), or that stopped staying when the planet moved (`keep_stable`).
        temporary (bool): Whether the planet's zone is wide enough for temporary captures (`irregular_population.holds_temporary_captures`): a flag, no object.
        belt_remainders (tuple): `(belt_order, kg)` for the part of `remainder_kg` that came from an inner belt (not stored: the belts' masses are, once `take_from_belts` has run).
    """
    moons: tuple = ()
    remainder_kg: float = 0.0
    above_floor: int = 0
    lost: int = 0
    temporary: bool = False
    belt_remainders: tuple = ()

    def to_dict(self):
        return {"moons": [moon.to_dict() for moon in self.moons], "remainder_kg": self.remainder_kg, "above_floor": self.above_floor, "lost": self.lost, "temporary": self.temporary}

    @classmethod
    def from_dict(cls, data):
        """The set a `to_dict` made; `None` or an empty dict is `NO_IRREGULARS`."""
        if not data:
            return NO_IRREGULARS
        return cls(tuple(IrregularMoon.from_dict(moon) for moon in data["moons"]), data["remainder_kg"], data["above_floor"], data["lost"], data["temporary"])


NO_IRREGULARS = IrregularSet()


def take_from_disks(disks, sets):
    """
    Takes the mass the planets caught out of the outer disks they caught it from (GEN.1085), so a belt's mass less what its planets hold is what the belt keeps: the stored bodies and the remainder of every
    `IrregularSet` in `sets` are shared among `disks` (`AsteroidBelt`s) by their primordial mass, and each disk's `mass_kg` falls by its share, never by more than its remainder (the mass its stored members do not
    hold), so the members it keeps are never short. The mass taken is added to the disk's `captured_mass_kg`, which is not stored: the planets' own sets say what they hold.

    Returns:
        float: The mass taken, kg.
    """
    caught = (sum(moon.mass_kg for held in sets for moon in held.moons if moon.belt_order is None)
              + sum(held.remainder_kg - sum(kg for _order, kg in held.belt_remainders) for held in sets))
    total = sum(disk.primordial_mass_earth for disk in disks)
    if not caught > 0.0 or not total > 0.0:
        return 0.0
    taken = 0.0
    for disk in disks:
        room = max(disk.mass_kg - sum(member["mass_kg"] for member in disk.members), 0.0)
        part = min(caught * disk.primordial_mass_earth / total, room)
        disk.mass_kg -= part
        disk.captured_mass_kg = disk.captured_mass_kg + part
        taken += part
    return taken


def take_from_belts(belts, sets):
    """
    Takes what the planets caught out of the inner belts they caught it from (GEN.1085, part 2): a belt's `mass_kg` falls by the mass of the bodies of it that planets hold and by its part of their remainders, never by more
    than the mass its stored members that are still in it do not hold. A stored member that a planet caught stays in the belt's list under its rank and designation, with the status `moon` (it is a moon now) and
    leaves the belt's mass. `belts` maps a belt's order (`belts_of`) to the belt. Returns the mass taken, kg.
    """
    taken = 0.0
    for order, belt in belts.items():
        moons = [moon for held in sets for moon in held.moons if moon.belt_order == order]
        remainder = sum(kg for held in sets for pair_order, kg in held.belt_remainders if pair_order == order)
        caught = sum(moon.mass_kg for moon in moons) + remainder
        if not caught > 0.0:
            continue
        for moon in moons:
            if moon.belt_rank is None:
                continue
            for member in belt.members:
                if member["rank"] == moon.belt_rank:
                    member["status"], member["soter_mu"] = body_status.MOON, None
        room = max(belt.mass_kg - sum(member["mass_kg"] for member in belt.members if member["status"] != body_status.MOON), 0.0)
        part = min(caught, room)
        belt.mass_kg -= part
        belt.captured_mass_kg = belt.captured_mass_kg + part
        taken += part
    return taken


def belts_of(bodies):
    """The belts among `bodies`, nearest first, as `{order: belt}`: a belt's order is its place among them, which a stored moon names it by."""
    return dict(enumerate(body for body in bodies if body.body_type == 'a'))


@dataclass(frozen=True)
class Source:
    """
    What the planets of a star can catch from.

    Attributes:
        reservoir_mass_earth (float): The primordial mass of the outer disks of the system, Earth masses (0 for none).
        inner_au, outer_au (float): The extent of those disks (0 and 0 for none).
        encounters (float): The giant-giant close encounters of the system (`capture_windows.Windows.encounters`).
        luminosity_w, aluminium_26_ratio (float): The star's, for the composition of what was caught.
        star_mass_kg, age_gy (float): The star's.
        belts (tuple): `(order, belt)` for each belt that is not an outer disk (`belts_of` gives the order): the inner belts a planet near enough can catch from (GEN.1085).
        taken (dict): `{order: set of ranks}` of the stored members planets have caught so far, so that a member is caught once; filled as the planets of the star are given their moons, in the order of the bodies.
    """
    reservoir_mass_earth: float
    inner_au: float
    outer_au: float
    encounters: float
    luminosity_w: float
    aluminium_26_ratio: float
    star_mass_kg: float
    age_gy: float
    belts: tuple = ()
    taken: dict = field(default_factory=dict, compare=False)


def source_of(bodies, star, windows, aluminium_26_ratio):
    """The `Source` of the planets among `bodies`: the outer-disk belts' primordial mass and extent (`belt.origin == "outer-disk"`), the other belts (the inner ones, each its own reservoir for a planet within reach),
    the windows of the system (`capture_windows.Windows`) and the star."""
    disks = [belt for belt in bodies if belt.body_type == 'a' and getattr(belt, "origin", None) == "outer-disk"]
    mass_earth = sum(belt.primordial_mass_earth for belt in disks)
    inner = min((belt.lower_limit for belt in disks), default=0.0)
    outer = max((belt.upper_limit for belt in disks), default=0.0)
    inner_belts = tuple((order, belt) for order, belt in belts_of(bodies).items() if getattr(belt, "origin", None) != "outer-disk" and (getattr(belt, "primordial_mass_earth", None) or 0.0) > 0.0)
    return Source(mass_earth, inner, outer, windows.encounters, star.luminosity, aluminium_26_ratio, star.mass, star.age, inner_belts)


def _compose(source, planet, window, mass_kg, stream):
    """The composition of a body of `window` from a stream of its own: its source distance, taxon, porosity and ice loss."""
    inner, outer = captured_composition.source_annulus(window, planet.distance, planet.mass, source.star_mass_kg, source.inner_au, source.outer_au)
    distance = captured_composition.draw_source_distance(inner, outer, stream.random())
    return captured_composition.compose(window, min(mass_kg, COMPOSITION_MASS_CAP_KG), distance, source.luminosity_w, source.aluminium_26_ratio,
                                        stream.random(), stream.random(), stream.random())


def _orbit(zone, planet, stream):
    """An orbit that stays in `zone`, drawn from `stream` (up to `ORBIT_TRIES`), or `None`."""
    with draw.bound(stream):
        for _ in range(ORBIT_TRIES):
            orbit = irregular_orbits.draw_orbit(zone, planet.distance)
            if irregular_zone.verdict(zone, orbit.semi_major_axis_au, orbit.eccentricity, orbit.inclination_deg) is None:
                return orbit
    return None


def for_planet(planet, source):
    """
    The `IrregularSet` of `planet`, a planet that orbits a star (not a moon, not a rogue), caught from `source`.

    The regular moons of the planet (`planet.moons`) set the inner edge of the zone. A planet with no outer disk to catch from and no gas window gets none.
    """
    if source.reservoir_mass_earth <= 0.0 and not source.belts:
        return NO_IRREGULARS
    outermost_moon_au = max((moon.distance for moon in planet.moons), default=0.0)
    zone = irregular_zone.zone(planet.distance, planet.mass, source.star_mass_kg, outermost_moon_au)
    if zone.empty:
        return NO_IRREGULARS
    temporary = irregular_population.holds_temporary_captures(zone.prograde_limit_au, planet.radius)
    gas_giant = (planet.envelope_mass_fraction or 0.0) >= moon_composition.GIANT_ENVELOPE_FRACTION
    pool = irregular_population.make_pool(zone.hill_radius_au, source.encounters, source.reservoir_mass_earth, gas_giant, getattr(planet, "moon_disk_lifetime_myr", 0.0) or 0.0)
    probability = irregular_population.intrinsic_probability(planet.mass, irregular_population.zone_radius_au(zone.hill_radius_au, planet.distance))
    age_yr = source.age_gy * 1.0e9
    key = f"irregulars:{source.star_mass_kg!r}:{planet.distance!r}:{planet.mass!r}"
    caught = pool.scale > 0.0 and (pool.encounter_fraction > 0.0 or pool.gas_fraction > 0.0)       # no window open, or no outer disk: no population to grind (its collision time is infinite)
    population = irregular_population.draw_irregulars(pool, probability, age_yr, draw.Stream(key)) if caught else irregular_population.Irregulars((), (), (), 0.0, 0)
    bodies = [(radius, mass, "gas" if gas else "encounter", None) for radius, mass, gas in zip(population.radii_km, population.masses_kg, population.gas_window)]
    chance = irregular_population.exchange_chance(RESERVOIR_BINARY_FRACTION, zone.hill_radius_au, gas_giant)
    exchange = draw.Stream(f"{key}:exchange")
    radius = irregular_population.draw_exchange(exchange.random(), exchange.random() * (1.0 - 1e-12), chance if source.reservoir_mass_earth > 0.0 else 0.0)       # a binary planetesimal comes from the outer disk
    if radius is not None:
        bodies.append((radius, irregular_population.mass_kg(radius), "exchange", None))
        bodies.sort(key=lambda body: -body[1])
    remainder_kg, above_floor, belt_remainders = population.remainder_kg, population.above_floor, []
    reach_au = ENCOUNTER_REACH_AU * zone.hill_radius_au / irregular_population.REFERENCE_HILL_AU
    caught_from_belts = False
    for order, belt in source.belts:
        if _gap_au(planet.distance, belt) > reach_au:
            continue
        belt_pool = irregular_population.make_pool(zone.hill_radius_au, source.encounters, belt.primordial_mass_earth)      # the encounter window only: no gas disk is left
        if not belt_pool.encounter_fraction > 0.0:
            continue
        caught = irregular_population.draw_irregulars(belt_pool, probability, age_yr, draw.Stream(f"{key}:belt{order}"))
        bodies += [(radius, mass, "encounter", order) for radius, mass in zip(caught.radii_km, caught.masses_kg)]
        remainder_kg, above_floor = remainder_kg + caught.remainder_kg, above_floor + caught.above_floor
        belt_remainders.append((order, caught.remainder_kg))
        caught_from_belts = True
    if caught_from_belts:
        bodies.sort(key=lambda body: -body[1])
    moons, lost, claimed = [], 0, {}
    for rank, (_radius, mass, window, order) in enumerate(bodies):
        stream = draw.Stream(f"{key}:{rank}" if order is None else f"{key}:belt{order}:{rank}")
        orbit = _orbit(zone, planet, stream)
        if orbit is None:
            lost += 1
            continue
        if order is not None:
            moons.append(_from_belt(source, claimed, dict(source.belts)[order], order, len(moons), window, mass, orbit, stream))
            continue
        captured = _compose(source, planet, window, mass, stream)
        makeup = captured.makeup
        moons.append(IrregularMoon(len(moons), window, mass, captured.radius_km, captured.density_kg_m3, captured.taxon, captured.source_distance_au, captured.ice_loss,
                                   makeup.rock, makeup.metal, makeup.water, makeup.ice, orbit))
    if not moons and not remainder_kg and not temporary:
        return NO_IRREGULARS
    return IrregularSet(tuple(moons), remainder_kg, above_floor, lost, temporary, tuple(belt_remainders))


def record(source, held):
    """Notes the stored members of belts that `held` (an `IrregularSet`) caught in `source.taken`, so that a later planet of the star cannot catch them again."""
    for moon in held.moons:
        if moon.belt_rank is not None:
            source.taken.setdefault(moon.belt_order, set()).add(moon.belt_rank)


def _gap_au(distance_au, belt):
    """How far the planet's orbit lies from the nearest edge of `belt`, AU (0 inside it)."""
    return max(belt.lower_limit - distance_au, distance_au - belt.upper_limit, 0.0)


def _from_belt(source, claimed, belt, order, rank, window, mass, orbit, stream):
    """
    The `IrregularMoon` of a body caught from an inner belt (GEN.1085, part 2). A body at least as heavy as the lightest stored member of the belt is one of its stored members (`belt.members`, heaviest first): the one
    nearest to it in mass of those no planet of the star has caught yet (`source.taken`, which `record` fills from the sets of the planets before it), whose mass, size, makeup and place it takes as they stand (the belt shows them), and whose rank it keeps. A lighter body is of the belt
    below its stored members: it is composed from the belt's own stretch like a body of the outer disk, and has no rank.
    """
    members = belt.members
    if members and mass >= min(member["mass_kg"] for member in members):
        caught = claimed.setdefault(order, set())
        free = [member for member in members if member["rank"] not in caught and member["rank"] not in source.taken.get(order, ())]
        if free:
            member = min(free, key=lambda one: abs(math.log(one["mass_kg"] / mass)))
            caught.add(member["rank"])
            return IrregularMoon(rank, window, member["mass_kg"], member["radius_km"], member["density_g_cm3"] * 1000.0, member["taxon"], member["distance_au"], 0.0,
                                 member["rock_fraction"], member["metal_fraction"], member["water_fraction"], member["ice_fraction"], orbit, order, member["rank"])
    distance = captured_composition.draw_source_distance(belt.lower_limit, belt.upper_limit, stream.random())
    captured = captured_composition.compose(window, min(mass, COMPOSITION_MASS_CAP_KG), distance, source.luminosity_w, source.aluminium_26_ratio, stream.random(), stream.random(), stream.random())
    makeup = captured.makeup
    return IrregularMoon(rank, window, mass, captured.radius_km, captured.density_kg_m3, captured.taxon, captured.source_distance_au, captured.ice_loss,
                         makeup.rock, makeup.metal, makeup.water, makeup.ice, orbit, order, None)



def keep_stable(planet):
    """
    Drops the irregular moons of `planet` that no longer stay in its zone (the planet's orbit, mass or regular moons changed: `physics.planets.drop_unstable_moons`), counting them as lost. Returns how many were dropped.
    """
    held = getattr(planet, "irregulars", NO_IRREGULARS)
    if not held.moons:
        return 0
    zone = irregular_zone.zone(planet.distance, planet.mass, planet.star.mass, max((moon.distance for moon in planet.moons), default=0.0))
    kept = tuple(moon for moon in held.moons
                 if irregular_zone.verdict(zone, moon.orbit.semi_major_axis_au, moon.orbit.eccentricity, moon.orbit.inclination_deg) is None)
    dropped = len(held.moons) - len(kept)
    if dropped:
        planet.irregulars = dataclasses.replace(held, moons=kept, lost=held.lost + dropped)
    return dropped
