# planetgen/generation/belt_birth.py

"""
The belts of a system, from its finished planets (GEN.1068, GEN.1156)
=====================================================================

A belt is a gap or an unfinished tail in the chain of planets the disk built (`physics/belt_origin.py`), not a slot the chain rolled for. This module reads the chain a star's planets
make once they are placed and gives back the `AsteroidBelt` of every stretch the disk left unfinished, each with the cause that made it, what the giants took from it and what is left
after the star's age of collisions:

1. The planets of the chain become `belt_origin.Body` (their distance, mass, whether they are giants, and whether a giant moved: it ends far from where it formed).
2. `belt_origin.build` names the stretches (gaps, the tail, the snow-line bump, the outer disk and, for a rare system, the rubble of a collision), traces each through the giants' clearing
   and a migrating giant's sweep, and gives the mass it formed with and the mass it has after them.
3. An outer disk swept by the outermost giant is split into the populations of `physics/scatter.py` (cold, resonant, scattered, hot), recorded on the belt with the sentence of what the
   giant did. The mass is not changed again: `belt_origin` already took the giants' share, and a second subtraction would count the same sweep twice.
   The mass the disk lost (`belt_origin`: cleared by the giant and scattered by its move) is then handed on (`hand_on`, GEN.1308, `physics/swept_hand_on.py`): the thrown-out share to the rogue builders
   (`disk_ejection`), the implanted share to the nearest inner belt, the captured share to the planets' irregular moons.
4. The stretch is cut to the orbits the star can hold (`orbit_floor_au`, `orbit_ceiling_au`), and the belt's mass after `age` of collisions (`belt_mass`, [R4785]) must be
   `MIN_PRESENT_MASS_EARTH` or more: a belt of dust the collisions have ground away is not a belt a survey would find.

The draws come from streams of the system's own, so the run's stream moves only for the belts that are kept (their composition).

Stand-in
--------
`MIN_PRESENT_MASS_EARTH`: the smallest mass a belt can have and be kept. It is the one number of the rule, tuned so that the share of main-sequence systems with an inner belt (the asteroid
belt, out to `INNER_BELT_SNOW_LINE_RATIO` snow lines) lands near 9 percent, inside the 4 to 14 percent that the surveys of GEN.1148 leave (Boss, GEN.1156: the inner belt counts there, the outer
and Kuiper belts are separate kinds). Without a floor 12 percent have one and the narrow belts near the star are ground to a few times 1e-10 Earth masses, a dust disk and not a belt. The
mass is that of warm dust a survey could see, a stand-in. [Stand-in BG14, tier Marginal educated guess from needed output, tweak L3.]
"""

from planetgen import tuning
from planetgen.generation.belt import AsteroidBelt, present_state
from planetgen.generation import remnant_orbits
from planetgen.physics import belt_origin, bulk_draw, scatter, swept_hand_on
from planetgen.physics import unit_registry as units
from planetgen.util import draw

_kg_to_msun = units.converter("kg", "Msun")
_kg_to_mearth = units.converter("kg", "Mearth")
_km_to_rearth = units.converter("km", "Rearth")

MIN_PRESENT_MASS_EARTH = 1.0e-9
"""float: The smallest mass, Earth masses, a belt can have after the star's age of collisions and be kept. [Stand-in BG14, tier Marginal educated guess from needed output, tweak L3.]"""

INNER_BELT_SNOW_LINE_RATIO = 1.2
"""float: A belt that ends within this many snow lines of the star is an inner belt (the asteroid belt kind); one that reaches further is an outer or Kuiper belt. [Stand-in BG15, tier Marginal educated guess from physics, tweak L0.]"""

MIGRATED_RATIO = 1.25
"""float: A giant ended this many times nearer (or further) than where it formed is one that moved (`belt_origin.Body.migrated`); the placement of a slot the planets were nudged to is less. [Stand-in BG16, tier Marginal educated guess from needed output, tweak L1.]"""


def is_inner(belt, snow_line_au):
    """Whether a belt is an inner one: it ends within `INNER_BELT_SNOW_LINE_RATIO` snow lines of the star."""
    return belt.upper_limit <= INNER_BELT_SNOW_LINE_RATIO * snow_line_au


def chain_of(planets):
    """The `belt_origin.Body` list of the planets in `planets` (a belt in the list is skipped), nearest the star first."""
    chain = []
    for planet in planets:
        if planet.body_type == 'a':
            continue
        formed = planet.formation_distance_au or planet.distance
        giant = planet.body_type == 'g'
        moved = giant and max(formed / planet.distance, planet.distance / formed) >= MIGRATED_RATIO
        chain.append(belt_origin.Body(planet.distance, planet.mass, planet.eccentricity or 0.0, giant, moved, formed, planet.radius))
    return sorted(chain, key=lambda body: body.distance_au)


def _stream(star, chain, label):
    return draw.Stream(f"belt-birth:{label}:{star.mass!r}:" + ":".join(f"{body.distance_au!r}:{body.mass_kg!r}" for body in chain))


def _debris(star, chain):
    """The `Debris` of a collision of large bodies, if the system's age brings one (a rare system), else an empty list."""
    if not chain:
        return []
    stream = _stream(star, chain, "debris")
    body = chain[min(int(stream.random() * len(chain)), len(chain) - 1)]
    found = belt_origin.collisional_debris(body.distance_au, body.eccentricity, star.age * 1.0e9, stream.random())
    return [found] if found is not None else []


def _fates(star_mass_kg, planets, origin, disk_edge_au, instability=0.0):
    """`(shares, sentence, giant_au, giant)` of an outer disk after the outermost giant (`scatter.scatter_shares`, the giant's distance and the giant itself), or `(None, None, None, None)` where the belt is not an outer disk or no giant has been."""
    giants = sorted((planet for planet in planets if planet.body_type == 'g'), key=lambda planet: planet.distance)
    if origin.origin != "outer-disk" or not giants:
        return None, None, None, None
    giant = giants[-1]
    moved_out = max(giant.distance - (giant.formation_distance_au or giant.distance), 0.0)
    shares = scatter.scatter_shares(_kg_to_msun(star_mass_kg), _kg_to_mearth(giant.mass), giant.distance, _km_to_rearth(giant.radius), moved_out, origin.inner_au, disk_edge_au, instability)
    return shares, scatter.history_line(shares, giant.distance, moved_out), giant.distance, giant


def build(system_config, star, planets, orbit_floor_au, orbit_ceiling_au, disk_edge_au, instability=0.0, disk_scatter=1.0):
    """
    The `AsteroidBelt` list of a system: the stretches the planets in `planets` leave unfinished, each as `belt_origin` traced it, cut to the orbits between `orbit_floor_au` and
    `orbit_ceiling_au` and kept when `MIN_PRESENT_MASS_EARTH` is left after the star's age. Nearest the star first. A system asked to hold a belt (`ASTEROID_BELT` True) that has none
    above the floor keeps the heaviest one it has, so the request is met wherever any stretch is left.

    Args:
        system_config (SystemConfig): The system's shared config.
        star (Star): The star the planets orbit.
        planets (list): The `Planet` list (a belt in it is skipped, and its stretch is not counted twice by the caller).
        orbit_floor_au, orbit_ceiling_au (float): The innermost and outermost orbits the star holds.
        disk_edge_au (float): The disk's outer edge, AU (infinite when the star's light no longer describes its disk: then `tuning.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER` snow lines).
        instability (float): 0 to 1, how much the system's giants scattered each other (`capture_windows.instability_level`), which excites a share of the cold outer disk (`scatter.scatter_shares`).
        disk_scatter (float): The factor the planets' disk was scattered by about the median of its star's mass (`planet_architecture.Architecture.disk_scatter`): the belts are made of the same disk (GEN.1261).
    """
    chain = chain_of(planets)
    star_mass_kg = remnant_orbits.formation_mass_kg(star)
    snow_line_au = star.snow_line_au
    edge_au = disk_edge_au if disk_edge_au < float("inf") else tuning.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER * snow_line_au
    origins = belt_origin.build(chain, star_mass_kg, snow_line_au, snow_line_au, tuning.PLANET_OUTER_GROWTH_AU, edge_au, u=_stream(star, chain, "survival").random(),
                                debris=_debris(star, chain), metal_ratio=bulk_draw.metal_ratio_to_solar(star), scatter=disk_scatter, instability=instability)
    found = []
    for origin in origins:
        lower, upper = max(origin.inner_au, orbit_floor_au), min(origin.outer_au, orbit_ceiling_au)
        if upper < lower * belt_origin.MIN_BELT_RATIO or not origin.mass_earth > 0.0:       # a stretch the giants cleared of everything is no belt
            continue
        _, mass_kg, _ = present_state(star, lower, lower, upper, origin.mass_earth)
        found.append((_kg_to_mearth(mass_kg), origin, lower, upper))
    kept = [item for item in found if item[0] >= MIN_PRESENT_MASS_EARTH]
    if not kept and found and system_config.ASTEROID_BELT is True:
        kept = [max(found, key=lambda item: item[0])]
    belts = []
    for _mass, origin, lower, upper in kept:
        shares, sentence, giant_au, giant = _fates(star_mass_kg, planets, origin, edge_au, instability)
        belt = AsteroidBelt(system_config, star, lower, lower, upper, origin=origin, scatter_fates=shares, scatter_giant_au=giant_au)
        if giant is not None:
            belt.scatter_giant_kg, belt.scatter_giant_radius_km = giant.mass, giant.radius
        if sentence:
            belt.origin_sentence = f"{belt.origin_sentence} {sentence[0].upper()}{sentence[1:]}."
        belts.append(belt)
    return belts


IMPLANT_SENTENCE = "A giant implanted {mass} of bodies from the outer disk into it."
"""str: What a belt that received an implant says of it (GEN.1308); `{mass}` is the mass in Earth masses, written by the caller."""


def hand_on(bodies, star):
    """
    Hands on the mass the outer disks of a star lost (GEN.1308, `physics/swept_hand_on.py`): for each outer-disk belt among `bodies` with a scatter step, its gone mass (what was cleared from it and
    scattered, `belt_origin` and the move of the giant) is shared among `ejected`, `implanted`, `captured` and `other` by the scatter step's rules (`swept_hand_on.hand_on`, the giant's Safronov number), and recorded in the disk's `swept_kg`. The implanted
    part goes to the inner belt (a belt that is not an outer disk) nearest to the disk's inner edge, which grows by it (`AsteroidBelt.receive`); a disk with no inner belt to receive it keeps the part in `swept_kg`
    and it is lost to the star's inner system. The disk's own mass is not changed: it is the mass left after the gone mass, which `belt_origin` already took.
    """
    belts = [body for body in bodies if getattr(body, "body_type", None) == 'a' and getattr(body, "origin", None) is not None]
    disks = [belt for belt in belts if belt.origin == "outer-disk" and belt.scatter_giant_kg is not None and belt.swept_kg is None]
    inner = [belt for belt in belts if belt.origin != "outer-disk"]
    for disk in disks:
        gone = units.convert(disk.cleared_mass_earth + disk.scattered_mass_earth, "Mearth", "kg")
        theta = scatter.safronov_number(units.convert(disk.scatter_giant_kg, "kg", "Mearth"), units.convert(disk.scatter_giant_radius_km, "km", "Rearth"), disk.scatter_giant_au, units.convert(star.mass, "kg", "Msun"))
        disk.swept_kg = swept_hand_on.hand_on(gone, theta)
        if disk.swept_kg is None or not disk.swept_kg["implanted"] > 0.0 or not inner:
            continue
        receiver = min(inner, key=lambda belt: abs(belt.upper_limit - disk.lower_limit))
        receiver.receive(disk.swept_kg["implanted"], IMPLANT_SENTENCE.format(mass=f"{units.convert(disk.swept_kg['implanted'], 'kg', 'Mearth'):.3g} Earth masses"))
