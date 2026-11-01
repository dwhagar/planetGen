# planetgen/physics/belt_origin.py

"""
Where a belt came from (GEN.1088, docs/design/systems/belt-origins.md sections 3 to 5)
=====================================================================================

Every belt and outer disc gets the cause that made it, a record of how its solids were formed and what took mass from them, and one sentence
for its page ("The leftover of a planet that never grew here ..."). The generator made the belt, so it knows. A pure law over the planet chain
the architecture (`planet_architecture`) built and the disc's solids; it reads no database, draws no random number (the two draws it needs come
in as `u` arguments) and moves nothing. GEN.1068 gives it the final chain in place of the 10 percent slot roll.

**Formation causes** (what made the belt; the stored `origin`):

- `gap`: two neighbours of the chain stand more than `GAP_SPACING_K` mutual Hill radii apart (the thinned chains cluster at 12 to 20), and the
  space between them, less five Hill radii of each, is wide enough: the solids there did not grow into a planet.
- `tail`: the chain stops short of `growth_limit_au`, where the disc has built less (`PLANET_OUTER_GROWTH_AU`), and the space past its last
  rocky planet is a belt (a last giant ends the chain, the outer disc takes over).
- `snow-line-bump`: solids pile up at the snow line as the system formed it (the formation snow line is stored); where that stretch is free of
  planets and is not already inside a gap or tail (which then take it as a second cause) it is its own belt.
- `outer-disc`: the disc past the last giant (or past the tail), out to `disc_edge_au`.
- `destroyed-planet`: the remains of a planet the system's history destroyed (the caller names it), `DESTROYED_RETAINED_SHARE` of its mass.
- `collisional-debris`: the rubble of a collision between large bodies, a rare event at `COLLISIONAL_DEBRIS_PER_GYR` of mass
  `COLLISIONAL_DEBRIS_MASS_EARTH`; `collisional_debris` turns the caller's draw into the event or none.

**Shaping causes** (what took mass afterwards; the stored `history`), applied in this order to each belt, solids weighted by the disc's own
surface density (`belt_mass.primordial_mass_earth`):

- `giant-stirring`: the chaotic zone of each giant, half-width `1.3 a mu^(2/7)` each side [R4727], and for a belt inside the giant the
  gaps at the resonances 4:1, 3:1, 5:2, 7:3 and 2:1 (`a_g (q/p)^(2/3)`), half-width `RESONANCE_HALF_WIDTH` of the resonance's distance for a
  Jupiter-mass giant and as the square root of the giant's mass ratio. What is in them is gone: `cleared_earth`. A belt within
  `STIR_ZONE_WIDTHS` chaotic half-widths of a giant is stirred whether or not mass is cleared, and says so.
- `instability-depletion`: a giant that moved (the architecture's hot Jupiter, `migrated`) through the belt's stretch swept that part of it: a share
  `INSTABILITY_SURVIVAL_RANGE` of the swept stretch survives, log-uniform in the draw `u`; the rest is `scattered_earth` (the scatter step, GEN.1089,
  decides where it went).

Mass is conserved: `primordial_mass_earth = mass_earth + cleared_earth + scattered_earth` (`conservation_error`). The mass later lost to
collisional grinding is not here: that is `belt_mass.belt_mass_earth`, which needs the age and the largest body.

Where the design says the stored origin is "the dominant" cause, this module stores the formation cause as `origin` (the one that made the
belt) and every shaping cause with the share of the mass it took in `history`, so the page can say both.

Use:

    chain = [Body(distance_au, mass_kg, eccentricity, giant, migrated, formation_au), ...]
    belts = build(chain, star_mass_kg, snow_line_au, formation_snow_line_au, growth_limit_au, disc_edge_au, u=0.5)
    for belt in belts: belt.origin, belt.mass_earth, belt.sentence
"""

import collections
import json
import math
import os

from planetgen import tuning
from planetgen.physics import belt_mass, clearance_ledger, elements, planet_architecture, unit_registry as units
from planetgen.physics.elements import DataError

_kg_to_mearth = units.converter("kg", "Mearth")

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "belt_origins.json")
"""str: The data file: the causes and the sentence of each."""

GAP_SPACING_K = 40.0
"""float: Two neighbours of the chain this many mutual Hill radii apart leave a gap (BG1)."""

MIN_BELT_RATIO = 1.1
"""float: The least outer-to-inner distance ratio of a belt (BG2)."""

MIN_BELT_MASS_EARTH = 1.0e-6
"""float: The least primordial mass of a belt, about 6e18 kg, a body of the individual tier of the ledger (BG3)."""

SNOW_BUMP_WIDTH = 0.15
"""float: The snow-line bump spans the formation snow line divided and multiplied by `1 + SNOW_BUMP_WIDTH` (BG4)."""

CHAOTIC_ZONE_COEFFICIENT = 1.3
"""float: The half-width of a planet's chaotic zone is this times `a mu^(2/7)` ([R4727], READ in [R3175, R3226], row BR2) (BG5)."""

RESONANCES = ((4, 1), (3, 1), (5, 2), (7, 3), (2, 1))
"""tuple: The period ratios `p:q` (giant : body) of the gaps a giant cuts in a belt inside its orbit, `a_g (q/p)^(2/3)` (BG6)."""

RESONANCE_HALF_WIDTH = 0.015
"""float: The half-width of a resonance gap as a share of its distance, for a giant of `JUPITER_MASS_RATIO` (BG7)."""

JUPITER_MASS_RATIO = 9.546e-4
"""float: Jupiter over the Sun, the giant `RESONANCE_HALF_WIDTH` is for."""

STIR_ZONE_WIDTHS = 3.0
"""float: A belt within this many chaotic half-widths of a giant, outside its chaotic zone, is stirred by it (BG8)."""

INSTABILITY_SURVIVAL_RANGE = (0.002, 0.01)
"""tuple: The share of a belt a giant that migrated through it leaves, the least and the most: a Grand Tack type event leaves a few thousandths to a
hundredth of the mass ([R4603], registry row BR4); a Nice type instability alone leaves 0.1 to 0.5 (BR5), which waits for the
architecture to record instabilities (GEN.1082) (BG9)."""

DESTROYED_RETAINED_SHARE = 0.5
"""float: The share of a destroyed planet's mass that stays in orbit as a belt (BG10)."""

DEBRIS_MIN_WIDTH = 0.05
"""float: The least half-width of the belt of debris round its distance, as a share of the distance (BG11)."""

COLLISIONAL_DEBRIS_PER_GYR = 0.003
COLLISIONAL_DEBRIS_MASS_EARTH = 1.0e-4
"""float: The chance per Gyr of a standing belt of rubble from a collision of large bodies, and its mass (BG12 and BG13)."""

Body = collections.namedtuple("Body", ("distance_au", "mass_kg", "eccentricity", "giant", "migrated", "formation_au"), defaults=(0.0, False, False, None))
"""One body of the planet chain: where it ends, its mass, its eccentricity, whether it is a giant, whether it moved in (the architecture's hot
Jupiter) and where it formed."""

Debris = collections.namedtuple("Debris", ("cause", "mass_earth", "distance_au", "eccentricity"), defaults=(0.0,))
"""The event behind a belt of rubble: `cause` is 'destroyed-planet' or 'collisional-debris', `mass_earth` the mass of the belt, `distance_au` the
distance of the event."""

Event = collections.namedtuple("Event", ("cause", "mass_earth", "giant_au", "text"))
"""One shaping cause of a belt: what it took (Earth masses), the giant that did it and the sentence."""

BeltOrigin = collections.namedtuple("BeltOrigin", (
    "origin", "co_causes", "inner_au", "outer_au", "formation_snow_line_au", "primordial_mass_earth", "mass_earth", "cleared_earth",
    "scattered_earth", "stirred", "history", "sentence"))
"""A belt and where it came from: the formation cause, the others that coincide with it, the stretch, the snow line at formation, the mass it
formed with, what is left, what the giants cleared, what an instability scattered, whether a giant stirs it, the shaping causes in order and
the sentence for its page."""


class BeltOriginError(DataError):
    """A cause file that is refused, or an input the module cannot read."""


def _load(path=DATA_FILE):
    with open(path, encoding="utf-8") as handle:
        document = json.load(handle)
    causes = document.get("causes")
    if not isinstance(causes, dict) or not causes:
        raise BeltOriginError("the file needs 'causes'")
    for code, cause in causes.items():
        if cause.get("kind") not in ("formation", "shaping") or not cause.get("label") or not cause.get("text") or not cause.get("fields"):
            raise BeltOriginError(f"cause {code!r}: needs a kind (formation or shaping), a label, its fields and a text")
        try:
            cause["text"].format(**{name: 1.0 for name in cause["fields"]})
            for extra in ("co", "alt"):
                if f"{extra}_text" in cause:
                    cause[f"{extra}_text"].format(**{name: 1.0 for name in cause[f"{extra}_fields"]})
        except (KeyError, IndexError, ValueError, TypeError) as error:
            raise BeltOriginError(f"cause {code!r}: the text does not match its fields ({error})") from error
    return causes


CAUSES = _load()
"""dict: Each cause's kind, label, fields and sentence."""

FORMATION_CAUSES = tuple(code for code, cause in CAUSES.items() if cause["kind"] == "formation")
SHAPING_CAUSES = tuple(code for code, cause in CAUSES.items() if cause["kind"] == "shaping")


# --- intervals --------------------------------------------------------------------------------------------------------------------------

def _clip(intervals, lower, upper):
    return [(max(a, lower), min(b, upper)) for a, b in intervals if min(b, upper) > max(a, lower)]


def _subtract(spans, cut):
    """The pieces of `spans` (a list of `(low, high)`) outside `cut` (one `(low, high)`)."""
    low, high = cut
    out = []
    for a, b in spans:
        if b <= low or a >= high:
            out.append((a, b))
            continue
        if a < low:
            out.append((a, low))
        if b > high:
            out.append((high, b))
    return out


# --- the chain --------------------------------------------------------------------------------------------------------------------------

def _hill_share(body, star_mass_kg):
    """`(m / 3M)^(1/3)` of a body: its Hill radius over its distance."""
    return (body.mass_kg / (3.0 * (star_mass_kg + body.mass_kg))) ** (1.0 / 3.0)


def keep_out(body, star_mass_kg):
    """The stretch of distance a belt may not take beside a body: its orbit and five Hill radii (`clearance_ledger.BELT_HILL_RADII`)."""
    share = _hill_share(body, star_mass_kg) * clearance_ledger.BELT_HILL_RADII
    return body.distance_au * (1.0 - body.eccentricity - share), body.distance_au * (1.0 + body.eccentricity + share)


def chaotic_half_width_au(body, star_mass_kg):
    """The half-width of the chaotic zone of a planet, each side of its orbit [R4727], in AU."""
    return CHAOTIC_ZONE_COEFFICIENT * body.distance_au * (body.mass_kg / star_mass_kg) ** (2.0 / 7.0)


def resonance_distances_au(body):
    """The distances inside a giant's orbit where the resonances of `RESONANCES` stand, `a (q/p)^(2/3)`, as `((p, q), au)`."""
    return [((p, q), body.distance_au * (q / p) ** (2.0 / 3.0)) for p, q in RESONANCES]


def _primordial(lower, upper, star_mass_kg, snow_line_au, metal_ratio, scatter):
    return belt_mass.primordial_mass_earth(lower, upper, star_mass_kg, snow_line_au, metal_ratio, scatter) if upper > lower else 0.0


# --- the formation causes ---------------------------------------------------------------------------------------------------------------

def _candidate(cause, lower, upper, **fields):
    return {"cause": cause, "inner_au": lower, "outer_au": upper, "co": [], **fields}


def _fits(lower, upper):
    return lower > 0.0 and upper >= lower * MIN_BELT_RATIO


def gaps(chain, star_mass_kg):
    """
    The belts of the gaps of the chain: for each pair of neighbours more than `GAP_SPACING_K` mutual Hill radii apart, the stretch between the
    keep-out of the inner and the keep-out of the outer, less `tuning.MIN_ASTEROID_BELT_SEPARATION` at its outer edge. Each is
    `{'cause': 'gap', 'inner_au', 'outer_au', 'inner_neighbour_au', 'outer_neighbour_au'}`.
    """
    found = []
    for near, far in zip(chain, chain[1:]):
        if planet_architecture.spacing_k(near.distance_au, near.mass_kg, far.distance_au, far.mass_kg, star_mass_kg) < GAP_SPACING_K:
            continue
        lower, upper = keep_out(near, star_mass_kg)[1], keep_out(far, star_mass_kg)[0] - tuning.MIN_ASTEROID_BELT_SEPARATION
        if _fits(lower, upper):
            found.append(_candidate("gap", lower, upper, inner_neighbour_au=near.distance_au, outer_neighbour_au=far.distance_au))
    return found


def tail(chain, star_mass_kg, growth_limit_au):
    """The unfinished tail: past the last rocky planet to `growth_limit_au`, only when no giant ends the chain (the outer disc takes that one)."""
    if not chain or chain[-1].giant:
        return []
    last = chain[-1]
    lower, upper = keep_out(last, star_mass_kg)[1], growth_limit_au
    return [_candidate("tail", lower, upper, inner_neighbour_au=last.distance_au)] if _fits(lower, upper) else []


def outer_disc(chain, star_mass_kg, growth_limit_au, disc_edge_au):
    """
    The outer disc: from the outer chaotic edge of the last giant and the keep-out of the last body, and from `growth_limit_au` too when the chain
    ends in a rocky planet (the tail has the space before it), out to `disc_edge_au`.
    """
    giants = [body for body in chain if body.giant]
    lower = growth_limit_au if not chain or not chain[-1].giant else 0.0
    if chain:
        lower = max(lower, keep_out(chain[-1], star_mass_kg)[1])
    if giants:
        lower = max(lower, giants[-1].distance_au + chaotic_half_width_au(giants[-1], star_mass_kg))
    neighbour = giants[-1].distance_au if giants else (chain[-1].distance_au if chain else growth_limit_au)
    return [_candidate("outer-disc", lower, disc_edge_au, inner_neighbour_au=neighbour)] if _fits(lower, disc_edge_au) else []


def snow_line_bump(chain, star_mass_kg, formation_snow_line_au, taken):
    """
    The snow-line bump: the stretch of `SNOW_BUMP_WIDTH` either side of the formation snow line that no planet's keep-out covers. If a belt of
    `taken` (the other formation candidates) already holds the snow line it gets the cause as a second one instead and nothing is returned.
    """
    for belt in taken:
        if belt["inner_au"] <= formation_snow_line_au <= belt["outer_au"]:
            belt["co"].append("snow-line-bump")
            return []
    spans = [(formation_snow_line_au / (1.0 + SNOW_BUMP_WIDTH), formation_snow_line_au * (1.0 + SNOW_BUMP_WIDTH))]
    for body in chain:
        spans = _subtract(spans, keep_out(body, star_mass_kg))
    spans = [span for span in spans if _fits(*span)]
    if not spans:
        return []
    lower, upper = max(spans, key=lambda span: span[1] / span[0])
    return [_candidate("snow-line-bump", lower, upper, snow_line_au=formation_snow_line_au)]


def destroyed_planet(mass_kg, distance_au, eccentricity=0.0):
    """The `Debris` of a planet of `mass_kg` the system's history destroyed at `distance_au`: `DESTROYED_RETAINED_SHARE` of its mass stays."""
    return Debris("destroyed-planet", _kg_to_mearth(DESTROYED_RETAINED_SHARE * mass_kg), distance_au, eccentricity)


def collisional_debris_chance(age_yr):
    """The chance that a system of `age_yr` has a standing belt of rubble from a collision of large bodies: `1 - exp(-rate t)`."""
    return 1.0 - math.exp(-COLLISIONAL_DEBRIS_PER_GYR * age_yr / 1.0e9)


def collisional_debris(distance_au, eccentricity, age_yr, u):
    """The `Debris` of a collision at `distance_au` if the draw `u` (uniform, 0 to 1) is under the chance for `age_yr`, else `None`."""
    if u < collisional_debris_chance(age_yr):
        return Debris("collisional-debris", COLLISIONAL_DEBRIS_MASS_EARTH, distance_au, eccentricity)
    return None


def _debris_candidate(debris):
    if debris.cause not in ("destroyed-planet", "collisional-debris"):
        raise BeltOriginError(f"debris cause {debris.cause!r} is not destroyed-planet or collisional-debris")
    width = max(debris.eccentricity, DEBRIS_MIN_WIDTH)
    return _candidate(debris.cause, debris.distance_au * (1.0 - width), debris.distance_au * (1.0 + width), event_au=debris.distance_au,
                      mass_earth=debris.mass_earth)


# --- the shaping causes -----------------------------------------------------------------------------------------------------------------

def cleared_spans(body, star_mass_kg, lower_au, upper_au):
    """The stretches of `[lower_au, upper_au]` a giant clears: its chaotic zone and, for the part inside its orbit, the resonance gaps."""
    half = chaotic_half_width_au(body, star_mass_kg)
    spans = [(body.distance_au - half, body.distance_au + half)]
    width = RESONANCE_HALF_WIDTH * math.sqrt(body.mass_kg / (star_mass_kg * JUPITER_MASS_RATIO))
    spans += [(distance * (1.0 - width), distance * (1.0 + width)) for _, distance in resonance_distances_au(body)]
    return _clip(sorted(spans), lower_au, upper_au)


def is_stirred(body, star_mass_kg, lower_au, upper_au):
    """Whether a giant stirs the stretch: it reaches within `STIR_ZONE_WIDTHS` chaotic half-widths of the giant's orbit."""
    reach = STIR_ZONE_WIDTHS * chaotic_half_width_au(body, star_mass_kg)
    return upper_au >= body.distance_au - reach and lower_au <= body.distance_au + reach


def moved_through(body, lower_au, upper_au, star_mass_kg):
    """Whether a giant that moved (`migrated`) passed through the stretch: between where it formed and where it ended, a chaotic half-width each side."""
    if not (body.giant and body.migrated and body.formation_au):
        return False
    half = chaotic_half_width_au(body, star_mass_kg)
    path = (min(body.formation_au, body.distance_au) - half, max(body.formation_au, body.distance_au) + half)
    return upper_au >= path[0] and lower_au <= path[1]


def instability_survival(u):
    """The share a moving giant leaves, log-uniform over `INSTABILITY_SURVIVAL_RANGE` in the draw `u` (0 to 1)."""
    low, high = INSTABILITY_SURVIVAL_RANGE
    return math.exp(math.log(low) + min(max(u, 0.0), 1.0) * (math.log(high) - math.log(low)))


# --- the belt ---------------------------------------------------------------------------------------------------------------------------

def _text(cause, **fields):
    return CAUSES[cause]["text"].format(**{name: fields[name] for name in CAUSES[cause]["fields"]})


def _trace(candidate, chain, star_mass_kg, snow_line_au, formation_snow_line_au, metal_ratio, scatter, u):
    lower, upper = candidate["inner_au"], candidate["outer_au"]

    def mass(low, high):
        return _primordial(low, high, star_mass_kg, snow_line_au, metal_ratio, scatter)

    primordial = candidate.get("mass_earth") or mass(lower, upper)
    scale = primordial / mass(lower, upper) if candidate.get("mass_earth") else 1.0
    remaining = [(lower, upper)]
    history, cleared, stirred = [], 0.0, False
    for body in sorted((b for b in chain if b.giant), key=lambda b: abs(math.log(b.distance_au / (lower * upper) ** 0.5))):
        taken = 0.0
        for span in cleared_spans(body, star_mass_kg, lower, upper):
            before = sum(mass(a, b) for a, b in remaining)
            remaining = _subtract(remaining, span)
            taken += before - sum(mass(a, b) for a, b in remaining)
        taken *= scale
        reaches = is_stirred(body, star_mass_kg, lower, upper)
        stirred = stirred or reaches
        if taken > 0.0 or reaches:
            cleared += taken
            text = (_text("giant-stirring", giant_au=body.distance_au, removed_percent=100.0 * taken / primordial) if taken > 0.0
                    else CAUSES["giant-stirring"]["alt_text"].format(giant_au=body.distance_au))
            history.append(Event("giant-stirring", taken, body.distance_au, text))
    scattered = 0.0
    movers = [body for body in chain if moved_through(body, lower, upper, star_mass_kg)]
    if movers:
        survival = instability_survival(u)
        half = chaotic_half_width_au(movers[0], star_mass_kg)
        path = (min(movers[0].formation_au, movers[0].distance_au) - half, max(movers[0].formation_au, movers[0].distance_au) + half)
        swept = sum(mass(a, b) for a, b in _clip(remaining, *path)) * scale
        scattered = swept * (1.0 - survival)
        history.append(Event("instability-depletion", scattered, movers[0].distance_au,
                             _text("instability-depletion", giant_from_au=movers[0].formation_au, giant_au=movers[0].distance_au,
                                   removed_percent=100.0 * scattered / primordial)))
    left = primordial - cleared - scattered
    fields = {key: value for key, value in candidate.items() if key not in ("cause", "co", "mass_earth")}
    fields.setdefault("snow_line_au", formation_snow_line_au)
    fields.setdefault("event_au", fields.get("event_au", 0.5 * (lower + upper)))
    fields.setdefault("inner_neighbour_au", lower)
    fields.setdefault("outer_neighbour_au", upper)
    parts = [_text(candidate["cause"], **fields)]
    parts += [CAUSES[co]["co_text"].format(**{name: fields[name] for name in CAUSES[co]["co_fields"]}) for co in candidate["co"]]
    sentence = " ".join(parts + [event.text for event in history])
    return BeltOrigin(candidate["cause"], tuple(candidate["co"]), lower, upper, formation_snow_line_au, primordial, left, cleared, scattered,
                      stirred, tuple(history), sentence)


def build(chain, star_mass_kg, snow_line_au, formation_snow_line_au, growth_limit_au, disc_edge_au, u=0.5, debris=(), metal_ratio=1.0,
          scatter=1.0):
    """
    Every belt of a system with the cause that made it: the gaps, the unfinished tail, the snow-line bump (or its co-cause), the outer disc and
    the `debris` events the caller names, each traced through the giants' clearing and the instability, in order of distance. Belts under
    `MIN_BELT_MASS_EARTH` of primordial mass are not belts.

    Args:
        chain (list): The `Body` list of the system's planets, in order of distance.
        star_mass_kg (float): The star's mass at formation.
        snow_line_au (float): The snow line now (the disc's surface density reads it).
        formation_snow_line_au (float): The snow line when the system formed (stored on the belt).
        growth_limit_au (float): Where the disc built less (`tuning.PLANET_OUTER_GROWTH_AU`).
        disc_edge_au (float): The outer edge of the disc.
        u (float): A uniform draw in (0, 1) for the survival share of a belt a giant swept through.
        debris (iterable): `Debris` events, from `destroyed_planet` and `collisional_debris`.

    Returns:
        list: `BeltOrigin`, nearest the star first.
    """
    chain = sorted(chain, key=lambda body: body.distance_au)
    candidates = gaps(chain, star_mass_kg) + tail(chain, star_mass_kg, growth_limit_au) + outer_disc(chain, star_mass_kg, growth_limit_au, disc_edge_au)
    candidates += snow_line_bump(chain, star_mass_kg, formation_snow_line_au, candidates)
    candidates += [_debris_candidate(event) for event in debris]
    belts = []
    for candidate in sorted(candidates, key=lambda c: c["inner_au"]):
        belt = _trace(candidate, chain, star_mass_kg, snow_line_au, formation_snow_line_au, metal_ratio, scatter, u)
        if belt.primordial_mass_earth >= MIN_BELT_MASS_EARTH:
            belts.append(belt)
    return belts


def conservation_error(belt):
    """The relative error of `primordial = left + cleared + scattered` (0 when the mass is conserved)."""
    total = belt.mass_earth + belt.cleared_earth + belt.scattered_earth
    return abs(total - belt.primordial_mass_earth) / belt.primordial_mass_earth
