# planetgen/physics/belt_origin.py

"""
Where a belt came from (GEN.1088, docs/design/systems/belt-origins.md sections 3 to 5)
=====================================================================================

Every belt and outer disk gets the cause that made it, a record of how its solids were formed and what took mass from them, and one sentence
for its page ("The leftover of a planet that never grew here ..."). The generator made the belt, so it knows. A pure law over the planet chain
the architecture (`planet_architecture`) built and the disk's solids; it reads no database, draws no random number (the two draws it needs come
in as `u` arguments) and moves nothing. GEN.1068 gives it the final chain in place of the 10 percent slot roll.

**Formation causes** (what made the belt; the stored `origin`):

- `gap`: two neighbours of the chain stand more than `GAP_SPACING_K` mutual Hill radii apart (the thinned chains cluster at 12 to 20), and the
  space between them, less five Hill radii of each, is wide enough: the solids there did not grow into a planet.
- `tail`: the chain stops short of `growth_limit_au`, where the disk has built less (`PLANET_OUTER_GROWTH_AU`), and the space past its last
  rocky planet is a belt (a last giant ends the chain, the outer disk takes over).
- `snow-line-bump`: solids pile up at the snow line as the system formed it (the formation snow line is stored); where that stretch is free of
  planets and is not already inside a gap or tail (which then take it as a second cause) it is its own belt.
- `outer-disk`: the disk past the last giant (or past the tail), out to `disk_edge_au`.
- `destroyed-planet`: the remains of a planet the system's history destroyed (the caller names it), `DESTROYED_RETAINED_SHARE` of its mass.
- `collisional-debris`: the rubble of a collision between large bodies, a rare event at `COLLISIONAL_DEBRIS_PER_GYR` of mass
  `COLLISIONAL_DEBRIS_MASS_EARTH`; `collisional_debris` turns the caller's draw into the event or none.

**Shaping causes** (what took mass afterwards; the stored `history`), applied in this order to each belt, solids weighted by the disk's own
surface density (`belt_mass.primordial_mass_earth`):

- `giant-stirring`: the chaotic zone of each giant, half-width `1.3 a mu^(2/7)` each side [R4727], and for a belt inside the giant the
  gaps at the resonances 4:1, 3:1, 5:2, 7:3 and 2:1 (`a_g (q/p)^(2/3)`), half-width `RESONANCE_HALF_WIDTH` of the resonance's distance for a
  Jupiter-mass giant and as the square root of the giant's mass ratio. What is in them is gone: `cleared_earth`. A belt within
  `STIR_ZONE_WIDTHS` chaotic half-widths of a giant is stirred whether or not mass is cleared, and says so.
- `instability-depletion`: a giant that moved (the architecture's hot Jupiter, `migrated`) through the belt's stretch swept that part of it: a share
  `INSTABILITY_SURVIVAL_RANGE` of the swept stretch survives, log-uniform in the draw `u`; the rest is `scattered_earth` (the scatter step, GEN.1089,
  decides where it went).

**What a stretch starts with** (GEN.1261): the solids of the disk in it, less the cold disk's drift loss (`planetesimal_share`: past the snow line only the share of solids that grew to planetesimals before the pebbles drifted in), and for a gap no more than the isolation mass at its center (`isolation_cap_earth`: more would have grown the planet). That is the `primordial_mass_earth` below.

Mass is conserved: `primordial_mass_earth = mass_earth + cleared_earth + scattered_earth` (`conservation_error`). The mass later lost to
collisional grinding is not here: that is `belt_mass.belt_mass_earth`, which needs the age and the largest body.

Where the design says the stored origin is "the dominant" cause, this module stores the formation cause as `origin` (the one that made the
belt) and every shaping cause with the share of the mass it took in `history`, so the page can say both.

Use:

    chain = [Body(distance_au, mass_kg, eccentricity, giant, migrated, formation_au), ...]
    belts = build(chain, star_mass_kg, snow_line_au, formation_snow_line_au, growth_limit_au, disk_edge_au, u=0.5)
    for belt in belts: belt.origin, belt.mass_earth, belt.sentence
"""

import collections
import json
import math
import os

from planetgen import tuning
from planetgen.physics import belt_mass, clearance_ledger, elements, formation, planet_architecture, planet_mass, scatter as scatter_rule, unit_registry as units
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

JUPITER_RADIUS_KM = 71492.0
"""float: The radius of a giant that has none stored, km (Jupiter's; a giant's radius changes little with its mass, row BG21)."""

DESTROYED_RETAINED_SHARE = 0.5
"""float: The share of a destroyed planet's mass that stays in orbit as a belt (BG10)."""

DEBRIS_MIN_WIDTH = 0.05
"""float: The least half-width of the belt of debris round its distance, as a share of the distance (BG11)."""

PEBBLE_CONVERSION_PER_ORBIT = 0.01
"""float: The share of the pebbles of a patch of the disk that the streaming instability turns into planetesimals each orbit, where it works (pebble-to-planetesimal simulations set the rate at 0.01 of the pebble surface density per orbit). [Stand-in BG17, tier Confident Unverified, tweak L0.]"""

PEBBLE_TO_GAS_SOLAR = 0.01
"""float: The solids-to-gas ratio of a disk of solar metal content, the interstellar value; it scales with the star's `metal_ratio`. [Stand-in BG18, tier Confident Unverified, tweak L0.]"""

PEBBLE_DRIFT_FACTOR = 1.1
"""float: The pebbles of a drift-limited disk fall in over `1 / (PEBBLE_DRIFT_FACTOR Z Omega)`: their Stokes number is the drift limit `0.55 Z / (eta)` of the drift-limited dust models
and the drift speed `2 St eta v_K`, so the pressure gradient `eta` drops out and the time is 1.1 inverse growth rates, whatever the radius. [Stand-in BG19, tier Marginal educated guess from physics, tweak L0.]"""


def planetesimal_share(metal_ratio=1.0):
    """
    The share of the solids of the cold disk that grew into planetesimals before the pebbles drifted in (GEN.1261). Pebbles turn into planetesimals at
    `PEBBLE_CONVERSION_PER_ORBIT` of their mass an orbit and fall in over `1 / (PEBBLE_DRIFT_FACTOR Z Omega)` (Z the solids-to-gas ratio, `PEBBLE_TO_GAS_SOLAR` times
    `metal_ratio`); the two compete, so the share that becomes planetesimals is `zeta / (zeta + PEBBLE_DRIFT_FACTOR Z)`, the same at every radius because both rates scale with
    the orbital frequency. The rest is the drift loss ("the majority of the dust reservoir is lost through drift", a pebble-drift disk model, read second-hand).
    """
    zeta = PEBBLE_CONVERSION_PER_ORBIT
    return zeta / (zeta + PEBBLE_DRIFT_FACTOR * PEBBLE_TO_GAS_SOLAR * max(metal_ratio, 0.0))


def isolation_cap_earth(lower, upper, star_mass_kg, snow_line_au, metal_ratio=1.0, scatter=1.0, separation=None):
    """
    The most solids a planetless stretch can hold (GEN.1261): the isolation mass a planet at its center would have taken from its feeding zone (`formation.isolation_mass_kg`,
    [R2430]) in the disk's own solids there. A stretch holding more would have grown the planet, so the chain would not have the gap.
    """
    center = math.sqrt(lower * upper)
    density = planet_mass.solids_surface_density_gcm2(center, star_mass_kg, snow_line_au, metal_ratio, scatter)
    return _kg_to_mearth(formation.isolation_mass_kg(center, density, star_mass_kg, tuning.OLIGARCH_SEPARATION_MEDIAN if separation is None else separation))


COLLISIONAL_DEBRIS_PER_GYR = 0.003
COLLISIONAL_DEBRIS_MASS_EARTH = 1.0e-4
"""float: The chance per Gyr of a standing belt of rubble from a collision of large bodies, and its mass (BG12 and BG13)."""

Body = collections.namedtuple("Body", ("distance_au", "mass_kg", "eccentricity", "giant", "migrated", "formation_au", "radius_km"), defaults=(0.0, False, False, None, None))
"""One body of the planet chain: where it ends, its mass, its eccentricity, whether it is a giant, whether it moved in (the architecture's hot
Jupiter), where it formed and its radius (km; `JUPITER_RADIUS_KM` for a giant that has none)."""

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
    """The unfinished tail: past the last rocky planet to `growth_limit_au`, only when no giant ends the chain (the outer disk takes that one)."""
    if not chain or chain[-1].giant:
        return []
    last = chain[-1]
    lower, upper = keep_out(last, star_mass_kg)[1], growth_limit_au
    return [_candidate("tail", lower, upper, inner_neighbour_au=last.distance_au)] if _fits(lower, upper) else []


def outer_disk(chain, star_mass_kg, growth_limit_au, disk_edge_au):
    """
    The outer disk: from the outer chaotic edge of the last giant and the keep-out of the last body, and from `growth_limit_au` too when the chain
    ends in a rocky planet (the tail has the space before it), out to `disk_edge_au`.
    """
    giants = [body for body in chain if body.giant]
    lower = growth_limit_au if not chain or not chain[-1].giant else 0.0
    if chain:
        lower = max(lower, keep_out(chain[-1], star_mass_kg)[1])
    if giants:
        lower = max(lower, giants[-1].distance_au + chaotic_half_width_au(giants[-1], star_mass_kg))
    neighbour = giants[-1].distance_au if giants else (chain[-1].distance_au if chain else growth_limit_au)
    return [_candidate("outer-disk", lower, disk_edge_au, inner_neighbour_au=neighbour)] if _fits(lower, disk_edge_au) else []


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

def scattering_reach_au(body, star_mass_kg):
    """How far out the giants of a system that scattered each other sweep, from the outermost one `body`: its apoapsis and the `STIR_ZONE_WIDTHS` chaotic half-widths beyond that it stirs (BG8, GEN.1341)."""
    return body.distance_au * (1.0 + body.eccentricity) + STIR_ZONE_WIDTHS * chaotic_half_width_au(body, star_mass_kg)


def scattering_loss_share(chain, star_mass_kg, lower_au, instability):
    """
    The share of the disk swept by the giants' scattering of each other that is thrown out (GEN.1341), for a stretch that starts at `lower_au`, as `(share, reach_au, giant)`, or
    `None` where nothing is: no instability, fewer than two giants, or a stretch that starts inside the outermost giant (the giants' own clearing takes it) or beyond the reach.

    The giants that scatter each other (`instability`, 0 to 1, `capture_windows.instability_level`) throw out the planetesimals they meet where the Safronov number of the outermost giant
    `Theta = q a / R` is large ([R5138]: above 1 the encounters end in ejection, below 1 in collision): the thrown share is `1 - exp(-Theta / THETA_EJECT)` of what lies
    between the giant and `scattering_reach_au`, scaled by the instability. The same rule the scatter step uses for the disk the giants sweep (`scatter.py`, rule 4), and with their own radii the
    survival is 0.35 for Jupiter, 0.51 for Saturn, 0.61 for Uranus and 0.39 for Neptune, about the 0.1 to 0.5 that instability runs leave (row BR5, [R5139]), with no number fitted to it. Row BG20.
    """
    giants = [body for body in chain if body.giant]
    if instability <= 0.0 or len(giants) < 2:
        return None
    outer = max(giants, key=lambda body: body.distance_au)
    reach = scattering_reach_au(outer, star_mass_kg)
    if lower_au < outer.distance_au or lower_au >= reach:
        return None
    radius_km = outer.radius_km or JUPITER_RADIUS_KM
    theta = scatter_rule.safronov_number(_kg_to_mearth(outer.mass_kg), radius_km / 6371.0, outer.distance_au, units.convert(star_mass_kg, "kg", "Msun"))
    return min(instability, 1.0) * (1.0 - math.exp(-theta / scatter_rule.THETA_EJECT)), reach, outer


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
    return bool(upper_au >= body.distance_au - reach and lower_au <= body.distance_au + reach)


def moved_through(body, lower_au, upper_au, star_mass_kg):
    """Whether a giant that moved (`migrated`) passed through the stretch: between where it formed and where it ended, a chaotic half-width each side."""
    if not (body.giant and body.migrated and body.formation_au):
        return False
    half = chaotic_half_width_au(body, star_mass_kg)
    path = (min(body.formation_au, body.distance_au) - half, max(body.formation_au, body.distance_au) + half)
    return bool(upper_au >= path[0] and lower_au <= path[1])


def instability_survival(u):
    """The share a moving giant leaves, log-uniform over `INSTABILITY_SURVIVAL_RANGE` in the draw `u` (0 to 1)."""
    low, high = INSTABILITY_SURVIVAL_RANGE
    return math.exp(math.log(low) + min(max(u, 0.0), 1.0) * (math.log(high) - math.log(low)))


# --- the belt ---------------------------------------------------------------------------------------------------------------------------

def _text(cause, **fields):
    return CAUSES[cause]["text"].format(**{name: fields[name] for name in CAUSES[cause]["fields"]})


def _trace(candidate, chain, star_mass_kg, snow_line_au, formation_snow_line_au, metal_ratio, scatter, u, instability=0.0):
    lower, upper = candidate["inner_au"], candidate["outer_au"]

    cold_share = planetesimal_share(metal_ratio)

    def mass(low, high):
        """The planetesimals of a stretch: the solids of the disk, less the cold disk's drift loss (past the snow line only `cold_share` of them grew to planetesimals)."""
        warm_high = min(high, formation_snow_line_au)
        warm = _primordial(low, warm_high, star_mass_kg, snow_line_au, metal_ratio, scatter) if warm_high > low else 0.0
        cold_low = max(low, formation_snow_line_au)
        cold = _primordial(cold_low, high, star_mass_kg, snow_line_au, metal_ratio, scatter) if high > cold_low else 0.0
        return warm + cold_share * cold

    whole = mass(lower, upper)
    primordial = candidate.get("mass_earth") or whole
    if candidate["cause"] == "gap" and not candidate.get("mass_earth"):
        primordial = min(primordial, isolation_cap_earth(lower, upper, star_mass_kg, snow_line_au, metal_ratio, scatter))
    scale = primordial / whole if whole > 0.0 else 1.0
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
    thrown = scattering_loss_share(chain, star_mass_kg, lower, instability)
    if thrown is not None and thrown[0] > 0.0:
        share, reach, giant = thrown
        swept = sum(mass(a, b) for a, b in _clip(remaining, lower, min(upper, reach))) * scale
        lost = min(swept * share, primordial - cleared - scattered)
        scattered += lost
        history.append(Event("planet-scattering", lost, giant.distance_au,
                             _text("planet-scattering", giant_au=giant.distance_au, reach_au=reach, removed_percent=100.0 * lost / primordial)))
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


def build(chain, star_mass_kg, snow_line_au, formation_snow_line_au, growth_limit_au, disk_edge_au, u=0.5, debris=(), metal_ratio=1.0,
          scatter=1.0, instability=0.0):
    """
    Every belt of a system with the cause that made it: the gaps, the unfinished tail, the snow-line bump (or its co-cause), the outer disk and
    the `debris` events the caller names, each traced through the giants' clearing and the instability, in order of distance. Belts under
    `MIN_BELT_MASS_EARTH` of primordial mass are not belts.

    Args:
        chain (list): The `Body` list of the system's planets, in order of distance.
        star_mass_kg (float): The star's mass at formation.
        snow_line_au (float): The snow line now (the disk's surface density reads it).
        formation_snow_line_au (float): The snow line when the system formed (stored on the belt).
        growth_limit_au (float): Where the disk built less (`tuning.PLANET_OUTER_GROWTH_AU`).
        disk_edge_au (float): The outer edge of the disk.
        u (float): A uniform draw in (0, 1) for the survival share of a belt a giant swept through.
        debris (iterable): `Debris` events, from `destroyed_planet` and `collisional_debris`.
        instability (float): 0 to 1, how much the giants scattered each other (`capture_windows.instability_level`): it throws out a share of the disk past the outermost giant (`scattering_loss_share`, GEN.1341).

    Returns:
        list: `BeltOrigin`, nearest the star first.
    """
    chain = sorted(chain, key=lambda body: body.distance_au)
    candidates = gaps(chain, star_mass_kg) + tail(chain, star_mass_kg, growth_limit_au) + outer_disk(chain, star_mass_kg, growth_limit_au, disk_edge_au)
    candidates += snow_line_bump(chain, star_mass_kg, formation_snow_line_au, candidates)
    candidates += [_debris_candidate(event) for event in debris]
    belts = []
    for candidate in sorted(candidates, key=lambda c: c["inner_au"]):
        belt = _trace(candidate, chain, star_mass_kg, snow_line_au, formation_snow_line_au, metal_ratio, scatter, u, instability)
        if belt.primordial_mass_earth >= MIN_BELT_MASS_EARTH:
            belts.append(belt)
    return belts


def conservation_error(belt):
    """The relative error of `primordial = left + cleared + scattered` (0 when the mass is conserved)."""
    total = belt.mass_earth + belt.cleared_earth + belt.scattered_earth
    return abs(total - belt.primordial_mass_earth) / belt.primordial_mass_earth
