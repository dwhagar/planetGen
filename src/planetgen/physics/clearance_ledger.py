# planetgen/physics/clearance_ledger.py

"""
The clearance ledger: where in a system a planet, a moon or a belt may sit, worked out once from the stars' tree and the bodies already placed
(GEN.1014, GEN.1006 section 4.6)
===================================================================================================================================================

Today three separate pieces answer "is there room here": the star tree's own orbits, a per-star `a_crit_au` the planet code reads, and the
floor, spacing and `min_orbit_distance` rules the placer applies one body at a time after the fact. The placer draws a distance and then pushes
it out until the rules hold, so a rule it forgets to repeat (a planet after a belt is spaced against the belt but not against the planet
before it, GEN.1023) lets a pair through. The ledger answers the question first: for a host (a star, or a close pair taken as one body) it
holds the **spans** of semi-major axis where a new body of a given mass and eccentricity may sit, and the placer draws from the spans
(`restricted_quantile`: the placer's own law cut to the spans, by inverting its cumulative distribution), so no fitting is needed after.

What the spans are made of, each rule one read:

- **The host's window** (`planet_zones.planet_zones`): from the engulfment or Roche floor the caller gives, and the Holman and Wiegert (1999)
  P-type limit of a close pair, out to the tightest S-type limit against the stars beside it at every level above (`ancestor_limits_au`;
  Holman and Wiegert, Eqs 1 and 3, equal to the paper's).
  The fits are for coplanar test particles with a mass ratio of 0.1 to 0.9 (S-type; 0.1 to 0.5 for P-type) and an eccentricity up to 0.8
  (0.7 for P-type): a window worked out outside that range is applied at its edge and `Ledger.extrapolated` says so, so a caller can keep a margin.
- **Each planet already placed** excludes the stretch round it where a new planet would stand closer than `MUTUAL_HILL_RADII_SEPARATION`
  mutual Hill radii, measured between the inner body's apoapsis and the outer one's periapsis (`planet_exclusion`; the same closed form
  `validation.mutual_min_distance_au` solves, now in both directions, so the new body may go inside an existing one too).
- **Each belt** excludes its annulus, widened on the planet's side by the planet's own five Hill radii and `MIN_ASTEROID_BELT_SEPARATION`
  (`belt_exclusion`). A belt adds no spacing between the planets either side of it: both of those exclusions stay.
- **A moon's window** is the same idea one level down (`moon_bounds_km`): outside the planet and the largest moon it could hold, inside the
  fraction of the planet's Hill radius that Domingos, Winter and Yokoyama (2006, MNRAS 373, 1227, read in the paper) fit to the restricted
  elliptic three-body problem, with the planet's and the moon's eccentricity terms and the retrograde limit.

The ledger is derived data: it can be rebuilt from the tree and the placed bodies, and nothing in it is stored. It imports no generation code, so
the independent verifier (`physics/clearance.py`) can hold the generator's placements against forms it re-derives itself; `tests/test_clearance_ledger.py`
does that for planets, moons and belts.

Distances are AU (km for a moon), masses kg (solar for the tree), as the callers' own.
"""

import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.physics import planet_zones, unit_registry as units

_msun_to_kg = units.converter("Msun", "kg")

DOMINGOS_PROGRADE = (0.4895, 1.0305, 0.2738)
"""tuple: `(a, b, c)` of the prograde fit `a_crit / R_H = a (1 - b e_planet - c e_moon)` (Domingos, Winter and Yokoyama 2006, MNRAS 373, 1227, read
in the paper: restricted elliptic three-body problem, planet orbits tested at 0.1 AU with `e_planet` 0 to 0.9). Source tier T1, tweak level L0."""
DOMINGOS_RETROGRADE = (0.9309, 1.0764, 0.9812)
"""tuple: The same fit for a retrograde moon (Domingos et al. 2006). Source tier T1, tweak level L0."""
DOMINGOS_PLANET_E_MAX = 0.9
"""float: The planet eccentricity the fit was made up to; past it the limit is applied at this edge and flagged `extrapolated`."""
HOLMAN_WIEGERT_S_MU = (0.1, 0.9)
"""tuple: The companion's mass fraction the S-type fit covers (Holman and Wiegert 1999)."""
HOLMAN_WIEGERT_P_MU = (0.1, 0.5)
"""tuple: The lighter side's mass fraction the P-type fit covers."""
HOLMAN_WIEGERT_S_E_MAX = 0.8
"""float: The eccentricity the S-type fit covers."""
HOLMAN_WIEGERT_P_E_MAX = 0.7
"""float: The eccentricity the P-type fit covers."""
KAPPA_MAX = 1.8
"""float: The cap on the spacing `MUTUAL_HILL_RADII_SEPARATION * ((m1 + m2) / 3M)^(1/3)` (the linear model has no finite answer at 2; the generator's own cap)."""
BELT_HILL_RADII = 5.0
"""float: How many of the planet's own Hill radii it keeps from a belt (the generator's `min_orbit_distance`)."""
BELT_SHARE_CAP = 0.9
"""float: The cap on the Hill share of a planet's distance that a belt exclusion allows (the generator's own)."""
ECCENTRICITY_DIVISOR_FLOOR = 0.1
"""float: The floor on `1 - e - kappa/2` so a very eccentric heavy pair still gets a finite exclusion (the generator's own)."""


@dataclass(frozen=True)
class Span:
    """One stretch of semi-major axis, `low < high`."""
    low: float
    high: float

    @property
    def width(self):
        return self.high - self.low


class SpanSet:
    """
    A set of disjoint spans, sorted. Immutable: every operation returns a new set. Empty spans are dropped and touching ones joined.
    """

    def __init__(self, spans=()):
        ordered = sorted((Span(float(low), float(high)) for low, high in _pairs(spans) if high > low), key=lambda span: span.low)
        merged = []
        for span in ordered:
            if merged and span.low <= merged[-1].high:
                if span.high > merged[-1].high:
                    merged[-1] = Span(merged[-1].low, span.high)
            else:
                merged.append(span)
        self.spans = tuple(merged)

    def __iter__(self):
        return iter(self.spans)

    def __len__(self):
        return len(self.spans)

    def __bool__(self):
        return bool(self.spans)

    def __repr__(self):
        return "SpanSet(" + ", ".join(f"({span.low:.6g}, {span.high:.6g})" for span in self.spans) + ")"

    def __eq__(self, other):
        return isinstance(other, SpanSet) and self.spans == other.spans

    @property
    def width(self):
        return sum(span.width for span in self.spans)

    def contains(self, value):
        """True if `value` lies in a span (the ends count)."""
        return any(span.low <= value <= span.high for span in self.spans)

    def covers(self, low, high):
        """True if one span holds the whole of `[low, high]`."""
        return any(span.low <= low and high <= span.high for span in self.spans)

    def subtract(self, low, high):
        """The set with the open stretch `(low, high)` taken out (a body may sit at the edge of an exclusion)."""
        if not high > low:
            return self
        out = []
        for span in self.spans:
            if span.high <= low or span.low >= high:
                out.append((span.low, span.high))
                continue
            if span.low < low:
                out.append((span.low, low))
            if span.high > high:
                out.append((high, span.high))
        return SpanSet(out)

    def intersect(self, low, high):
        """The set clipped to `[low, high]`."""
        return SpanSet((max(span.low, low), min(span.high, high)) for span in self.spans)

    def nearest(self, value):
        """The point of the set nearest `value`, or None for an empty set."""
        best = None
        for span in self.spans:
            point = min(max(value, span.low), span.high)
            if best is None or abs(point - value) < abs(best - value):
                best = point
        return best


def _pairs(spans):
    for span in spans:
        if isinstance(span, Span):
            yield span.low, span.high
        else:
            yield span


def restricted_quantile(u, spans, cdf, ppf):
    """
    The `u` quantile (0 to 1) of a law cut to `spans`: the exact conditional distribution, in one draw with no retry. The law's cumulative
    distribution `cdf` is taken at both ends of every span, `u` picks a point of the total probability the spans hold, and `ppf` (the law's
    inverse) turns it back into a distance.

    Args:
        u (float): A uniform draw in [0, 1).
        spans (SpanSet): Where the body may sit.
        cdf (callable): The law's cumulative distribution over the distance.
        ppf (callable): Its inverse, over [0, 1].

    Returns:
        float or None: The distance, or None when the spans hold none of the law's probability (no room: a counted event for the caller).
    """
    weights = [(cdf(span.high) - cdf(span.low), span) for span in spans]
    total = sum(weight for weight, _span in weights if weight > 0.0)
    if not total > 0.0:
        return None
    target = min(max(u, 0.0), 1.0) * total
    base = 0.0
    for weight, span in weights:
        if weight <= 0.0:
            continue
        if target <= base + weight or span is weights[-1][1]:
            return min(max(ppf(cdf(span.low) + (target - base)), span.low), span.high)
        base += weight
    return None


def log_uniform_quantile(u, spans):
    """`restricted_quantile` for a law uniform in the logarithm of the distance (the placer's default), over spans that start above 0."""
    return restricted_quantile(u, spans.intersect(1e-300, math.inf), math.log, math.exp)


# --- planets and belts ------------------------------------------------------------------------------------------------------------------

def spacing_kappa(mass_kg, other_mass_kg, central_mass_kg):
    """`MUTUAL_HILL_RADII_SEPARATION` mutual Hill radii as a fraction of the mean distance of two bodies of these masses about `central_mass_kg`."""
    return min(tuning.MUTUAL_HILL_RADII_SEPARATION * ((mass_kg + other_mass_kg) / (3.0 * central_mass_kg)) ** (1.0 / 3.0), KAPPA_MAX)


def planet_exclusion(distance_au, eccentricity, mass_kg, central_mass_kg, new_mass_kg, new_eccentricity=0.0):
    """
    The stretch of semi-major axis round a planet where a new planet of `new_mass_kg` and `new_eccentricity` would stand closer than
    `MUTUAL_HILL_RADII_SEPARATION` mutual Hill radii (Gladman 1993; Chambers, Wetherill and Boss 1996, covered through Obertas, Van
    Laerhoven and Tamayo 2017), measured between the inner body's apoapsis and the outer one's periapsis. With the mutual radius taken at the
    mean distance, the gap `x (1 - e_new) - d (1 + e) >= kappa (x + d) / 2` for a body outside gives `x >= d (1 + e + kappa/2) / (1 - e_new -
    kappa/2)` (the generator's closed form), and for one inside `x <= d (1 - e - kappa/2) / (1 + e_new + kappa/2)`.

    Returns:
        tuple: `(low_au, high_au)`, the open stretch a new planet may not sit in; `low_au` is 0 when the planet leaves no room inside it
            (its own periapsis is under half its spacing).
    """
    kappa = spacing_kappa(mass_kg, new_mass_kg, central_mass_kg)
    outer = distance_au * (1.0 + eccentricity + kappa / 2.0) / max(1.0 - new_eccentricity - kappa / 2.0, ECCENTRICITY_DIVISOR_FLOOR)
    room_inside = 1.0 - eccentricity - kappa / 2.0
    inner = distance_au * room_inside / (1.0 + new_eccentricity + kappa / 2.0) if room_inside > 0.0 else 0.0
    return inner, outer


def belt_exclusion(lower_au, upper_au, new_mass_kg, central_mass_kg, new_eccentricity=0.0):
    """
    The stretch of semi-major axis a new planet may not take because of a belt on `[lower_au, upper_au]`: its Hill sphere (five Hill radii, the
    generator's `min_orbit_distance`) must clear the belt, and past the belt's outer edge the `MIN_ASTEROID_BELT_SEPARATION` too. The Hill
    radius scales with the distance, so as in `validation.min_distance_past_belt_au` it is solved in closed form: outside `x (1 - e - 5c) >=
    upper + gap`, inside `x (1 + e + 5c) <= lower`, with `c = (m / 3M)^(1/3)`.

    Returns:
        tuple: `(low_au, high_au)`, the open stretch a new planet may not sit in.
    """
    share = min(BELT_HILL_RADII * (new_mass_kg / (3.0 * (central_mass_kg + new_mass_kg))) ** (1.0 / 3.0), BELT_SHARE_CAP)
    outer = (upper_au + tuning.MIN_ASTEROID_BELT_SEPARATION) / max(1.0 - share - new_eccentricity, ECCENTRICITY_DIVISOR_FLOOR)
    inner = lower_au / (1.0 + share + new_eccentricity)
    return inner, outer


# --- moons ------------------------------------------------------------------------------------------------------------------------------

def moon_hill_fraction(planet_eccentricity, moon_eccentricity=0.0, retrograde=False):
    """
    The largest semi-major axis a moon may have, in units of its planet's Hill radius (taken at the planet's semi-major axis), from the
    Domingos, Winter and Yokoyama (2006) fit `a (1 - b e_planet - c e_moon)`.

    Returns:
        tuple: `(fraction, extrapolated)`; the fraction is held at 0 or more, and `extrapolated` is True past the planet eccentricity the fit covers.
    """
    a, b, c = DOMINGOS_RETROGRADE if retrograde else DOMINGOS_PROGRADE
    extrapolated = planet_eccentricity > DOMINGOS_PLANET_E_MAX
    e_planet = min(max(planet_eccentricity, 0.0), DOMINGOS_PLANET_E_MAX)
    return max(a * (1.0 - b * e_planet - c * max(moon_eccentricity, 0.0)), 0.0), extrapolated


def moon_max_eccentricity(semi_major_axis_au, hill_radius_au, planet_eccentricity=0.0, retrograde=False):
    """
    The largest eccentricity a moon at `semi_major_axis_au` may have: the Domingos fit solved for the moon's term, `a_moon <= R_H a (1 - b e_planet
    - c e_moon)`. Zero for a moon already at the limit.

    Returns:
        float or None: `None` when the moon's eccentricity does not enter the limit (it has none to give) or the Hill radius is unknown.
    """
    if hill_radius_au is None or semi_major_axis_au <= 0.0:
        return None
    a, b, c = DOMINGOS_RETROGRADE if retrograde else DOMINGOS_PROGRADE
    e_planet = min(max(planet_eccentricity, 0.0), DOMINGOS_PLANET_E_MAX)
    return max((1.0 - b * e_planet - semi_major_axis_au / (a * hill_radius_au)) / c, 0.0)


def moon_bounds_km(planet_radius_km, hill_radius_km, planet_eccentricity, moon_eccentricity=0.0, retrograde=False, margin_km=0.0,
                   largest_moon_radius_km=0.0):
    """
    The edges `(low_km, high_km)` of where a moon may orbit its planet, km from the planet's centre: outside the planet, the largest moon it could
    hold and the caller's margin (the air's drag-free height), and inside the Domingos limit for the planet's and the moon's eccentricity (the fit
    is for the moon's semi-major axis, its own eccentricity term already making the limit lower). `low_km >= high_km` means no room.

    Returns:
        tuple: `(low_km, high_km, extrapolated)`.
    """
    fraction, extrapolated = moon_hill_fraction(planet_eccentricity, moon_eccentricity, retrograde)
    return planet_radius_km + largest_moon_radius_km + margin_km, hill_radius_km * fraction, extrapolated


# --- the ledger of one system ----------------------------------------------------------------------------------------------------------

@dataclass
class _Placed:
    distance_au: float
    eccentricity: float
    mass_kg: float


class Ledger:
    """
    The ledger of one system's stars: every host's window, and the bodies placed in each so far.

    Args:
        nodes (iterable): The tree's nodes (empty for a lone object), as `planet_zones.planet_zones` takes them.
        members (iterable): Its members.
        floors_au (dict, optional): `{member_no: AU}` inside which a member has removed every planet (its engulfment or Roche floor).
        disk_cap_au (float, optional): The disk size; `tuning.PLANET_ZONE_DISK_CAP_AU` when left out.
        ceilings_au (dict, optional): `{(kind, host_no): AU}` an outer edge tighter than the window's (the star's Hill sphere against the galaxy, the
            disk edge), taken with the window.
    """

    def __init__(self, nodes, members, floors_au=None, disk_cap_au=None, ceilings_au=None):
        self.nodes = list(nodes)
        self.members = list(members)
        self._zones = {(zone.kind, zone.host_no): zone
                       for zone in planet_zones.planet_zones(self.nodes, self.members, floors_au, disk_cap_au)}
        self._ceilings = dict(ceilings_au or {})
        self._planets = {}
        self._belts = {}

    # -- reading the windows

    def hosts(self):
        """The `(kind, host_no)` of every host that has a window, stars first."""
        return sorted(self._zones, key=lambda host: (host[0] != planet_zones.S_TYPE, host[1]))

    def zone(self, host):
        return self._zones[host]

    def zones(self):
        """Every window of the tree as `planet_zones.Zone`, the S-type ones first."""
        return [self._zones[host] for host in self.hosts()]

    def central_mass_kg(self, host):
        """The mass a body of `host` goes round, kg."""
        return _msun_to_kg(planet_zones.host_mass_solar(self._zones[host], self.nodes, self.members))

    def outer_au(self, host):
        """The outermost orbit the host holds: the tightest of its window's limit, the caller's ceiling and every S-type limit above it."""
        zone = self._zones[host]
        return min([zone.outer_au, self._ceilings.get(host, math.inf)] + planet_zones.ancestor_limits_au(zone, self.nodes, self.members))

    def stability_limit_au(self, host):
        """
        The tightest S-type limit above the host (Holman and Wiegert against the stars beside it at every level, or the host's Hill radius
        in the pair), or None for a host nothing hangs beside (a lone star, the root of the tree). Unlike `outer_au` it holds no disk size and
        no caller ceiling: it is the distance past which no orbit is stable, which `Star.a_crit_au` stores.
        """
        limits = planet_zones.ancestor_limits_au(self._zones[host], self.nodes, self.members)
        return min(limits) if limits else None

    def window(self, host):
        """The host's own window, with its ceiling, as a `SpanSet` (empty when the limits leave none)."""
        return SpanSet([(self._zones[host].inner_au, self.outer_au(host))])

    def extrapolated(self, host):
        """True if the host's window rests on a Holman and Wiegert fit outside the range it was made on (applied at the edge of that range)."""
        zone = self._zones[host]
        if zone.kind == planet_zones.S_TYPE:
            if not any(member.member_no == zone.host_no and member.parent_node_no is not None for member in self.members):
                return False
            node, beside = planet_zones.sibling_of(self.nodes, self.members, zone.host_no)
            own = {member.member_no: member.mass_solar for member in self.members}
            other = sum(own[number] for number in beside)
            return _mu_outside(other / (other + own[zone.host_no]), HOLMAN_WIEGERT_S_MU) or node.eccentricity > HOLMAN_WIEGERT_S_E_MAX
        node = next(node for node in self.nodes if node.node_no == zone.host_no)
        low, high = planet_zones.pair_sides_solar(self.nodes, self.members, zone.host_no)
        return _mu_outside(min(low, high) / (low + high), HOLMAN_WIEGERT_P_MU) or node.eccentricity > HOLMAN_WIEGERT_P_E_MAX

    # -- putting bodies in

    def place_planet(self, host, distance_au, mass_kg, eccentricity=0.0):
        """Records a planet in `host`; later `open_spans` calls leave room round it."""
        self._planets.setdefault(host, []).append(_Placed(distance_au, eccentricity, mass_kg))

    def place_belt(self, host, lower_au, upper_au):
        """Records a belt in `host`."""
        self._belts.setdefault(host, []).append((lower_au, upper_au))

    def planets(self, host):
        return list(self._planets.get(host, ()))

    def belts(self, host):
        return list(self._belts.get(host, ()))

    # -- asking for room

    def open_spans(self, host, mass_kg, eccentricity=0.0):
        """
        Where a new planet of `mass_kg` and `eccentricity` may sit in `host`: its window, less the exclusion round every planet and every belt
        already placed. Both neighbours of a gap are held (a belt between two planets does not release the planets from each other).
        """
        spans = self.window(host)
        central = self.central_mass_kg(host)
        for placed in self._planets.get(host, ()):
            low, high = planet_exclusion(placed.distance_au, placed.eccentricity, placed.mass_kg, central, mass_kg, eccentricity)
            spans = spans.subtract(low, high)
        for lower, upper in self._belts.get(host, ()):
            low, high = belt_exclusion(lower, upper, mass_kg, central, eccentricity)
            spans = spans.subtract(low, high)
        return spans

    def open_belt_spans(self, host, width_au=0.0):
        """
        Where a belt of `width_au` may go in `host`: its window, less the stretch that holds each planet's orbit and the planet's five Hill radii
        either side, and less each other belt with `MIN_ASTEROID_BELT_SEPARATION` round it. The lower edge of a span is where a belt's lower edge may
        start; a belt that must fit whole needs `Span.width >= width_au`.
        """
        spans = self.window(host)
        central = self.central_mass_kg(host)
        for placed in self._planets.get(host, ()):
            reach = BELT_HILL_RADII * placed.distance_au * (placed.mass_kg / (3.0 * (central + placed.mass_kg))) ** (1.0 / 3.0)
            spans = spans.subtract(placed.distance_au * (1.0 - placed.eccentricity) - reach, placed.distance_au * (1.0 + placed.eccentricity) + reach)
        for lower, upper in self._belts.get(host, ()):
            spans = spans.subtract(lower - tuning.MIN_ASTEROID_BELT_SEPARATION, upper + tuning.MIN_ASTEROID_BELT_SEPARATION)
        return SpanSet((span.low, span.high) for span in spans if span.width >= width_au)


def _mu_outside(mu, bounds):
    return not bounds[0] <= mu <= bounds[1]
