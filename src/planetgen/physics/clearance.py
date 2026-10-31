# planetgen/physics/clearance.py

"""
The independent overlap verifier (TEST.477, GEN.1006 section 4.7): `audit(rows)` reads only what a finished system
stores (its stars, the nodes of its tree of orbits, its planets, moons and belts) and lists every place two bodies overlap
or sit where the physics does not leave room for them.

It is the proof of "systems cannot overlap by construction", so it must not share a bug with the generator it checks. It
imports no placement, build or screening code (nothing from `planetgen.generation`, `physics.stability`, `physics.roche`,
`physics.system_tree`, `physics.kepler` or `physics.orbits`), only the unit-registry constants, and it re-derives every
limit from the published form. `tests/test_clearance.py` holds each form against the generator's own copy; a form that
drifts from it fails there.

What it checks, for every pair of bodies that the rows hold. A **hard** violation is a place where two bodies overlap or a limit
the generator itself promises is broken; an **advisory** is a physical limit the generator does not promise to keep today
(GEN.1006 parts P4 to P8 are what turn them into promises), counted so the work to come has its baseline.

1. **Stars (`roche`, hard)**: two single members of a node must clear each other at the pair's closest approach: neither
   fills more than `MAX_FILL_FACTOR` of its Roche lobe there (Eggleton 1983), except a pair the record marks as a merger
   product or a post-common-envelope pair.
2. **Subtrees (`overlap`, hard)**: for a node with a sub-node as a child, the pericentre must exceed the sum of its
   children's extents (a star's radius; a sub-node's greatest distance from its own barycentre to the surface of any star
   below it, at the node's apocentre): the bounding spheres of the two sides do not meet.
3. **Stability (`stability`, hard)**: for every node with a sub-node as a child, the triple of the sub-node's two sides and the
   other child must pass at least one of Mardling and Aarseth 2001, Eggleton and Kiseleva 1995 or Equation 4 of
   Vynatheya et al. 2022 (a triple the generator drew by any of them passes its own; one that fails all three is unstable
   whichever it used).
4. **Masses (`mass`, hard)**: each node's mass is the sum of its two children's.
5. **Planets**: the closest approach is outside the host star (`planet-in-star`, hard) and, if the density says so, outside its
   rigid-body Roche limit (`planet-roche`, advisory); the semi-major axis lies inside the S-type critical axis of its host
   toward each companion above it (`s-type`, hard, Holman and Wiegert 1999, equation 1) or outside the P-type one of the pair
   it goes round (`p-type`, hard, equation 3; both are fits for the semi-major axis of a test particle); two planets of one host
   do not cross (`planet-crossing`, hard) and are `MUTUAL_HILL_RADII_SEPARATION` mutual Hill radii apart (`planet-spacing`, hard).
6. **Moons**: inside 0.4895 of the planet's Hill radius (`moon-hill`, hard, Domingos, Winter and Yokoyama 2006), not touching the
   planet (`moon-in-planet`, hard) and outside its rigid-body Roche limit for the two densities (`moon-roche`, advisory).
7. **Belts**: the annulus does not hold a planet's orbit (`belt-planet`, hard) and lies inside the host's S-type critical axis
   (`belt-s-type`, hard).

Not checked yet, because the rows do not hold what they need: the history rule (a star that was a giant keeps its largest
radius only in the generator's memory, GEN.1006 4.3, which waits for the member record's `radius_max_au`), the Kozai-Lidov
screen and the 3+1 quadruple margin (GEN.1006 section 5). Every limit that is outside the range its fit was made on
(Holman and Wiegert's mass ratio 0.1 to 0.9 and eccentricity up to 0.8 for the S-type axis, 0.1 to 0.5 and 0.7 for the P-type one) is applied at the nearest edge of that range and the
violation says so (`extrapolated`).
"""

import math
from dataclasses import dataclass, field

from planetgen.physics import constants

AU_KM = constants.AU_TO_KM
SOLAR_MASS_KG = constants.SOLAR_MASS_TO_KG

MAX_FILL_FACTOR = 1.2
"""float: How far past its Roche lobe a star may reach at the pair's closest approach (GEN.1006 4.5; the generator's own cap)."""

MUTUAL_HILL_RADII_SEPARATION = 10.0
"""float: Neighbouring planets of one host are at least this many mutual Hill radii apart (Chambers, Wetherill and Boss 1996)."""

MOON_PROGRADE_HILL_FRACTION = 0.4895
"""float: A prograde moon is stable inside this fraction of its planet's Hill radius (Domingos, Winter and Yokoyama 2006)."""

RIGID_ROCHE_COEFFICIENT = 1.26
"""float: A rigid satellite breaks up inside `1.26 R (rho_host / rho_body)^(1/3)` of its host (Roche 1849, in Murray and Dermott 1999);
the fluid coefficient 2.44 would flag a rocky body the generator is entitled to place."""

LIMIT_TOLERANCE = 1.0e-4
"""float: How far past a limit a value may stand before it is a violation: the generator places bodies exactly on a limit, and the
rows hold single-precision numbers."""

MASS_TOLERANCE = 1.0e-3
"""float: How closely a node's stored mass equals its children's (the star rows store kilograms, the nodes solar masses)."""

NO_FILL_CHECK = ("merger-product", "post-common-envelope")
"""tuple: The `orbit_nodes.interaction` values whose pair was shaped by a past interaction, so the detached Roche rule does not hold."""

HW_S_MU_RANGE = (0.1, 0.9)
"""tuple: The mass ratios Holman and Wiegert 1999 fitted their S-type critical axis on (the companion's share of the pair)."""

HW_S_E_MAX = 0.8
"""float: The largest eccentricity their S-type fit covers."""

HW_P_MU_RANGE = (0.1, 0.5)
"""tuple: The mass ratios they fitted the P-type critical axis on (the lighter star's share of the pair)."""

HW_P_E_MAX = 0.7
"""float: The largest eccentricity their P-type fit covers."""


@dataclass(frozen=True)
class Violation:
    """One overlap or missing room: what `kind` of check failed, where (`where`: the rows involved), the measured value and the
    limit it broke, in the units named in `unit`, and whether the limit was outside the range its fit covers."""
    kind: str
    where: str
    measured: float
    limit: float
    unit: str = ""
    extrapolated: bool = False
    detail: str = field(default="", compare=False)
    hard: bool = True

    def __str__(self):
        edge = " (extrapolated)" if self.extrapolated else ""
        edge += "" if self.hard else " (advisory)"
        return f"{self.kind}{edge} at {self.where}: {self.measured:.6g} against {self.limit:.6g} {self.unit}{(' ' + self.detail) if self.detail else ''}"


# --- the forms, each from the paper -----------------------------------------------------------------------------------


def roche_lobe_fraction(mass_self, mass_other):
    """Eggleton 1983: the Roche-lobe radius of the first body over the separation, `0.49 q^(2/3) / (0.6 q^(2/3) + ln(1 + q^(1/3)))`."""
    q = mass_self / mass_other
    return 0.49 * q ** (2.0 / 3.0) / (0.6 * q ** (2.0 / 3.0) + math.log(1.0 + q ** (1.0 / 3.0)))


def mardling_aarseth_axis_ratio(q_out, e_out, inclination):
    """Mardling and Aarseth 2001: the smallest stable `a_out / a_in`, `q_out = m3 / (m1 + m2)`, mutual `inclination` in radians."""
    peri = 2.8 * ((1.0 + q_out) * (1.0 + e_out) / math.sqrt(1.0 - e_out)) ** 0.4 * (1.0 - 0.3 * inclination / math.pi)
    return peri / (1.0 - e_out)


def eggleton_kiseleva_axis_ratio(q_out, q_in, e_in, e_out):
    """Eggleton and Kiseleva 1995: the smallest stable `a_out / a_in`; `q_in = m2 / m1` of the inner pair."""
    q13 = q_out ** (1.0 / 3.0)
    y_crit = 1.0 + 3.7 / q13 - 2.2 / (1.0 + q13) + 1.4 / q_in ** (1.0 / 3.0) * (q13 - 1.0) / (q13 + 1.0)
    return y_crit * (1.0 + e_in) / (1.0 - e_out)


def vynatheya_axis_ratio(q_out, e_in, e_out, inclination):
    """Vynatheya, Hamers and Mardling 2022, equation 4 with the effective inner eccentricity of equation 3."""
    window = 1.0 - 5.0 / 3.0 * math.cos(inclination) ** 2
    e_max = math.sqrt(window) if window > 0.0 else 0.0
    e_eff = max(e_in, 0.5 * e_max ** 2)
    base = 2.4 * ((1.0 + q_out) / ((1.0 + e_eff) * math.sqrt(1.0 - e_out))) ** 0.4
    y_crit = base * (((1.0 - 0.2 * e_eff + e_out) / 8.0) * (math.cos(inclination) - 1.0) + 1.0)
    return y_crit * (1.0 + e_eff) / (1.0 - e_out)


def holman_wiegert_s_axis_ratio(mu, e):
    """Holman and Wiegert 1999, equation 1: the largest stable `a / a_b` of a test particle about the heavier star (`mu = m2 / (m1 + m2)`)."""
    return 0.464 - 0.380 * mu - 0.631 * e + 0.586 * mu * e + 0.150 * e ** 2 - 0.198 * mu * e ** 2


def holman_wiegert_p_axis_ratio(mu, e):
    """Holman and Wiegert 1999, equation 3: the smallest stable `a / a_b` of a test particle about the pair."""
    return 1.60 + 5.10 * e - 2.22 * e ** 2 + 4.12 * mu - 4.27 * e * mu - 5.09 * mu ** 2 + 4.61 * e ** 2 * mu ** 2


def _mutual_inclination(inc1, node1, inc2, node2):
    cos_i = math.cos(inc1) * math.cos(inc2) + math.sin(inc1) * math.sin(inc2) * math.cos(node1 - node2)
    return math.acos(max(-1.0, min(1.0, cos_i)))


def _clamped_fit_inputs(mu, e, mu_range=HW_S_MU_RANGE, e_max=HW_S_E_MAX):
    """`mu` and `e` held inside the Holman and Wiegert fit, and whether they had to be."""
    clamped_mu = min(max(mu, mu_range[0]), mu_range[1])
    clamped_e = min(max(e, 0.0), e_max)
    return clamped_mu, clamped_e, (clamped_mu != mu or clamped_e != e)


# --- the rows ---------------------------------------------------------------------------------------------------------


class _Tree:
    """The stars and nodes of one system as a tree: `children[node_no] = {slot: ("star", row) | ("node", row)}`."""

    def __init__(self, stars, nodes):
        self.stars = {row["id"]: row for row in stars}
        self.nodes = {row["node_no"]: row for row in nodes}
        self.children = {}
        for star in stars:
            if star.get("parent_node_no") is not None:
                self.children.setdefault(star["parent_node_no"], {})[star.get("slot") or 0] = ("star", star)
        for node in nodes:
            if node["parent_node_no"] is not None:
                self.children.setdefault(node["parent_node_no"], {})[node["slot"]] = ("node", node)
        self._mass = {}

    def kids(self, node_no):
        slots = self.children.get(node_no, {})
        return slots.get(0), slots.get(1)

    def mass_solar(self, kind_row):
        kind, row = kind_row
        if kind == "star":
            return row["mass_kg"] / SOLAR_MASS_KG
        return row["mass_solar"]

    def stars_below(self, node_no):
        out = []
        for kind, row in self.children.get(node_no, {}).values():
            out.extend([row] if kind == "star" else self.stars_below(row["node_no"]))
        return out

    def extent_au(self, kind_row):
        """The greatest distance from a side's own barycentre to the surface of any star below it, at the apocentres, AU."""
        kind, row = kind_row
        if kind == "star":
            return row["radius_km"] / AU_KM
        low, high = self.kids(row["node_no"])
        if low is None or high is None:
            return 0.0
        mass_low, mass_high = self.mass_solar(low), self.mass_solar(high)
        total = mass_low + mass_high
        apocentre = row["separation_km"] / AU_KM * (1.0 + row["eccentricity"])
        return max(apocentre * mass_high / total + self.extent_au(low), apocentre * mass_low / total + self.extent_au(high))

    def ancestors(self, star_id):
        """The nodes above a star, nearest first, each with the mass on the star's side and on the other side."""
        out = []
        star = self.stars[star_id]
        node_no, side = star.get("parent_node_no"), ("star", star)
        while node_no is not None:
            node = self.nodes[node_no]
            low, high = self.kids(node_no)
            mine, other = (low, high) if low is not None and low[1] is side[1] else (high, low)
            if mine is None or other is None:
                break
            out.append((node, self.mass_solar(mine), self.mass_solar(other)))
            side = ("node", node)
            node_no = node["parent_node_no"]
        return out


def _audit_nodes(tree):
    out = []
    for node in tree.nodes.values():
        where = f"node {node['node_no']}"
        low, high = tree.kids(node["node_no"])
        if low is None or high is None:
            out.append(Violation("mass", where, 0.0, 2.0, "children", detail="a node needs a child in slot 0 and in slot 1"))
            continue
        total = tree.mass_solar(low) + tree.mass_solar(high)
        if abs(total - node["mass_solar"]) > MASS_TOLERANCE * node["mass_solar"]:
            out.append(Violation("mass", where, node["mass_solar"], total, "Msun", detail="the node's mass is not its children's"))
        axis = node["separation_km"] / AU_KM
        peri = axis * (1.0 - node["eccentricity"])
        single_pair = low[0] == "star" and high[0] == "star"
        if single_pair:
            if node.get("interaction") in NO_FILL_CHECK:
                continue
            for mine, other in ((low, high), (high, low)):
                lobe = roche_lobe_fraction(tree.mass_solar(mine), tree.mass_solar(other)) * peri
                radius = mine[1]["radius_km"] / AU_KM
                if radius > MAX_FILL_FACTOR * lobe * (1.0 + LIMIT_TOLERANCE):
                    out.append(Violation("roche", f"{where}, star {mine[1]['id']}", radius / lobe, MAX_FILL_FACTOR, "of its lobe",
                                         detail=f"radius {radius:.5g} AU, lobe {lobe:.5g} AU at pericentre {peri:.5g} AU"))
            continue
        extent = tree.extent_au(low) + tree.extent_au(high)
        if peri < extent * (1.0 - LIMIT_TOLERANCE):
            out.append(Violation("overlap", where, peri, extent, "AU", detail="the pericentre is inside the sum of the two sides' extents"))
    return out


def _audit_stability(tree):
    out = []
    for outer in tree.nodes.values():
        low, high = tree.kids(outer["node_no"])
        if low is None or high is None:
            continue
        for inner_side, third_side in ((low, high), (high, low)):
            if inner_side[0] != "node":
                continue
            inner = inner_side[1]
            first, second = tree.kids(inner["node_no"])
            if first is None or second is None:
                continue
            m1, m2 = sorted((tree.mass_solar(first), tree.mass_solar(second)), reverse=True)
            m3 = tree.mass_solar(third_side)
            a_in, a_out = inner["separation_km"] / AU_KM, outer["separation_km"] / AU_KM
            e_in, e_out = inner["eccentricity"], outer["eccentricity"]
            inclination = _mutual_inclination(math.radians(inner["inclination_deg"]), math.radians(inner["ascending_node_deg"]),
                                              math.radians(outer["inclination_deg"]), math.radians(outer["ascending_node_deg"]))
            q_out, q_in = m3 / (m1 + m2), m2 / m1
            required = (mardling_aarseth_axis_ratio(q_out, e_out, inclination),
                        eggleton_kiseleva_axis_ratio(q_out, q_in, e_in, e_out),
                        vynatheya_axis_ratio(q_out, e_in, e_out, inclination))
            ratio = a_out / a_in
            if ratio < min(required) * (1.0 - LIMIT_TOLERANCE):
                out.append(Violation("stability", f"node {outer['node_no']} about node {inner['node_no']}", ratio, min(required),
                                     "a_out / a_in", detail="below Mardling and Aarseth, Eggleton and Kiseleva and Vynatheya et al."))
    return out


# --- planets, moons and belts -----------------------------------------------------------------------------------------


def _density(row):
    volume = 4.0 / 3.0 * math.pi * (row["radius_km"] * 1.0e3) ** 3
    return row["mass_kg"] / volume if volume > 0.0 and row.get("mass_kg") else None


def _host_of(row, tree):
    """What a planet goes round: `("star", star row)`, `("pair", node row)` (the pair's barycentre), or None."""
    star_id = row.get("star_id")
    if star_id is not None and star_id in tree.stars:
        return ("star", tree.stars[star_id])
    pairs = [node for node in tree.nodes.values()
             if all(kind == "star" and not row.get("dark") for kind, row in tree.children.get(node["node_no"], {}).values())
             and len(tree.children.get(node["node_no"], {})) == 2]
    lit = [row for row in tree.stars.values() if not row.get("dark")]
    if len(lit) == 1:
        return ("star", lit[0])
    if len(pairs) == 1:
        return ("pair", pairs[0])
    roots = [node for node in pairs if node["parent_node_no"] is None]
    return ("pair", roots[0]) if len(roots) == 1 else None


def _host_mass_solar(host, tree):
    kind, row = host
    if kind == "star":
        return row["mass_kg"] / SOLAR_MASS_KG
    return row["mass_solar"]


def _host_radius_km(host, tree):
    kind, row = host
    if kind == "star":
        return row["radius_km"]
    return max(star["radius_km"] for star in tree.stars_below(row["node_no"]))


def _s_limits(host, tree):
    """For a planet about a star, `(a_critical_au, extrapolated, where)` toward each companion above it."""
    out = []
    kind, row = host
    if kind != "star":
        return out
    for node, mine, other in tree.ancestors(row["id"]):
        mu = other / (mine + other)
        use_mu, use_e, extrapolated = _clamped_fit_inputs(mu, node["eccentricity"])
        a_b = node["separation_km"] / AU_KM
        out.append((holman_wiegert_s_axis_ratio(use_mu, use_e) * a_b, extrapolated, f"node {node['node_no']}"))
    return out


def _audit_planets(planets, tree):
    out = []
    by_host = {}
    for planet in planets:
        host = _host_of(planet, tree)
        where = f"planet {planet['id']}"
        if host is None:
            continue
        a_au = planet["distance_km"] / AU_KM
        e = planet.get("eccentricity") or 0.0
        peri = a_au * (1.0 - e)
        host_key = (host[0], host[1]["id"] if host[0] == "star" else host[1]["node_no"])
        by_host.setdefault(host_key, []).append((planet, a_au, e, host))
        host_radius = _host_radius_km(host, tree) / AU_KM
        if peri < host_radius:
            out.append(Violation("planet-in-star", where, peri, host_radius, "AU", detail="closest approach inside the host star"))
        rho_host = _density(host[1]) if host[0] == "star" else None
        rho_planet = _density(planet)
        if rho_host and rho_planet:
            roche = RIGID_ROCHE_COEFFICIENT * host_radius * (rho_host / rho_planet) ** (1.0 / 3.0)
            if peri < roche:
                out.append(Violation("planet-roche", where, peri, roche, "AU", hard=False, detail="closest approach inside the rigid-body Roche limit"))
        if host[0] == "star":
            for limit, extrapolated, node_where in _s_limits(host, tree):
                if a_au > limit * (1.0 + LIMIT_TOLERANCE):
                    out.append(Violation("s-type", f"{where} about star {host[1]['id']} toward {node_where}", a_au, limit, "AU", extrapolated))
        else:
            node = host[1]
            low, high = tree.kids(node["node_no"])
            mu_raw = min(tree.mass_solar(low), tree.mass_solar(high)) / node["mass_solar"]
            use_mu, use_e, extrapolated = _clamped_fit_inputs(mu_raw, node["eccentricity"], HW_P_MU_RANGE, HW_P_E_MAX)
            limit = holman_wiegert_p_axis_ratio(use_mu, use_e) * node["separation_km"] / AU_KM
            if a_au < limit * (1.0 - LIMIT_TOLERANCE):
                out.append(Violation("p-type", f"{where} about pair node {node['node_no']}", a_au, limit, "AU", extrapolated))
            for upper, extrapolated_up, node_where in _p_outer_limits(node, tree):
                if a_au > upper * (1.0 + LIMIT_TOLERANCE):
                    out.append(Violation("s-type", f"{where} about pair node {node['node_no']} toward {node_where}", a_au, upper, "AU", extrapolated_up))
    for members in by_host.values():
        members.sort(key=lambda item: item[1])
        for (inner, a_in, e_in, host), (outer, a_out, e_out, _host) in zip(members, members[1:]):
            where = f"planets {inner['id']} and {outer['id']}"
            if a_in * (1.0 + e_in) >= a_out * (1.0 - e_out):
                out.append(Violation("planet-crossing", where, a_in * (1.0 + e_in), a_out * (1.0 - e_out), "AU"))
            total = _host_mass_solar(host, tree) * SOLAR_MASS_KG
            if total > 0.0 and inner.get("mass_kg") and outer.get("mass_kg"):
                hill = ((inner["mass_kg"] + outer["mass_kg"]) / (3.0 * total)) ** (1.0 / 3.0) * 0.5 * (a_in + a_out)
                if a_out - a_in < MUTUAL_HILL_RADII_SEPARATION * hill * (1.0 - LIMIT_TOLERANCE):
                    out.append(Violation("planet-spacing", where, (a_out - a_in) / hill, MUTUAL_HILL_RADII_SEPARATION, "mutual Hill radii"))
    return out


def _p_outer_limits(node, tree):
    """The S-type critical axis of a pair (as one body) toward each companion above its node, for a planet about the pair."""
    out = []
    parent_no, child = node["parent_node_no"], node
    while parent_no is not None:
        parent = tree.nodes[parent_no]
        low, high = tree.kids(parent_no)
        mine, other = (low, high) if low[1] is child else (high, low)
        mass_mine, mass_other = tree.mass_solar(mine), tree.mass_solar(other)
        mu = mass_other / (mass_mine + mass_other)
        use_mu, use_e, extrapolated = _clamped_fit_inputs(mu, parent["eccentricity"])
        out.append((holman_wiegert_s_axis_ratio(use_mu, use_e) * parent["separation_km"] / AU_KM, extrapolated, f"node {parent_no}"))
        child, parent_no = parent, parent["parent_node_no"]
    return out


def _audit_moons(moons, planets, tree):
    out = []
    by_id = {planet["id"]: planet for planet in planets}
    for moon in moons:
        planet = by_id.get(moon["planet_id"])
        if planet is None or not planet.get("mass_kg"):
            continue
        host = _host_of(planet, tree)
        if host is None:
            continue
        where = f"moon {moon['id']}"
        host_mass = _host_mass_solar(host, tree) * SOLAR_MASS_KG
        hill_au = planet["distance_km"] / AU_KM * (planet["mass_kg"] / (3.0 * host_mass)) ** (1.0 / 3.0)
        a_moon = moon["distance_km"] / AU_KM
        e_moon = moon.get("eccentricity") or 0.0
        if a_moon > MOON_PROGRADE_HILL_FRACTION * hill_au * (1.0 + LIMIT_TOLERANCE):
            out.append(Violation("moon-hill", where, a_moon, MOON_PROGRADE_HILL_FRACTION * hill_au, "AU"))
        touching = (planet["radius_km"] + moon["radius_km"]) / AU_KM
        if a_moon * (1.0 - e_moon) < touching:
            out.append(Violation("moon-in-planet", where, a_moon * (1.0 - e_moon), touching, "AU", detail="closest approach inside the planet and moon radii"))
        rho_planet, rho_moon = _density(planet), _density(moon)
        if rho_planet and rho_moon:
            roche = RIGID_ROCHE_COEFFICIENT * planet["radius_km"] / AU_KM * (rho_planet / rho_moon) ** (1.0 / 3.0)
            if a_moon * (1.0 - e_moon) < roche:
                out.append(Violation("moon-roche", where, a_moon * (1.0 - e_moon), roche, "AU", hard=False))
    return out


def _audit_belts(belts, planets, tree):
    out = []
    for belt in belts:
        where = f"belt {belt['id']}"
        host = _host_of(belt, tree)
        lower, upper = belt["lower_limit_km"] / AU_KM, belt["upper_limit_km"] / AU_KM
        if lower >= upper:
            out.append(Violation("belt-planet", where, lower, upper, "AU", detail="an empty annulus"))
            continue
        for planet in planets:
            if _host_of(planet, tree) != host:
                continue
            a_au = planet["distance_km"] / AU_KM
            e = planet.get("eccentricity") or 0.0
            if a_au * (1.0 + e) > lower and a_au * (1.0 - e) < upper:
                out.append(Violation("belt-planet", f"{where} and planet {planet['id']}", a_au, upper, "AU", detail="the planet's orbit reaches into the belt"))
        if host is not None:
            for limit, extrapolated, node_where in _s_limits(host, tree):
                if upper > limit * (1.0 + LIMIT_TOLERANCE):
                    out.append(Violation("belt-s-type", f"{where} toward {node_where}", upper, limit, "AU", extrapolated))
    return out


def hard_violations(rows):
    """The hard violations of `audit(rows)`: bodies that overlap, or a limit the generator promises that is broken."""
    return [item for item in audit(rows) if item.hard]


def audit(rows):
    """
    Every violation in one finished system, from what it stores, advisories (`Violation.hard` false) included.

    Args:
        rows (dict): `stars` (`id`, `mass_kg`, `radius_km`, `parent_node_no`, `slot`; a dark member of the tree, GEN.1013, is a row
            too, with `dark` true: it weighs in the tree and takes up room, and no planet goes round it), `nodes` (the `orbit_nodes` columns),
            `planets` (`id`, `star_id`, `distance_km`, `eccentricity`, `radius_km`, `mass_kg`), `moons` (`id`, `planet_id`,
            `distance_km`, `eccentricity`, `radius_km`, `mass_kg`) and `belts` (`id`, `star_id`, `lower_limit_km`,
            `upper_limit_km`); the last three may be empty. A planet or belt with no `star_id` goes round the system's pair.

    Returns:
        list[Violation]: Empty for a system with room for every body.
    """
    tree = _Tree(rows["stars"], rows.get("nodes") or [])
    planets, moons, belts = rows.get("planets") or [], rows.get("moons") or [], rows.get("belts") or []
    found = _audit_nodes(tree) + _audit_stability(tree)
    found += _audit_planets(planets, tree) + _audit_moons(moons, planets, tree) + _audit_belts(belts, planets, tree)
    return found


def stored_rows(conn, system_id):
    """The rows `audit` reads for one stored system, read with plain SELECTs: the star rows' `parent_node_id` is the node number
    they hang from (`parent_node_no`), as `systemscene` reads it. The system's `dark_members` come after its stars, with negative ids
    and `dark` true."""
    stars = [dict(row, parent_node_no=row["parent_node_id"]) for row in conn.execute(
        "SELECT id, mass_kg, radius_km, parent_node_id, slot FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()]
    stars += [{"id": -1 - row["member_no"], "mass_kg": row["mass_solar"] * SOLAR_MASS_KG, "radius_km": (row["radius_au"] or 0.0) * AU_KM,
               "parent_node_no": row["parent_node_no"], "slot": row["slot"], "dark": True, "kind": row["kind"]}
              for row in conn.execute("SELECT member_no, parent_node_no, slot, kind, mass_solar, radius_au FROM dark_members"
                                      " WHERE star_system_id = ? ORDER BY member_no", (system_id,)).fetchall()]
    nodes = [dict(row) for row in conn.execute(
        "SELECT node_no, parent_node_no, slot, mass_solar, separation_km, eccentricity, inclination_deg, ascending_node_deg,"
        " arg_periapsis_deg, mean_anomaly_deg, period_years, kind, interaction, flags FROM orbit_nodes WHERE star_system_id = ?"
        " ORDER BY node_no", (system_id,)).fetchall()]
    planets = [dict(row) for row in conn.execute(
        "SELECT id, star_id, distance_km, eccentricity, radius_km, mass_kg FROM planets WHERE star_system_id = ? ORDER BY id",
        (system_id,)).fetchall()]
    moons = [dict(row) for row in conn.execute(
        "SELECT id, planet_id, distance_km, eccentricity, radius_km, mass_kg FROM moons WHERE star_system_id = ? ORDER BY id",
        (system_id,)).fetchall()]
    belts = [dict(row) for row in conn.execute(
        "SELECT id, star_id, lower_limit_km, upper_limit_km FROM asteroid_belts WHERE star_system_id = ? ORDER BY id",
        (system_id,)).fetchall()]
    return {"stars": stars, "nodes": nodes, "planets": planets, "moons": moons, "belts": belts}
