"""
GEN.1014: `planetgen.physics.clearance_ledger`, the ledger of where a planet, a moon or a belt may sit (GEN.1006 section 4.6). Four groups:
the span arithmetic and the truncated draw; each exclusion held against the generator's own closed form (so the ledger keeps the rule
the placer keeps today) and against the independent verifier's reading of the published one; the windows of a tree of stars; and a sweep
that places bodies only where the ledger says they may go and requires the verifier (`physics/clearance.py`) to find no hard violation.
"""

import math
import random
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation import validation
from planetgen.physics import clearance, constants, orbits, planet_zones
from planetgen.physics import clearance_ledger as ledger
from planetgen.physics.system_tree import OrbitNode, TreeMember

AU_KM = constants.AU_TO_KM
SOLAR = constants.SOLAR_MASS_TO_KG
EARTH = constants.EARTH_MASS_TO_KG


# --- spans --------------------------------------------------------------------------------------------------------------------------------

def test_spans_are_sorted_joined_and_empty_ones_dropped():
    spans = ledger.SpanSet([(5, 6), (1, 2), (1.5, 3), (7, 7), (9, 8)])
    assert [(s.low, s.high) for s in spans] == [(1, 3), (5, 6)]
    assert spans.width == pytest.approx(3.0) and not ledger.SpanSet()


def test_subtract_intersect_contains_and_nearest():
    spans = ledger.SpanSet([(1, 10)])
    cut = spans.subtract(2, 4).subtract(6, 7)
    assert [(s.low, s.high) for s in cut] == [(1, 2), (4, 6), (7, 10)]
    assert cut.contains(2) and not cut.contains(3) and cut.covers(4.5, 5.5) and not cut.covers(1.5, 4.5)
    assert [(s.low, s.high) for s in cut.intersect(5, 8)] == [(5, 6), (7, 8)]
    assert cut.nearest(3.0) in (2, 4) and cut.nearest(0.0) == 1 and cut.nearest(99.0) == 10 and ledger.SpanSet().nearest(1.0) is None
    assert spans.subtract(0, 20) == ledger.SpanSet()


def test_log_uniform_draw_stays_in_the_spans_and_keeps_their_share_of_the_law():
    spans = ledger.SpanSet([(1.0, 2.0), (8.0, 16.0)])
    rng = random.Random(1014)
    draws = [ledger.log_uniform_quantile(rng.random(), spans) for _ in range(20000)]
    assert all(spans.contains(x) for x in draws)
    assert ledger.log_uniform_quantile(0.0, spans) == pytest.approx(1.0) and ledger.log_uniform_quantile(1.0, spans) == pytest.approx(16.0)
    assert sum(x < 4.0 for x in draws) / len(draws) == pytest.approx(0.5, abs=0.02)   # equal log widths


def test_a_restricted_quantile_is_the_conditional_law_for_any_law():
    # an exponential law cut to [0, 1] and [3, 4]
    spans = ledger.SpanSet([(0.0, 1.0), (3.0, 4.0)])
    cdf, ppf = (lambda x: 1.0 - math.exp(-x)), (lambda p: -math.log(1.0 - p))
    low_mass, high_mass = cdf(1.0) - cdf(0.0), cdf(4.0) - cdf(3.0)
    rng = random.Random(3)
    draws = [ledger.restricted_quantile(rng.random(), spans, cdf, ppf) for _ in range(20000)]
    assert sum(x <= 1.0 for x in draws) / len(draws) == pytest.approx(low_mass / (low_mass + high_mass), abs=0.02)
    assert all(spans.contains(x) for x in draws)


def test_a_draw_with_no_probability_in_the_spans_says_so():
    assert ledger.restricted_quantile(0.5, ledger.SpanSet(), math.log, math.exp) is None
    assert ledger.log_uniform_quantile(0.5, ledger.SpanSet([(-3.0, -1.0)])) is None


# --- exclusions ---------------------------------------------------------------------------------------------------------------------------

def _stub(distance, mass, eccentricity, star_mass):
    return SimpleNamespace(distance=distance, mass=mass, eccentricity=eccentricity, star=SimpleNamespace(mass=star_mass))


@pytest.mark.parametrize("e_in,e_out,m_in,m_out", [(0.0, 0.0, 6e24, 6e24), (0.1, 0.05, 1.9e27, 3e25), (0.3, 0.4, 5e26, 5e26), (0.02, 0.6, 1e23, 2e28)])
def test_the_outer_edge_of_a_planet_exclusion_is_the_generators_closed_form(e_in, e_out, m_in, m_out):
    star = 1.0 * SOLAR
    _inner, outer = ledger.planet_exclusion(1.3, e_in, m_in, star, m_out, e_out)
    assert outer == pytest.approx(validation.mutual_min_distance_au(_stub(2.0, m_out, e_out, star), _stub(1.3, m_in, e_in, star)), rel=1e-12)


@pytest.mark.parametrize("e_in,e_out", [(0.0, 0.0), (0.1, 0.05), (0.3, 0.4)])
def test_both_edges_of_a_planet_exclusion_stand_ten_mutual_hill_radii_from_the_planet(e_in, e_out):
    star, m_a, m_b, d = 1.0 * SOLAR, 3e25, 8e26, 1.7
    inner, outer = ledger.planet_exclusion(d, e_in, m_a, star, m_b, e_out)
    assert 0.0 < inner < d < outer
    for x, (a_in, e_i, a_out, e_o) in ((outer, (d, e_in, outer, e_out)), (inner, (inner, e_out, d, e_in))):
        hill = orbits.mutual_hill_radius_au(m_a, m_b, a_in, a_out, star)
        gap = a_out * (1.0 - e_o) - a_in * (1.0 + e_i)
        assert gap / hill == pytest.approx(tuning.MUTUAL_HILL_RADII_SEPARATION, rel=1e-9)


def test_a_planet_too_heavy_to_leave_room_inside_leaves_none():
    inner, outer = ledger.planet_exclusion(1.0, 0.9, 1e28, 1.0 * SOLAR, 1e28, 0.0)
    assert inner == 0.0 and outer > 1.0


@pytest.mark.parametrize("mass,e", [(6e24, 0.0), (1.9e27, 0.2), (1e23, 0.5)])
def test_the_outer_edge_of_a_belt_exclusion_is_the_generators_closed_form(mass, e):
    star = 1.0 * SOLAR
    planet = _stub(0.0, mass, e, star)
    belt = SimpleNamespace(upper_limit=2.3)
    # the generator reads the planet's own Hill radius at the distance it sits at; the closed form divides by that distance, so any value agrees
    planet.distance = 4.0
    planet.hill_radius = orbits.calculate_hill_sphere(4.0 * constants.AU_TO_M, mass, star) / 1000.0
    _inner, outer = ledger.belt_exclusion(1.9, 2.3, mass, star, e)
    assert outer == pytest.approx(validation.min_distance_past_belt_au(planet, belt), rel=1e-9)


def test_a_planet_inside_a_belt_exclusion_keeps_five_hill_radii_from_the_belts_inner_edge():
    star, mass = 1.0 * SOLAR, 1.9e27
    inner, _outer = ledger.belt_exclusion(3.0, 3.4, mass, star)
    hill = inner * (mass / (3.0 * (star + mass))) ** (1.0 / 3.0)
    assert inner + 5.0 * hill == pytest.approx(3.0, rel=1e-9)


# --- moons --------------------------------------------------------------------------------------------------------------------------------

def test_the_domingos_limits_at_zero_eccentricity_are_the_papers_and_the_generators():
    assert ledger.moon_hill_fraction(0.0)[0] == pytest.approx(tuning.MOON_PROGRADE_STABLE_HILL_FRACTION) == pytest.approx(0.4895)
    assert ledger.moon_hill_fraction(0.0, retrograde=True)[0] == pytest.approx(0.9309)


def test_the_domingos_limits_fall_with_each_eccentricity_and_never_go_negative():
    assert ledger.moon_hill_fraction(0.0, 0.0)[0] > ledger.moon_hill_fraction(0.3, 0.0)[0] > ledger.moon_hill_fraction(0.3, 0.2)[0] > 0.0
    assert ledger.moon_hill_fraction(0.3, 0.0)[0] == pytest.approx(0.4895 * (1 - 1.0305 * 0.3))
    assert ledger.moon_hill_fraction(0.3, 0.2, retrograde=True)[0] == pytest.approx(0.9309 * (1 - 1.0764 * 0.3 - 0.9812 * 0.2))
    assert ledger.moon_hill_fraction(0.9, 0.9)[0] == 0.0
    assert ledger.moon_hill_fraction(0.9)[1] is False and ledger.moon_hill_fraction(0.95)[1] is True
    assert ledger.moon_hill_fraction(0.95)[0] == ledger.moon_hill_fraction(0.9)[0]


def test_a_moon_window_runs_from_the_planet_and_its_largest_moon_to_the_hill_limit():
    window, extrapolated = ledger.moon_window_km(6371.0, 1.5e6, 0.0, margin_km=300.0, largest_moon_radius_km=3000.0)
    assert not extrapolated and [(s.low, s.high) for s in window] == [(pytest.approx(9671.0), pytest.approx(0.4895 * 1.5e6))]
    assert ledger.moon_window_km(6371.0, 1.5e6, 0.0, retrograde=True)[0].width > window.width
    assert not ledger.moon_window_km(6371.0, 1.0e4, 0.0)[0]


# --- the windows of a tree ----------------------------------------------------------------------------------------------------------------

def _node(node_no, parent, slot, mass, axis_au, e=0.0, kind="close"):
    return OrbitNode(node_no, parent, slot, mass, axis_au * AU_KM, e, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, kind)


def _tree_of(kind):
    if kind == "single":
        return [], [TreeMember(0, None, 0, 1.0)]
    if kind == "wide":
        return [_node(0, None, 0, 2.0, 50.0, 0.2)], [TreeMember(0, 0, 0, 1.0), TreeMember(1, 0, 1, 1.0)]
    if kind == "close":
        return [_node(0, None, 0, 1.4, 0.1, 0.0)], [TreeMember(0, 0, 0, 0.7), TreeMember(1, 0, 1, 0.7)]
    if kind == "triple":      # a close pair with a wide third star
        return ([_node(1, None, 0, 2.0, 80.0, 0.3), _node(0, 1, 0, 1.0, 0.2, 0.1)],
                [TreeMember(0, 0, 0, 0.5), TreeMember(1, 0, 1, 0.5), TreeMember(2, 1, 1, 1.0)])
    if kind == "nested":      # a wide pair whose primary has a close companion that has a still closer one
        return ([_node(2, None, 0, 3.0, 400.0, 0.2), _node(1, 2, 0, 2.0, 12.0, 0.1), _node(0, 1, 0, 1.0, 0.05, 0.0)],
                [TreeMember(0, 0, 0, 0.5), TreeMember(1, 0, 1, 0.5), TreeMember(2, 1, 1, 1.0), TreeMember(3, 2, 1, 1.0)])
    raise ValueError(kind)


def test_a_lone_star_has_one_window_from_its_floor_to_the_disk_cap():
    nodes, members = _tree_of("single")
    book = ledger.Ledger(nodes, members, floors_au={0: 0.7})
    assert book.hosts() == [(planet_zones.S_TYPE, 0)]
    assert [(s.low, s.high) for s in book.window((planet_zones.S_TYPE, 0))] == [(0.7, tuning.PLANET_ZONE_DISK_CAP_AU)]
    assert book.central_mass_kg((planet_zones.S_TYPE, 0)) == pytest.approx(SOLAR)
    assert not book.extrapolated((planet_zones.S_TYPE, 0))


def test_a_close_pair_has_a_circumbinary_window_that_starts_outside_the_pair():
    nodes, members = _tree_of("close")
    book = ledger.Ledger(nodes, members)
    host = (planet_zones.P_TYPE, 0)
    assert host in book.hosts() and book.window(host).spans[0].low == pytest.approx(orbits.holman_wiegert_circumbinary_a_crit_au(0.1, 0.5, 0.0))
    assert book.central_mass_kg(host) == pytest.approx(1.4 * SOLAR)


def test_a_ceiling_the_caller_gives_cuts_the_window():
    nodes, members = _tree_of("single")
    book = ledger.Ledger(nodes, members, ceilings_au={(planet_zones.S_TYPE, 0): 35.0})
    assert book.window((planet_zones.S_TYPE, 0)).spans[0].high == 35.0


def test_a_stars_window_is_held_by_every_level_above_it():
    nodes, members = _tree_of("nested")
    book = ledger.Ledger(nodes, members)
    zones = {(z.kind, z.host_no): z for z in planet_zones.planet_zones(nodes, members)}
    limits = planet_zones.ancestor_limits_au(zones[(planet_zones.S_TYPE, 0)], nodes, members)
    assert len(limits) == 3 and book.window((planet_zones.S_TYPE, 0)).spans[0].high == pytest.approx(min(limits + [100.0]))
    assert book.window((planet_zones.S_TYPE, 0)).spans[0].high <= zones[(planet_zones.S_TYPE, 0)].outer_au


def test_a_window_outside_the_fit_range_says_so():
    nodes, members = _tree_of("wide")
    assert not ledger.Ledger(nodes, members).extrapolated((planet_zones.S_TYPE, 0))
    nodes, members = [_node(0, None, 0, 1.05, 10.0, 0.2)], [TreeMember(0, 0, 0, 1.0), TreeMember(1, 0, 1, 0.05)]
    assert ledger.Ledger(nodes, members).extrapolated((planet_zones.S_TYPE, 0))     # a mass fraction of 0.05
    nodes, members = [_node(0, None, 0, 2.0, 10.0, 0.85)], [TreeMember(0, 0, 0, 1.0), TreeMember(1, 0, 1, 1.0)]
    assert ledger.Ledger(nodes, members).extrapolated((planet_zones.S_TYPE, 0))    # an eccentricity of 0.85
    nodes, members = [_node(0, None, 0, 1.4, 0.1, 0.8)], [TreeMember(0, 0, 0, 0.7), TreeMember(1, 0, 1, 0.7)]
    assert ledger.Ledger(nodes, members).extrapolated((planet_zones.P_TYPE, 0))


# --- room for planets and belts -----------------------------------------------------------------------------------------------------------

def test_a_planet_leaves_no_room_close_to_itself_and_a_belt_none_inside_it():
    nodes, members = _tree_of("single")
    host = (planet_zones.S_TYPE, 0)
    book = ledger.Ledger(nodes, members, floors_au={0: 0.1})
    open_before = book.open_spans(host, 6e24)
    book.place_planet(host, 1.0, 6e24)
    after = book.open_spans(host, 6e24)
    assert open_before.width > after.width and not after.contains(1.0) and not after.contains(1.1) and after.contains(0.5) and after.contains(2.0)
    book.place_belt(host, 3.0, 3.4)
    withbelt = book.open_spans(host, 6e24)
    assert not withbelt.contains(3.2) and not withbelt.contains(3.45) and withbelt.contains(2.0) and withbelt.contains(6.0)
    belts = book.open_belt_spans(host, 0.2)
    assert not belts.contains(1.0) and not belts.contains(3.2) and belts.contains(2.0) and belts.contains(10.0)


def test_a_belt_between_two_planets_does_not_release_them_from_each_other():
    # the shape of GEN.1023: a giant, a belt, then a planet; the planet after the belt must stand ten mutual Hill radii from the giant too
    nodes, members = _tree_of("single")
    host = (planet_zones.S_TYPE, 0)
    giant, small = 1.9e27, 1.3e25
    book = ledger.Ledger(nodes, members, floors_au={0: 0.1})
    book.place_planet(host, 0.82, giant, 0.16)
    book.place_belt(host, 1.27, 1.63)
    first_outside = next(span.low for span in book.open_spans(host, small, 0.02) if span.low > 1.63)
    assert first_outside >= ledger.planet_exclusion(0.82, 0.16, giant, SOLAR, small, 0.02)[1] * (1.0 - 1e-12)
    hill = orbits.mutual_hill_radius_au(giant, small, 0.82, first_outside, SOLAR)
    gap = first_outside * (1 - 0.02) - 0.82 * (1 + 0.16)
    assert gap / hill >= tuning.MUTUAL_HILL_RADII_SEPARATION * (1 - 1e-9)


def test_no_room_is_an_empty_set():
    nodes, members = _tree_of("single")
    host = (planet_zones.S_TYPE, 0)
    book = ledger.Ledger(nodes, members, ceilings_au={host: 1.0}, floors_au={0: 0.9})
    book.place_planet(host, 0.95, 1e27)
    assert not book.open_spans(host, 1e27)
    assert ledger.log_uniform_quantile(0.5, book.open_spans(host, 1e27)) is None


# --- placing only where the ledger allows, and the verifier finds nothing ---------------------------------------------------------------------

def _star_row(star_id, mass_sol, node, slot):
    return {"id": star_id, "mass_kg": mass_sol * SOLAR, "radius_km": 0.005 * AU_KM, "parent_node_no": node, "slot": slot}


def _node_row(node):
    return {"node_no": node.node_no, "parent_node_no": node.parent_node_no, "slot": node.slot, "mass_solar": node.mass_solar,
            "separation_km": node.separation_km, "eccentricity": node.eccentricity, "inclination_deg": 0.0, "ascending_node_deg": 0.0,
            "arg_periapsis_deg": 0.0, "mean_anomaly_deg": 0.0, "period_years": 1.0, "kind": node.kind, "interaction": "detached", "flags": 0}


@pytest.mark.parametrize("kind", ["single", "wide", "close", "triple", "nested"])
def test_planets_and_belts_drawn_from_the_ledger_pass_the_independent_verifier(kind):
    nodes, members = _tree_of(kind)
    rng = random.Random(1014 + len(kind))
    for trial in range(40):
        book = ledger.Ledger(nodes, members, floors_au={m.member_no: 0.02 for m in members})
        rows = {"stars": [_star_row(m.member_no + 1, m.mass_solar, m.parent_node_no, m.slot) for m in members],
                "nodes": [_node_row(n) for n in nodes], "planets": [], "belts": []}
        for host in book.hosts():
            if host[0] == planet_zones.P_TYPE and sum(m.parent_node_no == host[1] for m in members) < 2:
                continue      # the stored rows name a circumbinary host only for a pair of two stars
            central = book.central_mass_kg(host)
            star_id = host[1] + 1 if host[0] == planet_zones.S_TYPE else None
            for _ in range(rng.randint(1, 7)):
                if rng.random() < 0.2:
                    spans = book.open_belt_spans(host, 0.3)
                    if not spans:
                        continue
                    low = ledger.log_uniform_quantile(rng.random(), ledger.SpanSet((s.low, s.high - 0.3) for s in spans))
                    if low is None:
                        continue
                    book.place_belt(host, low, low + 0.3)
                    rows["belts"].append({"id": len(rows["belts"]) + 1, "star_id": star_id, "lower_limit_km": low * AU_KM, "upper_limit_km": (low + 0.3) * AU_KM})
                    continue
                mass = EARTH * 10 ** rng.uniform(-1.0, 3.0)
                e = rng.uniform(0.0, 0.3)
                a = ledger.log_uniform_quantile(rng.random(), book.open_spans(host, mass, e))
                if a is None:
                    continue
                book.place_planet(host, a, mass, e)
                rows["planets"].append({"id": len(rows["planets"]) + 1, "star_id": star_id, "distance_km": a * AU_KM, "eccentricity": e,
                                        "radius_km": 6371.0 * (mass / EARTH) ** 0.3, "mass_kg": mass})
            assert central > 0.0
        found = [v for v in clearance.hard_violations(rows)]
        assert not found, f"{kind} trial {trial}: " + "; ".join(str(v) for v in found)
