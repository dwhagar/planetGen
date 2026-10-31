"""
TEST.477: `planetgen.physics.clearance`, the independent overlap verifier (GEN.1006 4.7). Three groups: each form the verifier
re-derives equals the generator's own copy of the same published formula (so it cannot drift from the paper while the
generator stays right, or the reverse); `audit` finds each kind of overlap in hand-made rows and nothing in a clear system;
and a seed sweep builds and stores systems of one to five stars and reads them back through `stored_rows`, requiring no
violation but the known ones, which are listed with their seeds (`KNOWN`) and expected to fail (strictly: a seed that
passes again must leave the list).
"""

import math
import os
import random

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.db import store
from planetgen.physics import clearance, constants, orbits, roche, stability
from planetgen.util import draw

AU_KM = constants.AU_TO_KM
SOLAR = constants.SOLAR_MASS_TO_KG


# --- the forms equal the generator's ----------------------------------------------------------------------------------

@pytest.mark.parametrize("masses", [(1.0, 1.0), (0.3, 2.0), (5.0, 0.05), (1.0, 100.0)])
def test_the_roche_lobe_is_the_generators(masses):
    assert clearance.roche_lobe_fraction(*masses) == pytest.approx(roche.roche_lobe_fraction(*masses), rel=1e-12)


def test_the_three_stability_forms_are_the_generators():
    rng = random.Random(477)
    for _ in range(200):
        q_in, q_out = rng.uniform(0.02, 1.0), 10 ** rng.uniform(-1.5, 1.5)
        e_in, e_out, i = rng.uniform(0.0, 0.8), rng.uniform(0.0, 0.8), rng.uniform(0.0, math.pi)
        assert clearance.mardling_aarseth_axis_ratio(q_out, e_out, i) == pytest.approx(stability.mardling_aarseth_axis_ratio(q_out, e_out, i), rel=1e-12)
        assert clearance.eggleton_kiseleva_axis_ratio(q_out, q_in, e_in, e_out) == pytest.approx(
            stability.eggleton_kiseleva_ycrit(q_out, q_in) * (1.0 + e_in) / (1.0 - e_out), rel=1e-12)
        assert clearance.vynatheya_axis_ratio(q_out, e_in, e_out, i) == pytest.approx(stability.eq4_axis_ratio(q_out, e_in, e_out, i), rel=1e-12)


def test_the_holman_wiegert_forms_are_the_generators():
    for mu in (0.1, 0.3, 0.5, 0.9):
        for e in (0.0, 0.2, 0.5, 0.8):
            assert clearance.holman_wiegert_s_axis_ratio(mu, e) * 7.0 == pytest.approx(
                orbits.holman_wiegert_critical_semimajor_axis(7.0, mu, e), rel=1e-12)
    for mu in (0.1, 0.3, 0.5):
        for e in (0.0, 0.2, 0.5, 0.7):
            assert clearance.holman_wiegert_p_axis_ratio(mu, e) * 7.0 == pytest.approx(
                orbits.holman_wiegert_circumbinary_a_crit_au(7.0, mu, e), rel=1e-12)


# --- audit on hand-made rows ------------------------------------------------------------------------------------------

def _star(star_id, mass_sol, radius_au, node, slot):
    return {"id": star_id, "mass_kg": mass_sol * SOLAR, "radius_km": radius_au * AU_KM, "parent_node_no": node, "slot": slot}


def _node(node_no, parent, slot, mass, axis_au, e=0.0, inclination=0.0, interaction="detached", kind="close"):
    return {"node_no": node_no, "parent_node_no": parent, "slot": slot, "mass_solar": mass, "separation_km": axis_au * AU_KM,
            "eccentricity": e, "inclination_deg": inclination, "ascending_node_deg": 0.0, "arg_periapsis_deg": 0.0,
            "mean_anomaly_deg": 0.0, "period_years": 1.0, "kind": kind, "interaction": interaction, "flags": 0}


def _pair(axis_au=1.0, e=0.0, radius_au=0.005):
    return {"stars": [_star(1, 1.0, radius_au, 0, 0), _star(2, 1.0, radius_au, 0, 1)], "nodes": [_node(0, None, 0, 2.0, axis_au, e)]}


def _triple(inner_au=1.0, outer_au=100.0, e_out=0.3, radius_au=0.005):
    return {"stars": [_star(1, 1.0, radius_au, 0, 0), _star(2, 1.0, radius_au, 0, 1), _star(3, 1.0, radius_au, 1, 1)],
            "nodes": [_node(1, None, 0, 3.0, outer_au, e_out), _node(0, 1, 0, 2.0, inner_au)]}


def _kinds(rows):
    return sorted({violation.kind for violation in clearance.audit(rows)})


def test_a_clear_pair_and_a_clear_triple_hold_no_violation():
    assert clearance.audit(_pair()) == []
    assert clearance.audit(_triple()) == []
    assert clearance.audit({"stars": [_star(1, 1.0, 0.005, None, 0)]}) == []


def test_two_stars_that_fill_their_lobes_at_closest_approach_overlap():
    assert _kinds(_pair(axis_au=0.02, radius_au=0.01)) == ["roche"]
    assert _kinds(_pair(axis_au=1.0, e=0.99, radius_au=0.005)) == ["roche"]


def test_a_merger_product_or_post_common_envelope_pair_is_not_held_to_the_detached_rule():
    rows = _pair(axis_au=0.02, radius_au=0.01)
    rows["nodes"][0]["interaction"] = "post-common-envelope"
    assert clearance.audit(rows) == []


def _with_a_dark_third(outer_au, e_out=0.3, mass=0.05, radius_au=0.0005):
    rows = _triple(inner_au=1.0, outer_au=outer_au, e_out=e_out)
    rows["stars"][2] = dict(_star(-3, mass, radius_au, 1, 1), dark=True, kind="brown-dwarf")
    rows["nodes"][0]["mass_solar"] = 2.0 + mass
    return rows


def test_a_dark_member_weighs_in_the_stability_of_the_tree_and_takes_up_room_in_it():
    assert clearance.audit(_with_a_dark_third(outer_au=100.0)) == []
    tight = _with_a_dark_third(outer_au=2.0, e_out=0.0, mass=1.0)
    assert "stability" in _kinds(tight)
    wide_mass = _with_a_dark_third(outer_au=100.0)
    wide_mass["nodes"][0]["mass_solar"] = 5.0
    assert "mass" in _kinds(wide_mass)


def test_no_planet_goes_round_a_dark_member_and_a_lone_star_with_one_is_the_planets_host():
    rows = {"stars": [_star(1, 1.0, 0.005, 0, 0), dict(_star(-2, 0.05, 0.0005, 0, 1), dark=True, kind="brown-dwarf")],
            "nodes": [_node(0, None, 0, 1.05, 30.0, 0.1, kind="wide")]}
    host = clearance._host_of({"star_id": None}, clearance._Tree(rows["stars"], rows["nodes"]))
    assert host[0] == "star" and host[1]["id"] == 1


def test_a_node_that_does_not_clear_its_children_holds_them_inside_each_other():
    assert "overlap" in _kinds(_triple(inner_au=1.0, outer_au=0.4, e_out=0.0))
    assert "overlap" not in _kinds(_triple(inner_au=1.0, outer_au=3.0, e_out=0.0))


def test_a_triple_too_close_to_be_stable_is_found():
    kinds = _kinds(_triple(inner_au=1.0, outer_au=2.0, e_out=0.0))
    assert "stability" in kinds
    assert "stability" not in _kinds(_triple(inner_au=1.0, outer_au=30.0, e_out=0.0))


def test_a_node_whose_mass_is_not_its_children_is_found():
    rows = _pair()
    rows["nodes"][0]["mass_solar"] = 3.0
    assert _kinds(rows) == ["mass"]


def _with_planets(rows, *planets, **more):
    return dict(rows, planets=[dict({"id": n + 1, "star_id": 1, "eccentricity": 0.0, "radius_km": 6371.0, "mass_kg": 6e24}, **p)
                               for n, p in enumerate(planets)], **more)


def test_a_planet_inside_its_star_or_its_roche_limit_is_found():
    rows = _with_planets(_pair(axis_au=300.0), {"distance_km": 0.001 * AU_KM})
    assert "planet-in-star" in _kinds(rows)
    assert "planet-in-star" not in _kinds(_with_planets(_pair(axis_au=300.0), {"distance_km": 1.0 * AU_KM}))


def test_a_planet_past_the_holman_wiegert_limit_of_its_star_is_found():
    wide = _pair(axis_au=10.0)           # equal masses, circular: critical axis 0.464 - 0.19 = 0.274 of 10 AU
    assert "s-type" in _kinds(_with_planets(wide, {"distance_km": 3.0 * AU_KM}))
    assert "s-type" not in _kinds(_with_planets(wide, {"distance_km": 2.0 * AU_KM}))


def test_a_circumbinary_planet_inside_the_p_type_limit_is_found():
    close = _pair(axis_au=0.2)           # equal masses, circular: critical axis 1.60 + 2.06 - 1.2725 = 2.3875 of 0.2 AU
    assert "p-type" in _kinds(_with_planets(close, {"star_id": None, "distance_km": 0.3 * AU_KM}))
    assert "p-type" not in _kinds(_with_planets(close, {"star_id": None, "distance_km": 1.0 * AU_KM}))


def test_two_planets_that_cross_or_crowd_each_other_are_found():
    base = _pair(axis_au=300.0)
    crossing = _with_planets(base, {"distance_km": 1.0 * AU_KM, "eccentricity": 0.5}, {"distance_km": 1.4 * AU_KM})
    assert "planet-crossing" in _kinds(crossing)
    crowded = _with_planets(base, {"distance_km": 1.0 * AU_KM}, {"distance_km": 1.05 * AU_KM})
    assert "planet-spacing" in _kinds(crowded)
    apart = _with_planets(base, {"distance_km": 1.0 * AU_KM}, {"distance_km": 2.5 * AU_KM})
    assert not {"planet-crossing", "planet-spacing"} & set(_kinds(apart))


def test_a_moon_outside_its_planets_hill_sphere_or_inside_its_roche_limit_is_found():
    rows = _with_planets(_pair(axis_au=300.0), {"distance_km": 1.0 * AU_KM, "mass_kg": 6e24})
    hill_km = 1.0 * AU_KM * (6e24 / (3 * 2.0 * SOLAR / 2.0)) ** (1.0 / 3.0)
    far = dict(rows, moons=[{"id": 1, "planet_id": 1, "distance_km": 0.6 * hill_km, "eccentricity": 0.0, "radius_km": 1737.0, "mass_kg": 7e22}])
    assert "moon-hill" in _kinds(far)
    inside = dict(rows, moons=[{"id": 1, "planet_id": 1, "distance_km": 9000.0, "eccentricity": 0.0, "radius_km": 1737.0, "mass_kg": 7e22}])
    assert "moon-in-planet" in _kinds(dict(rows, moons=[{"id": 1, "planet_id": 1, "distance_km": 7000.0, "eccentricity": 0.0, "radius_km": 1737.0, "mass_kg": 7e22}]))
    assert "moon-roche" in _kinds(inside)
    assert [item.hard for item in clearance.audit(inside) if item.kind == "moon-roche"] == [False] and clearance.hard_violations(inside) == []
    fine = dict(rows, moons=[{"id": 1, "planet_id": 1, "distance_km": 380000.0, "eccentricity": 0.0, "radius_km": 1737.0, "mass_kg": 7e22}])
    assert not {"moon-hill", "moon-roche", "moon-in-planet"} & set(_kinds(fine))


def test_a_belt_holding_a_planets_orbit_or_past_the_stability_limit_is_found():
    rows = _with_planets(_pair(axis_au=300.0), {"distance_km": 3.0 * AU_KM})
    belt = {"id": 1, "star_id": 1, "lower_limit_km": 2.5 * AU_KM, "upper_limit_km": 3.5 * AU_KM}
    assert "belt-planet" in _kinds(dict(rows, belts=[belt]))
    wide = dict(_with_planets(_pair(axis_au=10.0)), belts=[{"id": 1, "star_id": 1, "lower_limit_km": 2.5 * AU_KM, "upper_limit_km": 4.0 * AU_KM}])
    assert "belt-s-type" in _kinds(wide)


def test_a_limit_outside_the_fit_is_applied_at_its_edge_and_says_so():
    light = {"stars": [_star(1, 1.0, 0.005, 0, 0), _star(2, 0.05, 0.001, 0, 1)], "nodes": [_node(0, None, 0, 1.05, 10.0)]}
    found = clearance.audit(_with_planets(light, {"distance_km": 4.5 * AU_KM}))
    assert found and all(item.extrapolated for item in found if item.kind == "s-type")


# --- the sweep over generated, stored systems -------------------------------------------------------------------------

SCALE = int(os.environ.get("PLANETGEN_CLEARANCE_SCALE", "1"))
"""int: How many times the seeds of `SWEEP` (`PLANETGEN_CLEARANCE_SCALE=10` for the long sweep that finds the seeds of `KNOWN`)."""

SWEEP = {"single": (None, 12 * SCALE), "pair": (2, 12 * SCALE), "triple": (3, 24 * SCALE), "quadruple": (4, 16 * SCALE), "five": (5, 8 * SCALE)}
"""dict: Systems per kind of star count (`SystemConfig.STARS`, None for the generator's own mix) and how many seeds each."""

KNOWN = {}
"""dict: `(kind, seed) -> what the verifier finds` in the systems that are not clear today (GEN.1006 P8 is what removes them). A seed
here is always swept, and expected to fail: when it passes, leave the list."""


def _cases():
    out = []
    for kind, (_stars, count) in SWEEP.items():
        seeds = sorted(set(range(1, count + 1)) | {seed for (known_kind, seed) in KNOWN if known_kind == kind})
        for seed in seeds:
            marks = pytest.mark.xfail(strict=True, reason=KNOWN[(kind, seed)]) if (kind, seed) in KNOWN else ()
            out.append(pytest.param(kind, seed, id=f"{kind}-{seed}", marks=marks))
    return out


@pytest.mark.db
@pytest.mark.slow
@pytest.mark.parametrize("kind,seed", _cases())
def test_a_generated_system_holds_no_overlap(kind, seed, mysql_config):
    cfg = SystemConfig()
    cfg.STARS = SWEEP[kind][0]
    cfg.PLANETS = True
    with draw.bound(seed):
        system = StarSystem(system_config=cfg)
    system_id = store.save_system(system, cfg, mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        rows = clearance.stored_rows(conn, system_id)
    finally:
        conn.close()
    found = clearance.hard_violations(rows)
    assert not found, "\n".join(str(item) for item in found)
