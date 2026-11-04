"""
MAP.275: the look of a planet's stored ring (`web/maps/ringlook.py`): the colour and brightness follow the ice share by the GEN.670 rules, the opacity
follows the optical depth, and a moon digs a gap where its resonance falls in the ring.
"""

import pytest

from planetgen.web.maps import bodylook, ringlook

SATURN_MASS_KG = 5.683e26
SATURN_RADIUS_KM = 58232.0


def _planet(**ring):
    row = {"id": 1, "body_type": "g", "planet_class": "G", "radius_km": SATURN_RADIUS_KM, "mass_kg": SATURN_MASS_KG,
           "ring_kind": "broad", "ring_origin": "disrupted_moon", "ring_inner_km": 74500.0, "ring_outer_km": 140200.0,
           "ring_optical_depth": 0.5, "ring_ice_fraction": 0.99}
    row.update(ring)
    return row


def test_albedo_is_the_icy_rings_down_to_dust():
    assert ringlook.albedo(0.99) == ringlook.ALBEDO_ICE
    assert ringlook.albedo(0.95) == ringlook.ALBEDO_ICE
    assert ringlook.albedo(0.5) == pytest.approx(ringlook.ALBEDO_HALF)
    assert ringlook.albedo(0.0) == ringlook.ALBEDO_DUST
    assert ringlook.albedo(0.2) < ringlook.albedo(0.7)


def test_an_icy_thick_ring_reads_cream_yellow_and_a_thin_one_tan():
    """The B ring: blue about 0.65 of green, red equal to green. The C ring: blue 0.65, red 1.02 to 1.03 of green (research section 5)."""
    r, g, b = ringlook.rgb(0.99, 0.8)
    assert b / g == pytest.approx(0.65, abs=0.03) and r / g == pytest.approx(1.0, abs=0.02)
    r, g, b = ringlook.rgb(0.99, 0.05)
    assert r / g == pytest.approx(1.025, abs=0.02) and b / g == pytest.approx(0.65, abs=0.03)


def test_blue_slope_is_flat_for_an_icy_ring_whatever_its_dust():
    """The A and B rings hold 0.2 to 0.5 percent dust and are still at the full slope: it does not scale with the share."""
    assert ringlook.slopes(0.995, 0.5)[0] == ringlook.BLUE_SLOPE_ICY
    assert ringlook.slopes(0.985, 0.5)[0] == ringlook.BLUE_SLOPE_ICY
    assert ringlook.slopes(0.0, 0.5)[0] == ringlook.BLUE_SLOPE_DUST
    assert ringlook.slopes(0.7, 0.5)[0] < ringlook.BLUE_SLOPE_ICY


def test_a_dust_ring_is_dark_and_still_drawn_faintly():
    r, g, b = ringlook.rgb(0.0, 1e-8)
    assert max(r, g, b) < 40
    assert ringlook.alpha(1e-8) == ringlook.DUST_ALPHA_FLOOR
    assert ringlook.alpha(0.5) == pytest.approx(0.393, abs=0.001)
    assert ringlook.alpha(50.0) == ringlook.MAX_ALPHA


def test_a_moons_two_to_one_resonance_gap_falls_at_its_radius():
    """A Mimas-like moon (a = 185,500 km, 3.8e19 kg) digs the Cassini Division at 0.63 of its distance, 117,000 km."""
    moons = [{"distance_km": 185500.0, "mass_kg": 3.75e19}]
    gaps = ringlook.resonance_gaps(74500.0, 140200.0, SATURN_MASS_KG, moons)
    two_to_one = [g for g in gaps if abs(g[0] - (185500.0 * 0.5 ** (2 / 3) - 74500.0) / 65700.0) < 1e-3]
    assert len(two_to_one) == 1
    centre, width, depth = two_to_one[0]
    assert 0.03 < width < 0.1 and 0.3 <= depth <= 0.95


def test_no_gap_for_a_moon_whose_resonances_fall_outside_the_ring_or_a_tiny_one():
    assert ringlook.resonance_gaps(74500.0, 140200.0, SATURN_MASS_KG, [{"distance_km": 1.2e6, "mass_kg": 1e22}]) == []
    assert ringlook.resonance_gaps(74500.0, 140200.0, SATURN_MASS_KG, [{"distance_km": 185500.0, "mass_kg": 1e12}]) == []
    assert ringlook.resonance_gaps(74500.0, 140200.0, SATURN_MASS_KG, [{"distance_km": None, "mass_kg": 1e19}]) == []
    assert ringlook.resonance_gaps(140200.0, 74500.0, SATURN_MASS_KG, [{"distance_km": 185500.0, "mass_kg": 1e19}]) == []


def test_a_big_moon_clears_at_most_a_quarter_of_the_ring_per_resonance():
    gaps = ringlook.resonance_gaps(74500.0, 140200.0, SATURN_MASS_KG, [{"distance_km": 185500.0, "mass_kg": 1e24}])
    assert gaps and all(g[1] <= ringlook.GAP_MAX_SHARE for g in gaps)


def test_ring_look_gives_edges_in_radii_and_none_without_a_ring():
    look = ringlook.ring_look(_planet(moons=[{"distance_km": 185500.0, "mass_kg": 3.75e19}]))
    assert look["in"] == pytest.approx(74500.0 / SATURN_RADIUS_KM, abs=1e-3) and look["out"] == pytest.approx(140200.0 / SATURN_RADIUS_KM, abs=1e-3)
    assert look["kind"] == "broad" and look["gaps"] and len(look["rgb"]) == 3
    assert ringlook.ring_look(_planet(ring_kind=None)) is None
    assert ringlook.ring_look(_planet(ring_inner_km=None)) is None
    assert ringlook.ring_look(_planet(ring_outer_km=70000.0)) is None


def test_the_look_record_carries_the_ring_only_when_the_body_stores_one():
    assert bodylook.planet_look(_planet())["ring"]["kind"] == "broad"
    plain = _planet(ring_kind=None, ring_origin=None, ring_inner_km=None, ring_outer_km=None, ring_optical_depth=None, ring_ice_fraction=None)
    assert "ring" not in bodylook.planet_look(plain)
