# tests/test_giant_interior.py

"""
GEN.642: the interiors of gas giants, ice giants and brown dwarfs (`planetgen.physics.giant_interior`): the layers add up to
the body, hold exactly the stated heavy elements and natal helium, pressures rise downward, cold giants rain helium out of
the molecular envelope into a gradient and a shell, big giants get a dilute core and small ones a compact core, ice giants
get an ionic or superionic mantle over a rock core, brown dwarfs have no solid core, and the draw is deterministic.
"""

import math

import pytest

from planetgen.physics import giant_interior as gi
from planetgen.util import draw

ICES = {"water_liquid": 0.25, "nh3_ice": 0.02, "ch4_ice": 0.03}
GIANT_ICES = {"water_liquid": 0.6, "nh3_ice": 0.1, "ch4_ice": 0.1}

JUPITER = dict(mass_earth=317.8, radius_km=71492.0, heavy_mass_earth=42.0, envelope_z=0.04, helium_share=0.2745,
               temperature_1bar_k=166.0, ice_shares=ICES)
SATURN = dict(mass_earth=95.2, radius_km=60268.0, heavy_mass_earth=25.0, envelope_z=0.04, helium_share=0.2745,
              temperature_1bar_k=140.0, ice_shares=ICES)
NEPTUNE = dict(mass_earth=17.1, radius_km=24764.0, heavy_mass_earth=14.6, envelope_z=0.1, helium_share=0.2745,
               temperature_1bar_k=72.0, ice_shares=GIANT_ICES)
BROWN_DWARF = dict(mass_earth=40 * 317.8, radius_km=0.9 * 71492.0, heavy_mass_earth=40 * 317.8 * 0.015, envelope_z=0.015,
                   helium_share=0.2745, temperature_1bar_k=1000.0, is_brown_dwarf=True)


def interior(spec, seed=1):
    with draw.bound(seed):
        return gi.draw_interior(gi.GiantSubject(**spec))


def kinds(body):
    return [layer.kind for layer in body.layers]


@pytest.mark.parametrize("spec", [JUPITER, SATURN, NEPTUNE, BROWN_DWARF], ids=["jupiter", "saturn", "neptune", "brown dwarf"])
@pytest.mark.parametrize("seed", range(1, 9))
def test_every_draw_keeps_the_invariants(spec, seed):
    body = interior(spec, seed)
    body.validate(spec["heavy_mass_earth"], spec["helium_share"])
    assert math.isclose(math.fsum(layer.mass_earth for layer in body.layers), spec["mass_earth"], rel_tol=1e-9)
    shares = body.whole_body_shares()
    assert math.isclose(shares["metal"] + shares["rock"] + shares["water"] + shares["ice"] + shares["envelope"], 1.0,
                        rel_tol=1e-9)
    assert all(0.0 <= v <= 1.0 for v in shares.values())
    for layer in body.layers:
        assert 0.0 <= layer.bottom_radius_fraction <= layer.top_radius_fraction <= 1.0


def test_the_layers_run_from_the_surface_down_and_the_radii_close():
    body = interior(JUPITER)
    assert body.layers[0].top_radius_fraction == 1.0
    assert body.layers[-1].bottom_radius_fraction == pytest.approx(0.0, abs=1e-6)
    for upper, lower in zip(body.layers, body.layers[1:]):
        assert upper.bottom_radius_fraction == pytest.approx(lower.top_radius_fraction, abs=1e-6)


def test_the_draw_is_deterministic_and_the_seed_matters():
    assert interior(JUPITER, 4) == interior(JUPITER, 4)
    assert len({interior(JUPITER, seed).metallization_pressure_gpa for seed in range(1, 6)}) > 1


def test_a_jupiter_has_molecular_over_metallic_hydrogen_a_gradient_and_a_dilute_core():
    body = interior(JUPITER)
    assert kinds(body)[0] == "molecular_envelope" and body.layer("metallic_envelope")
    assert body.layer("molecular_envelope").state == "molecular_fluid"
    assert body.layer("metallic_envelope").state == "metallic_fluid"
    low, high = gi.METALLIZATION_PRESSURE_RANGE_GPA
    assert low <= body.metallization_pressure_gpa <= high
    assert body.layer("helium_gradient") is not None


def test_a_cold_giant_rains_helium_into_a_gradient_and_a_shell():
    saturn = interior(SATURN, 5)
    assert saturn.helium_depletion < 1.0
    assert saturn.layer("helium_gradient") and saturn.has_helium_shell
    molecular = saturn.layer("molecular_envelope")
    shell = saturn.layer("helium_shell")
    natal = SATURN["helium_share"]
    molecular_helium = molecular.components["helium"] / (molecular.mass_earth - molecular.heavy_mass_earth())
    assert molecular_helium < natal * 0.5
    assert shell.components["helium"] / shell.mass_earth > 0.5
    warm = interior(JUPITER, 5)
    assert not warm.has_helium_shell


def test_helium_depletion_is_one_when_warm_and_falls_when_cold():
    assert gi.helium_depletion(300.0) == 1.0
    assert gi.helium_depletion(100.0) == gi.HELIUM_RAIN_POINTS_K[0][1]
    assert gi.helium_depletion(150.0) < gi.helium_depletion(160.0) < 1.0


def test_a_big_giant_can_have_a_dilute_core_and_a_small_one_a_compact_core():
    big = [interior(JUPITER, seed) for seed in range(1, 12)]
    assert any(b.dilute_core_mass_earth > 0 for b in big)
    for b in big:
        assert b.dilute_core_mass_earth <= gi.DILUTE_CORE_MAX_MASS_FRACTION * b.total_mass_earth
    small = dict(mass_earth=25.0, radius_km=38000.0, heavy_mass_earth=10.0, envelope_z=0.05, helium_share=0.2745,
                 temperature_1bar_k=120.0, ice_shares=ICES)
    for seed in range(1, 12):  # a small giant draws a dilute core only for heavy elements that have nowhere else to go
        body = interior(small, seed)
        dilute = body.layer("dilute_core")
        if dilute:
            assert dilute.heavy_mass_earth() <= gi.DILUTE_CORE_MAX_HEAVY_FRACTION * dilute.mass_earth * (1 + 1e-9)
    cores = [interior(small, seed).core_mass_earth for seed in range(1, 12)]
    assert all(c == 0.0 or gi.COMPACT_CORE_RANGE_EARTH[0] <= c <= gi.COMPACT_CORE_RANGE_EARTH[1] for c in cores)


def test_the_core_layers_hold_the_rock_and_iron_and_the_envelope_the_ices_first():
    body = interior(JUPITER)
    shares = body.whole_body_shares()
    assert shares["metal"] + shares["rock"] > 0.0
    assert shares["envelope_z"] > 0.0


def test_an_ice_giant_has_a_superionic_or_ionic_mantle_over_a_rock_core():
    body = interior(NEPTUNE)
    assert body.kind == "ice_giant"
    assert kinds(body)[0] == "molecular_envelope" and "ice_mantle" in kinds(body) and body.layers[-1].kind == "rock_core"
    mantle = body.layer("ice_mantle")
    assert mantle.state in ("ionic_fluid", "superionic")
    assert body.layer("rock_core").state == "molten"
    assert body.layer("metallic_envelope") is None
    shares = body.whole_body_shares()
    assert shares["water"] > 0.0 and shares["rock"] > 0.0
    assert body.layer("rock_core").components.keys() <= {"iron", "rock"}


def test_a_brown_dwarf_has_no_solid_core_and_stays_metallic_below_the_transition():
    body = interior(BROWN_DWARF)
    assert body.kind == "brown_dwarf"
    assert body.core_mass_earth == 0.0 and body.dilute_core_mass_earth == 0.0
    assert set(kinds(body)) <= {"molecular_envelope", "metallic_envelope"}
    assert body.layers[-1].bottom_pressure_gpa > 1.0e5


def test_the_subject_refuses_nonsense():
    with pytest.raises(gi.InteriorError):
        gi.GiantSubject(**{**JUPITER, "heavy_mass_earth": 400.0})
    with pytest.raises(gi.InteriorError):
        gi.GiantSubject(**{**JUPITER, "helium_share": 1.2})
    with pytest.raises(gi.InteriorError):
        gi.GiantSubject(**{**JUPITER, "ice_shares": {"water_liquid": 0.8, "nh3_ice": 0.5}})
    with pytest.raises(gi.InteriorError):
        gi.GiantSubject(**{**JUPITER, "mass_earth": 0.0})


def test_validate_catches_broken_interiors():
    body = interior(JUPITER)
    first = body.layers[0]
    broken = gi.GiantInterior((gi.Layer(first.kind, first.state, first.top_pressure_gpa, first.bottom_pressure_gpa,
                                        first.mass_earth * 2, first.components),) + body.layers[1:],
                              body.total_mass_earth, body.kind, body.metallization_pressure_gpa, body.helium_depletion)
    with pytest.raises(gi.InteriorError):
        broken.validate()
    with pytest.raises(gi.InteriorError):
        body.validate(heavy_mass_earth=body.heavy_mass_earth() + 1.0)
    with pytest.raises(gi.InteriorError):
        body.validate(helium_share=0.5)


def test_the_polytrope_places_pressure_by_mass_and_the_thin_shell_formula_matches_the_limit():
    mass, radius = 1.9e27, 7.15e7
    profile = gi.Polytrope(1.0, mass, radius)
    r_centre, p_centre = profile.at_mass_inside(0.0)
    assert r_centre == pytest.approx(0.0, abs=1e-6)
    assert profile.at_mass_inside(1.0) == pytest.approx((1.0, 0.0), abs=1e-6)
    assert p_centre == pytest.approx(math.pi * 6.674e-11 * mass ** 2 / (8 * radius ** 4), rel=0.02)
    _, p_half = profile.at_mass_inside(0.5)
    _, p_tenth = profile.at_mass_inside(0.1)
    assert 0 < p_half < p_tenth < p_centre
    assert gi.thin_shell_base_pressure_pa(1.0e25, 1.0e24, 6.0e6) > 0.0


def test_the_compact_core_table_is_registered_with_open_rows():
    from planetgen.generation import possibility_table as pt
    table = pt.REGISTRY.get("giant-compact-core")
    assert [row.name for row in table.rows] == ["compact core", "no compact core"]


def test_the_description_has_a_sentence_for_every_layer():
    body = interior(JUPITER)
    sentences = gi.describe_interior(body)
    assert len(sentences) == len(body.layers)
    assert sentences[0].startswith("The molecular hydrogen and helium envelope (molecular fluid) holds")
    assert all(s.endswith(".") and any(word in s for word in ("hydrogen", "rock", "iron")) for s in sentences)
    ice = gi.describe_interior(interior(NEPTUNE))
    assert any("mantle of water, ammonia and methane" in s for s in ice)
