# tests/test_collision_outcome.py

"""
GEN.345: the outcome table of docs/design/collisions-and-mergers.md section 5.1 as one pure function. Every row has a
case, mass and momentum are conserved to a stated tolerance, the same input gives the same outcome in either order, and a
borderline case is flagged.
"""

import math

import pytest

from planetgen import tuning
from planetgen.physics import collision_outcome as co
from planetgen.physics import constants
from planetgen.physics import planets, stellar_abundances

EARTH, EARTH_R = constants.EARTH_MASS_TO_KG, 6.371e6
JUP = constants.JUPITER_MASS_TO_KG
SUN = constants.SOLAR_MASS_TO_KG
KMS = 1000.0
CONSERVED = 1e-12


def _earth(key, mass=1.0, speed_kms=0.0, **rest):
    radius = EARTH_R * mass ** (1.0 / 3.0)
    return co.Body(key, "T", mass * EARTH, radius, velocity_ms=(speed_kms * KMS, 0.0, 0.0), **rest)


def _giant(key, mass_jup, speed_kms=0.0, kind=None, z=stellar_abundances.PROTO_Z):
    mass = mass_jup * JUP
    kind = kind or ("B" if mass_jup >= 13.0 else "G")
    return co.Body(key, kind, mass, planets.giant_radius_km(mass) * 1000.0, velocity_ms=(speed_kms * KMS, 0.0, 0.0),
                   metallicity_z=z)


def _assert_conserved(a, b, outcome):
    mass, momentum = co.conservation_error(a, b, outcome)
    assert mass < CONSERVED and momentum < CONSERVED, (outcome.code, mass, momentum)


# --- The numbers ---------------------------------------------------------------------------------------

def test_the_catastrophic_energy_matches_the_documents_table():
    assert co.catastrophic_energy(EARTH, EARTH) == pytest.approx(6.4e7, rel=0.01)
    assert co.catastrophic_energy(EARTH, 0.1 * EARTH) == pytest.approx(1.1e8, rel=0.02)
    assert co.catastrophic_energy(EARTH, 0.01 * EARTH) == pytest.approx(6.4e8, rel=0.01)
    assert co.catastrophic_energy(JUP, JUP) == pytest.approx(3.0e9, rel=0.01)


def test_the_head_on_catastrophic_speed_of_two_earths_is_22_6_km_per_s():
    a, b = _earth("a"), _earth("b", 1.0)
    v_esc = co.escape_speed(2 * EARTH, 2 * EARTH_R)
    assert v_esc == pytest.approx(11.2e3, rel=0.01)
    v_cat = math.sqrt(8.0 * co.catastrophic_energy(EARTH, EARTH))  # Q_R = Q* for equal masses: mu / M = 1/4
    assert v_cat == pytest.approx(22.6e3, rel=0.01)
    below = co.resolve(a, co.Body("b", "T", b.mass_kg, b.radius_m, velocity_ms=(math.sqrt(v_cat ** 2 - v_esc ** 2) * 0.99, 0, 0)), 0.0, 0.0)
    above = co.resolve(a, co.Body("b", "T", b.mass_kg, b.radius_m, velocity_ms=(math.sqrt(v_cat ** 2 - v_esc ** 2) * 1.01, 0, 0)), 0.0, 0.0)
    assert below.numbers["x"] < 1.0 <= above.numbers["x"]
    assert below.code == "erosion" and above.code == "disruption"


def test_the_lens_fraction_runs_from_nothing_to_the_whole_body():
    assert co.lens_fraction(1.0, 0.0) == 0.0
    assert co.lens_fraction(1.0, 2.0) == 1.0
    assert co.lens_fraction(1.0, 1.0) == pytest.approx(0.5)
    assert co.lens_fraction(1.0, 3.0) == 1.0 and co.lens_fraction(1.0, -1.0) == 0.0


def test_the_largest_remnant_fraction_is_continuous():
    assert co.largest_remnant_fraction(0.0) == 1.0
    assert co.largest_remnant_fraction(1.0) == pytest.approx(0.5)
    assert co.largest_remnant_fraction(co.HERO_FADE_X) == pytest.approx(0.1)
    for edge in (1.0, co.HERO_FADE_X):
        assert co.largest_remnant_fraction(edge - 1e-9) == pytest.approx(co.largest_remnant_fraction(edge + 1e-9), rel=1e-6)
    assert co.largest_remnant_fraction(10.0) < co.largest_remnant_fraction(2.0)


# --- Two terrestrial bodies: every row -----------------------------------------------------------------------

def _tt(speed_kms, impact_over_sum, mass_b=1.0, **settings):
    a, b = _earth("a"), _earth("b", mass_b, speed_kms)
    return a, b, co.resolve(a, b, impact_over_sum * (a.radius_m + b.radius_m), 100.0, **settings)


def test_a_slow_head_on_pass_merges_into_one_body_with_the_radius_by_volume():
    a, b, o = _tt(1.0, 0.0)
    assert o.code == "merger" and o.destroyed == ("b",) and o.debris is None and o.lost_mass_kg == 0.0
    (s,) = o.survivors
    assert s.key == "a" and s.mass_kg == pytest.approx(2 * EARTH) and s.radius_m == pytest.approx(EARTH_R * 2 ** (1 / 3))
    assert s.velocity_ms == pytest.approx((0.5e3, 0.0, 0.0))
    _assert_conserved(a, b, o)


def test_a_slow_graze_merges():
    a, b, o = _tt(1.0, 0.9)
    assert o.code == "merger" and o.numbers["b"] > o.numbers["b_crit"]  # the speed rule takes it before the angle
    _assert_conserved(a, b, o)


def test_a_graze_between_the_merge_and_hit_and_run_speeds_graze_and_merges():
    v_esc = co.escape_speed(2 * EARTH, 2 * EARTH_R) / KMS
    speed = math.sqrt((1.35 * v_esc) ** 2 - v_esc ** 2)  # v_imp = 1.35 v_esc
    a, b, o = _tt(speed, 0.9)
    assert o.numbers["v_imp"] / (v_esc * KMS) == pytest.approx(1.35, rel=1e-6)
    assert o.code == "graze_and_merge" and len(o.survivors) == 1 and o.debris is None
    _assert_conserved(a, b, o)


def test_a_fast_graze_is_a_hit_and_run_that_leaves_both_alive_and_strips_a_field():
    a, b, o = _tt(20.0, 0.6)
    assert o.code == "hit_and_run" and o.destroyed == ()
    assert [s.key for s in o.survivors] == ["a", "b"]
    assert all(s.mass_kg < EARTH for s in o.survivors)
    assert o.debris is not None and o.debris.mass_kg >= tuning.COLLISION_FIELD_MIN_SHARE * 2 * EARTH
    assert o.debris.mass_kg == pytest.approx(2 * EARTH - sum(s.mass_kg for s in o.survivors))
    _assert_conserved(a, b, o)


def test_a_hit_and_run_that_strips_under_a_percent_is_only_logged():
    a, b, o = _tt(20.0, 0.97)
    assert o.code == "hit_and_run" and o.debris is None and o.lost_mass_kg > 0.0
    _assert_conserved(a, b, o)


def test_a_head_on_pass_at_moderate_speed_erodes_to_the_largest_remnant_and_a_field():
    a, b, o = _tt(10.0, 0.0)
    assert o.code == "erosion" and 0.1 <= o.numbers["x"] < 1.0
    (s,) = o.survivors
    assert s.mass_kg == pytest.approx(2 * EARTH * (1.0 - o.numbers["x"] / 2.0))
    assert o.debris.mass_kg == pytest.approx(2 * EARTH - s.mass_kg) and o.debris.hero_mass_kg == 0.0
    _assert_conserved(a, b, o)


def test_a_head_on_pass_past_the_catastrophic_energy_destroys_both_into_one_field_with_a_hero_body():
    a, b, o = _tt(40.0, 0.0)
    assert o.code == "disruption" and o.survivors == () and set(o.destroyed) == {"a", "b"}
    assert o.debris.mass_kg == pytest.approx(2 * EARTH)
    assert o.debris.hero_mass_kg == pytest.approx(co.largest_remnant_fraction(o.numbers["x"]) * 2 * EARTH)
    assert o.debris.composition == "stony" and o.debris.radius_m == pytest.approx(a.radius_m + b.radius_m)
    _assert_conserved(a, b, o)


def test_the_literal_rule_destroys_any_contact_and_leaves_no_hero():
    a, b, o = _tt(1.0, 0.0, always_destroy_terrestrial=True)
    assert o.code == "disruption" and o.survivors == () and o.debris.hero_mass_kg == 0.0
    assert o.debris.mass_kg == pytest.approx(2 * EARTH)
    _assert_conserved(a, b, o)


def test_a_field_of_unlike_parents_is_mixed_and_of_like_parents_is_theirs():
    a = _earth("a", composition="icy")
    like = co.resolve(a, _earth("b", 1.0, 40.0, composition="icy"), 0.0, 0.0)
    unlike = co.resolve(a, _earth("b", 1.0, 40.0, composition="carbonaceous"), 0.0, 0.0)
    assert like.debris.composition == "icy" and unlike.debris.composition == "mixed"


def test_a_small_body_on_a_large_one_at_speed_is_not_catastrophic():
    a, b = _earth("a"), _earth("b", 0.01, 30.0)
    o = co.resolve(a, b, 0.0, 0.0)
    assert o.numbers["x"] < 1.0 and o.code in ("merger", "erosion")
    _assert_conserved(a, b, o)


# --- Absorption, tides and the stars -------------------------------------------------------------------------

def test_a_terrestrial_body_hitting_a_giant_is_absorbed_and_the_giant_keeps_the_momentum():
    giant, rock = _giant("a", 1.0), _earth("b", 1.0, 30.0)
    o = co.resolve(giant, rock, 0.0, 0.0)
    assert o.code == "absorbed" and o.destroyed == ("b",) and o.debris is None and o.records == ("absorbed",)
    (s,) = o.survivors
    assert s.key == "a" and s.kind == "G" and s.mass_kg == pytest.approx(JUP + EARTH)
    assert s.radius_m == pytest.approx(planets.giant_radius_km(s.mass_kg) * 1000.0)
    _assert_conserved(giant, rock, o)


def test_a_terrestrial_body_hitting_a_star_is_absorbed_with_a_transient_flag():
    star = co.Body("a", "S", SUN, 6.957e8)
    rock = _earth("b", 1.0, 40.0)
    o = co.resolve(star, rock, 0.0, 0.0)
    assert o.code == "absorbed" and o.survivors[0].kind == "S" and "transient" in o.survivors[0].flags
    assert o.survivors[0].radius_m is None
    _assert_conserved(star, rock, o)


@pytest.mark.parametrize("compact, radius", [(1.4 * SUN, 1.2e4), (10 * SUN, 3.0e4)])
def test_a_rock_near_a_neutron_star_or_black_hole_is_tidally_disrupted_before_contact(compact, radius):
    x = co.Body("a", "X", compact, radius)
    rock = _earth("b", 1.0, 100.0)
    o = co.resolve(x, rock, 0.5 * (radius + rock.radius_m), 0.0)
    assert o.code == "tidal_disruption" and o.debris is None and o.destroyed == ("b",)
    assert o.survivors[0].kind == "X" and o.survivors[0].mass_kg == pytest.approx(compact + EARTH)
    _assert_conserved(x, rock, o)


def test_two_stars_merge_conserving_momentum_and_say_it_is_outside_the_table():
    a = co.Body("a", "S", SUN, 6.957e8, velocity_ms=(1e4, 0.0, 0.0))
    b = co.Body("b", "S", 0.5 * SUN, 4.0e8, velocity_ms=(-2e4, 0.0, 0.0))
    o = co.resolve(a, b, 0.0, 0.0)
    assert o.code == "stellar_merger" and o.survivors[0].mass_kg == pytest.approx(1.5 * SUN)
    assert any("outside the table" in f for f in o.survivors[0].flags)
    _assert_conserved(a, b, o)


# --- Giants, brown dwarfs and fusion -------------------------------------------------------------------------

def test_two_giants_merge_losing_a_small_share_as_gas_with_the_radius_from_the_relation():
    a, b = _giant("a", 1.0), _giant("b", 1.0, 30.0)
    o = co.resolve(a, b, 0.0, 0.0)
    assert o.code == "giant_merger" and o.debris is None and o.destroyed == ("b",)
    low, high = tuning.GIANT_MERGER_LOSS_RANGE
    assert low * 2 * JUP <= o.lost_mass_kg <= high * 2 * JUP
    (s,) = o.survivors
    assert s.mass_kg == pytest.approx(2 * JUP - o.lost_mass_kg) and s.kind == "G"
    assert s.radius_m == pytest.approx(planets.giant_radius_km(s.mass_kg) * 1000.0)
    assert {"rapidly rotating", "inflated", "hot"} <= set(s.flags)
    _assert_conserved(a, b, o)


def test_two_giants_can_make_a_brown_dwarf_and_it_is_recorded():
    a, b = _giant("a", 10.0), _giant("b", 10.0)
    o = co.resolve(a, b, 0.0, 0.0)
    assert o.survivors[0].kind == "B" and o.records == ("giant_merger", "brown_dwarf_formed")
    _assert_conserved(a, b, o)


def test_a_brown_dwarf_and_a_giant_make_a_brown_dwarf_without_a_new_record():
    a, b = _giant("a", 40.0), _giant("b", 5.0)
    o = co.resolve(a, b, 0.0, 0.0)
    assert o.survivors[0].kind == "B" and o.records == ("giant_merger",)


def test_two_massive_brown_dwarfs_light_a_star_with_no_planets_and_a_merger_flag():
    a, b = _giant("a", 70.0), _giant("b", 70.0)
    o = co.resolve(a, b, 0.0, 0.0)
    (s,) = o.survivors
    assert s.kind == "S" and s.radius_m is None
    assert o.records == ("giant_merger", "star_created")
    assert {"merger remnant", "no planets"} <= set(o.flags)
    assert s.mass_kg / SUN >= stellar_abundances.hydrogen_burning_limit_sol(stellar_abundances.PROTO_Z)
    _assert_conserved(a, b, o)


def test_the_fusion_line_moves_with_metallicity():
    pair = lambda z: co.resolve(_giant("a", 40.0, z=z), _giant("b", 40.0, z=z), 0.0, 0.0)
    solar = pair(stellar_abundances.PROTO_Z)
    metal_free = pair(1e-6)
    assert solar.survivors[0].mass_kg == pytest.approx(metal_free.survivors[0].mass_kg, rel=0.2)
    limit = lambda z: stellar_abundances.hydrogen_burning_limit_sol(z) * SUN
    assert limit(1e-6) > limit(stellar_abundances.PROTO_Z)


def test_a_result_within_three_percent_of_the_hydrogen_limit_is_flagged_borderline_and_the_floor_note_is_added():
    limit_sol = stellar_abundances.hydrogen_burning_limit_sol(stellar_abundances.PROTO_Z)
    # Two equal brown dwarfs whose merged mass, after the gas lost, lands at 1.01 times the limit.
    each = 40.0
    for _ in range(6):
        o = co.resolve(_giant("a", each, kind="B"), _giant("b", each, kind="B"), 0.0, 0.0)
        each *= 1.01 * limit_sol * SUN / o.survivors[0].mass_kg
    assert o.survivors[0].mass_kg / SUN == pytest.approx(1.01 * limit_sol, rel=1e-3)
    assert "borderline" in o.flags and o.records[-1] == "star_created"
    assert o.survivors[0].mass_kg / SUN < tuning.IMF_BREAKS_SOL[0] and "below the generator's floor" in o.flags
    far = co.resolve(_giant("a", 80.0), _giant("b", 80.0), 0.0, 0.0)
    assert "borderline" not in far.flags


# --- Determinism, order and refusals ----------------------------------------------------------------------

def test_the_same_input_gives_the_same_outcome_in_either_order():
    a, b = _earth("a", 1.0, 5.0), _earth("b", 0.6, -12.0)
    assert co.resolve(a, b, 3.0e6, 77.0) == co.resolve(b, a, 3.0e6, 77.0) == co.resolve(a, b, 3.0e6, 77.0)
    assert co.resolve(a, b, 3.0e6, 77.0).seed == "collision:a:b:77.0"


def test_every_row_conserves_mass_and_momentum_over_a_sweep_of_speeds_angles_and_ratios():
    seen = set()
    for ratio in (1.0, 0.3, 0.05):
        for speed in (0.5, 8.0, 15.0, 25.0, 60.0, 150.0):
            for angle in (0.0, 0.3, 0.55, 0.8, 0.95):
                a = _earth("a")
                b = co.Body("b", "T", ratio * EARTH, EARTH_R * ratio ** (1 / 3), velocity_ms=(speed * KMS, 2.0e3, -1.0e3),
                            position_m=(1.0e7, 2.0e7, 3.0e7))
                o = co.resolve(a, b, angle * (a.radius_m + b.radius_m), 5.0)
                seen.add(o.code)
                _assert_conserved(a, b, o)
    assert {"merger", "graze_and_merge", "hit_and_run", "erosion", "disruption"} <= seen


def test_the_refusals():
    a, b = _earth("a"), _earth("b")
    with pytest.raises(ValueError):
        co.resolve(a, co.Body("a", "T", EARTH, EARTH_R), 0.0, 0.0)
    with pytest.raises(ValueError):
        co.resolve(a, co.Body("b", "Q", EARTH, EARTH_R), 0.0, 0.0)
    with pytest.raises(ValueError):
        co.resolve(a, co.Body("b", "T", -1.0, EARTH_R), 0.0, 0.0)
    with pytest.raises(ValueError):
        co.resolve(a, b, a.radius_m + b.radius_m + 1.0, 0.0)
    with pytest.raises(ValueError):
        co.resolve(a, b, -1.0, 0.0)
    with pytest.raises(ValueError):
        co.catastrophic_energy(1.0, 2.0)
