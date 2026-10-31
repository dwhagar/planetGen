"""
GEN.1088: `planetgen.physics.belt_origin`, where a belt came from. Each formation cause fires on its own input and on no other, each shaping cause
takes the mass the rule says and no more, mass is conserved, and every belt has a sentence made of its own numbers. The Solar System is a test
target here (a Jupiter-like giant, a snow line, a Neptune-like last giant), never an input of the module.
"""

import json
import math

import pytest

from planetgen import tuning
from planetgen.physics import belt_mass, belt_origin as bo
from planetgen.physics.constants import EARTH_MASS_TO_KG as EARTH, SOLAR_MASS_TO_KG as SUN

JUPITER = 1.898e27
STAR = SUN
SNOW = 2.7
GROWTH = tuning.PLANET_OUTER_GROWTH_AU


def rocky(a, mass=5.97e24, e=0.0):
    return bo.Body(a, mass, e)


def giant(a, mass=JUPITER, e=0.05, migrated=False, formed=None):
    return bo.Body(a, mass, e, True, migrated, formed)


def build(chain, **kwargs):
    arguments = dict(star_mass_kg=STAR, snow_line_au=SNOW, formation_snow_line_au=SNOW, growth_limit_au=GROWTH, disc_edge_au=100.0)
    arguments.update(kwargs)
    return bo.build(chain, **arguments)


def causes(belts):
    return [belt.origin for belt in belts]


SOL_LIKE = [rocky(0.39, 3.3e23), rocky(0.72, 4.87e24), rocky(1.0), rocky(1.52, 6.4e23), giant(5.2), giant(9.5, 5.7e26), giant(19.2, 8.7e25),
            giant(30.0, 1.0e26, 0.01)]


# --- the cause file -------------------------------------------------------------------------------------------------------------------------

def test_the_cause_file_has_six_formation_causes_and_two_shaping_ones_with_sentences_that_format():
    assert set(bo.FORMATION_CAUSES) == {"gap", "tail", "snow-line-bump", "outer-disc", "collisional-debris", "destroyed-planet"}
    assert set(bo.SHAPING_CAUSES) == {"giant-stirring", "instability-depletion"}
    for code, cause in bo.CAUSES.items():
        assert cause["label"] and cause["text"].format(**{name: 1.5 for name in cause["fields"]}), code


def _write(tmp_path, causes):
    path = tmp_path / "belt_origins.json"
    path.write_text(json.dumps({"causes": causes}), encoding="utf-8")
    return path


@pytest.mark.parametrize("cause,match", [
    ({"kind": "formation", "label": "x", "fields": ["a"]}, "needs a kind"),
    ({"kind": "formation", "label": "x", "fields": ["a"], "text": ""}, "needs a kind"),
    ({"kind": "origin", "label": "x", "fields": ["a"], "text": "{a}"}, "needs a kind"),
    ({"kind": "formation", "label": "x", "fields": ["a"], "text": "{b}"}, "does not match its fields"),
    ({"kind": "formation", "label": "x", "fields": ["a"], "text": "{a:.3g}", "co_text": "{b}", "co_fields": ["a"]}, "does not match its fields"),
])
def test_the_loader_refuses_a_cause_without_a_kind_a_sentence_or_the_fields_its_sentence_uses(tmp_path, cause, match):
    with pytest.raises(bo.BeltOriginError, match=match):
        bo._load(_write(tmp_path, {"x": cause}))


def test_the_loader_refuses_a_file_without_causes(tmp_path):
    with pytest.raises(bo.BeltOriginError):
        bo._load(_write(tmp_path, {}))


# --- gap ------------------------------------------------------------------------------------------------------------------------------------

def test_two_neighbours_far_apart_leave_a_gap_belt_between_their_keep_outs():
    chain = [rocky(0.4), rocky(3.0)]
    (belt,) = [b for b in build(chain) if b.origin == "gap"]
    near, far = bo.keep_out(chain[0], STAR), bo.keep_out(chain[1], STAR)
    assert belt.inner_au == pytest.approx(near[1]) and belt.outer_au == pytest.approx(far[0] - tuning.MIN_ASTEROID_BELT_SEPARATION)
    assert belt.primordial_mass_earth == pytest.approx(belt_mass.primordial_mass_earth(belt.inner_au, belt.outer_au, STAR, SNOW))
    assert "0.4 AU to 3 AU" in belt.sentence


def test_neighbours_at_the_spacing_of_the_chain_leave_no_gap():
    assert [b for b in build([rocky(1.0), rocky(1.22)]) if b.origin == "gap"] == []


def test_a_gap_is_found_where_the_spacing_just_passes_the_limit_and_not_below_it():
    def spacing(a2):
        from planetgen.physics import planet_architecture
        return planet_architecture.spacing_k(1.0, 5.97e24, a2, 5.97e24, STAR)
    below = next(a for a in [1.0 + i * 0.01 for i in range(1, 400)] if spacing(a) >= bo.GAP_SPACING_K - 0.5)
    assert not bo.gaps([rocky(1.0), rocky(below)], STAR)
    assert bo.gaps([rocky(1.0), rocky(below * 1.5)], STAR)


def test_a_gap_too_narrow_after_the_keep_outs_is_not_a_belt():
    heavy = 1.0e28
    assert bo.gaps([bo.Body(1.0, heavy), bo.Body(1.6, heavy)], STAR) == []


# --- tail -----------------------------------------------------------------------------------------------------------------------------------

def test_a_chain_of_rocky_planets_that_stops_short_has_an_unfinished_tail_to_the_growth_limit():
    belts = build([rocky(0.5), rocky(1.0)])
    (belt,) = [b for b in belts if b.origin == "tail"]
    assert belt.inner_au == pytest.approx(bo.keep_out(rocky(1.0), STAR)[1]) and belt.outer_au == GROWTH
    assert "past 1 AU" in belt.sentence


def test_a_chain_that_ends_in_a_giant_has_no_tail_and_one_that_ends_past_the_growth_limit_has_none():
    assert "tail" not in causes(build([rocky(1.0), giant(5.2)]))
    assert bo.tail([rocky(GROWTH * 1.5)], STAR, GROWTH) == []
    assert bo.tail([], STAR, GROWTH) == []


# --- snow line ------------------------------------------------------------------------------------------------------------------------------

def test_solids_at_a_free_snow_line_are_a_belt_of_their_own():
    belts = build([rocky(1.0), giant(5.2)])
    (belt,) = [b for b in belts if b.origin == "snow-line-bump"]
    assert belt.inner_au == pytest.approx(SNOW / (1 + bo.SNOW_BUMP_WIDTH)) and belt.outer_au == pytest.approx(SNOW * (1 + bo.SNOW_BUMP_WIDTH))
    assert belt.formation_snow_line_au == SNOW and "2.7 AU" in belt.sentence


def test_a_snow_line_inside_a_gap_or_tail_is_a_second_cause_not_a_second_belt():
    belts = build([rocky(1.0)])
    assert "snow-line-bump" not in causes(belts)
    (tail,) = [b for b in belts if b.origin == "tail"]
    assert tail.co_causes == ("snow-line-bump",) and "snow line of the young system" in tail.sentence
    gap = [b for b in build([rocky(0.4), rocky(6.0)]) if b.origin == "gap"][0]
    assert gap.co_causes == ("snow-line-bump",)


def test_a_snow_line_under_a_planet_s_keep_out_makes_no_belt():
    assert bo.snow_line_bump([rocky(SNOW)], STAR, SNOW, []) == []


def test_the_snow_line_belt_keeps_the_snow_line_of_the_formation_not_of_now():
    (belt,) = [b for b in build([rocky(1.0), giant(9.0)], snow_line_au=4.0, formation_snow_line_au=2.2) if b.origin == "snow-line-bump"]
    assert belt.formation_snow_line_au == 2.2 and belt.inner_au == pytest.approx(2.2 / (1 + bo.SNOW_BUMP_WIDTH))


# --- outer disc -----------------------------------------------------------------------------------------------------------------------------

def test_the_outer_disc_starts_past_the_chaotic_zone_of_the_last_giant():
    chain = [rocky(1.0), giant(5.2), giant(30.0, 1.0e26, 0.01)]
    (disc,) = [b for b in build(chain) if b.origin == "outer-disc"]
    half = bo.chaotic_half_width_au(chain[-1], STAR)
    assert disc.inner_au == pytest.approx(max(30.0 + half, bo.keep_out(chain[-1], STAR)[1])) and disc.inner_au > 30.0 + half
    assert disc.outer_au == 100.0
    assert disc.stirred and "past 34.1 AU" in disc.sentence


def test_with_no_giant_the_outer_disc_takes_over_where_the_tail_ends_and_a_disc_with_no_room_is_none():
    belts = build([rocky(1.0)])
    (disc,) = [b for b in belts if b.origin == "outer-disc"]
    (tail,) = [b for b in belts if b.origin == "tail"]
    assert disc.inner_au == pytest.approx(tail.outer_au) and disc.inner_au == GROWTH
    assert bo.outer_disc([rocky(1.0)], STAR, GROWTH, 12.0) == []


def test_a_hot_jupiter_with_no_other_body_leaves_the_disc_past_it_as_the_outer_disc_and_sweeps_the_part_it_crossed():
    chain = [giant(0.05, migrated=True, formed=3.0)]
    (disc,) = [b for b in build(chain) if b.origin == "outer-disc"]
    assert disc.inner_au == pytest.approx(max(bo.keep_out(chain[0], STAR)[1], 0.05 + bo.chaotic_half_width_au(chain[0], STAR)))
    assert disc.history[-1].cause == "instability-depletion" and 0.0 < disc.scattered_earth < disc.primordial_mass_earth


# --- destroyed planets and collision rubble -------------------------------------------------------------------------------------------------

def test_a_destroyed_planet_leaves_half_its_mass_as_a_belt_round_its_orbit():
    debris = bo.destroyed_planet(EARTH, 2.0, 0.02)
    assert debris.cause == "destroyed-planet" and debris.mass_earth == pytest.approx(bo.DESTROYED_RETAINED_SHARE)
    (belt,) = [b for b in build([giant(8.0)], debris=[debris]) if b.origin == "destroyed-planet"]
    assert belt.primordial_mass_earth == pytest.approx(0.5)
    assert belt.inner_au == pytest.approx(2.0 * (1.0 - bo.DEBRIS_MIN_WIDTH)) and belt.outer_au == pytest.approx(2.0 * (1.0 + bo.DEBRIS_MIN_WIDTH))
    assert "destroyed near 2 AU" in belt.sentence
    wide = bo.destroyed_planet(EARTH, 2.0, 0.2)
    assert bo._debris_candidate(wide)["outer_au"] == pytest.approx(2.0 * 1.2)


def test_collision_rubble_follows_the_chance_for_the_age_and_the_draw():
    assert bo.collisional_debris_chance(0.0) == 0.0
    assert bo.collisional_debris_chance(1e9) == pytest.approx(1 - math.exp(-bo.COLLISIONAL_DEBRIS_PER_GYR))
    assert bo.collisional_debris_chance(4e9) > bo.collisional_debris_chance(1e9)
    chance = bo.collisional_debris_chance(4.5e9)
    assert bo.collisional_debris(2.0, 0.0, 4.5e9, chance * 1.01) is None
    event = bo.collisional_debris(2.0, 0.0, 4.5e9, chance * 0.99)
    assert event.cause == "collisional-debris" and event.mass_earth == bo.COLLISIONAL_DEBRIS_MASS_EARTH
    (belt,) = [b for b in build([giant(8.0)], debris=[event]) if b.origin == "collisional-debris"]
    assert belt.primordial_mass_earth == pytest.approx(bo.COLLISIONAL_DEBRIS_MASS_EARTH) and "collision" in belt.sentence


def test_a_debris_event_of_another_name_is_refused_and_a_dust_belt_is_dropped():
    with pytest.raises(bo.BeltOriginError):
        build([giant(8.0)], debris=[bo.Debris("comet", 1.0, 2.0)])
    assert [b for b in build([giant(8.0)], debris=[bo.Debris("collisional-debris", 1e-9, 2.0)]) if b.origin == "collisional-debris"] == []


# --- giants clear and stir ------------------------------------------------------------------------------------------------------------------

def test_the_chaotic_zone_and_the_kirkwood_gaps_of_a_jupiter_are_where_the_physics_puts_them():
    jupiter = giant(5.2, 1.898e27)
    assert bo.chaotic_half_width_au(jupiter, SUN) == pytest.approx(0.93, abs=0.02)
    distances = {ratio: distance for ratio, distance in bo.resonance_distances_au(jupiter)}
    assert distances[(3, 1)] == pytest.approx(2.50, abs=0.01) and distances[(5, 2)] == pytest.approx(2.82, abs=0.01)
    assert distances[(7, 3)] == pytest.approx(2.96, abs=0.01) and distances[(2, 1)] == pytest.approx(3.28, abs=0.01)
    assert distances[(4, 1)] == pytest.approx(2.06, abs=0.01)


def test_a_belt_inside_the_chaotic_zone_of_a_giant_is_cleared_whole_and_a_belt_far_away_is_not_touched():
    near = bo.Body(5.2, JUPITER, 0.0, True)
    spans = bo.cleared_spans(near, SUN, 4.8, 5.6)
    assert spans[0] == (4.8, 5.6)
    belts = build([giant(5.2), rocky(1.0)])
    outer = bo.cleared_spans(near, SUN, 40.0, 60.0)
    assert outer == []
    assert not bo.is_stirred(near, SUN, 40.0, 60.0) and bo.is_stirred(near, SUN, 6.0, 8.0)
    assert belts


def test_the_resonance_gaps_take_a_share_of_a_belt_inside_the_giant_and_none_outside_it():
    chain = [rocky(1.0), giant(5.2)]
    (belt,) = [b for b in build(chain) if b.origin == "snow-line-bump"]
    assert 0.1 < belt.cleared_earth / belt.primordial_mass_earth < 0.6 and belt.stirred
    assert belt.history[0].cause == "giant-stirring" and "5.2 AU cleared" in belt.sentence
    beyond = bo.cleared_spans(chain[1], STAR, 8.0, 12.0)
    assert [span for span in beyond if span[0] > 6.5] == []


def test_a_heavier_giant_clears_more_of_the_same_belt():
    def share(mass):
        (belt,) = [b for b in build([rocky(1.0), giant(5.2, mass)]) if b.inner_au <= SNOW <= b.outer_au]
        return belt.cleared_earth / belt.primordial_mass_earth
    assert share(JUPITER) > share(0.3 * JUPITER) > share(0.1 * JUPITER)


def test_a_belt_a_giant_only_stirs_loses_nothing_and_says_so():
    chain = [rocky(1.0), giant(5.2), giant(30.0, 1.0e26, 0.01)]
    (disc,) = [b for b in build(chain) if b.origin == "outer-disc"]
    assert disc.cleared_earth == 0.0 and disc.stirred and disc.history[0].text == "The giant at 30 AU stirs it."
    assert disc.mass_earth == disc.primordial_mass_earth


# --- a moving giant sweeps ------------------------------------------------------------------------------------------------------------------

def test_the_survival_share_runs_log_uniformly_over_the_range():
    low, high = bo.INSTABILITY_SURVIVAL_RANGE
    assert bo.instability_survival(0.0) == pytest.approx(low) and bo.instability_survival(1.0) == pytest.approx(high)
    assert bo.instability_survival(0.5) == pytest.approx(math.sqrt(low * high))
    assert bo.instability_survival(-1.0) == pytest.approx(low) and bo.instability_survival(9.0) == pytest.approx(high)


def test_only_a_giant_that_moved_and_crossed_the_belt_sweeps_it():
    mover = giant(0.05, migrated=True, formed=3.0)
    assert bo.moved_through(mover, 1.0, 2.0, STAR)
    assert not bo.moved_through(mover, 20.0, 30.0, STAR)
    assert not bo.moved_through(giant(0.05, formed=3.0), 1.0, 2.0, STAR), "it did not migrate"
    assert not bo.moved_through(rocky(0.05), 1.0, 2.0, STAR)
    assert not bo.moved_through(giant(0.05, migrated=True), 1.0, 2.0, STAR), "no record of where it formed"


def test_a_belt_a_moving_giant_crossed_loses_the_swept_part_to_the_scatter_and_keeps_the_rest():
    chain = [giant(0.05, migrated=True, formed=3.0), rocky(1.0)]
    (belt,) = [b for b in build(chain, u=0.0) if b.origin == "tail"]
    survival = bo.instability_survival(0.0)
    assert belt.scattered_earth > 0.0 and belt.history[-1].cause == "instability-depletion"
    swept = belt.scattered_earth / (1.0 - survival)
    assert 0.0 < swept < belt.primordial_mass_earth, "only the crossed stretch is swept"
    (kept,) = [b for b in build(chain, u=1.0) if b.origin == "tail"]
    assert kept.scattered_earth < belt.scattered_earth and "swept away" in belt.sentence
    (far,) = [b for b in build(chain, u=0.0) if b.origin == "outer-disc"]
    assert far.scattered_earth == 0.0


# --- mass is conserved ----------------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("chain", [
    SOL_LIKE,
    [rocky(0.4), rocky(3.0)],
    [giant(0.05, migrated=True, formed=3.0), rocky(1.0), giant(6.0, 3 * JUPITER)],
    [rocky(0.1), rocky(0.2), giant(1.0, 0.3 * JUPITER, 0.3), giant(4.0, 0.1 * JUPITER)],
])
@pytest.mark.parametrize("u", [0.0, 0.5, 1.0])
def test_every_belt_conserves_its_mass_and_has_a_sentence_made_of_its_numbers(chain, u):
    debris = [bo.destroyed_planet(2e24, 2.2, 0.1)]
    belts = build(chain, u=u, debris=debris)
    assert belts and [b.inner_au for b in belts] == sorted(b.inner_au for b in belts)
    for belt in belts:
        assert bo.conservation_error(belt) < 1e-9
        assert belt.mass_earth >= -1e-12 and belt.cleared_earth >= 0.0 and belt.scattered_earth >= 0.0
        assert belt.origin in bo.FORMATION_CAUSES and all(event.cause in bo.SHAPING_CAUSES for event in belt.history)
        assert "{" not in belt.sentence and belt.sentence.endswith(".") and belt.inner_au < belt.outer_au


def test_the_cleared_and_scattered_masses_are_the_sum_of_the_history():
    belts = build([giant(0.05, migrated=True, formed=3.0), rocky(1.0), giant(6.0)], u=0.3)
    for belt in belts:
        assert belt.cleared_earth + belt.scattered_earth == pytest.approx(sum(event.mass_earth for event in belt.history))


# --- the Solar System as a target -----------------------------------------------------------------------------------------------------------

def test_a_sol_like_chain_gets_a_main_belt_between_the_rocky_planets_and_jupiter_and_a_disc_past_the_last_giant():
    belts = build(SOL_LIKE)
    by_cause = {belt.origin: belt for belt in belts}
    main = by_cause["snow-line-bump"]
    assert 2.0 < main.inner_au < 2.7 < main.outer_au < 3.5 and main.stirred
    assert by_cause["outer-disc"].inner_au > 30.0
    assert 0.0 < main.cleared_earth < main.primordial_mass_earth


def test_the_same_input_gives_the_same_belts():
    assert build(SOL_LIKE, u=0.37) == build(SOL_LIKE, u=0.37)
    assert build(list(reversed(SOL_LIKE))) == build(SOL_LIKE)


def test_no_planets_at_all_still_gives_the_disc_a_belt_with_its_cause():
    belts = build([])
    assert [b.origin for b in belts] == ["snow-line-bump", "outer-disc"] or [b.origin for b in belts] == ["outer-disc"]
    assert all(b.sentence for b in belts)
