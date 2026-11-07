# tests/test_adaptation_traits.py

"""GEN.876: the adaptation traits of a body's life (`physics/adaptation_traits.py`, `physics/adaptation_draw.py`): the text, the grades and their factors, the
draw (a chance that never falls with the stress and never passes its cap, six draws always), the limits they widen (and the human scores they never touch),
and the life pass that stores them."""

import copy
import math

import pytest

from planetgen.generation import life
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import adaptation_draw, adaptation_traits as traits, habitability as hab, life_classes
from planetgen.physics.habitability import World
from planetgen.util import draw


def test_the_text_of_a_trait_set_reads_back_and_reads_as_words():
    found = {"sulfide": 3, "radiation": 2}
    text = traits.text_of(found)
    assert text == "radiation:2,sulfide:3", "the order of the stresses, not of the dict"
    assert traits.parse(text) == found
    assert traits.sentence(text) == "radiation (tolerant), sulfide (extreme)"
    assert traits.text_of({}) is None and traits.text_of({"radiation": 1}) is None
    assert traits.parse(None) == {} and traits.parse("") == {} and traits.sentence(None) is None


@pytest.mark.parametrize("text", ["radiation:1", "radiation:4", "uv:2", "radiation:2,radiation:3", "radiation"])
def test_a_text_that_is_not_a_trait_set_is_refused(text):
    with pytest.raises(ValueError):
        traits.parse(text)


def test_a_grade_never_loosens_a_limit_and_an_animal_never_beats_a_microbe():
    for stress in traits.STRESSES:
        animal = [traits.factor({stress: grade}, stress, True) for grade in traits.GRADES]
        assert animal[0] == 1.0 and animal == sorted(animal)
    micro = [traits.factor({"radiation": grade}, "radiation", False) for grade in traits.GRADES]
    assert micro == sorted(micro) and micro[0] == 1.0
    assert traits.factor({"radiation": 3}, "radiation", True) <= traits.factor({"radiation": 3}, "radiation", False)
    assert traits.factor({"radiation": 3}, "radiation", True) == traits.factor({"radiation": 2}, "radiation", True), "animals stop at tolerant"
    assert traits.factor({"sulfide": 3}, "sulfide", False) == 1.0, "no microbe limit on sulfide exists to widen"
    assert traits.factor(None, "radiation", True) == 1.0


def test_the_chance_of_a_trait_is_zero_at_no_stress_never_falls_and_never_passes_its_cap():
    for stress, onset in adaptation_draw.ONSET.items():
        assert adaptation_draw.chance(stress, 0.0) == 0.0 and adaptation_draw.chance(stress, onset) == 0.0
        amounts = [onset * 10 ** (step / 10) for step in range(0, 80)]
        chances = [adaptation_draw.chance(stress, amount) for amount in amounts]
        assert chances == sorted(chances) and max(chances) <= traits.P_CAP
        assert adaptation_draw.chance(stress, adaptation_draw.FULL[stress]) == pytest.approx(traits.P_CAP)


def test_a_body_with_no_stress_draws_no_trait_whatever_the_seed():
    calm = {"surface_dose_msv_yr": 1.0, "hydrosphere": "surface ocean", "atmospheric_pressure_pa": 1.0e5,
            "gas_partial_pa": {"o2": 2.1e4, "n2": 7.8e4, "co2": 40.0}}
    for seed in range(300):
        assert adaptation_draw.draw(calm, draw.Stream(seed)) == {}


def test_a_body_draws_six_uniforms_whatever_it_holds():
    """The draw count is fixed, so changing a trait's numbers moves no other draw of any body."""
    harsh = {"surface_dose_msv_yr": 1.0e9, "gas_partial_pa": {"co2": 9.0e4, "h2s": 5.0e3}, "atmospheric_pressure_pa": 1.0e5}
    for values in ({}, harsh):
        used, reference = draw.Stream("x"), draw.Stream("x")
        adaptation_draw.draw(values, used)
        for _ in range(2 * len(traits.STRESSES)):
            reference.random()
        assert used.random() == reference.random()


def test_the_harshest_world_draws_the_extreme_grade_more_often_than_a_mild_one():
    mild = {"surface_dose_msv_yr": 3.0e5}
    harsh = {"surface_dose_msv_yr": 1.0e9}
    def share(values, grade):
        found = [adaptation_draw.draw(values, draw.Stream(seed)).get("radiation") for seed in range(600)]
        return sum(1 for g in found if g == grade) / len(found)
    assert share(harsh, 3) > share(mild, 3)
    assert share(harsh, 2) + share(harsh, 3) > share(mild, 2) + share(mild, 3)
    assert share(harsh, 2) + share(harsh, 3) <= traits.P_CAP + 0.06


def _world(**kwargs):
    base = dict(pressure_kpa=101.3, temperature_c=15, gases={"O2": 0.21, "N2": 0.78})
    base.update(kwargs)
    return World(**base)


STRESSED = {
    "radiation": dict(surface_dose_sv_yr=80.0 * hab.SV_PER_GY),
    "sulfide": dict(gases={"O2": 0.21, "N2": 0.77, "H2S": 0.0004}),
    "carbon_dioxide": dict(gases={"O2": 0.21, "N2": 0.70, "CO2": 0.04}),
}


@pytest.mark.parametrize("stress", traits.STRESSES)
def test_a_trait_widens_its_own_limit_and_lowers_no_life_score(stress):
    plain = _world(**STRESSED[stress])
    tough = _world(**STRESSED[stress], traits={stress: 3})
    for score in (hab.l_rad, hab.l_rad_animal, hab.toxin_factor, hab.t_metazoa, hab.phi_bio, hab.phi_cpx):
        assert score(tough) >= score(plain) - 1e-12, score.__name__
    widened = {"radiation": hab.l_rad_animal, "sulfide": hab.toxin_factor, "carbon_dioxide": hab.toxin_factor}[stress]
    assert widened(tough) > widened(plain) + 0.05, "the stress is chosen so that the ordinary limit costs the life something"


@pytest.mark.parametrize("stress", traits.STRESSES)
def test_a_trait_never_changes_a_human_score(stress):
    plain = _world(**STRESSED[stress])
    tough = _world(**STRESSED[stress], traits={stress: 3})
    assert hab.phi4(plain) == hab.phi4(tough)


def test_a_trait_moves_only_the_limit_in_its_row():
    plain = _world(**STRESSED["radiation"])
    tough = _world(**STRESSED["radiation"], traits={"radiation": 3})
    assert (hab.l_solv(plain), hab.l_chem(plain), hab.l_ener(plain)) == (hab.l_solv(tough), hab.l_chem(tough), hab.l_ener(tough))
    assert hab.toxin_factor(plain) == hab.toxin_factor(tough)
    sulfur = _world(**STRESSED["sulfide"], traits={"sulfide": 2})
    assert hab.l_rad(sulfur) == hab.l_rad(_world(**STRESSED["sulfide"]))


def test_a_microbe_stands_a_dose_that_ends_an_animal():
    world = _world(surface_dose_sv_yr=300.0 * hab.SV_PER_GY)
    assert hab.l_rad(world) == 1.0 and hab.l_rad_animal(world) < 0.02
    extreme = _world(surface_dose_sv_yr=300.0 * hab.SV_PER_GY, traits={"radiation": 3})
    assert hab.l_rad_animal(extreme) == hab.l_rad_animal(_world(surface_dose_sv_yr=300.0 * hab.SV_PER_GY, traits={"radiation": 2})) > 0.9
    assert hab.animal_dose_share(_world(surface_dose_sv_yr=1.0e9)) == 0.0


# --- the life pass --------------------------------------------------------------------------------------------

SEEDS = range(5000, 5012)


def _bodies(system):
    for planet in system.all_planets:
        if planet.body_type == "a":
            continue
        yield from (planet, *planet.moons)


@pytest.fixture(scope="module")
def bodies():
    config = SystemConfig()
    config.PLANETS = True
    found = []
    for seed in SEEDS:
        with draw.bound(seed):
            found.extend(_bodies(StarSystem(system_config=config)))
    return found


def test_a_stored_trait_set_is_one_the_text_reads_and_only_a_living_body_has_it(bodies):
    for body in bodies:
        if body.life_class is None:
            assert body.life_traits is None
        elif body.life_traits is not None:
            assert traits.text_of(traits.parse(body.life_traits)) == body.life_traits


def _living(bodies):
    found = [body for body in bodies if body.life_class and body.body_type != "g"]
    if not found:
        pytest.skip("no body with life in these seeds")
    return found[0]


def test_the_life_pass_draws_the_traits_before_the_classes_are_scored(bodies, monkeypatch):
    body = copy.deepcopy(_living(bodies))
    seen = []
    real = life_classes.score_classes
    monkeypatch.setattr(adaptation_draw, "draw", lambda values, stream: {"radiation": 3})
    monkeypatch.setattr(life_classes, "score_classes",
                        lambda planet: seen.append(getattr(planet, "life_traits", "values")) or real(planet))
    life.apply_life_data(body, keep_stored=True)
    assert seen and seen[0] == "radiation:3"
    assert body.life_traits == "radiation:3"


def test_a_world_the_pass_leaves_lifeless_stores_no_trait(bodies, monkeypatch):
    body = copy.deepcopy(_living(bodies))
    monkeypatch.setattr(adaptation_draw, "draw", lambda values, stream: {"sulfide": 2})
    monkeypatch.setattr(life_classes, "dominant", lambda scores: None, raising=False)
    monkeypatch.setattr(life, "_choose", lambda planet, scores, keep: (None, False))
    life.apply_life_data(body, keep_stored=True)
    assert body.life_class is None and body.life_traits is None


def test_the_traits_of_a_body_come_from_a_stream_of_its_own(bodies):
    body = copy.deepcopy(_living(bodies))
    life.apply_life_data(body, keep_stored=True)
    first = body.life_traits
    run = draw.Stream(1)
    with draw.bound(run):
        life.apply_life_data(body, keep_stored=True)
    assert body.life_traits == first, "the same body, the same conditions, the same traits"
    assert run.random() == draw.Stream(1).random(), "an edit makes no draw of the run's"
