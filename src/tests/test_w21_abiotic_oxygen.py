# tests/test_w21_abiotic_oxygen.py

"""
TEST.491, the W21 audit of docs/design/systems/weird-worlds.md: an M dwarf world that lost its water leaves an oxygen air with no life, and the life
scores never read oxygen alone (Luger and Barnes 2015; Wordsworth and Pierrehumbert 2014).

Two links, each end to end over a sample: the air draw (`mixture_draw.draw_mixture`, the `j5` abiotic-oxygen row, the escape filter and the abiotic cap)
and the life scoring of a finished body (`life_classes.score_classes`). A whole generated system needs the stellar spectra dataset, which a test machine
does not always hold, so the chain is tested at these two joints; `generation/life._apply` stores the scores of a world with no class supported as the
classes its conditions support and no life.
"""

import copy


from planetgen.generation.possibility_table import Subject
from planetgen.physics import atmosphere, life_classes, mixture_draw as md
from planetgen.util import draw
from tests.test_atmosphere_mixture_evolved import ABIOTIC, BASE
from tests.test_life_classes import EARTH

SEEDS = range(240)
M_DWARF_ECOSPHERE = [{**ABIOTIC, "age_gy": age, "mass_earth": mass} for age in (1.5, 5.0, 9.0) for mass in (0.4, 1.0, 3.0)]


def _draws(changes, seeds=SEEDS):
    for seed in seeds:
        with draw.bound(seed):
            yield md.draw_mixture(Subject("terrestrial", {**BASE, **changes}))


def test_m_dwarf_ecosphere_worlds_can_hold_an_oxygen_air_and_it_is_the_abiotic_family_that_does():
    rich, families = [], set()
    for changes in M_DWARF_ECOSPHERE:
        for result in _draws(changes, range(80)):
            families.add(result.family)
            if result.mole_fractions.get("O2", 0.0) > 0.1:
                rich.append(result)
    assert rich, "no M dwarf ecosphere world kept an oxygen air"
    assert {result.family for result in rich} == {"abiotic-oxygen"}      # oxygen above a tenth comes from the water-loss family alone
    assert max(result.mole_fractions["O2"] for result in rich) > 0.5      # the several-oceans tail of Luger and Barnes
    assert all(result.pressure_bar > 0.0 for result in rich)


def test_no_other_host_gets_the_abiotic_oxygen_family_and_no_other_family_holds_more_than_a_trace_of_oxygen():
    for result in _draws({"host_type": "G", "age_gy": 5.0}, range(120)):
        assert result.family != "abiotic-oxygen"
    for result in _draws({"host_type": "K", "age_gy": 5.0}, range(120)):
        assert result.family != "abiotic-oxygen"
    for changes in M_DWARF_ECOSPHERE[:3]:
        for result in _draws(changes, range(120)):
            if result.family != "abiotic-oxygen":
                assert result.mole_fractions.get("O2", 0.0) < 0.05, result.family


def _gas(**partial_pa):
    return {atmosphere.GAS_FIELD: dict(partial_pa)}


def _supported(values):
    return {class_id for class_id, score in life_classes.score_classes(values).items() if score.score > 0}


def _dead_world(**over):
    """The Earth's mixture of gases on ground too hot for any solvent: no class can live there."""
    values = copy.deepcopy(EARTH)
    values.update(surface_temperature_k=470.0, hydrosphere="none", ocean_class=None, water_activity=0.0, water_mass_fraction=0.0)
    values.update(over)
    return values


def test_oxygen_alone_does_not_make_any_class_live_on_a_world_that_cannot_hold_life():
    oxygen = _dead_world(**_gas(O2=95000.0, N2=500.0, CO2=40.0, Ar=934.0))
    assert _supported(oxygen) == set()
    cold = copy.deepcopy(oxygen)
    cold.update(surface_temperature_k=60.0)
    assert _supported(cold) == set()
    assert life_classes.dominant(life_classes.score_classes(oxygen), draw) is None


def test_putting_oxygen_in_the_air_adds_only_the_classes_that_breathe_it():
    nitrogen = copy.deepcopy(EARTH)
    nitrogen.update(**_gas(N2=99300.0, CO2=40.0, Ar=934.0, H2O=1000.0))
    oxygen = copy.deepcopy(EARTH)
    oxygen.update(**_gas(O2=95000.0, N2=500.0, CO2=40.0, Ar=934.0, H2O=1000.0))
    gained = _supported(oxygen) - _supported(nitrogen)
    breathers = {row.id for row in life_classes.CLASSES if row.life["oxidant"] == "o2"}
    assert gained <= breathers
    assert _supported(nitrogen) - _supported(oxygen), "a pure-oxygen air should cost the anaerobes something"
    # and what the oxygen adds needs the water, warmth and light the other classes need too
    assert _supported(_dead_world(**_gas(O2=95000.0, N2=500.0, CO2=40.0, Ar=934.0))) == set()


def test_a_world_with_an_oxygen_air_and_no_class_is_scored_with_no_life():
    from planetgen.physics import life_scores
    assert life_scores.scores(0.0, 0.0, None) == (0, 0, 0, 0)
    assert life_scores.scores(0.0, 0.0, "multicellularity", grade=3) == (0, 0, 0, 0)
