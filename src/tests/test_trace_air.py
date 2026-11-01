# tests/test_trace_air.py

"""GEN.524 part 1: `physics.trace_air`, the trace gases of table m drawn into a finished air after the life stage."""

import copy

import pytest

from planetgen.generation import life
from planetgen.physics import atmosphere, biogenic_oxygen as bo, mixture_trace, trace_air
from planetgen.physics.mixture_base import MixtureError
from planetgen.util import draw

from tests.test_life_air import _worlds


@pytest.fixture(scope="module")
def airy_world():
    """A generated rocky world with a nitrogen, oxygen and carbon dioxide air of its own making stored on it."""
    for body in _worlds(range(0, 12)):
        if (body.atmospheric_pressure or 0) > 20000.0 and body.surface_temperature and body.greenhouse_factor is not None:
            body = copy.deepcopy(body)
            bo.store_air(body, {"N2": 78000.0, "O2": 20000.0, "CO2": 300.0, "Ar": 900.0})
            return body
    pytest.skip("no world with air in these seeds")


@pytest.fixture
def world(airy_world):
    body = copy.deepcopy(airy_world)
    body.life_stage = "multicellularity"
    return body


def _drawn(body, seeds=range(60)):
    """The gases of `TRACE_GASES` the air of `body` holds, over several draws on different bodies (the stream is the body's own)."""
    seen = set()
    for seed in seeds:
        other = copy.deepcopy(body)
        other.mass = body.mass * (1.0 + seed * 1e-3)
        seen |= set(trace_air.apply(other))
    return seen


def test_the_rows_are_the_design_flag_rows_and_hydrogen_chloride_and_fluoride_are_not_among_them():
    assert set(mixture_trace.WIRED_ROWS) == {"m3", "m4", "m11"} == set(mixture_trace.ABUNDANCES)
    assert set(trace_air.TRACE_GASES) == {"N2O", "PH3", "NO2", "CH3Cl"}
    assert not {"HCl", "HF"} & set(trace_air.TRACE_GASES)


def test_a_living_world_gets_trace_gases_and_the_air_keeps_its_total_pressure(world):
    before = atmosphere.gas_partials_pa(world)
    seen = _drawn(world)
    assert {"N2O", "PH3"} & seen
    assert not {"HCl", "HF"} & seen
    drawn = None
    for seed in range(60):
        other = copy.deepcopy(world)
        other.mass = world.mass * (1.0 + seed * 1e-3)
        if trace_air.apply(other):
            drawn = other
            break
    assert drawn is not None
    after = atmosphere.gas_partials_pa(drawn)
    assert sum(after.values()) == pytest.approx(sum(before.values()), rel=1e-12)
    assert sum(after.get(gas, 0.0) for gas in trace_air.TRACE_GASES) / sum(after.values()) < 1e-5
    assert all(after.get(gas, 0.0) > 0.0 for gas in trace_air.apply(copy.deepcopy(drawn)))


def test_the_amounts_sit_in_the_design_brackets(world):
    for seed in range(60):
        other = copy.deepcopy(world)
        other.mass = world.mass * (1.0 + seed * 1e-3)
        added = trace_air.apply(other)
        air = atmosphere.gas_partials_pa(other)
        total = sum(air.values())
        for formula in added:
            low, high, _ = next(spec for row in mixture_trace.ABUNDANCES.values() for name, spec in row.items() if name == formula)
            assert low * 0.999 <= air[formula] / total <= high * 1.001


def test_applying_twice_gives_the_same_air_and_a_lifeless_world_loses_them(world):
    for seed in range(60):
        other = copy.deepcopy(world)
        other.mass = world.mass * (1.0 + seed * 1e-3)
        if not trace_air.apply(other):
            continue
        once = atmosphere.gas_partials_pa(other)
        trace_air.apply(other)
        twice = atmosphere.gas_partials_pa(other)
        assert twice.keys() == once.keys()
        assert all(twice[gas] == pytest.approx(once[gas], rel=1e-9) for gas in once)
        other.life_stage = None
        assert trace_air.apply(other) == ()
        assert not set(trace_air.TRACE_GASES) & set(atmosphere.gas_partials_pa(other))
        return
    pytest.fail("no draw gave a trace gas")


def test_a_world_without_life_gets_none_because_the_design_flags_are_off(world):
    world.life_stage = None
    assert _drawn(world) == set()


def test_the_draw_is_made_on_the_body_s_own_stream_and_does_not_move_the_run(world):
    with draw.bound("run"):
        first = draw.random()
    with draw.bound("run"):
        trace_air.apply(copy.deepcopy(world))
        second = draw.random()
    assert first == second


def test_an_airless_body_and_a_giant_are_left_alone(world):
    airless = copy.deepcopy(world)
    airless.gas_partials_pa = {}
    assert trace_air.apply(airless) == ()
    giant = copy.deepcopy(world)
    giant.body_type = "g"
    assert trace_air.apply(giant) == ()


def test_the_subject_holds_every_quantity_table_m_reads(world):
    partials = atmosphere.gas_partials_pa(world)
    subject = trace_air.subject_of(world, partials)
    needed = set(mixture_trace.required_trace_inputs("terrestrial"))
    assert needed <= set(subject.quantities)
    assert subject.quantities["life_chemical"] is True
    assert subject.quantities["lightning"] is False


def test_apply_life_data_stores_the_trace_gases_after_the_life_stage(world):
    life.apply_life_data(world, keep_stored=True)
    air = atmosphere.gas_partials_pa(world)
    assert sum(air.values()) > 0.0
    assert not {"HCl", "HF"} & set(air)


def test_a_subject_missing_a_quantity_is_refused_by_name():
    with pytest.raises(MixtureError):
        mixture_trace.draw_trace_amounts(mixture_trace.Subject("terrestrial", {}))
