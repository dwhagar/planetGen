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
    pytest.fail("(TEST.737) no world with air in these seeds")


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


def test_the_stored_molar_mass_density_and_scale_height_follow_the_air_the_trace_gases_changed(world):
    """TEST.585: the trace gases move the mean molar mass at the eighth digit; the stored molar mass, density and scale height are those of the new air."""
    from planetgen.physics import atmosphere_state, constants
    drawn = None
    for seed in range(60):
        other = copy.deepcopy(world)
        other.mass = world.mass * (1.0 + seed * 1e-3)
        if trace_air.apply(other):
            drawn = other
            break
    assert drawn is not None
    state = atmosphere_state.air_state(atmosphere.gas_partials_pa(drawn), drawn.surface_temperature, drawn.gravity * constants.EARTH_GRAVITY)
    assert drawn.atm_molar_density == pytest.approx(state.mean_molar_mass_kg_mol, rel=1e-12)
    assert drawn.atm_density == pytest.approx(state.density_kg_m3, rel=1e-12)
    assert drawn.scale_height == pytest.approx(state.scale_height_m / 1000.0, rel=1e-12)


def test_a_hot_air_keeps_its_equilibrium_nitrogen_dioxide_and_nitrous_oxide():
    """GEN.1614: NO2 and N2O are gases of a hot air at their equilibrium amounts (far over the design brackets of a draw); `trace_air` takes off an air only what a draw can have given
    it (the top of the gas's bracket), so the equilibrium amount stays and the total pressure stays."""
    total = 4.0e7
    air = {"N2": 0.618 * total, "CO2": 0.38 * total, "NO2": 1.5e-3 * total, "N2O": 5.0e-4 * total}
    kept = trace_air._strip(air)
    assert sum(kept.values()) == pytest.approx(total, rel=1e-12)
    for formula, fraction in (("NO2", 1.5e-3), ("N2O", 5.0e-4)):
        assert kept[formula] / total == pytest.approx(fraction, abs=1.01 * trace_air.TRACE_CEILING[formula])      # off by the most a draw could have added, no more
        assert kept[formula] / total > 0.99 * fraction


def test_a_drawn_trace_gas_is_taken_off_whole():
    """What a draw gave an air is at most the top of its bracket, so `_strip` removes it whole and a second `apply` gives the same air (the module's rule)."""
    total = 1.0e5
    drawn = {"N2": 78000.0, "O2": 20000.0, "N2O": 0.9 * trace_air.TRACE_CEILING["N2O"] * total, "PH3": 1.0e-9 * total, "NO2": 5.0e-12 * total}
    kept = trace_air._strip(drawn)
    assert set(kept) == {"N2", "O2"}
    assert sum(kept.values()) == pytest.approx(sum(drawn.values()), rel=1e-12)


def test_the_ceiling_is_the_top_of_each_gas_s_design_bracket():
    for row in mixture_trace.WIRED_ROWS:
        for formula, spec in mixture_trace.ABUNDANCES[row].items():
            assert trace_air.TRACE_CEILING[formula] == spec[1]
    assert set(trace_air.TRACE_CEILING) == set(trace_air.TRACE_GASES)
