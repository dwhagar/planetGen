# tests/test_giant_structure.py

"""
GEN.733: the hydrostatic structure of a giant (`planetgen.physics.giant_structure`) on the fixture tables of
`tests/datasets_support.py` (an ideal gas held below a density cap, so the numbers are the fixture's; the physics of the real tables
is in `tests/test_datasets_real.py`): the fit closes the model on the radius it is given, pressure and mass inside fall and rise
together, the readers agree with each other, a denser composition needs a hotter adiabat, a radius no adiabat closes is said so,
and the missing tables are refused.
"""

import math

import pytest

from planetgen.datasets import fetch, hhe_eos
from planetgen.physics import constants, giant_structure as gs
from tests import datasets_support as support

EARTH = constants.EARTH_MASS_TO_KG
SATURN = dict(mass_kg=95.2 * EARTH, radius_m=58232e3, helium_share=0.27, envelope_z=0.04, heavy_mass_kg=20.0 * EARTH,
              temperature_1bar_k=134.0)
NEPTUNE = dict(mass_kg=17.1 * EARTH, radius_m=24764e3, helium_share=0.27, envelope_z=0.1, heavy_mass_kg=14.6 * EARTH,
               temperature_1bar_k=72.0, envelope_ice_fraction=0.5)


@pytest.fixture(autouse=True)
def tables(tmp_path, monkeypatch):
    support.reset_caches(tmp_path / "cache", monkeypatch)
    return support.install_hhe_eos(monkeypatch)


def fit(**changes):
    return gs.fit_adiabat(gs.StructureSubject(**{**SATURN, **changes}))


@pytest.mark.parametrize("spec", [SATURN, NEPTUNE], ids=["saturn", "neptune"])
def test_the_fit_closes_the_model_on_the_radius_it_is_given(spec):
    profile = gs.fit_adiabat(gs.StructureSubject(**spec))
    assert profile.closed and profile.radius_m == spec["radius_m"] and profile.mass_kg == spec["mass_kg"]
    assert gs.ADIABAT_RANGE_K[0] < profile.adiabat_1bar_k < gs.ADIABAT_RANGE_K[1]
    assert 3 <= profile.evaluations <= gs.MAX_EVALUATIONS


def test_pressure_rises_and_the_mass_and_radius_inside_fall_with_depth():
    profile = fit()
    pressures = [10 ** (5.0 + 0.25 * i) for i in range(1, 30)]
    pressures = [p for p in pressures if p < profile.central_pressure_pa]
    states = [profile.at_pressure(p) for p in pressures]
    assert all(a[0] > b[0] and a[1] > b[1] for a, b in zip(states, states[1:]))
    assert profile.at_pressure(gs.ONE_BAR_PA) == (1.0, 1.0) and profile.at_pressure(10 * profile.central_pressure_pa) == (0.0, 0.0)
    assert profile.central_pressure_pa > 3.0 / (8.0 * math.pi) * constants.G * SATURN["mass_kg"] ** 2 / SATURN["radius_m"] ** 4


def test_the_readers_are_inverses_of_each_other():
    profile = fit()
    for fraction in (0.9, 0.5, 0.1):
        radius, pressure = profile.at_mass_inside(fraction)
        assert profile.at_pressure(pressure)[1] == pytest.approx(fraction, rel=2e-3)
        assert profile.at_radius(radius)[0] == pytest.approx(fraction, rel=2e-3)
        assert profile.at_radius(radius)[1] == pytest.approx(pressure, rel=2e-3)
        assert profile.pressure_at_mass_inside_pa(fraction) == pressure
    assert profile.at_mass_inside(1.0) == (1.0, 0.0) and profile.at_mass_inside(0.0) == (0.0, profile.central_pressure_pa)
    assert profile.at_radius(1.0)[0] == 1.0 and profile.at_radius(0.0)[0] == 0.0


def test_a_denser_composition_needs_a_hotter_adiabat_to_fill_the_same_radius():
    base = fit().adiabat_1bar_k
    assert fit(helium_share=0.35).adiabat_1bar_k > base
    assert fit(envelope_z=0.12, heavy_mass_kg=25.0 * EARTH).adiabat_1bar_k > base
    assert fit(envelope_ice_fraction=0.8).adiabat_1bar_k < fit(envelope_ice_fraction=0.0).adiabat_1bar_k    # ice is less dense than rock and iron


def test_explicit_layers_and_the_heavy_mass_describe_the_same_body():
    core = (SATURN["heavy_mass_kg"] - SATURN["envelope_z"] * SATURN["mass_kg"]) / (1.0 - SATURN["envelope_z"])
    by_mass = fit()
    by_layers = fit(layers=(gs.Layer("core", core), gs.Layer("envelope")), heavy_mass_kg=0.0)
    assert by_layers.adiabat_1bar_k == pytest.approx(by_mass.adiabat_1bar_k, rel=1e-6)
    assert by_layers.central_pressure_pa == pytest.approx(by_mass.central_pressure_pa, rel=1e-6)


def test_a_mantle_between_the_core_and_the_envelope_is_read():
    layers = (gs.Layer("core", 3.0 * EARTH), gs.Layer("mantle", 6.0 * EARTH), gs.Layer("envelope"))
    profile = gs.fit_adiabat(gs.StructureSubject(**{**NEPTUNE, "layers": layers, "heavy_mass_kg": 0.0}))
    assert profile.closed
    inner = profile.at_mass_inside(3.0 * EARTH / NEPTUNE["mass_kg"])[1]
    outer = profile.at_mass_inside(9.0 * EARTH / NEPTUNE["mass_kg"])[1]
    assert gs.ONE_BAR_PA < outer < inner < profile.central_pressure_pa


def test_a_radius_no_adiabat_closes_is_taken_at_the_nearest_one_and_said_so():
    radius = SATURN["radius_m"]
    puffed = fit(radius_m=100 * radius)          # too light for its radius even on the hottest adiabat
    assert not puffed.closed and puffed.adiabat_1bar_k == pytest.approx(gs.ADIABAT_RANGE_K[1], rel=0.1)
    assert radius < puffed.radius_m < 100 * radius
    squeezed = fit(radius_m=0.1 * radius)        # too dense for its radius even on the coldest
    assert not squeezed.closed and squeezed.adiabat_1bar_k < 2 * gs.ADIABAT_RANGE_K[0]
    assert 0.1 * radius < squeezed.radius_m < radius


def test_a_body_with_no_heavy_core_is_one_envelope():
    profile = fit(heavy_mass_kg=SATURN["envelope_z"] * SATURN["mass_kg"])
    assert profile.closed
    assert profile.adiabat_1bar_k < fit().adiabat_1bar_k           # no dense core: a colder adiabat fills the radius


def test_nonsense_is_refused():
    with pytest.raises(gs.StructureError):
        gs.StructureSubject(**{**SATURN, "mass_kg": -1.0})
    with pytest.raises(gs.StructureError):
        gs.StructureSubject(**{**SATURN, "heavy_mass_kg": 2 * SATURN["mass_kg"]})
    with pytest.raises(gs.StructureError):
        gs.StructureSubject(**{**SATURN, "helium_share": 1.0})
    with pytest.raises(gs.StructureError):
        gs.StructureSubject(**{**SATURN, "layers": (gs.Layer("core", 1.0),)})        # the outermost layer is the envelope
    with pytest.raises(gs.StructureError):
        gs.Layer("crust", 1.0)
    with pytest.raises(gs.StructureError):
        fit(layers=(gs.Layer("core", 2 * SATURN["mass_kg"]), gs.Layer("envelope")))
    with pytest.raises(gs.StructureError):
        gs.heavy_volume_table(1.5, 0.3)


def test_the_heavy_elements_are_denser_with_depth_and_ice_is_lighter_than_rock():
    rock = gs.heavy_volume_table(0.0, 0.3)
    ice = gs.heavy_volume_table(1.0, 0.3)
    assert all(a >= b for a, b in zip(rock, rock[1:])) and all(i > r for i, r in zip(ice[:100], rock[:100]))
    assert gs.heavy_volume_table(0.5, 0.3)[10] == pytest.approx(0.5 * (ice[10] + rock[10]))


def test_without_the_tables_the_fit_says_so(tmp_path, monkeypatch):
    support.reset_caches(tmp_path / "empty", monkeypatch)
    monkeypatch.setattr(hhe_eos, "DATASET", hhe_eos.DATASET)
    with pytest.raises(fetch.DataMissing):
        fit()
