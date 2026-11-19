"""GEN.1767: the dry-surface carbonate buffer (`physics/carbonate_buffer.py`) and its place in the carbon balance of `surface_air._shares`."""

import dataclasses

import pytest

from planetgen.physics import carbonate_buffer as cb
from planetgen.physics import continental_weathering, sulfur_balance, surface_air as sa, surface_pressure, waterworld_carbon
from planetgen.util import draw

EARTH_CRUST_COLUMN_KG_M2 = 2.8e22 / 5.1e14


def test_the_pressure_of_the_urey_reaction_follows_the_data_set_and_the_sources_read():
    """[R1981: Table 2a] gives 88 bar at 740 K against the 92.1 bar Venus holds [R5313] and about 275 bar of pure CO2 near 800 K in the experiments of [R5314]."""
    assert cb.equilibrium_bar(740.0) == pytest.approx(92.1, rel=0.10)
    assert cb.equilibrium_bar(800.0) == pytest.approx(275.0, rel=0.15)
    assert cb.equilibrium_bar(500.0) == pytest.approx(0.081, rel=0.05)
    assert cb.equilibrium_bar(600.0) == pytest.approx(3.0, rel=0.05)
    temperatures = [400.0 + 25.0 * step for step in range(30)]
    pressures = [cb.equilibrium_bar(t) for t in temperatures]
    assert pressures == sorted(pressures) and len(set(pressures)) == len(pressures)       # the buffer pressure rises with the temperature


def test_the_crust_calcium_holds_about_the_carbonate_of_a_few_hundred_bar_and_the_sulfate_takes_its_share():
    earth = cb.capacity_bar(EARTH_CRUST_COLUMN_KG_M2, 0.0, 9.81)
    assert 300.0 < earth < 700.0
    calcium_mol_m2 = EARTH_CRUST_COLUMN_KG_M2 * sulfur_balance.CRUST_CALCIUM_WT / 0.040078
    assert cb.capacity_bar(EARTH_CRUST_COLUMN_KG_M2, 0.5 * calcium_mol_m2, 9.81) == pytest.approx(0.5 * earth, rel=1e-3)
    assert cb.capacity_bar(EARTH_CRUST_COLUMN_KG_M2, 2.0 * calcium_mol_m2, 9.81) == 0.0


def test_the_kept_pressure_is_the_reaction_while_the_calcium_lasts_and_the_surplus_after():
    pressure_600 = cb.equilibrium_bar(600.0)
    assert cb.kept_bar(200.0, 600.0, 500.0) == pytest.approx(pressure_600)                    # plenty of calcium: the air holds the pressure of the reaction
    assert cb.kept_bar(200.0, 600.0, 50.0) == pytest.approx(150.0)                           # 50 bar of carbonate at most: the rest stays in the air
    assert cb.kept_bar(200.0, 600.0, 0.0) == pytest.approx(200.0)                            # no calcium left: all the carbon stays
    assert cb.kept_bar(2.0, 700.0, 500.0) == pytest.approx(2.0)                              # a thin supply under the pressure of the reaction: nothing carbonates
    assert cb.kept_bar(200.0, 1000.0, 500.0) == pytest.approx(200.0)                         # hot enough that the reaction pressure is over any supply
    assert cb.kept_bar(200.0, cb.DRY_CARBONATION_MIN_K - 1.0, 500.0) == pytest.approx(200.0)  # under the kinetic floor a dry rock keeps its air as it is


def _dry_setting(released_j=1.0e31, **fields):
    fields = {"ocean_cover": 0.0, "water_kg": 0.0, **fields}
    return sa.Setting(gravity_m_s2=9.81, radius_m=6.371e6, mantle_kg=4.0e24, released_j=released_j, heat_flux_w_m2=0.0614, 
                      absorbed_flux_w_m2=240.0, released_at_seal_j=2.0e29, subject_of=None, **fields)


def _share(setting, temperature_k, weathers=False, carbon=1.0e-3):
    draw.set_run_seed(26)
    design = {}
    efficiencies = surface_pressure.draw_efficiencies(design)
    cycle = continental_weathering.draw_cycle(design)
    lifetimes = sulfur_balance.draw_lifetimes(design)
    record = {}
    water_world = sa.WaterWorld(False, waterworld_carbon.draw_waterworld({}), waterworld_carbon.draw_seafloor_law({}))
    share = sa._shares(setting, efficiencies, cycle, (lifetimes, 4.5 * sa.GY_S, 0.9), False, weathers, temperature_k, water_world, {}, record)({"C": carbon})["C"]
    return share, record.get("sink")


def test_a_dry_hot_body_holds_the_reaction_pressure_and_a_cool_or_hotter_one_keeps_its_carbon():
    setting = _dry_setting()
    unbuffered, none_ = _share(setting, 400.0)
    assert none_ is None                                                     # under the floor: no sink, the whole degassed carbon stays
    held, sink = _share(setting, 600.0)
    assert sink is not None and sink.solved_bar == pytest.approx(cb.equilibrium_bar(600.0))
    assert held == pytest.approx(unbuffered * sink.solved_bar / sink.supply_bar)
    assert held < unbuffered
    hot, none_ = _share(setting, 1200.0)
    assert none_ is None and hot == pytest.approx(unbuffered)               # the reaction pressure is over the supply: nothing carbonates
    assert sink.kept_at(700.0) > sink.solved_bar and sink.kept_at(300.0) == sink.supply_bar        # the sink answers at another temperature (the life stage asks again)


def test_a_crust_with_little_calcium_keeps_the_surplus_and_a_giant_or_a_weathering_body_is_not_the_buffers():
    little = _dry_setting(released_j=1.0e30)
    _, sink = _share(little, 600.0, carbon=1.0e-2)
    crust = sulfur_balance.crust_column_kg_m2(1.0e30, 6.371e6)
    assert sink.supply_bar == pytest.approx(106.7, rel=0.01) and cb.capacity_bar(crust, 0.0, 9.81) == pytest.approx(36.1, rel=0.01)
    assert sink.solved_bar == pytest.approx(sink.supply_bar - cb.capacity_bar(crust, 0.0, 9.81), rel=1e-3)     # a body that built a thin crust cannot carbonate all of it
    giant = dataclasses.replace(_dry_setting(), giant=True)
    assert _share(giant, 600.0)[1] is None
    weathering = _dry_setting(ocean_cover=0.7, water_kg=1.4e21)
    _, sink = _share(weathering, 600.0, weathers=True)
    assert sink is None or sink.kept_at(600.0) != cb.equilibrium_bar(600.0)  # an ocean body's sink is the weathering law's, not the buffer's
