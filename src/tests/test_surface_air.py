"""
The fixed point of the air of a body (GEN.411, GEN.422, `planetgen.physics.surface_air`): freeze-out, the water at a relative humidity, the
acceleration of a creep, the cages of the ice, and the stored air of generated systems.
"""

import dataclasses
import math
import types

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import aerosols, atmosphere, atmosphere_state, clathrate_brine, climate_branches, climate_draws, constants, continental_weathering, gases, giant_draw, greenhouse, greenhouse_table, orbit_elements, planet_type, radiation_field, saturation, sulfur_balance, surface_air as sa, surface_pressure
from planetgen import tuning
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy

EARTH_G = constants.EARTH_GRAVITY
EARTH_R_M = 6.371e6
EARTH_AREA = 4.0 * math.pi * EARTH_R_M ** 2


def _moles_for_bar(formula, bar, gravity=EARTH_G, area=EARTH_AREA):
    """The moles of one gas that make a hydrostatic column of `bar` of it."""
    return bar * 1.0e5 * area / gravity / gases.gas(formula).molar_mass


def _partials(left, gravity=EARTH_G, area=EARTH_AREA):
    weight = sum(n * gases.gas(f).molar_mass for f, n in left.items())
    pressure = weight * gravity / area
    total = sum(left.values())
    return {f: n / total * pressure for f, n in left.items()}


def test_a_gas_under_its_saturation_pressure_is_left_alone():
    moles = {"N2": _moles_for_bar("N2", 1.0), "CO2": _moles_for_bar("CO2", 0.01)}
    left, condensed, _ = sa.freeze_out(moles, EARTH_G, EARTH_AREA, 288.0)
    assert left == pytest.approx(moles, rel=1e-12)
    assert all(abs(n) < 1e-6 * moles["N2"] for n in condensed.values())


def test_a_gas_over_its_saturation_pressure_stands_on_the_curve_and_the_moles_are_kept():
    temperature = 150.0
    moles = {"N2": _moles_for_bar("N2", 0.5), "CO2": _moles_for_bar("CO2", 0.5)}
    left, condensed, _ = sa.freeze_out(moles, EARTH_G, EARTH_AREA, temperature)
    partial = _partials(left)
    assert partial["CO2"] == pytest.approx(saturation.pressure_pa("CO2", temperature), rel=1e-6)
    assert condensed["CO2"] > 0.0
    for formula in moles:
        assert left[formula] + condensed[formula] == pytest.approx(moles[formula], rel=1e-12)
    assert left["N2"] == pytest.approx(moles["N2"])


def test_the_condensing_of_one_gas_changes_the_partial_pressures_of_the_rest():
    """The column is well mixed, so a gas has the share `x` of the pressure of the whole column: when the heavy gas leaves, the mean
    molar mass and the pressure fall and the share of the nitrogen rises."""
    moles = {"N2": _moles_for_bar("N2", 0.5), "CO2": _moles_for_bar("CO2", 0.5)}
    before = _partials(moles)["N2"]
    left, _, _ = sa.freeze_out(moles, EARTH_G, EARTH_AREA, 150.0)
    after = _partials(left)["N2"]
    share = left["N2"] / sum(left.values())
    weight = sum(n * gases.gas(f).molar_mass for f, n in left.items())
    assert after == pytest.approx(share * weight * EARTH_G / EARTH_AREA, rel=1e-12)
    assert after != pytest.approx(before, rel=1e-3)


def test_water_over_an_ocean_stands_at_the_relative_humidity_of_its_saturation_pressure():
    temperature = 288.0
    moles = {"N2": _moles_for_bar("N2", 0.78), "O2": _moles_for_bar("O2", 0.21)}
    reservoir = _moles_for_bar("H2O", 10.0)  # more water than the air can hold
    left, condensed, remaining = sa.freeze_out(moles, EARTH_G, EARTH_AREA, temperature, {"H2O": (sa.RELATIVE_HUMIDITY, reservoir)})
    partial = _partials(left)
    assert partial["H2O"] == pytest.approx(sa.RELATIVE_HUMIDITY * saturation.pressure_pa("H2O", temperature), rel=1e-6)
    assert remaining["H2O"] + left["H2O"] == pytest.approx(reservoir, rel=1e-12)
    assert condensed["H2O"] < 0.0  # the reservoir fed the air


def test_a_small_reservoir_is_used_up_and_the_air_is_drier_than_the_humidity():
    temperature = 288.0
    moles = {"N2": _moles_for_bar("N2", 1.0)}
    reservoir = _moles_for_bar("H2O", 1.0e-4)
    left, _, remaining = sa.freeze_out(moles, EARTH_G, EARTH_AREA, temperature, {"H2O": (sa.RELATIVE_HUMIDITY, reservoir)})
    assert remaining["H2O"] == pytest.approx(0.0, abs=1e-9 * reservoir)
    assert left["H2O"] == pytest.approx(reservoir, rel=1e-9)
    assert _partials(left)["H2O"] < sa.RELATIVE_HUMIDITY * saturation.pressure_pa("H2O", temperature)


def test_a_body_with_less_than_half_the_earth_heat_flow_has_no_long_lived_geology():
    earth = surface_pressure.HEAT_FLUX_EARTH_W_M2
    assert not sa.long_lived_geology(0.49 * earth)
    assert sa.long_lived_geology(sa.LONG_LIVED_FLUX_FRACTION * earth)
    assert sa.long_lived_geology(earth)


def test_a_steady_creep_is_extrapolated_to_its_limit_and_a_swing_or_noise_is_not():
    """x -> 0.99 x + 0.1 has its fixed point at 10 and the plain damped step crawls toward it; the residuals fall by a steady ratio, so the
    step is the distance to the limit of that geometric series."""
    def residual(x):
        return 0.99 * x + 0.1 - x

    x, history = 1.0, []
    for _ in range(3):
        history.append(residual(x))
        x += sa.DAMPING * residual(x)
    step = sa._step(sa.DAMPING, residual(x), history)
    assert x + step == pytest.approx(10.0, rel=1e-3) or step == pytest.approx(sa.MAX_GAIN * sa.DAMPING * residual(x))
    assert abs(step) > abs(sa.DAMPING * residual(x))
    # a residual that changed sign is a swing, one that grew is a drift, one that fell by an unsteady ratio is noise: the damped step
    assert sa._step(sa.DAMPING, -0.5, [0.8, 0.4]) == pytest.approx(-0.5 * sa.DAMPING)
    assert sa._step(sa.DAMPING, 0.9, [0.4, 0.6]) == pytest.approx(0.9 * sa.DAMPING)
    assert sa._step(sa.DAMPING, 0.2, [1.0, 0.9]) == pytest.approx(0.2 * sa.DAMPING)
    # a first or second step has no history
    assert sa._step(sa.DAMPING, 0.5, []) == pytest.approx(0.5 * sa.DAMPING)
    assert sa._step(sa.DAMPING, 0.5, [0.7]) == pytest.approx(0.5 * sa.DAMPING)


LAND = sa.WaterWorld(sealed=False, draws=None, law=None)


def _setting(water=0.3, gravity=1.35):
    return sa.Setting(gravity_m_s2=gravity, radius_m=2.5e6, mantle_kg=1.0e22, released_j=0.0, heat_flux_w_m2=0.0, ocean_cover=0.0,
                      water_kg=1.0e18, absorbed_flux_w_m2=1.0, released_at_seal_j=0.0, subject_of=None, bulk_water_fraction=water)


def test_a_cold_icy_body_holds_its_guest_gases_in_the_ice_and_a_warm_or_dry_one_does_not():
    floor = clathrate_brine.dissociation_pressure_pa("CH4", 100.0)
    air = {"CH4": 10.0 * floor, "N2": 1.0e3, "Ar": 5.0e4}
    held = sa._hold_in_hydrate(air, 100.0, _setting(), 1.0e6)
    assert held["CH4"] == pytest.approx(floor)
    assert held["Ar"] == air["Ar"]  # argon forms no hydrate here
    assert sa._hold_in_hydrate(air, 290.0, _setting(), 1.0e6) == air  # no cages above the ice point
    assert sa._hold_in_hydrate(air, 100.0, _setting(water=1e-4), 1.0e6) == air  # a dry body has no ice skin
    thin = sa._hold_in_hydrate(air, 100.0, _setting(), 1.0)
    assert thin["CH4"] > held["CH4"]  # a thin skin holds less


def test_an_ocean_dissolves_sulfur_dioxide_and_leaves_the_nitrogen():
    setting = _setting()
    moles = {"N2": 1.0e18, "SO2": 1.0e15}
    left, dissolved = sa._dissolve_sulfur(moles, 290.0, setting)
    assert left["N2"] == moles["N2"]
    assert 0.0 < left["SO2"] < moles["SO2"]
    assert left["SO2"] + dissolved["SO2"] == pytest.approx(moles["SO2"])


def test_the_sulfur_the_air_keeps_is_what_the_sulfate_and_sulfide_sinks_leave():
    setting = sa.Setting(gravity_m_s2=9.81, radius_m=EARTH_R_M, mantle_kg=4.0e24, released_j=1.0e30, heat_flux_w_m2=0.0614, ocean_cover=0.7,
                         water_kg=1.4e21, absorbed_flux_w_m2=240.0, released_at_seal_j=2.0e29, subject_of=None)
    draw.set_run_seed(26)
    design = {}
    efficiencies = surface_pressure.draw_efficiencies(design)
    cycle = continental_weathering.draw_cycle(design)
    lifetimes = sulfur_balance.draw_lifetimes(design)
    inventory = {"C": 0.0, "S": 1.0e-3, "N": 1.0e-4}
    age_s = 4.5 * sa.GY_S
    wet = sa._shares(setting, efficiencies, cycle, (lifetimes, age_s, 0.9), True, True, 288.0, LAND, {})(inventory)
    dry_setting = dataclasses.replace(setting, water_kg=0.0)    # GEN.1562 N6: a layer of water over the full depth makes the sulfate sink wet whether or not an ocean stands
    dry = sa._shares(dry_setting, efficiencies, cycle, (lifetimes, age_s, 0.9), False, False, 288.0, LAND, {})(inventory)
    base = surface_pressure.degassed_share(setting.released_j, setting.mantle_kg) * efficiencies["S"]
    assert 0.0 < wet["S"] < base                        # the sinks take most of it
    assert wet["S"] == pytest.approx(base * lifetimes.wet_s / age_s, rel=1e-6)   # the steady-state burden of a day or so
    assert dry["S"] > wet["S"]                          # sulfation of the rock is slower than rain-out
    assert wet["N"] == dry["N"]                         # nitrogen has no sink here


def test_the_shares_the_surface_holds_start_from_the_heat_released_times_each_elements_degassing_efficiency():
    """TEST.538 (the falsifier of stand-ins AM28 and AM29): hydrogen is degassed at a drawn share of the mantle's inside the registered
    bracket, carbon, nitrogen and argon in full, sulfur at a half, and the air sees the heat released times each of them."""
    draw.set_run_seed(28)
    design = {}
    efficiencies = surface_pressure.draw_efficiencies(design)
    assert surface_pressure.DEGAS_EFFICIENCY_H[0] <= efficiencies["H"] <= surface_pressure.DEGAS_EFFICIENCY_H[1]
    assert design["degas_efficiency_h"] == efficiencies["H"]
    assert {element: efficiencies[element] for element in ("C", "Ar", "S")} == {"C": 1.0, "Ar": 1.0, "S": 0.5}
    assert surface_pressure.DEGAS_EFFICIENCY_N[0] <= efficiencies["N"] <= surface_pressure.DEGAS_EFFICIENCY_N[1]       # GEN.1584: nitrogen is drawn
    assert design["degas_efficiency_n"] == efficiencies["N"]
    cycle = continental_weathering.draw_cycle(design)
    lifetimes = sulfur_balance.draw_lifetimes(design)
    for released_j in (1.0e30, 0.0):
        setting = sa.Setting(gravity_m_s2=9.81, radius_m=EARTH_R_M, mantle_kg=4.0e24, released_j=released_j, heat_flux_w_m2=0.0614,
                             ocean_cover=0.7, water_kg=1.4e21, absorbed_flux_w_m2=240.0, released_at_seal_j=2.0e29, subject_of=None)
        shares = sa._shares(setting, efficiencies, cycle, (lifetimes, 4.5 * sa.GY_S, 0.9), False, False, 288.0, LAND, {})({})
        base = surface_pressure.degassed_share(released_j, setting.mantle_kg)
        assert (base > 0.0) == (released_j > 0.0)
        for element, efficiency in efficiencies.items():
            assert shares[element] == pytest.approx(base * efficiency)


def test_a_crust_too_small_for_the_sulfur_leaves_the_excess_in_the_air():
    setting = sa.Setting(gravity_m_s2=9.81, radius_m=EARTH_R_M, mantle_kg=4.0e24, released_j=1.0e24, heat_flux_w_m2=0.0614, ocean_cover=0.7,
                         water_kg=1.4e21, absorbed_flux_w_m2=240.0, released_at_seal_j=2.0e29, subject_of=None)
    draw.set_run_seed(27)
    design = {}
    efficiencies = surface_pressure.draw_efficiencies(design)
    cycle = continental_weathering.draw_cycle(design)
    lifetimes = sulfur_balance.draw_lifetimes(design)
    inventory = {"C": 0.0, "S": 1.0, "N": 0.0}         # an inventory far larger than the crust this heat built
    shares = sa._shares(setting, efficiencies, cycle, (lifetimes, 4.5 * sa.GY_S, 0.9), True, True, 288.0, LAND, {})(inventory)
    base = surface_pressure.degassed_share(setting.released_j, setting.mantle_kg) * efficiencies["S"]
    assert shares["S"] > 0.5 * base


def _ocean_setting(ocean_cover, water_kg, released_j=1.0e30, released_at_seal_j=2.0e29, **fields):
    return sa.Setting(gravity_m_s2=9.81, radius_m=EARTH_R_M, mantle_kg=4.0e24, released_j=released_j, heat_flux_w_m2=0.0614,
                      ocean_cover=ocean_cover, water_kg=water_kg, absorbed_flux_w_m2=240.0, released_at_seal_j=released_at_seal_j,
                      subject_of=None, **fields)


def _carbon_share(setting, water_world, temperature_k=288.0, carbon=1.0e-2):
    draw.set_run_seed(26)
    design = {}
    efficiencies = surface_pressure.draw_efficiencies(design)
    cycle = continental_weathering.draw_cycle(design)
    lifetimes = sulfur_balance.draw_lifetimes(design)
    return sa._shares(setting, efficiencies, cycle, (lifetimes, 4.5 * sa.GY_S, 0.9), True, True, temperature_k, water_world, {})({"C": carbon})["C"]


def test_the_land_weathering_works_on_shrinks_with_the_ice_cover_of_the_bodys_own_law():
    """GEN.1697: no rain on bare rock under a global ice sheet; a body with no ice law keeps its whole land."""
    law = climate_branches.IceAlbedoLaw(ice_free=0.3, ice_covered=0.62)
    held = _ocean_setting(0.7, 1.4e21)
    iced = dataclasses.replace(held, albedo_law=law)
    assert sa._weathering_land(held, 200.0) == pytest.approx(0.3)
    assert sa._weathering_land(iced, 310.0) == pytest.approx(0.3)
    assert sa._weathering_land(iced, 200.0) == 0.0
    middle = sa._weathering_land(iced, 0.5 * (law.free_above_k + law.covered_below_k))
    assert 0.0 < middle < 0.3 and middle == pytest.approx(0.3 * (1.0 - law.ice_fraction(0.5 * (law.free_above_k + law.covered_below_k))))


def test_a_frozen_world_with_the_ice_law_keeps_the_carbon_the_volcanoes_give_and_a_warm_one_does_not():
    """GEN.1697: weathering all but stops under the ice, so the carbon builds up until it melts it (Hoffman et al. 2017)."""
    law = climate_branches.IceAlbedoLaw(ice_free=0.3, ice_covered=0.62)
    held = _ocean_setting(0.7, 1.4e21)
    iced = dataclasses.replace(held, albedo_law=law, starlight_w_m2=1360.0)      # the energy limit of the runoff reads the starlight the ice law lets in
    open_ = sa.WaterWorld(False, sa.waterworld_carbon.draw_waterworld({}), sa.waterworld_carbon.draw_seafloor_law({}))
    assert _carbon_share(iced, open_, 260.0, 1.0) > _carbon_share(held, open_, 260.0, 1.0)      # (at 230 K the runoff is nil with or without the ice: the law has no rain there)
    assert _carbon_share(iced, open_, 230.0, 1.0) == pytest.approx(_carbon_share(held, open_, 230.0, 1.0), rel=1e-6)    # all of what the volcanoes give stays
    assert _carbon_share(iced, open_, 310.0, 1.0) == pytest.approx(_carbon_share(held, open_, 310.0, 1.0))


def test_the_water_world_is_sealed_only_when_its_ocean_floor_reaches_the_seal_pressure():
    """GEN.830: a shallow ocean (Earth's, 1.4e21 kg) leaves the seafloor at about 0.04 GPa, a body that is a tenth water at several GPa."""
    shallow = sa.water_world_of(_ocean_setting(0.7, 1.4e21), {})
    deep = sa.water_world_of(_ocean_setting(1.0, 4.0e23), {})
    assert not shallow.sealed and deep.sealed


def test_a_sealed_ocean_world_degasses_what_the_mantle_released_by_the_seal_time_and_not_after():
    setting = _ocean_setting(1.0, 4.0e23)
    draws, law = sa.waterworld_carbon.draw_waterworld({}), sa.waterworld_carbon.draw_seafloor_law({})
    sealed, open_ = sa.WaterWorld(True, draws, law), sa.WaterWorld(False, draws, law)
    assert surface_pressure.degassed_share(setting.released_at_seal_j, setting.mantle_kg) < surface_pressure.degassed_share(
        setting.released_j, setting.mantle_kg)
    assert 0.0 < _carbon_share(setting, sealed) <= 1.0 and 0.0 < _carbon_share(setting, open_) <= 1.0


def test_an_ocean_world_with_no_land_keeps_its_carbon_in_the_ocean_and_the_hotter_ocean_keeps_less():
    setting = _ocean_setting(1.0, 4.0e23)
    sealed = sa.WaterWorld(True, sa.waterworld_carbon.draw_waterworld({}), sa.waterworld_carbon.draw_seafloor_law({}))
    cool, hot = _carbon_share(setting, sealed, 280.0), _carbon_share(setting, sealed, 560.0)
    assert hot >= cool  # the solubility of carbon dioxide falls with the temperature, so the air holds more of it
    assert cool < 1.0   # an ocean takes some of it at all


def test_volcanoes_bring_up_no_more_hydrogen_than_the_mantle_holds_as_water():
    """GEN.837: recipes b2 and b3 deliver 0.1 to 0.2 percent of the planet's mass as water, a body of 0.03 percent water has none to give."""
    draw.set_run_seed(26)
    design = {}
    efficiencies = surface_pressure.draw_efficiencies(design)
    cycle = continental_weathering.draw_cycle(design)
    lifetimes = sulfur_balance.draw_lifetimes(design)
    inventory = {"H": 5.0, "C": 0.0, "S": 0.0, "N": 0.0}       # moles per 100 g of the mantle: 0.9 percent of the mantle's mass as water
    def share(mantle_water_kg):
        setting = _ocean_setting(0.7, 1.4e21, mantle_water_kg=mantle_water_kg)
        return sa._shares(setting, efficiencies, cycle, (lifetimes, 4.5 * sa.GY_S, 0.9), True, True, 288.0, LAND, {})(inventory)["H"]
    base = surface_pressure.degassed_share(1.0e30, 4.0e24) * efficiencies["H"]
    assert share(math.inf) == pytest.approx(base)                # no limit stated: the supply is as before
    assert share(0.0) == 0.0                                      # a dry mantle gives none
    some = 1.0e19                                                 # kg of water
    assert share(some) == pytest.approx(2.0 * tuning.MANTLE_OUTGAS_SHARE * some / 0.018015 / (5.0 * 4.0e24 / 0.1), rel=1e-3)
    assert 0.0 < share(some) < base


def _generated(seeds):
    for seed in seeds:
        with deterministic_entropy(seed):
            system = StarSystem(system_config=SystemConfig())
        for planet in system.planets:
            if hasattr(planet, "moons"):
                yield planet
                yield from planet.moons


def test_the_stored_partials_add_up_to_the_pressure_and_the_state_follows_from_them():
    bodies = 0
    for body in _generated(range(8)):
        if not body.atmospheric_pressure:
            assert atmosphere.gas_partials_pa(body) == {}
            continue
        bodies += 1
        air = atmosphere.gas_partials_pa(body)
        assert sum(air.values()) == pytest.approx(body.atmospheric_pressure, rel=1e-9)
        assert body.atm_density > 0.0 and body.scale_height > 0.0 and body.atm_molar_density > 0.0
        mean = sum(air[formula] * gases.gas(formula).molar_mass for formula in air) / body.atmospheric_pressure
        assert body.atm_molar_density == pytest.approx(mean, rel=1e-6)
        if body.body_type == "g":
            assert body.atmospheric_pressure == pytest.approx(sa.REFERENCE_PRESSURE_BAR * 1.0e5, rel=1e-9)
    assert bodies > 0


def test_every_generated_body_stores_the_density_scale_height_and_greenhouse_its_mixture_gives(monkeypatch):
    """GEN.423 and DB.67: the molar density, the density, the scale height and the greenhouse multiplier a body stores are the ones
    its stored partial pressures give at its surface temperature and gravity: the molar mass is `sum x_i M_i`, the density `P M / (R T)`,
    the scale height `R T / (M g)`, and the multiplier is the warming `atmosphere_state.greenhouse_warming` finds for the stored mixture
    at the starlight the body absorbs, with the sunlight its haze takes aloft (GEN.851), which is the one its stored temperature `T_e (1 + G)^(1/4)` is made
    of; the haze is read at the stored temperature, the fixed point's last step read it at the one before (under 0.1 K away; a giant's haze is the one of its single step), and a giant that radiates more than it absorbs
    (its internal flux) has a haze share of the radiated flux, which the test scales from the starlight it absorbs."""
    unsettled, real_solve = [], sa.solve

    def recording_solve(setting, start_temperature_k):
        result = real_solve(setting, start_temperature_k)
        if result is not None and result.partials_pa and not result.converged:
            unsettled.append(result.pressure_pa)         # a fixed point that has not settled in `MAX_STEPS` keeps its last step (`Air.converged`)
        return result

    monkeypatch.setattr(sa, "solve", recording_solve)
    bodies = skipped = 0
    for body in _generated(range(10)):
        air = atmosphere.gas_partials_pa(body)
        if not air or body.atmosphere == "None":
            assert body.greenhouse_factor == 0.0 or body.greenhouse_factor is None
            continue
        bodies += 1
        if any(body.atmospheric_pressure == pytest.approx(pressure, rel=1e-9) for pressure in unsettled):
            skipped += 1        # the test reads the settled airs: an unsettled one (a 2-cycle between a cool air and a hot one, GEN.1323 hydrogen on a hazy methane air) is a step off its own mixture
            continue
        gravity = body.gravity * constants.EARTH_GRAVITY
        state = atmosphere_state.air_state(air, body.surface_temperature, gravity)
        assert body.atm_molar_density == pytest.approx(state.mean_molar_mass_kg_mol, rel=1e-6)
        assert body.atm_density == pytest.approx(state.density_kg_m3, rel=1e-6)
        assert body.scale_height == pytest.approx(state.scale_height_m / 1000.0, rel=1e-6)
        assert body.bond_albedo is not None and 0.0 <= body.bond_albedo < 1.0
        absorbed = constants.STEFAN_BOLTZMANN_CONSTANT * body.surface_temperature ** 4 / (1.0 + body.greenhouse_factor)
        if body.atmospheric_pressure > greenhouse_table.PRESSURE_MAX_PA or body.surface_temperature >= sa.ionization.ceiling_k(air):
            continue            # the table holds its edge above 1000 bar and the fixed point stops at the ionization ceiling of the air: no reproduction there
        if getattr(body, "life_class", None):
            continue            # life changes the air after the temperature is solved (its methane and oxygen, GEN.855), and with it the haze the temperature was solved under
        haze_temperature = body.surface_temperature
        if not planet_type.has_ground(body):
            # a giant takes one step of the fixed point, so its haze is the one of the temperature the step starts from (the no-air one at the stored albedo)
            starlight = body.star.luminosity / (4.0 * math.pi * (body._star_distance_au * constants.AU_M) ** 2)
            haze_temperature = climate_draws.equilibrium_temperature_k(body.bond_albedo, starlight)
        haze = aerosols.haze_of(air, gravity, haze_temperature, aerosols.draw_haze_draws(aerosols.haze_stream_of(body)),
                                body.star.temperature, aerosols.far_uv_of(body.star, body._star_distance_au))
        if not planet_type.has_ground(body):
            # GEN.1030: a giant's stored G is against its starlight alone, so `T_e (1 + G)^(1/4)` with the albedo reproduces the temperature (GEN.412);
            # the air it has is warmed under the starlight and its own heat together
            radiated = absorbed + giant_draw.own_flux_w_m2(body)
            warming = atmosphere_state.greenhouse_warming(air, gravity, radiated, haze.absorbed_share * absorbed / radiated, haze.pressure_pa)
            assert body.surface_temperature == pytest.approx(warming.surface_temperature_k, rel=5e-3)
            continue
        warming = atmosphere_state.greenhouse_warming(air, gravity, absorbed, haze.absorbed_share, haze.pressure_pa)
        assert body.greenhouse_factor == pytest.approx(warming.factor, rel=1.5e-2, abs=2e-3)       # a thin cold air's G of a few 1e-4 moves T by under 0.05 percent and cannot be read back from the stored haze closer than that
    assert bodies > 5
    # the unsettled airs are few: more would mean the fixed point stopped settling, not that the test cannot read them. The share is a property of the
    # solver (a world near its boiling point swings between a cool air and a steam one, GEN.1344), not of one sample: over seeds 10 to 59 it is 10 of
    # 132 airs on main and 8 of 129 on Train 24 (GEN.993), and a draw of ten seeds ranges from none to three of 28, so a tenth was a luck bound, a fifth is three times the share
    assert skipped <= max(1, bodies // 5)


def test_a_system_replays_to_the_same_air():
    def airs():
        return [(b.atmospheric_pressure, b.surface_temperature, atmosphere.gas_partials_pa(b)) for b in _generated([5])]

    assert airs() == airs()


# --- the two bounds of the fixed point (GEN.646; stand-ins CC7 and CC8) ------------------------------------------------------------

def _earth_like_planet(seed=3):
    from planetgen.generation.planet import Planet
    from planetgen.generation.star import Star

    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    config.WEIRDNESS_DIAL = 0.0      # no oddity: W20 solves the air a second time from the other start, and these tests read one solve
    with draw.bound("surface-air-star"):
        star = Star(config)
    inner, outer = star.habitable_zone
    with draw.bound(seed):
        return Planet(config, star, star.habitable_zone, 0.5 * (inner + outer), mass=constants.EARTH_MASS_TO_KG, moon_count=0)


def test_a_runaway_greenhouse_stops_at_the_ionization_ceiling_and_the_body_is_named(monkeypatch, caplog):
    """GEN.1580: the surface stops where the main gas of the air is half ionized (`ionization.ceiling_k`), a limit no generated body should reach: the test suite makes
    a body there an error (`PLANETGEN_CEILING_IS_ERROR`, checked first), and in a run it is one normal-log line with the body, the uncapped value, the ceiling and the air.
    The ceiling is the guard of a body that can hold steam; the Earth-like body of the sample is under its runaway limit (GEN.1329 holds such a body under
    `MOIST_LIMIT_K`), so the test declares it past the limit."""
    import dataclasses
    import logging

    from planetgen.physics import air_final

    real = sa.atmosphere_state.greenhouse_warming
    monkeypatch.setattr(sa.atmosphere_state, "greenhouse_warming",
                        lambda *args: dataclasses.replace(real(*args), surface_temperature_k=1.0e6))
    monkeypatch.setattr(sa, "_steam_possible", lambda *args: True)
    with pytest.raises(air_final.CeilingReached, match="asked for 1000000 K"):                # the suite's switch: a body at the ceiling is an error
        _earth_like_planet()
    monkeypatch.delenv(air_final.CEILING_IS_ERROR_VARIABLE)
    caplog.clear()
    with caplog.at_level(logging.INFO, logger="planetgen"):
        planet = _earth_like_planet()
    assert planet.atmospheric_pressure > 0.0
    assert planet.surface_temperature == pytest.approx(sa.ionization.ceiling_k(atmosphere.gas_partials_pa(planet)))
    lines = [record.getMessage() for record in caplog.records if "Ionization ceiling" in record.getMessage()]
    assert len(lines) == 1 and "asked for 1000000 K" in lines[0] and "main gas" in lines[0]


def test_a_pressure_that_falls_to_a_vacuum_is_airless(monkeypatch):
    real = sa.mixture_draw.draw_mixture
    calls = []

    def fails_after_the_first(subject, **kwargs):
        calls.append(subject)
        if len(calls) > 1:
            raise sa.gas_equilibrium.EquilibriumError("no convergence")
        return real(subject, **kwargs)

    monkeypatch.setattr(sa.mixture_draw, "draw_mixture", fails_after_the_first)
    solved = []
    real_solve = sa.solve
    monkeypatch.setattr(sa, "solve", lambda *args: solved.append(real_solve(*args)) or solved[-1])
    planet = _earth_like_planet()                  # GEN.1580: at an Earth-like pressure a failure is not a vacuum: the body keeps the last air the equilibrium solved, unsettled
    assert planet.atmospheric_pressure > 0.0 and solved and not solved[0].converged
    assert len(calls) == 2 + sa.MAX_EQUILIBRIUM_STALLS                # the fixed point stepped back a few times first
    calls.clear()
    monkeypatch.setattr(sa, "VACUUM_PRESSURE_PA", 1.0e12)
    planet = _earth_like_planet()                  # with the vacuum edge above it, the same failure leaves an airless body
    assert planet.atmospheric_pressure == 0.0


def test_a_stalled_equilibrium_is_tried_again_from_the_middle_of_the_way_back(monkeypatch):
    real = sa.mixture_draw.draw_mixture
    calls = []

    def fails_once(subject, **kwargs):
        calls.append(subject)
        if len(calls) == 2:
            raise sa.gas_equilibrium.EquilibriumError("no convergence")
        return real(subject, **kwargs)

    monkeypatch.setattr(sa.mixture_draw, "draw_mixture", fails_once)
    planet = _earth_like_planet()
    assert len(calls) > 2 and planet.atmospheric_pressure > 0.0


# --- GEN.972: the ice-albedo law in the fixed point -------------------------------------------------------------------------------

def _iced_body(water_mass_fraction=0.2, body_type="t"):
    return types.SimpleNamespace(body_type=body_type, mass=6.0e24, star=types.SimpleNamespace(temperature=5770.0))


def test_a_body_with_surface_water_gets_an_ice_law_that_starts_at_its_drawn_albedo_and_reaches_the_bright_ice():
    law = sa.ice_albedo_law(_iced_body(), 0.31, surface_w=0.1)
    assert law.ice_free == 0.31 and law.ice_covered == pytest.approx(0.62)
    assert law.albedo(320.0) == 0.31 and law.albedo(200.0) == pytest.approx(0.62)


def test_a_bright_ice_never_covers_under_the_albedo_of_the_ice_free_body():
    law = sa.ice_albedo_law(_iced_body(), 0.8, surface_w=0.1)
    assert law.ice_covered == 0.8 and law.albedo(150.0) == 0.8


def test_a_giant_and_a_dry_body_keep_the_albedo_they_drew():
    assert sa.ice_albedo_law(_iced_body(body_type="g"), 0.3, surface_w=0.1) is None
    assert sa.ice_albedo_law(_iced_body(), 0.3, surface_w=0.0) is None
    assert sa.ice_albedo_law(_iced_body(), 0.3, surface_w=None) is None


def test_the_flux_of_a_setting_follows_the_law_at_the_temperature_and_is_held_without_one():
    held = _setting()
    assert held.absorbed_flux_at(100.0) == held.absorbed_flux_at(400.0) == held.absorbed_flux_w_m2
    law = sa.climate_branches.IceAlbedoLaw(ice_free=0.3, ice_covered=0.6)
    follows = dataclasses.replace(held, starlight_w_m2=1000.0, albedo_law=law)
    assert follows.absorbed_flux_at(320.0) == pytest.approx(0.7 * 250.0)
    assert follows.absorbed_flux_at(200.0) == pytest.approx(0.4 * 250.0)


def test_a_cold_world_with_the_law_is_brighter_and_colder_than_the_same_world_with_the_ice_free_albedo(monkeypatch):
    from tests.test_climate_sweep import _earth_twin, _hold_the_air_of_the_earth

    with monkeypatch.context() as held_air:
        _hold_the_air_of_the_earth(held_air)
        law_free = _earth_twin(150.0)
        held = dataclasses.replace(law_free, albedo_law=None, starlight_w_m2=0.0)
        with draw.bound("ice-law"):
            with_ice = sa.solve(law_free, 255.0)
        with draw.bound("ice-law"):
            without = sa.solve(held, 255.0)
    assert with_ice.albedo is not None and without.albedo is None
    assert with_ice.albedo > sa.climate_branches.ICE_FREE_ALBEDO
    assert with_ice.temperature_k < without.temperature_k


def test_the_albedo_and_greenhouse_a_generated_body_stores_give_back_its_temperature_with_the_ice_law():
    checked = 0
    for planet in _generated(range(8)):
        if planet.body_type != "t" or planet.is_moon or (planet.water_mass_fraction or 0) <= 1.0e-4:
            continue
        flux = radiation_field.starlight_w_m2(planet.star, orbit_elements.flux_distance_au(float(planet.distance), planet.eccentricity)) * (planet.light_factor or 1.0)       # GEN.1441: every star's light
        assert planet.surface_temperature == pytest.approx(climate_draws.surface_temperature_k(planet.bond_albedo, planet.greenhouse_factor, flux), rel=1e-6)
        checked += 1
    assert checked > 0


def test_an_airless_wet_body_settles_on_the_albedo_its_temperature_keeps():
    law = sa.climate_branches.IceAlbedoLaw(ice_free=0.3, ice_covered=0.62)
    cold, albedo_cold = sa.airless_state(law, 30.0, climate_draws.equilibrium_temperature_k(0.3, 30.0))
    assert albedo_cold == pytest.approx(0.62) and cold == pytest.approx(climate_draws.equilibrium_temperature_k(0.62, 30.0), rel=1e-12)
    hot_light = 3500.0                  # an ice-free equilibrium over 305 K
    warm, albedo_warm = sa.airless_state(law, hot_light, climate_draws.equilibrium_temperature_k(0.3, hot_light))
    assert albedo_warm == 0.3 and warm == pytest.approx(climate_draws.equilibrium_temperature_k(0.3, hot_light), rel=1e-12)


def test_the_escape_is_compared_to_the_pressure_the_air_settled_at_without_it(monkeypatch):
    """The loss acts on the air that stands. Compared to the pressure of the step before it (already what stays) it took its share again at every step,
    and the reducing mixtures of the first steps, which a settled air does not keep, took the whole air: an Earth-mass world in the habitable zone of a Sun-like
    star fell to a vacuum in most seeds."""
    real = sa.mixture_escape.apply
    columns = []

    def recording(subject, result, column_bar=None):
        columns.append(column_bar)
        return real(subject, result, column_bar)

    monkeypatch.setattr(sa.mixture_escape, "apply", recording)
    for seed in (2, 3, 5):                                           # (a seed whose air needs no escape step, or that leaves no air at all, is no sample of the order; 6 was one until GEN.1323)
        columns.clear()
        planet = _earth_like_planet(seed)
        assert planet.atmospheric_pressure > 0.0, seed
        assert columns and None not in columns, seed                 # no step before the air settled lost anything
        assert len(set(columns)) == 1, seed                          # and every step after it compares to the one settled pressure


def _escape_flags(monkeypatch, seed, max_steps=None):
    """The `escape` flag of every mixture draw of one solve (a seed's Earth-like planet), in step order."""
    flags = []
    real = sa.mixture_draw.draw_mixture

    def recording(subject, **kwargs):
        flags.append(bool(kwargs.get("escape")))
        return real(subject, **kwargs)

    monkeypatch.setattr(sa.mixture_draw, "draw_mixture", recording)
    if max_steps is not None:
        monkeypatch.setattr(sa, "MAX_STEPS", max_steps)
    _earth_like_planet(seed)
    return flags


def test_an_air_that_settles_without_the_escape_on_its_last_allowed_step_still_gets_the_pass_that_has_it_in(monkeypatch):
    """GEN.1594: the no-escape fixed point sets the escape on and goes on; if the loop's last step was the one that settled, the escape pass never ran and the air was called converged."""
    for seed in (1, 3, 7):
        flags = _escape_flags(monkeypatch, seed)
        if True not in flags:
            continue                                               # (an air with no light gas to lose has no escape pass)
        settled_at = flags.index(True)                             # the draws before it are the no-escape ones: this many steps
        monkeypatch.undo()
        again = _escape_flags(monkeypatch, seed, max_steps=settled_at)
        assert len(again) > settled_at and True in again, (seed, settled_at, len(again))
        return
    pytest.skip("no seed of the three needs an escape pass")


def test_the_carbon_sink_of_an_air_is_the_last_steps_and_none_when_the_last_step_does_not_weather(monkeypatch):
    """GEN.1594: the record was only written on steps that weather, so an air whose last step did not kept an earlier step's sink for the life stage to balance against."""
    planet = _earth_like_planet(3)
    assert getattr(planet, "_co2_sink", None) is not None            # the Earth-like air weathers
    calls = []
    real = sa._weathers

    def first_only(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs) if len(calls) == 1 else False

    monkeypatch.setattr(sa, "_weathers", first_only)
    again = _earth_like_planet(3)
    assert len(calls) > 1
    assert getattr(again, "_co2_sink", None) is None


# --- the late accretion in the supply (GEN.819) -------------------------------------------------------------------------------

EARTH_ORBITAL_SPEED_MS = 29780.0
VENUS_ORBITAL_SPEED_MS = 35020.0


def _late(setting, seed):
    draw.set_run_seed(seed)
    design = {}
    return sa._late_accretion(setting, design), design


def _n2_kg(late, setting):
    return late["N"] * setting.mantle_kg / 0.1 * 0.0140067      # moles of N atoms at 14.0067 g/mol, kg


def test_the_late_accretion_leaves_carbon_and_nitrogen_and_leaves_the_water_to_the_hydrosphere():
    setting = _ocean_setting(0.7, 1.4e21, orbital_speed_ms=EARTH_ORBITAL_SPEED_MS)
    late, design = _late(setting, 26)
    assert set(late) == {"C", "N"} and all(atoms > 0.0 for atoms in late.values())
    assert design["late_accreted_percent"] > 0.0 and 0.01 <= design["veneer_volatile_scale"] <= 1.0   # drawn and recorded for the replay
    assert _late(setting, 26)[0] == late                                                              # the same seed gives the same delivery


def test_a_giant_a_body_with_no_star_or_no_mantle_takes_no_late_accretion():
    setting = _ocean_setting(0.7, 1.4e21, orbital_speed_ms=EARTH_ORBITAL_SPEED_MS)
    assert sa._late_accretion(dataclasses.replace(setting, orbital_speed_ms=None), {}) == {}
    assert sa._late_accretion(dataclasses.replace(setting, mantle_kg=0.0), {}) == {}
    assert sa._late_accretion(dataclasses.replace(setting, giant=True), {}) == {}


def test_the_supply_adds_the_late_nitrogen_to_the_degassed_and_the_carbonate_cycle_acts_on_the_carbon_of_both():
    draw.set_run_seed(30)       # (GEN.1584: the nitrogen draw moved the stream; at this setting the cycle takes nearly all of the added carbon, and the sign of what is left, 1e-19, depends on the seed)
    design = {}
    efficiencies = surface_pressure.draw_efficiencies(design)
    cycle = continental_weathering.draw_cycle(design)
    lifetimes = sulfur_balance.draw_lifetimes(design)
    setting = _ocean_setting(0.7, 1.4e21)
    inventory = {"C": 1.0e-2, "N": 1.0e-3}

    def shares(late):
        return sa._shares(setting, efficiencies, cycle, (lifetimes, 4.5 * sa.GY_S, 0.9), True, True, 288.0, LAND, late)(inventory)
    before, after = shares({}), shares({"C": 2.0e-3, "N": 5.0e-4})
    assert after["N"] == pytest.approx(before["N"] + 5.0e-4 / 1.0e-3)
    assert after["C"] * inventory["C"] > before["C"] * inventory["C"]                                  # more carbon reaches the air...
    assert after["C"] * inventory["C"] < before["C"] * inventory["C"] + 2.0e-3                         # ...but the cycle takes a share of the added


def test_the_earth_twin_late_accretion_leaves_the_nitrogen_the_earth_has_and_venus_within_a_factor_of_five():
    """The bracket of the design note (docs/design/composition/impact-erosion-and-late-delivery.md section 4): the Earth's 3.9e18 kg of N2
    lies inside the scatter of the delivery (the middle within a factor of four, the degassed nitrogen adds the rest), Venus's 1.1e19 kg within a factor of five of its middle."""
    earth = _ocean_setting(0.7, 1.4e21, orbital_speed_ms=EARTH_ORBITAL_SPEED_MS)
    venus = dataclasses.replace(_ocean_setting(0.0, 0.0, orbital_speed_ms=VENUS_ORBITAL_SPEED_MS), mantle_kg=3.3e24, radius_m=6.052e6,
                                gravity_m_s2=8.87, absorbed_flux_w_m2=150.0)
    for setting, observed, factor in ((earth, 3.9e18, 4.0), (venus, 1.08e19, 5.0)):
        n2_kg = sorted(_n2_kg(_late(setting, seed)[0], setting) for seed in range(40))
        assert n2_kg[0] < observed < n2_kg[-1]                                                         # inside the spread
        assert observed / factor < n2_kg[len(n2_kg) // 2] < observed * factor
