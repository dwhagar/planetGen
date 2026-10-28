"""
The fixed point of the air of a body (GEN.411, GEN.422, `planetgen.physics.surface_air`): freeze-out, the water at a relative humidity, the
acceleration of a creep, the cages of the ice, and the stored air of generated systems.
"""

import math

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import atmosphere, atmosphere_state, clathrate_brine, constants, gases, saturation, sulfur_balance, surface_air as sa, surface_pressure
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


def _setting(water=0.3, gravity=1.35):
    return sa.Setting(gravity_m_s2=gravity, radius_m=2.5e6, mantle_kg=1.0e22, released_j=0.0, heat_flux_w_m2=0.0, ocean_cover=0.0,
                      water_kg=1.0e18, absorbed_flux_w_m2=1.0, subject_of=None, bulk_water_fraction=water)


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
                         water_kg=1.4e21, absorbed_flux_w_m2=240.0, subject_of=None)
    draw.set_run_seed(26)
    design = {}
    efficiencies = surface_pressure.draw_efficiencies(design)
    cycle = surface_pressure.draw_cycle(design)
    lifetimes = sulfur_balance.draw_lifetimes(design)
    inventory = {"C": 0.0, "S": 1.0e-3, "N": 1.0e-4}
    age_s = 4.5 * sa.GY_S
    wet = sa._shares(setting, efficiencies, cycle, (lifetimes, age_s, 0.9), True, True, 288.0)(inventory)
    dry = sa._shares(setting, efficiencies, cycle, (lifetimes, age_s, 0.9), False, False, 288.0)(inventory)
    base = surface_pressure.degassed_share(setting.released_j, setting.mantle_kg) * efficiencies["S"]
    assert 0.0 < wet["S"] < base                        # the sinks take most of it
    assert wet["S"] == pytest.approx(base * lifetimes.wet_s / age_s, rel=1e-6)   # the steady-state burden of a day or so
    assert dry["S"] > wet["S"]                          # sulfation of the rock is slower than rain-out
    assert wet["N"] == dry["N"]                         # nitrogen has no sink here


def test_a_crust_too_small_for_the_sulfur_leaves_the_excess_in_the_air():
    setting = sa.Setting(gravity_m_s2=9.81, radius_m=EARTH_R_M, mantle_kg=4.0e24, released_j=1.0e24, heat_flux_w_m2=0.0614, ocean_cover=0.7,
                         water_kg=1.4e21, absorbed_flux_w_m2=240.0, subject_of=None)
    draw.set_run_seed(27)
    design = {}
    efficiencies = surface_pressure.draw_efficiencies(design)
    cycle = surface_pressure.draw_cycle(design)
    lifetimes = sulfur_balance.draw_lifetimes(design)
    inventory = {"C": 0.0, "S": 1.0, "N": 0.0}         # an inventory far larger than the crust this heat built
    shares = sa._shares(setting, efficiencies, cycle, (lifetimes, 4.5 * sa.GY_S, 0.9), True, True, 288.0)(inventory)
    base = surface_pressure.degassed_share(setting.released_j, setting.mantle_kg) * efficiencies["S"]
    assert shares["S"] > 0.5 * base


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


def test_every_generated_body_stores_the_density_scale_height_and_greenhouse_its_mixture_gives():
    """GEN.423 and DB.67: the molar density, the density, the scale height and the greenhouse multiplier a body stores are the ones
    its stored partial pressures give at its surface temperature and gravity: the molar mass is `sum x_i M_i`, the density `P M / (R T)`,
    the scale height `R T / (M g)`, and the multiplier is the warming `atmosphere_state.greenhouse_warming` finds for the stored mixture
    at the starlight the body absorbs, which is the one its stored temperature `T_e (1 + G)^(1/4)` is made of."""
    bodies = 0
    for body in _generated(range(10)):
        air = atmosphere.gas_partials_pa(body)
        if not air or body.atmosphere == "None":
            assert body.greenhouse_factor == 0.0 or body.greenhouse_factor is None
            continue
        bodies += 1
        gravity = body.gravity * constants.EARTH_GRAVITY
        state = atmosphere_state.air_state(air, body.surface_temperature, gravity)
        assert body.atm_molar_density == pytest.approx(state.mean_molar_mass_kg_mol, rel=1e-6)
        assert body.atm_density == pytest.approx(state.density_kg_m3, rel=1e-6)
        assert body.scale_height == pytest.approx(state.scale_height_m / 1000.0, rel=1e-6)
        assert body.bond_albedo is not None and 0.0 <= body.bond_albedo < 1.0
        absorbed = constants.STEFAN_BOLTZMANN_CONSTANT * body.surface_temperature ** 4 / (1.0 + body.greenhouse_factor)
        warming = atmosphere_state.greenhouse_warming(air, gravity, absorbed)
        assert body.greenhouse_factor == pytest.approx(warming.factor, rel=1e-6, abs=1e-9)
    assert bodies > 5


def test_a_system_replays_to_the_same_air():
    def airs():
        return [(b.atmospheric_pressure, b.surface_temperature, atmosphere.gas_partials_pa(b)) for b in _generated([5])]

    assert airs() == airs()
