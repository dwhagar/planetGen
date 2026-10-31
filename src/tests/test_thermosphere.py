# tests/test_thermosphere.py

"""The thermosphere and thermal escape model (GEN.738, docs/design/composition/atmospheric-escape-model.md): Earth, Venus, Mars, Titan and the
mass loss of a hydrogen envelope against the nine planets of Kubyshkina et al. 2018."""

import math

import pytest

from planetgen.physics import constants, thermosphere as ts

ME = constants.EARTH_MASS_TO_KG
EARTH = dict(mass=ME, radius=6.371e6, teq=255.0)
MARS = dict(mass=0.107 * ME, radius=3.3895e6, teq=210.0)
VENUS = dict(mass=0.815 * ME, radius=6.0518e6, teq=232.0)
TITAN = dict(mass=1.3452e23, radius=2.575e6)

# Fig. 12 of Johnstone et al. 2018 (arXiv:1806.06897), read from the figure: neutral exobase temperature of the Earth against the CO2 mixing ratio
JOHNSTONE_EXOBASE = {4e-6: 1310.0, 4e-5: 1260.0, 4e-4: 1155.0, 4e-3: 985.0, 4e-2: 540.0, 4e-1: 290.0}

EARTH_AIR = {"N2": 0.78, "O2": 0.21, "Ar": 0.01, "CO2": 4e-4, "H2": 5e-7, "He": 5e-6}
VENUS_AIR = {"CO2": 0.965, "N2": 0.035, "H2O": 3e-5, "He": 1e-5}
MARS_AIR = {"CO2": 0.95, "N2": 0.026, "Ar": 0.019, "H2": 1e-5}
TITAN_AIR = {"N2": 0.95, "CH4": 0.05, "H2": 0.004}
TITAN_MESOPAUSE_K = 150.0
TITAN_XUV_EARTH = 0.011         # the Earth's 1 AU flux over 9.5 AU squared


def history(air, body, pressure_pa, age_gy, xuv_now, star=(1.0, constants.SOLAR_LUMINOSITY), mesopause=None):
    return ts.kept_shares(air, body["mass"], body["radius"], body["teq"], pressure_pa, age_gy, xuv_now, star[0], star[1], mesopause)


@pytest.mark.parametrize("co2", sorted(JOHNSTONE_EXOBASE))
def test_the_earth_exobase_follows_the_co2_cooling_of_the_upper_atmosphere_model(co2):
    temperature = ts.exobase_temperature_k(1.0, ts.mesopause_temperature_k(EARTH["teq"]), co2)
    assert temperature == pytest.approx(JOHNSTONE_EXOBASE[co2], rel=0.08)


def test_the_exobase_of_the_young_sun_is_four_thousand_kelvin_at_4_9_times_the_xuv():
    assert ts.exobase_temperature_k(4.9, ts.mesopause_temperature_k(EARTH["teq"]), 4e-4) == pytest.approx(4600.0, rel=0.06)


def test_the_exobase_rises_with_the_xuv_and_is_held_to_the_range_the_model_is_checked_in():
    meso = 199.0
    assert ts.exobase_temperature_k(0.0, meso, 4e-4) == meso
    assert ts.exobase_temperature_k(0.5, meso, 4e-4) < ts.exobase_temperature_k(1.0, meso, 4e-4) < ts.exobase_temperature_k(3.0, meso, 4e-4)
    assert ts.exobase_temperature_k(1000.0, meso, 4e-4) == ts.exobase_temperature_k(ts.XUV_FLUX_MAX_EARTH, meso, 4e-4)


def test_the_earth_exobase_altitude_follows_the_upper_atmosphere_model_now_and_at_the_young_suns_xuv():
    """Fig. 12 and 13 of Johnstone et al. 2018: 475 km at 1,155 K; about 2,000 to 2,500 km at 4,600 K (read from the figures)."""
    gravity = constants.G * EARTH["mass"] / EARTH["radius"] ** 2
    mean = ts.mean_mass_kg({"N2": 0.78, "O2": 0.21, "Ar": 0.01})
    meso = ts.mesopause_temperature_k(EARTH["teq"])
    now = (ts.exobase(EARTH["radius"], gravity, mean, 1155.0, meso)[0] - EARTH["radius"]) / 1e3
    young = (ts.exobase(EARTH["radius"], gravity, mean, 4600.0, meso)[0] - EARTH["radius"]) / 1e3
    assert now == pytest.approx(475.0, rel=0.15)
    assert 1800.0 < young < 2600.0


def test_mars_exobase_is_in_the_range_maven_measures_and_responds_as_it_does_to_the_xuv():
    meso = ts.mesopause_temperature_k(MARS["teq"])
    now = ts.exobase_temperature_k(0.43, meso, 0.95)
    assert 180.0 <= now <= 260.0
    # MAVEN: about 38 K per mW/m2 of dayside EUV, about 9.5 per mW/m2 of global mean; the model's slope is per Earth unit of 4.6 mW/m2
    slope = (ts.exobase_temperature_k(0.53, meso, 0.95) - ts.exobase_temperature_k(0.33, meso, 0.95)) / 0.2 / 4.6
    assert 4.0 < slope < 19.0


def test_venus_exobase_is_cold_because_its_air_is_carbon_dioxide():
    temperature = ts.exobase_temperature_k(1.9, ts.mesopause_temperature_k(VENUS["teq"]), 0.965)
    assert 199.0 <= temperature <= 320.0


def test_titan_exobase_and_hydrogen_flux_match_what_cassini_found():
    temperature = ts.exobase_temperature_k(TITAN_XUV_EARTH, TITAN_MESOPAUSE_K, 0.0)
    assert 150.0 <= temperature <= 190.0
    gravity = constants.G * TITAN["mass"] / TITAN["radius"] ** 2
    area = 4.0 * math.pi * TITAN["radius"] ** 2
    limit = ts.diffusion_limited_flux("H2", TITAN_AIR, gravity, TITAN_MESOPAUSE_K) * area
    assert limit == pytest.approx(9.2e27, rel=4.0)              # the escape rate of hydrogen of Cui et al. 2008 and Strobel 2010
    mean = ts.mean_mass_kg(TITAN_AIR)
    r_exo, density, _ = ts.exobase(TITAN["radius"], gravity, mean, temperature, TITAN_MESOPAUSE_K)
    mass_h2 = ts.gases.gas("H2").molar_mass / 6.02214076e23
    jeans = ts.jeans_flux(TITAN_AIR["H2"] * density, mass_h2, temperature, constants.G * TITAN["mass"], r_exo) * 4.0 * math.pi * r_exo ** 2
    assert jeans == pytest.approx(6e27, rel=4.0)
    assert jeans < limit * 10.0


def test_earth_keeps_its_nitrogen_oxygen_argon_and_carbon_dioxide_and_loses_its_hydrogen():
    kept, exobase = history(EARTH_AIR, EARTH, 1.013e5, 4.5, 1.0)
    assert exobase == pytest.approx(1155.0, rel=0.08)
    assert all(kept[gas] > 0.95 for gas in ("N2", "O2", "Ar", "CO2"))
    assert kept["H2"] < 1e-3


def test_venus_keeps_all_its_air_and_mars_keeps_its_carbon_dioxide_and_nitrogen_and_loses_hydrogen():
    venus, _ = history(VENUS_AIR, VENUS, 9.2e6, 4.5, 1.9)
    assert min(venus.values()) > 0.99
    mars, exobase = history(MARS_AIR, MARS, 636.0, 4.5, 0.43)
    assert 180.0 <= exobase <= 260.0
    assert mars["CO2"] > 0.99 and mars["N2"] > 0.99 and mars["H2"] < 1e-3


def test_the_share_kept_falls_with_the_xuv_the_star_puts_out_and_rises_with_the_gravity_of_the_body():
    air = {"N2": 0.9, "H2": 0.1}
    weak, _ = history(air, EARTH, 1.013e5, 4.5, 0.2)
    strong, _ = history(air, EARTH, 1.013e5, 4.5, 4.0)
    assert weak["H2"] >= strong["H2"]
    heavy, _ = history(air, dict(EARTH, mass=3.0 * ME, radius=1.4 * EARTH["radius"]), 1.013e5, 4.5, 1.0)
    light, _ = history(air, dict(EARTH, mass=0.3 * ME, radius=0.7 * EARTH["radius"]), 1.013e5, 4.5, 1.0)
    assert light["N2"] <= heavy["N2"] and light["H2"] <= heavy["H2"]


def test_a_small_warm_body_under_a_young_active_star_loses_the_heavy_gases_it_cannot_hold():
    moonlet = dict(mass=0.0327 * ME, radius=2.1e6, teq=300.0)
    kept, _ = history({"N2": 0.9, "H2O": 0.1}, moonlet, 2e4, 2.0, 10.0, star=(0.3, 0.02 * constants.SOLAR_LUMINOSITY))
    assert kept["N2"] < 0.5 and kept["H2O"] < 0.5


def test_an_older_body_has_lost_at_least_what_a_younger_one_has_and_a_gas_at_a_higher_pressure_is_less_eroded():
    air = {"N2": 0.9, "H2": 0.1}
    young, _ = history(air, MARS, 636.0, 1.0, 0.43)
    old, _ = history(air, MARS, 636.0, 4.5, 0.43)
    assert old["H2"] <= young["H2"]
    thin, _ = history({"N2": 0.9, "He": 0.1}, MARS, 100.0, 4.5, 3.0)
    thick, _ = history({"N2": 0.9, "He": 0.1}, MARS, 1e5, 4.5, 3.0)
    assert thick["He"] >= thin["He"]


def test_the_diffusion_limit_applies_only_to_a_gas_lighter_than_the_air():
    gravity = constants.G * EARTH["mass"] / EARTH["radius"] ** 2
    assert math.isinf(ts.diffusion_limited_flux("CO2", {"N2": 0.9, "CO2": 0.1}, gravity, 200.0))
    assert math.isinf(ts.diffusion_limited_flux("N2", {"N2": 1.0}, gravity, 200.0))
    assert ts.diffusion_limited_flux("H2", {"N2": 0.99, "H2": 0.01}, gravity, 200.0) > 0.0


def test_hydrogen_in_nitrogen_diffuses_at_the_measured_rate():
    # the binary diffusion coefficient of H2 in N2 at 273 K and 1 atm is 0.68 cm2/s (Marrero and Mason 1972)
    number_density = 101325.0 / (constants.BOLTZMANN * 273.0)
    coefficient = ts.binary_diffusion_parameter("H2", "N2", 273.0) / number_density
    assert coefficient * 1e4 == pytest.approx(0.68, rel=0.2)


def test_the_energy_limit_scales_as_the_flux_and_inversely_with_the_mass():
    base = ts.energy_limited_mass_loss_kg_s(4.6e-3, 7e6, ME)
    assert ts.energy_limited_mass_loss_kg_s(9.2e-3, 7e6, ME) == pytest.approx(2.0 * base)
    assert ts.energy_limited_mass_loss_kg_s(4.6e-3, 7e6, 2.0 * ME) == pytest.approx(0.5 * base)


# Table 2 of Kubyshkina et al. 2018 (arXiv:1810.06920): planet, escape parameter, radius (Earth), distance (AU), XUV flux (erg/cm2/s), published rate (g/s)
KUBYSHKINA_TABLE_2 = [
    ("HD209458b", 90.0, 15.45, 0.047, 1086.0, 9.6e9), ("GJ436b", 58.0, 4.25, 0.02887, 1760.0, 2.3e9),
    ("Kepler-11b", 18.0, 1.97, 0.091, 278.0, 1.9e9), ("HD189733b", 179.0, 12.74, 0.03, 24778.0, 4.5e9),
    ("GJ3470b", 37.0, 4.18, 0.03557, 1868.0, 1.6e10), ("HD149026b", 61.0, 8.04, 0.04288, 6886.0, 4.5e10),
    ("HAT-P-11b", 48.5, 4.72, 0.053, 3236.0, 1.3e10), ("55 Cnc e", 16.0, 1.99, 0.01544, 570.0, 4.9e10),
    ("HD97658b", 34.0, 2.24, 0.08, 955.0, 1.8e9),
]


@pytest.mark.parametrize("name,escape_parameter,radius,distance,flux,published", KUBYSHKINA_TABLE_2)
def test_the_hydro_based_rate_reproduces_the_published_rate_of_each_planet_of_the_papers_table(name, escape_parameter, radius, distance, flux, published):
    teq = 1000.0
    mass = escape_parameter * constants.BOLTZMANN * teq * radius * 6.371e6 / (constants.G * ts._M_HYDROGEN_KG) / ME
    assert ts.hydro_based_mass_loss_g_s(mass, radius, teq, distance, flux) == pytest.approx(published, rel=1.5)       # within a factor of 2.5


def test_the_hydro_based_rate_grows_with_the_flux_and_falls_as_the_planet_grows_heavier():
    base = ts.hydro_based_mass_loss_g_s(5.0, 2.0, 800.0, 0.05, 1000.0)
    assert ts.hydro_based_mass_loss_g_s(5.0, 2.0, 800.0, 0.05, 3000.0) > base
    assert ts.hydro_based_mass_loss_g_s(20.0, 2.0, 800.0, 0.05, 1000.0) < base


def test_the_xuv_flux_comes_from_the_stored_value_else_from_the_star_and_is_zero_without_one():
    class Subject:
        def __init__(self, quantities):
            self.quantities = quantities

        def get(self, name, default=None):
            return self.quantities.get(name, default)
    assert ts.xuv_flux_earth(Subject({"xuv_flux_earth": 2.5})) == 2.5
    assert ts.xuv_flux_earth(Subject({})) == 0.0
    sun = Subject({"host_mass_solar": 1.0, "age_gy": 4.6, "luminosity_solar": 1.0, "insolation": 1.0})
    assert ts.xuv_flux_earth(sun) == pytest.approx(1.0, rel=0.5)
    far = Subject({"host_mass_solar": 1.0, "age_gy": 4.6, "luminosity_solar": 1.0, "insolation": 0.25})
    assert ts.xuv_flux_earth(far) == pytest.approx(ts.xuv_flux_earth(sun) * 0.25, rel=1e-6)


# --- GEN.841: the exobase search is bounded and the caches change no value ---------------------------------------------------

@pytest.mark.parametrize("args, expected", [
    ((6.371e6, 9.81, 4.65e-26, 1100.0, 200.0), (6792504.786904568, 132120832223323.83, 37844.14551331694)),
    ((3.4e6, 3.7, 7.3e-26, 400.0, 150.0), (3681795.8243822805, 208540182732463.34, 23976.19458507194)),
    ((7e7, 25.0, 3.3e-27, 3000.0, 100.0), (76238023.20909084, 8395996362979.754, 595521.9349601399)),
    ((1e6, 0.5, 7.3e-26, 5000.0, 30.0), (1002836.95, 21407076186867.75, 233567.62765516146)),   # reaches MAX_EXOBASE_RADII, the longest search
])
def test_the_exobase_search_gives_the_values_it_gave_before_it_was_bounded(args, expected):
    assert ts.exobase(*args) == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), 0.0, -3.0])
@pytest.mark.parametrize("slot", range(5))
def test_an_exobase_of_a_non_finite_or_non_positive_input_is_an_error_and_not_a_search_that_never_ends(slot, bad):
    """A NaN temperature made every comparison false and the height NaN, so the search ran for ever (GEN.841)."""
    args = [6.371e6, 9.81, 4.65e-26, 1100.0, 200.0]
    args[slot] = bad
    with pytest.raises(ValueError, match="finite positive"):
        ts.exobase(*args)


def test_the_exobase_search_stops_at_its_step_bound_and_returns_the_height_it_reached(monkeypatch):
    monkeypatch.setattr(ts, "MAX_EXOBASE_STEPS", 3)
    radius, density, scale = ts.exobase(6.371e6, 9.81, 4.65e-26, 1100.0, 200.0)
    assert 6.371e6 < radius < 6.371e6 * 1.01     # three steps of 0.25 scale heights up from the surface, and not the 6.79e6 m it closes at
    assert math.isfinite(density) and density > 0.0 and math.isfinite(scale) and scale > 0.0


def test_no_air_in_the_swept_range_comes_near_the_step_bound(monkeypatch):
    """The corners of the sweep that has the longest searches give the same exobase with a tenth of the bound as with the bound."""
    corners = [(6.4e6, 3.0, 4.65e-26, 5000.0, 100.0), (7e7, 0.5, 7.3e-26, 5000.0, 100.0), (1e6, 30.0, 2.7e-26, 5000.0, 100.0),
               (3e7, 0.5, 7.3e-26, 5000.0, 100.0), (6.4e6, 0.5, 7.3e-26, 1000.0, 30.0)]
    full = [ts.exobase(*c) for c in corners]
    monkeypatch.setattr(ts, "MAX_EXOBASE_STEPS", ts.MAX_EXOBASE_STEPS // 10)
    assert [ts.exobase(*c) for c in corners] == full


def test_the_cached_diffusion_parameter_and_the_precomputed_mean_change_no_value():
    uncached = ts.binary_diffusion_parameter.__wrapped__
    for a, b, t in [("H2", "N2", 150.0), ("He", "CO2", 400.0), ("H2O", "Ar", 250.0), ("CH4", "SO2", 300.0)]:
        assert ts.binary_diffusion_parameter(a, b, t) == uncached(a, b, t)
        assert ts.binary_diffusion_parameter(a, b, t) == ts.binary_diffusion_parameter(a, b, t)
    air = {"N2": 0.8, "H2": 0.15, "CO2": 0.05}
    mean = ts.mean_mass_kg(air)
    assert ts.diffusion_limited_flux("H2", air, 9.81, 200.0, mean) == ts.diffusion_limited_flux("H2", air, 9.81, 200.0)


def test_a_starved_xuv_flux_loses_next_to_nothing_instead_of_overflowing():
    """A planet lit by a black hole's faint glow gets an XUV flux near 1e-60: the boil-off test is made in logs, so no exp overflows."""
    rate = ts.hydro_based_mass_loss_g_s(20.0, 5.0, 300.0, 0.002, 1.0e-60)
    assert 0.0 <= rate < 1.0
