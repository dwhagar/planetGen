# tests/test_brown_dwarf.py

"""GEN.252, GEN.549: a brown dwarf as a star: its mass, how bright and hot it is at an age (the cooling track of Baraffe
et al. 2015, continued at [R0606]'s old-age exponent), the L, T or Y type that temperature makes
([R2374: Table 13]), and where those reach a page (a free-floating brown dwarf, a brown dwarf member of a
multi-star tree).

The 20 pc census's densities ([R2374: Table 15]) are test targets that show the cooling laws give a
plausible spread of temperatures. None of them is an input of the generator."""


import pathlib

import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.physics import brown_dwarf, brown_dwarf_types, constants, substellar_ceiling
from planetgen.util import draw
from planetgen.web import scatter_members_view
from tests import render_support


def test_the_masses_run_from_the_deuterium_limit_to_the_hydrogen_limit():
    low, high = brown_dwarf.mass_range_solar()
    assert low == pytest.approx(13 * constants.JUPITER_MASS_TO_KG / constants.SOLAR_MASS_TO_KG)
    assert high == pytest.approx(0.08)      # where the stars' mass function begins (GEN.936), 83.8 Jupiter masses
    with draw.bound(252):
        masses = [brown_dwarf.sample_mass_solar() for _ in range(4000)]
    assert low <= min(masses) and max(masses) <= high


def test_the_mass_function_is_the_censuss_power_law():
    """dN/dM ~ M^-0.6: the share below a mass M is (M^0.4 - low^0.4) / (high^0.4 - low^0.4)."""
    low, high = brown_dwarf.mass_range_solar()
    with draw.bound(252):
        masses = [brown_dwarf.sample_mass_solar() for _ in range(20000)]
    power = 1.0 - tuning.BROWN_DWARF_MASS_FUNCTION_SLOPE
    for mass in (0.02, 0.04, 0.06):
        expected = (mass ** power - low ** power) / (high ** power - low ** power)
        assert sum(m < mass for m in masses) / len(masses) == pytest.approx(expected, abs=0.015)


def test_the_track_gives_back_the_models_values_at_their_masses_and_ages():
    """[R0254]'s tabulated points: `(mass, age Gy, Teff K, radius Rsun)`."""
    for mass, age, temperature, radius in ((0.07, 10.0, 1626, 0.09), (0.05, 0.3, 2063, 0.109), (0.08, 3.0, 2343, 0.099),
                                           (0.04, 0.1, 2311, 0.13), (0.06, 1.0, 1584, 0.094)):
        assert brown_dwarf.effective_temperature_k(mass, age) == pytest.approx(temperature)
        assert brown_dwarf.radius_sol(mass, age) == pytest.approx(radius)


def test_the_temperature_falls_with_age_and_grows_with_mass():
    for mass in (0.0124, 0.03, 0.05, 0.065, 0.07, 0.0764):
        temperatures = [brown_dwarf.effective_temperature_k(mass, age) for age in (0.02, 0.1, 0.3, 1.0, 3.0, 10.0)]
        assert temperatures == sorted(temperatures, reverse=True), mass
    for age in (0.1, 1.0, 5.0, 10.0):
        temperatures = [brown_dwarf.effective_temperature_k(mass, age) for mass in (0.0124, 0.03, 0.05, 0.065, 0.0764)]
        assert temperatures == sorted(temperatures), age


def test_a_track_that_ends_before_10_gy_cools_at_burrows_exponent_and_one_that_reaches_it_holds():
    """[R0606], eq. 2: `Teff` falls as `t^-0.32` for an old brown dwarf. The 0.05 solar mass row ends at 0.625
    Gy; the 0.075 row reaches 10 Gy and holds the hydrogen-burning plateau."""
    assert brown_dwarf.effective_temperature_k(0.05, 5.0) == pytest.approx(1555.0 * (5.0 / 0.625) ** -0.32)
    assert brown_dwarf.effective_temperature_k(0.075, 13.0) == pytest.approx(2037.0)
    assert brown_dwarf.radius_sol(0.05, 5.0) == pytest.approx(0.099)


def test_the_track_is_continuous_in_mass_where_the_table_stops_holding_a_mass_for_longer():
    """GEN.1137: 0.0699 solar masses at 5 Gy was 1,220 K and 0.07 was 1,627 K, because the lighter mass cooled on at Burrows'
    exponent from its last age and the heavier one held the plateau of its table. The cooling past a body's last age now goes
    from the power law at the lightest mass of that age's table to the next heavier mass's own change (so the plateau is
    approached, not jumped to), and two masses a few percent apart give temperatures and radii within a few percent."""
    for age in (1.5, 5.0, 10.0, 13.0):
        for lighter, heavier in ((0.0699, 0.07), (0.0599, 0.06), (0.0499, 0.05), (0.0399, 0.04), (0.0299, 0.03)):
            assert brown_dwarf.effective_temperature_k(lighter, age) == pytest.approx(brown_dwarf.effective_temperature_k(heavier, age), rel=0.02), (lighter, age)
            assert brown_dwarf.radius_sol(lighter, age) == pytest.approx(brown_dwarf.radius_sol(heavier, age), rel=0.02), (lighter, age)
    # Across the whole table no mass step of a thousandth of a solar mass moves the temperature by more than 15 %.
    masses = [0.0125 + 0.001 * i for i in range(150)]
    for age in (0.5, 3.0, 10.0):
        temperatures = [brown_dwarf.effective_temperature_k(mass, age) for mass in masses]
        assert max(abs(b / a - 1.0) for a, b in zip(temperatures, temperatures[1:])) < 0.15, age


def test_a_young_brown_dwarf_near_the_limit_is_a_late_m_dwarf_no_hotter_than_the_peak_of_the_limit_mass():
    """The cap is the BHAC15 temperature of the hydrogen-burning limit (GEN.494): 2,950 K at the peak, no longer the 2,800 K of the rogue clip."""
    for mass in (0.065, 0.0764):
        temperature = brown_dwarf.effective_temperature_k(mass, 0.01)
        assert 2800.0 < temperature <= substellar_ceiling.HOTTEST_K
        assert brown_dwarf.spectral_type(temperature)[0] == "M"
    assert brown_dwarf.effective_temperature_k(0.0764, 0.01) == substellar_ceiling.HOTTEST_K


def test_no_brown_dwarf_is_hotter_than_the_limit_mass_of_its_age():
    """Below 0.07 solar masses a body cools for ever, so at any age it is no hotter than the limit mass at that age (2,050 K at 1 Gy, 1,630 K at 5 Gy)."""
    for age in (0.01, 0.1, 1.0, 5.0, 10.0, 13.0):
        for mass in (0.0124, 0.03, 0.05, 0.065, 0.0699):
            assert brown_dwarf.effective_temperature_k(mass, age) <= substellar_ceiling.hottest_k(age) + 1e-9
    assert brown_dwarf.effective_temperature_k(0.0699, 0.01) == pytest.approx(substellar_ceiling.ceiling_k(0.01), rel=0.01)  # the track of a body just under the limit mass is the limit mass's while it is young


def test_the_radius_is_about_jupiters_when_old_and_larger_when_young():
    jupiter = constants.JUPITER_RADIUS_KM * constants.KM_TO_M_FACTOR / constants.SOLAR_RADIUS_M
    for mass in (0.03, 0.05, 0.07, 0.0764):
        assert brown_dwarf.radius_sol(mass, 0.05) > brown_dwarf.radius_sol(mass, 1.0) >= brown_dwarf.radius_sol(mass, 10.0)
    for mass in (0.03, 0.05, 0.07, 0.0764):
        assert 0.7 * jupiter < brown_dwarf.radius_sol(mass, 10.0) < 1.5 * jupiter
    # The table loses a mass when it cools out of the models' range, the lightest (the deuterium limit) after 0.01 Gy, and the
    # radius it is held at is the young one (1.7 Jupiter radii at 13 Jupiter masses): the models give no old radius there.
    assert 0.7 * jupiter < brown_dwarf.radius_sol(0.0124, 10.0) < 1.8 * jupiter
    old = [brown_dwarf.radius_sol(mass, 8.0) for mass in (0.03, 0.05, 0.065, 0.075)]
    assert old == sorted(old, reverse=True)


def test_the_luminosity_is_stefan_boltzmann_of_the_tracks_radius_and_temperature():
    for mass, age in ((0.03, 3.0), (0.05, 1.0), (0.07, 8.0), (0.0764, 0.5)):
        found = brown_dwarf.state(mass, age)
        assert found["luminosity_sol"] == pytest.approx(brown_dwarf.luminosity_sol(mass, age))
        assert found["luminosity_sol"] == pytest.approx(
            found["radius_sol"] ** 2 * (found["temperature_k"] / tuning.SUN_EFFECTIVE_TEMPERATURE_K) ** 4)


def test_a_mass_outside_the_models_has_no_track():
    for mass in (0.005, 2.0):
        with pytest.raises(ValueError):
            brown_dwarf.effective_temperature_k(mass, 1.0)


def test_the_types_start_where_kirkpatrick_et_al_put_them():
    assert brown_dwarf_types.teff_of_type(0) == pytest.approx(2237.5)
    assert brown_dwarf_types.teff_of_type(10) == pytest.approx(1254.8, abs=0.5)
    assert brown_dwarf_types.teff_of_type(20) == pytest.approx(460.6, abs=0.5)
    # Y0 is inside the 371 to 510 K they measured.
    assert 371.0 < brown_dwarf_types.teff_of_type(20) < 510.0
    # The pieces meet where one range ends and the next begins.
    for number in (8.75, 14.75):
        below, above = brown_dwarf_types.teff_of_type(number - 1e-9), brown_dwarf_types.teff_of_type(number + 1e-9)
        assert below == pytest.approx(above, abs=0.5)


def test_the_map_read_backwards_gives_the_number_back():
    for number in (0.5, 4.0, 8.0, 10.0, 14.0, 18.0, 20.0, 21.5):
        assert brown_dwarf_types.number_of_teff(brown_dwarf_types.teff_of_type(number)) == pytest.approx(number, abs=1e-6)
    assert brown_dwarf_types.number_of_teff(2300.0) is None
    assert brown_dwarf_types.number_of_teff(100.0) == brown_dwarf_types.MAP.pieces[-1].s_max


@pytest.mark.parametrize("temperature, letters", [
    (3000.0, "M"), (2800.0, "ML"), (1600.0, "L"), (1200.0, "LT"), (900.0, "T"), (500.0, "TY"), (200.0, "Y"),
])
def test_the_type_is_drawn_about_the_mean_of_the_map_and_stays_with_the_temperature(temperature, letters):
    assert brown_dwarf.spectral_type(temperature)[0] in letters
    assert brown_dwarf.spectral_type(temperature) == brown_dwarf.spectral_type(temperature)


def test_dwarfs_of_one_temperature_spread_over_types_by_the_maps_scatter():
    """[R2374]: 134 K among L dwarfs, 79 K among T and Y. At 1600 K (L5) the types spread over about
    134 / 100 K per subtype, a few subtypes; the flat L/T transition spreads them over many more."""
    def spread(center):
        numbers = [brown_dwarf_types.drawn_number(center, draw.Stream(f"spread:{center}:{k}"))[0] for k in range(2000)]
        numbers = [n for n in numbers if n is not None]
        mean = sum(numbers) / len(numbers)
        return (sum((n - mean) ** 2 for n in numbers) / len(numbers)) ** 0.5

    assert 0.7 < spread(1600.0) < 1.6                                  # 134 K at about 100 K per subtype near L5
    assert spread(1220.0) > 2.0                                         # the nearly flat L/T transition
    assert 0.4 < spread(800.0) < 0.9                                    # 79 K at about 130 K per subtype near T8


def test_a_cooler_brown_dwarf_has_a_later_type_on_average():
    means = []
    for temperature in range(2000, 250, -250):
        numbers = [brown_dwarf_types.drawn_number(float(temperature), draw.Stream(f"mean:{temperature}:{k}"))[0]
                   for k in range(400)]
        numbers = [n for n in numbers if n is not None]
        means.append(sum(numbers) / len(numbers))
    assert means == sorted(means)


def test_the_generator_draws_the_type_through_the_map():
    """The wiring guard (GEN.760): `brown_dwarf.spectral_type` calls `brown_dwarf_types.drawn_number`, and nothing reads
    the old fits of the mean relation."""
    source = pathlib.Path(brown_dwarf.__file__).read_text(encoding="utf-8")
    assert "brown_dwarf_types.drawn_number(" in source
    assert not hasattr(tuning, "BROWN_DWARF_SPECTRAL_TYPE_FITS")
    assert not hasattr(brown_dwarf, "type_temperature_k") and not hasattr(brown_dwarf, "spectral_number")


def test_the_temperatures_of_a_census_spread_as_the_20_pc_census_does():
    """Table 15: the density per 150 K bin (1e-3 per pc3) at 600 to 1200 K is 2.80, 1.99, 1.72, 1.11 (complete there).
    The cooling laws, a constant birth rate over 10 Gy and the mass function give those shares within a factor 2."""
    low, high = brown_dwarf.mass_range_solar()
    census = {600: 2.80, 750: 1.99, 900: 1.72, 1050: 1.11, 1200: 1.95, 1350: 0.94, 1500: 0.81}
    counts = dict.fromkeys(census, 0)
    with draw.bound(252):
        for _ in range(60000):
            temperature = brown_dwarf.effective_temperature_k(brown_dwarf.sample_mass_solar(), draw.uniform(0.0, 10.0))
            for edge in counts:
                if edge <= temperature < edge + 150:
                    counts[edge] += 1
    shares = {edge: n / sum(counts.values()) for edge, n in counts.items()}
    observed = {edge: n / sum(census.values()) for edge, n in census.items()}
    for edge in (600, 750, 900, 1050):
        assert 0.5 <= shares[edge] / observed[edge] <= 2.0, (edge, shares[edge], observed[edge])


def test_the_hot_l_dwarfs_reach_the_censuss_share():
    census = {1500: 0.81, 1650: 0.78, 1800: 0.50, 1950: 0.72, 1200: 1.95, 1350: 0.94}
    counts = dict.fromkeys(census, 0)
    with draw.bound(252):
        for _ in range(60000):
            temperature = brown_dwarf.effective_temperature_k(brown_dwarf.sample_mass_solar(), draw.uniform(0.0, 10.0))
            for edge in counts:
                if edge <= temperature < edge + 150:
                    counts[edge] += 1
    for edge in (1500, 1650, 1800, 1950):
        share = counts[edge] / sum(counts.values())
        assert 0.5 <= share / (census[edge] / sum(census.values())) <= 2.0


def test_the_free_floating_brown_dwarf_has_the_laws_mass_radius_and_a_type():
    with draw.bound(252):
        dwarfs = [RoguePlanet(SystemConfig(), mass_bin="brown-dwarf") for _ in range(30)]
    for dwarf in dwarfs:
        solar = dwarf.mass_kg / constants.SOLAR_MASS_TO_KG
        assert dwarf.radius_km == pytest.approx(brown_dwarf.radius_sol(solar, dwarf.age_gy) * constants.SOLAR_RADIUS_M
                                                / constants.KM_TO_M_FACTOR)
        assert dwarf.effective_temperature_k == pytest.approx(brown_dwarf.effective_temperature_k(solar, dwarf.age_gy),
                                                              rel=1e-6)
        assert dwarf.spectral_type == brown_dwarf.spectral_type(dwarf.effective_temperature_k)
        section = render_support.phenomenon_section(dwarf)
        assert f"{dwarf.spectral_type} brown dwarf" in render_support.prose_text(section)
        assert f"type {dwarf.spectral_type}" in render_support.facts_text(section)
    assert RoguePlanet(SystemConfig(), mass_bin="terrestrial").spectral_type is None


def test_a_brown_dwarf_member_has_a_type_and_the_radius_the_laws_give():
    with draw.bound(252):
        leaf = multistar._star_leaf(0.05, 4.0e9, draw)
    found = brown_dwarf.state(0.05, 4.0)
    assert leaf.kind == "brown-dwarf" and leaf.state == found["spectral_type"]
    assert leaf.radius_au == pytest.approx(found["radius_sol"] * constants.SOLAR_RADIUS_M / constants.AU_M)
    assert leaf.state[0] in "LTY"


def test_the_members_page_names_a_brown_dwarfs_type():
    tree = {"members": [{"member_no": 0, "kind": "main-sequence", "state": "V", "mass_solar": 0.6},
                        {"member_no": 1, "kind": "brown-dwarf", "state": "T5", "mass_solar": 0.05},
                        {"member_no": 2, "kind": "brown-dwarf", "state": None, "mass_solar": 0.05}],
            "nodes": []}
    labels = [row[1] for row in scatter_members_view.members_view(tree)["members"]]
    assert labels == ["Main-sequence star", "Brown dwarf (T5)", "Brown dwarf"]


def _ages(population=None, age_law=None, count=60):
    config = SystemConfig()
    config.POPULATION = population
    config.AGE_LAW = age_law
    with draw.bound(553):
        return [RoguePlanet(config, mass_bin="brown-dwarf").age_gy for _ in range(count)]


@pytest.mark.parametrize("population", ["young", "intermediate", "old", "bulge", "halo"])
def test_a_free_brown_dwarf_takes_the_age_window_of_its_population(population):
    low, high = tuning.STELLAR_POPULATION_AGE_RANGES_GY[population]
    ages = _ages(population)
    assert all(low <= age <= high for age in ages)
    assert max(ages) - min(ages) > 0.3 * (high - low)


def test_a_free_brown_dwarf_made_without_a_place_has_the_whole_disks_window():
    low, high = tuning.STAR_FORMATION_AGE_RANGE_GY
    ages = _ages()
    assert all(low <= age <= high for age in ages) and max(ages) > 0.8 * high


def test_a_free_brown_dwarfs_age_follows_the_populations_age_law():
    from planetgen.galaxy.age_law import AgeLaw
    window = tuning.STELLAR_POPULATION_AGE_RANGES_GY["old"]
    front = AgeLaw(ages=(window[0], window[1]), rates=(1.0, 0.0))     # all the rate at the young edge
    ages = _ages("old", front, count=200)
    assert sum(ages) / len(ages) < window[0] + 0.4 * (window[1] - window[0])
    flat = _ages("old", None, count=200)
    assert sum(flat) / len(flat) > window[0] + 0.4 * (window[1] - window[0])


def test_a_galaxy_sectors_brown_dwarfs_take_the_populations_of_its_place(monkeypatch):
    """GEN.553: the sector fill gives each brown dwarf a population from the place's mix, and so an age inside that
    population's window; the other kinds are built from the bare config as before."""
    import math
    from types import SimpleNamespace
    from planetgen.galaxy import profile
    from planetgen.galaxy.density import build_galaxy_shape
    from planetgen.galaxy.skeleton import expected_system_count_at_density_1
    from planetgen.generation import bright_stars, run_sector
    from planetgen.physics.units import ly_to_pc

    shape = build_galaxy_shape(
        disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
        arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)
    edge_pc = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
    here = profile.profile_at((8128.0, 0.0, 0.0), shape, edge_pc, expected_system_count_at_density_1(
        tuning.DEFAULT_SECTOR_EDGE_LY))
    fill = bright_stars.FillContext(here)
    made = []

    def build(config, _distance):
        dwarf = RoguePlanet(config, mass_bin="brown-dwarf")
        made.append((config.POPULATION, dwarf.age_gy))
        return dwarf

    sector = SimpleNamespace(entries=[SimpleNamespace(star_system=SimpleNamespace(stars=[0] * 400))], stage_base=None)
    monkeypatch.setattr(run_sector, "SECTOR_PHENOMENON_KINDS", (("brown-dwarf", "rogue-planet", build),))

    def not_placed(*args, **kwargs):
        raise ValueError("no room")      # the draw is what is read; placing the dwarf is not needed

    sector.add_phenomenon = not_placed
    monkeypatch.setattr(run_sector, "add_star_hosted_nebulae", lambda *args, **kwargs: [])
    with draw.bound(553):
        run_sector.generate_sector_phenomena(sector, None, fill=fill, scattered=False)
    assert len(made) >= 30
    for population, age in made:
        low, high = tuning.STELLAR_POPULATION_AGE_RANGES_GY[population]
        assert low <= age <= high
    assert len({population for population, _age in made}) > 1
