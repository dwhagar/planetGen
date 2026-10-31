# tests/test_biogenic_oxygen.py

"""GEN.855: `physics.biogenic_oxygen`, the oxygen oxygenic life puts in the air, and its wiring in `generation.life`."""

import pytest

from planetgen import tuning
from planetgen.generation import evolution
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import atmosphere, biogenic_oxygen as bo, gases, life_classes

from tests.fuzz_support import deterministic_entropy


def _worlds(seeds):
    for seed in seeds:
        config = SystemConfig()
        config.STAR_TYPE = "G2V"
        config.PLANETS = True
        config.HABITABLE_WORLD = True
        with deterministic_entropy(seed):
            system = StarSystem(system_config=config)
        for planet in system.planets:
            if planet.body_type != "g" and planet.body_type != "a":
                yield planet


OXYGENIC_WORLDS_WANTED = 3
"""int: How many worlds of oxygenic life the sample must hold for the tests that read them."""

FIRST_SEEDS = 12
"""int: The seeds always drawn; more follow until `OXYGENIC_WORLDS_WANTED` oxygenic worlds are held (one seed in eight gives one)."""

LAST_SEED = 80
"""int: The seed after which the sample stops growing (the tests then fail on the missing oxygenic worlds)."""


def _makes_oxygen(world):
    """Whether the world has oxygenic life that can build the oxygen (a delay): a request for a habitable world only raises the odds (GEN.1130), so
    a sample of them holds worlds that cannot, and the tests read the ones that can."""
    return bool(world.life_class_id and evolution.releases_oxygen(world.life_class_id)) and bo.profile_of(world).delay_gyr is not None


@pytest.fixture(scope="module")
def worlds():
    held = list(_worlds(range(0, FIRST_SEEDS)))
    seed = FIRST_SEEDS
    while seed < LAST_SEED and sum(1 for w in held if _makes_oxygen(w)) < OXYGENIC_WORLDS_WANTED:
        held += list(_worlds([seed]))
        seed += 1
    return held


def test_earth_numbers_balance_at_the_oxygen_of_earths_air():
    assert bo.WEATHERING_PER_YR * bo.EARTH_O2_MOL == pytest.approx(bo.EARTH_BURIAL_MOL_YR - bo.EARTH_REDUCTANT_MOL_YR)
    assert bo.EARTH_STORE_MOL / (bo.EARTH_BURIAL_MOL_YR - bo.EARTH_REDUCTANT_MOL_YR) / 1e9 == pytest.approx(0.714, abs=0.01)


def test_a_more_reduced_mantle_puts_out_more_reductants_inside_the_range():
    assert bo.redox_factor(None) == 1.0
    assert bo.redox_factor(bo.EARTH_MANTLE_DELTA_IW) == pytest.approx(1.0)
    assert bo.redox_factor(bo.EARTH_MANTLE_DELTA_IW - 1.0) == pytest.approx(10 ** bo.REDOX_DECADE_PER_UNIT)
    low, high = bo.REDOX_FACTOR_RANGE
    assert bo.redox_factor(-50.0) == high and bo.redox_factor(50.0) == low


def test_the_sea_makes_nothing_outside_its_liquid_range_and_most_near_the_reference():
    cold, _, _, hot = bo.OCEAN_SPAN_K
    assert bo.temperature_factor(None) == 0.0
    assert bo.temperature_factor(cold) == 0.0 and bo.temperature_factor(hot) == 0.0
    assert bo.temperature_factor(bo.REFERENCE_TEMPERATURE_K) == pytest.approx(1.0)
    assert 0.0 < bo.temperature_factor(295.0)


def test_burning_takes_the_oxygen_the_reduced_gases_want_and_leaves_their_oxides():
    moles = {"N2": 10.0, "CH4": 1.0, "H2": 2.0, "CO": 1.0, "O2": 0.5}
    burnt, taken = bo.burn(moles)
    assert taken == pytest.approx(1.0 * 2.0 + 2.0 * 0.5 + 1.0 * 0.5)
    assert set(burnt) == {"N2", "CO2", "O2"} and burnt["CO2"] == pytest.approx(2.0) and burnt["O2"] == 0.5


def test_the_fire_limit_is_the_share_of_the_molecules_the_mixture_allows():
    from planetgen.physics import mixture_oxygen
    moles = {"N2": 100.0}
    limit = bo._fire_limit_mol(moles)
    assert limit / (limit + 100.0) == pytest.approx(mixture_oxygen.BIOTIC_O2_MAX)


def test_the_oxygen_comes_after_the_delay_and_not_before():
    profile = bo.Profile(2.0, 1.0, 5.0, 3.0, 2.0e4, 4.0)
    assert bo.oxygen_mol(profile, 3.99) == 0.0 and bo.oxygen_mol(profile, 4.0) == 3.0
    assert bo.oxygen_mol(bo.NONE, 100.0) == 0.0


def test_a_world_without_a_sunlit_sea_makes_no_oxygen(worlds):
    seen = 0
    for world in worlds:
        if getattr(world, "hydrosphere", None) not in bo.LIVING_HYDROSPHERES:
            assert bo.productivity(world, 5.1e14) == 0.0
            seen += 1
    assert seen > 3


def test_only_oxygenic_life_holds_oxygen_it_made_and_only_in_a_sunlit_sea(worlds):
    seen = 0
    for world in worlds:
        if not world.life_class_id:
            continue
        o2_kpa = atmosphere.partial_kpa(world, "o2")
        if o2_kpa > 0.2 and not evolution.releases_oxygen(world.life_class_id):
            # an abiotic oxygen air is allowed, but an oxygenic world is the only one that grows its own
            assert world.life_age_complex_cells_gyr is not None
        if evolution.releases_oxygen(world.life_class_id):
            assert world.hydrosphere in bo.LIVING_HYDROSPHERES, (world.name, world.hydrosphere)
            seen += 1
    assert seen >= 1


def test_the_stored_air_stays_one_air_after_the_oxygen_goes_in(worlds):
    for world in worlds:
        if world.life_class_id and evolution.releases_oxygen(world.life_class_id) and atmosphere.partial_kpa(world, "o2") > 1.0:
            partials = atmosphere.gas_partials_pa(world)
            assert sum(partials.values()) == pytest.approx(world.atmospheric_pressure, rel=2e-3), world.name
            assert not (set(partials) & set(bo.REDUCED_GASES)) or partials["O2"] < 1e-3 * sum(partials.values())
            assert all(gases.gas(f) for f in partials)
            return
    pytest.skip("no oxygen made in these seeds")


def test_an_oxygenic_class_eats_light():
    oxygenic = [row for row in life_classes.CLASSES if row.releases_oxygen]
    assert oxygenic
    for row in oxygenic:
        assert row.light and row.metabolism == "oxygenic", row.id


def test_an_oxygenic_world_past_its_delay_holds_oxygen_and_one_before_it_holds_none(worlds):
    """GEN.873: the air follows the life's age. A world whose photosynthesis is older than the delay holds the oxygen; one younger holds what
    its air had before (none made by life)."""
    past = before = 0
    for world in worlds:
        if not (world.life_class_id and evolution.releases_oxygen(world.life_class_id)):
            continue
        profile = bo.profile_of(world)
        age = world.star.age - world.life_age_photosynthesis_gyr
        if profile.delay_gyr is None:
            continue
        o2 = atmosphere.partial_kpa(world, "o2")
        if age >= profile.delay_gyr:
            assert o2 > 0.5, (world.name, o2, age, profile.delay_gyr)
            past += 1
        else:
            assert o2 < 0.2, (world.name, o2, age, profile.delay_gyr)
            before += 1
    assert past + before >= 1
