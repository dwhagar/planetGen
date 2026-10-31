"""
GEN.236: the sector profile's kind factors (`SectorProfile.kind_factors`) and the planets a metal-poor star makes. The rate
per star of each rare kind stays what it was; the place scales it (supernova remnants by the gas, molecular clouds and
nebulae by the gas and the young stars, rogue planets and comets by the planets the place's [Fe/H] allows), and a
metal-poor star's disk makes fewer planets and almost no giants.
"""

import math
import types
from types import SimpleNamespace

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import nebula_field, profile, superbubble_field
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.skeleton import expected_system_count_at_density_1
from planetgen.generation import nebula_contents, phenomenon_scatter, run_sector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, metallicity, planet_architecture
from planetgen.physics.units import ly_to_pc
from planetgen.util import draw

pytestmark = pytest.mark.usefixtures("large_spiral_fallback")

SHAPE = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
    arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)
EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
E_VALUE = expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY)
SUN = (8128.0, 0.0, 0.0)
INNER_DISK = (3000.0, 0.0, 0.0)
HALO = (20000.0, 3000.0, 15000.0)
ABOVE_THE_DISK = (8128.0, 0.0, 1500.0)
EVERYWHERE_ONE = ("planetary-nebula", "brown-dwarf")


def _at(position):
    return profile.profile_at(position, SHAPE, EDGE_PC, E_VALUE)


# --- The occurrence functions ---------------------------------------------------------------------------------------


def test_planet_occurrence_is_one_at_and_above_the_floor_and_falls_below_it():
    floor = metallicity.floor_fe_h()
    assert metallicity.planet_occurrence_factor(floor) == 1.0
    assert metallicity.planet_occurrence_factor(0.4) == 1.0
    assert metallicity.planet_occurrence_factor(floor - 1.0) == pytest.approx(
        10.0 ** -tuning.PLANET_OCCURRENCE_METALLICITY_EXPONENT)
    values = [metallicity.planet_occurrence_factor(floor - 0.2 * step) for step in range(10)]
    assert values == sorted(values, reverse=True)


def test_giant_occurrence_is_one_at_and_above_the_floor_and_falls_faster_than_planets():
    floor = metallicity.floor_fe_h()
    assert metallicity.giant_occurrence_factor(floor) == 1.0
    assert metallicity.giant_occurrence_factor(0.3) == 1.0
    for fe_h in (-0.8, -1.2, -2.0):
        assert metallicity.giant_occurrence_factor(fe_h) < metallicity.planet_occurrence_factor(fe_h) < 1.0
    assert metallicity.giant_occurrence_factor(-1.6) < 0.01


def test_the_expected_factor_of_a_population_is_the_gauss_hermite_mean_of_its_distribution():
    peaks = ((1.0, -0.4, 0.15),)
    assert metallicity.expected_factor(lambda fe_h: fe_h, peaks) == pytest.approx(-0.4)
    assert metallicity.expected_factor(lambda fe_h: (fe_h + 0.4) ** 2, peaks) == pytest.approx(0.15 ** 2)
    two = ((0.25, -1.5, 0.1), (0.75, -0.5, 0.1))
    assert metallicity.expected_factor(lambda fe_h: fe_h, two) == pytest.approx(0.25 * -1.5 + 0.75 * -0.5)
    assert metallicity.expected_factor(lambda fe_h: 1.0, ()) == 1.0


def test_the_expected_factor_clamps_a_draw_as_the_sampler_does():
    low, high = tuning.METALLICITY_RANGE_FE_H
    assert metallicity.expected_factor(lambda fe_h: fe_h, ((1.0, 5.0, 0.01),)) == pytest.approx(high)
    assert metallicity.expected_factor(lambda fe_h: fe_h, ((1.0, -9.0, 0.01),)) == pytest.approx(low)


# --- The profile's factors ------------------------------------------------------------------------------------------


@pytest.mark.parametrize("position", [SUN, INNER_DISK, HALO, ABOVE_THE_DISK, (0.0, 0.0, 0.0)])
def test_planetary_nebulae_and_brown_dwarfs_are_not_scaled_by_place(position):
    factors = _at(position).kind_factors
    for kind in EVERYWHERE_ONE:
        assert factors[kind] == 1.0


@pytest.mark.parametrize("position", [SUN, INNER_DISK, HALO, ABOVE_THE_DISK, (0.0, 0.0, 0.0)])
def test_the_black_hole_factor_is_the_failed_share_the_places_stars_make(position):
    """GEN.563: not 1 but the share of the place's massive stars that failed to explode and made their remnant: below
    the galaxy's metal-free share, and the halo's metal-poor stars make more than the disk's."""
    factor = _at(position).kind_factors["black-hole"]
    assert 0.0 < factor <= SHAPE.structure.black_hole_metal_free_share / (1.0 - 0.00637) * 1.001


@pytest.mark.parametrize("position", [SUN, INNER_DISK, HALO, ABOVE_THE_DISK, (0.0, 0.0, 0.0)])
def test_the_neutron_star_factor_is_what_the_black_holes_leave_less_the_unbound_share(position):
    """GEN.554: the massive stars that did not fail to explode, bound: below what the solar share leaves less the
    galaxy's unbound share."""
    structure = SHAPE.structure
    factor = _at(position).kind_factors["neutron-star"]
    assert 0.0 < factor <= (1.0 - structure.black_hole_solar_share) * (1.0 - structure.neutron_star_unbound_share) / (1.0 - 0.00637) * 1.001


def test_a_metal_rich_place_makes_more_neutron_stars_than_a_metal_poor_one():
    assert _at(SUN).kind_factors["neutron-star"] > _at(HALO).kind_factors["neutron-star"]


def test_a_metal_poor_place_makes_more_black_holes_than_a_metal_rich_one():
    assert _at(HALO).kind_factors["black-hole"] > _at(SUN).kind_factors["black-hole"]
    assert _at(SUN).kind_factors["black-hole"] >= _at(INNER_DISK).kind_factors["black-hole"]


@pytest.mark.parametrize("position", [SUN, INNER_DISK, HALO, ABOVE_THE_DISK, (0.0, 0.0, 0.0)])
def test_every_factor_is_between_zero_and_one(position):
    factors = _at(position).kind_factors
    assert factors
    assert all(0.0 <= value <= 1.0 for value in factors.values()), factors


def test_a_profile_is_still_a_pure_function_with_its_factors():
    assert _at(SUN).kind_factors == _at(SUN).kind_factors
    assert _at(SUN) == _at(SUN)


def test_the_supernova_remnant_factor_is_the_gas_fraction():
    here = _at(SUN)
    assert here.kind_factors["supernova-remnant"] == pytest.approx(here.gas.cold_fraction)
    assert 0.0 < here.kind_factors["supernova-remnant"] <= 1.0
    assert _at(INNER_DISK).kind_factors["supernova-remnant"] == 1.0      # the inner disk is gas-rich: capped at 1


@pytest.mark.parametrize("position", [HALO, ABOVE_THE_DISK])
def test_no_supernova_remnant_cloud_or_nebula_factor_in_the_halo_or_above_the_disk(position):
    factors = _at(position).kind_factors
    assert factors["supernova-remnant"] == pytest.approx(0.0, abs=1e-6)
    assert factors["molecular-cloud"] == pytest.approx(0.0, abs=1e-6)
    assert factors["nebula"] == pytest.approx(0.0, abs=1e-6)


def test_a_superbubble_zeroes_the_clouds_and_nebulae_but_not_the_remnants(monkeypatch):
    open_sky = _at(SUN)
    assert open_sky.kind_factors["nebula"] > 0.0
    bubble = SimpleNamespace(ambient_density_cm3=0.003)
    monkeypatch.setattr(superbubble_field, "containing", lambda seed, shape, position: bubble)
    inside = profile.profile_at(SUN, SHAPE, EDGE_PC, E_VALUE, galaxy_seed=bytes(range(16)))
    assert inside.gas.superbubble is bubble
    assert inside.kind_factors["molecular-cloud"] == 0.0
    assert inside.kind_factors["nebula"] == 0.0
    assert inside.kind_factors["supernova-remnant"] == open_sky.kind_factors["supernova-remnant"]


def test_the_cloud_factor_is_the_cold_gas_times_the_young_share_against_the_solar_circles():
    here = _at(SUN)
    reference = profile._reference_young_share(SHAPE)
    expected = here.gas.cold_fraction * min(1.0, here.populations["young"].share / reference)
    assert here.kind_factors["molecular-cloud"] == pytest.approx(expected)
    assert here.kind_factors["nebula"] == pytest.approx(here.gas.cold_fraction)
    assert 0.0 < here.kind_factors["molecular-cloud"] < here.kind_factors["nebula"]


def test_rogue_planets_and_comets_follow_the_planet_occurrence_of_the_places_metallicity():
    for position in (SUN, INNER_DISK, HALO):
        factors = _at(position).kind_factors
        assert factors["rogue-planet"] == factors["comet"]
    assert _at(SUN).kind_factors["rogue-planet"] > 0.9
    assert _at(HALO).kind_factors["rogue-planet"] < 0.1
    here = _at(HALO)
    expected = metallicity.expected_factor(metallicity.planet_occurrence_factor, here.populations["halo"].fe_h)
    assert here.kind_factors["rogue-planet"] == pytest.approx(expected, rel=0.02)


# --- The vectorised gas factor ----------------------------------------------------------------------------------------


def test_the_array_gas_factor_matches_the_one_point_form():
    points = [SUN, INNER_DISK, HALO, ABOVE_THE_DISK, (0.0, 0.0, 0.0), (-9000.0, 2000.0, 500.0), (5000.0, 3000.0, 100.0)]
    x, y, z = (np.array([p[axis] for p in points]) for axis in range(3))
    array = nebula_field.gas_factor_array(x, y, z, SHAPE)
    scalar = [nebula_field.gas_factor(p, SHAPE) for p in points]
    assert array == pytest.approx(scalar, rel=1e-9, abs=1e-12)


def test_the_array_gas_factor_broadcasts_a_band_height_against_a_ring_of_bins():
    x, y = np.array([[8128.0, 0.0], [0.0, 8128.0]]), np.array([[0.0, 8128.0], [8128.0, 0.0]])
    z = np.array([[0.0], [900.0]])
    out = nebula_field.gas_factor_array(x, y, z, SHAPE)
    assert out.shape == (2, 2)
    assert out[0, 0] == pytest.approx(nebula_field.gas_factor((8128.0, 0.0, 0.0), SHAPE))
    assert out[1, 0] == pytest.approx(nebula_field.gas_factor((0.0, 8128.0, 900.0), SHAPE))
    assert out[1, 1] == pytest.approx(nebula_field.gas_factor((8128.0, 0.0, 900.0), SHAPE))


# --- The scatter's supernova remnants ---------------------------------------------------------------------------------


def test_the_scatters_gas_weight_is_zero_in_the_halo_and_the_gas_fraction_in_the_disk():
    x, y, z = np.array([HALO[0], SUN[0]]), np.array([HALO[1], SUN[1]]), np.array([HALO[2], SUN[2]])
    weight = phenomenon_scatter._gas_weight(x, y, z, SHAPE)
    assert weight[0] == pytest.approx(0.0, abs=1e-6)
    assert weight[1] == pytest.approx(_at(SUN).kind_factors["supernova-remnant"])


def test_only_the_supernova_remnants_are_gas_weighted_in_the_scatter():
    assert phenomenon_scatter.GAS_WEIGHTED_KINDS == ("supernova-remnant",)


def test_the_scatters_expectation_for_a_remnant_carries_the_gas_and_a_black_hole_does_not():
    assert phenomenon_scatter._expectation_factor(SHAPE, "supernova-remnant", None, SUN[0], 5000.0) < 1e-6
    disk_remnant = phenomenon_scatter._expectation_factor(SHAPE, "supernova-remnant", None, SUN[0], 0.0)
    assert disk_remnant == pytest.approx(
        min(1.0, nebula_field.gas_factor(SUN, SHAPE))
        * phenomenon_scatter.remnant_distribution.vertical_factor("supernova-remnant", 0.0, SHAPE.disk_scale_height_pc)
        * phenomenon_scatter.remnant_distribution.radial_factor("supernova-remnant", SUN[0]))


# --- The fill -----------------------------------------------------------------------------------------------------------


def _seeded_sector(count=7):
    from planetgen.galaxy.sector import SpaceSector
    sector = SpaceSector("Kind Factor Test", edge_ly=100.0)
    for _ in range(count):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = False
        sector.add_system(StarSystem(system_config=cfg), system_config=cfg)
    return sector


def test_a_fill_multiplies_each_kinds_mean_by_its_factor_and_leaves_a_kind_without_one_alone(monkeypatch):
    sector = _seeded_sector()
    stars = run_sector.sector_star_count(sector)
    means = {}
    current = {}

    def fake(mean):
        means[current["kind"]] = mean
        return 0

    monkeypatch.setattr(run_sector, "_sample_poisson_count", fake)
    original = run_sector.fill_stages.stage

    def tracking_stage(sector_, name, sub=None):
        current["kind"] = sub
        return original(sector_, name, sub)

    monkeypatch.setattr(run_sector.fill_stages, "stage", tracking_stage)
    factors = {"rogue-planet": 0.25, "comet": 0.0, "supernova-remnant": 0.5}
    run_sector.generate_sector_phenomena(sector, SimpleNamespace(), kind_factors=factors)
    rate = tuning.phenomenon_rate_per_star
    assert means["rogue-planet"] == pytest.approx(rate("rogue-planet") * stars * 0.25)
    assert means["comet"] == 0.0
    assert means["brown-dwarf"] == pytest.approx(rate("brown-dwarf") * stars)


def test_a_fill_with_no_factors_draws_the_same_as_one_with_every_factor_one():
    def kinds(factors):
        with draw.bound("kind-factors-same"):
            sector = _seeded_sector(5)
            entries = run_sector.generate_sector_phenomena(sector, SimpleNamespace(), kind_factors=factors)
        return sorted(entry.phenomenon_type for entry in entries)
    assert kinds(None) == kinds({"rogue-planet": 1.0, "comet": 1.0, "nebula": 1.0})


def test_a_factor_of_zero_rolls_no_rogue_planet_or_comet(monkeypatch):
    monkeypatch.setitem(tuning.PHENOMENON_RATE_SCALE, "brown-dwarf", 0.0)    # saved as a rogue planet of its own bin
    with draw.bound("kind-factors-zero"):
        sector = _seeded_sector(7)
        entries = run_sector.generate_sector_phenomena(
            sector, SimpleNamespace(), kind_factors={"rogue-planet": 0.0, "comet": 0.0})
    assert not [e for e in entries if e.phenomenon_type in ("rogue-planet", "comet")]


def test_the_hosted_nebula_chance_is_scaled_by_the_nebula_factor(monkeypatch):
    class Star:
        type = "O1V"
        age = 0.0005
        temperature = 30000.0

        def __init__(self, number):
            self.mass = 2.0e31 * (1.0 + number * 1e-6)  # each star has its own roll (GEN.320: drawn from the star's own seed)

    class Entry:
        made = 0

        def __init__(self):
            Entry.made += 1
            star = Star(Entry.made)
            self.star_system = types.SimpleNamespace(stars=[star], placed_star=star, name="x")
            self.position = (0.0, 0.0, 0.0)

    class Sector:
        name = "s"

        def __init__(self, count):
            self.entries = [Entry() for _ in range(count)]
            self.added = 0

        def add_phenomenon(self, nebula, kind, position=None):
            self.added += 1
            return nebula

    monkeypatch.setattr(nebula_contents, "Nebula", lambda config, nebula_class=None, host=None, age_myr=None: types.SimpleNamespace(name="n", nebula_class="C"))
    monkeypatch.setattr(nebula_contents, "choose_weighted_class", lambda classes, label: classes[0])
    rates = {}
    for factor in (1.0, 0.4, 0.0):
        draw.set_run_seed(236)
        sector = Sector(4000)
        run_sector.add_star_hosted_nebulae(sector, types.SimpleNamespace(), chance_factor=factor)
        rates[factor] = sector.added / 4000
    assert abs(rates[1.0] - 0.25) < 0.03
    assert abs(rates[0.4] - 0.10) < 0.03
    assert rates[0.0] == 0.0


# --- The planet cull --------------------------------------------------------------------------------------------------


def _architectures(fe_h, count=60):
    """The architectures of `count` Sun-like stars of [Fe/H] `fe_h` (Z / Z0 = 10^fe_h), each from its own stream."""
    return [planet_architecture.build(constants.SOLAR_MASS_TO_KG, 2.34, 10.0 ** fe_h, fe_h, 30.0, stream=draw.Stream(f"cull:{fe_h}:{i}"))
            for i in range(count)]


@pytest.mark.parametrize("fe_h", [-1.0, -1.6, -2.4])
def test_a_metal_poor_star_carries_fewer_planets_than_a_solar_one(fe_h):
    poor = sum(len(a.planets) for a in _architectures(fe_h))
    solar = sum(len(a.planets) for a in _architectures(0.0))
    assert poor < 0.6 * solar


@pytest.mark.parametrize("fe_h", [-1.0, -1.6, -2.4])
def test_a_metal_poor_star_makes_almost_no_giants(fe_h):
    giants = sum(p.runaway for a in _architectures(fe_h, 100) for p in a.planets)
    solar_giants = sum(p.runaway for a in _architectures(0.0, 100) for p in a.planets)
    assert giants <= 0.1 * solar_giants


def test_no_halo_system_carries_more_planets_than_a_disk_one_does_on_average():
    counts = []
    for index in range(20):
        with draw.bound(f"halo-cull:{index}"):
            cfg = SystemConfig()
            cfg.STAR_TYPE = "G2V"
            cfg.POPULATION = "halo"
            cfg.BINARY_SYSTEM = False
            system = StarSystem(system_config=cfg)
        counts.append(len(system.planets))
    assert sum(counts) / len(counts) < 0.5 * 6.9      # a disk G dwarf holds about 6.9


def test_a_close_binary_proxy_takes_its_primarys_metallicity():
    from planetgen.generation.binary import BinaryStarProxy
    from planetgen.generation.star import Star
    config = SystemConfig()
    primary, secondary = Star(config, star_type="G2V", metallicity_fe_h=-1.4), Star(config, star_type="M5V", metallicity_fe_h=-0.2)
    proxy = BinaryStarProxy(config, primary, secondary)
    assert proxy._primary in (primary, secondary)      # the heavier of the two
    assert proxy.metallicity_fe_h == proxy._primary.metallicity_fe_h
    assert proxy.metallicity_fe_h in (-1.4, -0.2)


# --- The Type Ia share and the old share (GEN.319) --------------------------------------------------------------------


def _populations(young_share):
    """Population profiles that hold `young_share` of the stars as young and the rest as old."""
    return {name: types.SimpleNamespace(
        share=young_share if name == "young" else (1.0 - young_share if name == "old" else 0.0), fe_h=((1.0, 0.0, 0.1),))
        for name in profile.POPULATIONS}


def test_the_type_ia_share_is_the_solar_circles_chance_at_the_solar_circles_mix():
    reference = profile._reference_young_share(SHAPE)
    assert profile.type_ia_share_at(_populations(reference), SHAPE) == pytest.approx(
        tuning.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE)


def test_the_type_ia_share_falls_as_the_young_share_rises_and_is_one_with_none():
    reference = profile._reference_young_share(SHAPE)
    shares = [profile.type_ia_share_at(_populations(reference * factor), SHAPE) for factor in (0.0, 0.25, 1.0, 3.0)]
    assert shares[0] == 1.0
    assert shares == sorted(shares, reverse=True) and shares[-1] > 0.0
    assert profile.type_ia_share_at(_populations(1.0), SHAPE) == 0.0


@pytest.mark.parametrize("position", [HALO, ABOVE_THE_DISK])
def test_a_place_without_young_stars_has_only_type_ia_remnants(position):
    here = _at(position)
    assert here.populations["young"].share < 1e-3
    assert here.type_ia_share > 0.99


def test_the_galactic_centre_with_a_sliver_of_young_stars_is_nearly_all_type_ia():
    assert _at((0.0, 0.0, 0.0)).type_ia_share > 0.8


def test_the_disk_has_both_kinds_and_the_report_carries_the_share():
    from planetgen.galaxy import profile_report
    for position in (SUN, INNER_DISK):
        here = _at(position)
        assert 0.0 < here.type_ia_share < 1.0
        assert profile_report.profile_report(here)["type_ia_share"] == here.type_ia_share
    assert _at(SUN) == _at(SUN)


def test_the_planetary_nebula_factor_is_one_in_the_disk_within_a_few_percent_and_in_every_older_place():
    for position in (SUN, INNER_DISK, (5000.0, 5000.0, 0.0), (0.0, 12000.0, 0.0)):
        assert 0.95 < _at(position).kind_factors["planetary-nebula"] <= 1.0
    for position in (HALO, ABOVE_THE_DISK):
        assert _at(position).kind_factors["planetary-nebula"] == 1.0
    assert profile.kind_factors_at(_populations(0.9), 1.0, 1.0, SHAPE)["planetary-nebula"] < 0.12


def test_the_disks_other_factors_are_what_they_were():
    here = _at(SUN)
    assert here.kind_factors["supernova-remnant"] == pytest.approx(here.gas.cold_fraction)
    assert here.kind_factors["brown-dwarf"] == 1.0


def _remnants(share, count=300):
    made = []
    with draw.bound(f"type-ia-{share}"):
        for _ in range(count):
            config = SystemConfig()
            config.TYPE_IA_SHARE = share
            made.append(run_sector.SupernovaRemnant(config))
    return made


def test_a_remnant_is_type_ia_with_the_places_share():
    assert {remnant.progenitor_type for remnant in _remnants(1.0)} == {"Type Ia"}
    assert {remnant.progenitor_type for remnant in _remnants(0.0)} == {"core-collapse"}
    solar = sum(remnant.progenitor_type == "Type Ia" for remnant in _remnants(None)) / 300
    assert solar == pytest.approx(tuning.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE, abs=0.08)
    assert sum(remnant.progenitor_type == "Type Ia" for remnant in _remnants(0.8)) / 300 == pytest.approx(0.8, abs=0.08)


def test_a_halo_remnant_has_no_neutron_star_or_black_hole():
    assert all(remnant.compact_remnant is None for remnant in _remnants(_at(HALO).type_ia_share))


def test_the_scattered_remnants_and_the_fill_take_the_places_share():
    class Sector:
        def add_phenomenon(self, phenomenon, kind, position=None):
            return phenomenon

    args = SimpleNamespace()
    row = {"kind": "supernova-remnant", "subtype": None, "mass_solar": None}
    for share, expected in ((1.0, "Type Ia"), (0.0, "core-collapse")):
        with draw.bound(f"scattered-{share}"):
            remnants = [run_sector._build_scattered(Sector(), args, row, (0.0, 0.0, 0.0), 1000.0, None, share)
                        for _ in range(20)]
        assert {remnant.progenitor_type for remnant in remnants} == {expected}
