"""
Magnetars as a scatter class of neutron stars (GEN.278, docs/design/systems/later-phenomena-scatter-rules.md section 3): about a
hundred rows over a Milky Way-sized galaxy, each on the young population and not in the halo, whatever the scatter's
mass cut; the fill never draws one; a row builds into a magnetar of the row's mass.
"""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import density, mass_chain
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.skeleton import expected_system_count_at_density_1
from planetgen.generation import phenomenon_scatter as scatter
from planetgen.generation import run_galaxy, run_sector
from planetgen.generation.phenomena import registry
from planetgen.generation.phenomena.compact_remnant import NeutronStar
from planetgen.generation.system import SystemConfig
from planetgen.galaxy.sector import SpaceSector
from tests.toy_shape import toy_shape

MAGNETAR = ("neutron-star", scatter.MAGNETAR_SUBTYPE)
MILKY_WAY = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
    bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
    masses=(5.0e10, mass_chain.median_halo_mass(5.0e10)),
)
TOY = toy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
EDGE_PC = tuning.DEFAULT_SECTOR_EDGE_PC


def _young_and_total_integrals(shape):
    """`(young, every component)` of the model's stars, by integrating `population_densities_array` over a grid (the
    arms averaged over 24 angles) and `integrated_masses` with the halo: the same units, times `k_norm`."""
    length, height = shape.disk_scale_length_pc, shape.disk_scale_height_pc
    r = (np.arange(240) + 0.5) * (12.0 * length / 240)
    z = (np.arange(160) + 0.5) * (8.0 * height / 160)
    angle = (np.arange(24) + 0.5) * (2.0 * math.pi / 24)
    x = (r[:, None] * np.cos(angle)[None, :])[None, :, :]
    y = (r[:, None] * np.sin(angle)[None, :])[None, :, :]
    young = density.population_densities_array(x, y, z[:, None, None], shape)["young"].mean(axis=2)
    dr, dz = 12.0 * length / 240, 8.0 * height / 160
    young_total = float((young * 2.0 * math.pi * r[None, :] * dr).sum() * 2.0 * dz)
    masses = density.integrated_masses(shape)
    every = (sum(masses.values()) + density.halo_mass(shape)) * shape.k_norm
    return young_total, every


def test_the_magnetars_are_the_last_scatter_class_and_follow_the_young_stars():
    assert scatter.SCATTER_CLASSES[-1] == MAGNETAR
    assert registry.scatter_class_gate(*MAGNETAR) == "young stars"
    assert registry.scatter_class_gate("neutron-star", None) == "stellar density"
    assert registry.scatter_class_gate("black-hole", "stellar") == "stellar density"


def test_the_cut_does_not_touch_the_magnetars():
    for cut in (1.0, 8.0, 14.0, 20.0):
        assert scatter.share_above(*MAGNETAR, cut) == 1.0
        assert scatter.mass_range(*MAGNETAR, cut, True) == tuning.NEUTRON_STAR_MASS_RANGE_SOLAR
        assert scatter.mass_range(*MAGNETAR, cut, False) is None
    # So the fill, which draws what the cut leaves below it, draws none.
    center = (8200.0, 0.0, 0.0)
    for seed in range(40):
        draws, _rng = scatter.below_cut_draws((9, 0, 0), center, TOY, 1.0e6, 14.0, seed)
        assert MAGNETAR not in [(kind, subtype) for kind, subtype, _range, _count in draws]


def test_the_rate_is_a_ten_millionth_of_the_neutron_stars_per_young_star():
    """Of the research neutron star rate (a stand-in: the magnetars follow the young stars, not the metals, GEN.554)."""
    assert scatter.class_rate_per_star(*MAGNETAR) == pytest.approx(
        tuning.phenomenon_rate_per_star("neutron-star") * 1.0e-7 / tuning.MAGNETAR_YOUNG_STAR_SHARE)


def test_the_young_share_is_the_young_populations_share_of_the_stars():
    young, every = _young_and_total_integrals(MILKY_WAY)
    assert tuning.MAGNETAR_YOUNG_STAR_SHARE == pytest.approx(young / every, rel=0.05)


def test_a_milky_way_holds_about_a_hundred_magnetars():
    young, _every = _young_and_total_integrals(MILKY_WAY)
    e_value = expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY)
    expected = scatter.class_rate_per_star(*MAGNETAR) * e_value * young / EDGE_PC ** 3
    assert 30.0 <= expected <= 250.0     # about 100 active, 30 to 600 [R3346]
    # Against the neutron stars the same arithmetic gives the 1e9 of the Galaxy.
    stars = e_value * (young / tuning.MAGNETAR_YOUNG_STAR_SHARE) / EDGE_PC ** 3
    assert 3.0e8 <= scatter.class_rate_per_star("neutron-star", None) * stars <= 3.0e9


def _toy_rows(subtype_wanted, min_mass_solar=14.0, boost=1.0):
    extents = [(shell, min(shell, 3)) for shell in range(1, 9)]
    rows = []
    for shell, high in extents:
        rows += [row for row in scatter.scatter_shell(TOY, shell, high, EDGE_PC, 600.0 * boost, 5,
                                                      min_mass_solar=min_mass_solar) if row[3] == "neutron-star"
                 and row[4] == subtype_wanted]
    return rows


def test_a_scatter_places_magnetars_whatever_the_cut_and_no_other_neutron_stars_at_the_default(monkeypatch):
    monkeypatch.setattr(tuning, "MAGNETAR_PER_NEUTRON_STAR", 1.0e-1)     # enough of them for a toy galaxy
    rows = _toy_rows(scatter.MAGNETAR_SUBTYPE)
    assert len(rows) > 20
    assert _toy_rows(None) == []         # the default cut of 14 solar masses scatters no ordinary neutron star
    column = {name: i for i, name in enumerate(scatter.PHENOMENON_SCATTER_COLUMNS)}
    for row in rows:
        low, high = tuning.NEUTRON_STAR_MASS_RANGE_SOLAR
        assert low <= row[column["mass_solar"]] <= high
        assert row[column["velocity_x_kms"]] is None


def test_magnetars_sit_in_the_thin_young_layer_not_the_thick_one_remnants_have(monkeypatch):
    monkeypatch.setattr(tuning, "MAGNETAR_PER_NEUTRON_STAR", 1.0e-1)
    magnetars = _toy_rows(scatter.MAGNETAR_SUBTYPE, min_mass_solar=1.0)
    ordinary = _toy_rows(None, min_mass_solar=1.0)
    z = scatter.PHENOMENON_SCATTER_COLUMNS.index("position_z_mpc")
    assert len(magnetars) > 20 and len(ordinary) > 20
    magnetar_height = np.mean([abs(row[z]) for row in magnetars]) / 1000.0
    ordinary_height = np.mean([abs(row[z]) for row in ordinary]) / 1000.0
    assert magnetar_height < 0.95 * ordinary_height
    assert max(abs(row[z]) for row in magnetars) / 1000.0 < 20.0 * TOY.disk_scale_height_pc * 50 / 300


def test_the_young_population_has_nothing_in_the_halo_or_far_above_the_plane():
    """The weight a magnetar band is drawn by is the young density: none where only halo and thick-disk stars are."""
    x = np.array([[8200.0, 3000.0], [0.0, 15000.0]])
    y = np.zeros_like(x)
    z = np.array([[20000.0], [60000.0]])
    populations = density.population_densities_array(x, y, z, MILKY_WAY)
    assert (populations["halo"] > 0.0).all()
    assert (populations["young"] < 1.0e-60 * populations["halo"]).all()
    assert (populations["cluster"] == 0.0).all()
    in_plane = density.population_densities_array(np.array([8200.0]), np.array([0.0]), 0.0, MILKY_WAY)
    assert in_plane["young"][0] > 0.0


def test_a_magnetar_row_builds_a_magnetar_of_its_mass(mysql_config):
    args = run_galaxy._default_generation_args(config=mysql_config)
    row = {"kind": "neutron-star", "subtype": scatter.MAGNETAR_SUBTYPE, "seed": 0, "mass_solar": 1.62}
    for seed in range(12):
        row["seed"] = seed
        entry = run_sector._seeded_build(seed, lambda: run_sector._build_scattered(
            SpaceSector(name="Probe"), args, row, (0.0, 0.0, 0.0), 1000.0, 14.0))
        star = entry.phenomenon
        assert star.pulsar_type == "magnetar" and star.mass_solar == pytest.approx(1.62)
        low, high = tuning.MAGNETAR_AGE_RANGE_GY
        assert low <= star.age <= high
        assert star.magnetic_field_gauss >= tuning.MAGNETAR_MAGNETIC_FIELD_GAUSS_RANGE[0]


def test_a_free_neutron_star_is_never_young_enough_to_be_a_magnetar(monkeypatch):
    monkeypatch.setattr(tuning, "MAGNETAR_BIRTH_FRACTION", 1.0)
    assert tuning.COMPACT_REMNANT_AGE_RANGE_GY[0] * 1.0e9 >= tuning.MAGNETAR_ACTIVE_LIFETIME_YR
    config = SystemConfig()
    assert all(NeutronStar(config).pulsar_type != "magnetar" for _ in range(300))
    # A remnant's core carries the remnant's age, and keeps the birth lottery.
    assert NeutronStar(config, age_gy=1.0e-6).pulsar_type == "magnetar"
