"""
The planets of a white dwarf and the disk of a planetary nebula's close pair (GEN.318,
docs/design/systems/exotic-environments-planets-and-compact-binaries.md 1.2 and 1.3): the engulfment radius by the progenitor's
mass, the widening of the survivors' orbits by the mass the star lost, and no planets inside a post-AGB pair's disk.
"""

import math

import pytest

from planetgen import tuning
from planetgen.generation import remnant_orbits
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, mass_loss, planets as planet_physics, remnant_mass
from planetgen.util import draw

WHITE_DWARF_TYPE = "O7VII"


def _system(seed, *, binary, planets=True, **config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = WHITE_DWARF_TYPE
    cfg.BINARY_SYSTEM = binary
    cfg.WIDE_BINARY = False
    if planets:
        cfg.NUM_ORBITS = 5
    for name, value in config.items():
        setattr(cfg, name, value)
    with draw.bound(draw.Stream(f"white-dwarf-planets:{seed}")):
        return StarSystem(system_config=cfg)


def test_the_engulfment_radius_runs_through_the_jovian_points_of_mustill_and_villaver_2012():
    for progenitor, au in tuning.WD_ENGULFMENT_POINTS_AU:
        assert mass_loss.engulfment_radius_au(progenitor) == pytest.approx(au)
    assert mass_loss.engulfment_radius_au(1.0) == pytest.approx(2.7)
    assert mass_loss.engulfment_radius_au(0.8) == pytest.approx(2.7)
    assert mass_loss.engulfment_radius_au(1.75) == pytest.approx(3.25)
    assert mass_loss.engulfment_radius_au(3.75) == pytest.approx(6.9)
    assert mass_loss.engulfment_radius_au(8.0) == pytest.approx(10.0 + 6.2 / 2.5 * 3.0)


@pytest.mark.parametrize("initial, final, expected", [(1.0, 0.54, 1.85), (2.0, 0.65, 3.1), (5.0, 0.9, 5.6),
                                                      (8.0, 1.3, 6.2)])
def test_an_orbit_grows_by_the_designs_mass_ratio(initial, final, expected):
    sun = constants.SOLAR_MASS_TO_KG
    assert mass_loss.expansion_factor(initial * sun, final * sun) == pytest.approx(expected, rel=0.02)


def test_a_planet_in_the_mass_loss_keeps_its_share_of_the_total():
    sun = constants.SOLAR_MASS_TO_KG
    jupiter = 1.898e27
    factor = mass_loss.expansion_factor(2.0 * sun, 0.65 * sun, jupiter)
    # a (M_star + M_p) is the same before and after
    assert factor * (0.65 * sun + jupiter) == pytest.approx(2.0 * sun + jupiter)
    assert factor < mass_loss.expansion_factor(2.0 * sun, 0.65 * sun)


def test_a_progenitor_is_the_models_or_the_remnant_relation_inverted_inside_its_range():
    assert mass_loss.progenitor_mass_sol(0.6, initial_mass_sol=2.2) == 2.2
    assert mass_loss.progenitor_mass_sol(remnant_mass.white_dwarf_mass_sol(3.0)) == pytest.approx(3.0)
    low, high = tuning.WD_PROGENITOR_MASS_RANGE_SOL
    assert mass_loss.progenitor_mass_sol(0.45) == low
    assert mass_loss.progenitor_mass_sol(1.4) == high


def test_a_white_dwarfs_survivors_lie_past_its_progenitors_engulfment_radius_widened_by_the_mass_it_lost():
    found = 0
    for seed in range(12):
        system = _system(seed, binary=False)
        star = system.star
        progenitor = mass_loss.progenitor_mass_sol(star.mass / constants.SOLAR_MASS_TO_KG, star.initial_mass_sol)
        factor = mass_loss.expansion_factor(progenitor * constants.SOLAR_MASS_TO_KG, star.mass)
        assert factor > 1.0
        for body in system.planets:
            inner_edge = body.lower_limit if body.body_type == "a" else body.distance
            assert inner_edge >= mass_loss.engulfment_radius_au(progenitor) * factor * 0.99
            if body.body_type != "a":
                assert body.period == pytest.approx(planet_physics.calculate_orbital_period_years(body.distance, star.mass),
                                                    rel=1e-6)
                found += 1
    assert found


def test_a_star_that_lost_nothing_keeps_its_orbits():
    for seed in range(6):
        with draw.bound(draw.Stream(f"main-sequence:{seed}")):
            cfg = SystemConfig()
            cfg.STAR_TYPE = "G2V"
            cfg.NUM_ORBITS = 4
            system = StarSystem(system_config=cfg)
        before = [(body.distance, getattr(body, "period", None)) for body in system.planets]
        remnant_orbits.widen_for_mass_loss(system.planets, system.star)
        assert [(body.distance, getattr(body, "period", None)) for body in system.planets] == before


def test_a_close_pair_with_a_hot_white_dwarf_holds_a_disk_and_no_planets():
    pairs = 0
    for seed in range(40):
        system = _system(seed, binary=True)
        if system.binary_type != "close":
            continue
        pairs += 1
        assert remnant_orbits.is_post_agb_pair(system.star)
        assert system.planets == [] and system.secondary_planets == []
        primary = system.primary_star
        low, high = tuning.POST_AGB_DISK_MASS_RANGE_SOLAR
        assert primary.disk_kind == remnant_orbits.POST_AGB_DISK_KIND
        assert low <= primary.disk_mass_solar <= high
        assert 0.0 < primary.disk_inner_au <= primary.disk_outer_au
        assert system.secondary_star.disk_kind is None
    assert pairs >= 5


def test_a_single_white_dwarf_has_no_disk_and_planets_may_stay():
    kept = 0
    for seed in range(10):
        system = _system(seed, binary=False)
        assert not remnant_orbits.is_post_agb_pair(system.star)
        assert system.star.disk_kind is None
        kept += len(system.planets)
    assert kept > 0


def test_a_cool_white_dwarfs_pair_is_not_a_planetary_nebulas():
    with draw.bound(draw.Stream("cool-white-dwarf")):
        cfg = SystemConfig()
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = False
        system = StarSystem(system_config=cfg)
    assert system.binary_type == "close"
    cooled = [m for m in system.star.stars if m.yerkes_class in ("VII", "D")
              and m.temperature < tuning.PLANETARY_NEBULA_CENTRAL_MIN_TEMPERATURE_K]
    if cooled:
        assert not remnant_orbits.is_post_agb_pair(system.star)


def test_the_disk_descriptor_survives_a_save_and_load_of_the_system():
    for seed in range(60):
        system = _system(seed, binary=True)
        if system.binary_type == "close":
            break
    else:
        pytest.fail("no close pair in 60 seeds")
    loaded = StarSystem.from_dict(system.to_dict())
    primary = loaded.primary_star
    assert primary.disk_kind == remnant_orbits.POST_AGB_DISK_KIND
    assert primary.disk_mass_solar == pytest.approx(system.primary_star.disk_mass_solar)
    assert not math.isnan(primary.disk_inner_au)


def test_the_star_page_reads_a_post_agb_disk_as_a_circumbinary_one():
    from types import SimpleNamespace
    from planetgen.describe import build
    star = SimpleNamespace(disk_kind="post-agb", disk_mass_solar=0.03, disk_inner_au=0.4, disk_outer_au=120.0)
    text = str(build._star_disk_text(star))
    assert "Circumbinary post-AGB disk of" in text and "no bodies in it" in text
