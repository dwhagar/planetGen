# tests/test_swept_hand_on.py

"""
GEN.1308: the mass an outer disk lost is handed on to where the scatter step says it went (`physics/swept_hand_on.py`, `generation/belt_birth.hand_on`, `generation/disk_ejection.py`,
`irregular_moons.take_from_disks`), and the mass ledger adds up exactly.
"""

import math
import types

import pytest

from planetgen.generation import belt as belt_module, belt_birth, disk_ejection, irregular_moons
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import scatter, swept_hand_on, unit_registry as units
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy

EARTH_KG = units.convert(1.0, "Mearth", "kg")

def test_the_gone_mass_is_shared_by_the_ejection_rule_and_the_swept_split_and_adds_up():
    parts = swept_hand_on.hand_on(1.0e24, 10.0)
    assert sum(parts.values()) == pytest.approx(1.0e24, rel=1e-12)
    ejected = 1.0 - math.exp(-10.0 / scatter.THETA_EJECT)
    split = dict(scatter.SWEPT_SPLIT)
    assert parts["ejected"] == pytest.approx(1.0e24 * ejected)
    assert parts["implanted"] == pytest.approx(1.0e24 * (1.0 - ejected) * split["implanted_inner_belt"]) and parts["captured"] == pytest.approx(1.0e24 * (1.0 - ejected) * split["captured_moons"])
    assert parts["other"] == pytest.approx(1.0e24 * (1.0 - ejected) * (1.0 - split["implanted_inner_belt"] - split["captured_moons"]))
    assert swept_hand_on.hand_on(1.0e24, 100.0)["ejected"] > parts["ejected"]            # a stronger giant throws out more


def test_a_sweep_of_nothing_or_a_giant_with_no_safronov_number_hands_nothing_on():
    assert swept_hand_on.hand_on(0.0, 10.0) is None and swept_hand_on.hand_on(1.0e24, None) is None and swept_hand_on.hand_on(1.0e24, 0.0) is None


def test_the_rogue_bodies_of_an_ejected_mass_follow_the_size_law():
    stream = draw.Stream("test:disk-ejection")
    assert disk_ejection.rogue_bodies(0.0, stream) == (0, 0.0)
    count, total = disk_ejection.rogue_bodies(5.0 * EARTH_KG, stream)
    assert count >= 0 and total <= 5.0 * EARTH_KG and (count == 0) == (total == 0.0)
    small, _ = disk_ejection.rogue_bodies(1.0e20, draw.Stream("test:disk-ejection:small"))
    assert small == 0                              # a disk of small bodies throws out no planet


def test_a_disk_loss_is_a_thrown_planets_record_with_the_giant_as_its_scatterer(monkeypatch):
    monkeypatch.setattr(disk_ejection, "rogue_bodies", lambda ejected_kg, stream: (2, 0.4 * EARTH_KG))
    star = types.SimpleNamespace(mass=1.989e30, snow_line_au=2.7)
    belt = types.SimpleNamespace(body_type='a', swept_kg={"ejected": 3.0 * EARTH_KG}, scatter_giant_kg=1.9e27, scatter_giant_radius_km=7.0e4, scatter_giant_au=5.2, distance=40.0, lower_limit=30.0)
    [loss] = disk_ejection.disk_losses([belt], star)
    assert loss.kind == "disk" and loss.count == 2 and loss.mass_kg == pytest.approx(0.4 * EARTH_KG) and loss.ejected == 1.0
    assert loss.perturber_au == 5.2 and loss.perturber_mass_kg == 1.9e27 and loss.ice_rich and loss.time_yr == disk_ejection.DISK_EJECTION_TIME_YR and loss.theta > 0.0 and loss.hill_radius_au > 0.0
    assert disk_ejection.disk_losses([types.SimpleNamespace(body_type='a', swept_kg=None)], star) == []


def _system(seed):
    with deterministic_entropy(1308):
        draw.set_run_seed(seed)
        cfg = SystemConfig()
        cfg.STAR_TYPE, cfg.PLANETS, cfg.MOONS, cfg.BINARY_SYSTEM = "G2V", True, False, False
        return StarSystem(system_config=cfg)


def test_the_mass_an_outer_disk_lost_is_all_accounted_for_in_a_built_system():
    system = _system(51)             # a system whose outer disk a giant swept (an instability scattered 98 Earth masses of it)
    belts = [body for body in system.planets if body.body_type == 'a']
    [disk] = [belt for belt in belts if belt.swept_kg]
    gone = (disk.cleared_mass_earth + disk.scattered_mass_earth) * EARTH_KG
    assert sum(disk.swept_kg.values()) == pytest.approx(gone, rel=1e-9)
    assert disk.swept_kg["ejected"] > disk.swept_kg["implanted"] > 0.0 and disk.swept_kg["captured"] > 0.0
    assert disk.mass_kg / EARTH_KG <= disk.primordial_mass_earth - disk.cleared_mass_earth - disk.scattered_mass_earth + 1e-12       # the disk's own mass is what is left, not changed again
    received = sum(belt.implanted_kg for belt in belts if belt.origin != "outer-disk")
    assert received == pytest.approx(disk.swept_kg["implanted"], rel=1e-9) or received == 0.0       # the implant goes whole to the inner belt nearest the disk, or nowhere where the system has none


def test_an_inner_belt_that_receives_an_implant_grows_by_it_and_says_so():
    for seed in range(1, 60):                  # the first system with exactly one belt that is not the outer disk (the draws behind a seed move with every change to the generator)
        system = _system(seed)
        inners = [body for body in system.planets if body.body_type == 'a' and body.origin != "outer-disk"]
        if len(inners) == 1:
            break
    [inner] = inners
    star = belt_module._STARS[inner]
    before = inner.primordial_mass_earth
    plain = belt_module.present_state(star, inner.distance, inner.lower_limit, inner.upper_limit, before - inner.cleared_mass_earth - inner.scattered_mass_earth)[1]
    assert inner.receive(1.0e23, "A giant implanted it.")
    assert inner.primordial_mass_earth == pytest.approx(before + 1.0e23 / EARTH_KG) and inner.implanted_kg == 1.0e23
    assert inner.mass_kg >= plain and inner.origin_sentence.endswith("A giant implanted it.")
    assert inner.mass_kg / EARTH_KG <= inner.primordial_mass_earth - inner.cleared_mass_earth - inner.scattered_mass_earth + 1e-12
    assert not inner.receive(0.0)


def test_hand_on_gives_the_implant_to_the_inner_belt_nearest_the_disks_inner_edge():
    received = []
    near = types.SimpleNamespace(body_type='a', origin="gap", upper_limit=3.0, receive=lambda mass, sentence: received.append(("near", mass)) or True)
    far = types.SimpleNamespace(body_type='a', origin="tail", upper_limit=1.0, receive=lambda mass, sentence: received.append(("far", mass)) or True)
    disk = types.SimpleNamespace(body_type='a', origin="outer-disk", lower_limit=20.0, scatter_giant_kg=1.9e27, scatter_giant_radius_km=7.0e4, scatter_giant_au=5.2, cleared_mass_earth=0.0, scattered_mass_earth=10.0,
                                 swept_kg=None)
    belt_birth.hand_on([far, near, disk], types.SimpleNamespace(mass=1.989e30))
    assert received == [("near", pytest.approx(disk.swept_kg["implanted"]))] and disk.swept_kg["implanted"] > 0.0
    belt_birth.hand_on([far, near, disk], types.SimpleNamespace(mass=1.989e30))                      # once only
    assert len(received) == 1


def _disk(captured_kg):
    return types.SimpleNamespace(primordial_mass_earth=10.0, mass_kg=1.0e26, members=[], captured_mass_kg=0.0, swept_kg={"captured": captured_kg})


def test_the_planets_catch_from_the_share_the_scatter_step_handed_on_before_the_disk():
    disk = _disk(4.0e21)
    assert irregular_moons.take_from_disks([disk], [irregular_moons.IrregularSet((), 1.0e21, 0, 0, False)]) == 0.0      # all of it from the share: the disk's own mass does not fall
    assert disk.mass_kg == 1.0e26 and disk.captured_mass_kg == pytest.approx(1.0e21)
    disk = _disk(4.0e21)
    taken = irregular_moons.take_from_disks([disk], [irregular_moons.IrregularSet((), 1.0e22, 0, 0, False)])            # more than the share: the rest comes out of the disk
    assert taken == pytest.approx(6.0e21) and disk.mass_kg == pytest.approx(1.0e26 - 6.0e21) and disk.captured_mass_kg == pytest.approx(1.0e22)
