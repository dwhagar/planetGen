# tests/test_gen_debris_disk.py

"""
GEN.787: a main-sequence star may hold a cold debris disk, at the share the Herschel DEBRIS survey found for its kind
[R4081], as a ring whose radius follows the blackbody law for dust of a drawn temperature.
"""

import math

import pytest

from planetgen import tuning
from planetgen.describe import sheet
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import constants, debris_disk
from planetgen.util import draw


def _share(mass_sol, trials=4000):
    found = sum(debris_disk.draw_disk(mass_sol, 1.0, draw.Stream(f"t:{mass_sol}:{i}")) is not None for i in range(trials))
    return found / trials


@pytest.mark.parametrize("mass,expected", [(0.3, 0.022), (0.7, 0.130), (0.95, 0.143), (1.2, 0.239), (2.0, 0.24), (8.0, 0.24)])
def test_the_share_of_stars_with_a_disk_is_the_surveys_for_their_kind(mass, expected):
    assert debris_disk.incidence(mass) == expected
    assert _share(mass) == pytest.approx(expected, abs=0.03)


def test_the_disk_is_a_ring_at_the_blackbody_radius_of_its_dust_scaled_out():
    assert debris_disk.blackbody_radius_au(278.3, 1.0) == pytest.approx(1.0)
    assert debris_disk.blackbody_radius_au(139.15, 4.0) == pytest.approx(8.0)       # 4 times further for half the temperature, 2 for 4 L
    disks = [d for d in (debris_disk.draw_disk(1.0, 1.0, draw.Stream(f"ring:{i}")) for i in range(2000)) if d]
    assert len(disks) > 200
    for disk in disks:
        assert disk["disk_kind"] == "debris" and disk["disk_mass_solar"] is None
        assert disk["disk_inner_au"] < disk["disk_outer_au"]
        assert disk["disk_outer_au"] / disk["disk_inner_au"] == pytest.approx(1.5 / 0.7)
    belts = sorted(math.sqrt(d["disk_inner_au"] * d["disk_outer_au"]) / math.sqrt(0.7 * 1.5) for d in disks)
    assert belts[len(belts) // 2] == pytest.approx(1.74 * debris_disk.blackbody_radius_au(58.0, 1.0), rel=0.1)


def test_a_brighter_star_has_its_belt_further_out():
    near = debris_disk.draw_disk(1.0, 1.0, _always(0.3))
    far = debris_disk.draw_disk(2.0, 25.0, _always(0.3))
    assert far["disk_inner_au"] == pytest.approx(5.0 * near["disk_inner_au"])


class _always:
    """A stream whose first draw is under every incidence and whose Gaussian is the same each time."""
    def __init__(self, value):
        self.value = value

    def random(self):
        return 0.0

    def gauss(self, mu, sigma):
        return 0.0


def test_the_same_star_always_gets_the_same_disk():
    def build():
        return debris_disk.draw_disk(1.0, 1.0, draw.Stream("debris:1.0:4.6"))
    assert build() == build()


def _star(yerkes, mass_sol=1.0, age=4.6):
    params = {"type": "G2V", "yerkes_class": yerkes, "mass_kg": mass_sol * constants.SOLAR_MASS_TO_KG, "radius_km": 6.957e5,
              "temperature_k": 5772.0, "luminosity_w": constants.SOLAR_LUMINOSITY, "age_gy": age, "lifespan_gy": 10.0,
              "initial_mass_sol": mass_sol, "phase_end_age_gy": 10.0, "metallicity_fe_h": 0.0}
    return Star.from_params(params, SystemConfig(), name="Test")


def test_a_main_sequence_star_may_hold_a_debris_disk_and_a_giant_never_does():
    dwarfs = [_star("V", age=1.0 + i * 0.01) for i in range(300)]
    held = [s for s in dwarfs if s.disk_kind]
    assert 0.05 < len(held) / len(dwarfs) < 0.30
    assert all(s.disk_kind == "debris" and s.disk_mass_solar is None and s.disk_inner_au > 0 for s in held)
    assert not any(_star("III", age=1.0 + i * 0.01).disk_kind for i in range(100))


def test_the_sheet_rule_says_which_stars_have_disk_rows():
    star = {"disk_kind": "debris"}
    assert sheet.applies("star", "G2V", "disk_kind", star) and sheet.applies("star", "G2V", "disk_outer_au", star)
    assert not sheet.applies("star", "G2V", "disk_mass_solar", star)
    assert sheet.applies("star", "M5V", "disk_mass_solar", {"disk_kind": "class-ii"})
    assert not sheet.applies("star", "G2V", "disk_kind", {"disk_kind": None})
    assert ("star", "disk_kind") not in sheet.OPEN_GAPS and ("star", "disk_mass_solar") not in sheet.OPEN_GAPS


def test_a_star_with_a_debris_disk_survives_its_dictionary():
    star = next(s for s in (_star("V", age=1.0 + i * 0.01) for i in range(300)) if s.disk_kind)
    rebuilt = Star.from_dict(star.to_dict(), SystemConfig())
    assert (rebuilt.disk_kind, rebuilt.disk_mass_solar, rebuilt.disk_inner_au, rebuilt.disk_outer_au) == (
        "debris", None, star.disk_inner_au, star.disk_outer_au)
