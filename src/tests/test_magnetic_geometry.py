# tests/test_magnetic_geometry.py

"""
GEN.1350: the direction of a magnetic dipole (planets, moons, rogue planets, neutron stars) and the polar dipole field of a star
(`physics/magnetic_geometry.py`), drawn from streams of the body's own.
"""

import math
import statistics

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import NeutronStar
from planetgen.generation.system import StarSystem
from planetgen.physics import chemically_peculiar, magnetic_geometry as mg

KEYS = [f"body-{n}" for n in range(3000)]


def _tilts(kind):
    return [mg.of_body(kind, key)[0] for key in KEYS]


def test_an_axial_dipole_lies_near_the_spin_axis_and_a_body_of_our_own_system_is_inside_the_draw():
    tilts = _tilts("axial")
    assert all(0.0 <= tilt <= 90.0 for tilt in tilts)
    assert statistics.median(tilts) == pytest.approx(7.0, abs=0.6)
    assert sum(1 for tilt in tilts if tilt >= 11.0) / len(tilts) == pytest.approx(0.18, abs=0.03)       # Earth's 11 degrees is neither rare nor the rule
    assert sum(1 for tilt in tilts if tilt <= 4.0) / len(tilts) > 0.1                                      # Ganymede's 4 degrees and Saturn's 0 are common enough


@pytest.mark.parametrize("kind", ["multipolar", "ice-giant"])
def test_a_multipolar_dipole_and_an_ice_giants_point_any_way(kind):
    tilts = _tilts(kind)
    assert statistics.median(tilts) == pytest.approx(60.0, abs=3.0)     # cos(tilt) uniform, folded to 0 to 90 degrees
    assert max(tilts) > 85.0 and min(tilts) < 5.0
    assert sum(1 for tilt in tilts if 40.0 <= tilt <= 65.0) / len(tilts) > 0.25      # Uranus (59) and Neptune (47) lie where the draw puts a third of them


def test_the_offset_class_follows_the_kind_of_dynamo():
    assert {mg.of_body(kind, "k")[1] for kind in ("ice-giant", "multipolar", "mercury")} == {"offset"}
    assert mg.of_body("axial", "k")[1] == "centred"


def test_a_body_gets_the_same_dipole_every_time_and_no_other_body_shares_it():
    assert mg.of_body("axial", "same") == mg.of_body("axial", "same")
    assert len({mg.of_body("axial", key)[0] for key in KEYS[:50]}) == 50


def _inclinations(pulsar_type, age_yr):
    return [mg.neutron_star_inclination_deg(pulsar_type, age_yr, key) for key in KEYS]


def test_a_neutron_stars_axes_align_with_age_and_a_recycled_pulsar_is_isotropic_again():
    young, middle, old = _inclinations("young", 1.0e4), _inclinations("young", 2.0e7), _inclinations("non-pulsing", 5.0e9)
    assert statistics.median(young) == pytest.approx(60.0, abs=3.0)
    assert statistics.median(middle) < statistics.median(young) - 15.0
    assert max(old) < 0.01
    recycled = _inclinations("millisecond", 5.0e9)
    assert statistics.median(recycled) == pytest.approx(60.0, abs=3.0) and all(0.0 <= value <= 90.0 for value in recycled)


def test_a_generated_neutron_star_stores_its_inclination():
    star = NeutronStar(SystemConfig())
    assert 0.0 <= star.magnetic_inclination_deg <= 90.0
    assert "magnetic_inclination_deg" in star.SERIALIZABLE_FIELDS


def test_the_rossby_number_of_the_sun_is_near_two():
    assert mg.rossby_number(1.0, 25.4 * 24.0) == pytest.approx(1.76, abs=0.05)
    assert mg.convective_turnover_days(0.5) > mg.convective_turnover_days(1.0)


def _field(mass_sol, period_hours, age_gy, temperature_k=5772.0, yerkes="V"):
    mass_kg = mass_sol * 1.989e30
    return mg.star_field_gauss(yerkes, mass_sol, mass_kg, 696000.0 * mass_sol ** 0.8, temperature_k, age_gy, period_hours)


def test_a_cool_dwarfs_field_follows_the_rossby_number_about_the_suns():
    suns = [_field(1.0, 25.4 * 24.0, 4.0 + n * 1.0e-3) for n in range(400)]
    assert statistics.median(suns) == pytest.approx(mg.B_SUN_POLAR_GAUSS, rel=0.25)
    assert statistics.stdev(math.log10(value) for value in suns) == pytest.approx(mg.FIELD_SCATTER_DEX, abs=0.05)
    young = statistics.median(_field(1.0, 4.0 * 24.0, 0.3 + n * 1.0e-4) for n in range(200))
    assert young / statistics.median(suns) == pytest.approx((25.4 / 4.0) ** 1.4, rel=0.35)      # a fast young Sun has a much stronger field
    assert max(_field(0.2, 12.0, 0.1 + n * 1.0e-4, temperature_k=3200.0) for n in range(200)) <= mg.MAX_FIELD_GAUSS


def test_a_fully_convective_m_dwarf_that_spins_fast_has_a_field_of_hundreds_of_gauss():
    field = statistics.median(_field(0.2, 24.0, 1.0 + n * 1.0e-4, temperature_k=3200.0) for n in range(200))
    assert 100.0 < field <= mg.MAX_FIELD_GAUSS


@pytest.mark.parametrize("args", [
    dict(mass_sol=1.0, period_hours=None, age_gy=4.0),                                         # no rotation period
    dict(mass_sol=1.0, period_hours=600.0, age_gy=4.0, yerkes="III"),                         # a giant
    dict(mass_sol=1.0, period_hours=600.0, age_gy=4.0, yerkes="D"),                           # a white dwarf
    dict(mass_sol=0.04, period_hours=600.0, age_gy=4.0, temperature_k=2000.0),                # below the range of the turnover time
])
def test_a_star_the_sources_do_not_cover_has_no_field(args):
    assert _field(**args) is None


def test_an_ap_dwarf_has_the_polar_field_of_its_kind_and_an_am_dwarf_none():
    found = {}
    for n in range(4000):
        mass_sol, temperature_k, age_gy = 2.0 + n * 1.0e-4, 8000.0 + n * 0.5, 0.5 + n * 1.0e-4
        peculiar = chemically_peculiar.of_star("V", mass_sol * 1.989e30, 1.7 * 696000.0, temperature_k, age_gy)
        if peculiar is not None and peculiar.kind not in found:
            found[peculiar.kind] = (peculiar, _field(mass_sol, 100.0, age_gy, temperature_k=temperature_k))
    ap, ap_field = found[chemically_peculiar.AP]
    assert ap_field == ap.field_gauss and 300.0 <= ap_field <= 40000.0
    if chemically_peculiar.AM in found:
        assert found[chemically_peculiar.AM][1] is None


def test_every_generated_body_with_a_dipole_has_a_direction_and_one_without_has_none():
    seen = 0
    for _ in range(6):
        cfg = SystemConfig()
        cfg.MOONS = True
        cfg.MAX_PLANETS = True
        system = StarSystem(system_config=cfg)
        for planet in system.planets:
            if planet.body_type == "a":
                continue
            for body in (planet, *planet.moons):
                if (body.magnetic_moment_a_m2 or 0.0) > 0.0:
                    seen += 1
                    assert 0.0 <= body.dipole_tilt_deg <= 90.0 and body.dipole_offset_class in ("centred", "offset")
                else:
                    assert body.dipole_tilt_deg is None and body.dipole_offset_class is None
        for star in system.stars:
            assert star.surface_dipole_field_gauss is None or 0.0 < star.surface_dipole_field_gauss <= 40000.0
    assert seen
