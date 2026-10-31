"""
GEN.809 (wiring): a giant's magnetic moment and dipole class are the dynamo of its interior's metallic layer (`magnetism.giant_field`,
`giant_dynamo.giant_field`) in place of the three anchors of Neptune, Saturn and Jupiter, which an ice giant keeps until its layers carry radii.
"""

import functools

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import dynamo, giant_draw, giant_interior, magnetism
from planetgen.util import draw


@functools.lru_cache(maxsize=None)
def _giants(count=30):
    found = {"gas_giant": [], "ice_giant": []}
    for seed in range(count):
        with draw.bound(seed):
            system = StarSystem(SystemConfig())
        for planet in system.planets:
            if planet.body_type == "g" and planet.interior_kind in found:
                found[planet.interior_kind].append(planet)
    return found


def _first(kind):
    giants = _giants()[kind]
    if not giants:
        pytest.skip(f"no {kind} in the sample systems")
    return giants[0]


def test_a_gas_giants_moment_is_its_dynamos():
    planet = _first("gas_giant")
    found = magnetism.giant_field(planet)
    assert found is not None and found.region.layer == "metallic_envelope"
    assert planet.magnetic_moment_a_m2 == pytest.approx(found.field.moment_a_m2, rel=1e-12)
    assert found.power_w == pytest.approx(dynamo.giant_power_w(
        giant_draw.own_flux_w_m2(planet) * 4.0 * 3.141592653589793 * (planet.radius * 1000.0) ** 2), rel=1e-9)


def test_a_slow_spin_makes_a_gas_giants_field_multipolar_and_a_fast_one_a_dipole():
    planet = _first("gas_giant")
    planet.rotation_period_hours = 400.0
    magnetism.generate_field(planet)
    assert planet.dipole_class == "multipolar"
    planet.rotation_period_hours = 8.0
    magnetism.generate_field(planet)
    assert planet.dipole_class in ("earth-like", "strong")


def test_an_ice_giant_with_no_placed_layer_keeps_the_anchor_moment():
    planet = _first("ice_giant")
    assert magnetism.giant_field(planet) is None
    mass_earth = planet.mass / 5.972e24
    expected = magnetism._giant_moment(mass_earth)
    assert planet.magnetic_moment_a_m2 == pytest.approx(expected, rel=2.0)    # log-normal 0.4 dex about the anchor law
    assert planet.dipole_class != "none"


def test_a_giant_without_a_stored_interior_has_no_dynamo_field():
    planet = _first("gas_giant")
    kept = planet.interior_kind
    try:
        planet.interior_kind = None
        assert giant_interior.interior_of(planet) is None and magnetism.giant_field(planet) is None
    finally:
        planet.interior_kind = kept
