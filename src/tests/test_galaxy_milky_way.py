# tests/test_galaxy_milky_way.py

"""
GEN.118 The galaxy bulge was about 40 times too light, so edge-on views
showed no bulge, and GEN.119 the density model had no thick disk
(docs/TODO.md).

`planetgen plan`'s default galaxy shape against the Milky Way's published
structure: Bland-Hawthorn & Gerhard 2016 (ARA&A 54:529, "BHG16") for the
masses, scale lengths and heights; Dwek et al. 1995 (ApJ 445:716) for the
bulge's boxy shape, fitted to the COBE/DIRBE edge-on image; Wegg & Gerhard
2013 (MNRAS 435:1874) for the bar's angle; McKee et al. 2015 (ApJ 814:13)
for the stars' surface density at the Sun. Each tolerance is the published
uncertainty, or the spread between the sources.
"""

import argparse
import math

import pytest

from planetgen import tuning
from planetgen.cli.generate import add_plan_arguments
from planetgen.galaxy import density
from planetgen.physics import mathcheck


def _default_shape():
    args = argparse.ArgumentParser()
    add_plan_arguments(args)
    a = args.parse_args([])
    amplitude, level = density.arm_terms(a.arm_density, a.interarm_density)
    return density.build_galaxy_shape(
        a.disk_scale_length_pc, a.disk_scale_height_pc, a.bulge_scale_radius_pc, a.bulge_amplitude,
        a.arm_count, math.radians(a.pitch_angle_deg), amplitude, arm_level=level, core_amplitude=a.core_density,
    )


SHAPE = _default_shape()
TERMS = density.model_terms(SHAPE)
SOLAR_RADIUS_PC = tuning.GALAXY_SOLAR_RADIUS_TO_SCALE_LENGTH * SHAPE.disk_scale_length_pc


def _at(radius, angle, z=0.0):
    return (radius * math.cos(angle), radius * math.sin(angle), z)


def _ring_mean(radius, z, steps=360):
    """The density averaged around a ring (arms and bar averaged out)."""
    return sum(density.relative_density(_at(radius, 2 * math.pi * i / steps, z), SHAPE)
               for i in range(steps)) / steps


def _local_scale(profile, a, b):
    """The exponential scale of `profile` between `a` and `b`."""
    return (b - a) / math.log(profile(a) / profile(b))


def test_the_math_check_uses_planetgen_plans_defaults():
    assert mathcheck._default_galaxy_shape() == SHAPE


def test_the_bulge_holds_a_milky_way_share_of_the_stars():
    # BHG16: 1.4-1.7e10 of 5e10 M_sun of stars are in the bulge.
    masses = density.component_masses(SHAPE)
    assert 1.4 / 5.0 <= masses["bulge"] / sum(masses.values()) <= 1.7 / 5.0


def test_the_thick_disk_holds_a_milky_way_share_of_the_disk():
    # BHG16: thick disk 6 +- 3 e9 against thin disk 3.5 +- 1 e10 M_sun.
    masses = density.component_masses(SHAPE)
    assert 3.0 / 45.0 <= masses["thick_disk"] / masses["thin_disk"] <= 9.0 / 25.0


def test_the_disks_have_the_milky_ways_scale_lengths_face_on():
    # Face on, the disk's surface density between 5 and 12 kpc falls off
    # with the thin disk's 2.6 +- 0.5 kpc (BHG16), the thick disk's
    # shorter 2.0 kpc pulling it in only a little.
    def surface(radius):
        return sum(_ring_mean(radius, z, steps=36) for z in range(0, 6000, 25))

    assert _local_scale(surface, 5000.0, 12000.0) == pytest.approx(2600.0, abs=500.0)


def test_the_disks_have_the_milky_ways_scale_heights_at_the_sun():
    # The thin disk falls off with 300 +- 50 pc and the thick disk with
    # 900 +- 180 pc (BHG16), and far above the plane the thick disk is
    # all that's left.
    sun = _at(SOLAR_RADIUS_PC, 1.0)

    def part(index):
        return lambda z: density._components((sun[0], sun[1], z), SHAPE)[index]

    assert _local_scale(part(1), 1000.0, 2000.0) == pytest.approx(300.0, abs=50.0)
    assert _local_scale(part(2), 3000.0, 4500.0) == pytest.approx(900.0, abs=180.0)
    assert _local_scale(lambda z: _ring_mean(SOLAR_RADIUS_PC, z, steps=36), 3000.0, 4500.0) == pytest.approx(
        900.0, abs=180.0)


def test_the_thick_disk_is_four_percent_of_the_plane_at_the_sun():
    # BHG16: 0.04 +- 0.02.
    _bulge, thin, thick, _arm = density._components(_at(SOLAR_RADIUS_PC, 1.0), SHAPE)
    assert thick / thin == pytest.approx(0.04, abs=0.02)


def _bulge_half_density_pc(direction):
    """How far the bulge reaches along a unit vector in the bar's own
    frame before its density halves."""
    cos_a, sin_a = TERMS["bar_cos"], TERMS["bar_sin"]

    def bulge(distance):
        u, v, w = (distance * c for c in direction)
        return density._components((u * cos_a - v * sin_a, u * sin_a + v * cos_a, w), SHAPE)[0]

    peak, low, high = bulge(0.0), 0.0, 20000.0
    for _ in range(60):
        mid = (low + high) / 2
        low, high = (mid, high) if bulge(mid) > peak / 2 else (low, mid)
    return low


def test_the_bulge_is_cobes_boxy_bar():
    # Dwek et al. 1995 G2: exp(-r_s^2/2) halves at r_s = sqrt(2 ln 2), so
    # at 1.177 x (1.58, 0.62, 0.43) kpc along, across and above the bar.
    reach = math.sqrt(2 * math.log(2))
    assert _bulge_half_density_pc((1, 0, 0)) == pytest.approx(reach * 1580.0, rel=0.01)
    assert _bulge_half_density_pc((0, 1, 0)) == pytest.approx(reach * 620.0, rel=0.01)
    assert _bulge_half_density_pc((0, 0, 1)) == pytest.approx(reach * 430.0, rel=0.01)
    # Boxy, not an ellipsoid: on the diagonal it reaches farther than an
    # ellipsoid with the same axes would.
    diagonal = (math.sqrt(0.5), 0.0, math.sqrt(0.5))
    ellipsoid = reach / math.sqrt(0.5 / 1580.0 ** 2 + 0.5 / 430.0 ** 2)
    assert _bulge_half_density_pc(diagonal) > 1.02 * ellipsoid


def test_the_bar_leads_the_sun_by_27_degrees():
    # Wegg & Gerhard 2013: 27 +- 2 degrees, its near end at positive
    # galactic longitude: ahead of the Sun in the direction of rotation
    # (counterclockwise from galactic north, increasing theta).
    angles = [2 * math.pi * i / 3600 for i in range(3600)]
    densest = max(angles, key=lambda a: density._components(_at(1500.0, a), SHAPE)[0])
    lead = math.degrees((densest - TERMS["solar_angle_rad"]) % math.pi)
    assert lead == pytest.approx(27.0, abs=2.0)


def test_the_bars_central_density_matches_its_mass_and_the_suns_disk():
    # The bulge's central density over the disk's surface density at the
    # Sun. Published: 1.55e10 M_sun (BHG16) in Dwek's G2 shape, whose
    # volume is x0 y0 z0 times 20.65, gives 1.78 M_sun/pc^3 at the center;
    # the stars' surface density at the Sun is about 33 M_sun/pc^2
    # (McKee et al. 2015). So 0.054 per parsec, with the bulge mass's
    # +-10% and the surface density's +-10%.
    center = density._components((0.0, 0.0, 0.0), SHAPE)[0]
    sun = _at(SOLAR_RADIUS_PC, 0.3)

    def disk(z):
        _bulge, thin, thick, _arm = density._components((sun[0], sun[1], z), SHAPE)
        return thin + thick

    surface = 2 * sum(disk(z + 2.5) for z in range(0, 12000, 5)) * 5
    published = 1.55e10 / (1580.0 * 620.0 * 430.0 * 20.65) / 33.0
    assert center / surface == pytest.approx(published, rel=0.2)


def test_the_center_is_bulge_dominated():
    shares = density.population_densities((0.0, 0.0, 0.0), SHAPE)
    assert shares["bulge"] > 0.7 * sum(shares.values())
    above = density.population_densities((0.0, 0.0, 800.0), SHAPE)
    assert above["bulge"] > 0.5 * sum(above.values())


def _edge_on_column(x, z):
    """Edge-on light at sky position `(x, z)`, seen from the Sun's side:
    the density summed along the Sun-center line."""
    angle = TERMS["solar_angle_rad"]
    toward = (math.cos(angle), math.sin(angle))
    across = (-toward[1], toward[0])
    return sum(density.relative_density(
        (x * across[0] + y * toward[0], x * across[1] + y * toward[1], z), SHAPE)
        for y in range(-20000, 20001, 50))


def test_the_bulge_stands_out_edge_on():
    # 1 kpc above the plane, edge on, the center outshines the disk 4 kpc
    # out. The old 200 pc bulge gave 2.5 (a flat disk); McMillan 2017's
    # Milky Way model gives 3.8.
    assert _edge_on_column(0.0, 1000.0) > 3.8 * _edge_on_column(4000.0, 1000.0)
    # And 600 pc up, most of the center's light is the bulge's.
    angle = TERMS["solar_angle_rad"]
    bulge = sum(density._components(_at(y, angle, 600.0), SHAPE)[0] for y in range(-20000, 20001, 50))
    assert SHAPE.k_norm * bulge > 0.5 * _edge_on_column(0.0, 600.0)
