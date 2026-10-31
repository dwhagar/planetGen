# tests/test_galaxy_potential.py

"""
GEN.338: the galactic potential as one pure module. Every mass and scale comes from the generated galaxy, so the tests check the shape and the internal
consistency of the curve and the integrator-facing properties, not the Milky Way's numbers.
"""

import math

import numpy as np
import pytest

from planetgen.galaxy import enclosed_mass, mass_chain, shape_draw
from planetgen.galaxy import potential as P

G = P.G_PC_KMS2_PER_MSUN


@pytest.fixture(scope="module")
def shape():
    return shape_draw.default_shape(mass_chain.FALLBACK_SEED)


@pytest.fixture(scope="module")
def galaxy(shape):
    return P.from_shape(shape, mass_chain.FALLBACK_SEED)


def test_the_constants_are_the_codes_own():
    assert G == pytest.approx(4.30091e-3, rel=1e-4)
    assert P.KMS_TO_PC_PER_MYR == pytest.approx(1.0227, rel=1e-4)


def test_the_rotation_curve_rises_then_stays_flat_or_falls_gently(galaxy):
    radii = np.array([100.0, 300.0, 1000.0, 3000.0, 8000.0, 16000.0, 40000.0])
    speed = P.circular_speed_kms(galaxy, radii)
    assert speed[0] < speed[1] < speed[2]
    assert speed.max() > 0 and speed[3:].min() > 0.7 * speed.max() and speed[3:].max() < 1.15 * speed.max()


def test_the_circular_speed_is_the_enclosed_mass_for_the_spherical_parts():
    pot = P.Potential(spheres=(P.Hernquist("bulge", 1.0e10, 500.0),), halo=P.Halo(5.0e11, 10.0, 200000.0), black_hole_solar=4.0e6)
    for r in (50.0, 500.0, 5000.0, 50000.0):
        sphere = 1.0e10 * r * r / (r + 500.0) ** 2
        scale = pot.nfw_scale_pc
        halo = pot.nfw_characteristic_solar * (math.log1p(r / scale) - (r / scale) / (1.0 + r / scale))
        hole = 4.0e6 * r ** 3 / (r * r + 1.0) ** 1.5
        assert P.circular_speed_kms(pot, r) == pytest.approx(math.sqrt(G * (sphere + halo + hole) / r), rel=1e-10)


def _example():
    return P.Potential(spheres=(P.Hernquist("bulge", 1.0e10, 600.0), P.Hernquist("halo", 1.0e9, 8000.0)),
                       disks=(P.Disk("thin_disk", 5.0e10, 3000.0, 300.0), P.Disk("thick_disk", 8.0e9, 2500.0, 900.0)),
                       halo=P.Halo(8.0e11, 12.0, 220000.0), black_hole_solar=4.0e6)


@pytest.mark.parametrize("point", [(8000.0, 0.0, 0.0), (5000.0, 3000.0, 1100.0), (-700.0, 400.0, -250.0), (20.0, -10.0, 5.0), (30000.0, 20000.0, -15000.0)])
def test_the_pull_is_minus_the_gradient_of_the_potential(point):
    pot = _example()
    step = 1.0e-3 * max(1.0, math.sqrt(sum(c * c for c in point)))
    gradient = []
    for axis in range(3):
        up, down = list(point), list(point)
        up[axis] += step
        down[axis] -= step
        gradient.append((P.potential(pot, *up) - P.potential(pot, *down)) / (2.0 * step))
    pull = P.acceleration(pot, *point)
    scale = math.sqrt(sum(c * c for c in pull))
    for measured, expected in zip(pull, gradient):
        assert measured == pytest.approx(-expected, abs=1e-5 * scale)


def test_a_miyamoto_nagai_disk_has_the_density_the_poisson_equation_gives():
    mass, a, b = 5.0e10, 3000.0, 300.0
    pot = P.Potential(disks=(P.Disk("thin_disk", mass, a, b),))
    for radius, height in ((8000.0, 0.0), (8000.0, 400.0), (3000.0, 1100.0)):
        step = 2.0

        def second(f, axis):
            point = [radius, 0.0, height]
            up, down = list(point), list(point)
            up[axis] += step
            down[axis] -= step
            return (f(*up) - 2.0 * f(*point) + f(*down)) / step ** 2

        potential_at = lambda x, y, z: P.potential(pot, x, y, z)
        # Phi depends on R = sqrt(x^2 + y^2): the Laplacian at (R, 0, z) is d2/dx2 + d2/dy2 + d2/dz2 with y evaluated through R.
        laplacian = second(potential_at, 0) + (P.potential(pot, radius, step, height) + P.potential(pot, radius, -step, height) - 2 * P.potential(pot, radius, 0, height)) / step ** 2 + second(potential_at, 2)
        root = math.sqrt(height ** 2 + b ** 2)
        exact = (b * b * mass / (4.0 * math.pi)) * (a * radius ** 2 + (a + 3.0 * root) * (a + root) ** 2) / ((radius ** 2 + (a + root) ** 2) ** 2.5 * root ** 3)
        assert laplacian / (4.0 * math.pi * G) == pytest.approx(exact, rel=2e-3)


def test_the_pull_is_finite_and_zero_at_the_centre(galaxy):
    assert P.acceleration(galaxy, 0.0, 0.0, 0.0) == (0.0, 0.0, 0.0)
    assert math.isfinite(P.potential(galaxy, 0.0, 0.0, 0.0))
    near = P.acceleration(galaxy, 1.0e-9, 0.0, 0.0)
    assert all(math.isfinite(c) for c in near)
    assert P.potential(P.Potential(halo=P.Halo(1.0e12, 10.0, 200000.0)), 0.0, 0.0, 0.0) == pytest.approx(-G * 1.0e12 / (200000.0 / 10.0)
                                                                                                        / (math.log1p(10.0) - 10.0 / 11.0))


def test_a_galaxy_without_a_component_is_legal_and_adds_nothing():
    only_disk = P.Potential(disks=(P.Disk("thin_disk", 5.0e10, 3000.0, 300.0), P.Disk("thick_disk", 0.0, 2500.0, 900.0)))
    assert P.acceleration(only_disk, 8000.0, 0.0, 0.0) == P.acceleration(P.Potential(disks=(only_disk.disks[0],)), 8000.0, 0.0, 0.0)
    nothing = P.Potential()
    assert P.acceleration(nothing, 100.0, 5.0, 2.0) == (0.0, 0.0, 0.0) and P.potential(nothing, 1.0, 2.0, 3.0) == 0.0
    with pytest.raises(ValueError):
        P.Potential(spheres=(P.Hernquist("bulge", -1.0, 500.0),))
    with pytest.raises(ValueError):
        P.Potential(halo=P.Halo(1.0, 0.0, 1.0))


def test_an_arrays_answer_is_the_floats(galaxy):
    x = np.array([8000.0, 100.0, -3000.0])
    z = np.array([0.0, 50.0, 700.0])
    ax, ay, az = P.acceleration(galaxy, x, 0.5 * x, z)
    for i in range(3):
        single = P.acceleration(galaxy, float(x[i]), 0.5 * float(x[i]), float(z[i]))
        assert (ax[i], ay[i], az[i]) == pytest.approx(single, rel=1e-12)
    assert P.potential(galaxy, 1000.0, 2000.0, 100.0) == P.potential(galaxy, 1000.0, 2000.0, 100.0)


def test_the_vertical_force_pulls_towards_the_plane_and_grows_with_the_disk_mass(galaxy):
    assert P.vertical_force(galaxy, 8000.0, 1100.0) > 0.0 > P.vertical_force(galaxy, 8000.0, -1100.0)
    heavier = P.Potential(galaxy.spheres, tuple(d._replace(mass_solar=2.0 * d.mass_solar) for d in galaxy.disks), galaxy.halo, galaxy.black_hole_solar)
    assert P.vertical_force(heavier, 8000.0, 1100.0) > P.vertical_force(galaxy, 8000.0, 1100.0)


def _leapfrog(pot, position, velocity, step_myr, steps):
    dt = step_myr * P.KMS_TO_PC_PER_MYR            # in pc / (km/s): a Myr is 1.0227 of them
    x, v = np.array(position, float), np.array(velocity, float)
    a = np.array(P.acceleration(pot, *x))
    track = []
    for _ in range(steps):
        v = v + 0.5 * dt * a
        x = x + dt * v
        a = np.array(P.acceleration(pot, *x))
        v = v + 0.5 * dt * a
        track.append((x.copy(), v.copy()))
    return track


def test_a_circular_orbit_keeps_its_radius_and_angular_momentum_over_five_gyr(galaxy):
    radius = 8000.0
    speed = P.circular_speed_kms(galaxy, radius)
    track = _leapfrog(galaxy, (radius, 0.0, 0.0), (0.0, speed, 0.0), 1.0, 5000)
    radii = np.array([math.hypot(x[0], x[1]) for x, _ in track])
    assert np.abs(radii - radius).max() < 1.0
    lz = np.array([x[0] * v[1] - x[1] * v[0] for x, v in track])
    assert np.abs(lz / lz[0] - 1.0).max() < 1e-10
    energy = np.array([0.5 * float(v @ v) + P.potential(galaxy, *x) for x, v in track[::500]])
    assert abs(energy[-1] - energy[0]) / abs(energy[0]) < 1e-3 * 5.0
