# src/tests/test_stellar_table.py

"""
The tabulated stars on top of the analytic potential (GEN.1285, TEST.627; galaxy/stellar_table.py).

The table is held against closed forms (a Hernquist sphere and a Miyamoto-Nagai disk given to it as the "model", so the whole potential is the table's), against its
own potential's gradient (the pull and the potential are one field), and, for the real density model, against an independent ring sum by complete elliptic
integrals (the method of reports/map383-bench/gen1285b.py) for the circular speed of the stars.
"""

import dataclasses
import math

import numpy as np
import pytest
from scipy.special import ellipk

from planetgen.galaxy import density, potential as P, stellar_table as T
from tests.audit import representation

G = P.G_PC_KMS2_PER_MSUN


def _zero(x, y, z):
    return 0.0 * (x + y + z)


def _hernquist(mass, scale):
    def density_fn(x, y, z):
        r = np.sqrt(x * x + y * y + z * z)
        return mass * scale / (2.0 * math.pi * r * (r + scale) ** 3)
    return density_fn


def _miyamoto_nagai(mass, a, b):
    def density_fn(x, y, z):
        cyl2 = x * x + y * y
        s = np.sqrt(z * z + b * b)
        return (b * b * mass / (4.0 * math.pi)) * (a * cyl2 + (a + 3.0 * s) * (a + s) ** 2) / ((cyl2 + (a + s) ** 2) ** 2.5 * s ** 3)
    return density_fn


@pytest.fixture(scope="module")
def sphere_table():
    return T.build(_hernquist(1.0e10, 300.0), _zero, G)


@pytest.fixture(scope="module")
def disk_table():
    return T.build(_miyamoto_nagai(1.0e10, 3000.0, 500.0), _zero, G)


@pytest.fixture(scope="module")
def galaxy():
    return P.from_shape(representation.SHAPE, None)


def test_a_residual_of_nothing_is_no_field():
    table = T.build(_hernquist(1.0e10, 300.0), _hernquist(1.0e10, 300.0), G)
    assert table.potential(100.0, 20.0, -30.0) == pytest.approx(0.0, abs=1.0e-9)
    assert max(abs(v) for v in table.pull(100.0, 20.0, -30.0)) < 1.0e-9
    assert table.density(100.0, 20.0, -30.0) == pytest.approx(0.0, abs=1.0e-12)


@pytest.mark.parametrize("radius", (3.0, 30.0, 300.0, 3000.0, 3.0e4, 2.0e5))
def test_the_table_holds_a_hernquist_sphere_in_potential_and_pull(sphere_table, radius):
    x, y, z = 0.6 * radius, 0.0, 0.8 * radius
    assert sphere_table.potential(x, y, z) / (-G * 1.0e10 / (radius + 300.0)) == pytest.approx(1.0, abs=0.002)
    pull = math.sqrt(sum(v * v for v in sphere_table.pull(x, y, z)))
    assert pull / (G * 1.0e10 / (radius + 300.0) ** 2) == pytest.approx(1.0, abs=0.003)


@pytest.mark.parametrize("radius, height", ((500.0, 0.0), (3000.0, 0.0), (3000.0, 400.0), (8000.0, 0.0), (8000.0, 1500.0), (2.0e4, 0.0), (2.0e4, 3000.0), (100.0, 50.0)))
def test_the_table_holds_a_miyamoto_nagai_disk_in_potential(disk_table, radius, height):
    exact = -G * 1.0e10 / math.sqrt(radius ** 2 + (3000.0 + math.sqrt(height ** 2 + 500.0 ** 2)) ** 2)
    assert disk_table.potential(radius, 0.0, height) / exact == pytest.approx(1.0, abs=0.005)


def test_the_pull_is_the_gradient_of_the_potential(disk_table):
    for radius, height, phase in ((300.0, 50.0, 0.3), (1500.0, 400.0, 2.0), (8128.0, 1500.0, 0.3), (15000.0, 3000.0, 2.0), (2000.0, -1500.0, 0.3), (30.0, -20.0, 2.0)):
        x, y, z = radius * math.cos(phase), radius * math.sin(phase), height
        step = 0.003 * math.hypot(radius, height)
        pull = np.array(disk_table.pull(x, y, z))
        gradient = np.array([(disk_table.potential(*(np.array([x, y, z]) + step * axis)) - disk_table.potential(*(np.array([x, y, z]) - step * axis))) / (2.0 * step)
                             for axis in np.eye(3)])
        assert np.linalg.norm(pull + gradient) / np.linalg.norm(pull) < 0.03


def test_the_pull_mirrors_across_the_plane_and_vanishes_at_the_center(disk_table):
    above = disk_table.pull(900.0, 300.0, 250.0)
    below = disk_table.pull(900.0, 300.0, -250.0)
    assert (above[0], above[1], above[2]) == pytest.approx((below[0], below[1], -below[2]), rel=1.0e-9)
    assert disk_table.pull(0.0, 0.0, 0.0) == (0.0, 0.0, 0.0)
    assert disk_table.potential(0.0, 0.0, 0.0) == pytest.approx(disk_table.potential(0.0, 0.0, 1.0e-6), rel=1.0e-6)


def test_beyond_the_grid_the_table_falls_off_without_a_jump(sphere_table):
    inside, outside = sphere_table.potential(T.OUTER_PC * 0.999, 0.0, 0.0), sphere_table.potential(T.OUTER_PC * 1.001, 0.0, 0.0)
    assert abs(outside - inside) < 0.01 * abs(inside) + 1.0e-12
    assert abs(sphere_table.pull(T.INNER_PC * 0.5, 0.0, 0.0)[0]) < abs(sphere_table.pull(T.INNER_PC * 2.0, 0.0, 0.0)[0])


def test_the_same_model_gives_the_same_table():
    one = T.build(_hernquist(1.0e10, 300.0), _zero, G)
    two = T.build(_hernquist(1.0e10, 300.0), _zero, G)
    assert one.digest == two.digest and repr(one) == f"StellarTable({one.digest})"
    assert T.build(_hernquist(2.0e10, 300.0), _zero, G).digest != one.digest


def test_without_a_correction_the_potential_is_the_analytic_one(galaxy):
    bare = dataclasses.replace(galaxy, correction=None)
    assert galaxy.correction is not None
    assert P.potential(bare, 100.0, 50.0, 20.0) == pytest.approx(P.potential(galaxy, 100.0, 50.0, 20.0) - galaxy.correction.potential(100.0, 50.0, 20.0), rel=1.0e-12)
    ax, ay, az = P.acceleration(bare, 100.0, 50.0, 20.0)
    cx, cy, cz = galaxy.correction.pull(100.0, 50.0, 20.0)
    assert P.acceleration(galaxy, 100.0, 50.0, 20.0) == pytest.approx((ax + cx, ay + cy, az + cz), rel=1.0e-12)


def test_the_stellar_density_follows_the_density_model(galaxy):
    shape = representation.SHAPE
    scale = 1.0 / (shape.k_norm * density.model_terms(shape)["mass_per_solar"])
    angles = np.linspace(0.0, 2.0 * math.pi, 32, endpoint=False)
    for radius, height in ((10.0, 0.0), (300.0, 0.0), (1500.0, 0.0), (8128.0, 0.0), (8128.0, 300.0), (8128.0, 2000.0), (16000.0, 0.0), (26000.0, 0.0), (1000.0, 300.0)):
        model = float(density.total_density_array(radius * np.cos(angles), radius * np.sin(angles), np.full_like(angles, height), shape).mean()) * scale
        assert float(P.stellar_density(galaxy, radius, 0.0, height)) / model == pytest.approx(1.0, abs=0.04)


def _ring_speed(shape, radii):
    """The circular speed of the density model's stars alone, km/s, by summing rings of the azimuth-averaged density (complete elliptic integrals)."""
    terms = density.model_terms(shape)
    edge_r = np.geomspace(0.5, 60000.0, 121)
    edge_z = np.concatenate(([0.0], np.geomspace(0.5, 60000.0, 121)))
    ring = np.sqrt((edge_r[:-1] ** 2 + edge_r[1:] ** 2) / 2.0)
    height = (edge_z[:-1] + edge_z[1:]) / 2.0
    volume = (math.pi * (edge_r[1:] ** 2 - edge_r[:-1] ** 2))[:, None] * (2.0 * np.diff(edge_z))[None, :]
    angle = (np.arange(24) + 0.5) * (2.0 * math.pi / 24)
    rho = density.total_density_array(ring[:, None, None] * np.cos(angle), ring[:, None, None] * np.sin(angle), height[None, :, None], shape).mean(axis=2)
    mass = rho / shape.k_norm / terms["mass_per_solar"] * volume
    rg, zg = np.meshgrid(ring, height, indexing="ij")

    def phi(radius):
        k2 = np.minimum(4.0 * radius * rg / ((radius + rg) ** 2 + zg ** 2), 1.0 - 1.0e-12)
        return -G * np.sum(mass * 2.0 * ellipk(k2) / (math.pi * np.sqrt((radius + rg) ** 2 + zg ** 2)))

    return [math.sqrt(r * (phi(r * 1.01) - phi(r * 0.99)) / (0.02 * r)) for r in radii]


def test_the_stars_of_the_real_model_turn_at_the_speed_an_independent_ring_sum_gives(galaxy):
    radii = (1000.0, 2000.0, 4000.0, 8128.0, 12000.0, 20000.0)
    stars_only = dataclasses.replace(galaxy, halo=P.Halo(0.0, 1.0, 1.0), black_hole_solar=0.0)
    for radius, reference in zip(radii, _ring_speed(representation.SHAPE, radii)):
        assert float(P.circular_speed_kms(stars_only, radius)) / reference == pytest.approx(1.0, abs=0.04)


def test_the_analytic_stars_alone_miss_the_model_by_the_amount_the_item_records(galaxy):
    """The reason for GEN.1285: the Miyamoto-Nagai and Hernquist stand-ins undershoot the inner and solar-circle rotation of the density model."""
    bare = dataclasses.replace(galaxy, correction=None, halo=P.Halo(0.0, 1.0, 1.0), black_hole_solar=0.0)
    stars_only = dataclasses.replace(galaxy, halo=P.Halo(0.0, 1.0, 1.0), black_hole_solar=0.0)
    assert float(P.circular_speed_kms(bare, 8128.0)) < 0.9 * float(P.circular_speed_kms(stars_only, 8128.0))
    assert float(P.circular_speed_kms(bare, 300.0)) > 1.3 * float(P.circular_speed_kms(stars_only, 300.0))


def test_a_query_of_many_points_is_one_vector_call(galaxy):
    points = np.random.default_rng(1).normal(size=(3, 50000)) * 5000.0
    ax, ay, az = P.acceleration(galaxy, *points)
    assert ax.shape == ay.shape == az.shape == (50000,) and np.all(np.isfinite(ax)) and np.all(np.isfinite(az))


# --- TEST.627: where the table is held to the ring sum, the speed at the Sun, the cache, and the stand-in registry ---------------------------------------------------------

TOLERANCE = 0.04
"""float: The agreement with the independent ring sum the table states: four percent of the circular speed of the stars."""


def test_the_table_agrees_with_the_ring_sum_at_the_bulge_the_solar_circle_the_disk_edge_and_beyond_16_kpc(galaxy):
    """The four places the item names: the bulge (500 pc and 1.5 kpc), the solar circle, the edge of the disk (three scale lengths of the thin disk out and the
    thick disk's) and past 16 kpc (22 and 30 kpc), each inside `TOLERANCE`."""
    shape = representation.SHAPE
    thin = next(d for d in galaxy.disks if d.name == "thin_disk") if any(d.name == "thin_disk" for d in galaxy.disks) else galaxy.disks[0]
    places = {
        "bulge": (500.0, 1500.0), "solar circle": (8128.0,),
        "disk edge": (3.0 * thin.scale_length_pc, 4.5 * thin.scale_length_pc), "beyond 16 kpc": (22000.0, 30000.0),
    }
    radii = tuple(r for group in places.values() for r in group)
    stars_only = dataclasses.replace(galaxy, halo=P.Halo(0.0, 1.0, 1.0), black_hole_solar=0.0)
    for radius, reference in zip(radii, _ring_speed(shape, radii)):
        got = float(P.circular_speed_kms(stars_only, radius))
        assert got / reference == pytest.approx(1.0, abs=TOLERANCE), (radius, got, reference)


def test_the_whole_galaxy_turns_at_the_sun_at_the_speed_the_model_gives(galaxy):
    """The Sun's circular speed of the test galaxy (`representation.SHAPE`) with the stars, the dark halo, and the central hole is 215.8 km/s; the analytic stars
    alone gave 198.7. The stars-only part of that is held to the ring sum above; the halo and the hole are closed forms, so the whole is the sum in quadrature."""
    whole = float(P.circular_speed_kms(galaxy, 8128.0))
    assert whole == pytest.approx(215.8, abs=3.0)
    assert float(P.circular_speed_kms(dataclasses.replace(galaxy, correction=None), 8128.0)) == pytest.approx(198.7, abs=3.0)
    stars_only = dataclasses.replace(galaxy, halo=P.Halo(0.0, 1.0, 1.0), black_hole_solar=0.0)
    dark_only = dataclasses.replace(galaxy, spheres=(), disks=(), correction=None, black_hole_solar=0.0)
    assert math.hypot(float(P.circular_speed_kms(stars_only, 8128.0)), float(P.circular_speed_kms(dark_only, 8128.0))) == pytest.approx(whole, rel=0.01)


def test_the_table_is_built_once_per_galaxy_and_again_only_when_the_galaxy_changes(monkeypatch):
    """The cache (`for_galaxy`) keeps one `Potential` per shape and halo seed: asking again builds nothing, and a different galaxy builds a new table."""
    built = []
    real = T.build

    def counting(*args, **kwargs):
        built.append(1)
        return real(*args, **kwargs)

    monkeypatch.setattr(T, "build", counting)
    P._for_galaxy.cache_clear()
    try:
        shape = representation.SHAPE
        first = P.for_galaxy(shape, None)
        assert len(built) == 1
        assert P.for_galaxy(shape, None) is first
        assert P.for_galaxy(shape._replace(components=representation.WITH_COMPONENTS.components), None) is first, "star clusters and satellites do not change the galaxy"
        assert len(built) == 1
        other = P.for_galaxy(shape._replace(k_norm=shape.k_norm * 1.1), None)
        assert other is not first and len(built) == 2, "a different galaxy is a new table"
    finally:
        P._for_galaxy.cache_clear()


def test_the_registry_records_that_the_analytic_disks_and_spheres_are_the_baseline_since_the_table():
    """GP1 and GP2 stay in the stand-in registry (the closed forms are still the baseline the table adds to) and say what the table changed; GP3 (the hole's softening)
    is untouched. The item's wish that they be gone is not what the code is: the closed forms remain."""
    import os

    path = os.path.join(os.path.dirname(__file__), "..", "..", "docs", "design", "research", "gap-filling-plan.md")
    with open(path, encoding="utf-8") as handle:
        rows = {line.split("|")[1].strip(): line for line in handle if line.startswith("| GP")}
    assert {"GP1", "GP2", "GP3"} <= set(rows)
    assert "since GEN.1285" in rows["GP1"] and "stellar_table" in rows["GP1"]
    assert "since GEN.1285" in rows["GP2"] and "stellar_table" in rows["GP2"]
    assert "since GEN.1285" not in rows["GP3"]
