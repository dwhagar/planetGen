# tests/test_nebula_contour.py

"""
GEN.1267 and MAP.423: a nebula is the part of its sphere where the galaxy's density, read for it, is at least the form level. The generator's
containment test (`nebula_density.inside`, `containment.contour_test`) and the map's drawing (`nebula_volume.sample_group`) read one density; a
hot-star nebula anchors its own density on its sphere, so it has a shape in the outer disk, where the model's young-star density alone is under the
form level.
"""

import math

import numpy as np
import pytest

from planetgen.galaxy import containment, density, nebula_blocks, nebula_density as nd, nebula_field, nebula_volume
from planetgen.physics.units import ly_to_pc

SEED = bytes(range(16))
FAMILIES = nebula_blocks.FAMILIES


@pytest.fixture(scope="module")
def shape():
    return density.build_galaxy_shape(
        disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
        arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4, nuclear_amplitude=1.0)


def _ball(center, radius, n, seed=5):
    """`n` points in the cube of side 2.4 radii round `center`."""
    rng = np.random.default_rng(seed)
    return np.asarray(center) + rng.uniform(-1.2, 1.2, (n, 3)) * radius


@pytest.mark.parametrize("radius_kpc", [2.0, 4.0, 6.0, 8.0, 10.0, 12.0])
def test_a_hot_star_nebula_has_a_shape_wherever_it_is_in_the_disk(shape, radius_kpc):
    """MAP.423: the model's young density is under the form level from about 6 kpc out, so a nebula read from the model alone is empty near the
    Sun and beyond; anchored on its own sphere it fills a good part of it at every radius."""
    center, radius = (radius_kpc * 1000.0, 0.0, 0.0), 30.0
    points = _ball(center, radius, 3000)
    held = nd.inside(SEED, shape, points, center, radius, "young")
    in_sphere = np.linalg.norm(points - np.asarray(center), axis=1) < radius
    share = held.sum() / max(1, in_sphere.sum())
    assert 0.15 < share < 0.85, share
    assert not held[~in_sphere].any()


def test_without_the_anchor_the_young_density_leaves_the_solar_circle_empty(shape):
    """The defect the anchor fixes: on the solar circle the model alone never reaches the form level."""
    points = _ball((8000.0, 0.0, 0.0), 30.0, 3000)
    assert (nd.density_at(SEED, shape, points, "young", nd.CONTAIN_DETAIL) >= nd.FORM_LEVEL).mean() < 0.02


def test_the_contour_is_not_a_sphere_and_is_the_same_every_time(shape):
    center, radius = (9000.0, 0.0, 0.0), 40.0
    points = _ball(center, radius, 4000)
    first = nd.inside(SEED, shape, points, center, radius, "young")
    assert (first == nd.inside(SEED, shape, points, center, radius, "young")).all()
    sphere = np.linalg.norm(points - np.asarray(center), axis=1) < radius
    assert (sphere & ~first).sum() > 0.1 * sphere.sum() and first.sum() > 0       # lumpy: part of the sphere is outside, part inside


def test_a_dark_cloud_is_held_by_the_cloud_field_density(shape):
    """A field cloud's contour is the gas density the map draws it from; it needs no anchor."""
    found = []
    for i in range(100, 124):
        for j in range(-3, 4):
            for k in range(-1, 2):
                found += [(c, ly_to_pc(n.radius_ly)) for n, c in nebula_field.cell_clouds(SEED, shape, (i, j, k)) if ly_to_pc(n.radius_ly) >= 20.0]
    assert found
    center, radius = found[0]
    points = _ball(center, radius, 2000)
    held = nd.inside(SEED, shape, points, center, radius, "gas")
    in_sphere = np.linalg.norm(points - np.asarray(center), axis=1) < radius
    assert held.any() and not held[~in_sphere].any()
    assert (held == (in_sphere & (nd.density_at(SEED, shape, points, "gas", nd.CONTAIN_DETAIL) >= nd.FORM_LEVEL))).mean() > 0.999


def test_a_star_placed_inside_a_nebula_lies_inside_the_outline_the_map_draws(shape):
    """The map's points of a hot-star nebula, read at the containment test's octaves, are the points the test holds (inside the window where the
    drawing is at full strength)."""
    center, radius = (8000.0, 0.0, 0.0), 40.0
    nebula = {"x": center[0], "y": center[1], "z": center[2], "radius_pc": radius, "family": "emission"}
    drawn = nebula_volume.sample_group([nebula], None, lambda p, f, d, a: nd.density_at(SEED, shape, p, f, d, a), nd.FORM_LEVEL, nd.CONTAIN_DETAIL,
                                       FAMILIES, level=3)
    assert len(drawn.ix) > 100
    points = np.stack([drawn.ix, drawn.iy, drawn.iz], axis=1) * drawn.spacing_pc
    core = np.linalg.norm(points - np.asarray(center), axis=1) < nebula_volume.WINDOW_FROM * radius
    strong = drawn.strength >= 128
    held = nd.inside(SEED, shape, points, center, radius, "young")
    assert core.any() and (strong & core).any()
    assert held[strong & core].mean() > 0.99, "a point drawn at half strength or more is inside the contour"


def test_the_map_draws_a_hot_star_nebula_in_the_outer_disk(shape):
    """MAP.423: a 40 pc emission nebula at 10 kpc drew no point at all; it now draws a body that is full strength in its core."""
    nebula = {"x": 10000.0, "y": 0.0, "z": 0.0, "radius_pc": 40.0, "family": "emission"}
    drawn = nebula_volume.sample_group([nebula], None, lambda p, f, d, a: nd.density_at(SEED, shape, p, f, d, a), nd.FORM_LEVEL, nd.MAX_DETAIL,
                                       FAMILIES)
    assert len(drawn.ix) > 100 and int(drawn.strength.max()) >= 250 and int((drawn.strength >= 128).sum()) > 40


def test_a_container_holds_a_point_only_inside_its_contour(shape):
    center, radius = (9000.0, 0.0, 0.0), 40.0
    seen = []
    inside = containment.contour_test(lambda: (SEED, shape), "D", center, radius)
    points = _ball(center, radius, 300)
    answers = [inside(tuple(p)) for p in points]
    assert any(answers) and not all(answers)
    container = {"column": "inside_nebula_id", "id": 1, "center": center, "radius_pc": radius, "inside": inside}
    for point, answer in zip(points, answers):
        found = containment.innermost_container(tuple(point), [container])
        near = np.linalg.norm(point - np.asarray(center)) <= radius
        assert (found is container) == (answer and near)
    del seen


def test_a_database_never_planned_holds_the_whole_sphere():
    inside = containment.contour_test(lambda: None, "D", (0.0, 0.0, 0.0), 10.0)
    assert inside((1.0, 2.0, 3.0))


def test_the_form_level_comparison_is_rounded(shape):
    """The comparison is made on the density rounded to `CONTAIN_ROUND`, so a last-bit difference cannot move a star across the outline."""
    center, radius = (9000.0, 0.0, 0.0), 40.0
    points = _ball(center, radius, 500)
    value = nd.density_at(SEED, shape, points, "young", nd.CONTAIN_DETAIL, (np.asarray(center), radius))
    bumped = (value.astype(np.float64) + 4e-6)
    rounded = np.rint(bumped / nd.CONTAIN_ROUND) >= round(nd.FORM_LEVEL / nd.CONTAIN_ROUND)
    original = np.rint(value.astype(np.float64) / nd.CONTAIN_ROUND) >= round(nd.FORM_LEVEL / nd.CONTAIN_ROUND)
    assert (rounded != original).mean() < 0.01
