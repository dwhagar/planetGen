# tests/test_nebula_density.py

"""
MAP.365 (the pure side): `galaxy.nebula_density.density_at`, the density a nebula is drawn from on the Galaxy Map. A pure
function of position, seed and shape in 0..1; the same point agrees in every call (so nodes and sectors meet without a
seam); octaves refine and do not move the coarse form; a cloud is not a circle; a superbubble clears the gas; the object
ID is never read.
"""

import math

import numpy as np
import pytest
from scipy import ndimage

from planetgen.galaxy import density, nebula_density as nd, nebula_field, superbubble_field

SEED = bytes(range(16))
OTHER_SEED = bytes(range(16, 32))
GRID = 64


@pytest.fixture(scope="module")
def shape():
    return density.build_galaxy_shape(
        disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
        arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4, nuclear_amplitude=1.0)


@pytest.fixture(scope="module")
def clouds(shape):
    """Kept clouds of 20 pc or more, a few thousand parsecs out on the plane: rows `(x, y, z, radius_pc)`."""
    found = []
    for i in range(100, 130):
        for j in range(-5, 6):
            for k in range(-1, 2):
                found += [tuple(row) for row in nd._cell_spheres(SEED, shape, (i, j, k)) if row[3] >= 20.0]
    assert len(found) >= 14
    return found


def _slice(cloud, scale=1.8, n=GRID):
    x, y, z, radius = cloud
    axis = np.linspace(-scale * radius, scale * radius, n)
    gx, gy = np.meshgrid(axis, axis, indexing="ij")
    return np.stack([x + gx, y + gy, np.full_like(gx, z)], axis=-1)


def _inside(shape, cloud, detail, field="gas"):
    return nd.density_at(SEED, shape, _slice(cloud), field, detail) >= nd.FORM_LEVEL


def _radial_spread(mask):
    """The standard deviation over the mean of the distance from the shape's centroid to its boundary pixels: 0 for a circle."""
    boundary = np.argwhere(mask & ~ndimage.binary_erosion(mask))
    if not len(boundary):
        return 0.0
    distance = np.linalg.norm(boundary - np.argwhere(mask).mean(axis=0), axis=1)
    return float(distance.std() / distance.mean())


# --- the contract --------------------------------------------------------------------------------------------------

def test_the_values_are_float32_in_0_to_1_in_the_shape_of_the_points(shape, clouds):
    points = _slice(clouds[0])
    for field in nd.FIELDS:
        for detail in (0, 3, nd.MAX_DETAIL, nd.MAX_DETAIL + 5, -2):
            values = nd.density_at(SEED, shape, points, field, detail)
            assert values.dtype == np.float32 and values.shape == points.shape[:-1]
            assert values.min() >= 0.0 and values.max() <= 1.0 and np.isfinite(values).all()
    assert nd.density_at(SEED, shape, [5000.0, 0.0, 0.0], "gas", 2).shape == ()
    assert nd.density_at(SEED, shape, np.zeros((0, 3)), "gas", 2).shape == (0,)


def test_a_bad_field_or_points_are_refused(shape):
    with pytest.raises(ValueError, match="field is one of"):
        nd.density_at(SEED, shape, [[0.0, 0.0, 0.0]], "dust", 0)
    with pytest.raises(ValueError, match="last axis of 3"):
        nd.density_at(SEED, shape, [[0.0, 0.0]], "gas", 0)


def test_the_wavelengths_halve_and_the_detail_is_clamped(shape, clouds):
    assert nd.wavelength_pc(0) == nd.W0_PC and nd.wavelength_pc(3) == nd.W0_PC / 8
    points = _slice(clouds[0])
    assert np.array_equal(nd.density_at(SEED, shape, points, "gas", nd.MAX_DETAIL), nd.density_at(SEED, shape, points, "gas", 50))
    assert np.array_equal(nd.density_at(SEED, shape, points, "gas", 0), nd.density_at(SEED, shape, points, "gas", -1))


# --- one point, one value: no seam between nodes ------------------------------------------------------------------

@pytest.mark.parametrize("field", nd.FIELDS)
def test_a_point_has_the_same_value_whatever_else_is_asked_with_it(shape, clouds, field):
    points = _slice(clouds[0]).reshape(-1, 3)
    whole = nd.density_at(SEED, shape, points, field, 4)
    pick = np.random.default_rng(1).choice(len(points), 40, replace=False)
    assert np.array_equal(nd.density_at(SEED, shape, points[pick], field, 4), whole[pick])
    assert np.array_equal(nd.density_at(SEED, shape, points[::-1], field, 4)[::-1], whole)


def test_two_grids_that_meet_on_a_border_agree_with_one_grid_over_both(shape, clouds):
    """A sector border through the middle of a cloud: the nodes on each side ask for their own grid, and the values on both sides
    (and on the border plane they share) are the ones one grid over the whole cloud gives, so no seam shows."""
    whole_points = _slice(clouds[0], n=64)
    whole = nd.density_at(SEED, shape, whole_points, "gas", 5)
    left, right = whole_points[:33], whole_points[32:]          # both hold row 32, the border
    assert np.array_equal(nd.density_at(SEED, shape, left, "gas", 5), whole[:33])
    assert np.array_equal(nd.density_at(SEED, shape, right, "gas", 5), whole[32:])


def test_the_same_galaxy_gives_the_same_values_and_another_seed_another_shape(shape, clouds):
    points = _slice(clouds[0])
    first = nd.density_at(SEED, shape, points, "gas", 5)
    assert np.array_equal(first, nd.density_at(SEED, shape, points, "gas", 5))
    other = nd.density_at(OTHER_SEED, shape, points, "gas", 5)
    assert not np.array_equal(first, other)


def test_the_object_id_is_never_read(shape, clouds, monkeypatch):
    def refuse(*args, **kwargs):
        raise AssertionError("the shape read an object ID")
    for name in ("cloud_object_id", "birth_sector", "centroid_mpc"):
        monkeypatch.setattr(nebula_field, name, refuse)
    nd._cell_spheres.cache_clear()
    nd.density_at(SEED, shape, _slice(clouds[0]), "gas", 4)


# --- the look Boss asked for ---------------------------------------------------------------------------------------

def test_a_clouds_core_is_dense_and_the_edge_is_not_a_hard_step(shape, clouds):
    for cloud in clouds[:6]:
        values = nd.density_at(SEED, shape, _slice(cloud), "gas", 4)
        centre = values[GRID // 2 - 3:GRID // 2 + 3, GRID // 2 - 3:GRID // 2 + 3]
        assert centre.mean() > 0.6                                                  # a dense core
        edge = np.abs(np.diff(values, axis=0))
        assert edge.max() < 0.5 and 0.05 < values.mean() < 0.95                     # soft falloff, not a mask


def test_the_outline_of_a_cloud_is_not_a_circle(shape, clouds):
    spreads = [_radial_spread(_inside(shape, cloud, 4)) for cloud in clouds[:12]]
    assert sum(s > 0.15 for s in spreads) >= 8, spreads                           # a circle's is under 0.02
    assert max(_radial_spread(_inside(shape, cloud, 0)) for cloud in clouds[:12]) > 0.15      # bulbous even at the coarsest level


def test_more_octaves_add_filaments_without_moving_the_coarse_form(shape, clouds):
    overlap = np.zeros((nd.MAX_DETAIL, len(clouds[:12])))
    for c, cloud in enumerate(clouds[:12]):
        masks = [_inside(shape, cloud, d) for d in range(nd.MAX_DETAIL + 1)]
        for d in range(nd.MAX_DETAIL):
            overlap[d, c] = (masks[d] & masks[d + 1]).sum() / (masks[d] | masks[d + 1]).sum()
    assert overlap.mean(axis=1)[0] > 0.6 and (overlap.mean(axis=1)[1:3] > 0.75).all() and (overlap.mean(axis=1)[3:] > 0.9).all()
    assert overlap.min() > 0.45, overlap


def test_adding_an_octave_changes_a_value_by_no_more_than_its_amplitude(shape, clouds):
    points = _slice(clouds[1])
    for detail in range(nd.MAX_DETAIL):
        change = np.abs(nd.density_at(SEED, shape, points, "gas", detail + 1) - nd.density_at(SEED, shape, points, "gas", detail))
        assert change.max() <= nd.octave_bound(detail) + 1e-6
    assert nd.octave_bound(1) == pytest.approx(nd.octave_bound(0) * nd.PERSISTENCE)


def test_finer_detail_has_more_edge_than_coarse(shape, clouds):
    """Filaments: the boundary of the inside set is longer at 6 octaves than at 1 (the same area, more edge)."""
    longer = 0
    for cloud in clouds[:12]:
        coarse, fine = _inside(shape, cloud, 1), _inside(shape, cloud, 6)
        edge = lambda m: int((m & ~ndimage.binary_erosion(m)).sum())
        longer += edge(fine) > edge(coarse)
    assert longer >= 9


# --- what the density is made from ---------------------------------------------------------------------------------

def test_a_placed_cloud_always_has_a_core_and_empty_space_has_none(shape, clouds):
    for cloud in clouds[:12]:
        x, y, z, _r = cloud
        assert nd.density_at(SEED, shape, [x, y, z], "gas", 0) >= nd.FORM_LEVEL
    far = np.array([[0.0, 0.0, 20000.0], [30000.0, 0.0, 0.0], [0.0, 0.0, -5000.0]])
    assert (nd.density_at(SEED, shape, far, "gas", 6) < nd.FORM_LEVEL).all()


def test_a_superbubble_clears_the_gas_inside_it(shape):
    bubbles = superbubble_field.all_bubbles(SEED, shape, 9000.0, 200.0)
    assert bubbles
    bubble = max(bubbles, key=lambda b: b.radius_pc)
    centre = np.array(bubble.center_pc)
    inside = centre + np.array([[0.0, 0.0, 0.0], [bubble.radius_pc * 0.3, 0.0, 0.0], [0.0, bubble.radius_pc * 0.5, 0.0]])
    assert (nd.density_at(SEED, shape, inside, "gas", 0) == 0.0).all()


def test_the_model_gas_shapes_the_outline_the_young_field_follows_the_young_stars(shape):
    inner, solar, outer = (np.array([[r, 0.0, 0.0]]) for r in (1500.0, 8000.0, 25000.0))
    young = [float(nd.density_at(SEED, shape, p, "young", 0)[0]) for p in (inner, solar, outer)]
    assert young[0] > young[1] > young[2] >= 0.0
    plane, above = np.array([[8000.0, 0.0, 0.0]]), np.array([[8000.0, 0.0, 1500.0]])
    assert nd.density_at(SEED, shape, plane, "young", 0)[0] > nd.density_at(SEED, shape, above, "young", 0)[0]
    assert (nd.density_at(SEED, shape, np.array([[8000.0, 0.0, 0.0]]) + np.zeros((5, 3)), "young", 6) <= 1.0).all()


def test_the_constants_are_consistent():
    assert 0.0 < nd.FORM_LEVEL < 1.0 and nd.YOUNG_CEILING < 1.0
    assert nd.GAS_BACKGROUND * 1.0 + nd.octave_bound(0) + nd.octave_bound(1) < 1.5       # filaments can reach out, not fill the sky
    assert nd.W0_PC / 2 ** (nd.MAX_DETAIL - 1) < 1.0 + 1e-9 + 1.0
