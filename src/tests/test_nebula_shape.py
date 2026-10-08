"""
GEN.75: a nebula's shape -- seeded metaballs, domain-warped by noise, read
at an isovalue (`planetgen.galaxy.nebula_shape`): the same seed gives the
same shape, it round-trips through JSON, it is bulbous and irregular (not a
sphere), `contains` and the marching-cubes mesh agree, and the surface stays
inside the bounding sphere of radius 1.
"""

import json
import random

import numpy as np
import pytest

from planetgen.galaxy import nebula_shape as ns

SEEDS = range(8)


def _shape(seed):
    return ns.draw_shape(random.Random(seed))


def _unit_ball_points(count=4000, seed=1):
    pts = np.random.default_rng(seed).uniform(-1, 1, (count * 3, 3))
    return pts[np.linalg.norm(pts, axis=1) <= 1][:count]


def test_the_same_seed_draws_the_same_shape_and_others_differ():
    assert _shape(3).to_dict() == _shape(3).to_dict()
    assert _shape(3).to_dict() != _shape(4).to_dict()
    # 4 to 8 centres, each with a radius.
    for seed in SEEDS:
        shape = _shape(seed)
        assert 4 <= len(shape.centres) <= 8 and len(shape.radii) == len(shape.centres)


def test_a_shape_round_trips_through_json():
    shape = _shape(2)
    again = ns.NebulaShape.from_dict(json.loads(json.dumps(shape.to_dict())))
    pts = _unit_ball_points(500)
    assert again.to_dict() == shape.to_dict()
    assert (again.contains_many(pts) == shape.contains_many(pts)).all()
    with pytest.raises(ValueError):
        ns.NebulaShape.from_dict({**shape.to_dict(), "version": 99})


def test_the_noise_is_deterministic_and_stays_in_range():
    noise = ns.GradientNoise(7)
    pts = np.random.default_rng(0).uniform(-20, 20, (2000, 3))
    values = noise(pts)
    assert np.allclose(values, ns.GradientNoise(7)(pts))
    assert not np.allclose(values, ns.GradientNoise(8)(pts))
    assert -1.01 <= values.min() and values.max() <= 1.01
    assert abs(float(noise(np.array([0.0, 0.0, 0.0])))) < 1e-9  # a lattice point


@pytest.mark.parametrize("seed", SEEDS)
def test_the_shape_is_a_bulbous_region_inside_the_unit_sphere(seed):
    shape = _shape(seed)
    pts = _unit_ball_points()
    inside = shape.contains_many(pts)
    fill = inside.mean()
    # Something is there, and it is far from filling its bounding sphere
    # (a sphere would be 1.0).
    assert 0.02 < fill < 0.6, fill
    # Nothing past the bounding sphere counts as inside.
    outside = np.random.default_rng(2).uniform(-1.5, 1.5, (3000, 3))
    outside = outside[np.linalg.norm(outside, axis=1) > 1.0]
    assert not shape.contains_many(outside).any()
    assert not shape.contains([2.0, 0.0, 0.0])


@pytest.mark.parametrize("seed", SEEDS)
def test_the_mesh_lies_on_the_boundary_of_what_contains_says(seed):
    shape = _shape(seed)
    verts, faces = shape.mesh("full")
    assert len(verts) > 100 and len(faces) > 100
    assert faces.min() >= 0 and faces.max() < len(verts)
    assert np.linalg.norm(verts, axis=1).max() <= 1.0
    # Every vertex is on the isosurface, give or take the grid's error
    # (the field is steep at the edge, so the error is in value, not place).
    error = np.abs(shape.field(verts * shape.scale) - shape.iso)
    assert np.median(error) < 0.05 and error.max() < 0.35
    low_verts, low_faces = shape.mesh("low")
    assert 0 < len(low_faces) < len(faces)


def test_the_low_poly_mesh_is_coarser_than_the_full_one_and_not_a_sphere():
    shape = _shape(1)
    low = shape.mesh("low")[1]
    full = shape.mesh("full")[1]
    assert len(low) * 3 < len(full)
    verts, _ = shape.mesh("full")
    radii = np.linalg.norm(verts - verts.mean(axis=0), axis=1)
    assert radii.std() / radii.mean() > 0.1, "an irregular outline, not a sphere"
