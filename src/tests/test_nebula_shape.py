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

from planetgen.util import draw
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


# --- Storage and the API (schema v56) ----------------------------------------------------


def test_shape_columns_round_trip():
    shape = _shape(4)
    scalars, balls = ns.to_columns(shape)
    assert set(scalars) == set(ns.SCALAR_COLUMNS)
    again = ns.from_columns(scalars, list(reversed(balls)))
    assert again.to_dict() == shape.to_dict()


def test_a_nebulas_shape_comes_from_its_properties_not_the_random_state():
    args = ("M", 120.0, 300.0, 15.0, 2.0, "H2")
    state = draw.getstate()
    first = ns.shape_for_nebula(*args)
    assert draw.getstate() == state, "drawing a shape must not move the shared random state"
    assert first.to_dict() == ns.shape_for_nebula(*args).to_dict()
    assert first.to_dict() != ns.shape_for_nebula("M", 121.0, 300.0, 15.0, 2.0, "H2").to_dict()


def _nebula_in_db(mysql_config, **attrs):
    from planetgen.db import store
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.phenomena.nebula import Nebula

    nebula = Nebula(SystemConfig(), name="Shapely")
    for name, value in attrs.items():
        setattr(nebula, name, value)
    conn = store.get_connection(mysql_config)
    try:
        nebula_id = store.insert_nebula(conn, nebula, placement={
            "center_x_pc": 5.0, "center_y_pc": 1.0, "center_z_pc": 2.0, "galactic_radius_pc": 5.5})
        conn.commit()
    finally:
        conn.close()
    return nebula_id, nebula


def test_a_saved_nebula_stores_its_shape(mysql_config):
    from planetgen.db import query, store

    nebula_id, nebula = _nebula_in_db(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        row = conn.execute("SELECT shape_scale, shape_warp_octaves FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
        balls = conn.execute("SELECT COUNT(*) AS n FROM nebula_shape_balls WHERE nebula_id = ?", (nebula_id,)).fetchone()
        _row, loaded = query.nebula_shape(conn, nebula_id)
        # Deleting the nebula takes its balls with it.
        conn.execute("DELETE FROM nebulae WHERE id = ?", (nebula_id,))
        left = conn.execute("SELECT COUNT(*) AS n FROM nebula_shape_balls WHERE nebula_id = ?", (nebula_id,)).fetchone()
        conn.commit()
        with pytest.raises(ValueError):
            query.nebula_shape(conn, nebula_id)
    finally:
        conn.close()
    assert row["shape_scale"] == pytest.approx(nebula.get_shape().scale)
    assert 4 <= balls["n"] <= 8 and left["n"] == 0
    assert loaded.to_dict() == nebula.get_shape().to_dict()


def test_a_nebula_saved_before_shapes_gets_the_shape_a_fresh_save_would(mysql_config):
    from planetgen.db import query, store

    nebula_id, nebula = _nebula_in_db(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM nebula_shape_balls WHERE nebula_id = ?", (nebula_id,))
        conn.execute("UPDATE nebulae SET shape_scale = NULL WHERE id = ?", (nebula_id,))
        conn.commit()
        _row, loaded = query.nebula_shape(conn, nebula_id)
    finally:
        conn.close()
    assert loaded.to_dict() == nebula.get_shape().to_dict()
