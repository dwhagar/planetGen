# tests/test_component_field.py

"""
GEN.160: the globular clusters' density in the sector gate (`galaxy/component_field.py`). A cluster's cell density is
the King profile integrated over the cell's cube, so a cluster's cells add up to its mass over 0.4 and the four
exemplars of docs/design/globular-clusters.md section 3 land where its table puts them; the density reaches
`relative_density`, `population_densities`, the profile and the skeleton's bound, and no cell the gate would let in is
left outside the outline.
"""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import component_field, components, profile, skeleton, spherical
from planetgen.galaxy import density as galaxyDensity
from planetgen.physics import king

EDGE = 4.0
E1 = 6.3
"""float: Systems a cell holds at density 1, about what the default edge gives."""

EXEMPLARS = {
    # name: (mass, r_h, W0, central cell range (low, high), cells with at least one system)
    "typical": (2.0e5, 5.0, 5, (3.0e4, 8.4e4), 965),
    "47 Tuc": (8.0e5, 6.0, 9, (1.1e5, 5.5e5), 5848),
    "sparse": (5.0e4, 8.0, 3, (3.4e3, 5.1e3), 1099),
}


def _cluster(mass, half, w0, center=(0.0, 0.0, 0.0), name="Test"):
    return components.Component(
        kind="cluster", name=name, source="generated", center_pc=center, velocity_kms=(0.0, 0.0, 0.0), epoch_unix=0.0,
        mass_solar=mass, half_mass_radius_pc=half, profile="king", age_gy=11.0, metallicity_fe_h=-1.5, seed=1,
        king_w0=w0, tidal_radius_pc=king.tidal_radius_pc(w0, half))


def _lattice(field, rt, shift=0.0):
    """Expected systems of every cell on a lattice of cells covering the tidal sphere, the cluster `shift` pc from a
    cell's middle on each axis."""
    n = int(rt / EDGE) + 3
    axis = range(-n, n + 1)
    return np.array([[[field.relative_density((EDGE * i - shift, EDGE * j - shift, EDGE * k - shift)) * E1
                       for k in axis] for j in axis] for i in axis])


@pytest.mark.parametrize("name", EXEMPLARS)
def test_the_cells_of_a_cluster_add_up_to_its_systems(name):
    mass, half, w0, _central, _count = EXEMPLARS[name]
    field = component_field.ComponentField([_cluster(mass, half, w0)], EDGE, E1)
    cells = _lattice(field, king.tidal_radius_pc(w0, half))
    assert cells.sum() == pytest.approx(king.expected_systems(mass), rel=0.03)


@pytest.mark.parametrize("name", EXEMPLARS)
def test_the_exemplars_land_where_the_design_notes_table_puts_them(name):
    mass, half, w0, (low, high), count = EXEMPLARS[name]
    field = component_field.ComponentField([_cluster(mass, half, w0)], EDGE, E1)
    rt = king.tidal_radius_pc(w0, half)
    middle, corner = _lattice(field, rt), _lattice(field, rt, shift=EDGE / 2)
    assert middle.max() == pytest.approx(high, rel=0.12)       # centred in a cell
    assert corner.max() == pytest.approx(low, rel=0.12)        # on a cell corner
    for cells in (middle, corner):
        assert int((cells >= 1.0).sum()) == pytest.approx(count, rel=0.12)


def test_a_cell_far_from_every_cluster_holds_none_and_costs_a_lookup():
    field = component_field.ComponentField([_cluster(2.0e5, 5.0, 5)], EDGE, E1)
    assert field.relative_density((5000.0, 0.0, 0.0)) == 0.0
    assert field.relative_density((0.0, 0.0, 200.0)) == 0.0
    assert field.relative_density((0.0, 0.0, 0.0)) > 1e3


def test_the_density_follows_the_kings_table():
    """Far from the core a cell is small against the profile, so its density is the table's at its distance."""
    mass, half, w0 = 2.0e5, 5.0, 5
    field = component_field.ComponentField([_cluster(mass, half, w0)], EDGE, E1)
    for radius in (14.0, 18.0, 22.0):
        got = field.relative_density((radius, 0.0, 0.0)) * E1 / EDGE ** 3
        assert got == pytest.approx(king.density_per_pc3(w0, mass, half, radius) / tuning.CLUSTER_MEAN_SYSTEM_MASS_SOL,
                                    rel=0.05)
    assert field.relative_density((king.tidal_radius_pc(w0, half) + 10.0, 0.0, 0.0)) == 0.0


def test_clusters_add_and_a_shape_without_any_is_untouched():
    one = _cluster(2.0e5, 5.0, 5, center=(30.0, 0.0, 0.0), name="A")
    two = _cluster(2.0e5, 5.0, 5, center=(-30.0, 0.0, 0.0), name="B")
    both = component_field.ComponentField([one, two], EDGE, E1)
    alone = component_field.ComponentField([one], EDGE, E1)
    other = component_field.ComponentField([two], EDGE, E1)
    for point in ((30.0, 0.0, 0.0), (0.0, 0.0, 0.0), (-33.0, 2.0, 1.0)):
        assert both.relative_density(point) == pytest.approx(alone.relative_density(point) + other.relative_density(point))
    assert len(both) == 2 and both.names == ("A", "B")


def _shape(found=()):
    shape = galaxyDensity.build_galaxy_shape(
        disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=700.0, bulge_amplitude=2.0,
        arm_count=2, pitch_angle_rad=math.radians(12.0), arm_amplitude=0.3)
    if found:
        shape = shape._replace(components=component_field.ComponentField(found, EDGE, E1))
    return shape


def test_the_cluster_reaches_the_density_model_and_the_profile():
    plain = _shape()
    shape = _shape([_cluster(2.0e5, 5.0, 5, center=(3000.0, 0.0, 1500.0))])
    core = (3000.0, 0.0, 1500.0)
    assert plain.components is None
    assert galaxyDensity.relative_density(core, shape) > galaxyDensity.relative_density(core, plain) + 1e3
    mix = galaxyDensity.population_densities(core, shape)
    assert mix["cluster"] > 1e3
    assert sum(mix.values()) == pytest.approx(galaxyDensity.relative_density(core, shape), rel=1e-9)
    assert galaxyDensity.predicted_star_count(core, shape, E1) == pytest.approx(
        E1 * galaxyDensity.relative_density(core, shape))
    where = profile.profile_at(core, shape, EDGE, E1)
    assert where.region == "cluster" and where.populations["cluster"].share > 0.99
    assert where.expected_systems == pytest.approx(E1 * where.density)
    away = (5000.0, 0.0, 0.0)
    assert galaxyDensity.relative_density(away, shape) == pytest.approx(galaxyDensity.relative_density(away, plain))
    assert galaxyDensity.population_densities(away, shape)["cluster"] == 0.0


def test_the_one_point_densities_sum_to_the_total_with_clusters():
    shape = _shape([_cluster(2.0e5, 5.0, 5, center=(3000.0, 0.0, 1500.0))])
    for point in ((3000.0, 0.0, 1500.0), (3004.0, 0.0, 1500.0), (3000.0, 20.0, 1490.0), (100.0, 0.0, 0.0)):
        assert sum(galaxyDensity.population_densities(point, shape).values()) == pytest.approx(
            galaxyDensity.relative_density(point, shape), rel=1e-9)


def test_a_shape_with_components_compares_and_hashes_by_its_field():
    found = [_cluster(2.0e5, 5.0, 5)]
    assert _shape() == _shape()
    assert _shape(found) != _shape(found)          # another field object
    shape = _shape(found)
    assert shape == shape and hash(shape) == hash(shape)
    assert "components" not in galaxyDensity.shape_with_terms(shape)


def test_the_bound_is_never_below_a_cells_density():
    found = [_cluster(2.0e5, 5.0, 5, center=(3000.0, 0.0, 1500.0))]
    shape = _shape(found)
    rng = np.random.default_rng(3)
    for _ in range(300):
        point = np.array([3000.0, 0.0, 1500.0]) + rng.normal(0.0, 12.0, 3)
        bound = skeleton.bound_relative_density_at(shape, math.hypot(point[0], point[1]), point[2])
        assert bound >= galaxyDensity.relative_density(tuple(point), shape) * (1 - 1e-9)


def test_no_cell_above_the_gate_is_left_outside_the_outline():
    """A cluster 12 kpc off the plane is no part of the run of bands the disk qualifies: every cell the gate lets
    in is still a candidate, by the cluster's own bands."""
    center = (3000.0, 0.0, 12000.0)
    shape = _shape([_cluster(2.0e5, 5.0, 5, center=center)])
    threshold = 1.0 / E1
    extents, _outer, _confirmed = skeleton.build_shell_extents(shape, EDGE, threshold, max_shell=3400)
    high = dict(extents)
    plain, _o, _c = skeleton.build_shell_extents(_shape(), EDGE, threshold, max_shell=3400)
    assert dict(plain) != high                    # the cluster changed the outline
    cells = spherical.cells_touching_sphere(center, king.tidal_radius_pc(5, 5.0) + 6.0, EDGE)
    gated = 0
    for cell in cells:
        address = tuple(cell[:3])
        shell, band = address[:2]
        position = spherical.sector_position_pc(*address, EDGE)
        if galaxyDensity.relative_density(position, shape) >= threshold:
            gated += 1
            assert shell in high and abs(band) <= high[shell], address
    assert gated > 100


def test_a_cluster_outside_the_walked_shells_does_not_extend_the_outline():
    far = _cluster(2.0e5, 5.0, 5, center=(60000.0, 0.0, 0.0))
    threshold = 1.0 / E1
    with_far, outer, _ = skeleton.build_shell_extents(_shape([far]), EDGE, threshold, max_shell=20000)
    without, outer_plain, _ = skeleton.build_shell_extents(_shape(), EDGE, threshold, max_shell=20000)
    assert outer == outer_plain and with_far == without


def test_the_field_is_made_once_from_equal_rows():
    from planetgen.db import store
    found = tuple(_cluster(2.0e5, 5.0, 5, name=f"C{i}", center=(100.0 * i, 0.0, 0.0)) for i in range(3))
    assert store._field_of(found, EDGE, E1) is store._field_of(tuple(found), EDGE, E1)


def test_an_unchanged_table_is_not_read_again(monkeypatch):
    """GEN.244: a galaxy holds tens of thousands of open clusters, and every page asks for the field; the rows are read
    only when the table's signature changes."""
    from planetgen.db import store

    class Answer:
        def __init__(self, row):
            self.row = row

        def fetchone(self):
            return self.row

    signature = {"db": "d", "n": 3, "first_id": 1, "last_id": 3, "mass": 6.0e5}

    class Conn:
        def execute(self, sql, params=()):
            return Answer({"1": 1} if "information_schema" in sql else dict(signature))

    reads = []
    found = [_cluster(2.0e5, 5.0, 5, name=f"C{i}", center=(100.0 * i, 0.0, 0.0)) for i in range(3)]
    monkeypatch.setattr(store, "get_galaxy_components", lambda conn, kind=None: reads.append(1) or list(found))
    monkeypatch.setattr(store, "_FIELD_BY_SIGNATURE", store.collections.OrderedDict())
    first = store._component_field(Conn(), EDGE, E1)
    assert first is not None and store._component_field(Conn(), EDGE, E1) is first and len(reads) == 1
    signature["last_id"] = 9           # the table was written again
    store._component_field(Conn(), EDGE, E1)
    assert len(reads) == 2
