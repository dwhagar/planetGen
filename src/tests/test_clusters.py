# tests/test_clusters.py

"""
GEN.159: the generated globular cluster system (`generation/clusters.py`) is a pure function of the galaxy seed, has the
count it is given, and each cluster is a well-formed King component whose edge fits inside its Jacobi radius.
"""

import math
import statistics

import pytest

from planetgen import tuning
from planetgen.galaxy import components, enclosed_mass, shape_draw
from planetgen.generation import clusters
from planetgen.physics import king

SEED = bytes(range(16))
COUNT = 157


@pytest.fixture(scope="module")
def system():
    return clusters.generate_clusters(SEED, COUNT, epoch_unix=123.0)


def test_the_count_is_the_one_asked_for(system):
    assert len(system) == COUNT
    assert len(clusters.generate_clusters(SEED, count=12)) == 12


def test_the_same_seed_gives_the_same_clusters_and_another_gives_others(system):
    assert system == clusters.generate_clusters(SEED, COUNT, epoch_unix=123.0)
    other = clusters.generate_clusters(bytes(reversed(SEED)), COUNT, epoch_unix=123.0)
    assert [c.center_pc for c in system] != [c.center_pc for c in other]
    assert clusters.generate_clusters(7, COUNT) == clusters.generate_clusters(7, COUNT)


def test_a_larger_count_keeps_the_first_clusters(system):
    assert clusters.generate_clusters(SEED, count=200, epoch_unix=123.0)[:157] == system


def test_every_cluster_is_a_well_formed_king_component(system):
    low_mass, high_mass = tuning.CLUSTER_MASS_RANGE_SOL
    low_r, high_r = tuning.CLUSTER_HALF_MASS_RADIUS_RANGE_PC
    low_w0, high_w0 = tuning.KING_W0_RANGE
    low_fe, high_fe = tuning.METALLICITY_RANGE_FE_H
    names = set()
    for c in system:
        assert (c.kind, c.source, c.profile, c.epoch_unix) == ("cluster", "generated", "king", 123.0)
        assert low_mass <= c.mass_solar <= high_mass
        assert low_r <= c.half_mass_radius_pc <= high_r
        assert low_w0 <= c.king_w0 <= high_w0 and c.king_w0 == int(c.king_w0)
        assert c.tidal_radius_pc == pytest.approx(king.tidal_radius_pc(c.king_w0, c.half_mass_radius_pc))
        assert low_fe <= c.metallicity_fe_h <= high_fe
        assert tuning.CLUSTER_AGE_RANGE_GY[0] <= c.age_gy <= tuning.CLUSTER_AGE_RANGE_GY[1]
        assert 0 <= c.seed < 2 ** 63
        assert c.id is None and c.parent_component_id is None and c.catalogue_id is None
        names.add(c.name)
    assert len(names) == len(system)


def test_a_clusters_edge_fits_inside_its_jacobi_radius_unless_it_is_already_at_the_lowest_potential(system):
    for c in system:
        distance = math.sqrt(sum(v * v for v in c.center_pc))
        jacobi = king.jacobi_radius_pc(c.mass_solar, distance, enclosed_mass.enclosed_mass_solar(shape_draw.default_shape(SEED), SEED, max(distance, 1.0)))
        assert c.tidal_radius_pc <= jacobi * (1 + 1e-9) or c.king_w0 == tuning.KING_W0_RANGE[0], c


def test_only_the_top_potential_is_core_collapsed(system):
    collapsed = [c for c in system if c.core_collapsed]
    assert collapsed and all(c.king_w0 == tuning.KING_W0_RANGE[1] for c in collapsed)


def test_the_system_has_the_designs_shape():
    big = clusters.generate_clusters(SEED, count=4000)
    assert statistics.median(c.mass_solar for c in big) == pytest.approx(tuning.CLUSTER_MASS_PEAK_SOL, rel=0.1)
    assert statistics.median(c.half_mass_radius_pc for c in big) == pytest.approx(5.0, rel=0.4)
    poor = [c for c in big if c.metallicity_fe_h < -1.0]
    assert len(poor) / len(big) == pytest.approx(tuning.CLUSTER_FE_H_GROUPS[0][0], abs=0.04)
    distance = [math.sqrt(sum(v * v for v in c.center_pc)) for c in big]
    assert max(distance) <= tuning.CLUSTER_HALO_RADIUS_RANGE_PC[1] * 1.001  # inside the cell key's reach, 262,144 pc
    assert sum(c.expected_systems for c in clusters.generate_clusters(SEED, COUNT)) == pytest.approx(1.8e8, rel=0.5)


def test_a_component_checks_its_own_fields():
    base = dict(kind="cluster", name="X", source="generated", center_pc=(0, 0, 0), velocity_kms=(0, 0, 0), epoch_unix=0.0,
                mass_solar=1e5, half_mass_radius_pc=5.0, age_gy=11.0, metallicity_fe_h=-1.5, seed=1)
    components.Component(profile="king", king_w0=7, tidal_radius_pc=40.0, **base)
    components.Component(profile="plummer", tidal_radius_pc=40.0, **base)
    with pytest.raises(ValueError):
        components.Component(profile="king", **base)
    with pytest.raises(ValueError):
        components.Component(profile="plummer", **base)
    with pytest.raises(ValueError):
        components.Component(profile="plummer", king_w0=7, tidal_radius_pc=40.0, **base)
    with pytest.raises(ValueError):
        components.Component(profile="king", king_w0=7, tidal_radius_pc=40.0, **{**base, "kind": "galaxy"})
    with pytest.raises(ValueError):
        components.Component(profile="king", king_w0=7, tidal_radius_pc=40.0, **{**base, "source": "guess"})
    with pytest.raises(ValueError):
        components.Component(profile="nfw", tidal_radius_pc=40.0, **base)
