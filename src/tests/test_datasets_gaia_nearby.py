# tests/test_datasets_gaia_nearby.py

"""GEN.1534: the Gaia DR3 sample within 25 pc (`datasets/gaia_nearby.py`): setup downloads it, only a test and the fit read it, and its reader merges the stars of a system and finds the nearest neighbour."""

import math

import pytest

from planetgen.datasets import catalogue, fetch, gaia_nearby


def test_the_sample_is_fetched_by_setup_but_bulk_generation_does_not_need_it():
    assert gaia_nearby.DATASET in catalogue.ALL and not gaia_nearby.DATASET.required
    piece = gaia_nearby.DATASET.pieces[0]
    assert piece.url.startswith("https://gea.esac.esa.int/tap-server/") and "ORDER+BY+source_id" in piece.url and piece.size == 766879
    assert gaia_nearby.DATASET.fill == "SEP1"


def test_stars_closer_than_the_merge_distance_are_one_system_at_their_mean_place():
    stars = [(0.0, 0.0, 0.0), (0.05, 0.0, 0.0), (0.1, 0.0, 0.0), (5.0, 0.0, 0.0)]       # a chain through the middle one, and a loner
    centres = gaia_nearby.systems(stars)
    assert len(centres) == 2
    assert sorted(round(float(c[0]), 3) for c in centres) == [0.05, 5.0]


def test_the_nearest_neighbour_is_asked_only_of_systems_inside_the_centre_distance_and_found_among_all():
    centres = [(1.0, 0.0, 0.0), (3.0, 0.0, 0.0), (20.0, 0.0, 0.0), (21.5, 0.0, 0.0)]
    near = gaia_nearby.nearest_neighbour_distances_pc(centres)       # only the two inside 15 pc ask
    assert [round(float(d), 6) for d in near] == [2.0, 2.0]
    far = gaia_nearby.nearest_neighbour_distances_pc(centres, centre_pc=22.0)
    assert len(far) == 4 and round(float(far[0]), 6) == 1.5


def test_the_density_is_the_count_over_the_sphere_of_the_sample():
    assert gaia_nearby.system_density_per_pc3([(0.0, 0.0, 0.0)] * 10) == pytest.approx(10 / (4.0 / 3.0 * math.pi * 25.0 ** 3))


def test_the_real_sample_holds_the_counts_the_fit_recalls():
    if fetch.locate(gaia_nearby.DATASET).missing:
        pytest.skip("the Gaia sample is not downloaded (python3 -m planetgen.cli.datasets fetch)")
    assert len(gaia_nearby.positions_pc()) == 5323
    centres = gaia_nearby.systems()
    assert 4800 <= len(centres) <= 4900
    assert 0.07 <= gaia_nearby.system_density_per_pc3(centres) <= 0.078
    nearest = gaia_nearby.nearest_neighbour_distances_pc(centres)
    assert 900 <= len(nearest) <= 1000 and 1.2 <= float(nearest[len(nearest) // 2]) <= 1.45
