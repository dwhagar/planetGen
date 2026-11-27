# tests/test_gen161_cluster_floor.py

"""
GEN.161: a globular cluster places its own giants and its horizontal branch at plan time whatever the galaxy's luminosity floor is.
The scatter's default floor (18,000 solar luminosities) leaves a cluster entirely to its sectors' fill; a cluster's floor is the faint edge of
the horizontal branch (`tuning.CLUSTER_BRIGHT_FLOOR_SOL`, 35), taken in one place (`star_population.population_floor`) so the pass, the
fill's share and the fill's luminosity cap agree and no star is drawn twice.
"""

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import component_field, profile
from planetgen.generation import bright_stars, component_pass, star_population
from planetgen.util import draw
from tests import test_gen_component_pass as passes
from tests.test_satellites import E1, EDGE, SHAPE
from tests import test_satellites as satellite_tests

DEFAULT_FLOOR = tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL
COLUMNS = store.BRIGHT_STAR_COLUMNS
SOLAR_LUMINOSITY_W = 3.828e26


def test_the_floor_of_a_cluster_is_the_faint_edge_of_the_horizontal_branch_and_only_a_clusters():
    assert tuning.CLUSTER_BRIGHT_FLOOR_SOL == tuning.HORIZONTAL_BRANCH_LUMINOSITY_RANGE_SOL[0] == 35.0
    assert star_population.population_floor("cluster", DEFAULT_FLOOR) == 35.0
    assert star_population.population_floor("cluster", 20.0) == 20.0, "a floor already lower than the cluster's is kept"
    for other in ("old", "young", "satellite", "halo", None):
        assert star_population.population_floor(other, DEFAULT_FLOOR) == DEFAULT_FLOOR
    assert star_population.population_floor("cluster", None) is None


def test_the_share_a_cluster_keeps_is_all_of_its_stars_over_the_floor_whatever_the_galaxys_floor_and_a_band_adds_none():
    shares = bright_stars.band_fractions(DEFAULT_FLOOR)
    assert shares["cluster"] == pytest.approx(star_population.bright_star_fraction(35.0, "cluster"))
    assert 1.0e-4 < shares["cluster"] < 5.0e-3
    assert shares["old"] == pytest.approx(star_population.bright_star_fraction(DEFAULT_FLOOR, "old"))
    assert shares["old"] < shares["cluster"] / 1000.0
    assert bright_stars.band_fractions(1000.0, DEFAULT_FLOOR)["cluster"] == 0.0, "the top pass placed every cluster star over 35 already"
    assert shares["cluster"] == pytest.approx(bright_stars.band_fractions(100.0)["cluster"]), "a lower galaxy floor above the cluster's changes nothing"


def test_a_default_floor_galaxy_places_each_clusters_giants_and_its_whole_horizontal_branch():
    cluster = passes._cluster(mass=2.0e5, half=5.0)
    rows = [dict(zip(COLUMNS, row)) for row in passes._stars([cluster], floor=DEFAULT_FLOOR)]
    expected = component_pass.systems_of(cluster) * bright_stars.band_fractions(DEFAULT_FLOOR)["cluster"]
    assert len(rows) > 150 and abs(len(rows) - expected) < 5.0 * expected ** 0.5
    assert all(row["luminosity_w"] / SOLAR_LUMINOSITY_W >= 35.0 - 1.0e-9 for row in rows)
    branch = [row for row in rows if "Horizontal" in row["star_type"]]
    giants = [row for row in rows if "Giant" in row["star_type"]]
    assert len(branch) > 30 and len(giants) > 30, "both the branch and the giant branches are there"
    assert max(row["temperature_k"] for row in branch) > 7500.0, "the colour spread reaches the blue horizontal branch (metal-poor)"
    assert {row["population"] for row in rows} == {"cluster"}


def test_the_horizontal_branch_count_follows_the_branchs_lifetime():
    """The branch's count is the model's own: its stars are the cluster's stars on the branch, so the share is `HORIZONTAL_BRANCH_DURATION_GY` over the age window (row PS2)."""
    cluster = passes._cluster(mass=2.0e5, half=5.0)
    rows = [dict(zip(COLUMNS, row)) for row in passes._stars([cluster], floor=DEFAULT_FLOOR)]
    branch = sum(1 for row in rows if "Horizontal" in row["star_type"])
    assert 0.2 < branch / len(rows) < 0.7


def test_the_fill_takes_what_the_pass_placed_off_a_clusters_count_and_caps_its_own_stars_under_the_floor():
    found = satellite_tests._placed()
    cluster = passes._cluster(center=(8000.0, 0.0, 0.0), mass=3.0e6, half=4.0)
    shape = SHAPE._replace(components=component_field.ComponentField(list(found) + [cluster], EDGE, E1))
    where = profile.profile_at(cluster.center_pc, shape, EDGE, E1)
    fill = bright_stars.FillContext(where, min_luminosity_sol=DEFAULT_FLOOR)
    assert fill.densities["cluster"] > 0.0
    share = fill.bright_share()
    wanted = sum(density * star_population.bright_star_fraction(DEFAULT_FLOOR, name) for name, density in fill.densities.items()) / sum(fill.densities.values())
    assert share == pytest.approx(wanted), "the fill counts the same share the pass placed"
    assert share > 100.0 * star_population.bright_star_fraction(DEFAULT_FLOOR, "old"), "the cluster's giants and branch count as placed, not only the disk's few"
    capped = {}
    for _ in range(60):
        config = fill.apply(_config(), draw.Stream(f"fill-{_}"))
        capped.setdefault(config.POPULATION, set()).add(config.MAX_STAR_LUMINOSITY_SOL)
    assert capped["cluster"] == {35.0}
    for name, caps in capped.items():
        if name != "cluster":
            assert caps == {DEFAULT_FLOOR}, name


def _config():
    from planetgen.generation.config import SystemConfig
    return SystemConfig()
