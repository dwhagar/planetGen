# tests/test_datasets_real.py

"""
The downloaded datasets of GEN.548 against numbers measured in the sky (GEN.590), on the real files. They are skipped on a
machine whose cache folder (`python3 -m planetgen.cli.datasets status`) does not hold every dataset of the catalogue (TEST.312: a
machine that fetched before a dataset was added, such as the ozone cross sections, lacks it), which is every test machine that has
not run `update.sh`; the offline fixtures of `test_datasets_readers.py` and `test_stellar_composition.py` cover the code.
"""

import pytest

from planetgen.datasets import catalogue, fetch, isochrones, neutron_star_matter, white_dwarf_models
from planetgen.physics import stellar_composition as sc

pytestmark = [pytest.mark.real_datasets, pytest.mark.skipif(
    any(fetch.locate(dataset).missing for dataset in catalogue.ALL),
    reason="the downloaded datasets are not all in this machine's cache folder")]


def test_the_main_sequence_of_a_sun_like_star_lives_about_ten_gigayears():
    point = isochrones.main_sequence_point(1.0, 4.6, 0.0)
    assert point is not None and 10 ** point[0] == pytest.approx(1.1, abs=0.15)
    assert isochrones.main_sequence_lifetime_gy(1.0, 0.0) == pytest.approx(10.0, rel=0.1)
    assert isochrones.main_sequence_lifetime_gy(2.0, 0.0) < isochrones.main_sequence_lifetime_gy(1.0, 0.0) / 4


def test_a_metal_poor_star_is_hotter_than_a_metal_rich_one_of_the_same_mass_and_age():
    poor, rich = isochrones.main_sequence_point(0.9, 5.0, -1.0), isochrones.main_sequence_point(0.9, 5.0, 0.0)
    assert poor[1] > rich[1]


def test_the_sequences_give_sirius_b_its_radius():
    """Sirius B: 1.018 solar masses, 25,000 K, 0.0080 solar radii (Joyce et al. 2018)."""
    assert white_dwarf_models.radius_solar(1.018, 25000.0) == pytest.approx(0.0080, rel=0.06)


def test_a_heavier_white_dwarf_is_smaller_and_a_thick_layer_one_larger():
    assert white_dwarf_models.radius_solar(1.2, 12000.0) < white_dwarf_models.radius_solar(0.6, 12000.0)
    assert white_dwarf_models.radius_solar(0.6, 12000.0, "thick") > white_dwarf_models.radius_solar(0.6, 12000.0, "thin")


def test_every_mass_model_of_the_outer_crust_starts_with_iron_56_and_ends_near_the_drip_density():
    for model in neutron_star_matter.crust_models():
        layers = neutron_star_matter.outer_crust(model)
        assert (layers[0].z, layers[0].a) == (26, 56), model
        drip = neutron_star_matter.mass_density_g_cm3(layers[-1].upper_density_fm3)
        assert 3e11 < drip < 6e11, (model, drip)


def test_a_neutron_star_from_the_real_tables_has_an_iron_surface():
    for n in range(20):
        star = sc.star_makeup(sc.NEUTRON_STAR, None, 0.1, 1.4, f"ns{n}")
        assert star.detail["outer_crust"][0] == "Fe-56" and star.detail["crust_model"] in neutron_star_matter.crust_models()


def test_the_catalogue_reports_nothing_missing():
    assert catalogue.missing() == []
