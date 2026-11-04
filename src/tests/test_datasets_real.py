# tests/test_datasets_real.py

"""
The downloaded datasets of GEN.548 against numbers measured in the sky (GEN.590), on the real files. They are skipped on a
machine whose cache folder (`python3 -m planetgen.cli.datasets status`) does not hold every dataset of the catalogue (TEST.312: a
machine that fetched before a dataset was added, such as the ozone cross sections, lacks it), which is every test machine that has
not run `update.sh`; the offline fixtures of `test_datasets_readers.py` and `test_stellar_composition.py` cover the code.
"""

import pytest

from planetgen.datasets import catalogue, fetch, hhe_eos, isochrones, neutron_star_matter, white_dwarf_models
from planetgen.physics import constants, giant_structure, stellar_composition as sc

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
    """Sirius B: 1.018 solar masses, 25,000 K, 0.0080 solar radii [R2208]."""
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


def test_hydrogen_helium_at_a_megabar_is_a_little_under_one_gram_per_cubic_centimetre():
    """[R0708] at 100 GPa and 5,000 K, Y = 0.275 (Jupiter's mid-envelope is near 1 g/cm3)."""
    assert 600.0 < hhe_eos.density_kg_m3(5000.0, 100e9) < 1100.0
    assert hhe_eos.density_kg_m3(5000.0, 100e9, 0.5) > hhe_eos.density_kg_m3(5000.0, 100e9, 0.1)      # helium is the denser fluid
    assert hhe_eos.entropy_mj_kg_k(8000.0, 100e9) > hhe_eos.entropy_mj_kg_k(4000.0, 100e9)


EARTH = constants.EARTH_MASS_TO_KG


def _fit(mass_earth, radius_km, heavy_earth, z, t1bar, y=0.27, ice=0.0, core_ice=0.0):
    return giant_structure.fit_adiabat(giant_structure.StructureSubject(
        mass_earth * EARTH, radius_km * 1e3, y, envelope_z=z, heavy_mass_kg=heavy_earth * EARTH, temperature_1bar_k=t1bar,
        envelope_ice_fraction=ice, core_ice_fraction=core_ice))


def test_the_ice_giants_close_on_adiabats_near_their_measured_one_bar_temperatures():
    """Uranus (76 K at one bar, [R2699]) and Neptune (72 K): the model has the radius and finds the adiabat itself, so
    for a heavy-element mass and an ice share inside the range of the interior models the adiabat is within a third of what was
    measured (the tables and the radius decide it; the shares are the only choice, and it moves the answer by tens of kelvin)."""
    uranus = _fit(14.54, 25362.0, 12.5, 0.1, 76.0, ice=0.8, core_ice=0.8)
    neptune = _fit(17.15, 24622.0, 14.5, 0.1, 72.0, ice=0.8, core_ice=0.6)
    assert uranus.closed and neptune.closed
    assert 50.0 < uranus.adiabat_1bar_k < 100.0 and 40.0 < neptune.adiabat_1bar_k < 100.0
    assert 1.0e11 < uranus.central_pressure_pa < 4.0e12


def test_a_gas_giant_closes_and_a_lighter_one_is_less_deep():
    jupiter = _fit(317.8, 69911.0, 25.0, 0.04, 166.0)
    saturn = _fit(95.16, 58232.0, 20.0, 0.04, 134.0)
    assert jupiter.closed and saturn.closed and 3.0e12 < jupiter.central_pressure_pa < 1.5e13
    assert saturn.central_pressure_pa < jupiter.central_pressure_pa


def test_a_brown_dwarf_is_degenerate_at_the_centre():
    mass = 40 * 317.8
    dwarf = _fit(mass, 0.9 * 71492.0, 0.015 * mass, 0.015, 1000.0)
    assert dwarf.closed and dwarf.central_pressure_pa > 1.0e13


def test_the_catalogue_reports_nothing_missing():
    assert catalogue.missing() == []
