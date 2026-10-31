# tests/test_satellites.py

"""
GEN.237: the Milky Way's satellite galaxies as galaxy components. The exponential and Plummer profile tables
(`physics/component_profiles.py`), the table of real satellites and its embedding in the galaxy frame
(`generation/satellites.py`), and what the satellites do in the sector gate and the sector profile (`ComponentField`
with the `satellite` kind, `galaxy/profile.py`, the fill's [Fe/H]).
"""

import dataclasses
import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import component_field, profile, skeleton, spherical
from planetgen.galaxy import density as galaxyDensity
from planetgen.generation import bright_stars
from tests import reference_satellites as satellites
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import component_profiles, metallicity
from planetgen.physics.units import ly_to_pc
from planetgen.galaxy.skeleton import expected_system_count_at_density_1
from planetgen.util import draw

SHAPE = galaxyDensity.build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
    arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)
EDGE = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
E1 = expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY)
SOLAR_RADIUS = tuning.REFERENCE_RADIUS_SCALE_LENGTHS * SHAPE.disk_scale_length_pc
SOLAR_ANGLE = galaxyDensity.model_terms(SHAPE)["reference_angle_rad"]


def _placed():
    return satellites.place_reference_satellites(7, SOLAR_RADIUS, SOLAR_ANGLE, epoch_unix=tuning.SCATTER_PHASE_ORIGIN_UNIX)


def _by_name(found):
    return {c.name: c for c in found}


def _field():
    return component_field.ComponentField(_placed(), EDGE, E1)


@pytest.mark.parametrize("name", ["exponential", "plummer"])
def test_a_profile_table_holds_the_whole_mass_and_puts_half_of_it_at_the_half_mass_radius(name):
    table = component_profiles.profile_table(name)
    assert table.mass_fraction[0] == 0.0 and table.mass_fraction[-1] == 1.0
    assert np.all(np.diff(table.mass_fraction) >= 0.0)
    assert table.mass_within_x(table.half_mass_x) == pytest.approx(0.5, abs=1e-3)
    assert table.density[0] == 1.0 and np.all(np.diff(table.density) <= 0.0)
    assert table.density_at_x(table.tidal_x * 1.01) == 0.0
    # The cut is where the design says: six half-mass radii of the uncut profile, so the table's own half-mass radius
    # is a little under the uncut one.
    assert table.tidal_x / table.half_mass_x == pytest.approx(tuning.SATELLITE_TIDAL_OVER_HALF_MASS, rel=0.05)


def test_an_exponential_profile_integrates_to_its_closed_form():
    table = component_profiles.profile_table("exponential")
    x = np.linspace(0.0, table.tidal_x, 200001)
    integral = np.trapezoid(x ** 2 * np.exp(-x), x)
    assert table.mass_over_central == pytest.approx(integral, rel=1e-4)
    assert table.half_mass_x == pytest.approx(2.674, rel=0.01)


def test_a_plummer_profile_matches_the_half_mass_radius_of_the_uncut_sphere():
    table = component_profiles.profile_table("plummer")
    assert table.half_mass_x == pytest.approx(1.305, rel=0.03)


def test_a_component_is_cut_at_the_radius_its_table_says():
    for row in satellites.MILKY_WAY_SATELLITES:
        table = component_profiles.profile_table(row.profile)
        cut = component_profiles.cut_radius_pc(row.profile, row.half_light_radius_pc)
        assert cut == pytest.approx(row.half_light_radius_pc / table.half_mass_x * table.tidal_x)
    with pytest.raises(ValueError):
        component_profiles.profile_table("nfw")


def test_the_table_of_satellites_holds_the_clouds_at_the_distances_the_design_gives():
    found = _by_name(_placed())
    assert len(found) == len(satellites.MILKY_WAY_SATELLITES) == 12
    # The counts follow the stars' mass over the mean system mass, which the primary mass function (GEN.929) made 5% heavier
    # than the hand value: the design's 5.09e9, 1.07e9, 2.09e8 and 4.15e7 times 0.951.
    # GEN.366: the distances and masses are McConnachie (2012) Tables 2 and 4 (51, 64, 26 and 147 kpc; 1.5e9, 4.6e8, 2.1e7 and 2.0e7 solar masses), so the counts
    # are the design's scaled by the mass ratio: 4.84e9 x 0.6, 1.02e9 x 0.885, 1.99e8 x 0.21.
    for name, distance_kpc, systems in (("Large Magellanic Cloud", 51.0, 2.90e9), ("Small Magellanic Cloud", 64.0, 9.02e8),
                                        ("Sagittarius dwarf", 26.0, 4.18e7), ("Fornax dwarf", 147.0, 3.95e7)):
        component = found[name]
        sun = np.array(satellites.sun_frame(SOLAR_RADIUS, SOLAR_ANGLE)[0])
        assert np.linalg.norm(np.array(component.center_pc) - sun) == pytest.approx(distance_kpc * 1000.0, rel=1e-9)
        assert component.expected_systems == pytest.approx(systems, rel=0.03)
    assert {c.kind for c in found.values()} == {"satellite"} and {c.source for c in found.values()} == {"catalogue"}


def test_a_direction_on_the_sky_is_a_direction_in_the_galaxy():
    sun = satellites.sun_frame(SOLAR_RADIUS, SOLAR_ANGLE)[0]
    assert math.hypot(sun[0], sun[1]) == pytest.approx(SOLAR_RADIUS) and sun[2] == 0.0
    centre = satellites.satellite_position_pc(0.0, 0.0, SOLAR_RADIUS, SOLAR_RADIUS, SOLAR_ANGLE)
    assert centre == pytest.approx((0.0, 0.0, 0.0), abs=1e-6)                  # longitude 0 looks at the galactic centre
    anti = satellites.satellite_position_pc(180.0, 0.0, 1000.0, SOLAR_RADIUS, SOLAR_ANGLE)
    assert math.hypot(anti[0], anti[1]) == pytest.approx(SOLAR_RADIUS + 1000.0)  # and 180 away from it
    north = satellites.satellite_position_pc(0.0, 90.0, 500.0, SOLAR_RADIUS, SOLAR_ANGLE)
    assert north[2] == pytest.approx(-500.0)                                  # the pole points to -z in a proper frame
    # Longitude 90 lies the way the disk turns: counter-clockwise seen from +z, ahead of the Sun in azimuth.
    sun = satellites.sun_frame(SOLAR_RADIUS, SOLAR_ANGLE)[0]
    ahead = satellites.satellite_position_pc(90.0, 0.0, 1000.0, SOLAR_RADIUS, SOLAR_ANGLE)
    assert np.array(ahead) - np.array(sun) == pytest.approx(
        1000.0 * np.array([-math.sin(SOLAR_ANGLE), math.cos(SOLAR_ANGLE), 0.0]), abs=1e-6)


def test_the_satellites_scale_with_the_galaxy_and_keep_their_size():
    toy = satellites.place_reference_satellites(7, SOLAR_RADIUS / 10.0, SOLAR_ANGLE, scale=0.1)
    full = _by_name(_placed())
    for component in toy:
        twin = full[component.name]
        assert np.allclose(component.center_pc, np.array(twin.center_pc) / 10.0)
        assert (component.mass_solar, component.half_mass_radius_pc) == (twin.mass_solar, twin.half_mass_radius_pc)


def test_the_seeds_follow_the_galaxy_seed_and_the_name():
    first, again, other = _placed(), _placed(), satellites.place_reference_satellites(8, SOLAR_RADIUS, SOLAR_ANGLE)
    assert [c.seed for c in first] == [c.seed for c in again]
    assert [c.seed for c in first] != [c.seed for c in other]
    assert len({c.seed for c in first}) == len(first)


def test_a_clouds_profile_integrates_to_its_systems_through_the_field():
    """The density the gate sees, summed over the satellite, is its mass over 0.4: integrate it along a ray."""
    field = _field()
    cloud = _by_name(_placed())["Small Magellanic Cloud"]
    cut = cloud.tidal_radius_pc
    radii = np.geomspace(5.0, cut, 400)
    # relative_density is systems per cell over E1; a cell is EDGE^3, so per pc^3 it is rho * E1 / EDGE^3.
    per_pc3 = np.array([field.relative_density((cloud.center_pc[0] + r, cloud.center_pc[1], cloud.center_pc[2]), "satellite")
                        for r in radii]) * E1 / EDGE ** 3
    # The Large Cloud lies 20 kpc away and reaches into none of the SMC's tidal sphere.
    inner = np.trapezoid(4.0 * math.pi * radii ** 2 * per_pc3, radii)
    inner += per_pc3[0] * 4.0 / 3.0 * math.pi * radii[0] ** 3
    assert inner == pytest.approx(cloud.expected_systems, rel=0.02)


def test_a_satellites_cells_show_up_as_the_satellite_population_with_its_own_metallicity():
    found = _placed()
    shape = SHAPE._replace(components=component_field.ComponentField(found, EDGE, E1))
    cloud = _by_name(found)["Large Magellanic Cloud"]
    where = profile.profile_at(cloud.center_pc, shape, EDGE, E1)
    parts = ("satellite_young", "satellite_intermediate", "satellite")
    assert where.region == "satellite" and sum(where.populations[name].share for name in parts) > 0.99
    assert sum(where.mix[name] for name in parts) > 1.0 and where.mix["cluster"] == 0.0
    assert sum(where.mix.values()) == pytest.approx(where.density, rel=1e-9)
    for name in parts:      # the Cloud's young and old stars share the Cloud's metallicity
        ((weight, mean, sigma),) = where.populations[name].fe_h
        assert (mean, sigma) == (cloud.metallicity_fe_h, tuning.COMPONENT_FE_H_SIGMA["satellite"]) and weight > 0.0
    assert where.populations["satellite"].age_gy == tuning.STELLAR_POPULATION_AGE_RANGES_GY["satellite"]
    assert where.expected_systems == pytest.approx(E1 * where.density)
    # A cell at the Sun has no satellite and takes the generic peaks if one is asked for.
    sun = profile.profile_at((SOLAR_RADIUS, 0.0, 0.0), shape, EDGE, E1)
    assert sun.mix["satellite"] == 0.0
    assert sun.populations["satellite"].fe_h == tuning.METALLICITY_SATELLITE_PEAKS


def test_a_cell_where_two_satellites_overlap_has_a_peak_for_each():
    a, b = _placed()[3:5]                                   # Fornax and Sculptor, moved to one place
    b = dataclasses.replace(b, center_pc=a.center_pc)
    field = component_field.ComponentField([a, b], EDGE, E1)
    peaks = field.fe_h_peaks(a.center_pc, "satellite")
    assert sorted(mean for _w, mean, _s in peaks) == sorted([a.metallicity_fe_h, b.metallicity_fe_h])
    assert field.fe_h_peaks(a.center_pc, "cluster") == ()


def test_the_field_sums_kinds_and_leaves_the_arrays_alone():
    found = _placed()
    field = component_field.ComponentField(found, EDGE, E1)
    cloud = _by_name(found)["Large Magellanic Cloud"].center_pc
    densities = field.densities(cloud)
    assert set(densities) == set(component_field.POPULATIONS_IN_THE_GATE)
    assert densities["cluster"] == 0.0
    assert sum(densities[name] for name in component_field.SATELLITE_POPULATIONS) == pytest.approx(
        field.relative_density(cloud))
    shape = SHAPE._replace(components=field)
    arrays = galaxyDensity.population_densities_array(np.array([cloud[0]]), np.array([cloud[1]]), np.array([cloud[2]]), shape)
    assert all(arrays[name][0] == 0.0 for name in component_field.POPULATIONS_IN_THE_GATE)     # the component pass is GEN.240
    assert galaxyDensity.population_densities(cloud, shape)["satellite"] > 0.5


def test_the_bound_and_the_outline_take_a_satellite_in():
    found = _placed()
    shape = SHAPE._replace(components=component_field.ComponentField(found, EDGE, E1))
    for name in ("Small Magellanic Cloud", "Sagittarius dwarf"):
        centre = _by_name(found)[name].center_pc
        r_cyl = math.hypot(centre[0], centre[1])
        part = shape.components.relative_density(centre)
        assert shape.components.bound(r_cyl, centre[2]) >= part * (1 - 1e-9) > 0.0
        added = skeleton.bound_relative_density_at(shape, r_cyl, centre[2]) - skeleton.bound_relative_density_at(
            SHAPE, r_cyl, centre[2])
        assert added == pytest.approx(shape.components.bound(r_cyl, centre[2]))
        shell, band, _wedge = spherical.sector_address_at(centre, EDGE)
        alone = component_field.ComponentField([_by_name(found)[name]], EDGE, E1)
        if name == "Small Magellanic Cloud":
            assert abs(alone.high_band(shell, spherical.band_at, EDGE)) >= abs(band)
        else:                                  # GEN.627: the Sagittarius dwarf's peak is under one system a cell, so it adds none
            assert alone._density_reach[0] == 0.0 and alone.high_band(shell, spherical.band_at, EDGE) == -1
    assert shape.components.high_band(10, spherical.band_at, EDGE) == -1       # a shell 10 edges out reaches none


def test_a_cloud_star_has_the_clouds_metallicity_and_a_cluster_or_satellite_cannot_be_sampled_without_peaks():
    stream = lambda: draw.Stream("test-satellite-fe-h")
    with pytest.raises(ValueError):
        metallicity.sample_fe_h("satellite", 1.0e5, stream())
    with pytest.raises(ValueError):
        metallicity.sample_fe_h("cluster", 1.0e5, stream())
    peaks = ((1.0, -0.4, 0.3),)
    drawn = [metallicity.sample_fe_h("satellite", 1.0e5, draw.Stream(f"t{i}"), peaks=peaks) for i in range(2000)]
    assert np.mean(drawn) == pytest.approx(-0.4, abs=0.03) and np.std(drawn) == pytest.approx(0.3, abs=0.03)
    two = ((1.0, -2.0, 0.1), (3.0, 0.0, 0.1))
    drawn = np.array([metallicity.sample_fe_h("satellite", 1.0e5, draw.Stream(f"u{i}"), peaks=two) for i in range(2000)])
    assert (drawn < -1.0).mean() == pytest.approx(0.25, abs=0.04)


def test_the_fill_gives_a_satellite_system_the_peaks_of_its_place_and_the_star_draws_from_them():
    found = _placed()
    shape = SHAPE._replace(components=component_field.ComponentField(found, EDGE, E1))
    cloud = _by_name(found)["Small Magellanic Cloud"]
    where = profile.profile_at(cloud.center_pc, shape, EDGE, E1)
    config = bright_stars.FillContext(where).apply(SystemConfig(), draw.Stream("test-fill"))
    assert config.POPULATION == "satellite"
    assert config.FE_H_PEAKS == where.populations["satellite"].fe_h
    drawn = [Star(config, galactic_center_dist_ly=1.6e5).metallicity_fe_h for _ in range(40)]
    assert np.mean(drawn) == pytest.approx(cloud.metallicity_fe_h, abs=0.2)
    assert bright_stars.star_fe_h({"metallicity_fe_h": -1.0}) == -1.0      # the star carries the birth [Fe/H] its component drew
    elsewhere = bright_stars.FillContext(profile.profile_at((SOLAR_RADIUS, 0.0, 0.0), shape, EDGE, E1))
    assert elsewhere.apply(SystemConfig(), draw.Stream("test-fill")).FE_H_PEAKS is None


@pytest.mark.fresh_database
def test_satellites_round_trip_through_the_store(mysql_config):
    found = _placed()
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, found, kind="satellite")
        conn.commit()
        back = store.get_galaxy_components(conn, kind="satellite")
        assert [c.name for c in back] == [c.name for c in found]
        exact = {c.name: c for c in back}
        for original in found:
            stored = exact[original.name]
            assert (stored.profile, stored.king_w0, stored.catalogue_id, stored.source) == (
                original.profile, None, original.name, "catalogue")
            assert stored.tidal_radius_pc == pytest.approx(original.tidal_radius_pc)
            assert stored.center_pc == pytest.approx(original.center_pc)
    finally:
        conn.rollback()
        conn.close()


def test_a_field_made_of_clusters_and_satellites_is_one_object_for_equal_rows():
    found = tuple(_placed())
    assert store._field_of(found, EDGE, E1) is store._field_of(tuple(found), EDGE, E1)
