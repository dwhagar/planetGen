"""VIEW.43: the dust of a galaxy from its own gas (`planetgen/physics/dust_model.py`)."""

import collections
import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import disk_history
from planetgen.physics import dust_model as dm, sky_photometry as sky

MW_MASS, MW_LENGTH = 5.0e10, 2600.0
SMALL_MASS, SMALL_LENGTH = 1.6e8, 481.0
SUN_PC = 8100.0


def _profile(mass, length):
    thin = mass * 0.8
    history = disk_history.history(thin, length, disk_history.disk_speed_kms(mass))
    return dm.profile_from_history(history, length)


@pytest.fixture(scope="module")
def milky_way():
    return _profile(MW_MASS, MW_LENGTH)


@pytest.fixture(scope="module")
def small():
    return _profile(SMALL_MASS, SMALL_LENGTH)


# --- the parts -----------------------------------------------------------------------------------------------------------

def test_the_constants_are_the_ones_read_or_argued():
    assert dm.HYDROGEN_ATOMS_PER_AV_CM2 == 2.21e21                      # [R5123: Eq. 1]
    assert dm.hydrogen_mass_fraction() == pytest.approx(0.713, abs=0.005)
    assert dm._MSUN_PC2_IN_HYDROGEN_PER_CM2 == pytest.approx(1.248e20, rel=2e-3)
    assert dm._G_PC_KMS2_PER_MSUN == pytest.approx(4.30e-3, rel=1e-3)


def test_the_dust_goes_with_the_metals_and_is_the_solar_ratio_at_solar():
    assert dm.dust_to_gas_ratio(0.0) == pytest.approx(1.0)
    assert dm.dust_to_gas_ratio(-0.3) == pytest.approx(10 ** -0.3)      # linear in the metals above the break of the law
    assert dm.dust_to_gas_ratio(0.2) == pytest.approx(10 ** 0.2)
    many = dm.dust_to_gas_ratio(np.array([-0.3, 0.0, 0.2]))
    assert many.shape == (3,) and many[1] == pytest.approx(1.0)
    assert dm.dust_to_gas_ratio(-1.0) < 10 ** -1.0                       # under the break the law steepens: poorer in dust than the metals alone


def test_the_gas_is_solar_at_the_reference_radius_and_poorer_outside_it():
    assert dm.gas_oh_dex(8200.0, 8200.0) == 0.0
    assert dm.gas_oh_dex(18200.0, 8200.0) == pytest.approx(-0.34)
    assert dm.gas_oh_dex(3200.0, 8200.0) > 0.0


def test_the_layer_is_a_gas_sheet_in_its_own_gravity_and_no_thinner_than_the_floor():
    # sigma^2 / (pi G Sigma): 9 km/s in 53 solar masses per pc2 is 113 pc
    assert dm.layer_height_pc(53.1) == pytest.approx(113.0, rel=0.02)
    assert dm.layer_height_pc(1000.0) == dm.LAYER_FLOOR_PC
    assert dm.layer_height_pc(0.0) > 1.0e6                                # no gravity, no layer: the sheet is unbounded, never a division by zero
    assert dm.layer_height_pc(25.0, dispersion_kms=18.0) == pytest.approx(4.0 * dm.layer_height_pc(25.0), rel=1e-9)


def test_the_column_follows_the_gas_the_visible_share_and_the_metals():
    base = dm.column_av(11.8, 0.0)
    assert base == pytest.approx(0.24, rel=0.05)                         # the measured column through the layer at the Sun, [R5125]
    assert dm.column_av(23.6, 0.0) == pytest.approx(2.0 * base)
    assert dm.column_av(11.8, -0.3) == pytest.approx(0.5 * base, rel=0.01)
    assert dm.column_av(11.8, 0.0, visible_share=1.0) == pytest.approx(2.0 * base)
    assert dm.column_av(0.0, 0.0) == 0.0


# --- the Milky Way against what is measured -------------------------------------------------------------------------------

def test_the_milky_way_has_the_measured_column_and_layer_at_the_sun(milky_way):
    assert float(milky_way.column_at(SUN_PC)) == pytest.approx(0.24, abs=0.03)       # 0.12 through a half layer, [R5125]
    assert 100.0 < float(milky_way.height_at(SUN_PC)) < 150.0                          # 125 to 134 pc, a single sech^2 layer
    assert 0.7 < float(milky_way.plane_av_per_kpc(SUN_PC)) < 1.3                      # 0.85 +- 0.20, [R5125]; 0.8 to 1.0, [R5124]


def test_the_layer_flares_and_the_plane_dims_less_outward(milky_way):
    radii = np.array([6000.0, 8000.0, 10000.0, 12500.0, 14000.0])
    heights = milky_way.height_at(radii)
    assert np.all(np.diff(heights) > 0)
    assert 300.0 < float(milky_way.height_at(12500.0)) < 700.0                         # thin 340 pc and thick 1,000 pc at 12 to 13 kpc, [R5124]
    plane = milky_way.plane_av_per_kpc(radii)
    assert np.all(np.diff(plane) < 0)
    assert float(plane[3]) == pytest.approx(0.14, abs=0.1)                             # thin plus thick x 3.1 at 12 to 13 kpc, [R5124]


def test_the_plane_from_the_sun_to_the_center_dims_by_tens_of_magnitudes_not_two(milky_way):
    radii = np.linspace(50.0, SUN_PC, 4000)
    total = np.trapezoid(milky_way.av_per_pc(radii, 0.0, 0.0), radii)
    assert 25.0 < total < 60.0                                                         # A_V 31 to 59 in the literature, 38 to 42 for [R5127]
    old = np.trapezoid(1.0e-3 * np.exp(-(8000.0 - radii) / sky.DUST_SCALE_LENGTH_PC), radii)
    assert old < 3.0                                                                   # the fixed law gave 2.3


def test_beyond_the_table_there_is_no_dust(milky_way):
    assert float(milky_way.av_per_pc(milky_way.radius_pc[-1] * 1.5, 0.0, 0.0)) == 0.0
    assert float(milky_way.column_at(0.0)) > 0.0                                       # inside the first ring the first ring answers


# --- other galaxies --------------------------------------------------------------------------------------------------------

def test_a_small_galaxy_is_nearly_clear_where_the_fixed_law_gave_nine_and_a_half_magnitudes_a_kpc(small, milky_way):
    site = 3.15 * SMALL_LENGTH
    assert 0.05 < float(small.plane_av_per_kpc(site)) < 0.3
    assert float(small.plane_av_per_kpc(site)) < 0.3 * float(milky_way.plane_av_per_kpc(SUN_PC))
    assert float(small.height_at(site)) > 3.0 * float(milky_way.height_at(SUN_PC))     # a thicker layer of thinner gas
    old = math.exp((8000.0 - 2885.0) / sky.DUST_SCALE_LENGTH_PC)
    assert old > 9.0


def test_the_extinction_per_kpc_falls_with_the_disk(small, milky_way):
    middle = _profile(1.0e10, 1800.0)
    at = lambda p, length: float(p.plane_av_per_kpc(3.15 * length))                    # noqa: E731
    assert at(small, SMALL_LENGTH) < at(middle, 1800.0) < at(milky_way, MW_LENGTH)


# --- the layer integrates to the column ------------------------------------------------------------------------------------

def test_the_extinction_through_the_layer_is_the_column(milky_way):
    z = np.linspace(-5000.0, 5000.0, 200001)
    for radius in (3000.0, SUN_PC, 12000.0):
        through = np.trapezoid(milky_way.av_per_pc(radius, 0.0, z), z)
        assert through == pytest.approx(float(milky_way.column_at(radius)), rel=2e-3)


def test_height_follows_sech_squared(milky_way):
    h = float(milky_way.height_at(SUN_PC))
    mid = float(milky_way.av_per_pc(SUN_PC, 0.0, 0.0))
    assert float(milky_way.av_per_pc(SUN_PC, 0.0, h)) == pytest.approx(mid / math.cosh(1.0) ** 2)
    assert float(milky_way.av_per_pc(0.0, SUN_PC, 0.0)) == pytest.approx(mid)         # only the radius matters
    arrays = milky_way.av_per_pc(np.array([SUN_PC, 9000.0]), np.zeros(2), np.array([0.0, 50.0]))
    assert arrays.shape == (2,)


# --- the profile as a value ------------------------------------------------------------------------------------------------

def test_a_profile_is_a_value_with_a_key(milky_way):
    again = _profile(MW_MASS, MW_LENGTH)
    assert again == milky_way and hash(again) == hash(milky_way) and again.key() == milky_way.key()
    other = _profile(1.0e10, 1800.0)
    assert other != milky_way and other.key() != milky_way.key()
    assert len(milky_way.key()) == 16


@pytest.mark.parametrize("radius, column, height", [
    ((1.0,), (1.0,), (1.0,)),
    ((1.0, 1.0), (1.0, 1.0), (1.0, 1.0)),
    ((1.0, 2.0), (1.0, -1.0), (1.0, 1.0)),
    ((1.0, 2.0), (1.0, 1.0), (1.0, 0.0)),
    ((1.0, 2.0), (1.0,), (1.0, 1.0)),
])
def test_a_table_that_is_not_two_increasing_radii_with_a_positive_height_and_a_column_raises(radius, column, height):
    with pytest.raises(ValueError):
        dm.DustProfile(radius, column, height)


# --- how the sky reads it ---------------------------------------------------------------------------------------------------

def test_the_sky_law_takes_the_profile_and_no_dust_stays_no_dust(milky_way):
    assert float(sky.dust_av_per_pc(SUN_PC, 0.0, 0.0, profile=milky_way)) == pytest.approx(float(milky_way.av_per_pc(SUN_PC, 0.0, 0.0)))
    assert float(sky.dust_av_per_pc(SUN_PC, 0.0, 0.0, profile=milky_way, scale_length_pc=1.0)) == pytest.approx(float(milky_way.av_per_pc(SUN_PC, 0.0, 0.0)))
    zero = sky.dust_av_per_pc(np.array([SUN_PC, 9000.0]), np.zeros(2), np.zeros(2), profile=milky_way, av_per_kpc_at_sun=0)
    assert zero.shape == (2,) and not zero.any()
    assert float(sky.dust_av_per_pc(8000.0, 0.0, 0.0)) == pytest.approx(1.0e-3)       # without a profile it is the fixed law


def test_a_line_of_sight_in_the_plane_near_the_sun_dims_by_about_a_magnitude_a_kpc(milky_way):
    a_v = sky.galactic_extinction_av((SUN_PC, 0.0, 0.0), (0.0, 1.0, 0.0), 1000.0, profile=milky_way)
    assert 0.7 < a_v < 1.4
    up = sky.galactic_extinction_av((SUN_PC, 0.0, 0.0), (0.0, 0.0, 1.0), 5000.0, profile=milky_way)
    assert up == pytest.approx(0.5 * float(milky_way.column_at(SUN_PC)), rel=0.02)    # a half layer straight up
    free = sky.galactic_extinction_av((SUN_PC, 0.0, 0.0), (0.0, 1.0, 0.0), 1000.0, profile=milky_way, av_per_kpc_at_sun=0)
    assert free == 0.0


# --- the shape ---------------------------------------------------------------------------------------------------------------

def test_a_profile_for_a_shape_reads_the_shapes_gas_and_scale_length(monkeypatch):
    Shape = collections.namedtuple("Shape", "disk_scale_length_pc k_norm")
    seen = []

    def history_of(shape):
        seen.append(shape)
        return disk_history.history(MW_MASS * 0.8, shape.disk_scale_length_pc, disk_history.disk_speed_kms(MW_MASS))

    monkeypatch.setattr(disk_history, "history_of", history_of)
    dm._profile_for_shape.cache_clear()
    one = dm.profile_for_shape(Shape(MW_LENGTH, 3.0))
    again = dm.profile_for_shape(Shape(MW_LENGTH, 7.0))                                # the normalization does not change the gas
    assert one == again and len(seen) == 1
    assert one == _profile(MW_MASS, MW_LENGTH)
    assert one.radius_pc[0] < tuning.REFERENCE_RADIUS_SCALE_LENGTHS * MW_LENGTH < one.radius_pc[-1]
    dm._profile_for_shape.cache_clear()
