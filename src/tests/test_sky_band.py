"""VIEW.27: the glow band of the unresolved stars (`planetgen/galaxy/sky_band.py`): the cube map, the cache key and cell, the
subtraction of the resolved stars, and the ray-march against an analytic density."""

import math

import numpy as np
import pytest

from planetgen.galaxy import sky_band as band
from planetgen.galaxy.density import build_galaxy_shape

NO_DUST = {"av_per_kpc_at_sun": 0.0}
SUN = (8128.0, 0.0, 0.0)


def uniform(level):
    """A density function with the same relative density everywhere."""
    return lambda x, y, z, shape: np.full(np.broadcast(x, y, z).shape, float(level))


def calls_counter():
    calls = []

    def density_at(x, y, z, shape):
        calls.append(1)
        return np.ones(np.broadcast(x, y, z).shape)

    return density_at, calls


@pytest.fixture(scope="module")
def milky_way_like():
    return build_galaxy_shape(
        disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
        arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)


# --- the cube map -----------------------------------------------------------------------------------------------------

def test_face_directions_are_unit_vectors_that_cover_the_sphere():
    d = band.face_directions(8)
    assert d.shape == (6, 8, 8, 3)
    assert np.allclose(np.linalg.norm(d, axis=-1), 1.0)
    flat = d.reshape(-1, 3)
    for axis in range(3):
        assert flat[:, axis].max() > 0.9 and flat[:, axis].min() < -0.9
    assert np.allclose(flat.mean(axis=0), 0.0, atol=1e-12)


def test_every_pixel_looks_itself_up():
    d = band.face_directions(8)
    index = np.arange(6 * 8 * 8).reshape(6, 8, 8)
    assert (band.face_lookup(index, d) == index).all()


def test_a_face_needs_a_pixel():
    with pytest.raises(ValueError):
        band.face_directions(0)


# --- the share of the light that is not drawn -------------------------------------------------------------------------

@pytest.mark.parametrize("cut, share", [(4.81, 0.07), (2.31, 0.35), (-0.19, 0.81)])
def test_the_share_of_faint_light_is_the_local_measure(cut, share):
    assert band.faint_share(cut) == pytest.approx(share)


def test_the_share_runs_from_none_to_all_and_never_falls_as_the_cut_gets_brighter():
    cuts = np.linspace(-15.0, 15.0, 301)
    shares = band.faint_share(cuts)
    assert shares[0] == 1.0 and shares[-1] == 0.0
    assert (np.diff(shares) <= 1e-12).all()


# --- the steps and the cell -------------------------------------------------------------------------------------------

def test_steps_cover_the_ray_and_are_fine_near_the_observer():
    s, ds = band.step_distances(400, 30000.0)
    assert ds.sum() == pytest.approx(30000.0)
    assert (np.diff(s) > 0).all() and s[0] > 0
    assert ds[0] < 1.0 < ds[-1]


def test_the_cell_holds_the_observer_and_its_center_is_in_the_middle():
    assert band.cell_of((0.0, 0.0, 0.0)) == (0, 0, 0)
    assert band.cell_of((49.9, 50.0, -0.1)) == (0, 1, -1)
    assert band.cell_center((0, 1, -1)) == (25.0, 75.0, -25.0)
    assert band.cell_of(band.cell_center((7, -3, 2))) == (7, -3, 2)


def test_the_key_is_the_same_in_a_cell_and_changes_with_everything_the_image_depends_on():
    base = band.cache_key((1, 2, 3), "galaxy-1")
    assert band.cache_key((1, 2, 3), "galaxy-1") == base
    assert len(base) == 32
    changes = [
        band.cache_key((1, 2, 4), "galaxy-1"), band.cache_key((1, 2, 3), "galaxy-2"),
        band.cache_key((1, 2, 3), "galaxy-1", limiting_magnitude=6.5), band.cache_key((1, 2, 3), "galaxy-1", pixels=16),
        band.cache_key((1, 2, 3), "galaxy-1", steps=200), band.cache_key((1, 2, 3), "galaxy-1", max_distance_pc=20000.0),
        band.cache_key((1, 2, 3), "galaxy-1", dust={"scale_length_pc": 3000.0}),
    ]
    assert len({base, *changes}) == len(changes) + 1


def test_the_key_follows_the_tables_of_a_dust_profile():
    from planetgen.physics import dust_model
    thin = dust_model.DustProfile((0.0, 5000.0, 10000.0), (0.2, 0.1, 0.0), (100.0, 150.0, 300.0))
    same = dust_model.DustProfile((0.0, 5000.0, 10000.0), (0.2, 0.1, 0.0), (100.0, 150.0, 300.0))
    other = dust_model.DustProfile((0.0, 5000.0, 10000.0), (0.3, 0.1, 0.0), (100.0, 150.0, 300.0))
    keys = [band.cache_key((1, 2, 3), "g", dust={"profile": p}) for p in (thin, same, other)]
    assert keys[0] == keys[1] != keys[2]
    assert keys[0] != band.cache_key((1, 2, 3), "g")


def test_the_key_follows_the_constants_of_the_image(monkeypatch):
    base = band.cache_key((0, 0, 0), "g")
    monkeypatch.setattr(band, "LOCAL_V_LUMINOSITY_DENSITY", 0.05)
    assert band.cache_key((0, 0, 0), "g") != base
    monkeypatch.undo()
    monkeypatch.setattr(band, "MODEL_VERSION", band.MODEL_VERSION + 1)
    assert band.cache_key((0, 0, 0), "g") != base
    monkeypatch.undo()
    monkeypatch.setattr(band, "FAINT_LIGHT_TABLE", ((-10.0, 1.0), (12.0, 0.0)))
    assert band.cache_key((0, 0, 0), "g") != base


# --- the ray-march against a density with a known answer ---------------------------------------------------------------

def test_uniform_density_without_dust_and_with_every_star_in_the_glow_is_the_analytic_flux():
    # f_faint = 1 (a limit so bright that nothing is drawn): J = j0 * rho * max_distance.
    smax = 5000.0
    b = band.surface_brightness(SUN, None, limiting_magnitude=-100.0, pixels=4, steps=200, max_distance_pc=smax,
                                dust=NO_DUST, density_at=uniform(2.0))
    flux = band.LOCAL_V_LUMINOSITY_DENSITY * 2.0 * smax
    expected = band.SB_ZERO_POINT - 2.5 * math.log10(band.OMEGA_ARCSEC2 * flux)
    assert b.shape == (6, 4, 4)
    assert np.allclose(b, expected, atol=1e-4)


def test_the_subtraction_leaves_exactly_the_faint_share_of_the_light():
    smax, m_lim = 20000.0, 7.5
    b = band.surface_brightness(SUN, None, limiting_magnitude=m_lim, pixels=2, steps=400, max_distance_pc=smax,
                                dust=NO_DUST, density_at=uniform(1.0))
    s = np.linspace(0.0, smax, 4_000_001)[1:]
    share = band.faint_share(m_lim - 5.0 * np.log10(s / 10.0))
    flux = band.LOCAL_V_LUMINOSITY_DENSITY * np.sum(share) * (s[1] - s[0])
    expected = band.SB_ZERO_POINT - 2.5 * math.log10(band.OMEGA_ARCSEC2 * flux)
    assert np.allclose(b, expected, atol=0.01)


def test_a_brighter_limit_leaves_more_light_in_the_glow_and_a_fainter_one_less():
    levels = [band.surface_brightness(SUN, None, limiting_magnitude=m, pixels=2, steps=300, dust=NO_DUST,
                                      density_at=uniform(1.0))[0, 0, 0] for m in (4.5, 6.0, 7.5, 9.0, 12.0)]
    assert all(a < b for a, b in zip(levels, levels[1:]))             # a lower surface brightness number is brighter


def test_a_limit_that_draws_everything_leaves_a_dark_sky():
    b = band.surface_brightness(SUN, None, limiting_magnitude=60.0, pixels=2, steps=200, dust=NO_DUST,
                                density_at=uniform(1.0))
    assert (b == band.DARK_MAG).all()


def test_an_empty_galaxy_is_dark():
    b = band.surface_brightness(SUN, None, pixels=2, steps=50, density_at=uniform(0.0))
    assert (b == band.DARK_MAG).all()


def test_dust_dims_the_band():
    clear = band.surface_brightness(SUN, None, pixels=2, steps=300, dust=NO_DUST, density_at=uniform(1.0))
    dusty = band.surface_brightness(SUN, None, pixels=2, steps=300, density_at=uniform(1.0))
    assert (dusty >= clear).all()
    plane = band.face_lookup(dusty, np.array([0.0, 1.0, 0.0]))
    pole = band.face_lookup(dusty, np.array([0.0, 0.0, 1.0]))
    assert plane - band.face_lookup(clear, np.array([0.0, 1.0, 0.0])) > pole - band.face_lookup(clear, np.array([0.0, 0.0, 1.0]))


# --- the band of a galaxy shaped like the Milky Way --------------------------------------------------------------------

def test_the_band_of_a_milky_way_like_galaxy_is_brightest_towards_the_center_and_in_the_plane(milky_way_like):
    b = band.surface_brightness(SUN, milky_way_like)
    center = band.face_lookup(b, np.array([-1.0, 0.0, 0.04]))
    anticenter = band.face_lookup(b, np.array([1.0, 0.0, 0.0]))
    pole = band.face_lookup(b, np.array([0.0, 0.0, 1.0]))
    assert center < anticenter < pole
    assert 21.0 < center < 23.5                 # the sky-view.md 3.1 figures for the Milky Way are 21.6 to 23.0
    assert 22.0 < anticenter < 24.5
    assert 22.5 < pole < 25.5
    assert 20.0 < b.min() and b.max() < 26.5


def test_moving_twenty_light_years_changes_the_band_by_a_few_hundredths_of_a_magnitude(milky_way_like):
    here = band.surface_brightness(SUN, milky_way_like, pixels=4, steps=300)
    moved = band.surface_brightness((SUN[0] + 6.13, SUN[1], SUN[2]), milky_way_like, pixels=4, steps=300)
    assert np.abs(here - moved).max() < 0.02


# --- the cache --------------------------------------------------------------------------------------------------------

def test_a_cell_is_worked_out_once_and_shared_by_the_systems_in_it():
    density_at, calls = calls_counter()
    store = band.MemoryBandStore()
    kwargs = dict(pixels=2, steps=20, density_at=density_at)
    first = band.band_for((10.0, 10.0, 10.0), None, "g", store, **kwargs)
    worked = len(calls)
    assert worked > 0
    second = band.band_for((40.0, 3.0, 49.0), None, "g", store, **kwargs)
    assert len(calls) == worked and np.array_equal(first, second)
    band.band_for((60.0, 10.0, 10.0), None, "g", store, **kwargs)                # the next cell
    assert len(calls) > worked
    after = len(calls)
    band.band_for((10.0, 10.0, 10.0), None, "other galaxy", store, **kwargs)     # same cell, another galaxy
    assert len(calls) > after


def test_the_band_of_a_cell_is_worked_from_its_center():
    store = band.MemoryBandStore()
    stored = band.band_for((10.0, 10.0, 10.0), None, "g", store, pixels=2, steps=50, density_at=uniform(1.0))
    direct = band.surface_brightness((25.0, 25.0, 25.0), None, pixels=2, steps=50, density_at=uniform(1.0))
    assert np.array_equal(stored, direct)


def test_files_keep_a_band_between_runs_and_a_bad_file_is_a_miss(tmp_path):
    store = band.FileBandStore(tmp_path / "bands")
    assert store.get("missing") is None
    image = np.arange(24, dtype=np.float32).reshape(6, 2, 2)
    store.put("k", image)
    assert np.array_equal(band.FileBandStore(tmp_path / "bands").get("k"), image)
    assert [p.name for p in (tmp_path / "bands").iterdir()] == ["k.npy"]
    (tmp_path / "bands" / "bad.npy").write_bytes(b"not an array")
    assert store.get("bad") is None


def test_a_file_store_serves_band_for_and_recomputes_a_damaged_file(tmp_path):
    density_at, calls = calls_counter()
    store = band.FileBandStore(tmp_path)
    first = band.band_for((5.0, 5.0, 5.0), None, "g", store, pixels=2, steps=20, density_at=density_at)
    worked = len(calls)
    again = band.band_for((6.0, 7.0, 8.0), None, "g", band.FileBandStore(tmp_path), pixels=2, steps=20,
                          density_at=density_at)
    assert len(calls) == worked and np.array_equal(first, again)
    for path in tmp_path.iterdir():
        path.write_bytes(b"torn")
    third = band.band_for((5.0, 5.0, 5.0), None, "g", store, pixels=2, steps=20, density_at=density_at)
    assert len(calls) > worked and np.array_equal(first, third)


# --- the glow image the page samples (VIEW.13, `skyband.js`) ----------------------------------------------------------------

def test_the_page_unit_is_ten_to_the_minus_point_four_of_the_magnitude_over_twenty_two_and_dark_is_zero():
    values = band.values_of(np.array([22.0, 24.5, 19.5, band.DARK_MAG]))
    assert values == pytest.approx([1.0, 0.1, 10.0, 0.0])
    assert band.SB_REF == 22.0 and (band.EQUIRECT_WIDTH, band.EQUIRECT_HEIGHT) == (256, 128)


def test_the_tint_has_a_luminance_of_one_and_is_warm():
    red, green, blue = band.tint()
    assert 0.2126 * red + 0.7152 * green + 0.0722 * blue == pytest.approx(1.0)
    assert red > green > blue > 0.0


def test_the_equirect_image_is_the_band_in_the_directions_the_page_looks_up():
    sb = np.full((6, 8, 8), 22.0)
    sb[4] = 20.0                                                     # the +z face (the north pole), brighter
    image = band.equirect(sb, 64, 32)
    assert image.shape == (32, 64, 3) and image.dtype == np.float32
    north, equator, south = image[0, :, 1], image[16, :, 1], image[31, :, 1]
    assert north.min() > equator.max()          # row 0 is the north pole
    assert south.mean() == pytest.approx(equator.mean(), rel=0.05)
    assert np.allclose(equator, equator[0], rtol=0.02)                                         # the equator is the same all the way round
    luminance = image @ np.array([0.2126, 0.7152, 0.0722])
    assert luminance[16].mean() == pytest.approx(1.0, rel=0.02)


def test_the_equirect_texel_at_a_longitude_looks_that_way():
    values = np.full((6, 8, 8), 30.0)
    values[0] = 20.0                                                 # the +x face only
    image = band.equirect(values, 64, 32)
    row = image[16, :, 1]
    assert row[31] == row[32] == row.max() and row[0] < row[31] / 100.0 and row[16] < row[31] / 100.0 and row[48] < row[31] / 100.0   # column 32 is longitude 0


def test_half_float_bytes_are_little_endian_and_six_for_each_texel():
    image = band.equirect(np.full((6, 8, 8), 22.0), 256, 128)
    raw = band.to_half_bytes(image)
    assert len(raw) == 256 * 128 * 6
    back = np.frombuffer(raw, dtype="<f2").reshape(128, 256, 3)
    assert np.allclose(back, image, rtol=2e-3)
