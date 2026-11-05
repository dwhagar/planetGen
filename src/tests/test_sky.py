"""VIEW.23: the sky catalogue of a system from stored stars (`planetgen/galaxy/sky.py`): direction, apparent magnitude, colour, light-time
position, the stand-in velocity, the own-system exclusion and the report."""

import math

import numpy as np
import pytest

from planetgen.galaxy import sky
from planetgen.physics import photometry, sky_photometry
from planetgen.physics.constants import SPEED_OF_LIGHT_M_S

NO_DUST = {"av_per_kpc_at_sun": 0.0}
PC_M = 3.0856775814913673e16
YEAR_S = 365.25 * 86400.0
OBSERVER = sky.Observer((8128.0, 0.0, 0.0), 0.0, (0.0, 0.0, 0.0), -1)


def flat_speed(radius):
    return np.full(np.shape(radius), 220.0)


def stars(rows):
    """`SkyStars` from dicts; each takes its defaults from a Sun at rest 10 pc along +x."""
    default = dict(uid="s", tier="charted", system_id=-1, position=(8138.0, 0.0, 0.0), epoch=0.0, velocity=(0.0, 0.0, 0.0), luminosity=1.0,
                   temperature=5772.0, age=np.nan, phase_end=np.nan)
    full = [{**default, **row} for row in rows]
    return sky.SkyStars(
        np.array([r["uid"] for r in full], dtype=object), np.array([r["tier"] for r in full], dtype=object),
        np.array([r["system_id"] for r in full], dtype=int), np.array([r["position"] for r in full], dtype=float).reshape(-1, 3),
        np.array([r["epoch"] for r in full], dtype=float), np.array([r["velocity"] for r in full], dtype=float).reshape(-1, 3),
        np.array([r["luminosity"] for r in full], dtype=float), np.array([r["temperature"] for r in full], dtype=float),
        np.array([r["age"] for r in full], dtype=float), np.array([r["phase_end"] for r in full], dtype=float))


def make(rows, **kwargs):
    kwargs.setdefault("dust", NO_DUST)
    kwargs.setdefault("circular_speed_kms", flat_speed)
    return sky.sky_from_stars(kwargs.pop("observer", OBSERVER), stars(rows), kwargs.pop("at_unix", 0.0), **kwargs)


# --- where and how bright --------------------------------------------------------------------------------------------

def test_a_sun_at_rest_ten_parsecs_away_is_where_it_is_and_magnitude_4_82():
    result = make([{}])
    assert len(result.uid) == 1
    assert result.direction[0] == pytest.approx((1.0, 0.0, 0.0), abs=1e-9)
    assert result.magnitude[0] == pytest.approx(4.82, abs=0.02)
    assert result.distance_pc[0] == pytest.approx(10.0, rel=1e-9)
    assert result.light_time_years[0] == pytest.approx(10.0 * PC_M / SPEED_OF_LIGHT_M_S / YEAR_S, rel=1e-9)


def test_a_moving_star_is_seen_where_it_was_when_its_light_left():
    # 100 pc along +x, moving at 300 km/s along +y: the light left tau ago, so the star is seen back along its track.
    result = make([{"position": (8228.0, 0.0, 0.0), "velocity": (0.0, 300.0, 0.0), "luminosity": 100.0}], at_unix=1.0e9)
    d, tau = result.distance_pc[0] * PC_M, result.light_time_years[0] * YEAR_S
    moved = np.array([100.0 * PC_M, 300e3 * 1.0e9, 0.0])
    apparent = moved - np.array([0.0, 300e3, 0.0]) * tau
    assert d / SPEED_OF_LIGHT_M_S == pytest.approx(tau, rel=1e-9)
    assert np.allclose(result.direction[0], apparent / np.linalg.norm(apparent), atol=1e-9)
    assert result.distance_pc[0] * PC_M == pytest.approx(np.linalg.norm(apparent), rel=1e-9)


def test_an_observer_that_moves_sees_the_star_where_the_observer_has_got_to():
    # An observer that moves +y at 20 km/s for 1e9 s sees a star at rest 10 pc along +x displaced by -20e3 * 1e9 m along y.
    moving = sky.Observer(OBSERVER.position_pc, 0.0, (0.0, 20.0, 0.0), -1)
    result = make([{}], observer=moving, at_unix=1.0e9)
    offset = np.array([10.0 * PC_M, -20e3 * 1.0e9, 0.0])
    assert np.allclose(result.direction[0], offset / np.linalg.norm(offset), atol=1e-9)


def test_stars_with_other_epochs_are_brought_to_the_skys_time():
    # The same star stored at two epochs gives the same sky.
    one = make([{"position": (8138.0, 0.0, 0.0), "velocity": (0.0, 50.0, 0.0), "epoch": 0.0}], at_unix=2.0e9)
    p = 50e3 * 1.0e9 / PC_M
    two = make([{"position": (8138.0, p, 0.0), "velocity": (0.0, 50.0, 0.0), "epoch": 1.0e9}], at_unix=2.0e9)
    assert np.allclose(one.direction, two.direction, atol=1e-9)
    assert one.magnitude[0] == pytest.approx(two.magnitude[0], abs=1e-6)


def test_dust_dims_a_star_in_the_plane():
    far = {"position": (8128.0 + 2000.0, 0.0, 0.0), "luminosity": 1000.0}
    clear = make([far], limit=30.0)
    dusty = make([far], limit=30.0, dust={})
    assert dusty.magnitude[0] - clear.magnitude[0] == pytest.approx(
        sky_photometry.galactic_extinction_av(OBSERVER.position_pc, (1.0, 0.0, 0.0), 2000.0), abs=1e-6)
    assert dusty.magnitude[0] > clear.magnitude[0] + 1.0


# --- who is kept -----------------------------------------------------------------------------------------------------

def test_stars_fainter_than_the_limit_are_left_out_and_the_brightest_come_first():
    result = make([{"uid": "a"}, {"uid": "b", "luminosity": 1000.0}, {"uid": "c", "luminosity": 0.001}], limit=7.5)
    assert list(result.uid) == ["b", "a"]
    assert result.magnitude[0] < result.magnitude[1] <= 7.5
    assert make([{"uid": "c", "luminosity": 0.001}], limit=7.5).report["visible"] == 0
    assert len(make([{"uid": "c", "luminosity": 0.001}], limit=14.0).uid) == 1


def test_the_observers_own_system_is_left_out_with_its_companions():
    observer = sky.Observer(OBSERVER.position_pc, 0.0, (0.0, 0.0, 0.0), 7)
    result = make([{"uid": "own", "system_id": 7}, {"uid": "companion", "system_id": 7, "position": (8128.0, 0.0, 0.001)},
                   {"uid": "other", "system_id": 8}], observer=observer)
    assert list(result.uid) == ["other"]
    assert result.report["excluded_own_system"] == 2


def test_a_star_with_no_luminosity_or_temperature_is_counted_and_left_out():
    result = make([{"uid": "a", "luminosity": np.nan}, {"uid": "b", "temperature": np.nan}, {"uid": "c"}])
    assert list(result.uid) == ["c"] and result.report["no_photometry"] == 2


def test_a_star_that_had_not_formed_when_its_light_left_is_left_out():
    tau_years = 10.0 * PC_M / SPEED_OF_LIGHT_M_S / YEAR_S
    result = make([{"uid": "young", "age": tau_years / 2.0}, {"uid": "old", "age": tau_years * 2.0}, {"uid": "unknown"}])
    assert sorted(result.uid) == ["old", "unknown"]
    assert result.report["not_formed_yet"] == 1


def test_a_star_whose_phase_ended_within_the_light_time_is_counted():
    tau_years = 10.0 * PC_M / SPEED_OF_LIGHT_M_S / YEAR_S
    result = make([{"uid": "a", "age": 1.0e6, "phase_end": 1.0e6 - tau_years / 2.0}, {"uid": "b", "age": 1.0e6, "phase_end": 2.0e6}])
    assert len(result.uid) == 2 and result.report["other_phase"] == 1


def test_an_empty_catalogue_is_an_empty_sky():
    result = sky.sky_from_stars(OBSERVER, sky.empty_stars(), 0.0)
    assert len(result.uid) == 0 and result.report["visible"] == 0 and result.direction.shape == (0, 3)


# --- the report ------------------------------------------------------------------------------------------------------

def test_the_report_says_which_tiers_were_reached():
    result = make([{"uid": "a", "tier": "charted"}, {"uid": "b", "tier": "bright_star", "position": (8328.0, 0.0, 0.0), "luminosity": 1.0e5},
                   {"uid": "c", "tier": "cluster", "luminosity": 1e-6}])
    tiers = result.report["tiers"]
    assert set(tiers) == {"charted", "bright_star", "cluster"}
    assert tiers["charted"]["visible"] == 1 and tiers["bright_star"]["visible"] == 1 and tiers["cluster"]["visible"] == 0
    assert tiers["bright_star"]["farthest_pc"] == pytest.approx(200.0, rel=1e-6)
    assert tiers["charted"]["faintest_magnitude"] == pytest.approx(result.magnitude[result.tier == "charted"][0])
    assert tiers["cluster"]["faintest_magnitude"] is None
    assert result.report["stars_in"] == 3 and result.report["visible"] == 2 and result.report["limit"] == 7.5


# --- the stand-in velocity -------------------------------------------------------------------------------------------

def test_an_unbuilt_star_gets_the_same_velocity_every_time_and_another_star_another():
    p = np.array([[8128.0, 0.0, 0.0]] * 2)
    a = sky.unbuilt_velocity_kms(p, ["x", "y"], "g", flat_speed)
    b = sky.unbuilt_velocity_kms(p, ["x", "y"], "g", flat_speed)
    assert np.array_equal(a, b) and not np.allclose(a[0], a[1])
    assert not np.allclose(sky.unbuilt_velocity_kms(p[:1], ["x"], "other galaxy", flat_speed), a[:1])


def test_unbuilt_velocities_go_round_the_galaxy_counterclockwise_at_the_circular_speed():
    n = 3000
    uid = [f"u{i}" for i in range(n)]
    for position, tangent in [((8000.0, 0.0, 0.0), (0.0, 1.0, 0.0)), ((0.0, 5000.0, 0.0), (-1.0, 0.0, 0.0)), ((-3000.0, 0.0, 0.0), (0.0, -1.0, 0.0))]:
        v = sky.unbuilt_velocity_kms(np.tile(position, (n, 1)), uid, "g", flat_speed)
        assert v.mean(axis=0) @ np.array(tangent) == pytest.approx(220.0, abs=1.5)
        assert np.abs(v.mean(axis=0) - 220.0 * np.array(tangent)).max() < 1.5
    v = sky.unbuilt_velocity_kms(np.tile((8000.0, 0.0, 0.0), (n, 1)), uid, "g", flat_speed)
    assert v[:, 0].std() == pytest.approx(35.0, rel=0.06)
    assert v[:, 1].std() == pytest.approx(25.0, rel=0.06)
    assert v[:, 2].std() == pytest.approx(20.0, rel=0.06)


def test_a_star_with_no_stored_velocity_uses_the_stand_in_and_the_report_counts_it():
    row = {"velocity": (np.nan, np.nan, np.nan), "position": (8128.0, 100.0, 0.0), "uid": "scatter-1"}
    result = make([row])
    assert result.report["unbuilt_velocity"] == 1
    assert make([row]).direction.tolist() == result.direction.tolist()


def test_the_sky_of_a_star_on_the_galactic_axis_does_not_divide_by_zero():
    v = sky.unbuilt_velocity_kms(np.array([[0.0, 0.0, 10.0]]), ["a"], "g", flat_speed)
    assert np.all(np.isfinite(v))


# --- colour ----------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("temperature", [2500.0, 3500.0, 5772.0, 9000.0, 20000.0, 45000.0])
def test_the_colour_table_agrees_with_the_direct_colour(temperature):
    assert np.allclose(sky.star_colours(np.array([temperature]))[0], photometry.star_colour(temperature).linear, atol=0.01)


def test_hot_stars_are_blue_and_cool_ones_red():
    hot, cool = sky.star_colours(np.array([30000.0, 3000.0]))
    assert hot[2] > hot[0] and cool[0] > cool[2]
    assert max(hot) == pytest.approx(1.0, abs=0.01) and max(cool) == pytest.approx(1.0, abs=0.01)


# --- the shells and the source ---------------------------------------------------------------------------------------

def test_the_floor_is_the_luminosity_that_just_reaches_the_limit_at_that_distance():
    for distance in (1.0, 10.0, 300.0, 5000.0):
        floor = float(sky.luminosity_floor_solar(distance, 7.5))
        m = photometry.SUN_ABSOLUTE_BOLOMETRIC_MAG - 2.5 * math.log10(floor) + 5.0 * math.log10(distance / 10.0)
        assert m == pytest.approx(7.5, abs=1e-9)
    assert sky.luminosity_floor_solar(100.0, 7.5) > sky.luminosity_floor_solar(10.0, 7.5)
    assert sky.luminosity_floor_solar(100.0, 6.5) > sky.luminosity_floor_solar(100.0, 7.5)


def test_shells_grow_to_the_maximum():
    radii = sky.shell_radii_pc(40000.0)
    assert radii[0] == sky.FIRST_SHELL_PC and radii[-1] == 40000.0
    assert all(b > a for a, b in zip(radii, radii[1:]))


def test_asking_by_shells_finds_every_star_the_whole_catalogue_shows():
    rng = np.random.default_rng(3)
    n = 4000
    directions = rng.normal(size=(n, 3))
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    distance = 10.0 ** rng.uniform(0.0, 3.3, n)
    luminosity = 10.0 ** rng.uniform(-3.0, 4.0, n)
    rows = [{"uid": f"s{i}", "position": tuple(np.array(OBSERVER.position_pc) + directions[i] * distance[i]), "luminosity": float(luminosity[i]),
             "temperature": float(rng.uniform(3000.0, 12000.0))} for i in range(n)]
    everything = stars(rows)
    calls = []

    def source(centre, inner, outer, floor):
        calls.append((inner, outer, floor))
        d = np.linalg.norm(everything.position_pc - np.array(centre), axis=1)
        keep = (d >= inner) & (d < outer) & (everything.luminosity_sol >= floor)
        return sky.SkyStars(*[getattr(everything, name)[keep] for name in sky.SkyStars._fields])

    asked = sky.sky_from(OBSERVER, source, 0.0, 7.5, circular_speed_kms=flat_speed, dust=NO_DUST, max_distance_pc=4000.0)
    whole = sky.sky_from_stars(OBSERVER, everything, 0.0, 7.5, circular_speed_kms=flat_speed, dust=NO_DUST)
    assert len(whole.uid) > 50
    assert sorted(asked.uid) == sorted(whole.uid)
    assert asked.report["shells"] == sky.shell_radii_pc(4000.0)
    assert [c[0] for c in calls] == [0.0] + sky.shell_radii_pc(4000.0)[:-1]
    assert all(b[2] >= a[2] for a, b in zip(calls, calls[1:]))


# --- the sphere query of VIEW.22 ---------------------------------------------------------------------------------------

def fake_sphere(rows):
    """What `db.sky_sphere.SphereStars` carries, by attribute name, from `(source, ident, system, cluster, x, y, z, ...)` dicts."""
    from types import SimpleNamespace
    default = dict(source=0, ident=1, system_id=1, cluster_id=0, x=8138.0, y=0.0, z=0.0, v=(0.0, 0.0, 0.0), epoch=0.0, lum=1.0, temp=5772.0, age=np.nan,
                   life=np.inf)
    full = [{**default, **r} for r in rows]
    col = lambda name, dtype=float: np.array([r[name] for r in full], dtype=dtype)
    return SimpleNamespace(
        source=col("source", np.int8), ident=col("ident", np.int64), system_id=col("system_id", np.int64), cluster_id=col("cluster_id", np.int64),
        x_pc=col("x"), y_pc=col("y"), z_pc=col("z"), distance_pc=np.array([np.linalg.norm([r["x"] - 8128.0, r["y"], r["z"]]) for r in full]),
        vx_kms=np.array([r["v"][0] for r in full]), vy_kms=np.array([r["v"][1] for r in full]), vz_kms=np.array([r["v"][2] for r in full]),
        epoch_unix=col("epoch"), luminosity_sol=col("lum"), temperature_k=col("temp"), age_gy=col("age"), lifespan_gy=col("life"))


def test_a_sphere_becomes_stars_with_tiers_uids_and_years():
    sphere = fake_sphere([{"source": 0, "ident": 5, "system_id": 9, "age": 2.0, "life": 10.0}, {"source": 1, "ident": 6, "cluster_id": 3},
                          {"source": 1, "ident": 7}, {"source": 2, "ident": 8, "system_id": 0, "temp": np.nan}])
    converted = sky.stars_from_sphere(sphere)
    assert list(converted.tier) == ["charted", "cluster", "bright_star", "phenomenon"]
    assert list(converted.uid) == ["0:5", "1:6", "1:7", "2:8"]
    assert list(converted.system_id) == [9, -1, -1, -1]
    assert converted.age_years[0] == 2.0e9 and converted.phase_end_age_years[0] == 1.0e10 and np.isinf(converted.phase_end_age_years[1])
    assert converted.position_pc.shape == (4, 3) and converted.velocity_kms.shape == (4, 3)


def test_a_star_with_no_epoch_holds_at_the_skys_time():
    nan = float("nan")
    result = make([{"uid": "bright", "velocity": (nan, nan, nan), "epoch": nan}], at_unix=5.0e9)
    assert result.report["unknown_epoch"] == 1 and result.report["unbuilt_velocity"] == 1
    assert result.distance_pc[0] == pytest.approx(10.0, rel=0.01)


def test_the_sphere_source_gives_a_shell_and_the_sky_comes_out_of_it():
    sphere = fake_sphere([{"ident": 1, "x": 8128.0 + 5.0}, {"ident": 2, "x": 8128.0 + 30.0, "lum": 50.0},
                          {"ident": 3, "x": 8128.0 + 400.0, "lum": 2000.0}])
    asked = []

    def query(centre, radius, floor):
        asked.append((radius, floor))
        keep = (sphere.distance_pc <= radius) & (sphere.luminosity_sol >= floor)
        from types import SimpleNamespace
        return SimpleNamespace(**{k: v[keep] for k, v in vars(sphere).items()})

    source = sky.sphere_source(query)
    shell = source(OBSERVER.position_pc, 8.0, 16.0, 0.0)
    assert list(shell.uid) == []
    shell = source(OBSERVER.position_pc, 16.0, 32.0, 0.0)
    assert list(shell.uid) == ["0:2"]
    result = sky.sky_from(OBSERVER, source, 0.0, 7.5, circular_speed_kms=flat_speed, dust=NO_DUST, max_distance_pc=500.0)
    assert sorted(result.uid) == ["0:1", "0:2", "0:3"]
    assert result.report["tiers"]["charted"]["visible"] == 3


# --- what the packed list needs (VIEW.29) -----------------------------------------------------------------------------

def test_the_sky_carries_the_unretarded_offset_the_velocity_and_the_dimmed_absolute_magnitude():
    row = {"position": (8228.0, 0.0, 0.0), "velocity": (0.0, 300.0, 0.0), "luminosity": 100.0, "temperature": 9000.0}
    result = make([row], at_unix=1.0e9, dust={})
    advanced = np.array([100.0 * PC_M, 300e3 * 1.0e9, 0.0]) / PC_M
    assert np.allclose(result.offset_pc[0], advanced, rtol=1e-12)
    assert np.allclose(result.velocity_kms[0], (0.0, 300.0, 0.0)) and result.temperature_k[0] == 9000.0
    d = result.distance_pc[0]
    assert result.dimmed_absolute_magnitude[0] + 5.0 * math.log10(d / 10.0) == pytest.approx(result.magnitude[0], abs=1e-12)
    clear = make([row], at_unix=1.0e9)
    assert result.dimmed_absolute_magnitude[0] > clear.dimmed_absolute_magnitude[0]        # the dust of the line of sight is in it
    empty = sky.sky_from_stars(OBSERVER, sky.empty_stars(), 0.0)
    assert empty.offset_pc.shape == (0, 3) and empty.dimmed_absolute_magnitude.shape == (0,)


# --- the database entry point (VIEW.23 with VIEW.22) -------------------------------------------------------------------

def test_built_to_pc_is_the_last_shell_complete_down_to_its_floor_without_a_gap():
    inf = float("inf")
    assert sky.built_to_pc([(0.0, 8.0, 1e-5, 1e-6), (8.0, 16.0, 2e-5, 2e-5), (16.0, 32.0, 8e-5, 1e-3), (32.0, 64.0, 3e-4, 1e-9)]) == 16.0
    assert sky.built_to_pc([(0.0, 8.0, 1e-5, 1e-3)]) == 0.0
    assert sky.built_to_pc([(0.0, 8.0, 1e-5, inf)]) == 0.0
    assert sky.built_to_pc([]) == 0.0


def test_the_sphere_source_reports_the_levels_of_each_shell():
    from types import SimpleNamespace
    sphere = fake_sphere([{"ident": 1, "x": 8128.0 + 5.0}])
    sphere.levels = SimpleNamespace(complete_above_sol=0.5)
    levels = []
    source = sky.sphere_source(lambda centre, radius, floor: sphere, levels)
    source(OBSERVER.position_pc, 0.0, 8.0, 0.001)
    bare = fake_sphere([{"ident": 2}])                    # a query with no levels is complete above nothing it can say
    sky.sphere_source(lambda c, r, f: bare, levels)(OBSERVER.position_pc, 8.0, 16.0, 0.01)
    assert levels == [(0.0, 8.0, 0.001, 0.5), (8.0, 16.0, 0.01, float("inf"))]


def test_the_database_entry_point_asks_the_sphere_query_on_its_connection(monkeypatch):
    from types import SimpleNamespace
    sphere = fake_sphere([{"ident": 1, "x": 8128.0 + 5.0}, {"ident": 2, "x": 8128.0 + 30.0, "lum": 50.0}])
    sphere.levels = SimpleNamespace(complete_above_sol=0.0)
    seen = []

    def stars_in_sphere(conn, centre, radius, floor):
        seen.append((conn, radius, floor))
        keep = (sphere.distance_pc <= radius) & (sphere.luminosity_sol >= floor)
        picked = SimpleNamespace(**{k: v[keep] for k, v in vars(sphere).items() if k != "levels"})
        picked.levels = sphere.levels
        return picked

    monkeypatch.setattr(sky.sky_sphere, "stars_in_sphere", stars_in_sphere)
    result = sky.sky_from_database("the connection", OBSERVER, 0.0, 7.5, circular_speed_kms=flat_speed, dust=NO_DUST, max_distance_pc=500.0)
    assert sorted(result.uid) == ["0:1", "0:2"]
    assert all(conn == "the connection" for conn, *_ in seen) and len(seen) == len(sky.shell_radii_pc(500.0))
    assert result.report["built_to_pc"] == 500.0 and result.report["complete_above_sol"] == 0.0
