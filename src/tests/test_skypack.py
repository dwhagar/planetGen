"""VIEW.29: the packed sky list (`planetgen/web/lib/skypack.py`). static/skypack.js reads the same bytes (tests/js/skypack.test.mjs, fed by
`test_js_unit.sky_pack`)."""

import math

import numpy as np
import pytest

from planetgen.galaxy import sky
from planetgen.web.lib import skypack
from tests.test_sky import OBSERVER, make

C_KMS = 299792.458
PC_KM = 3.0856775814913673e13


def a_sky(count=60, seed=5):
    rng = np.random.default_rng(seed)
    tiers = ["charted", "bright_star", "phenomenon", "cluster"]
    rows = []
    for i in range(count):
        direction = rng.normal(size=3)
        direction /= np.linalg.norm(direction)
        rows.append({"uid": f"s{i}", "tier": tiers[i % 4], "position": tuple(np.array(OBSERVER.position_pc) + direction * 10.0 ** rng.uniform(0.5, 2.7)),
                     "velocity": tuple(rng.normal(0.0, 40.0, 3)), "luminosity": 10.0 ** rng.uniform(0.0, 4.0), "temperature": rng.uniform(3000.0, 25000.0)})
    return make(rows, limit=9.0, dust={}), rows


def test_a_sky_comes_back_close_to_what_went_in():
    result, _ = a_sky()
    assert len(result.uid) > 40
    observer_velocity = (4.0, -7.0, 2.0)
    back = skypack.unpack_sky(skypack.pack_sky(result, 1.5e9, 9.0, observer_velocity, meta={"observer": "x"}))
    assert back["count"] == len(result.uid) and back["epoch_unix"] == 1.5e9 and back["limit"] == pytest.approx(9.0)
    assert back["observer_velocity_kms"] == observer_velocity and back["meta"] == {"observer": "x"}
    # float32 offsets: a part in 1e7; the relative velocity to the step of its code; the magnitude to 1/1024; the temperature to 1.6 percent
    assert np.allclose(back["offset_pc"], result.offset_pc, rtol=2e-7, atol=1e-9)
    assert np.abs(back["velocity_kms"] - (result.velocity_kms - np.array(observer_velocity))).max() <= back["velocity_scale"] / 2 + 1e-9
    assert np.abs(back["magnitude"] - result.dimmed_absolute_magnitude).max() <= 0.5 / skypack.MAGNITUDE_STEPS_PER_MAG + 1e-12
    assert np.allclose(back["temperature_k"], result.temperature_k, rtol=0.02)
    assert list(back["index"]) == list(range(len(result.uid)))
    names = {v: k for k, v in skypack.KIND_CODES.items()}
    assert [names[k] for k in back["kind"]] == list(result.tier)


def test_a_star_is_26_bytes_and_a_gzipped_list_is_about_half_that():
    import gzip
    result, _ = a_sky(2000)
    raw = skypack.pack_sky(result, 0.0, 9.0)
    assert skypack.BYTES_PER_STAR == 26
    n = len(result.uid)
    assert len(raw) == 52 + 2 + n * 26                                  # the header, the empty JSON and the columns
    assert len(gzip.compress(raw)) < 0.9 * len(raw)


def test_the_card_can_put_a_star_where_the_server_sees_it():
    """apparent = offset + relative velocity * t - velocity * tau, tau = |apparent| / c: from the packed values it gives the direction
    `sky_from_stars` computed."""
    result, _ = a_sky(30)
    observer_velocity = np.array([4.0, -7.0, 2.0])
    back = skypack.unpack_sky(skypack.pack_sky(result, 0.0, 9.0, observer_velocity))
    v = back["velocity_kms"] + observer_velocity                           # the star's own, km/s
    tau = np.zeros(len(v))                                                  # seconds; t = 0, the offsets hold at the sky's time
    for _ in range(4):
        apparent = back["offset_pc"] - v * tau[:, None] / PC_KM
        tau = np.linalg.norm(apparent, axis=1) * PC_KM / C_KMS
    seen = apparent / np.linalg.norm(apparent, axis=1)[:, None]
    assert np.allclose(seen, result.direction, atol=2e-5)
    magnitude = back["magnitude"] + 5.0 * np.log10(np.linalg.norm(apparent, axis=1) / 10.0)
    assert np.allclose(magnitude, result.magnitude, atol=2e-3)


def test_an_empty_sky_packs_and_reads():
    empty = sky.sky_from_stars(OBSERVER, sky.empty_stars(), 0.0)
    raw = skypack.pack_sky(empty, 5.0, 7.5)
    back = skypack.unpack_sky(raw)
    assert back["count"] == 0 and back["offset_pc"].shape == (0, 3) and skypack.is_whole(raw)


def test_a_torn_or_foreign_list_is_not_whole():
    result, _ = a_sky(10)
    raw = skypack.pack_sky(result, 0.0, 9.0)
    assert skypack.is_whole(raw)
    assert not skypack.is_whole(raw[:-1]) and not skypack.is_whole(raw + b"\0") and not skypack.is_whole(None)
    assert not skypack.is_whole(b"PGT1" + raw[4:])
    with pytest.raises(skypack.SkyPackError):
        skypack.unpack_sky(raw[:-5])


def test_a_hot_and_a_cool_star_get_their_temperatures_and_no_temperature_is_none():
    result, _ = a_sky(8)
    cool = result._replace(temperature_k=np.where(np.arange(len(result.uid)) == 0, np.nan, result.temperature_k))
    back = skypack.unpack_sky(skypack.pack_sky(cool, 0.0, 9.0))
    assert math.isnan(back["temperature_k"][0]) and not math.isnan(back["temperature_k"][1])


def test_what_cannot_be_packed_is_refused():
    result, _ = a_sky(4)
    with pytest.raises(skypack.SkyPackError):
        skypack.pack_sky(result._replace(tier=np.array(["moon"] * len(result.uid), dtype=object)), 0.0, 9.0)
    with pytest.raises(skypack.SkyPackError):
        skypack.pack_sky(result._replace(dimmed_absolute_magnitude=result.dimmed_absolute_magnitude + 100.0), 0.0, 9.0)
    with pytest.raises(skypack.SkyPackError):
        skypack.pack_sky(result, 0.0, 9.0, indices=[2 ** 32] * len(result.uid))
