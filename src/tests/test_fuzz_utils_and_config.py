# tests/test_fuzz_utils_and_config.py

"""
Property-based / brute-force tests for the package's plumbing and pure
math: `stellarObjects/utils.py` (unit conversions, formatters, orbital
helpers, samplers), `appconfig.py` (hostile `config.json` contents),
`config.py` (`SystemConfig` round trips), `serialization.py`,
`progressFile.py`, `log.py`'s credential redaction, and sanity of every
constant in `physical_constants.py`/`program_constants.py`.

Name-generation helpers from `utils.py` live in `test_fuzz_names.py`.
See `fuzz_support.py` for the `ci`/`deep` profiles.
"""

import argparse
import json
import logging
import math
import os
import re
import sys
import types
from unittest import mock

import pytest
from hypothesis import assume, example, given, settings
from hypothesis import strategies as st

from stellarObjects import (
    _db, appconfig, log, physical_constants, program_constants, progressFile, serialization, utils,
)
from stellarObjects.config import SERIALIZABLE_FIELDS, SystemConfig
from tests.fuzz_support import any_float, finite, hostile_text, non_finite, scaled

positive = st.floats(min_value=1e-9, max_value=1e12, allow_nan=False, allow_infinity=False)
angle = st.floats(min_value=-1e4, max_value=1e4, allow_nan=False, allow_infinity=False)
normal_float = finite.filter(lambda x: x == 0 or abs(x) >= 1e-300)


def _close(a, b, rel=1e-9, abs_=1e-12):
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_)


def _config(markdown):
    config = SystemConfig()
    config.MARKDOWN = markdown
    return config


# ---------------------------------------------------------------------------
# utils: unit conversions
# ---------------------------------------------------------------------------

_ROUND_TRIPS = [
    (utils.ly_to_milliparsecs, utils.milliparsecs_to_ly),
    (utils.mpc_to_pc, utils.pc_to_mpc),
    (utils.pc_to_ly, utils.ly_to_pc),
    (utils.ly_to_au, utils.au_to_ly),
]


@given(x=normal_float)  # (a subnormal like 5e-324 legitimately underflows to 0 under /1000)
def test_unit_conversions_round_trip_and_preserve_sign(x):
    for forward, back in _ROUND_TRIPS:
        y = forward(x)
        assert math.isfinite(y)
        assert _close(back(y), x, rel=1e-12)
        assert (y > 0) == (x > 0) and (y == 0) == (x == 0)


@given(x=non_finite)
def test_unit_conversions_propagate_non_finite_without_raising(x):
    for forward, back in _ROUND_TRIPS:
        y = forward(x)
        assert math.isnan(y) if math.isnan(x) else y == x
        back(y)


def test_unit_conversions_agree_with_each_other():
    # ly -> mpc -> pc must equal ly -> pc directly.
    for ly in (1e-6, 1.0, 3.26156, 26_000.0, 1e9):
        assert _close(utils.mpc_to_pc(utils.ly_to_milliparsecs(ly)), utils.ly_to_pc(ly), rel=1e-9)
    assert _close(utils.pc_to_ly(1.0), 3.2616, rel=1e-3)


# ---------------------------------------------------------------------------
# utils: number formatting
# ---------------------------------------------------------------------------

_SCI_MD = re.compile(r"^(-?\d+(?:\.\d+)?) × 10<sup>(-?\d+)</sup>$")
_SCI_WIKI = re.compile(r"^\{\{Exp\|(-?\d+(?:\.\d+)?)\|(-?\d+)\}\}$")


@given(x=normal_float, precision=st.integers(0, 8), markdown=st.booleans())
def test_to_scientific_notation_parses_back_to_the_value(x, precision, markdown):
    out = utils.to_scientific_notation(_config(markdown), x, precision)
    if x == 0:
        assert out == "0"
        return
    match = (_SCI_MD if markdown else _SCI_WIKI).match(out)
    assert match, out
    coefficient, exponent = float(match.group(1)), int(match.group(2))
    assert exponent in (math.floor(math.log10(abs(x))), math.floor(math.log10(abs(x))) + 1)
    assert 1 <= abs(coefficient) < 10, out
    assert math.isclose(coefficient * 10.0 ** exponent, x, rel_tol=0.51 * 10 ** -precision + 1e-12)


@given(x=finite.filter(lambda v: v != 0), precision=st.integers(0, 6), markdown=st.booleans())
@example(x=9.999, precision=2, markdown=True)          # was "10.00 x 10^0"
@example(x=99.996, precision=2, markdown=True)
@example(x=0.0099999, precision=2, markdown=False)
@example(x=-9.9999e30, precision=2, markdown=True)
@example(x=1e-320, precision=2, markdown=True)          # subnormal, was "10.02 x 10^-321"
def test_to_scientific_notation_coefficient_is_normalized(x, precision, markdown):
    out = utils.to_scientific_notation(_config(markdown), x, precision)
    coefficient = float((_SCI_MD if markdown else _SCI_WIKI).match(out).group(1))
    assert 1 <= abs(coefficient) < 10, out


@given(x=any_float, markdown=st.booleans())
@example(x=math.nan, markdown=True)
@example(x=math.inf, markdown=True)
@example(x=-math.inf, markdown=False)
@example(x=5e-324, markdown=True)   # was ZeroDivisionError
def test_to_scientific_notation_never_raises(x, markdown):
    assert isinstance(utils.to_scientific_notation(_config(markdown), x), str)


@given(age=finite, precision=st.integers(0, 6))
def test_format_age_string_picks_the_right_unit(age, precision):
    out = utils.format_age_string(age, precision)
    assert out.endswith("Billion Years" if age >= 1.0 else "Million Years")


@given(age=non_finite)
def test_format_age_string_non_finite_does_not_raise(age):
    assert isinstance(utils.format_age_string(age), str)


@given(value=normal_float, threshold=finite, digits=st.integers(0, 6),
       sci=st.one_of(st.none(), st.integers(0, 6)), markdown=st.booleans())
def test_format_length_km_never_raises_on_finite_values(value, threshold, digits, sci, markdown):
    out = utils.format_length_km(_config(markdown), value, threshold, digits, sci)
    assert out.endswith(" km")
    if value <= threshold and "\u00d7" not in out:
        assert float(out[:-3].replace(",", "")) == round(value, digits)
    if value <= threshold and abs(round(value, digits)) >= 1e4:
        assert "\u00d7 10" in out  # UX.20: 5+ whole digits go scientific


@given(value=normal_float.filter(lambda v: v != 0), markdown=st.booleans(), low=st.integers(0, 8))
def test_format_relative_to_sol_picks_percent_or_multiplier(value, markdown, low):
    sol = physical_constants.SOLAR_MASS_TO_KG
    out = utils.format_relative_to_sol(_config(markdown), value, sol, "kg", low)
    ratio = value / sol
    if ratio < program_constants.PERCENT_SOL_THRESHOLD_HIGH:
        assert out.endswith("% of Sol)")
    else:
        assert out.endswith("× Sol)")


_PERIOD_TEXT = re.compile(r"^-?[0-9.,]+(?: \u00d7 10[\u207b\u2070\u00b9\u00b2\u00b3\u2074-\u2079]+)? "
                          r"(\u00b5s|ms|s|minutes?|hours?|days?|years?|ky|My|Gy)$")


@given(years=st.floats(min_value=0, max_value=1e15))
def test_format_period_years_is_one_number_and_one_unit(years):
    text = utils.format_period_years(years)
    match = _PERIOD_TEXT.match(text)
    assert match, (years, text)
    unit = match.group(1)
    if unit.rstrip("s") in ("minute", "hour", "day", "year"):
        # Singular exactly when the number shown is 1.
        assert (text.split(" ")[0] == "1") == (not unit.endswith("s")), text


@given(years=st.floats(min_value=-1e9, max_value=0))
def test_format_period_years_non_positive_does_not_raise(years):
    assert isinstance(utils.format_period_years(years), str)


@given(sentences=st.lists(hostile_text, max_size=6))
def test_to_paragraph_joins_with_single_spaces(sentences):
    assert utils.to_paragraph(sentences) == " ".join(sentences)


@given(props=st.dictionaries(st.text(alphabet="abcdefghij_", min_size=1, max_size=12), hostile_text, max_size=8),
       template=hostile_text, header=st.one_of(st.none(), hostile_text),
       key_map=st.one_of(st.none(), st.dictionaries(st.text(max_size=5), hostile_text, max_size=3)))
def test_properties_to_string_shapes(props, template, header, key_map):
    wiki = utils.properties_to_string(_config(False), props, template, header, key_map)
    assert wiki.startswith("{{" + template + "\n") and wiki.endswith("\n}}")
    for key, value in props.items():
        assert f"\n|{key}={value}\n" in wiki
    md = utils.properties_to_string(_config(True), props, template, header, key_map)
    body = md.split("| Property | Value |\n|---|---|", 1)
    assert len(body) == 2
    if header:
        assert md.startswith(header)
    # The rendered page must stay inert whatever the values hold.
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "html", "lib"))
    import mdconvert
    assert "<script" not in mdconvert.markdown_to_html(md).lower()


# ---------------------------------------------------------------------------
# utils: physics helpers
# ---------------------------------------------------------------------------

@given(low=finite, high=finite, mode=any_float, spread=st.floats(0.5, 10))
def test_sample_bounded_bell_stays_in_bounds(low, high, mode, spread):
    value = utils.sample_bounded_bell(low, high, mode, spread, max_attempts=50)
    if high <= low:
        assert value == low
    else:
        assert low <= value <= high


@given(lum=st.floats(min_value=0, max_value=1e33))
def test_habitable_zone_and_snow_line_are_ordered_and_scale_with_sqrt(lum):
    inner, outer = utils.calculate_habitable_zone(lum)
    assert 0 <= inner <= outer
    snow = utils.snow_line_au(lum)
    assert snow >= 0
    if lum > 0:
        inner4, outer4 = utils.calculate_habitable_zone(4 * lum)
        assert _close(inner4, 2 * inner) and _close(outer4, 2 * outer)
        assert _close(utils.snow_line_au(4 * lum), 2 * snow)


def test_solar_reference_points():
    inner, outer = utils.calculate_habitable_zone(physical_constants.SOLAR_LUMINOSITY)
    assert 0.9 < inner < 1.0 < outer < 1.5
    assert _close(utils.disk_surface_density_scale(physical_constants.SOLAR_MASS_TO_KG), 1.0)
    assert _close(utils.snow_line_au(physical_constants.SOLAR_LUMINOSITY),
                  physical_constants.SNOW_LINE_AU_AT_1_LSUN)


@given(d=positive, m=positive, big=positive)
def test_hill_sphere_scaling(d, m, big):
    r = utils.calculate_hill_sphere(d, m, big)
    assert r > 0
    assert _close(utils.calculate_hill_sphere(2 * d, m, big), 2 * r)
    assert _close(utils.calculate_hill_sphere(d, 8 * m, big), 2 * r)


@given(sep=positive, mu=any_float, e=any_float)
def test_holman_wiegert_is_clamped_and_positive(sep, mu, e):
    assume(not math.isnan(mu) and not math.isnan(e))
    a = utils.holman_wiegert_critical_semimajor_axis(sep, mu, e)
    mu_lo, mu_hi = physical_constants.HOLMAN_WIEGERT_MU_RANGE
    e_lo, e_hi = physical_constants.HOLMAN_WIEGERT_ECCENTRICITY_RANGE
    clamped = utils.holman_wiegert_critical_semimajor_axis(sep, min(max(mu, mu_lo), mu_hi), min(max(e, e_lo), e_hi))
    assert a == clamped
    assert 0 < a < sep


@given(m1=positive, m2=positive, d1=positive, d2=positive, big=positive)
def test_mutual_hill_radius_units_agree_and_are_symmetric(m1, m2, d1, d2, big):
    au = utils.mutual_hill_radius_au(m1, m2, d1, d2, big)
    m = utils.mutual_hill_radius_m(d1 * physical_constants.AU_TO_M, d2 * physical_constants.AU_TO_M, m1, m2, big)
    assert _close(m, au * physical_constants.AU_TO_M, rel=1e-9)
    assert _close(au, utils.mutual_hill_radius_au(m2, m1, d2, d1, big))
    assert au > 0


def test_wide_binary_samplers_stay_in_range():
    lo, hi = program_constants.WIDE_BINARY_SEPARATION_MIN_AU, program_constants.WIDE_BINARY_SEPARATION_MAX_AU
    for _ in range(scaled(60) * 20):
        assert lo <= utils.sample_wide_binary_separation_au() <= hi
        assert 0 <= utils.sample_wide_binary_eccentricity() <= program_constants.WIDE_BINARY_ECCENTRICITY_MAX


@given(d=positive, snow=positive, scale=st.floats(0, 1e3))
def test_mmsn_density_jumps_by_the_ice_factor_at_the_snow_line(d, snow, scale):
    density = utils.mmsn_surface_density_gcm2(d, snow, scale)
    bare = physical_constants.MMSN_SOLID_SURFACE_DENSITY_SOL_GCM2 * scale * d ** physical_constants.MMSN_SURFACE_DENSITY_EXPONENT
    expected = bare * (physical_constants.SNOW_LINE_ICE_BOOST_FACTOR if d >= snow else 1)
    assert _close(density, expected)
    assert density >= 0


@given(d=st.floats(1e-3, 1e4), sigma=st.floats(1e-6, 1e4), star=st.floats(1e27, 1e33))
def test_isolation_mass_scaling(d, sigma, star):
    m = utils.isolation_mass_kg(d, sigma, star)
    assert m > 0 and math.isfinite(m)
    assert _close(utils.isolation_mass_kg(d, 4 * sigma, star), 8 * m, rel=1e-9)
    assert _close(utils.isolation_mass_kg(d, sigma, 4 * star), m / 2, rel=1e-9)


@given(distance=st.one_of(st.floats(max_value=0, allow_nan=False), st.floats(1e-6, 1e7)))
def test_galactic_orbit(distance):
    speed, period = utils.calculate_galactic_orbit(distance)
    if distance <= 0:
        assert (speed, period) == (0.0, 0.0)
        return
    assert 0 < speed <= physical_constants.GALACTIC_ROTATION_FLAT_VELOCITY_KMS
    assert period > 0 and math.isfinite(period)
    assert isinstance(utils.format_galactic_orbit(speed, period), str)


@given(distance=st.one_of(st.none(), st.floats(1, 1e6)), phase=st.one_of(st.none(), angle))
def test_generate_galactic_orbit_fields(distance, phase):
    speed, period, phase_out, interval = utils.generate_galactic_orbit_fields(distance, phase)
    assert speed > 0 and period > 0 and interval > 0
    if phase is None:
        assert 0 <= phase_out <= 360
    else:
        assert phase_out == phase


@given(d=positive, period=positive)
def test_circular_orbital_speed_and_update_interval(d, period):
    v = utils.circular_orbital_speed_kms(d, period)
    assert v > 0 and _close(utils.circular_orbital_speed_kms(2 * d, period), 2 * v)
    interval = utils.minimum_update_interval_years(period)
    assert interval > 0 and _close(utils.minimum_update_interval_years(2 * period), 2 * interval)
    # Advancing by the interval always moves the phase by at least one ulp near 360.
    assert interval / period * 360 >= math.ulp(360.0) * (1 - 1e-12)


@given(r=st.floats(0, 1e6), inc=angle, node=angle, phase=angle)
def test_orbital_position_lies_on_the_orbit(r, inc, node, phase):
    x, y, z = utils.orbital_position_au(r, inc, node, phase)
    assert _close(math.sqrt(x * x + y * y + z * z), r, rel=1e-9, abs_=1e-9)
    _x0, _y0, z0 = utils.orbital_position_au(r, 0.0, node, phase)
    assert z0 == 0.0 or abs(z0) < 1e-9 * max(r, 1)


@given(parent=positive, children=st.lists(st.tuples(st.floats(0, 1e6), finite, finite, finite), max_size=6))
def test_reflex_offset_is_a_sum_of_two_body_terms(parent, children):
    ox, oy, oz = utils.calculate_reflex_offset(parent, children)
    ex = ey = ez = 0.0
    for m, x, y, z in children:
        mu = m / (parent + m)
        ex -= mu * x
        ey -= mu * y
        ez -= mu * z
    assert (ox, oy, oz) == (ex, ey, ez)
    assert utils.calculate_reflex_offset(parent, []) == (0.0, 0.0, 0.0)


@given(radius=st.floats(0, 1e6), density=st.floats(0, 100))
def test_calculate_object_mass_with_given_density(radius, density):
    volume, mass = utils.calculate_object_mass("M", radius, {}, {}, object_density=density)
    assert _close(volume, 4 / 3 * math.pi * radius ** 3)
    assert _close(mass, volume * physical_constants.KM_TO_M_FACTOR ** 3 * density * 1000)
    assert mass >= 0


@pytest.mark.parametrize("planet_class", sorted(program_constants.PLANET_CLASSES))
def test_calculate_object_mass_random_density_within_class_range(planet_class):
    kind = program_constants.PLANET_CLASSES[planet_class]["type"]
    lo, hi = physical_constants.PLANET_DENSITY[kind]
    volume, mass = utils.calculate_object_mass(planet_class, 1000.0, program_constants.PLANET_CLASSES,
                                               physical_constants.PLANET_DENSITY)
    density = mass / (volume * physical_constants.KM_TO_M_FACTOR ** 3 * 1000)
    assert lo - 1e-9 <= density <= hi + 1e-9


@given(letter=st.characters(), yerkes=st.sampled_from(["V", "III", "Ia", "VII", "IV", "sd", ""]),
       age=st.floats(0, 20), lifespan=st.one_of(st.just(math.inf), st.floats(0, 1e4)),
       proxy=st.booleans())
def test_star_class_helpers_on_fake_stars(letter, yerkes, age, lifespan, proxy):
    star = types.SimpleNamespace(type=letter + "2" + yerkes, yerkes_class=yerkes, age=age, lifespan=lifespan)
    target = types.SimpleNamespace(_primary=star) if proxy else star
    assert utils.get_star_spectral_class(target) == letter.upper()  # (may be 2 chars, e.g. "ß" -> "SS")
    profile = utils.get_star_evolutionary_profile(target)
    base = program_constants.STAR_EVOLUTION.get(letter.upper(), {})
    if not base:
        assert profile == {}
    elif yerkes == "V":
        assert profile is base
    else:
        scales = profile["supported_evolutionary_scales"]
        assert scales and set(scales) <= {"fast", "normal", "slow"}


# ---------------------------------------------------------------------------
# appconfig
# ---------------------------------------------------------------------------

json_scalar = st.one_of(st.none(), st.booleans(), st.integers(-2**63, 2**63), st.floats(allow_nan=False),
                        hostile_text)
json_value = st.recursive(json_scalar, lambda inner: st.one_of(st.lists(inner, max_size=3),
                                                              st.dictionaries(st.text(max_size=8), inner,
                                                                              max_size=3)),
                          max_leaves=10)
config_keys = st.sampled_from(sorted(appconfig.DEFAULT_CONFIG) + ["unknown", "mysql", "wiki"])
config_object = st.dictionaries(config_keys, st.one_of(json_value, st.dictionaries(
    st.sampled_from(["host", "port", "password", "wikijs", "dir", "x"]), json_value, max_size=3)), max_size=6)


def _load_from(tmp_dir, raw_text):
    path = os.path.join(tmp_dir, "config.json")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(raw_text)
    with mock.patch.object(appconfig, "CONFIG_PATH", path):
        return appconfig.load_config()


@settings(max_examples=scaled(60))
@given(overrides=config_object)
def test_load_config_merges_any_object_without_touching_defaults(tmp_path_factory, overrides):
    import copy
    before = copy.deepcopy(appconfig.DEFAULT_CONFIG)
    merged = _load_from(str(tmp_path_factory.mktemp("cfg")), json.dumps(overrides))
    assert appconfig.DEFAULT_CONFIG == before
    assert set(before) <= set(merged)
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(before.get(key), dict):
            for sub in before[key]:
                assert sub in merged[key]
        else:
            assert merged[key] == value
    # A second load is independent of the first result.
    merged["mysql"] = "mutated"
    assert _load_from(str(tmp_path_factory.mktemp("cfg")), json.dumps(overrides)).get("mysql") != "mutated" \
        or overrides.get("mysql") == "mutated"


@pytest.mark.parametrize("raw", ["[]", "null", "1", '"text"', "true", "", "{", "{\"a\": }", "﻿{}",
                                 "NaN"])
def test_load_config_on_broken_or_non_object_json_fails_loudly_not_silently(tmp_path, raw):
    """A broken config.json must never silently load as something that
    isn't the full default shape: either raise, or return every key."""
    with pytest.raises(ValueError):  # json.JSONDecodeError is a ValueError too
        _load_from(str(tmp_path), raw)


@pytest.mark.parametrize("raw", ["[]", "null", "1", '"text"', "true"])
def test_load_config_names_the_problem_for_non_object_json(tmp_path, raw):
    with pytest.raises(ValueError, match="must contain a JSON object"):
        _load_from(str(tmp_path), raw)


def _configure_with_config(tmp_path, raw):
    path = tmp_path / "config.json"
    path.write_text(raw, encoding="utf-8")
    with mock.patch.dict(os.environ, {"PLANETGEN_LOG_FILE": ""}):
        os.environ.pop("PLANETGEN_DEBUG", None)
        os.environ.pop("PLANETGEN_LOG_FILE", None)
        try:
            with mock.patch.object(appconfig, "CONFIG_PATH", str(path)):
                log.configure(log.NORMAL)
        finally:
            with mock.patch.object(appconfig, "CONFIG_PATH", str(tmp_path / "missing.json")):
                log.configure(log.NORMAL)


@pytest.mark.parametrize("raw", ["[]", "{", '{"debug": {"x": 1}, "log_file": "/nonexistent-dir/x.log"}',
                                 '{"debug": false, "log_file": 5}', '"just a string"'])
def test_configure_survives_a_hostile_config_file(tmp_path, raw):
    _configure_with_config(tmp_path, raw)
    assert not log.debug_log_active()


@given(log_file=json_value.filter(lambda v: bool(v) and not isinstance(v, str)))
@example(log_file=["not", "a", "path"])   # was TypeError from WatchedFileHandler
@example(log_file=5)
@example(log_file={"a": 1})
@settings(max_examples=scaled(15))
def test_configure_survives_a_non_string_log_file(tmp_path_factory, log_file):
    _configure_with_config(tmp_path_factory.mktemp("cfg"), json.dumps({"debug": True, "log_file": log_file}))
    assert not log.debug_log_active()


@given(env=st.text(alphabet=st.characters(blacklist_categories=("Cs",), blacklist_characters="\x00"), max_size=10),
       config_debug=json_value)
def test_debug_enabled_env_wins_and_returns_bool(env, config_debug):
    with mock.patch.dict(os.environ, {"PLANETGEN_DEBUG": env}):
        result = appconfig.debug_enabled({"debug": config_debug})
    assert result is (env.strip().lower() not in appconfig._FALSE_STRINGS)
    os.environ.pop("PLANETGEN_DEBUG", None)
    with mock.patch.dict(os.environ, {}):
        os.environ.pop("PLANETGEN_DEBUG", None)
        assert isinstance(appconfig.debug_enabled({"debug": config_debug}), bool)
        assert appconfig.debug_enabled({}) is False


@given(value=st.sampled_from(appconfig._FALSE_STRINGS), case=st.sampled_from([str.upper, str.lower, str.title]),
       pad=st.sampled_from(["", " ", "\t", " \n"]))
@example(value="false", case=str.lower, pad="")   # was True: bool("false")
@example(value="0", case=str.lower, pad="")
@example(value="off", case=str.lower, pad="")
@example(value="no", case=str.title, pad="")
def test_debug_enabled_treats_false_strings_in_config_like_the_env_var(value, case, pad):
    with mock.patch.dict(os.environ, {}):
        os.environ.pop("PLANETGEN_DEBUG", None)
        assert appconfig.debug_enabled({"debug": pad + case(value) + pad}) is False
        assert appconfig.debug_enabled({"debug": "yes"}) is True
        assert appconfig.debug_enabled({"debug": True}) is True


@given(env=st.one_of(st.none(), st.just(""), st.text(alphabet="abc/._-", min_size=1, max_size=20)),
       configured=json_value)
def test_log_file_path_precedence(env, configured):
    with mock.patch.dict(os.environ, {}):
        os.environ.pop("PLANETGEN_LOG_FILE", None)
        if env is not None:
            os.environ["PLANETGEN_LOG_FILE"] = env
        out = appconfig.log_file_path({"log_file": configured})
    if env:
        assert out == env
    elif configured:
        assert out == configured
    else:
        assert out == appconfig.DEFAULT_CONFIG["log_file"]


# ---------------------------------------------------------------------------
# config.SystemConfig / serialization
# ---------------------------------------------------------------------------

@given(values=st.fixed_dictionaries({}, optional={f.lower(): json_value for f in SERIALIZABLE_FIELDS}),
       extra=st.dictionaries(st.text(max_size=10), json_value, max_size=3))
def test_system_config_round_trip_with_arbitrary_values(values, extra):
    data = {**{k: v for k, v in extra.items() if k not in values}, **values}
    config = SystemConfig.from_dict(data)
    defaults = SystemConfig()
    out = config.to_dict()
    assert set(out) == {f.lower() for f in SERIALIZABLE_FIELDS}
    for field in SERIALIZABLE_FIELDS:
        key = field.lower()
        assert out[key] == (data[key] if key in data else getattr(defaults, field))
    # Unknown keys never become attributes; run-time bookkeeping untouched.
    for key in extra:
        if key not in out and key.isidentifier():
            assert not hasattr(config, key) or hasattr(defaults, key)
    assert config.system_flavor_count == 0 and config.recent_flavor_texts == []
    assert SystemConfig.from_dict(json.loads(json.dumps(out))).to_dict() == json.loads(json.dumps(out))


@given(upper=st.dictionaries(st.sampled_from(SERIALIZABLE_FIELDS), json_value, min_size=1))
def test_system_config_ignores_uppercase_keys(upper):
    assert SystemConfig.from_dict(upper).to_dict() == SystemConfig().to_dict()


@given(fields=st.lists(st.from_regex(r"[A-Za-z_][A-Za-z0-9_]{0,10}", fullmatch=True), unique_by=str.lower,
                       max_size=8),
       data=st.data())
def test_fields_round_trip_on_any_object(fields, data):
    values = {f: data.draw(json_value) for f in fields}
    source = types.SimpleNamespace(**values)
    as_dict = serialization.fields_to_dict(source, fields)
    assert as_dict == {f.lower(): v for f, v in values.items()}
    target = types.SimpleNamespace(untouched=1)
    serialization.fields_from_dict(target, as_dict, fields)
    assert all(getattr(target, f) == values[f] for f in fields)
    assert target.untouched == 1
    partial = types.SimpleNamespace(**{f: "keep" for f in fields})
    serialization.fields_from_dict(partial, {}, fields)
    assert all(getattr(partial, f) == "keep" for f in fields)


def test_fields_to_dict_raises_on_a_missing_attribute():
    with pytest.raises(AttributeError):
        serialization.fields_to_dict(types.SimpleNamespace(), ["missing"])


def _serializable_classes():
    import importlib
    import pkgutil
    import stellarObjects
    found = {"config.SystemConfig": SERIALIZABLE_FIELDS}
    for info in pkgutil.iter_modules(stellarObjects.__path__):
        module = importlib.import_module(f"stellarObjects.{info.name}")
        for name, obj in vars(module).items():
            if isinstance(obj, type) and obj.__module__ == module.__name__ and "SERIALIZABLE_FIELDS" in vars(obj):
                found[f"{info.name}.{name}"] = obj.SERIALIZABLE_FIELDS
    return found


def test_every_serializable_fields_list_is_collision_free_when_lowercased():
    classes = _serializable_classes()
    assert len(classes) >= 10
    for label, fields in classes.items():
        lowered = [f.lower() for f in fields]
        assert len(set(lowered)) == len(lowered), (label, sorted(f for f in lowered if lowered.count(f) > 1))
        assert all(f.isidentifier() for f in fields), label


# ---------------------------------------------------------------------------
# progressFile
# ---------------------------------------------------------------------------

@pytest.fixture
def progress_env(tmp_path):
    path = tmp_path / "progress.json"
    with mock.patch.dict(os.environ, {progressFile.ENV_VAR: str(path)}), \
            mock.patch.object(progressFile, "_last_write", 0.0):
        yield path


@given(completed=st.one_of(st.integers(-10**18, 10**18), finite),
       total=st.one_of(st.none(), st.integers(0, 10**18), finite),
       description=st.one_of(st.none(), hostile_text))
def test_progress_report_writes_valid_json(tmp_path_factory, completed, total, description):
    path = tmp_path_factory.mktemp("progress") / "progress.json"
    with mock.patch.dict(os.environ, {progressFile.ENV_VAR: str(path)}):
        progressFile.report(completed, total, description, force=True)
    body = json.loads(path.read_text(encoding="utf-8"))
    assert body["completed"] == completed and body["total"] == total and body["description"] == description
    assert isinstance(body["updated_at"], float)
    assert [p.name for p in path.parent.iterdir()] == ["progress.json"]


def test_progress_report_is_throttled_unless_forced(progress_env):
    progressFile.report(1, 10, "x", force=True)
    first = progress_env.read_text()
    progressFile.report(2, 10, "x")
    assert progress_env.read_text() == first
    progressFile.report(3, 10, "x", force=True)
    assert json.loads(progress_env.read_text())["completed"] == 3


def test_progress_report_without_env_is_a_no_op(tmp_path):
    with mock.patch.dict(os.environ, {}):
        os.environ.pop(progressFile.ENV_VAR, None)
        progressFile.report(1, force=True)
        with mock.patch.dict(os.environ, {progressFile.ENV_VAR: ""}):
            progressFile.report(1, force=True)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("target", ["missing-dir/progress.json", "\x01weird/../p.json"])
def test_progress_report_never_raises_on_bad_paths(tmp_path, target):
    with mock.patch.dict(os.environ, {progressFile.ENV_VAR: str(tmp_path / target)}):
        progressFile.report(1, 2, "x", force=True)


def test_progress_report_cleans_up_its_temp_file_when_the_rename_fails(tmp_path):
    """Regression: with the target path a directory, os.replace failed and
    the `.progress-*` temp file was left behind."""
    target = tmp_path / "progress.json"
    target.mkdir()
    with mock.patch.dict(os.environ, {progressFile.ENV_VAR: str(target)}):
        progressFile.report(1, 2, "x", force=True)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["progress.json"]


@given(value=st.one_of(st.builds(object), st.sets(st.integers(), max_size=2), st.binary(max_size=4),
                       st.just(math.nan), st.just(math.inf)),
       description=st.one_of(st.none(), st.builds(object), hostile_text))
@example(value=object(), description=None)   # was TypeError from json.dump, plus a leaked temp file
@example(value={1, 2}, description=None)
@example(value=b"bytes", description=None)
def test_progress_report_never_raises_on_unserializable_values(tmp_path_factory, value, description):
    tmp = tmp_path_factory.mktemp("progress")
    with mock.patch.dict(os.environ, {progressFile.ENV_VAR: str(tmp / "p.json")}):
        progressFile.report(value, description=description, force=True)
    assert not [p for p in tmp.iterdir() if p.name.startswith(".progress-")]


# ---------------------------------------------------------------------------
# log.py / _db.py credential redaction
# ---------------------------------------------------------------------------

secret_value = st.text(alphabet="abcdefghijklmnopqrstuvwxyz0123456789!@#", min_size=12, max_size=24)
benign_arg = st.text(alphabet="abcdefghijklmnop-=_/.", max_size=12).filter(lambda a: "password" not in a.lower())


@given(secret=secret_value, before=st.lists(benign_arg, max_size=4), after=st.lists(benign_arg, max_size=4),
       flag=st.sampled_from(["--mysql-password", "--MYSQL-PASSWORD", "--password", "-password", "--new-password"]),
       joined=st.booleans())
def test_redacted_argv_never_shows_a_password_value(secret, before, after, flag, joined):
    assume(not any(secret in a for a in before + after))
    argv = ["generate.py", *before, *([f"{flag}={secret}"] if joined else [flag, secret]), *after]
    with mock.patch.object(sys, "argv", argv):
        shown = log._redacted_argv()
    assert len(shown) == len(argv)
    assert secret not in " ".join(shown)
    assert "<withheld>" in " ".join(shown)


@given(argv=st.lists(hostile_text, max_size=8))
def test_redacted_argv_never_raises(argv):
    with mock.patch.object(sys, "argv", argv):
        assert len(log._redacted_argv()) == len(argv)


@given(secret=secret_value, cut=st.integers(len("--mysql-p"), len("--mysql-password")), joined=st.booleans())
@example(secret="hunter2secret", cut=len("--mysql-pass"), joined=False)   # was logged in the clear
def test_redacted_argv_covers_argparse_abbreviations(secret, cut, joined):
    option = "--mysql-password"[:cut]
    argv = ["generate.py", *([f"{option}={secret}"] if joined else [option, secret])]
    parser = argparse.ArgumentParser()
    _db.add_mysql_connection_args(parser)
    if cut > len("--mysql-p"):  # "--mysql-p" alone is ambiguous with --mysql-port
        assert parser.parse_args(argv[1:]).mysql_password == secret  # argparse really takes it
    with mock.patch.object(sys, "argv", argv):
        assert secret not in " ".join(log._redacted_argv())


def test_redacted_argv_leaves_other_options_alone():
    argv = ["generate.py", "--mysql-host", "db", "--mysql-port", "3306", "--mysql-user", "u", "--name", "X"]
    with mock.patch.object(sys, "argv", argv):
        assert log._redacted_argv() == argv


@given(secret=secret_value, column=st.sampled_from(["password_hash", "token_hash", "key_hash", "PASSWORD",
                                                    "api_token", "secret_key"]),
       padding=st.text(alphabet=" \n\t", max_size=3))
def test_sql_log_line_withholds_params_of_credential_statements(secret, column, padding):
    sql = f"UPDATE admin_users SET{padding}{column} = ? WHERE id = ?"
    line = _db._sql_for_log(sql, (secret, 1))
    assert secret not in line and "<withheld: credentials>" in line
    assert "\n" not in line


@given(sql=hostile_text, params=st.lists(hostile_text, max_size=4))
def test_sql_log_line_is_one_bounded_line(sql, params):
    line = _db._sql_for_log(sql, tuple(params))
    assert "\n" not in line.split(" | params=")[0]
    assert len(line) <= _db._SQL_LOG_LIMIT + 40


def test_debug_log_file_never_contains_command_line_password(tmp_path):
    log_path = tmp_path / "planetgen.log"
    secret = "S3cr3t-Value-For-Redaction"
    argv = ["generate.py", "--mysql-password", secret, f"--mysql-password={secret}x"]
    with mock.patch.dict(os.environ, {"PLANETGEN_DEBUG": "1", "PLANETGEN_LOG_FILE": str(log_path)}), \
            mock.patch.object(sys, "argv", argv):
        try:
            log.configure(log.NORMAL)
            log.trace("SQL line: %s", _db._sql_for_log("UPDATE admin_users SET password_hash = ?", (secret,)))
            for handler in logging.getLogger("planetgen").handlers:
                handler.flush()
        finally:
            os.environ.pop("PLANETGEN_DEBUG", None)
            log.configure(log.NORMAL)
    text = log_path.read_text(encoding="utf-8")
    assert "Debug log opened" in text and "<withheld>" in text
    assert secret not in text


# ---------------------------------------------------------------------------
# physical_constants / program_constants sanity
# ---------------------------------------------------------------------------

def _public_constants(module):
    return {k: v for k, v in vars(module).items() if k.isupper() and not k.startswith("_")}


def _walk(value, path):
    """Yields (path, value) for every leaf and every numeric 2-sequence."""
    if isinstance(value, dict):
        for key, item in value.items():
            yield from _walk(item, f"{path}[{key!r}]")
    elif isinstance(value, (list, tuple)):
        if len(value) == 2 and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in value):
            yield path, tuple(value)
        for i, item in enumerate(value):
            yield from _walk(item, f"{path}[{i}]")
    else:
        yield path, value


_EXPONENT_OR_SIGNED = re.compile(r"EXPONENT|ROUND_|DEFAULT_M$|_OFFSET|_DELTA|_SLOPE|_INTERCEPT")


@pytest.mark.parametrize("module", [physical_constants, program_constants], ids=lambda m: m.__name__)
def test_every_numeric_constant_is_finite_and_every_pair_ordered(module):
    for name, value in _public_constants(module).items():
        for path, leaf in _walk(value, name):
            if isinstance(leaf, tuple):
                assert leaf[0] <= leaf[1], (path, leaf)
            elif isinstance(leaf, float):
                assert math.isfinite(leaf), (path, leaf)
            elif isinstance(leaf, str):
                assert leaf == leaf.strip() and leaf, (path, leaf)


def test_top_level_scalar_constants_are_positive_unless_signed_by_nature():
    for module in (physical_constants, program_constants):
        for name, value in _public_constants(module).items():
            if isinstance(value, (int, float)) and not isinstance(value, bool) and not _EXPONENT_OR_SIGNED.search(name):
                assert value > 0, (module.__name__, name, value)


@pytest.mark.parametrize("module", [physical_constants, program_constants], ids=lambda m: m.__name__)
def test_named_min_max_pairs_are_ordered(module):
    constants = _public_constants(module)
    pairs = 0
    for name, value in constants.items():
        for lo, hi in (("_MIN", "_MAX"), ("MIN_", "MAX_"), ("_LOW", "_HIGH"), ("_INNER", "_OUTER")):
            if lo in name and name.replace(lo, hi) in constants:
                other = constants[name.replace(lo, hi)]
                if isinstance(value, (int, float)) and isinstance(other, (int, float)):
                    assert value <= other, (name, value, other)
                    pairs += 1
    if module is program_constants:
        assert pairs >= 5


def test_conversion_constants_are_mutually_consistent():
    pc = physical_constants
    assert _close(pc.LY_TO_AU * pc.AU_TO_LY, 1.0, rel=1e-9)
    assert _close(pc.AU_TO_KM * pc.KM_TO_M_FACTOR, pc.AU_TO_M, rel=1e-9)
    assert _close(pc.LY_TO_M, pc.LY_TO_AU * pc.AU_TO_M, rel=1e-4)  # both rounded to 4 significant digits
    assert _close(pc.AU_PER_PARSEC, 1000 * pc.AU_PER_MILLIPARSEC, rel=1e-9)
    assert 3.15e7 < pc.SECONDS_PER_YEAR < 3.16e7
    lo, hi = pc.HOLMAN_WIEGERT_MU_RANGE
    assert 0 <= lo <= hi <= 1
    lo, hi = pc.HOLMAN_WIEGERT_ECCENTRICITY_RANGE
    assert 0 <= lo <= hi < 1


def test_probability_tables_are_well_formed():
    pg = program_constants
    for name in ("SPECTRAL_PROBABILITIES_LARGE_STAR", "SPECTRAL_PROBABILITIES_NORMAL", "PLANET_CLASS_PROBABILITIES",
                 "PHENOMENON_DENSITY_PC3", "PHENOMENON_RATE_SCALE"):
        weights = getattr(pg, name)
        assert weights and all(w >= 0 for w in weights.values()), name
        assert sum(weights.values()) > 0, name
    for letter, chance in pg.BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS.items():
        assert 0 <= chance <= 1, letter
    assert 0 <= pg.WIDE_BINARY_ECCENTRICITY_MAX < 1
    assert 0 <= pg.WIDE_BINARY_DEFAULT_CHANCE <= 1


def test_lookup_tables_reference_each_other_consistently():
    pg = program_constants
    pc = physical_constants
    assert set(pg.PLANET_CLASS_PROBABILITIES) <= set(pg.PLANET_CLASSES)
    for letter, spec in pg.PLANET_CLASSES.items():
        assert spec["type"] in pc.PLANET_DENSITY, (letter, spec["type"])
    spectral = set(pc.TEMP_RANGES)
    for table in (pc.SPECTRAL_LUMINOSITY_RANGES, pc.SPECTRAL_MASS_RANGES):
        assert set(table) == spectral
    assert set(pg.SPECTRAL_PROBABILITIES_NORMAL) <= spectral
    assert set(pg.SPECTRAL_PROBABILITIES_LARGE_STAR) <= spectral
    assert set(pg.BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS) <= spectral
    for scale in ("fast", "normal", "slow"):
        assert pg.EVOLUTIONARY_TIMELINES[scale]["technological_civilization"] > 0
    fast, normal, slow = (pg.EVOLUTIONARY_TIMELINES[s]["technological_civilization"] for s in ("fast", "normal", "slow"))
    assert fast <= normal <= slow
