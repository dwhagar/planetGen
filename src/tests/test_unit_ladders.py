# tests/test_unit_ladders.py

"""
The shared unit ladders (UX.23): `planetgen.util.ladders` holds the tables,
`planetgen.util.format.format_quantity` is the one algorithm, and
`html/static/unitladder.js` is its mirror in the browser, fed by the
generated `unitladders.js`. These tests hold the table itself (sorted, named,
every factor against astropy), the properties of the algorithm (a value is
shown within rounding of itself, a larger value never picks a smaller rung,
nothing raises or prints "nan"), and that the browser gives the same text as
Python for thousands of values per ladder, ties and rung boundaries included
(needs Node, skipped without it).
"""

import json
import math
import os
import random
import re
import shutil
import subprocess
import sys

import pytest
from hypothesis import given
from hypothesis import strategies as st

from planetgen.physics import constants
from planetgen.util import ladders
from planetgen.util.format import format_quantity, pick_rung
from tests import fuzz_support  # noqa: F401  (loads the hypothesis profile)

_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
_STATIC = os.path.join(_ROOT, "src", "html", "static")

_SUPERSCRIPTS = str.maketrans("⁻⁰¹²³⁴⁵⁶⁷⁸⁹", "-0123456789")
_TEXT = re.compile(r"^(-?[\d,]+(?:\.\d+)?)(?: × 10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+))? ([^()]+?)(?: \(.*\))?$")

magnitudes = st.builds(
    lambda mantissa, exponent, sign: sign * mantissa * 10.0 ** exponent,
    st.floats(min_value=1.0, max_value=9.999, allow_nan=False),
    st.integers(min_value=-14, max_value=32),
    st.sampled_from((1, -1)),
)


CASES = [(ladder, system) for ladder in ladders.LADDERS for system in ladders.SYSTEMS]
"""list: Every ladder under every unit preference (UX.78)."""


def _case_id(case):
    return f"{case[0].quantity}-{case[1]}"


def _primary(ladder, text):
    """The (value in the ladder's base unit, rung) a formatted `text` reads back as."""
    match = _TEXT.match(text)
    assert match, text
    number = float(match.group(1).replace(",", ""))
    if match.group(2):
        number *= 10.0 ** int(match.group(2).translate(_SUPERSCRIPTS))
    unit = match.group(3)
    for rung in ladder.rungs:
        if unit in (rung.symbol, rung.singular):
            return number * rung.factor + rung.offset, rung
    raise AssertionError(f"{unit!r} is on no rung of {ladder.quantity}")


@pytest.mark.parametrize("ladder, system", CASES, ids=[_case_id(case) for case in CASES])
class TestTable:
    def test_rungs_are_sorted_and_named(self, ladder, system):
        ladder = ladder.for_system(system)
        factors = [rung.factor for rung in ladder.rungs]
        starts = [rung.factor if rung.start is None else rung.start for rung in ladder.rungs]
        # A rung may be a smaller unit than the one below it when it starts higher (gauss after tesla).
        assert len(set(factors)) == len(factors)
        assert starts == sorted(starts)
        assert len({rung.symbol for rung in ladder.rungs}) == len(ladder.rungs)
        assert all(rung.name for rung in ladder.rungs)

    def test_the_lowest_rung_takes_everything_below_it(self, ladder, system):
        ladder = ladder.for_system(system)
        assert ladder.rungs[0].start == 0.0

    def test_asides_are_complete(self, ladder, system):
        ladder = ladder.for_system(system)
        carriers = [rung for rung in ladder.rungs if rung.with_aside]
        assert bool(carriers) == bool(ladder.asides)
        if ladder.asides:
            sizes = [size for size, _ in ladder.asides]
            assert sizes == sorted(sizes, reverse=True) and sizes[-1] == 0.0
            for _, symbol in ladder.asides:
                assert ladder.rung(symbol)

    def test_zero_names_a_rung(self, ladder, system):
        ladder = ladder.for_system(system)
        if ladder.zero is not None:
            assert ladder.rung(ladder.zero)


def _astropy_factor(symbol):
    """How many SI base units one `symbol` is, by astropy; None when it has no unit of that name."""
    u = pytest.importorskip("astropy.units")
    from astropy import constants as const
    named = {"km": u.km, "AU": u.au, "mpc": u.mpc, "cpc": 0.01 * u.pc, "ly": u.lyr, "pc": u.pc, "kpc": u.kpc,
             "Mpc": u.Mpc, "Gpc": u.Gpc, "km/h": u.km / u.hr, "km/s": u.km / u.s, "Mm/s": u.Mm / u.s,
             "c": const.c, "µs": u.us, "ms": u.ms, "s": u.s, "minutes": u.min, "hours": u.hr, "days": u.day,
             "years": u.year, "kyr": u.kyr, "Myr": u.Myr, "Gyr": u.Gyr,
             "nK": 1e-9 * u.K, "µK": u.uK, "mK": u.mK, "K": u.K,
             "nPa": 1e-9 * u.Pa, "µPa": u.uPa, "mPa": u.mPa, "Pa": u.Pa, "kPa": u.kPa, "MPa": u.MPa, "GPa": u.GPa,
             "TPa": 1e12 * u.Pa, "m/s²": u.m / u.s ** 2, "bar": u.bar, "mbar": u.mbar,
             "g": const.g0,
             # The IAU defines no lunar mass: 7.342e22 kg is project-nominal (constants.MOON_MASS_KG).
             "W": u.W, "kW": u.kW, "MW": u.MW, "GW": u.GW, "TW": 1e12 * u.W, "PW": 1e15 * u.W, "EW": 1e18 * u.W,
             "L☉": const.L_sun, "nW/m²": 1e-9 * u.W / u.m ** 2, "µW/m²": 1e-6 * u.W / u.m ** 2,
             "mW/m²": 1e-3 * u.W / u.m ** 2, "W/m²": u.W / u.m ** 2, "kW/m²": 1e3 * u.W / u.m ** 2,
             "MW/m²": 1e6 * u.W / u.m ** 2,
             "particles/cm³": 1e6 / u.m ** 3, "nT": u.nT, "µT": u.uT, "mT": u.mT, "T": u.T,
             "G": 1e-4 * u.T,  # astropy's gauss is a CGS unit, not convertible without an equivalency
             "ft": u.imperial.ft, "mi": u.imperial.mi, "lb": u.imperial.lb, "mph": u.imperial.mi / u.h,
             "mi/s": u.imperial.mi / u.s, "ft/s²": u.imperial.ft / u.s ** 2,
             "Gm": u.Gm, "Tm": u.Tm, "Pm": u.Pm, "Em": u.Em, "Zm": 1e21 * u.m, "Ym": 1e24 * u.m,
             "ZW": 1e21 * u.W, "YW": 1e24 * u.W, "RW": 1e27 * u.W, "QW": 1e30 * u.W,
             "°F": 5 / 9 * u.K,  # the size of a degree; the zero is checked against astropy's temperature scales below
             "psi": 4.4482216152605 / 0.0254 ** 2 * u.Pa,  # lbf (exact, in N) per square inch; astropy's imperial.psi is rounded
             "kg": u.kg, "lunar masses": 7.342e22 * u.kg, "Earth masses": const.M_earth,
             "Jupiter masses": const.M_jup, "solar masses": const.M_sun}
    quantity = named.get(symbol)
    if quantity is None:
        return None
    if not isinstance(quantity, u.Quantity):
        quantity = 1 * quantity
    return float(quantity.si.value)


@pytest.mark.parametrize("ladder", ladders.LADDERS, ids=lambda ladder: ladder.quantity)
def test_every_factor_agrees_with_astropy(ladder):
    for rung in ladder.rungs:
        expected = _astropy_factor(rung.symbol)
        assert expected is not None, f"add {rung.symbol} to the astropy table in this test"
        assert rung.factor == pytest.approx(expected, rel=1e-12), rung.symbol
    # No atm (astropy has none), psi (its imperial.psi is rounded, 1.4e-8 off), °C or °F (affine).
    for companion in ladder.companions:
        expected = _astropy_factor(companion.symbol)
        if expected is not None:
            assert companion.factor == pytest.approx(expected, rel=1e-12), companion.symbol


def test_the_temperature_companions_agree_with_astropys_scales():
    u = pytest.importorskip("astropy.units")
    for kelvin in (0.0, 100.0, 273.15, 288.0, 5778.0):
        celsius, fahrenheit = (ladders.TEMPERATURE.companions[0], ladders.TEMPERATURE.companions[1])
        assert (kelvin - celsius.offset) / celsius.factor == pytest.approx(
            (kelvin * u.K).to(u.deg_C, u.temperature()).value, abs=1e-9)
        assert (kelvin - fahrenheit.offset) / fahrenheit.factor == pytest.approx(
            (kelvin * u.K).to(u.imperial.deg_F, u.temperature()).value, abs=1e-9)


@pytest.mark.parametrize("kilograms, text", [
    (1_000.0, "1,000 kg"),
    (7.0e20, "7.00 × 10²⁰ kg"),
    (1.3e22, "0.177 lunar masses"),
    (7.35e22 * 4.49, "4.49 lunar masses"),
    (0.107 * constants.EARTH_MASS_TO_KG, "0.107 Earth masses"),
    (17.1 * constants.EARTH_MASS_TO_KG, "17.1 Earth masses"),
    (0.299 * constants.JUPITER_MASS_TO_KG, "0.299 Jupiter masses"),
    (0.122 * constants.SOLAR_MASS_TO_KG, "0.122 solar masses"),
    (4.3e6 * constants.SOLAR_MASS_TO_KG, "4.30 × 10⁶ solar masses"),
    (constants.EARTH_MASS_TO_KG, "1 Earth mass"),
])
def test_mass_switches_units_at_the_decided_points(kilograms, text):
    """The design note's prototype outputs (docs/design/units-and-number-formatting.md)."""
    assert format_quantity(ladders.MASS, kilograms) == text


@pytest.mark.parametrize("ladder, value, text", [
    (ladders.POWER, 9.01e-31, "9.01 × 10⁻³¹ W"),
    (ladders.POWER, 47e12, "47 TW"),
    (ladders.POWER, 335e15, "335 PW"),
    (ladders.POWER, constants.SOLAR_LUMINOSITY, "1 L☉"),
    (ladders.POWER, 0.00155 * constants.SOLAR_LUMINOSITY, "0.00155 L☉"),
    (ladders.POWER, 2.6e13 * constants.SOLAR_LUMINOSITY, "2.60 × 10¹³ L☉"),
    (ladders.POWER, 1e-7 * constants.SOLAR_LUMINOSITY, "38.3 EW"),
    (ladders.FLUX, 0.087, "87 mW/m²"),
    (ladders.FLUX, 1361.0, "1.36 kW/m²"),
    (ladders.NUMBER_DENSITY, 3.5e6, "3.5 particles/cm³"),
    (ladders.NUMBER_DENSITY, 2.5e15, "2.50 × 10⁹ particles/cm³"),
    (ladders.MAGNETIC_FIELD, 5e-5, "50 µT"),
    (ladders.MAGNETIC_FIELD, 0.3, "300 mT"),
    (ladders.MAGNETIC_FIELD, 1e8, "1.00 × 10¹² G"),
])
def test_other_ladders_switch_units_at_the_decided_points(ladder, value, text):
    assert format_quantity(ladder, value) == text


def test_the_mass_rungs_start_where_boss_decided():
    """UX.22, Boss 2026-10-10: lunar masses from 0.01 lunar mass, Earth masses from 0.1, Jupiter masses
    from 0.1 M_J, solar masses from the 0.075 M☉ hydrogen-burning limit."""
    starts = {rung.symbol: rung.start / rung.factor for rung in ladders.MASS.rungs
              if "auto" in rung.systems and rung.symbol != "kg"}
    assert starts == pytest.approx({"lunar masses": 0.01, "Earth masses": 0.1, "Jupiter masses": 0.1,
                                    "solar masses": 0.075})


def test_atm_and_psi_are_the_exact_definitions():
    atm, psi = ladders.PRESSURE.companions[:2]
    assert atm.factor == 101_325
    assert psi.factor == pytest.approx(4.4482216152605 / 0.0254 ** 2, rel=1e-15)  # lbf (exact, in N) per square inch


@pytest.mark.parametrize("ladder, system", CASES, ids=[_case_id(case) for case in CASES])
class TestProperties:
    @given(value=magnitudes)
    def test_the_text_reads_back_within_rounding(self, ladder, system, value):
        if value == 0:
            return
        shaped = ladder.for_system(system)
        read, rung = _primary(shaped, format_quantity(ladder, value, system=system))
        # Rounded to three figures in the rung's own unit, so a °F figure near its zero is a coarse kelvin value.
        assert abs(read - value) <= 0.0051 * abs(value - rung.offset) + 1e-9 * abs(value)
        assert rung is pick_rung(shaped, value)

    @given(first=magnitudes, second=magnitudes)
    def test_a_larger_size_never_picks_a_smaller_rung(self, ladder, system, first, second):
        ladder = ladder.for_system(system)
        low, high = sorted((abs(first), abs(second)))
        assert ladder.rungs.index(pick_rung(ladder, low)) <= ladder.rungs.index(pick_rung(ladder, high))

    @given(value=st.one_of(magnitudes, st.floats(), st.integers(), st.none(), st.sampled_from((0, 0.0, -0.0, 1e300, 10 ** 400))))
    def test_nothing_raises_and_only_non_numbers_are_a_dash(self, ladder, system, value):
        try:
            text = format_quantity(ladder, value, system=system)
        except OverflowError:
            assert isinstance(value, int)  # an int past float range is the caller's to guard (fmt.dash_unless_finite)
            return
        assert text
        if value is None or (isinstance(value, float) and math.isnan(value)):
            assert text == "–"
        else:
            assert "nan" not in text

    def test_each_rung_boundary_picks_that_rung(self, ladder, system):
        ladder = ladder.for_system(system)
        for rung in ladder.rungs[1:]:
            start = rung.factor if rung.start is None else rung.start
            assert pick_rung(ladder, start) is rung
            assert pick_rung(ladder, start * (1 - 1e-13)) is rung
            assert pick_rung(ladder, start * 0.999) is not rung


def _check_values(ladder):
    """A seeded batch for `ladder`: log-uniform values of both signs, every
    rung boundary and its neighbours, and the half-way cases that round
    differently in two languages (9.995, 99.95, 999.5 of every rung)."""
    rng = random.Random(f"unit-ladders-{ladder.quantity}")
    values = [rng.choice((1, -1)) * 10 ** rng.uniform(-12, 30) for _ in range(4000)]
    for rung in ladder.rungs:
        for factor in (1, 1 - 1e-12, 1 + 1e-12, 1 - 1e-13, 0.999, 9.995, 99.95, 999.5, 9999.5, 1e6 + 0.5, 1e6 - 0.5, 2.5):
            values.append(rung.factor * factor)
            values.append(-rung.factor * factor)
        values.append((rung.factor if rung.start is None else rung.start) * (1 + 1e-9))
    values += [0.0, -0.0, 1e300, 1e-300, 5e-324, 1e-5, 4e-4, 9.99e-4, 1e-3]
    return values


@pytest.mark.skipif(shutil.which("node") is None, reason="Node isn't installed")
def test_the_browser_copy_gives_the_same_text_as_python():
    batches = {f"{ladder.quantity}|{system}": _check_values(ladder) for ladder, system in CASES}
    script = (
        "const m = await import(process.argv[1]);"
        "const batches = JSON.parse((await import('node:fs')).readFileSync(0, 'utf8'));"
        "const out = {};"
        "for (const [key, values] of Object.entries(batches)) {"
        "  const [name, system] = key.split('|');"
        "  out[key] = values.map(v => m.formatQuantity(m.LADDERS[name], v, { system }));"
        "}"
        "console.log(JSON.stringify(out));"
    )
    module_url = "file://" + os.path.abspath(os.path.join(_STATIC, "unitladder.js"))
    result = subprocess.run(
        ["node", "--input-type=module", "-e", script, module_url],
        input=json.dumps(batches), check=True, capture_output=True, text=True,
    )
    browser = json.loads(result.stdout)
    for ladder, system in CASES:
        key = f"{ladder.quantity}|{system}"
        wrong = [(value, format_quantity(ladder, value, system=system), shown)
                 for value, shown in zip(batches[key], browser[key])
                 if format_quantity(ladder, value, system=system) != shown]
        assert not wrong, f"{key}: {len(wrong)} differ, first {wrong[:5]}"


@pytest.mark.parametrize("ladder, value, auto, metric, customary", [
    (ladders.LENGTH, 100.0, "0.1 km", "0.1 km", "328 ft"),
    (ladders.LENGTH, 384.4e6, "384,400 km", "384,400 km", "238,855 mi"),
    (ladders.LENGTH, 1.496e11, "1 AU", "150 Gm", "1 AU"),
    (ladders.LENGTH, 9.46073e15, "30.7 cpc (1 ly)", "9.46 Pm", "30.7 cpc (1 ly)"),
    (ladders.SPEED, 30e3, "30 km/s", "30 km/s", "18.6 mi/s"),
    (ladders.SPEED, 1.0, "3.6 km/h", "3.6 km/h", "2.24 mph"),
    (ladders.TEMPERATURE, 288.0, "288 K (15 °C, 59 °F)", "288 K (15 °C)", "58.7 °F (15 °C, 288 K)"),
    (ladders.TEMPERATURE, 0.5, "500 mK (-273 °C, -459 °F)", "500 mK (-273 °C)", "-459 °F (-273 °C, 0.5 K)"),
    (ladders.PRESSURE, 610.0, "610 Pa (0.00602 atm, 0.0885 psi, 6.1 mbar)", "610 Pa (6.1 mbar)",
     "0.0885 psi (0.00602 atm)"),
    (ladders.PRESSURE, 101325.0, "101 kPa (1 atm, 14.7 psi)", "101 kPa", "14.7 psi (1 atm)"),
    (ladders.GRAVITY, 3.71, "3.71 m/s² (0.378 g)", "3.71 m/s²", "12.2 ft/s² (0.378 g)"),
    (ladders.MASS, constants.EARTH_MASS_TO_KG, "1 Earth mass", "5.97 × 10²⁴ kg", "1.32 × 10²⁵ lb"),
    (ladders.MASS, 100.0, "100 kg", "100 kg", "220 lb"),
    (ladders.POWER, constants.SOLAR_LUMINOSITY, "1 L☉", "383 YW", "1 L☉"),
    (ladders.MAGNETIC_FIELD, 1e8, "1.00 × 10¹² G", "1.00 × 10⁸ T", "1.00 × 10¹² G"),
])
def test_each_unit_preference_reads_as_decided(ladder, value, auto, metric, customary):
    """UX.78: Automatic is the ladders as they were; Metric only has no AU, ly, pc, Earth masses, L☉, atm,
    psi or °F; Customary puts ft, mi, mph, lb, °F, psi and ft/s² first and keeps the astronomical units."""
    assert format_quantity(ladder, value) == auto
    assert format_quantity(ladder, value, system="auto") == auto
    assert format_quantity(ladder, value, system="metric") == metric
    assert format_quantity(ladder, value, system="customary") == customary


def test_an_unknown_unit_preference_is_an_error():
    with pytest.raises(ValueError):
        format_quantity(ladders.LENGTH, 1.0, system="imperial")


def test_unitladders_js_is_up_to_date():
    result = subprocess.run([sys.executable, os.path.join(_ROOT, "scripts", "build_unit_ladders.py"), "--check"],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout


def test_a_value_under_a_thousandth_of_its_unit_is_scientific_not_a_row_of_zeros():
    assert format_quantity(ladders.LENGTH, 0.4) == "4.00 × 10⁻⁴ km"
    assert format_quantity(ladders.LENGTH, 0.5) == "5.00 × 10⁻⁴ km"
    assert format_quantity(ladders.LENGTH, 1.0) == "0.001 km"
    assert format_quantity(ladders.TIME, 2.5e-13) == "2.50 × 10⁻⁷ µs"
    assert format_quantity(ladders.LENGTH, 0) == "0 km"
    assert format_quantity(ladders.LENGTH, 1e-300).endswith("× 10⁻³⁰³ km")


@pytest.mark.skipif(shutil.which("node") is None, reason="Node isn't installed")
def test_the_spoken_form_spells_the_unit_out():
    script = (
        "const m = await import(process.argv[1]);"
        "console.log(JSON.stringify([m.formatDistanceKmSpoken(149597870.7), m.formatDistanceKmSpoken(384400),"
        " m.formatDistanceKmSpoken(3.0857e13), m.formatDistanceKmSpoken(null)]));"
    )
    module_url = "file://" + os.path.abspath(os.path.join(_STATIC, "unitladder.js"))
    out = subprocess.run(["node", "--input-type=module", "-e", script, module_url],
                         check=True, capture_output=True, text=True).stdout
    assert json.loads(out) == ["1 astronomical units", "384,400 kilometers", "1 parsecs (3.26 light-years)", "–"]


@pytest.mark.skipif(shutil.which("node") is None, reason="Node isn't installed")
def test_the_browser_body_radius_is_scientific_km_like_python():
    from planetgen.util.format import format_body_radius_km
    radii = [6371, 1737.4, 696_000, 0.5, 2.5e9]
    script = (
        "const m = await import(process.argv[1]);"
        "const values = JSON.parse((await import('node:fs')).readFileSync(0, 'utf8'));"
        "console.log(JSON.stringify(values.map(v => m.formatBodyRadiusKm(v))));"
    )
    module_url = "file://" + os.path.abspath(os.path.join(_STATIC, "unitladder.js"))
    out = subprocess.run(["node", "--input-type=module", "-e", script, module_url], input=json.dumps(radii),
                         check=True, capture_output=True, text=True).stdout
    python = [re.sub(r"<sup>(-?\d+)</sup>", lambda m: m.group(1).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")),
                     format_body_radius_km(km)) for km in radii]
    assert json.loads(out) == python
