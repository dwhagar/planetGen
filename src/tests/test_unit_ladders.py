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
            return number * rung.factor, rung
    raise AssertionError(f"{unit!r} is on no rung of {ladder.quantity}")


@pytest.mark.parametrize("ladder", ladders.LADDERS, ids=lambda ladder: ladder.quantity)
class TestTable:
    def test_rungs_are_sorted_and_named(self, ladder):
        factors = [rung.factor for rung in ladder.rungs]
        starts = [rung.factor if rung.start is None else rung.start for rung in ladder.rungs]
        assert factors == sorted(factors) and len(set(factors)) == len(factors)
        assert starts == sorted(starts)
        assert len({rung.symbol for rung in ladder.rungs}) == len(ladder.rungs)
        assert all(rung.name for rung in ladder.rungs)

    def test_the_lowest_rung_takes_everything_below_it(self, ladder):
        assert ladder.rungs[0].start == 0.0

    def test_asides_are_complete(self, ladder):
        carriers = [rung for rung in ladder.rungs if rung.with_aside]
        assert bool(carriers) == bool(ladder.asides)
        if ladder.asides:
            sizes = [size for size, _ in ladder.asides]
            assert sizes == sorted(sizes, reverse=True) and sizes[-1] == 0.0
            for _, symbol in ladder.asides:
                assert ladder.rung(symbol)

    def test_zero_names_a_rung(self, ladder):
        if ladder.zero is not None:
            assert ladder.rung(ladder.zero)


def _astropy_factor(symbol):
    """How many SI base units one `symbol` is, by astropy; None when it has no unit of that name."""
    u = pytest.importorskip("astropy.units")
    from astropy import constants as const
    named = {"km": u.km, "AU": u.au, "mpc": u.mpc, "cpc": 0.01 * u.pc, "ly": u.lyr, "pc": u.pc, "kpc": u.kpc,
             "Mpc": u.Mpc, "Gpc": u.Gpc, "km/h": u.km / u.hr, "km/s": u.km / u.s, "Mm/s": u.Mm / u.s,
             "c": const.c, "µs": u.us, "ms": u.ms, "s": u.s, "minutes": u.min, "hours": u.hr, "days": u.day,
             "years": u.year, "ky": u.kyr, "My": u.Myr, "Gy": u.Gyr}
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


@pytest.mark.parametrize("ladder", ladders.LADDERS, ids=lambda ladder: ladder.quantity)
class TestProperties:
    @given(value=magnitudes)
    def test_the_text_reads_back_within_rounding(self, ladder, value):
        if value == 0:
            return
        read, rung = _primary(ladder, format_quantity(ladder, value))
        assert read == pytest.approx(value, rel=0.0051)
        assert rung is pick_rung(ladder, value)

    @given(first=magnitudes, second=magnitudes)
    def test_a_larger_size_never_picks_a_smaller_rung(self, ladder, first, second):
        low, high = sorted((abs(first), abs(second)))
        assert ladder.rungs.index(pick_rung(ladder, low)) <= ladder.rungs.index(pick_rung(ladder, high))

    @given(value=st.one_of(magnitudes, st.floats(), st.integers(), st.none(), st.sampled_from((0, 0.0, -0.0, 1e300, 10 ** 400))))
    def test_nothing_raises_and_only_non_numbers_are_a_dash(self, ladder, value):
        try:
            text = format_quantity(ladder, value)
        except OverflowError:
            assert isinstance(value, int)  # an int past float range is the caller's to guard (fmt.dash_unless_finite)
            return
        assert text
        if value is None or (isinstance(value, float) and math.isnan(value)):
            assert text == "–"
        else:
            assert "nan" not in text

    def test_each_rung_boundary_picks_that_rung(self, ladder):
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
    values += [0.0, -0.0, 1e300]
    return values


@pytest.mark.skipif(shutil.which("node") is None, reason="Node isn't installed")
def test_the_browser_copy_gives_the_same_text_as_python():
    batches = {ladder.quantity: _check_values(ladder) for ladder in ladders.LADDERS}
    script = (
        "const m = await import(process.argv[1]);"
        "const batches = JSON.parse((await import('node:fs')).readFileSync(0, 'utf8'));"
        "const out = {};"
        "for (const [name, values] of Object.entries(batches)) {"
        "  out[name] = values.map(v => m.formatQuantity(m.LADDERS[name], v));"
        "}"
        "console.log(JSON.stringify(out));"
    )
    module_url = "file://" + os.path.abspath(os.path.join(_STATIC, "unitladder.js"))
    result = subprocess.run(
        ["node", "--input-type=module", "-e", script, module_url],
        input=json.dumps(batches), check=True, capture_output=True, text=True,
    )
    browser = json.loads(result.stdout)
    for ladder in ladders.LADDERS:
        wrong = [(value, format_quantity(ladder, value), shown)
                 for value, shown in zip(batches[ladder.quantity], browser[ladder.quantity])
                 if format_quantity(ladder, value) != shown]
        assert not wrong, f"{ladder.quantity}: {len(wrong)} differ, first {wrong[:5]}"


def test_unitladders_js_is_up_to_date():
    result = subprocess.run([sys.executable, os.path.join(_ROOT, "scripts", "build_unit_ladders.py"), "--check"],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout
