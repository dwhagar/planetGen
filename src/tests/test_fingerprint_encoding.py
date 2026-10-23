# tests/test_fingerprint_encoding.py

"""
GEN.304: the canonical encoder every fingerprint hashes (`planetgen.util.fingerprint_encoding`): the same value gives
the same text on any machine, a 1-ulp change in a float leaves the digest unchanged and a change in the ninth digit
changes it, and the structural digest ignores float noise altogether.
"""

import decimal
import math

import pytest

from planetgen.util import fingerprint_encoding as encoding


def test_a_value_has_one_text():
    value = {"b": [1, 2.5, None, True, False], "a": "x", "c": b"\x01\xff"}
    assert encoding.encode(value) == '{"a":"x","b":[1,2.50000000e+00,null,true,false],"c":b"01ff"}'
    assert encoding.encode({"a": 1, "b": 2}) == encoding.encode({"b": 2, "a": 1})
    assert encoding.encode([1, 2]) != encoding.encode([2, 1])
    assert encoding.encode((1, 2)) == encoding.encode([1, 2])


def test_the_types_do_not_blur():
    texts = [encoding.encode(value) for value in (1, 1.0, "1", True, None, b"1", [1], {"1": 1})]
    assert len(set(texts)) == len(texts)


@pytest.mark.parametrize("value, text", [
    (0.1 + 0.2, "3.00000000e-01"),
    (-0.0, "0.00000000e+00"),
    (0.0, "0.00000000e+00"),
    (1234567890.123, "1.23456789e+09"),
    (1.5e-300, "1.50000000e-300"),
    (6.02214076e23, "6.02214076e+23"),
    (float("nan"), "nan"),
    (float("inf"), "inf"),
    (float("-inf"), "-inf"),
    (decimal.Decimal("1.500"), "1.50000000e+00"),
    (decimal.Decimal("-0"), "0.00000000e+00"),
])
def test_floats_are_nine_significant_digits(value, text):
    assert encoding.float_token(value) == text


def test_integers_are_exact():
    assert encoding.int_token(2 ** 53) == str(2 ** 53)
    assert encoding.int_token(2 ** 53 + 1) == hex(2 ** 53 + 1)
    assert encoding.int_token(-(2 ** 53 + 1)) == "-" + hex(2 ** 53 + 1)
    assert encoding.encode(2 ** 95 + 12345) == hex(2 ** 95 + 12345)
    assert encoding.encode(2 ** 53 + 1) != encoding.encode(2 ** 53 + 2)


def test_strings_are_not_normalised():
    composed, decomposed = "é", "é"
    assert encoding.encode(composed) != encoding.encode(decomposed)
    assert encoding.encode(composed) == '"\\u00e9"'
    assert encoding.encode('a"b\\c\n') == '"a\\"b\\\\c\\n"'
    assert encoding.encode("\U0001f600") == '"\\U0001f600"'
    assert encoding.encode("\ud800") == '"\\ud800"'


def test_a_key_must_be_a_string_and_a_value_a_known_type():
    with pytest.raises(TypeError):
        encoding.encode({1: 2})
    with pytest.raises(TypeError):
        encoding.encode({"a": {1, 2}})


def test_a_last_place_difference_leaves_the_digest_unchanged():
    value = 3.0800000123456
    noisy = math.nextafter(value, math.inf)
    assert noisy != value
    for kind in ("full", "structural", "floats"):
        left = encoding.digests({"mass": [value, 1.0]})
        right = encoding.digests({"mass": [noisy, 1.0]})
        assert getattr(left, kind) == getattr(right, kind)


def test_a_change_in_the_ninth_digit_changes_the_float_digest_only():
    left = encoding.digests({"mass": 3.08000001, "n": 4})
    right = encoding.digests({"mass": 3.08000002, "n": 4})
    assert left.full != right.full and left.floats != right.floats
    assert left.structural == right.structural


def test_a_change_of_structure_changes_the_structural_digest():
    base = encoding.digests({"a": 1.0, "b": "x"})
    for other in ({"a": 1.0, "b": "y"}, {"a": 1.0, "b": "x", "c": None}, {"a": [1.0], "b": "x"}, {"a": 1, "b": "x"}):
        assert encoding.digests(other).structural != base.structural


def test_the_digests_are_pinned():
    """Golden vectors: the same value gives the same digest on Python 3.12 and 3.13 and on every machine; a change
    here is a new format (`FORMAT_TAG`), never an accident."""
    value = {"name": "Sol é", "n": 2 ** 60, "xs": [1.0, 0.1, -2.5e-9, 123456789.123456789], "raw": b"\x00\x01"}
    assert encoding.encode(value) == (
        '{"n":0x1000000000000000,"name":"Sol \\u00e9","raw":b"0001",'
        '"xs":[1.00000000e+00,1.00000000e-01,-2.50000000e-09,1.23456789e+08]}')
    digests = encoding.digests(value)
    assert digests.full == PINNED["full"]
    assert digests.structural == PINNED["structural"]
    assert digests.floats == PINNED["floats"]


PINNED = {
    "full": "8833503189a72b57b21d25b667090769bb8e4188c1ce35dedc411ffb571e0441",
    "structural": "fc7700eceed9279c9fc9d69488264521838b584223f340d9aa69517ba76e30ba",
    "floats": "810090f55af4efd9b5c022553c58b9660410e0b67ce925e56bbbbcb7e3dc5e59",
}


def test_the_tag_and_the_kind_are_part_of_the_digest():
    assert encoding.text_digest("full", "x") != encoding.text_digest("structural", "x")
    with pytest.raises(ValueError):
        encoding.text_digest("other", "x")


def test_rows_with_no_order_digest_alike():
    rows = [{"t": "a", "v": 1.0}, {"t": "b", "v": 2.0}, {"t": "a", "v": 0.5}]
    assert encoding.row_set_digests(rows) == encoding.row_set_digests(list(reversed(rows)))
    noisy = [{"t": "a", "v": math.nextafter(1.0, 2.0)}, {"t": "b", "v": 2.0}, {"t": "a", "v": 0.5}]
    assert encoding.row_set_digests(noisy) == encoding.row_set_digests(rows)
    assert encoding.row_set_digests(rows[:2]) != encoding.row_set_digests(rows)
