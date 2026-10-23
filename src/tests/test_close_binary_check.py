# tests/test_close_binary_check.py -- TEST.196: the generated close-binary shares against the published ranges.

import pytest

from tests.audit import close_binary
from tests.audit.close_binary import Bin

TOTAL = {}


@pytest.fixture(scope="module")
def bins():
    return close_binary.run_check()


def test_every_bin_with_a_published_range_is_inside_it(bins):
    found = close_binary.flags(bins)
    assert not found, "\n" + close_binary.format_report(bins, found)


def test_the_report_prints_every_bin_and_the_g_dwarf_level(bins):
    text = close_binary.format_report(bins)
    for name, _mass, _letter in close_binary.MASS_BINS:
        assert name in text
    assert "G dwarf level at [Fe/H] 0" in text and "never fails" in text
    assert text.count("\n") > len(close_binary.FE_H_BINS)


def test_the_share_inside_10_au_falls_with_metallicity_for_a_solar_type_star(bins):
    shares = [bins[("G dwarf", fe_h)].share for fe_h in close_binary.FE_H_BINS]
    assert shares == sorted(shares, reverse=True)


def test_the_check_flags_a_bin_outside_its_range():
    ok = {key: Bin(key[0], 1.0, key[1], 1000, 200, 100, 100) for key in
          ((n, f) for n, _m, _l in close_binary.MASS_BINS for f in close_binary.FE_H_BINS)}
    assert any(flag.row.startswith("G dwarf ratio") for flag in close_binary.flags(ok))      # a flat G dwarf line
    broken = dict(ok)
    broken[("O star", 0.0)] = Bin("O star", 28.0, 0.0, 1000, 20, 20, 0)                       # 2%: under the O span
    assert any(flag.row == "O share inside 10 AU" for flag in close_binary.flags(broken))
    broken[("mid B", 0.0)] = Bin("mid B", 7.0, 0.0, 1000, 20, 20, 0)
    assert any(flag.row == "mid B share inside 10 AU" for flag in close_binary.flags(broken))
    skewed = dict(ok)
    skewed[("G dwarf", 0.5)] = Bin("G dwarf", 1.0, 0.5, 1000, 200, 100, 600)                 # an outer share that moved
    assert any(flag.row == "G dwarf outer share" for flag in close_binary.flags(skewed))
