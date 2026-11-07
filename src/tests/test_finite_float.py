# tests/test_finite_float.py

"""SEC.41: `checks.finite_float` is the one place a request's text becomes a float, and it refuses what is not finite."""

import pytest

from planetgen.util.checks import finite_float


def test_a_finite_number_in_any_spelling_is_a_float():
    assert finite_float("1.5") == 1.5 and finite_float(" -2 ") == -2.0 and finite_float("1e3") == 1000.0 and finite_float(4) == 4.0


@pytest.mark.parametrize("raw", ["nan", "NaN", "-nan", "inf", "-inf", "+Infinity", "-infinity"])
def test_nan_and_infinity_are_refused_as_a_value_error(raw):
    with pytest.raises(ValueError, match="finite"):
        finite_float(raw)


@pytest.mark.parametrize("raw", ["", "soon", "1,5", None])
def test_text_that_is_no_number_is_a_value_or_type_error(raw):
    with pytest.raises((ValueError, TypeError)):
        finite_float(raw)
