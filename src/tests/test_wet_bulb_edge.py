# tests/test_wet_bulb_edge.py

"""
GEN.1586: the wet-bulb line the oxygen-water classes M, O and R read is the sourced uncompensable limit of `protection_limits` (31 C, [R4521]), not a margin under it. The classes read the
wet bulb of the hottest day (the daytime peak of the hottest band), so an edge below the limit puts a body whose air a person can live in on the edge: Earth's own median is 28 C.
"""

from planetgen.physics import class_fit, protection_limits


def _wet_bulb_lines(letter):
    return [c for cls in class_fit.CLASSES if cls.letter == letter and cls.kind == "t" for c in cls.when if c["quantity"] == "wet_bulb_c"]


def test_m_o_and_r_read_the_sourced_wet_bulb_limit():
    limit = protection_limits.limit("wet_bulb_c")
    for letter in ("M", "O", "R"):
        [line] = _wet_bulb_lines(letter)
        assert line.get("lt", line.get("min")) == limit, f"class {letter}: the wet-bulb line {line} is not the sourced limit {limit} C"


def test_an_earth_like_wet_bulb_is_inside_class_m():
    [line] = _wet_bulb_lines("M")
    assert 28.0 < line["lt"], "Earth's median hottest-day wet bulb (28 C) must not sit on M's edge"
