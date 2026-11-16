# tests/test_white_dwarf_pollution_rate.py

"""
GEN.1231: a polluted white dwarf has a rate at which rock falls on it. The flag stays the drawn one (25 to 50 percent of the cool ones); the rate is the belt's own when the mass its star threw inward is
enough to pollute a white dwarf, else the typical rate (`belt_history.TYPICAL_POLLUTION_RATE_G_S`, a stand-in). The mass a belt's star threw inward is stored with the belt (`asteroid_belts.thrown_inward_kg`).
"""

from planetgen.physics import belt_history


def test_a_belt_that_threw_enough_inward_gives_its_own_rate():
    thrown_kg = 1.0e21
    rate = belt_history.rate_on_white_dwarf(thrown_kg)
    assert rate == belt_history.pollution_rate_g_s(thrown_kg)
    assert belt_history.is_polluted(rate) and rate != belt_history.TYPICAL_POLLUTION_RATE_G_S


def test_no_belt_or_a_belt_that_gives_too_small_a_rate_leaves_the_typical_rate():
    assert belt_history.rate_on_white_dwarf(None) == belt_history.TYPICAL_POLLUTION_RATE_G_S
    assert belt_history.rate_on_white_dwarf(0.0) == belt_history.TYPICAL_POLLUTION_RATE_G_S
    tiny = 1.0e10                       # kg: over the delivery time this is far under the lowest rate seen on a polluted white dwarf
    assert not belt_history.is_polluted(belt_history.pollution_rate_g_s(tiny))
    assert belt_history.rate_on_white_dwarf(tiny) == belt_history.TYPICAL_POLLUTION_RATE_G_S


def test_the_typical_rate_is_a_rate_that_pollutes():
    assert belt_history.is_polluted(belt_history.TYPICAL_POLLUTION_RATE_G_S)
