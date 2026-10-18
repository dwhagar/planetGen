"""
The two rounded unit factors are exact (GEN.572): the light-years in a parsec and the astronomical units in a light-year
that `physics/hii.py`, `generation/phenomena/asteroid_field.py` and `generation/phenomena_plausibility.py` used to type
as 3.26156 and 63241.077 are now the ratios of the definitions in `physics/constants.py`.
"""

import math

from planetgen.generation import phenomena_plausibility
from planetgen.generation.phenomena import asteroid_field
from planetgen.physics import constants as pc, hii


def test_the_three_sites_read_the_exact_factors():
    assert hii.LY_PER_PC == pc.PARSEC_M / pc.LIGHTYEAR_M
    assert asteroid_field.AU_PER_LY == pc.LIGHTYEAR_M / pc.AU_M
    assert phenomena_plausibility._AU_PER_LY == pc.LIGHTYEAR_M / pc.AU_M


def test_the_factors_are_not_the_old_rounded_values():
    assert math.isclose(hii.LY_PER_PC, 3.2615637771674333, rel_tol=1e-15)
    assert math.isclose(asteroid_field.AU_PER_LY, 63241.07708426628, rel_tol=1e-15)
    assert hii.LY_PER_PC != 3.26156 and asteroid_field.AU_PER_LY != 63241.077


def test_the_factors_agree_with_the_unit_helpers():
    from planetgen.physics import units
    assert math.isclose(hii.LY_PER_PC, units.pc_to_ly(1.0), rel_tol=1e-15)
