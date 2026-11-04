# tests/test_hycean_candidate.py

"""The Hycean box as a derived flag (GEN.963)."""

import math

from planetgen.physics import constants, hycean_candidate as hc, hydrosphere, planet_classes


def test_the_density_of_known_planets():
    assert math.isclose(hc.density_g_cm3(1.0, 1.0), 5.51, abs_tol=0.02)         # the Earth
    assert math.isclose(hc.density_g_cm3(17.15, 3.88), 1.64, abs_tol=0.05)      # Neptune


def test_a_sub_neptune_with_a_temperate_envelope_is_a_candidate():
    # 4.8 Earth masses, 2.1 Earth radii, 380 K (the size and temperature of the temperate sub-Neptunes the box was drawn around; density 2.9)
    assert hc.density_g_cm3(4.8, 2.1) < 3.0
    assert hc.hycean_candidate("T", 4.8, 2.1, 380.0)


def test_each_range_excludes():
    assert not hc.in_box(0.5, 1.5, 300.0)                 # too light
    assert not hc.in_box(12.0, 2.5, 300.0)                # too heavy
    assert not hc.in_box(5.0, 1.0, 300.0)                 # too small (a dense rocky world: 20 g/cc)
    assert not hc.in_box(8.2, 2.74, 300.0)                # GJ 1214 b's size: over 2.6 Earth radii
    assert not hc.in_box(5.0, 2.0, 700.0)                 # too hot
    assert not hc.in_box(5.0, 2.0, 100.0)                 # too cold
    assert not hc.in_box(9.0, 1.45, 300.0)                # inside the mass, radius and temperature ranges but a density of 16 g/cc
    assert not hc.in_box(2.0, 2.6, 300.0)                 # a density of 0.6 g/cc


def test_the_flag_sits_on_the_classes_t_i_and_v_only():
    assert hc.hycean_candidate("V", 3.5, 1.9, 300.0)
    assert not hc.hycean_candidate("S", 3.5, 1.9, 300.0)
    assert not hc.hycean_candidate("P", 3.5, 1.9, 300.0)


def test_the_flag_agrees_with_the_classes_the_code_has():
    letters = {row["letter"] for row in planet_classes.DOCUMENT["rows"]}
    assert set(hc.HYCEAN_CLASS_LETTERS) <= letters
    # a gas dwarf of 5 Earth masses is class T by mass alone: the flag adds nothing to the class
    values = {"body_type": "g", "mass_kg": 5.0 * constants.EARTH_MASS_TO_KG}
    assert planet_classes.classify_values(values).letter == "T"
    assert hydrosphere.HYCEAN_MIN_H2_SHARE == 0.5            # the state the flag is compared with is the hydrosphere's
