"""
GEN.1525 (row 5 of the GEN.945 audit): the remnant mass ranges of `tuning` held against the data the repository holds. The neutron star ceiling is read against the maximum masses of
the 21 cold-matter equations of state (`physics/neutron_star_structure.py`, the integrated table committed beside it); the white dwarf range against the masses of the cooling
sequences (`datasets/white_dwarf_models`, skipped where the archive is not downloaded).
"""

import statistics

import pytest

from planetgen import tuning
from planetgen.datasets import fetch, neutron_star_matter, white_dwarf_models
from planetgen.physics import neutron_star_structure


def _maxima():
    return [neutron_star_structure.maximum_mass_solar(index) for index in range(1, neutron_star_matter.EOS_COUNT + 1)]


def test_the_neutron_star_ceiling_is_the_middle_of_the_maximum_masses_of_the_equations_of_state():
    maxima = _maxima()
    ceiling = tuning.NEUTRON_STAR_MASS_RANGE_SOLAR[1]
    assert min(maxima) < ceiling < max(maxima)                              # some equations of state support it and some do not: the ones that do are drawn for a star that heavy
    assert ceiling == pytest.approx(statistics.median(maxima), rel=0.02)    # 2.2 against the median 2.22 of the 21


def test_a_neutron_star_of_the_ceiling_has_an_equation_of_state_that_supports_it_and_the_floor_has_all_of_them():
    low, high = tuning.NEUTRON_STAR_MASS_RANGE_SOLAR
    assert neutron_star_structure.supporting(high) and len(neutron_star_structure.supporting(low)) == neutron_star_matter.EOS_COUNT


def test_the_white_dwarf_range_lies_inside_the_masses_of_the_cooling_sequences_to_the_scaled_oxygen_neon_edge():
    try:
        masses = white_dwarf_models.masses_sol()
    except fetch.DataMissing:
        pytest.skip("the white dwarf sequences are not downloaded")
    low, high = tuning.WD_MASS_RANGE_SOL
    assert min(masses) <= low and high <= max(masses) * 1.05       # the sequences end at 1.30; the 1.35 edge is read from the 1.30 sequence (registry WL7)
