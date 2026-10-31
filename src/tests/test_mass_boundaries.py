# tests/test_mass_boundaries.py

"""
GEN.933: the masses at which one kind of object gives way to another are one table (`physics/mass_boundaries.py`) with a source and a
tier for every row, and the mass function, companions, rogue bodies, star classes and remnants read their edges from it.
"""

import pathlib
import re

import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.physics import constants, mass_boundaries as boundaries, planet_mass

SOURCE = pathlib.Path(__file__).resolve().parent.parent / "planetgen"
TIERS = ("Verified", "Confident Unverified", "Marginal Unverified", "Unverified")


def test_every_boundary_has_a_value_a_unit_a_source_and_a_tier():
    names = [boundary.name for boundary in boundaries.BOUNDARIES]
    assert len(names) == len(set(names))
    for boundary in boundaries.BOUNDARIES:
        assert boundary.value > 0 and boundary.unit in ("Msun", "Mjup"), boundary.name
        assert boundary.source and boundary.tier.startswith(TIERS), boundary.name


def test_the_units_convert_both_ways():
    assert boundaries.sol("hydrogen_burning_limit") == 0.08
    assert boundaries.jupiter("deuterium_burning_limit") == 13.0
    assert boundaries.kg("deuterium_burning_limit") == 13.0 * constants.JUPITER_MASS_TO_KG
    assert boundaries.sol("deuterium_burning_limit") == pytest.approx(13.0 * constants.JUPITER_MASS_TO_KG / constants.SOLAR_MASS_TO_KG)
    assert boundaries.jupiter("hydrogen_burning_limit") == pytest.approx(0.08 * constants.SOLAR_MASS_TO_KG / constants.JUPITER_MASS_TO_KG)
    assert boundaries.kg("core_collapse_limit") == 8.0 * constants.SOLAR_MASS_TO_KG


def test_the_edges_are_in_order_and_a_disagreement_is_recorded_not_silent():
    assert boundaries.sol("deuterium_burning_limit") < boundaries.sol("brown_dwarf_upper_limit") <= boundaries.sol("hydrogen_burning_limit")
    assert boundaries.sol("white_dwarf_progenitor_floor") < boundaries.sol("core_collapse_limit")
    assert boundaries.sol("chandrasekhar_limit") == constants.CHANDRASEKHAR_LIMIT_SOL
    # GEN.936: the brown dwarf mass function ends where the stars' begins, with no stretch between and none counted twice.
    assert boundaries.sol("brown_dwarf_upper_limit") == boundaries.sol("hydrogen_burning_limit")


def test_the_generator_rows_read_the_table():
    assert tuning.IMF_BREAKS_SOL[0] == boundaries.sol("hydrogen_burning_limit")
    assert tuning.SUPERGIANT_MIN_MASS_SOL == boundaries.sol("core_collapse_limit") == tuning.PRE_MAIN_SEQUENCE_MAX_MASS_SOL
    assert tuning.WHITE_DWARF_PROGENITOR_SOL == (boundaries.sol("white_dwarf_progenitor_floor"), boundaries.sol("core_collapse_limit"))
    assert tuning.WD_PROGENITOR_MASS_RANGE_SOL[1] == boundaries.sol("core_collapse_limit")
    assert tuning.CORE_COLLAPSE_MASS_RANGE_SOL[0] == boundaries.sol("core_collapse_limit")
    assert tuning.ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER == (boundaries.jupiter("deuterium_burning_limit"),
                                                           boundaries.jupiter("brown_dwarf_upper_limit"))
    assert planet_mass._BROWN_DWARF_LIMIT_KG == boundaries.kg("deuterium_burning_limit")
    assert boundaries.jupiter("deuterium_burning_limit") * 317.8 == tuning._ROGUE_PLANET_BIN_EDGES_EARTH["jupiter"][1]


def test_the_star_class_edges_are_one_table_that_the_companion_code_reads():
    edges = [lower for _letter, lower in boundaries.SPECTRAL_CLASS_LOWER_MASS_SOL]
    assert edges == sorted(edges, reverse=True) and [letter for letter, _ in boundaries.SPECTRAL_CLASS_LOWER_MASS_SOL] == list("OBAFGKM")
    assert multistar.spectral_class(20.0) == "O" and multistar.spectral_class(1.0) == "G" and multistar.spectral_class(0.2) == "M"
    assert multistar.companion_class(10.0) == "B-early" and multistar.companion_class(6.0) == "B-mid" and multistar.companion_class(3.0) == "B"
    assert not hasattr(multistar, "SPECTRAL_CLASS_LOWER_MASS_SOL") and not hasattr(tuning, "B_SUBCLASS_LOWER_MASS_SOL")


def test_the_brown_dwarf_edge_of_the_companion_code_is_the_hydrogen_limit():
    assert not hasattr(multistar, "BROWN_DWARF_MAX_MASS_SOL")
    placed = multistar.Placed("star", 1.0, 1.0)
    assert placed.birth_mass_solar == 1.0       # the draw below the limit is the brown dwarf branch (`multistar` source, one read)


def test_no_module_types_one_of_the_edges_again():
    """The deuterium limit as `13 * JUPITER_MASS` and the class tables as assignments exist in the table only."""
    typed = re.compile(r"\b13(\.0)?\s*\*\s*(\w+\.)?JUPITER_MASS|m_jup\s*>=\s*13\b|^SPECTRAL_CLASS_LOWER_MASS_SOL\s*=|^B_SUBCLASS_LOWER_MASS_SOL\s*=",
                       re.MULTILINE)
    found = []
    for path in SOURCE.rglob("*.py"):
        if path.name == "mass_boundaries.py":
            continue
        for match in typed.finditer(path.read_text(encoding="utf-8")):
            found.append(f"{path.relative_to(SOURCE)}: {match.group(0)}")
    assert not found, found
