"""
GEN.1557: the words of the gas dwarf (T) and ice giant (I) classes agree with the mix their bodies store: mostly rock and metal under an envelope for a gas dwarf, an interior of rock, metal and ices
for an ice giant, and the day-counter line for a body with no surface says the rock lies deep under the envelope.
"""

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import planet_class_text, planet_classes
from planetgen.util import draw
from planetgen.web.lib import day_counter
from planetgen.physics import day_night


def words_of(letter):
    rows = [row for row in planet_classes.CLASSES if row.letter == letter]
    return [row.phrase for row in rows] + [row.sentence for row in rows] + [planet_class_text.description(letter)]


def test_a_gas_dwarf_is_not_described_as_a_small_core_under_an_envelope():
    for text in words_of("T"):
        assert "small core" not in text, text
    assert all("rock" in text for text in words_of("T"))


def test_an_ice_giant_is_not_described_as_mostly_ices():
    for text in words_of("I"):
        assert "mostly water" not in text and "mostly ices" not in text, text
    assert all("rock" in text for text in words_of("I"))


def test_the_no_surface_line_says_any_rock_lies_under_the_envelope():
    notes = next(value for value in vars(day_counter).values() if isinstance(value, dict) and day_night.NO_SURFACE in value)
    line = notes[day_night.NO_SURFACE]
    assert "rock lies deep" in line


def test_the_stored_gas_dwarfs_are_mostly_rock_and_metal_under_a_small_share_of_envelope():
    """The claim of the class T sentence, on the bodies a generated system stores (seed 2001 has several)."""
    config = SystemConfig()
    config.STAR_TYPE, config.PLANETS, config.COMETS, config.BINARY_SYSTEM = "G2V", True, False, False
    with draw.bound(2001):
        system = StarSystem(system_config=config)
    dwarfs = [p for p in system.planets if getattr(p, "planet_class", "") == "T"]
    assert dwarfs, "the seed no longer makes a gas dwarf; pick another"
    for body in dwarfs:
        solids = (body.rock_mass_fraction or 0.0) + (body.metal_mass_fraction or 0.0)
        assert solids > 0.5 > (body.envelope_mass_fraction or 0.0), (solids, body.envelope_mass_fraction)
